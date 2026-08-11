function [rmsVal, info] = oceanRMS(grid, latDeg, lonDeg, ocean, opts)
%OCEANRMS Open-ocean RMS of a field - the standard GRACE noise metric.
%
%   RMSVAL = shLowLevel.oceanRMS(GRID, LATDEG, LONDEG, OCEAN) returns the
%   area-weighted RMS of GRID over the OPEN ocean: ocean points more than
%   MinDistanceKm (default 1000) from the nearest non-ocean point.
%
%   Far from land a GRACE/GRACE-FO field should contain almost no real
%   signal, so what remains there is error. This makes the open-ocean RMS
%   the standard noise metric of the field - processing centres quote it
%   to compare filters and solutions (Dahle et al. 2025), and it is the
%   natural value for shLowLevel.leakageCorrect's NoiseLevel, which needs
%   the noise in the units of the field being corrected.
%
%   OCEAN is USER-SUPPLIED, like Love numbers: this toolbox does not
%   assume a coastline. Base MATLAB's 'coastlines' data set is not
%   present in every installation, and a wrong ocean mask silently
%   changes the number, so it has to be a decision the caller makes.
%
%   Inputs
%     grid       (nLat x nLon) double   the field to measure
%     latDeg     (1 x nLat | nLat x 1) double  geocentric latitudes [deg]
%     lonDeg     (1 x nLon | nLon x 1) double  longitudes [deg]
%     ocean      (nLat x nLon) logical, or a function handle
%                f(latDeg, lonDeg) -> logical/[0,1] evaluated on the grid.
%                true = ocean
%
%   Options
%     MinDistanceKm (1000)  erode the ocean mask by this great-circle
%             distance, i.e. keep only points farther than this from any
%             non-ocean point. 1000 km is the common convention; 0 uses
%             the whole ocean and will be contaminated by leakage from
%             the continents
%     Weighted (true)  weight by cos(latitude). Grid cells are not equal
%             in area, and an unweighted RMS over a lat/lon grid
%             over-counts the polar rows. On white noise the two agree,
%             which is exactly why the difference is easy to miss until
%             the field has structure
%     Subsample (1)  use every n-th boundary point for the distance
%             computation. The erosion costs O(nOcean x nLand); on fine
%             grids raise this, at the cost of a distance error of order
%             n grid steps
%
%   Outputs
%     rmsVal     (1,1) double  open-ocean RMS, in GRID's units. NaN when
%                the erosion leaves no points
%     info       (1,1) struct  fields: nPixels (1,1 double, points used),
%                fractionUsed (1,1 double, of the whole grid),
%                mask (nLat x nLon logical, the eroded region - plot it
%                once before trusting the number), minDistanceKm,
%                weighted (1,1 logical)
%
%   Example
%     isOcean = @(la, lo) ~inpolygon(lo, la, landLon, landLat);
%     noise = shLowLevel.oceanRMS(ewh, -89:89, 0:359, isOcean);
%     m = shLowLevel.leakageCorrect(ewh, -89:89, 0:359, ...
%             Filter = "gauss300", Mask = basin, NoiseLevel = noise);
%
%   See also shLowLevel.leakageCorrect, shLowLevel.spatialStats.
%
%   Developed by Matthias Weigelt with the help of Claude (Opus 5),
%   2026-08-11 (v3.5.1).
arguments
    grid double
    latDeg double
    lonDeg double
    ocean
    opts.MinDistanceKm (1,1) double {mustBeNonnegative} = 1000
    opts.Weighted (1,1) logical = true
    opts.Subsample (1,1) double {mustBeInteger, mustBePositive} = 1
end
latDeg = latDeg(:).';
lonDeg = lonDeg(:).';
nLat = numel(latDeg);
nLon = numel(lonDeg);
assert(isequal(size(grid), [nLat, nLon]), 'shLowLevel:oceanRMS:badSize', ...
    'grid must be nLat x nLon (%d x %d), got %d x %d.', ...
    nLat, nLon, size(grid, 1), size(grid, 2));
[LO, LA] = meshgrid(lonDeg, latDeg);
if isa(ocean, 'function_handle')
    isOc = logical(ocean(LA, LO));
else
    assert(isequal(size(ocean), [nLat, nLon]), ...
        'shLowLevel:oceanRMS:badMask', ...
        'ocean must be a function handle or an nLat x nLon mask.');
    isOc = logical(ocean);
end
assert(any(isOc(:)), 'shLowLevel:oceanRMS:emptyOcean', ...
    'The ocean mask selects no points.');

keep = erodeMask(isOc, LA, LO, opts.MinDistanceKm, opts.Subsample);
n = nnz(keep);
w = cosd(LA);
if n == 0
    rmsVal = NaN;
else
    v = grid(keep);
    if opts.Weighted
        ww = w(keep);
        rmsVal = sqrt(sum(ww .* v.^2) / sum(ww));
    else
        rmsVal = sqrt(mean(v.^2));
    end
end
info = struct('nPixels', n, 'fractionUsed', n / numel(grid), ...
    'mask', keep, 'minDistanceKm', opts.MinDistanceKm, ...
    'weighted', opts.Weighted);
end

% ------------------------------------------------------------- helpers
function keep = erodeMask(isOc, LA, LO, dkm, sub)
%ERODEMASK Ocean points farther than dkm from the nearest non-ocean point.
keep = isOc;
if dkm <= 0 || all(isOc(:))
    return
end
R = 6371;
la = deg2rad(LA(isOc));
lo = deg2rad(LO(isOc));
nla = deg2rad(LA(~isOc));
nlo = deg2rad(LO(~isOc));
nla = nla(1:sub:end);
nlo = nlo(1:sub:end);
dmin = inf(numel(la), 1);
% chunk over the boundary points: the full nOcean x nLand matrix is
% large enough to matter on fine grids
for c = 1:500:numel(nla)
    e = min(c + 499, numel(nla));
    cosPsi = sin(la) .* sin(nla(c:e)).' + ...
        cos(la) .* cos(nla(c:e)).' .* cos(lo - nlo(c:e).');
    d = R * acos(min(1, max(-1, cosPsi)));
    dmin = min(dmin, min(d, [], 2));
end
keep(isOc) = dmin > dkm;
end
