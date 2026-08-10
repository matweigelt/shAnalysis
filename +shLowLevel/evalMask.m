function [mask, grid] = evalMask(idx, region, opts)
%EVALMASK Evaluate a region definition on the toolbox quadrature grid.
%
%   [MASK, GRID] = shLowLevel.evalMask(IDX, REGION) evaluates REGION on the
%   Gauss-Legendre quadrature grid of shLowLevel.synthesisMatrix(IDX) and
%   returns MASK (Ngrid x 1, in [0,1]) plus the GRID struct. Shared by
%   shLowLevel.basinKernel and shLowLevel.slepianBasis.
%
%   REGION forms:
%     function handle  m = f(latDeg, lonDeg), vectorized, values in [0,1]
%     K x 2 double     polygon vertices [latDeg lonDeg], lon in [0, 360);
%                      point-in-polygon via inpolygon (lon treated
%                      planar - split polygons crossing lon = 0 yourself)
%     Ngrid x 1 double mask already on the quadrature grid (validated)
%
%   Options
%     BufferKm (0)  >0 grows, <0 shrinks the region by great-circle
%                   distance [km] to the region boundary points, computed
%                   exactly in latitude chunks (O(N*Nregion); fine to
%                   Lmax ~ 60-96, cost noted in the doc page)
%     R (6378136.3) sphere radius [m] for the km conversion
%     OverSample (2)  see arguments block
%   Outputs
%     mask  (Ngrid,1) double in [0,1]
%     grid  struct from shLowLevel.synthesisMatrix (latDeg, lonDeg per ring)
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    idx (1,1) struct
    region
    opts.BufferKm (1,1) double = 0
    opts.R (1,1) double = 6378136.3
    opts.OverSample (1,1) double {mustBeInteger, mustBePositive} = 2
end
% OverSample refines the quadrature grid (v2.2.1): the Gauss rule is
% exact for the band-limited Y products at any NLat >= Lmax+1, but a
% region INDICATOR is not band-limited - its staircase quantization on
% the default Lmax+1 rings reaches ~18% of a small cap's area. Factor 2
% is the accuracy/cost default; callers must pass the SAME factor to
% their own synthesisMatrix call (basinKernel/slepianBasis do).
[~, ~, grid] = shLowLevel.synthesisMatrix(idx, ...
    NLat = opts.OverSample * (idx.Lmax + 1), ...
    NLon = opts.OverSample * (2 * idx.Lmax + 2));
latRow = repelem(grid.latDeg(:), numel(grid.lonDeg));
lonRow = repmat(grid.lonDeg(:), numel(grid.latDeg), 1);
N = numel(latRow);
if isa(region, 'function_handle')
    mask = double(region(latRow, lonRow));
elseif isnumeric(region) && size(region, 2) == 2 && size(region, 1) >= 3
    mask = double(inpolygon(lonRow, latRow, region(:,2), region(:,1)));
elseif isnumeric(region) && numel(region) == N
    mask = double(region(:));
else
    error('shLowLevel:evalMask:badRegion', ...
        'REGION must be a function handle, a Kx2 [lat lon] polygon, or an Ngrid x 1 mask.');
end
mask = min(max(mask, 0), 1);

if opts.BufferKm ~= 0
    inR = mask > 0.5;
    src = find(xor(opts.BufferKm > 0, ~inR));    % grow: from inside pts
    if opts.BufferKm < 0, src = find(~inR); end  % shrink: dist to outside
    if opts.BufferKm > 0, src = find(inR); end
    dLim = abs(opts.BufferKm) * 1e3 / opts.R;    % central angle [rad]
    la1 = deg2rad(latRow); lo1 = deg2rad(lonRow);
    la2 = deg2rad(latRow(src))'; lo2 = deg2rad(lonRow(src))';
    hit = false(N, 1);
    chunk = 2000;
    for i0 = 1:chunk:N
        ii = i0:min(i0+chunk-1, N);
        cosd_ = sin(la1(ii)).*sin(la2) + cos(la1(ii)).*cos(la2).*cos(lo1(ii) - lo2);
        hit(ii) = any(acos(min(max(cosd_, -1), 1)) <= dLim, 2);
    end
    if opts.BufferKm > 0
        mask = double(inR | hit);                % grow by dLim
    else
        mask = double(inR & ~hit);               % shrink by dLim
    end
end
end
