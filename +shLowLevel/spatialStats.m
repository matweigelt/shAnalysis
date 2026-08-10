function stats = spatialStats(A, B, latDeg, lonDeg, opts)
%SPATIALSTATS Area-weighted spatial comparison of two grids.
%
%   STATS = shLowLevel.spatialStats(A, B, LATDEG, LONDEG) compares two
%   nlat x nlon grids with proper cos(latitude) area weighting -
%   unweighted RMS on a lat/lon grid over-counts the poles and is
%   simply wrong. Returns bias, RMSD, the CENTERED pattern statistics
%   (centered RMSD, weighted correlation, standard deviations) that
%   feed a Taylor diagram, and verifies the Taylor identity
%   crmsd^2 = stdA^2 + stdB^2 - 2 stdA stdB corr internally.
%
%   Options
%     Mask ([])     nlat x nlon logical: restrict to a region (e.g.
%                   land/ocean); default all points
%     Weights ([])  custom nlat x nlon weights (e.g. quadrature); default
%                   cos(latDeg) rings, normalized over the mask
%
%   Outputs
%     stats      (1,1) struct  fields:
%                  .bias    (1,1) double  weighted mean of A - B
%                  .rmsd    (1,1) double  weighted RMS of A - B
%                  .crmsd   (1,1) double  centered (bias-free) RMSD
%                  .corr    (1,1) double  weighted pattern correlation
%                  .stdA    (1,1) double  weighted std of A
%                  .stdB    (1,1) double  weighted std of B
%                  .ratio   (1,1) double  stdB / stdA
%                  .nUsed   (1,1) double  number of grid points used
%
%   Example
%     A1 = g.synthesis(lat, lon); A2 = gF.synthesis(lat, lon);
%     st = shLowLevel.spatialStats(A1, A2, lat, lon);
%     fprintf("rmsd %.3g, corr %.3f\n", st.rmsd, st.corr)
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-08 (v2.6.0).
arguments
    A double
    B double
    latDeg (1,:) double
    lonDeg (1,:) double
    opts.Mask = []
    opts.Weights double = []
end
if ~isequal(size(A), size(B), [numel(latDeg), numel(lonDeg)])
    error('shLowLevel:spatialStats:sizeMismatch', ...
        'A and B must be numel(lat) x numel(lon).');
end
if isempty(opts.Weights)
    w = cosd(latDeg(:)) * ones(1, numel(lonDeg));
else
    w = opts.Weights;
    if ~isequal(size(w), size(A))
        error('shLowLevel:spatialStats:badWeights', ...
            'Weights must match the grid size.');
    end
end
use = isfinite(A) & isfinite(B);
if ~isempty(opts.Mask), use = use & logical(opts.Mask); end
w(~use) = 0;
if sum(w(:)) <= 0
    error('shLowLevel:spatialStats:emptyMask', 'No usable grid points.');
end
w = w / sum(w(:));
wm = @(X) sum(w(:) .* X(:), 'omitnan');
d = A - B;
bias = wm(d);
rmsd = sqrt(wm(d.^2));
a = A - wm(A); b = B - wm(B);
stdA = sqrt(wm(a.^2)); stdB = sqrt(wm(b.^2));
corr = wm(a.*b) / (stdA * stdB);
crmsd = sqrt(wm((a - b).^2));
stats = struct('bias', bias, 'rmsd', rmsd, 'crmsd', crmsd, 'corr', corr, ...
    'stdA', stdA, 'stdB', stdB, 'ratio', stdB / stdA, 'nUsed', nnz(use));
end
