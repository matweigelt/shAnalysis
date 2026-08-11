function [k, info] = gridScaling(model, latDeg, lonDeg, opts)
%GRIDSCALING Per-pixel leakage/attenuation scaling factors from a model.
%
%   [K, INFO] = shLowLevel.gridScaling(MODEL, LATDEG, LONDEG, Filter = ...)
%   computes the classic gridded scaling factors: push a MODEL field
%   series through the SAME processing chain the data saw and regress,
%   per pixel, the true value on the filtered one,
%
%       k(i,j) = sum_t( true(i,j,t) * filt(i,j,t) )
%                / sum_t( filt(i,j,t)^2 )
%
%   Apply as CORRECTED = K .* FILTEREDDATA. This is the least-squares
%   gain per pixel, so it restores amplitude in the least-squares sense
%   over the model's time behaviour, not pointwise at any single epoch.
%
%   The factors depend on the model's SPATIAL PATTERN - that is the
%   method's standing caveat, not an implementation detail. A model with
%   a different pattern gives different factors, and where the model
%   carries no signal the ratio is two numerical zeros. Those pixels are
%   returned as NaN rather than as a large number that would silently
%   multiply your data; INFO reports how many were dropped, so a model
%   that does not cover your region of interest is visible instead of
%   quietly producing nonsense.
%
%   Inputs
%     model      (nLat x nLon x T) double  UNFILTERED model fields, one
%                slice per epoch, in the same quantity as the data
%     latDeg     (1 x nLat | nLat x 1) double  geocentric latitudes [deg]
%     lonDeg     (1 x nLon | nLon x 1) double  longitudes [deg]
%
%   Options
%     Filter ("gauss300")  the chain the DATA saw: "none" |
%             "gauss<radiusKm>" | "DDKn" | a W struct from
%             shLowLevel.readDDK / shLowLevel.designFilter
%     Nmax (NaN)  expansion degree (NaN: floor((nLat - 1) / 2))
%     MinSignal (1e-3)  pixels whose summed filtered power is below this
%             FRACTION of the strongest pixel are returned as NaN
%     Clip ([])  (1 x 2) double  clamp k to [lo hi] ([]: no clamping).
%             Factors far from 1 usually mean the model does not
%             describe that pixel, not that the signal is huge
%     Quiet (false)  suppress the coverage summary
%
%   Outputs
%     k          (nLat x nLon) double  scaling factors; NaN where the
%                model carries too little signal to define one
%     info       (1,1) struct  fields: covered (1,1 double, fraction of
%                pixels with a finite k), kMedian / kMin / kMax (1,1
%                double, over the finite pixels), clipped (1,1 double,
%                number of clamped pixels), nmax (1,1 double),
%                filter (1,1 string), epochs (1,1 double)
%
%   Python-validated in tools/dev/validate_leakage.py: amplitude
%   invariance (k depends on the pattern, not the model's scale, to
%   1e-15), error reduction on a field whose mixture differs from the
%   model's, and NaN coverage where the model is silent.
%
%   Example
%     k = shLowLevel.gridScaling(modelEWH, -89:89, 0:359, ...
%             Filter = "gauss300");
%     corrected = k .* filteredEWH;      % NaN outside the model's reach
%
%   See also shLowLevel.leakageCorrect, shLowLevel.basinScaling.
%
%   Developed by Matthias Weigelt with the help of Claude (Opus 5),
%   2026-08-11 (v3.2.0).
arguments
    model double
    latDeg double
    lonDeg double
    opts.Filter = "gauss300"
    opts.Nmax (1,1) double = NaN
    opts.MinSignal (1,1) double {mustBeNonnegative} = 1e-3
    opts.Clip double = []
    opts.Quiet (1,1) logical = false
end
latDeg = latDeg(:).';
lonDeg = lonDeg(:).';
nLat = numel(latDeg);
nLon = numel(lonDeg);
assert(size(model, 1) == nLat && size(model, 2) == nLon, ...
    'shLowLevel:gridScaling:badSize', ...
    'model must be nLat x nLon x T (%d x %d x T), got %d x %d.', ...
    nLat, nLon, size(model, 1), size(model, 2));
if ~isempty(opts.Clip)
    assert(numel(opts.Clip) == 2 && opts.Clip(1) < opts.Clip(2), ...
        'shLowLevel:gridScaling:badClip', 'Clip must be [lo hi].');
end
T = size(model, 3);
nmax = opts.Nmax;
if ~isfinite(nmax)
    nmax = floor((nLat - 1) / 2);
end
assert(nmax >= 2, 'shLowLevel:gridScaling:badNmax', ...
    'The grid supports nmax = %d; give at least 3 latitudes.', nmax);
wFilt = resolveFilter(opts.Filter, nmax);

num = zeros(nLat, nLon);
den = zeros(nLat, nLon);
for t = 1:T
    tru = model(:, :, t);
    flt = applyChain(tru, latDeg, lonDeg, nmax, wFilt);
    num = num + tru .* flt;
    den = den + flt .* flt;
end
k = nan(nLat, nLon);
dmax = max(den(:));
if dmax > 0
    strong = den >= opts.MinSignal * dmax;
    k(strong) = num(strong) ./ den(strong);
end
nClip = 0;
if ~isempty(opts.Clip)
    lo = opts.Clip(1); hi = opts.Clip(2);
    bad = isfinite(k) & (k < lo | k > hi);
    nClip = nnz(bad);
    k(k < lo) = lo;
    k(k > hi) = hi;
end
fin = isfinite(k);
info = struct('covered', nnz(fin) / numel(k), ...
    'kMedian', median(k(fin), 'omitnan'), ...
    'kMin', min(k(fin)), 'kMax', max(k(fin)), ...
    'clipped', nClip, 'nmax', nmax, ...
    'filter', string(filterName(opts.Filter)), 'epochs', T);
if isempty(info.kMedian), info.kMedian = NaN; end
if isempty(info.kMin), info.kMin = NaN; end
if isempty(info.kMax), info.kMax = NaN; end
if ~opts.Quiet
    fprintf(['gridScaling: %.1f%% of pixels covered, k median %.3f ' ...
             '(range %.3f..%.3f), %d epochs\n'], 100 * info.covered, ...
            info.kMedian, info.kMin, info.kMax, T);
end
end

% ------------------------------------------------------------- helpers
function out = applyChain(field, latDeg, lonDeg, nmax, wFilt)
[C, S] = shLowLevel.shAnalysisGrid(field, latDeg, lonDeg, nmax);
if ~isempty(wFilt)
    if isstruct(wFilt)
        [C, S] = shLowLevel.applyDDK(C, S, wFilt);
    else
        C = C .* wFilt(:);
        S = S .* wFilt(:);
    end
end
out = shLowLevel.shSynthesis(C, S, 1, 1, latDeg, lonDeg);
end

function w = resolveFilter(F, nmax)
if isstruct(F)
    w = F;
    return
end
assert(isstring(F) || ischar(F), 'shLowLevel:gridScaling:badFilter', ...
    'Filter must be a string or a W struct.');
F = string(F);
if F == "none"
    w = [];
elseif startsWith(F, "gauss")
    r = double(extractAfter(F, "gauss"));
    assert(isfinite(r) && r > 0, 'shLowLevel:gridScaling:badFilter', ...
        'Gaussian filter must be "gauss<radiusKm>" (got "%s").', F);
    w = shLowLevel.shGaussianWeights(nmax, r);
elseif startsWith(F, "DDK")
    w = shLowLevel.readDDK(char(F), Nmax = nmax);
else
    error('shLowLevel:gridScaling:badFilter', ...
        'Filter must be "none", "gauss<km>", "DDKn", or a W struct.');
end
end

function s = filterName(F)
if isstruct(F)
    s = F.name;
else
    s = char(string(F));
end
end
