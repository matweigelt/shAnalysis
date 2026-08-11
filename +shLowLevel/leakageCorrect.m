function [mCorr, info] = leakageCorrect(grid, latDeg, lonDeg, opts)
%LEAKAGECORRECT Iterative forward modelling of filter leakage/attenuation.
%
%   [MCORR, INFO] = shLowLevel.leakageCorrect(GRID, LATDEG, LONDEG,
%   Filter = ...) recovers the mass field that, pushed through the SAME
%   processing chain the data saw, reproduces the observed filtered
%   GRID. The chain F is analysis -> filter -> synthesis, and the
%   correction is the fixed-point iteration
%
%       m_(k+1) = m_k + Gain * (GRID - F(m_k)),      m_0 = GRID
%
%   which converges to a field whose filtered image matches the
%   observation. Filtering a GRACE solution removes real signal and
%   spreads it across basin boundaries; forward modelling puts it back
%   without needing an external model of the signal itself - only the
%   filter has to be known. That is its advantage over scaling factors
%   (shLowLevel.gridScaling), which need a model whose spatial pattern
%   you must trust.
%
%   Mask = LOGICAL confines the solution to pixels where mass is known
%   to live (a basin, a land mask, an ice sheet). This is the standard
%   and much better conditioned variant: the iteration then only has to
%   explain the observation with mass inside the region, so leakage into
%   the surroundings is removed rather than redistributed.
%
%   Inputs
%     grid       (nLat x nLon) double   observed FILTERED field, in the
%                units you want back (EWH [m], mass, ...)
%     latDeg     (1 x nLat | nLat x 1) double  geocentric latitudes [deg]
%     lonDeg     (1 x nLon | nLon x 1) double  longitudes [deg]
%
%   Options
%     Filter ("gauss300")  the chain the DATA saw: "none" |
%             "gauss<radiusKm>" | "DDKn" | a W struct from
%             shLowLevel.readDDK / shLowLevel.designFilter. Getting this
%             wrong is the one way to make the correction meaningless -
%             it must be the filter actually applied to the data
%     Mask ([])  (nLat x nLon) logical  confine the solution to these
%             pixels ([]: unconstrained, global)
%     Nmax (NaN)  expansion degree of the intermediate analysis
%             (NaN: floor((nLat - 1) / 2), the highest degree the grid
%             supports)
%     Gain (1)  relaxation factor. Values up to about 3 accelerate
%             convergence; large values DIVERGE (validated: 5 diverges
%             on the reference problem), which is detected rather than
%             returned
%     MaxIter (50)  iteration cap
%     Tol (1e-4)  stop when max|residual| / max|grid| falls below this
%     GM (3.986004415e14)  reference constants used to convert the grid
%             to coefficients and back. They CANCEL in the chain, so the
%             result does not depend on them; they exist so an unusual
%             pair can be matched to the data
%     R (6378136.3)  reference radius [m], see GM
%     Quiet (false)  suppress the per-run summary
%
%   Outputs
%     mCorr      (nLat x nLon) double   leakage-corrected field, same
%                units and grid as GRID
%     info       (1,1) struct  fields: iterations (1,1 double),
%                residual (1,1 double, final relative residual),
%                history (1,K double, the residual per iteration - plot
%                it, a rising curve means Gain is too large),
%                converged (1,1 logical), nmax (1,1 double),
%                filter (1,1 string), masked (1,1 logical)
%
%   The correction is only as good as the filter description and the
%   mask. With an exact mask on a synthetic disc it recovers the truth
%   to better than 1e-3; with no mask it converges more slowly and
%   overshoots at the edges. Python-validated in
%   tools/dev/validate_leakage.py (fixed point, zero-in/zero-out,
%   monotone residual, divergence bound).
%
%   Example
%     ewh = g.synthesis(-89:89, 0:359, quantity = "ewh", kn = kn);
%     [m, info] = shLowLevel.leakageCorrect(ewh, -89:89, 0:359, ...
%         Filter = "gauss300", Mask = basinMask);
%     fprintf("%d iterations, residual %.2e\n", info.iterations, info.residual);
%
%   See also shLowLevel.gridScaling, shLowLevel.basinDeconvolve,
%   shLowLevel.basinScaling.
%
%   Developed by Matthias Weigelt with the help of Claude (Opus 5),
%   2026-08-11 (v3.2.0).
arguments
    grid double
    latDeg double
    lonDeg double
    opts.Filter = "gauss300"
    opts.Mask = []
    opts.Nmax (1,1) double = NaN
    opts.Gain (1,1) double {mustBePositive} = 1
    opts.MaxIter (1,1) double {mustBeInteger, mustBePositive} = 50
    opts.Tol (1,1) double {mustBePositive} = 1e-4
    opts.GM (1,1) double = 3.986004415e14
    opts.R (1,1) double = 6378136.3
    opts.Quiet (1,1) logical = false
end
latDeg = latDeg(:).';
lonDeg = lonDeg(:).';
assert(isequal(size(grid), [numel(latDeg), numel(lonDeg)]), ...
    'shLowLevel:leakageCorrect:badSize', ...
    'grid must be nLat x nLon (%d x %d), got %d x %d.', ...
    numel(latDeg), numel(lonDeg), size(grid, 1), size(grid, 2));
nmax = opts.Nmax;
if ~isfinite(nmax)
    nmax = floor((numel(latDeg) - 1) / 2);
end
assert(nmax >= 2, 'shLowLevel:leakageCorrect:badNmax', ...
    'The grid supports nmax = %d; give at least 3 latitudes.', nmax);
mask = opts.Mask;
useMask = ~isempty(mask);
if useMask
    assert(isequal(size(mask), size(grid)), ...
        'shLowLevel:leakageCorrect:badMask', ...
        'Mask must have the same size as grid.');
    mask = logical(mask);
    assert(any(mask(:)), 'shLowLevel:leakageCorrect:emptyMask', ...
        'Mask selects no pixels - nothing to solve for.');
end
wFilt = resolveFilter(opts.Filter, nmax);

scale = max(abs(grid(:)));
if scale == 0                            % nothing to correct
    mCorr = grid;
    info = struct('iterations', 0, 'residual', 0, 'history', [], ...
        'converged', true, 'nmax', nmax, ...
        'filter', string(filterName(opts.Filter)), 'masked', useMask);
    return
end

m = grid;
if useMask, m(~mask) = 0; end
hist = zeros(1, opts.MaxIter);
converged = false;
for k = 1:opts.MaxIter
    r = grid - applyChain(m, latDeg, lonDeg, nmax, wFilt, opts.GM, opts.R);
    m = m + opts.Gain * r;
    if useMask, m(~mask) = 0; end
    hist(k) = max(abs(r(:))) / scale;
    % a growing residual means Gain is past the stability bound; say so
    % instead of returning a field that looks like a result
    if ~isfinite(hist(k)) || (k > 2 && hist(k) > 10 * hist(1))
        error('shLowLevel:leakageCorrect:diverged', ...
            ['The iteration diverged at step %d (relative residual ' ...
             '%.3g, started at %.3g). Gain = %g is too large for this ' ...
             'filter; 1 is safe and values above ~3 are not.'], ...
            k, hist(k), hist(1), opts.Gain);
    end
    if hist(k) < opts.Tol
        converged = true;
        break
    end
end
hist = hist(1:k);
mCorr = m;
info = struct('iterations', k, 'residual', hist(end), 'history', hist, ...
    'converged', converged, 'nmax', nmax, ...
    'filter', string(filterName(opts.Filter)), 'masked', useMask);
if ~opts.Quiet
    if converged
        fprintf(['leakageCorrect: converged in %d iterations ' ...
                 '(residual %.2e)\n'], k, hist(end));
    else
        fprintf(['leakageCorrect: STOPPED at MaxIter = %d with ' ...
                 'residual %.2e > Tol = %.1e - raise MaxIter or Gain\n'], ...
                k, hist(end), opts.Tol);
    end
end
end

% ------------------------------------------------------------- helpers
function out = applyChain(field, latDeg, lonDeg, nmax, wFilt, GM, R)
%APPLYCHAIN analysis -> filter -> synthesis, the operator the data saw.
%   GM and R must be the SAME on both sides: shAnalysisGrid converts a
%   quantity into Stokes coefficients using them, and shSynthesis
%   converts back. They therefore cancel exactly and the result does not
%   depend on their values - but mismatching them (or leaving synthesis
%   at GM = R = 1) scales the field by ~1e-7 and silently ruins the
%   iteration.
[C, S] = shLowLevel.shAnalysisGrid(field, latDeg, lonDeg, nmax, ...
    GM = GM, R = R);
if ~isempty(wFilt)
    if isstruct(wFilt)
        [C, S] = shLowLevel.applyDDK(C, S, wFilt);
    else
        C = C .* wFilt(:);
        S = S .* wFilt(:);
    end
end
out = shLowLevel.shSynthesis(C, S, GM, R, latDeg, lonDeg);
end

function w = resolveFilter(F, nmax)
%RESOLVEFILTER "none" | "gauss<km>" | "DDKn" | W struct -> weights or W.
if isstruct(F)
    w = F;
    return
end
assert(isstring(F) || ischar(F), 'shLowLevel:leakageCorrect:badFilter', ...
    'Filter must be a string or a W struct.');
F = string(F);
if F == "none"
    w = [];
elseif startsWith(F, "gauss")
    r = double(extractAfter(F, "gauss"));
    assert(isfinite(r) && r > 0, 'shLowLevel:leakageCorrect:badFilter', ...
        'Gaussian filter must be "gauss<radiusKm>" (got "%s").', F);
    w = shLowLevel.shGaussianWeights(nmax, r);
elseif startsWith(F, "DDK")
    w = shLowLevel.readDDK(char(F), Nmax = nmax);
else
    error('shLowLevel:leakageCorrect:badFilter', ...
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
