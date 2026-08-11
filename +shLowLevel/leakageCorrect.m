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
%   THE MASK MUST COVER EVERY REGION THAT CAN HOLD MASS, NOT ONLY THE
%   ONE YOU ARE MEASURING. This is the mistake that costs most: a mask
%   drawn around the target alone tells the iteration that mass exists
%   nowhere else, so signal from neighbouring sources is forced into the
%   target and the result is biased HIGH. Include the neighbours in the
%   mask and integrate the corrected field over the target afterwards.
%   Measured against the published GravIS Greenland series (COST-G RL02,
%   -231.1 Gt/yr over the matching span): a Greenland-only mask gives
%   -259.8 Gt/yr (+12%), while a mask that also admits the Canadian
%   Arctic, Iceland and Svalbard gives -242.5 Gt/yr (+5%) - with the
%   neighbours absorbing -112 Gt/yr that the first mask had nowhere to
%   put. The union mask also converges faster.
%
%   Inputs
%     grid       (nLat x nLon) double   observed FILTERED field, in the
%                units you want back (EWH [m], mass, ...)
%     latDeg     (1 x nLat | nLat x 1) double  geocentric latitudes [deg]
%     lonDeg     (1 x nLon | nLon x 1) double  longitudes [deg]
%
%   STOPPING IS REGULARISATION. This is an ill-posed inverse problem and
%   the iteration SEMICONVERGES: the error against the truth falls,
%   reaches a minimum, and then RISES again as the iteration begins to
%   fit noise. Meanwhile the residual keeps shrinking, so "iterate until
%   nothing changes" is precisely the wrong instruction - on the
%   reference problem in tools/dev/validate_stopping.py the final
%   solution is 361x worse than the best one, and a step-size tolerance
%   of 1e-3 stops 89x past the optimum (1e-4 never triggers at all).
%
%   Give NOISELEVEL and the run stops on the discrepancy principle
%   (Morozov): as soon as the residual reaches the noise level of the
%   data, because fitting the data more closely than its own noise is
%   fitting noise. On the reference problem that lands within 2% of the
%   best attainable solution. Without NoiseLevel the iteration is
%   UNREGULARISED, the answer depends on MaxIter, and a warning says so.
%
%   Options
%     NoiseLevel (0)  RMS noise of GRID, in GRID's units. The standard
%             estimate is the RMS of the field over the open ocean far
%             from any signal (> 1000 km from the coast), which is what
%             GRACE processing centres use as a noise metric. 0 disables
%             the discrepancy principle and falls back to Tol
%     ResidualRegion ([])  where the discrepancy principle measures the
%             residual. [] uses MASK when one is given, and the whole
%             grid otherwise. This matters: with a mask the model
%             describes ONE region while the data contains every other
%             mass signal on the globe, so a global residual never falls
%             to the noise level and the principle never fires. On the
%             GravIS Greenland case the residual was 0.0085 globally
%             against 0.0008 near the mask, with a noise level of 0.0024
%     Tau (1.2)  safety factor in the discrepancy principle; stop when
%             rms(residual) <= Tau * NoiseLevel. Values of 1.2 to 1.5
%             land nearest the optimum, 1.0 slightly overfits
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
%     MaxIter (200)  iteration cap. Small regions under a strong
%             filter converge slowly - a 6-degree disc under a 500 km
%             Gaussian needs about 260 iterations at Gain = 1 and about
%             170 at Gain = 2. The FIELD is already accurate long before
%             the step criterion is met (0.997 of the truth after 80 of
%             those 260), so a run that stops at MaxIter is usually
%             usable; check info.step to see how close it got
%     Tol (1e-4)  fallback stopping rule, used ONLY when NoiseLevel is 0.
%             Stop when the relative CHANGE OF THE SOLUTION between
%             two iterations falls below this. Not the residual: with a
%             mask the problem is inconsistent (no field confined to the
%             region reproduces an observation that has energy outside
%             it), so the residual plateaus at a nonzero floor while the
%             solution is converged
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
%                residual (1,1 double, final relative residual - with a
%                mask this floors above Tol and that is correct),
%                history (1,K double, the residual per iteration - plot
%                it, a rising curve means Gain is too large),
%                step (1,K double, the relative change of the solution
%                per iteration; this is what Tol tests),
%                stoppedBy (1,1 string: "discrepancy" | "step" |
%                "maxIter" | "zeroField" - always check this, it says
%                whether the result was regularised),
%                residualRMS (1,1 double, over ResidualRegion),
%                residualRMSGlobal (1,1 double, over the whole grid - a
%                large gap between the two means unmodelled signal
%                elsewhere) and noiseLevel (1,1 double, as given),
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
    opts.NoiseLevel (1,1) double {mustBeNonnegative} = 0
    opts.ResidualRegion = []
    opts.Tau (1,1) double {mustBePositive} = 1.2
    opts.MaxIter (1,1) double {mustBeInteger, mustBePositive} = 200
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
% Where the discrepancy principle looks. The residual must be measured
% WHERE THE MODEL IS RESPONSIBLE. With a mask the model describes one
% region while the data contains every other mass signal on the globe -
% Antarctica, Alaska, land hydrology - so a GLOBAL residual can never
% fall to the noise level however long the iteration runs, and the
% principle silently never fires. Measured on the GravIS Greenland case:
% residual 0.0085 globally against 0.0008 near the mask and a noise
% level of 0.0024.
resReg = opts.ResidualRegion;
if isempty(resReg)
    resReg = mask;                     % [] when unmasked -> global
else
    assert(isequal(size(resReg), size(grid)), ...
        'shLowLevel:leakageCorrect:badResidualRegion', ...
        'ResidualRegion must have the same size as grid.');
    resReg = logical(resReg);
    assert(any(resReg(:)), 'shLowLevel:leakageCorrect:emptyResidualRegion', ...
        'ResidualRegion selects no pixels.');
end
if isempty(resReg)
    resSel = true(size(grid));
else
    resSel = resReg;
end

scale = max(abs(grid(:)));
if scale == 0                            % nothing to correct
    mCorr = grid;
    info = struct('iterations', 0, 'residual', 0, 'history', [], ...
        'step', [], 'converged', true, 'stoppedBy', "zeroField", ...
        'residualRMS', 0, 'residualRMSGlobal', 0, ...
        'noiseLevel', opts.NoiseLevel, 'nmax', nmax, ...
        'filter', string(filterName(opts.Filter)), 'masked', useMask);
    return
end

m = grid;
if useMask, m(~mask) = 0; end
stoppedBy = "maxIter";
hist = zeros(1, opts.MaxIter);
step = zeros(1, opts.MaxIter);
converged = false;
for k = 1:opts.MaxIter
    r = grid - applyChain(m, latDeg, lonDeg, nmax, wFilt, opts.GM, opts.R);
    prev = m;
    m = m + opts.Gain * r;
    if useMask, m(~mask) = 0; end
    hist(k) = max(abs(r(:))) / scale;
    rmsR = sqrt(mean(r(resSel).^2));
    % DISCREPANCY PRINCIPLE (Morozov): once the residual has reached the
    % noise level of the data, further iterations fit noise. This is the
    % regularisation; see the note on semiconvergence in the help.
    if opts.NoiseLevel > 0 && rmsR <= opts.Tau * opts.NoiseLevel
        stoppedBy = "discrepancy";
        converged = true;
        break
    end
    % Convergence is judged on the CHANGE IN THE SOLUTION, not on the
    % residual. With a mask the problem is inconsistent - no field
    % confined to the region reproduces an observation that has energy
    % outside it - so the residual plateaus at a nonzero floor while the
    % solution is perfectly converged. Stopping on the residual would
    % report failure on exactly the case the mask exists for.
    step(k) = max(abs(m(:) - prev(:))) / max(max(abs(m(:))), eps);
    % a growing residual means Gain is past the stability bound; say so
    % instead of returning a field that looks like a result
    if ~isfinite(hist(k)) || (k > 2 && hist(k) > 10 * hist(1))
        error('shLowLevel:leakageCorrect:diverged', ...
            ['The iteration diverged at step %d (relative residual ' ...
             '%.3g, started at %.3g). Gain = %g is too large for this ' ...
             'filter; 1 is safe and values above ~3 are not.'], ...
            k, hist(k), hist(1), opts.Gain);
    end
    if opts.NoiseLevel == 0 && step(k) < opts.Tol
        stoppedBy = "step";
        converged = true;
        break
    end
end
hist = hist(1:k);
step = step(1:k);
mCorr = m;
info = struct('iterations', k, 'residual', hist(end), 'history', hist, ...
    'step', step, 'converged', converged, 'stoppedBy', stoppedBy, ...
    'residualRMS', rmsOver(grid - applyChain(m, latDeg, lonDeg, ...
        nmax, wFilt, opts.GM, opts.R), resSel), ...
    'residualRMSGlobal', rmsOver(grid - applyChain(m, latDeg, lonDeg, ...
        nmax, wFilt, opts.GM, opts.R), true(size(grid))), ...
    'noiseLevel', opts.NoiseLevel, 'nmax', nmax, ...
    'filter', string(filterName(opts.Filter)), 'masked', useMask);
if opts.NoiseLevel == 0
    warning('shLowLevel:leakageCorrect:unregularised', ...
        ['No NoiseLevel given, so the iteration is UNREGULARISED and ' ...
         'the result depends on where it stops (%d iterations here). ' ...
         'This is an ill-posed problem: the residual keeps falling ' ...
         'long after the solution has begun to degrade. Pass ' ...
         'NoiseLevel= (e.g. the open-ocean RMS of your field) to stop ' ...
         'on the discrepancy principle instead.'], k);
end
if ~opts.Quiet
    if converged
        fprintf(['leakageCorrect: converged in %d iterations ' ...
                 '(step %.2e, residual %.2e)\n'], k, step(end), hist(end));
    else
        fprintf(['leakageCorrect: STOPPED at MaxIter = %d with step ' ...
                 '%.2e > Tol = %.1e - raise MaxIter or Gain\n'], ...
                k, step(end), opts.Tol);
    end
end
end

% ------------------------------------------------------------- helpers
function v = rmsOver(F, sel)
v = sqrt(mean(F(sel).^2));
end

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
