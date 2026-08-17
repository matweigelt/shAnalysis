function [ts, rep] = kalmanChain(obsIn, opts)
%KALMANCHAIN Kalman-smoothed SH series in one call (Kurtenbach/Kvas).
%
%   [TS, REP] = shLowLevel.kalmanChain(OBSIN, ModelSeries=TSM) is the
%   single point of access to the Kalman/VAR module: it turns noisy
%   per-epoch gravity information plus a geophysical model series into
%   a temporally smoothed shSeries, encoding the correctly ordered
%   chain of Kurtenbach (2012) and Kvas (2019):
%     1. reduce the model series TSM by its own climatology (mean,
%        trend, annual, semi-annual) - the residual approximates the
%        short-term process (Kurtenbach eq. 3.81)
%     2. estimate the VAR(Order) process model from the residuals
%        (shLowLevel.estimateVAR; Order=1 is Kurtenbach's B/Q)
%     3. build one observation record per requested epoch: a solution
%        with noise covariance, a SINEX normal equation, or a gap
%     4. forward Kalman filter with stationary initialization
%        (shLowLevel.kalmanFilter)
%     5. RTS smoother (shLowLevel.rtsSmoother) - with step 4 exactly
%        the joint adjustment over all epochs (Kvas Sec. 2.3)
%     6. repack into an shSeries with formal sigmas; in solution mode
%        the observation climatology removed in step 3 is restored, so
%        TS is directly comparable to the input series
%
%   Observation modes (by class of OBSIN):
%     shSeries      solution mode: each epoch's coefficient vector is
%                   the observation l = x + v. The series is reduced by
%                   its own climatology first (Climatology=true); the
%                   noise covariance comes from NoiseCov or, if empty,
%                   from the series' formal sigmas (errors loudly if
%                   neither exists - no silent noise assumptions).
%     string/char   neq mode: a folder of SINEX files carrying
%                   +SOLUTION/NORMAL_EQUATION_MATRIX blocks (e.g. ITSG,
%                   shLowLevel.fetchITSGSINEX). Each file contributes
%                   N and b = N * x_est at its epoch. THE NEQS ARE
%                   TAKEN AS-IS: they refer to the provider's
%                   background models, so TS is a RESIDUAL series with
%                   respect to that background - restoring it (GAX,
%                   static field) is the caller's responsibility, as
%                   everywhere else in the toolbox.
%
%   Inputs
%     obsIn  (1 x 1) shSeries | (1 x 1) string  observations: series
%            (solution mode) or SINEX folder (neq mode)
%
%   Options
%     ModelSeries (required)  (1 x 1) shSeries  process realization at
%            the target sampling and nmax (truncate first if needed);
%            e.g. GAX/AOD1B, ESA ESM or a hydrology model as SH series
%     Order (1)         (1 x 1) VAR order p (see estimateVAR)
%     Shrink (0)        (1 x 1) diagonal loading (see estimateVAR)
%     CondFun ([])      covariance conditioning (see estimateVAR)
%     Epochs ([])       (K x 1) double  target epoch grid [decimal
%            years]. [] uses the observation epochs as-is. A grid entry
%            with no observation within Tolerance becomes a gap epoch:
%            prediction only - this is how daily interpolation across
%            missing days works
%     Tolerance (0.5/365.25)  (1 x 1) max |epoch difference| [yr] to
%            match an observation to a grid entry
%     NoiseCov ([])     solution mode: (P x P) double full R, (P x 1)
%            double diagonal, or (P x T) double per-epoch diagonals;
%            [] derives diagonals from the series' formal sigmas
%     Climatology (true)  solution mode: reduce OBSIN by its own
%            climatology before filtering and restore it afterwards
%     Smoother (true)   run the RTS pass (false: filter only, e.g. for
%            NRT-style one-way runs a la Kvas Ch. 3; then the internal
%            filter stores only diagonals and memory stays flat)
%     QC ("none")       "none" | "flag" | "reject": innovation-based
%            per-epoch quality control (Kvas Sec. 3.3; see
%            kalmanFilter). "reject" keeps blunder epochs out of the
%            recursion - they become prediction-only and are counted
%            in REP.nRejected
%     QCAlpha (1e-3)    false-alarm level of the test
%     StoreCov ("auto") where the filter keeps the full covariances the
%            smoother needs: "auto" = "full" (RAM) with the smoother,
%            "diag" without; "matfile" streams them through a v7.3 MAT
%            on disk instead - RAM stays flat for long daily runs,
%            results identical to the last bit. The file is deleted
%            after smoothing (it is chain-internal); REP.memGB then
%            reports the DISK footprint
%     Pattern ("*.snx*")  neq mode: file pattern inside the folder
%
%   Outputs
%     ts   (1 x 1) shSeries  smoothed series on the epoch grid, formal
%          1-sigma stacks filled, every step in the history
%     rep  (1 x 1) struct  report: mode, order, nEpochs, nGaps,
%          qcStat (1 x K double, NaN when QC="none"), nRejected,
%          specRadius, meanContribution (mean data share per epoch,
%          Kurtenbach Sec. 3.3.2; NaN in gaps), epochs, gap (logical),
%          model (the estimateVAR struct, for reuse), memGB (covariance
%          storage actually used)
%
%   Example  % daily-type smoothing of a noisy series with a GAX prior
%     tsm = shSeries.fromFolder(gaxFolder).truncate(40);
%     [ts, rep] = shLowLevel.kalmanChain(tsObs, ModelSeries=tsm, Order=1);
%
%   Memory: the smoother stores two full (pP x pP) covariances per
%   epoch: n40/p1/T365 ~ 16 GB. For long daily series run in windows
%   or set Smoother=false; REP.memGB reports the actual footprint.
%
%   Reference: Kurtenbach, DGK C-683 (2012); Kvas, TU Graz PhD thesis
%   (2019). Numerics pre-validated in Python
%   (tools/dev/validate_kalman.py).
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-17, 20:10 UTC.

arguments
    obsIn
    opts.ModelSeries (1,1) shSeries
    opts.Order (1,1) double {mustBeInteger, mustBePositive} = 1
    opts.Shrink (1,1) double {mustBeNonnegative} = 0
    opts.CondFun = []
    opts.Epochs (:,1) double = []
    opts.Tolerance (1,1) double {mustBePositive} = 0.5/365.25
    opts.NoiseCov double = []
    opts.Climatology (1,1) logical = true
    opts.Smoother (1,1) logical = true
    opts.QC (1,1) string {mustBeMember(opts.QC, ["none","flag","reject"])} = "none"
    opts.QCAlpha (1,1) double {mustBeInRange(opts.QCAlpha, 0, 1, "exclusive")} = 1e-3
    opts.StoreCov (1,1) string {mustBeMember(opts.StoreCov, ["auto","full","matfile"])} = "auto"
    opts.Pattern (1,1) string = "*.snx*"
end

tsm = opts.ModelSeries;
nmax = tsm.nmax;
idx = shLowLevel.shIndex(nmax);
P = idx.P;

% ---- 1./2. process model from the model series residuals
[~, resid] = tsm.climatology;
X = zeros(P, resid.nEpochs);
for k = 1:resid.nEpochs
    g = resid.at(k);
    X(:, k) = shLowLevel.vecFromCS(g.C, g.S, idx);
end
model = shLowLevel.estimateVAR(X, Order=opts.Order, ...
    Shrink=opts.Shrink, CondFun=opts.CondFun);

% ---- 3. observation records
solMode = isa(obsIn, 'shSeries');
if solMode
    tsObs = obsIn;
    if tsObs.nmax ~= nmax
        error('shLowLevel:kalmanChain:nmaxMismatch', ...
            'Observation nmax %d ~= model nmax %d - truncate first.', ...
            tsObs.nmax, nmax);
    end
    if opts.Climatology
        [climObs, tsObs] = tsObs.climatology;
    else
        climObs = [];
    end
    obsEpochs = tsObs.epochs(:);
    getObs = @(k) solutionRecord(tsObs, k, idx, opts.NoiseCov);
elseif isstring(obsIn) || ischar(obsIn)
    files = dir(fullfile(char(obsIn), char(opts.Pattern)));
    if isempty(files)
        error('shLowLevel:kalmanChain:noFiles', ...
            'No files matching %s in %s.', opts.Pattern, obsIn);
    end
    snx = cell(1, numel(files));
    obsEpochs = zeros(numel(files), 1);
    for k = 1:numel(files)
        snx{k} = shLowLevel.readSINEX( ...
            fullfile(files(k).folder, files(k).name), Index=idx);
        if ~strcmp(snx{k}.kind, 'NEQ')
            error('shLowLevel:kalmanChain:notNEQ', ...
                '%s carries no normal-equation block (kind = ''%s'').', ...
                files(k).name, snx{k}.kind);
        end
        obsEpochs(k) = snx{k}.epoch;
    end
    [obsEpochs, ord] = sort(obsEpochs);
    snx = snx(ord);
    climObs = [];
    getObs = @(k) neqRecord(snx{k});
else
    error('shLowLevel:kalmanChain:badInput', ...
        'OBSIN must be an shSeries or a SINEX folder path.');
end

if isempty(opts.Epochs)
    epGrid = obsEpochs;
else
    epGrid = sort(opts.Epochs(:));
end
T = numel(epGrid);
obs = repmat(struct('l', [], 'R', [], 'N', [], 'b', []), 1, T);
for t = 1:T
    [d, k] = min(abs(obsEpochs - epGrid(t)));
    if d <= opts.Tolerance
        obs(t) = getObs(k);
    end
end

% ---- 4./5. filter and smooth
if opts.StoreCov == "auto"
    store = "full"; if ~opts.Smoother, store = "diag"; end
else
    store = opts.StoreCov;
end
filt = shLowLevel.kalmanFilter(model, obs, StoreCov=store, ...
    QC=opts.QC, QCAlpha=opts.QCAlpha);
if opts.Smoother
    smo = shLowLevel.rtsSmoother(filt);
    xs = smo.xs; sig = smo.sig;
    if store == "matfile" && isfile(filt.covFile)
        delete(filt.covFile);                  % chain-internal temp
    end
else
    xs = filt.xf(1:P, :); sig = sqrt(max(filt.dPf(1:P, :), 0));
    if store == "matfile" && isfile(filt.covFile)
        delete(filt.covFile);                  % never read: no leak
    end
end

% ---- 6. repack (restore the observation climatology on the grid first)
L1 = nmax + 1;
Cs = zeros(L1, L1, T); Ss = zeros(L1, L1, T);
sC = zeros(L1, L1, T); sS = zeros(L1, L1, T);
for t = 1:T
    [Cs(:,:,t), Ss(:,:,t)] = shLowLevel.csFromVec(xs(:, t), idx);
    [sC(:,:,t), sS(:,:,t)] = shLowLevel.csFromVec(sig(:, t), idx);
    if solMode && opts.Climatology
        gc = climObs.eval(epGrid(t));
        Cs(:,:,t) = Cs(:,:,t) + gc.C;
        Ss(:,:,t) = Ss(:,:,t) + gc.S;
    end
end
ts = shSeries(Cs, Ss=Ss, Epochs=epGrid, SigmaCs=sC, SigmaSs=sS, ...
    ProductType="KALMAN", History=sprintf( ...
    "kalmanChain: %s mode, VAR(%d), %d epochs (%d gaps), smoother=%d", ...
    tern(solMode, "solution", "neq"), opts.Order, T, nnz(filt.gap), ...
    opts.Smoother));

memGB = 8 * numel(epGrid) * (model.order * P)^2 * 2 * double(opts.Smoother) / 1e9;
rep = struct('mode', tern(solMode, "solution", "neq"), ...
    'order', opts.Order, 'nEpochs', T, 'nGaps', nnz(filt.gap), ...
    'qcStat', filt.qcStat, 'nRejected', nnz(filt.qcReject), ...
    'specRadius', model.specRadius, ...
    'meanContribution', mean(filt.contrib, 1), ...
    'epochs', epGrid, 'gap', filt.gap, 'model', model, 'memGB', memGB);
end

% ---------------------------------------------------------------- local
function o = solutionRecord(tsObs, k, idx, noiseCov)
g = tsObs.at(k);
o = struct('l', shLowLevel.vecFromCS(g.C, g.S, idx), 'R', [], 'N', [], 'b', []);
if isempty(noiseCov)
    if isempty(g.sigmaC) || all(g.sigmaC(:) == 0)
        error('shLowLevel:kalmanChain:noNoise', ...
            ['Epoch %d has no formal sigmas and NoiseCov is empty - ' ...
             'pass NoiseCov explicitly (no silent noise assumptions).'], k);
    end
    o.R = shLowLevel.vecFromCS(g.sigmaC, g.sigmaS, idx).^2;
elseif isvector(noiseCov)
    o.R = noiseCov(:);
elseif size(noiseCov, 2) == size(noiseCov, 1)
    o.R = noiseCov;
else
    o.R = noiseCov(:, k);
end
end

function o = neqRecord(s)
o = struct('l', [], 'R', [], 'N', s.M, 'b', s.M * s.x);
end

function out = tern(c, a, b)
if c, out = a; else, out = b; end
end
