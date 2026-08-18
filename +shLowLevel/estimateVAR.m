function model = estimateVAR(X, opts)
%ESTIMATEVAR Empirical VAR(p) process model from an SH state series.
%
%   MODEL = shLowLevel.estimateVAR(X) estimates the stochastic process
%   model of Kurtenbach et al. (2012) from a centred state realization:
%   the vector autoregressive model
%       x_t = Phi_1 x_{t-1} + ... + Phi_p x_{t-p} + w,  w ~ N(0, Q)
%   determined by the Yule-Walker equations from the empirical
%   (cross-)covariances Sigma(h) = 1/(T-h) sum x_i x_{i-h}' of the
%   series (Kvas 2019, Sec. 2.4). For Order = 1 this is EXACTLY the
%   closed form of Kurtenbach (2012), eqs. (3.84)-(3.85):
%       B = Sigma(1) Sigma(0)^-1,  Q = Sigma(0) - B Sigma(0) B'
%   (unit-tested against it). The realization X is typically a
%   geophysical model series (AOD1B/GAX, ESA ESM, hydrology) reduced
%   by mean, trend and (semi-)annual cycle - shSeries.climatology does
%   that reduction; the state ordering is shLowLevel.shIndex.
%
%   The model drives shLowLevel.kalmanFilter / rtsSmoother: with the
%   stationary initialization used there, filtering + smoothing is
%   identical to the joint least-squares adjustment over all epochs
%   (Kvas 2019, Sec. 2.3).
%
%   Inputs
%     X  (P x T) double  centred state realization, one shIndex-ordered
%        coefficient vector per column; T epochs at the target sampling
%        (daily model output for daily solutions)
%
%   Options
%     Order (1)     (1 x 1) VAR order p >= 1. Kurtenbach: 1. Higher
%                   orders capture longer temporal correlation (Kvas)
%                   at p-fold state dimension in the filter.
%     Shrink (0)    (1 x 1) diagonal loading factor: Sigma(0) +
%                   Shrink * trace(Sigma(0))/P * I. Stabilizes the
%                   Yule-Walker solve when T is not >> P.
%     Structure ("full")  coupling structure of the estimate:
%            "full"       dense Yule-Walker solve (Kurtenbach/Kvas) -
%                         needs T >> P; the right choice for daily
%                         sampling with conditioned covariances
%            "diagonal"   independent per-coefficient AR(p): Phi_k and
%                         Q diagonal - no tuning, immune to short
%                         samples
%            "orderblock" block-diagonal per (order m, C/S) block:
%                         couples coefficients within an order, none
%                         across - the structured middle ground;
%                         requires Blocks
%            MEASURED on the live ITSG monthly series (n40, T = 257,
%            36-month holdout, one-step RMS; tools/dev job 2026-08-18):
%            climatology-only 4.24e-12; "full" with weak shrink WORSE
%            (4.61e-12, overfit at T ~ 0.13 P); "diagonal" 3.35e-12;
%            "orderblock" ~ 3.36e-12 stationary. For monthly series
%            prefer "diagonal"; "full" remains the default for the
%            daily regime the Kalman module targets.
%     Blocks ({})  cell of index vectors partitioning 1:P, required
%            for Structure="orderblock" - e.g. from shIndex:
%            one block per (m, C/S) pair
%     CondFun ([])  function handle C -> C applied to every Sigma(h)
%                   before the solve, for covariance conditioning a la
%                   Kvas Sec. 2.4: land/ocean block masking and the
%                   distance-dependent correlation taper exp(-psi/psi0)
%                   (his eq. 2.120), built by the caller in the spatial
%                   domain. [] leaves the covariances untouched.
%
%   Outputs
%     model  (1 x 1) struct  process model:
%            Phi (1 x p cell of P x P double) coefficient matrices,
%            Q (P x P double) symmetric PSD process-noise covariance,
%            Sigma0 (P x P double) empirical auto-covariance Sigma(0)
%            (the stationary state covariance and the filter's initial
%            covariance), order (1 x 1), P (1 x 1), specRadius (1 x 1)
%            spectral radius of the companion matrix (< 1 = stable; a
%            warning is raised otherwise - filtering an unstable model
%            diverges over long gaps)
%
%   Example
%     [~, resid] = tsModel.climatology;              % reduce long-term parts
%     idx = shLowLevel.shIndex(tsModel.nmax);
%     X = zeros(idx.P, resid.nEpochs);
%     for k = 1:resid.nEpochs
%         g = resid.at(k);
%         X(:, k) = shLowLevel.vecFromCS(g.C, g.S, idx);
%     end
%     model = shLowLevel.estimateVAR(X, Order=1);
%
%   Reference: Kurtenbach, DGK C-683 (2012), Sec. 3.2; Kvas, TU Graz
%   PhD thesis (2019), Secs. 2.3-2.4; Luetkepohl (2005).
%   Numerics pre-validated in Python (tools/dev/validate_kalman.py).
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-17, 20:10 UTC.

arguments
    X (:,:) double {mustBeNonempty}
    opts.Order (1,1) double {mustBeInteger, mustBePositive} = 1
    opts.Shrink (1,1) double {mustBeNonnegative} = 0
    opts.CondFun = []
    opts.Structure (1,1) string {mustBeMember(opts.Structure, ["full","diagonal","orderblock"])} = "full"
    opts.Blocks cell = {}
end

[P, T] = size(X);
p = opts.Order;
if T <= p + 1
    error('shLowLevel:estimateVAR:tooShort', ...
        'Need T > Order + 1 epochs (got T = %d, Order = %d).', T, p);
end
if T < 3 * P && opts.Shrink == 0 && opts.Structure == "full"
    warning('shLowLevel:estimateVAR:shortSeries', ...
        ['T = %d epochs for P = %d parameters: the empirical covariance ' ...
         'is poorly conditioned. Consider Shrink > 0 or CondFun.'], T, P);
end

% ---- empirical (cross-)covariances Sigma(0..p) (Kurtenbach 3.82/3.83)
Sig = cell(1, p + 1);
Sig{1} = (X * X.') / T;
for h = 1:p
    Sig{h+1} = (X(:, 1+h:end) * X(:, 1:end-h).') / (T - h);
end
if opts.Shrink > 0
    Sig{1} = Sig{1} + opts.Shrink * trace(Sig{1}) / P * eye(P);
end
if ~isempty(opts.CondFun)
    for h = 1:p + 1
        Sig{h} = opts.CondFun(Sig{h});
    end
end

% ---- Yule-Walker solve, per structure. "full" is the dense original;
%      "diagonal"/"orderblock" run the SAME solve on 1 x 1 / block
%      submatrices of the Sigma(h), which is exact when the true
%      process is (block-)decoupled and a strong regularizer when it
%      is not (measured on the live monthly series - see help).
switch opts.Structure
    case "full"
        blocks = {1:P};
    case "diagonal"
        blocks = num2cell(1:P);
    case "orderblock"
        if isempty(opts.Blocks)
            error('shLowLevel:estimateVAR:needBlocks', ...
                'Structure="orderblock" requires Blocks (e.g. one index vector per (m, C/S) pair from shIndex).');
        end
        cover = sort(cell2mat(cellfun(@(b) b(:).', opts.Blocks, ...
            'UniformOutput', false)));
        if ~isequal(cover, 1:P)
            error('shLowLevel:estimateVAR:badBlocks', ...
                'Blocks must partition 1:%d exactly.', P);
        end
        blocks = opts.Blocks;
end
Phi = cell(1, p);
for k = 1:p, Phi{k} = zeros(P); end
Q = zeros(P);
for b = 1:numel(blocks)
    ib = blocks{b}(:).';
    nb = numel(ib);
    S = zeros(p * nb);
    for i = 1:p
        for j = 1:p
            h = i - j;
            if h >= 0
                blk = Sig{h+1}(ib, ib);
            else
                blk = Sig{-h+1}(ib, ib).';
            end
            S((i-1)*nb+1:i*nb, (j-1)*nb+1:j*nb) = blk;
        end
    end
    G = zeros(nb, p * nb);
    for k = 1:p
        G(:, (k-1)*nb+1:k*nb) = Sig{k+1}(ib, ib);
    end
    PhiB = (S.' \ G.').';
    Qb = Sig{1}(ib, ib);
    for k = 1:p
        Phi{k}(ib, ib) = PhiB(:, (k-1)*nb+1:k*nb);
        Qb = Qb - PhiB(:, (k-1)*nb+1:k*nb) * Sig{k+1}(ib, ib).';
    end
    Q(ib, ib) = Qb;
end
Q = (Q + Q.') / 2;
PhiAll = [Phi{:}];

% ---- stability diagnostic on the companion matrix
B = zeros(p * P);
B(1:P, :) = PhiAll;
if p > 1
    B(P+1:end, 1:(p-1)*P) = eye((p-1)*P);
end
specRadius = max(abs(eig(B)));
if specRadius >= 1
    warning('shLowLevel:estimateVAR:unstable', ...
        ['Companion spectral radius %.4f >= 1: the estimated process is ' ...
         'not stationary. Filtering will diverge over long gaps - reduce ' ...
         'the state dimension, increase Shrink, or condition with CondFun.'], ...
        specRadius);
end

Sig0 = (Sig{1} + Sig{1}.') / 2;
if opts.Structure ~= "full"
    M = zeros(P);
    for b = 1:numel(blocks), M(blocks{b}, blocks{b}) = 1; end
    Sig0 = Sig0 .* M;    % structure-consistent: B Sig0 B' + Q ~ Sig0
end
model = struct('Phi', {Phi}, 'Q', Q, 'Sigma0', Sig0, ...
    'order', p, 'P', P, 'specRadius', specRadius);
end
