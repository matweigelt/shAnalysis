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
end

[P, T] = size(X);
p = opts.Order;
if T <= p + 1
    error('shLowLevel:estimateVAR:tooShort', ...
        'Need T > Order + 1 epochs (got T = %d, Order = %d).', T, p);
end
if T < 3 * P && opts.Shrink == 0
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

% ---- Yule-Walker: [Sigma(1)..Sigma(p)] = [Phi_1..Phi_p] * S,
%      S block-Toeplitz with S{i,j} = Sigma(i-j), Sigma(-h) = Sigma(h)'
S = zeros(p * P);
for i = 1:p
    for j = 1:p
        h = i - j;
        if h >= 0
            blk = Sig{h+1};
        else
            blk = Sig{-h+1}.';
        end
        S((i-1)*P+1:i*P, (j-1)*P+1:j*P) = blk;
    end
end
G = [Sig{2:p+1}];                       % (P x p*P)
PhiAll = (S.' \ G.').';                 % solve, no explicit inverse
Phi = cell(1, p);
for k = 1:p
    Phi{k} = PhiAll(:, (k-1)*P+1:k*P);
end

% ---- process noise Q = Sigma(0) - sum Phi_k Sigma(k)'  (Luetkepohl)
Q = Sig{1};
for k = 1:p
    Q = Q - Phi{k} * Sig{k+1}.';
end
Q = (Q + Q.') / 2;

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

model = struct('Phi', {Phi}, 'Q', Q, 'Sigma0', (Sig{1} + Sig{1}.') / 2, ...
    'order', p, 'P', P, 'specRadius', specRadius);
end
