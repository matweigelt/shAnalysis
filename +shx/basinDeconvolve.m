function [avgHat, out] = basinDeconvolve(B, op, opts)
%BASINDECONVOLVE Unbiased basin averages by exact operator deconvolution.
%   [avgHat, out] = shx.basinDeconvolve(B, op) removes filter attenuation
%   and inter-basin leakage from basin averages, exploiting that W_t is a
%   KNOWN linear operator (no hydrology-model scale factors).
%
%   B (P x K): SH vectors of the basin kernels (band-limited indicator
%   functions in idx ordering). Model: the residual field within the basin
%   set is x = B*c. Then
%       B' * xf_res = (B' * W_t * B) * c   =>   c = A_t \ (B' * xf_res)
%   which is exact in the noise-free case and unbiased under noise
%   (validated to machine precision / MC in the tests and the Python
%   prototype). The deterministic part is added back via B'*model (it was
%   never filtered).
%
%   Outputs (K x T):
%     avgHat        deconvolved basin averages  ((B'B)*c + B'model)/diag(B'B)
%     out.avgNaive  naive averages of the filtered field (attenuated+leaky)
%     out.c         span coefficients
%     out.attn      diag(A_t)./diag(B'B) attenuation factors (K x T)
%     out.condA     cond(A_t) per month
%
%   Options:
%     Ridge (0)  Tikhonov parameter for ill-conditioned basin sets:
%                c = (A'A + Ridge^2 I) \ (A' v)
%   v2.5: out.sigma now includes the OLS parameter uncertainty of the
%   restored deterministic part (leverage x residual variance; MC-
%   validated), and the constrained-operator noise covariance is exact.
%   Outputs
%     avg        (K x T) double   deconvolved basin averages
%     out        struct: c (K x T) deconvolved coefficients, sigma (K x T) 1-sigma incl. deterministic term (v2.5), avgNaive, attn, condA
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    B double
    op (1,1) struct
    opts.Ridge (1,1) double {mustBeNonnegative} = 0
end

K = size(B, 2);
T = numel(op.tYears);
BtB  = B' * B;
dBtB = diag(BtB);

% Gv = V'B and Hu = Ut*B give A_t = Gv' diag(g) Hu and the deconvolved-
% average noise covariance  cov(c) = s_t * A^-1 (Gv' diag(g.^2) Gv) A^-T
% (derived from W N W' = V diag(g^2) V', N = V V'; Monte-Carlo validated,
% 0.7% at 4e4 samples). Block layout assembles the same quantities per block.
if isfield(op, 'layout') && strcmp(op.layout, 'blocks')
    Gv = zeros(size(B)); Hu = zeros(size(B));
    for k = 1:numel(op.blocks)
        b = op.blocks(k);
        Gv(b.rows, :) = b.V'  * B(b.rows, :);
        Hu(b.rows, :) = b.Ut  * B(b.rows, :);
    end
    lamAll = op.lam;
else
    Gv = op.V'  * B;
    Hu = op.Ut  * B;
    lamAll = op.lam;
end

% diagnose non-finite operator content precisely rather than letting NaN
% propagate silently into the solves (v2.1)
chk = struct('Gv', Gv, 'Hu', Hu, 'lam', lamAll, 's', op.s, ...
    'Xfres', op.Xfres, 'model', op.model);
fn = fieldnames(chk);
for q = 1:numel(fn)
    val = chk.(fn{q});
    if ~all(isfinite(val(:)))
        error('shx:basinDeconvolve:nonFiniteOperator', ...
            'Non-finite entries in %s (%d of %d): the filter operator is corrupt.', ...
            fn{q}, nnz(~isfinite(val(:))), numel(val));
    end
end

constrained = isfield(op, 'Ac') && ~isempty(op.Ac);
c = zeros(K, T); attn = zeros(K, T); condA = zeros(T, 1);
sigmaAvg = zeros(K, T);
isBlocks = isfield(op, 'layout') && strcmp(op.layout, 'blocks');
for t = 1:T
    if isBlocks && isfield(op, 'sBlocks')
        % per-row gain and noise factor honoring per-order-band VCE (v2.2)
        g = zeros(numel(lamAll), 1); sVec = zeros(numel(lamAll), 1);
        for kB = 1:numel(op.blocks)
            rB = op.blocks(kB).rows;
            g(rB) = op.blocks(kB).lam ./ (op.blocks(kB).lam + op.sBlocks(kB, t));
            sVec(rB) = op.sBlocks(kB, t);
        end
    else
        g = lamAll ./ (lamAll + op.s(t));
        sVec = repmat(op.s(t), numel(lamAll), 1);
    end
    if constrained
        % constraint correction alters W: A and (below, v2.5) the noise
        % covariance Q are built via the exact stored operator
        A = B' * shx.opApply(op, B, t);
    else
        A = Gv' * (g .* Hu);                     % == B' W_t B
    end
    if ~all(isfinite(A(:)))
        error('shx:basinDeconvolve:nonFiniteA', ...
            'Non-finite A at month %d (s(t)=%.6g, g range [%.3g, %.3g]).', ...
            t, op.s(t), min(g), max(g));
    end
    if norm(A) == 0
        error('shx:basinDeconvolve:zeroA', ...
            ['B''W_t B vanished at month %d (total damping: s(t)=%.6g, ' ...
             'max g=%.3g) - the deconvolution is undefined there.'], ...
            t, op.s(t), max(g));
    end
    v  = B' * op.Xfres(:, t);
    if constrained
        % EXACT noise covariance for the constrained W (v2.5):
        % Q = B' W (s_t N) W' B with W' B from the stored operator and
        % N = V V' (full layout; constraints never use the block path)
        WtB = shx.opApply(op, B, t, 'transp');
        Yv  = op.V' * WtB;
        Q = op.s(t) * (Yv' * Yv);
    else
        Q = Gv' * ((sVec .* g.^2) .* Gv);        % B' W N_t W' B (banded-aware)
    end
    if opts.Ridge > 0
        % Ridge is RELATIVE to the spectral norm of A (v2.1): an absolute
        % ridge is silently absorbed by floating-point addition whenever
        % ||A'A|| >> ridge^2 (eps(1e4) = 1.8e-12), leaving the system
        % exactly singular. lamR = Ridge*||A||_2 is scale-invariant; the
        % covariance uses the same regularized inverse.
        lamR = opts.Ridge * norm(A);
        Ai = (A'*A + lamR^2 * eye(K)) \ A';
        c(:, t) = Ai * v;
        covC = Ai * Q * Ai';
    else
        c(:, t) = A \ v;
        covC = (A \ Q) / A';
    end
    covAvg = (BtB * covC * BtB') ./ (dBtB * dBtB');
    sigmaAvg(:, t) = sqrt(max(diag(covAvg), 0));
    attn(:, t) = diag(A) ./ dBtB;
    condA(t)   = cond(A);
end

detPart = B' * op.model;                       % K x T, unfiltered
avgHat  = (BtB * c + detPart) ./ dBtB;

% deterministic-part uncertainty (v2.5): the restored model carries the
% OLS parameter error of the per-coefficient fit. With independent
% coefficient fits, Var(B'model(:,t)) = h_t * (B.^2)' * resVar where
% h_t = a_t'(A'A)^-1 a_t is the design leverage carried on the operator.
% Monte-Carlo validated (3000 runs, K=2, T=48): emp/pred 1.184 without
% the term (1.26 at seasonal leverage peaks), 1.003 with it (max 1.08).
% Pointwise per-epoch sigma; the det-part errors are correlated across
% epochs through the shared fit (documented).
if isfield(op, 'detLeverage') && ~isempty(op.detLeverage)
    varDetBase = (B.^2)' * op.detResVar;             % K x 1
    sigmaDet = sqrt(max(varDetBase * op.detLeverage, 0)) ./ dBtB;
    sigmaAvg = sqrt(sigmaAvg.^2 + sigmaDet.^2);
end

out.avgNaive = (B' * (op.Xfres + op.model)) ./ dBtB;
out.c = c; out.attn = attn; out.condA = condA;
% 1-sigma uncertainty of avgHat: noise cov(c) propagated through the
% basin-average map (BtB * c)./dBtB, PLUS the deterministic-fit term
% (v2.5). With constraints the noise covariance uses the exact
% constrained operator (v2.5).
out.sigma = sigmaAvg;
end
