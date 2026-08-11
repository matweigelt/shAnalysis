function [Xf, op, info] = tvANSFilter(X, tYears, idx, opts)
%TVANSFILTER Time-variable anisotropic Wiener filtering of a SH series.
%
%   [XF, OP, INFO] = shLowLevel.tvANSFilter(X, TYEARS, IDX) runs the tvANS chain
%   on the coefficient series X (P x T, IDX ordering, shLowLevel.shIndex):
%
%     1. deterministic model fit (bias/trend/annual/semi-annual); the
%        stochastic filter only ever sees the residuals
%     2. noise covariance shape N (empirical per-order/parity blocks, or
%        a released full covariance via opts.NoiseCov ([]))
%     3. per-month VCE factors s(t):  N_t = s(t) * N
%     4. iterative data-driven signal covariance S
%     5. one generalized eigendecomposition  S*U = N*U*diag(lam),
%        U'*N*U = I, giving EXACTLY
%            W_t = S (S + s_t N)^(-1) = V * diag(lam./(lam+s_t)) * U'
%        with V = inv(U'): one O(P^3) factorization for the whole series,
%        O(P^2) per month (equivalence unit-tested)
%     6. optional hard linear constraints (opts.Constraints ([]), P x q, e.g.
%        an ocean-mass kernel): Ac'*xf = Ac'*x to machine precision
%     7. deterministic model added back unfiltered
%
%   Inputs
%     X            (P,T) double    coefficient vectors per epoch
%     tYears       (T,1) double    decimal years
%     idx          struct          from shLowLevel.shIndex
%     opts.NoiseCov     (P,P) double  full covariance (else empirical)
%     opts.Constraints  (P,q) double  constraint kernels
%     opts.SignalMode   'isotropic' (default) | 'inhomogeneous'
%     opts.NIterSignal  (1,1) double  signal-covariance iterations, 3
%     opts.Robust       (1,1) logical robust deterministic fit, false
%     opts.Shrinkage    (1,1) double  empirical noise-cov shrinkage, 0.1
%     opts.VCEMinDegree (1,1) double  default round(2/3 * idx.Lmax)
%     opts.Blocks       'auto' (default) | 'on' | 'off'
%         Block-diagonal fast path: with SignalMode='isotropic', empirical
%         noise covariance and no constraints, both S and N are block-
%         diagonal per (order, C/S, degree parity), so the generalized
%         eigendecomposition decouples into O(Lmax^2) small blocks. This
%         removes the O(P^3) wall and the dense P x P memory footprint,
%         making Lmax = 96-120 tractable. 'auto' engages it exactly when
%         the conditions hold; results are IDENTICAL to the full path
%         (validated: 1.8e-15). 'on' errors if conditions are violated.
%   Outputs
%     Xf     (P,T) double   filtered series (deterministic part restored)
%     op     struct         filter operator: Ut, V, lam, s, Ac, SA, M,
%                           idx, tYears, model, Xres, Xfres
%                           (use with shLowLevel.opApply / shLowLevel.basinDeconvolve /
%                           shLowLevel.resolutionMap)
%     info   struct         noise/signal/VCE diagnostics; additionally
%                           info.sigmaXfres (P x T): 1-sigma posterior
%                           uncertainty of the filtered residual per
%                           coefficient/month, diag((I-W_t)S) =
%                           (V.^2)*(s_t*lam./(lam+s_t)) (validated to
%                           3.5e-16 against the direct matrix product).
%                           With constraints (v2.5) this is the EXACT
%                           diagonal of (W-I)S(W-I)' + s*W*N*W' - the
%                           earlier unconstrained formula underestimated
%                           along constrained directions.
%
%   Claude (Fable 5), 2026-08-07.
%
%   Options
%     VCEBands ([]) (1,:) double  order-band edges, e.g. [0 16 33 61]:
%         one monthly VCE factor per order band instead of a single
%         s(t) for the whole field. Requires the block path
%         (Blocks="on"/"auto"); the edges must be ascending and cover
%         0..Lmax. [] (default) = one factor per month (v2.2).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    X double
    tYears (:,1) double
    idx (1,1) struct
    opts.NoiseCov double = []
    opts.Constraints double = []
    opts.SignalMode {mustBeMember(opts.SignalMode, ...
        {'isotropic','inhomogeneous'})} = 'isotropic'
    opts.NIterSignal (1,1) double = 3
    opts.Robust (1,1) logical = false
    opts.Shrinkage (1,1) double = 0.1
    opts.VCEMinDegree (1,1) double = round(2/3 * idx.Lmax)
    opts.VCEBands (1,:) double = []  % order-band edges, e.g. [0 16 33 61]:
                                     % per-band monthly VCE factors
                                     % (block path only, v2.2)
    opts.Blocks {mustBeMember(opts.Blocks, {'auto','on','off'})} = 'auto'
end

P = idx.P; T = numel(tYears);
assert(isequal(size(X), [P T]), 'shLowLevel:tvANSFilter:badSize', ...
    'X must be P x T in idx ordering.');

% 1. deterministic separation
[model, Xres, ~, Afit, ~, detResVar] = shLowLevel.fitDeterministicModel( ...
    X, tYears, Robust = opts.Robust);
% leverage of the deterministic fit, h_t = a_t'(A'A)^-1 a_t: carried on
% the operator so basinDeconvolve can propagate the parameter
% uncertainty of the restored deterministic part (v2.5)
detLev = sum((Afit / (Afit' * Afit)) .* Afit, 2)';        % 1 x T

% block path eligibility (v2.2.2: an EXTERNAL NoiseCov is allowed when
% it is block-diagonal in the (order, C/S, parity) partition - verified
% by buildNoiseCov below; Blocks='auto' falls back to the full path for
% non-block-diagonal covariances, Blocks='on' errors loudly)
canBlock = strcmp(opts.SignalMode, 'isotropic') ...
    && isempty(opts.Constraints);
switch opts.Blocks
    case 'on'
        assert(canBlock, 'shLowLevel:tvANSFilter:blocksUnavailable', ...
            ['Blocks=''on'' requires SignalMode=''isotropic'' ' ...
             'and no constraints.']);
        useBlocks = true;
    case 'off'
        useBlocks = false;
    otherwise
        useBlocks = canBlock;
end

% 2. noise covariance shape
if isempty(opts.NoiseCov)
    if useBlocks
        [N, infoN] = shLowLevel.buildNoiseCov(Xres, idx, ...
            Shrinkage = opts.Shrinkage, Assemble = 'blocks');
    else
        [N, infoN] = shLowLevel.buildNoiseCov(Xres, idx, Shrinkage = opts.Shrinkage);
    end
else
    if useBlocks
        try
            [N, infoN] = shLowLevel.buildNoiseCov(Xres, idx, Mode = 'full', ...
                FullCov = opts.NoiseCov, Assemble = 'blocks');
        catch ME
            if strcmp(ME.identifier, 'shLowLevel:buildNoiseCov:notBlockDiagonal') ...
                    && strcmp(opts.Blocks, 'auto')
                useBlocks = false;               % quiet fallback for 'auto'
                [N, infoN] = shLowLevel.buildNoiseCov(Xres, idx, Mode = 'full', ...
                    FullCov = opts.NoiseCov);
            else
                rethrow(ME);
            end
        end
    else
        [N, infoN] = shLowLevel.buildNoiseCov(Xres, idx, Mode = 'full', ...
            FullCov = opts.NoiseCov);
    end
end

% 3. monthly VCE factors
s = shLowLevel.vceRescale(N, Xres, idx, MinDegree = opts.VCEMinDegree);
assert(all(isfinite(s)), 'shLowLevel:tvANSFilter:nonFiniteVCE', ...
    'Non-finite VCE factors at months: %s', mat2str(find(~isfinite(s))'));
banded = ~isempty(opts.VCEBands);
if banded
    assert(useBlocks, 'shLowLevel:tvANSFilter:bandsNeedBlocks', ...
        ['VCEBands requires the block path (Blocks="on"/"auto"): band-wise ' ...
         'noise scaling breaks the single global eigendecomposition.']);
    edges = sort(opts.VCEBands(:)');
    assert(numel(edges) >= 2 && edges(1) <= 0 && edges(end) > idx.Lmax, ...
        'shLowLevel:tvANSFilter:badBands', ...
        'VCEBands must be ascending order-edges covering 0..%d (e.g. [0 %d %d]).', ...
        idx.Lmax, ceil(idx.Lmax/2), idx.Lmax + 1);
    nBand = numel(edges) - 1;
    sBand = zeros(nBand, numel(s));
    for q = 1:nBand
        inB = idx.m >= edges(q) & idx.m < edges(q+1);
        if any(inB & (idx.n >= opts.VCEMinDegree))
            sBand(q, :) = shLowLevel.vceRescale(N, Xres, idx, ...
                MinDegree = opts.VCEMinDegree, Rows = inB);
        else
            sBand(q, :) = s;                    % band empty above MinDegree
        end
    end
    assert(all(isfinite(sBand(:))), 'shLowLevel:tvANSFilter:nonFiniteVCE', ...
        'Non-finite band VCE factors.');
end

% 4. signal covariance
[S, infoS] = shLowLevel.buildSignalCov(Xres, N, idx, ...
    Mode = opts.SignalMode, NIter = opts.NIterSignal);

% 5. generalized eigendecomposition (one-time cost)
if useBlocks
    % per-(order, C/S, parity) blocks: identical to the full path
    % (Python-validated 1.8e-15), O(sum p_b^3) instead of O(P^3)
    nb = numel(N.blocks);
    blocks = struct('rows', cell(1, nb), 'Ut', [], 'V', [], 'lam', []);
    for k = 1:nb
        r  = N.blocks{k};
        Sb = full(S(r, r));
        [Ub, Db] = eig(Sb, N.mats{k}, 'chol');
        blocks(k).rows = r;
        blocks(k).lam  = max(real(diag(Db)), 0);
        blocks(k).Ut   = Ub';
        blocks(k).V    = inv(Ub'); %#ok<MINV>
        if ~all(isfinite(blocks(k).lam)) || ~all(isfinite(Ub(:))) ...
                || ~all(isfinite(blocks(k).V(:)))
            error('shLowLevel:tvANSFilter:nonFiniteEig', ...
                ['Non-finite eigenfactors in block %d (rows %d..%d, ' ...
                 'n=%d..%d, m=%d, cs=%d): lam range [%.3g, %.3g], ' ...
                 'diag(N) range [%.3g, %.3g], diag(S) range [%.3g, %.3g].'], ...
                k, r(1), r(end), min(idx.n(r)), max(idx.n(r)), ...
                idx.m(r(1)), idx.cs(r(1)), ...
                min(blocks(k).lam), max(blocks(k).lam), ...
                min(diag(N.mats{k})), max(diag(N.mats{k})), ...
                min(diag(Sb)), max(diag(Sb)));
        end
    end
    lam = zeros(P, 1);
    for k = 1:nb, lam(blocks(k).rows) = blocks(k).lam; end
    Ac = []; SA = []; M = [];
    % per-block monthly factors: the band factor of the block's order when
    % banded, else the global factor - downstream consumers (opApply,
    % posterior sigmas, basinDeconvolve) use sBlocks uniformly (v2.2)
    nb2 = numel(blocks);
    sBlocks = repmat(s(:)', nb2, 1);
    if banded
        for k2 = 1:nb2
            mB = idx.m(blocks(k2).rows(1));
            qB = find(mB >= edges(1:end-1) & mB < edges(2:end), 1);
            sBlocks(k2, :) = sBand(qB, :);
        end
    end
    op = struct('layout', 'blocks', 'blocks', blocks, 'Ut', [], 'V', [], ...
        'lam', lam, 's', s, 'sBlocks', sBlocks, 'Ac', Ac, 'SA', SA, ...
        'M', M, 'idx', idx, ...
        'tYears', tYears, 'model', model, 'Xres', Xres, 'Xfres', [], ...
        'detLeverage', detLev, 'detResVar', detResVar(:));
else
    [U, D] = eig(S, N, 'chol');
    lam = max(real(diag(D)), 0);
    Ut  = U';
    V   = inv(Ut);                          %#ok<MINV> reused every month

    % 6. constraints (precompute)
    Ac = opts.Constraints;
    if ~isempty(Ac)
        SA = S * Ac;
        M  = Ac' * SA;
    else
        SA = []; M = [];
    end
    op = struct('layout', 'full', 'blocks', [], 'Ut', Ut, 'V', V, ...
        'lam', lam, 's', s, 'Ac', Ac, 'SA', SA, 'M', M, 'idx', idx, ...
        'tYears', tYears, 'model', model, 'Xres', Xres, 'Xfres', [], ...
        'detLeverage', detLev, 'detResVar', detResVar(:));
end

% 7. filter each month, add model back
Xfres = zeros(P, T);
for t = 1:T
    Xfres(:, t) = shLowLevel.opApply(op, Xres(:, t), t);
end
op.Xfres = Xfres;
Xf = Xfres + model;

% posterior 1-sigma of the filtered residuals: diag((I-W_t)S)
%   = (V.^2) * (s_t*lam./(lam+s_t))   (Python-validated, 3.5e-16)
sigmaXfres = zeros(P, T);
if strcmp(op.layout, 'blocks')
    for k = 1:numel(op.blocks)
        V2 = op.blocks(k).V.^2;
        lb = op.blocks(k).lam;
        for t = 1:T
            skt = op.sBlocks(k, t);              % = s(t) unless banded
            sigmaXfres(op.blocks(k).rows, t) = ...
                V2 * (skt * lb ./ (lb + skt));
        end
    end
else
    V2 = op.V.^2;
    if isempty(Ac)
        for t = 1:T
            sigmaXfres(:, t) = V2 * (s(t) * lam ./ (lam + s(t)));
        end
    else
        % EXACT posterior for the CONSTRAINED filter (v2.5). With
        %   W_t = W0 + S*Ac*M^-1*Ac'*(I - W0),  M = Ac'*S*Ac,
        % the error covariance is
        %   Cov = (W-I) S (W-I)' + s_t W N W'
        % which collapses to (I-W0)S only for the UNCONSTRAINED optimum.
        % In the eig basis (S = V diag(lam) V', N = V V', W0 = V diag(g) U',
        % g = lam/(lam+s)) the diagonal costs O(P^2 q) per month.
        % Python-validated against the brute-force covariance (2e-15);
        % the constrained-direction identity Ac'*Cov*Ac = s_t*Ac'*N*Ac
        % holds exactly. NOTE the honest correction: the previous
        % unconstrained formula UNDERESTIMATED the constrained-filter
        % error (validation ratios down to ~0.72 along constrained
        % directions), not "upper bound" as documented in v2.2-v2.4.
        Cq  = Ac' * op.V;                      % q x P
        CqT = Cq';                             % P x q
        Fq  = (lam .* CqT) / op.M;             % P x q (t-independent)
        VF  = op.V * Fq;                       % P x q (t-independent)
        for t = 1:T
            st = s(t);
            g  = lam ./ (lam + st);
            w1 = (1 - g).^2 .* lam;            % D*Lam*D diagonal
            w2 = g .* (1 - g);                 % G*E' row scale
            A2 = Cq * (w1 .* CqT);             % q x q
            A4 = Cq * ((1 - g).^2 .* CqT);     % q x q
            d  = V2 * (w1 + st * g.^2) ...
               - 2 * sum((op.V * (w1 .* CqT)) .* VF, 2) ...
               + sum((VF * A2) .* VF, 2) ...
               + 2 * st * sum((op.V * (w2 .* CqT)) .* VF, 2) ...
               + st * sum((VF * A4) .* VF, 2);
            sigmaXfres(:, t) = d;
        end
    end
end
sigmaXfres = sqrt(max(sigmaXfres, 0));
if ~isempty(Ac)
    infoConstr = ['posterior sigma exact for the constrained filter ' ...
        '(v2.5: (W-I)S(W-I)'' + s*WNW'' diagonal)'];
else
    infoConstr = '';
end

info.noise = infoN; info.signal = infoS; info.vce = s;
info.lam = lam;
info.sigmaXfres = sigmaXfres;
info.sigmaNote = infoConstr;
end
