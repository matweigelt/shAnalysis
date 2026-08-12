function [N, info] = buildNoiseCov(Xres, idx, opts)
%BUILDNOISECOV Monthly-solution noise covariance (shape matrix).
%   [N, info] = shLowLevel.buildNoiseCov(Xres, idx) estimates the noise
%   covariance from the residual coefficient series Xres (P x T, residuals
%   about the deterministic model). Two modes:
%
%   'empirical' (default): per-order, per-parity empirical covariance
%       across time, assembled block-diagonal. This reproduces the GRACE
%       error structure responsible for striping (correlations between
%       same-parity degrees of the same order) without any prior model.
%       Cross-order and cross-parity terms are exactly zero. Shrinkage (0.1)
%       toward the diagonal regularizes small-sample blocks.
%       Caveat: at low degrees the residual variability is signal-
%       dominated, so N is conservative there (the Wiener weight errs
%       toward less filtering of strong signal - acceptable). For the
%       rigorous path use mode 'full'.
%
%   'full': pass a released full covariance (e.g. ITSG / COST-G) via
%       opts.FullCov ([]), already reordered to idx ordering. Returned as-is
%       (symmetrized). Combine with shLowLevel.vceRescale for monthly scaling.
%
%   Options:
%     Mode       ('empirical') | 'full'
%     FullCov    ([])   P x P matrix for mode 'full'
%     Shrinkage  (0.1)  gamma in C <- (1-g)*C + g*diag(diag(C))
%
%   info: struct with block bookkeeping.
%   Inputs
%     Xres (P x T) residual coefficient vectors (shIndex ordering)
%   Options
%     Mode ("empirical")  ("empirical" | "full") estimate from Xres or
%         wrap a user-supplied dense covariance
%     FullCov ([])     (P x P) dense covariance for Mode="full"
%     Shrinkage (0.1)  (1 x 1) shrink toward the diagonal, 0..1
%     Assemble ("full") ("full" | "blocks") container layout
%
%   Outputs
%     N          (P x P) double or block struct   noise covariance (order/parity block-diagonal for Assemble ('full')='blocks')
%     info       struct: mode, shrinkage, blocks metadata
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    Xres double
    idx (1,1) struct
    opts.Mode {mustBeMember(opts.Mode, {'empirical','full'})} = 'empirical'
    opts.FullCov double = []
    opts.Shrinkage (1,1) double {mustBeInRange(opts.Shrinkage, 0, 1)} = 0.1
    opts.Assemble {mustBeMember(opts.Assemble, {'full','blocks'})} = 'full'
end

if strcmp(opts.Mode, 'full')
    assert(~isempty(opts.FullCov) && isequal(size(opts.FullCov), [idx.P idx.P]), ...
        'opts.FullCov must be P x P in idx ordering.');
    Nf = (opts.FullCov + opts.FullCov') / 2;
    if strcmp(opts.Assemble, 'full')
        N = Nf;
        info.mode = 'full';
        return
    end
    % Assemble='blocks' (v2.2.2): slice a user-supplied covariance into
    % the (order, C/S, parity) block container - enables the fast block
    % path of tvANSFilter with EXTERNAL noise models (diagonal N, DDK-
    % style block covariances, blockwise SINEX approximations). The
    % matrix must actually BE block-diagonal in this partition; anything
    % else errors loudly instead of silently discarding correlations.
    rowsC = {}; matsC = {};
    for m = 0:idx.Lmax
        for cs = 0:double(m > 0)
            for par = 0:1
                nn = max(m, idx.minDegree):idx.Lmax;
                nn = nn(mod(nn, 2) == par);
                if isempty(nn), continue; end
                rows = squeeze(idx.pos(nn+1, m+1, cs+1));
                rows = rows(rows > 0);
                if isempty(rows), continue; end
                rowsC{end+1} = rows(:); %#ok<AGROW>
                matsC{end+1} = Nf(rows, rows); %#ok<AGROW>
            end
        end
    end
    onMass = 0;
    for k = 1:numel(rowsC)
        onMass = onMass + sum(matsC{k}(:).^2);
    end
    offMass = sum(Nf(:).^2) - onMass;
    assert(offMass <= 1e-24 * max(sum(Nf(:).^2), realmin), ...
        'shLowLevel:buildNoiseCov:notBlockDiagonal', ...
        ['The supplied NoiseCov is not block-diagonal in the ' ...
         '(order, C/S, parity) partition (off-block share %.2e of the ' ...
         'Frobenius mass). Use Blocks="off" for full covariances.'], ...
        offMass / max(sum(Nf(:).^2), realmin));
    dAll = zeros(idx.P, 1);
    for k = 1:numel(rowsC), dAll(rowsC{k}) = diag(matsC{k}); end
    N = struct('layout', 'blocks', 'blocks', {rowsC}, 'mats', {matsC}, ...
        'diagN', dAll, 'P', idx.P);
    info.mode = 'full-blocks';
    info.nBlocks = numel(rowsC);
    return
end

T = size(Xres, 2);
assert(size(Xres,1) == idx.P, 'Xres must be P x T.');

rowsC = {}; matsC = {};
for m = 0:idx.Lmax
    for cs = 0:double(m > 0)
        for par = 0:1
            nn = max(m, idx.minDegree):idx.Lmax;
            nn = nn(mod(nn, 2) == par);
            if isempty(nn), continue; end
            rows = squeeze(idx.pos(nn+1, m+1, cs+1));
            rows = rows(rows > 0);
            if isempty(rows), continue; end
            C = cov(Xres(rows, :)');                    % time as samples
            if isscalar(C), C = var(Xres(rows, :)); end
            C = (1 - opts.Shrinkage)*C + opts.Shrinkage*diag(diag(C));
            C = (C + C') / 2;
            rowsC{end+1} = rows(:); %#ok<AGROW>
            matsC{end+1} = C;       %#ok<AGROW>
        end
    end
end
nBlocks = numel(rowsC);
dAll = zeros(idx.P, 1);
for k = 1:nBlocks, dAll(rowsC{k}) = diag(matsC{k}); end
% PD guard (v2.1 fix): the floor must reference the smallest POSITIVE
% variance. Coefficient rows with EXACTLY constant residuals (variance 0)
% are legitimate - e.g. identically-zero sectoral S in synthetic data, or
% coefficients fixed by convention - and min(dAll) = 0 would collapse the
% old 1e-10*min guard to zero, leaving N singular and eig(S, N, 'chol')
% returning Inf eigenvalues (diagnosed via shLowLevel:tvANSFilter:nonFiniteEig
% on real MATLAB runs; reproduced and fix validated in the Python port).
% With the floor, zero-variance rows get lam = S/floor >> 1, i.e. gain
% ~ 1: they pass through the filter unchanged - the correct
% no-noise-information behavior for rows whose data is exactly constant.
dpos = dAll(dAll > 0);
if isempty(dpos)
    error('shLowLevel:buildNoiseCov:allZeroVariance', ...
        'Every coefficient residual series is exactly constant - no noise to estimate.');
end
floorVal = 1e-10 * min(dpos);
for k = 1:nBlocks
    matsC{k} = matsC{k} + floorVal * eye(numel(rowsC{k}));
end
dAll = dAll + floorVal;

if strcmp(opts.Assemble, 'blocks')
    % block form: struct instead of a dense P x P matrix (memory-safe for
    % high Lmax; consumed by tvANSFilter/vceRescale/buildSignalCov)
    N = struct('layout', 'blocks', 'blocks', {rowsC}, 'mats', {matsC}, ...
        'diagN', dAll, 'P', idx.P);
else
    N = zeros(idx.P);
    for k = 1:nBlocks
        N(rowsC{k}, rowsC{k}) = matsC{k};
    end
    N = (N + N') / 2;
end

info.mode = 'empirical';
info.nBlocks = nBlocks;
info.T = T;
info.shrinkage = opts.Shrinkage;
end
