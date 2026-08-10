function [k, info] = basinScaling(op, b, tsModel, opts)
%BASINSCALING Leakage/attenuation scaling factors by forward modelling.
%
%   [K, INFO] = shLowLevel.basinScaling(OP, B, TSMODEL) estimates the classic
%   basin scaling (gain) factor: push a MODEL series (hydrology model,
%   mascon-derived field, synthetic pattern) through the SAME filter
%   operator OP that processed the data, and regress true against
%   filtered basin averages:
%       k = sum(aTrue .* aFilt) / sum(aFilt.^2)
%   Scaled data series: k * (b' xFiltered)/(b'b). Complements
%   shLowLevel.basinDeconvolve: deconvolution inverts the operator on a kernel
%   set; scaling corrects a single series using external signal
%   structure - the standard literature approach, with its standard
%   caveat (K depends on the model's spatial pattern; INFO reports the
%   regression sigma and per-month factors so the assumption is
%   checkable).
%
%   Inputs
%     op       operator, one of
%                struct  from shLowLevel.tvANSFilter (applied via shLowLevel.opApply)
%                (P,P) double  static matrix operator (Gaussian, DDK,
%                        fan as a matrix in idx ordering)
%                handle  @(x, t) -> filtered x
%     b        (P,1)      basin kernel (shLowLevel.basinKernel, idx order)
%     tsModel  shSeries   model series matched by nearest epoch within
%                         0.05 yr against opts.tYears / op.tYears
%   Options
%     idx                 required for matrix/handle operators
%                         (shLowLevel.shIndex); structs carry their own
%     tYears              required for matrix/handle operators
%     PerMonth (false)    also return monthly factors aTrue./aFilt
%   Outputs
%     k     (1,1) double  regression scaling factor
%     info  struct: sigmaK (regression 1-sigma), kMonthly (T x 1, if
%           PerMonth), matched (T x 1 logical), aTrue, aFilt
%
%   Claude (Fable 5), 2026-08-07 (v2.4).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    op
    b (:,1) double
    tsModel (1,1) shSeries
    opts.idx struct = struct()
    opts.tYears (1,:) double = []
    opts.PerMonth (1,1) logical = false
end
if isstruct(op) && isfield(op, 'idx')
    idx = op.idx; tY = op.tYears;
    applyOp = @(x, t) shLowLevel.opApply(op, x, t);
elseif isnumeric(op)
    assert(isfield(opts.idx, 'P') && ~isempty(opts.tYears), ...
        'shLowLevel:basinScaling:needIdx', ...
        'Matrix operators need idx= and tYears=.');
    idx = opts.idx; tY = opts.tYears;
    assert(isequal(size(op), [idx.P idx.P]), 'shLowLevel:basinScaling:badOp', ...
        'Matrix operator must be P x P.');
    applyOp = @(x, t) op * x;
elseif isa(op, 'function_handle')
    assert(isfield(opts.idx, 'P') && ~isempty(opts.tYears), ...
        'shLowLevel:basinScaling:needIdx', ...
        'Handle operators need idx= and tYears=.');
    idx = opts.idx; tY = opts.tYears;
    applyOp = op;
else
    error('shLowLevel:basinScaling:badOp', ...
        'Operator must be a tvANS struct, a P x P matrix, or a handle.');
end
assert(numel(b) == idx.P, 'shLowLevel:basinScaling:badKernel', ...
    'Kernel must be P x 1 in idx ordering (P = %d).', idx.P);
T = numel(tY);
n1 = min(tsModel.nmax + 1, idx.Lmax + 1);
aTrue = nan(T, 1); aFilt = nan(T, 1);
matched = false(T, 1);
bb = b' * b;
for t = 1:T
    [dt, j] = min(abs(tsModel.epochs - tY(t)));
    if dt > 0.05, continue; end
    matched(t) = true;
    Cm = zeros(idx.Lmax + 1); Sm = Cm;
    Cm(1:n1, 1:n1) = tsModel.Cs(1:n1, 1:n1, j);
    Sm(1:n1, 1:n1) = tsModel.Ss(1:n1, 1:n1, j);
    x = shLowLevel.vecFromCS(Cm, Sm, idx);
    xf = applyOp(x, t);
    aTrue(t) = (b' * x) / bb;
    aFilt(t) = (b' * xf) / bb;
end
assert(nnz(matched) >= 3, 'shLowLevel:basinScaling:tooFewMatches', ...
    'Only %d model epochs matched tYears within 0.05 yr.', nnz(matched));
u = aFilt(matched); v = aTrue(matched);
k = (u' * v) / (u' * u);
res = v - k * u;
sigmaK = sqrt((res' * res) / max(numel(u) - 1, 1) / (u' * u));
info = struct('sigmaK', sigmaK, 'matched', matched, ...
    'aTrue', aTrue, 'aFilt', aFilt, 'nMatched', nnz(matched));
if opts.PerMonth
    info.kMonthly = aTrue ./ aFilt;
end
end
