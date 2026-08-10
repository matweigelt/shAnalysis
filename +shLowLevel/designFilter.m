function [W, info] = designFilter(SigmaC, SigmaS, opts)
%DESIGNFILTER Build a DDK-class anisotropic filter from user covariances.
%
%   W = shLowLevel.designFilter(SIGMAC, SIGMAS) constructs the Wiener-
%   type decorrelation filter W = (N + Alpha * inv(S))^-1 * N per
%   order/parity block - the same construction as the published DDK
%   filters, but matched to YOUR solution's error information instead of
%   a fixed release. N is the normal (inverse noise covariance) matrix
%   built from the formal sigmas (diagonal) or from a full covariance
%   (Noise=); S is the signal covariance from a Kaula-type degree rule
%   or your own model (Signal= / DegVar=). The result uses the exact
%   block format of shLowLevel.readDDK, so it drops into applyDDK,
%   g.filter/ts.filter("ddk", W = ...) and every downstream tool
%   unchanged. Python-validated: block eigenvalues in (0, 1) for SPD
%   inputs, identity limit for Alpha -> 0, pinned 2 x 2 reference.
%
%   Inputs
%     SigmaC  (nmax+1 x nmax+1) double  formal sigmas of the cosine
%             coefficients, C(n+1, m+1) layout (e.g. g.sigmaC)
%     SigmaS  (nmax+1 x nmax+1) double  same for the sine coefficients
%
%   Options
%     Alpha (1)        regularization weight: larger = stronger
%                      smoothing (DDK1-like), smaller = weaker (DDK7+)
%     Kaula (1e-5)     signal RMS per coefficient K / n^2 (ignored when
%                      DegVar or Signal is given)
%     DegVar ([])      (nmax+1 x 1) double  per-coefficient signal
%                      VARIANCE per degree (overrides Kaula)
%     Noise ([])       (P x P) double  full noise covariance in
%                      shLowLevel.shIndex ordering (overrides sigmas)
%     Signal ([])      (P x P) double  full signal covariance, same
%                      ordering (overrides Kaula/DegVar)
%     Idx ([])         shLowLevel.shIndex struct, required with
%                      Noise/Signal
%     MinDegree (2)    first filtered degree; below it W passes through
%     Name ("designFilter")  stored in W.name
%
%   Outputs
%     W          (1,1) struct  fields: nmax (1,1 double), name (char),
%                blocks (1,K struct: m (1,1 double), cs ('c'|'s'),
%                n (1,B double, degrees), M (B x B double)) - identical
%                layout to shLowLevel.readDDK
%     info       (1,1) struct  fields: alpha (1,1 double), mode (string:
%                "diagonal" | "fullcov"), gainRange (1,2 double,
%                min/max eigenvalue over all blocks)
%
%   Example
%     W = shLowLevel.designFilter(g.sigmaC, g.sigmaS, Alpha = 1);
%     gF = g.applyDDK(W);                % apply like any DDK
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-10 (v3.1.0).
arguments
    SigmaC double
    SigmaS double
    opts.Alpha (1,1) double {mustBePositive} = 1
    opts.Kaula (1,1) double {mustBePositive} = 1e-5
    opts.DegVar double = []
    opts.Noise double = []
    opts.Signal double = []
    opts.Idx = []
    opts.MinDegree (1,1) double {mustBeInteger, mustBeNonnegative} = 2
    opts.Name (1,1) string = "designFilter"
end
fullcov = ~isempty(opts.Noise) || ~isempty(opts.Signal);
if fullcov && isempty(opts.Idx)
    error('shLowLevel:designFilter:needIdx', ...
        'Noise=/Signal= require the matching Idx= (shLowLevel.shIndex).');
end
if ~isequal(size(SigmaC), size(SigmaS))
    error('shLowLevel:designFilter:sizeMismatch', ...
        'SigmaC and SigmaS must share the same size.');
end
nmax = size(SigmaC, 1) - 1;
n0 = max(opts.MinDegree, 0);
% per-degree signal variance
if ~isempty(opts.DegVar)
    dv = opts.DegVar(:);
    if numel(dv) ~= nmax + 1
        error('shLowLevel:designFilter:badDegVar', ...
            'DegVar must have nmax+1 entries (degrees 0..nmax).');
    end
else
    nn = (0:nmax)';
    dv = (opts.Kaula ./ max(nn, 1).^2).^2;      % per-coefficient variance
end
blocks = struct('m', {}, 'cs', {}, 'n', {}, 'M', {});
evmin = inf; evmax = -inf;
for m = 0:nmax
    for cs = 'cs'
        if cs == 's' && m == 0, continue, end
        nlo = max(n0, m);                       % degrees >= max(m, MinDegree)
        if m == 0, nlo = max(nlo, 1); end       % degree 0 never filtered
        nvec = nlo:nmax;
        if isempty(nvec), continue, end
        B = numel(nvec);
        if fullcov
            sel = pickIdx(opts.Idx, nvec, m, cs);
            Nblk = inv(opts.Noise(sel, sel));
            if isempty(opts.Signal)
                Sblk = diag(dv(nvec + 1));
            else
                Sblk = opts.Signal(sel, sel);
            end
        else
            if cs == 'c', sg = diag(SigmaC(nvec + 1, m + 1));
            else, sg = diag(SigmaS(nvec + 1, m + 1)); end
            if any(diag(sg) <= 0)
                error('shLowLevel:designFilter:badSigma', ...
                    'All sigmas of filtered degrees must be positive (m = %d).', m);
            end
            Nblk = diag(1 ./ diag(sg).^2);
            Sblk = diag(dv(nvec + 1));
        end
        M = (Nblk + opts.Alpha * inv(Sblk)) \ Nblk;             %#ok<MINV>
        ev = eig((M + M') / 2);
        evmin = min(evmin, min(ev)); evmax = max(evmax, max(ev));
        blocks(end+1) = struct('m', m, 'cs', cs, 'n', nvec, 'M', M); %#ok<AGROW>
    end
end
W = struct('nmax', nmax, 'name', char(opts.Name), 'blocks', blocks);
info = struct('alpha', opts.Alpha, ...
    'mode', string(ternary(fullcov, "fullcov", "diagonal")), ...
    'gainRange', [evmin, evmax]);
end

function sel = pickIdx(idx, nvec, m, cs)
% linear indices of (n in nvec, m, cs) in shLowLevel.shIndex ordering
sel = zeros(numel(nvec), 1);
for k = 1:numel(nvec)
    hit = find(idx.n == nvec(k) & idx.m == m & idx.cs == (cs == 's'), 1);
    assert(~isempty(hit), 'shLowLevel:designFilter:idxMiss', ...
        'Idx has no entry for n=%d, m=%d, %s.', nvec(k), m, cs);
    sel(k) = hit;
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
