function [circ, noiseR, info] = eofSeparate(R, w, opts)
%EOFSEPARATE Separate coherent modes from noise in a residual stack.
%
%   [CIRC, NOISER, INFO] = shLowLevel.eofSeparate(R, W) splits a
%   time-centred residual stack into a coherent low-mode part (for
%   ocean applications: residual circulation - interannual dynamics
%   left after the trend + seasonal fit) and a noise remainder, using
%   an area-weighted EOF decomposition with the North et al. (1982)
%   rule: leading modes are kept while their eigenvalue is separated
%   from the next one by more than its sampling uncertainty
%   lambda * sqrt(2/T).
%
%   Physical limitation, stated up front: modes whose variance lies at
%   or below the noise floor (the Marchenko-Pastur bulk edge, roughly
%   sigma^2 * (1 + sqrt(Q/T))^2 for Q pixels and T epochs) are NOT
%   separable by ANY eigenvalue criterion - the Python pre-validation
%   demonstrated exactly this failure mode before the amplitudes were
%   dimensioned above the bulk. Absence of separated modes means "not
%   detectable here", not "no circulation exists".
%
%   Inputs
%     R  (Q x T) double  residual stack, one row per (masked) pixel -
%        e.g. the per-pixel residuals after oceanChain's trend +
%        seasonal fit
%     w  (Q x 1) double  positive area weights of the pixels
%
%   Options
%     NKeep ([])    (1 x 1) fixed number of modes; [] applies the
%                   North rule
%     MaxModes (10) (1 x 1) hard cap on kept modes
%
%   Outputs
%     circ   (Q x T) double  coherent (circulation) part
%     noiseR (Q x T) double  time-centred residual minus CIRC
%     info   (1 x 1) struct  lam (variance per mode), nKeep,
%            northSeparated (logical, per leading mode),
%            varExplained (fraction, kept modes), rmsNoise
%            (area-weighted RMS of NOISER)
%
%   Example
%     % after [out, rep] = oceanChain(...): rebuild the pixel residuals
%     % R (Q x T), weights wA, then
%     [circ, nz, inf1] = shLowLevel.eofSeparate(R, wA);
%     inf1.nKeep              % modes the North rule accepts
%     inf1.rmsNoise           % the de-circulated noise level
%
%   Error identifiers
%     shLowLevel:eofSeparate:badSize  R and w are inconsistent
%
%   Validated against a Python (numpy) reference: two planted modes
%   above the noise bulk are recovered (nKeep 2, field correlation
%   0.97) and the noise RMS is reproduced to 0.5%.
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.11.0).
arguments
    R (:,:) double
    w (:,1) double {mustBePositive}
    opts.NKeep double {mustBeScalarOrEmpty} = []
    opts.MaxModes (1,1) double {mustBePositive, mustBeInteger} = 10
end
[Q, T] = size(R);
if numel(w) ~= Q
    error('shLowLevel:eofSeparate:badSize', ...
        'R has %d rows but w has %d entries.', Q, numel(w));
end
Rc = R - mean(R, 2);
W = sqrt(w);
[U, S, V] = svd(W .* Rc, 'econ');
s = diag(S);
lam = s.^2 / T;
% North et al. (1982): separated while the eigenvalue gap exceeds the
% sampling uncertainty of the leading eigenvalue
dl = lam * sqrt(2 / T);
nsep = lam(1:end-1) - lam(2:end) > dl(1:end-1);
if isempty(opts.NKeep)
    n = 0;
    while n < numel(nsep) && nsep(n + 1)
        n = n + 1;
    end
else
    n = opts.NKeep;
end
n = min([n, opts.MaxModes, numel(s)]);
if n > 0
    circ = (U(:, 1:n) ./ W) * S(1:n, 1:n) * V(:, 1:n)';
else
    circ = zeros(Q, T);
end
noiseR = Rc - circ;
wRep = repmat(w, 1, T);
rmsNoise = sqrt(sum(wRep .* noiseR.^2, 'all') / sum(wRep, 'all'));
info = struct('lam', lam, 'nKeep', n, 'northSeparated', nsep, ...
    'varExplained', sum(lam(1:n)) / sum(lam), 'rmsNoise', rmsNoise);
end
