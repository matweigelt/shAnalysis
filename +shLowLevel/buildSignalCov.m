function [S, info] = buildSignalCov(Xres, N, idx, opts)
%BUILDSIGNALCOV Iterative, data-driven signal covariance.
%   [S, info] = shLowLevel.buildSignalCov(Xres, N, idx) estimates the signal
%   covariance from the residual series itself - no hydrology prior.
%
%   Iteration (default 3 passes):
%     1. W = S (S + N)^-1 from the current S
%     2. filter the series: Xf = W * Xres
%     3. re-estimate per-degree signal variances from Xf
%   Initialization uses the raw degree variances (noise-inflated at high
%   degrees; the first filter pass removes that and the iteration
%   converges quickly).
%
%   Modes:
%     'isotropic'      S = diag of per-coefficient variances sigma_n^2/(2n+1)
%     'inhomogeneous'  additionally modulates the covariance by a spatial
%        signal-variance map d(theta,lambda) estimated from the filtered
%        fields:  S = G * S_iso * G'  with  G = A * diag(sqrt(d/mean(d))) * Y,
%        where Y/A are the exact synthesis/analysis operators. For a
%        constant map, G = I and the isotropic case is recovered exactly
%        (unit-tested). This yields Klees-style non-stationarity - less
%        smoothing where the data show strong signal (continents), more
%        over the quiet ocean.
%
%   Options:
%     Mode        ('isotropic') | 'inhomogeneous'
%     NIter       (3)
%     FloorRel    (1e-8)  relative variance floor
%     MapSmooth   (true)  cap the variance map dynamic range to [0.1, 10]
%
%   info: c per iteration (per-coefficient variances), variance map (mode 2).
%   Outputs
%     S          (P x P) double   signal covariance in idx ordering
%     info       struct: mode, iterations, degree-variance model
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    Xres double
    N   % (P,P) double, or block struct from buildNoiseCov(Assemble='blocks')
    idx (1,1) struct
    opts.Mode {mustBeMember(opts.Mode, {'isotropic','inhomogeneous'})} = 'isotropic'
    opts.NIter (1,1) double {mustBeInteger, mustBePositive} = 3
    opts.FloorRel (1,1) double = 1e-8
    opts.MapSmooth (1,1) logical = true
end

P = idx.P;
assert(size(Xres,1) == P, 'Xres must be P x T.');

degVarPerCoef = @(X) accumVar(X, idx);

v = degVarPerCoef(Xres);                     % init: signal + noise
v = max(v, opts.FloorRel * max(v));
blockN = isstruct(N);
if blockN
    S = [];                                  % never densified in block mode
else
    S = diag(v);
end
info.c = v;

Xf = Xres;
if blockN
    assert(strcmp(opts.Mode, 'isotropic'), ...
        'shLowLevel:buildSignalCov:blocksIsotropicOnly', ...
        'Block-form N supports Mode=''isotropic'' only.');
end
for it = 1:opts.NIter
    if blockN
        Xf = zeros(size(Xres));
        for k = 1:numel(N.blocks)
            r  = N.blocks{k};
            Sb = diag(v(r));
            Xf(r, :) = (Sb / (Sb + N.mats{k})) * Xres(r, :);
        end
    else
        W  = S / (S + N);
        Xf = W * Xres;
    end
    v  = degVarPerCoef(Xf);
    v  = max(v, opts.FloorRel * max(v));
    if ~blockN, S = diag(v); end
    info.c(:, it+1) = v;
end
if blockN
    S = spdiags(v, 0, P, P);                 % sparse diagonal: memory-safe
end

if strcmp(opts.Mode, 'inhomogeneous')
    [Y, w] = shLowLevel.synthesisMatrix(idx);
    F = Y * Xf;                              % Ngrid x T filtered fields
    d = mean(F.^2, 2);
    d = d / (w' * d);                        % weighted global mean = 1
    if opts.MapSmooth
        d = min(max(d, 0.1), 10);
        d = d / (w' * d);
    end
    G = Y' * ((w .* sqrt(d)) .* Y);          % A * diag(sqrt(d)) * Y
    S = G * S * G';
    S = (S + S') / 2;
    info.varMap = d;
end
end

function v = accumVar(X, idx)
% per-coefficient variance = degree variance / (2n+1), averaged over time
sig2n = accumarray(idx.n + 1, mean(X.^2, 2) .* 1, [idx.Lmax + 1, 1]);
% mean over time of sum over (m, cs) per degree:
% accumarray above sums mean squared coefficients per degree -> sigma_n^2
v = sig2n(idx.n + 1) ./ (2*idx.n + 1);
end
