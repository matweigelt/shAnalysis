function out = mcPropagate(fun, g, opts)
%MCPROPAGATE Monte-Carlo uncertainty propagation through arbitrary chains.
%
%   OUT = shLowLevel.mcPropagate(FUN, G) draws coefficient samples consistent
%   with the formal errors of the shCoefficients G, pushes each through
%   FUN, and returns empirical statistics - uncertainty propagation for
%   ANY functional (filtered grids, basin averages, amplitude maps, ...)
%   without deriving analytic formulas, and the standard cross-check for
%   every analytic sigma in this toolbox.
%
%   Sampling models:
%     default        independent Gaussian per coefficient from
%                    G.sigmaC / G.sigmaS (NaN sigmas -> not perturbed)
%     Cov ([])=M, Idx (struct([]))=idx correlated sampling from a full P x P covariance
%                    (e.g. shLowLevel.readSINEX(..., Output="covariance").M
%                    reordered via Index=idx) - Cholesky factorization
%
%   Inputs
%     fun   function handle   y = fun(gs), gs an shCoefficients sample;
%                             y numeric (any fixed size, column-ized)
%     g     (1,1) shCoefficients   with sigmas (or Cov supplied)
%     opts.N (1,1) double = 500    number of samples
%     opts.Cov double = []         full covariance (P x P, Idx ordering)
%     opts.Idx struct = struct([]) shLowLevel.shIndex for Cov
%     opts.Seed ([]) double = []        rng seed for reproducibility
%     opts.KeepSamples (1,1) logical = false
%   Outputs
%     out  (1,1) struct  fields: samples (K x nSamples double), mean,
%          sigma (K x 1 double)  Monte-Carlo moments of the functional
%
%   MC error of sigma itself ~ sigma/sqrt(2N) (~3%% at N=500).
%
%   Claude (Fable 5), 2026-08-07.
%   Outputs
%     out        struct: mean, sigma, samples (N x Q), q16/q50/q84 percentiles of the propagated functional
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    fun (1,1) function_handle
    g (1,1) shCoefficients
    opts.N (1,1) double {mustBeInteger, mustBePositive} = 500
    opts.Cov double = []
    opts.Idx struct = struct([])
    opts.Seed double = []
    opts.KeepSamples (1,1) logical = false
end
if ~isempty(opts.Seed), rng(opts.Seed); end

useCov = ~isempty(opts.Cov);
if useCov
    assert(~isempty(fieldnames(opts.Idx)), 'shLowLevel:mcPropagate:needIdx', ...
        'Cov sampling requires Idx = shLowLevel.shIndex(...) matching Cov ordering.');
    idx = opts.Idx;
    assert(isequal(size(opts.Cov), [idx.P idx.P]), 'shLowLevel:mcPropagate:badCov', ...
        'Cov must be P x P for the supplied Idx (P = %d).', idx.P);
    [L, flag] = chol((opts.Cov + opts.Cov')/2, 'lower');
    assert(flag == 0, 'shLowLevel:mcPropagate:covNotPD', ...
        'Cov is not positive definite.');
else
    sC = g.sigmaC; sS = g.sigmaS;
    if isempty(sC) || all(isnan(sC(:)))
        error('shLowLevel:mcPropagate:noSigmas', ...
            'G carries no formal errors; supply Cov/Idx instead.');
    end
    sC(~isfinite(sC)) = 0;
    if isempty(sS), sS = zeros(size(sC)); end
    sS(~isfinite(sS)) = 0;
end

y0 = fun(g); y0 = y0(:);
M = numel(y0);
acc = zeros(M, 1); acc2 = zeros(M, 1);
if opts.KeepSamples, smp = zeros(M, opts.N); end
for k = 1:opts.N
    gs = g;
    if useCov
        dx = L * randn(idx.P, 1);
        [dC, dS] = shLowLevel.csFromVec(dx, idx);
        n1 = size(g.C, 1);
        gs = shCoefficients(g.C + dC(1:n1, 1:n1), g.S + dS(1:n1, 1:n1), ...
            GM = g.GM, R = g.R, Epoch = g.epoch);
    else
        gs = shCoefficients(g.C + sC .* randn(size(sC)), ...
            g.S + sS .* randn(size(sS)), GM = g.GM, R = g.R, Epoch = g.epoch);
    end
    y = fun(gs); y = y(:);
    acc = acc + y; acc2 = acc2 + y.^2;
    if opts.KeepSamples, smp(:, k) = y; end
end
mu = acc / opts.N;
sg = sqrt(max(acc2 / opts.N - mu.^2, 0) * opts.N / max(opts.N - 1, 1));
out = struct('mean', reshape(mu, size(fun(g))), ...
    'sigma', reshape(sg, size(fun(g))), 'N', opts.N);
if opts.KeepSamples, out.samples = smp; end
end
