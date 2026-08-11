function [k, info] = sensitivityKernel(idx, region, opts)
%SENSITIVITYKERNEL Basin kernel with an explicit leakage/noise trade-off.
%
%   [K, INFO] = shLowLevel.sensitivityKernel(IDX, REGION) builds an
%   averaging kernel for REGION that balances two things a plain
%   indicator kernel cannot: signal leaking in from OUTSIDE the region,
%   and GRACE error propagating into the average.
%
%   shLowLevel.basinKernel returns the band-limited exact indicator. It
%   recovers the region perfectly from perfect data - but real data are
%   noisy, and the indicator has energy at every degree, so it amplifies
%   the high-degree error without limit. The usual remedy is to smooth,
%   which trades noise against leakage by picking a filter radius and
%   hoping. This function makes the trade explicit, minimising
%
%       J(k) = ||k - kExact||^2_M  +  Alpha * k' N k
%
%   over the kernel, where the first term measures leakage against the
%   far-field weighting M and the second the propagated noise. The
%   solution is closed form,
%
%       k = (M + Alpha * N)^-1 * M * kExact,
%
%   so no iteration and no stopping criterion is involved. Alpha = 0
%   returns the exact indicator; large Alpha shrinks the kernel towards
%   zero. This is the construction behind the ESA CCI and GravIS gridded
%   ice products (Swenson & Wahr 2002; Groh & Horwath 2021; Doehne et
%   al. 2023).
%
%   Inputs
%     idx        (1,1) struct  from shLowLevel.shIndex
%     region     the region: a function handle f(latDeg,lonDeg) -> [0,1],
%                a K x 2 polygon [latDeg lonDeg], or a mask - anything
%                shLowLevel.basinKernel accepts
%
%   Options
%     Alpha (1)  the trade-off weight. Larger = smoother kernel, less
%             noise, more leakage. There is no universally right value:
%             sweep it and look at INFO.leakage against INFO.noise (an
%             L-curve), then take the corner. Report which you used
%     Noise ([])  (P x P) double  error covariance in IDX ordering, or a
%             (P x 1) vector of formal sigmas, or [] to use a degree-
%             dependent default that mimics GRACE error growth. A real
%             covariance is much better: it is what makes the kernel
%             TAILORED rather than merely smoothed
%     FarField ([])  (P x 1) double  per-coefficient weight for the
%             leakage term ([]: 1/(n+1), which emphasises the low
%             degrees where a compact basin's leakage lives). The
%             benefit over a plain Gaussian depends strongly on this
%             choice - 2% to 32% less leakage at matched noise across
%             plausible weightings - so it is worth thinking about
%             rather than accepting
%     BufferKm (0), TaperKm (0), R (6378136.3), OverSample (2)
%             passed to shLowLevel.basinKernel for the exact indicator
%
%   Outputs
%     k          (P,1) double  the kernel, IDX ordering, ready for
%                shSeries.basinAverage or b' * x like any other kernel
%     info       (1,1) struct  fields: leakage (1,1 double,
%                ||k - kExact||_M), noise (1,1 double, sqrt(k' N k)),
%                alpha (1,1 double), kExact (P,1 double, the indicator
%                it started from), areaFraction (1,1 double), gain
%                (1,1 double, k' kExact / (kExact' kExact) - how much of
%                the region's own signal survives; well below 1 means
%                Alpha is too large)
%
%   Validated in tools/dev/validate_senskernel.py: the trade-off is
%   monotone in Alpha (so it is a real dial), the two limits are exact,
%   an L-curve corner exists, and at matched noise the tailored kernel
%   leaks less than a Gaussian.
%
%   Example
%     idx = shLowLevel.shIndex(60);
%     for a = logspace(-2, 4, 20)                 % sweep, then choose
%         [~, s(end+1)] = shLowLevel.sensitivityKernel(idx, basin, ...
%             Alpha = a, Noise = sigmaVec);
%     end
%     loglog([s.leakage], [s.noise])              % the L-curve
%
%   See also shLowLevel.basinKernel, shLowLevel.basinDeconvolve,
%   shLowLevel.leakageCorrect.
%
%   Developed by Matthias Weigelt with the help of Claude (Opus 5),
%   2026-08-11 (v3.7.0).
arguments
    idx (1,1) struct
    region
    opts.Alpha (1,1) double {mustBeNonnegative} = 1
    opts.Noise double = []
    opts.FarField double = []
    opts.BufferKm (1,1) double = 0
    opts.TaperKm (1,1) double {mustBeNonnegative} = 0
    opts.R (1,1) double = 6378136.3
    opts.OverSample (1,1) double {mustBeInteger, mustBePositive} = 2
end
[kExact, bInfo] = shLowLevel.basinKernel(idx, region, ...
    BufferKm = opts.BufferKm, TaperKm = opts.TaperKm, R = opts.R, ...
    OverSample = opts.OverSample);
kExact = kExact(:);
P = numel(kExact);

% far-field weighting: diagonal, one weight per coefficient
if isempty(opts.FarField)
    mVec = 1 ./ (idx.n(:) + 1);
else
    mVec = opts.FarField(:);
    assert(numel(mVec) == P, 'shLowLevel:sensitivityKernel:badFarField', ...
        'FarField must have %d entries (idx.P).', P);
end
assert(all(mVec >= 0), 'shLowLevel:sensitivityKernel:badFarField', ...
    'FarField weights must be nonnegative.');

% noise: full covariance, sigma vector, or a degree-dependent default
Nfull = [];
if isempty(opts.Noise)
    % GRACE-like error growth; the SHAPE is what matters here, since
    % Alpha absorbs any overall scale
    nVec = (1 + (idx.n(:) / 8).^3).^2;
elseif isvector(opts.Noise)
    assert(numel(opts.Noise) == P, 'shLowLevel:sensitivityKernel:badNoise', ...
        'A Noise vector must have %d entries (idx.P).', P);
    nVec = opts.Noise(:).^2;
else
    assert(isequal(size(opts.Noise), [P P]), ...
        'shLowLevel:sensitivityKernel:badNoise', ...
        'A Noise matrix must be %d x %d (idx.P).', P, P);
    Nfull = (opts.Noise + opts.Noise.') / 2;      % symmetrise
    nVec = diag(Nfull);
end

% k = (M + Alpha N)^-1 M kExact. With diagonal M and N this is one
% divide per coefficient; a full covariance needs the solve.
if isempty(Nfull)
    k = (mVec .* kExact) ./ (mVec + opts.Alpha * nVec);
    noise = sqrt(sum(nVec .* k.^2));
else
    A = diag(mVec) + opts.Alpha * Nfull;
    k = A \ (mVec .* kExact);
    noise = sqrt(max(k' * Nfull * k, 0));
end
d = k - kExact;
leakage = sqrt(sum(mVec .* d.^2));
den = kExact' * kExact;
if den > 0
    gain = (k' * kExact) / den;
else
    gain = NaN;
end
info = struct('leakage', leakage, 'noise', noise, 'alpha', opts.Alpha, ...
    'kExact', kExact, 'areaFraction', bInfo.areaFraction, 'gain', gain);
end
