function [b, info] = basinKernel(idx, region, opts)
%BASINKERNEL Band-limited basin kernel(s) from a region definition.
%
%   [B, INFO] = shLowLevel.basinKernel(IDX, REGION) builds the SH coefficient
%   vector of the band-limited indicator function of REGION by exact
%   Gauss-Legendre quadrature,  b = Y' * (w .* mask),  ready for
%   shSeries.basinAverage / shLowLevel.basinDeconvolve (P x 1, IDX ordering).
%   Replaces the error-prone manual Y'*(w.*mask) construction.
%
%   REGION: function handle f(latDeg,lonDeg) -> [0,1], K x 2 polygon
%   [latDeg lonDeg], or an Ngrid x 1 mask (see shLowLevel.evalMask).
%
%   Options
%     BufferKm (0)   grow (>0) / shrink (<0) the region by great-circle
%                    distance before quadrature - counters leakage by
%                    matching the kernel to the filter's resolution
%     TaperKm  (0)   >0: multiply the kernel spectrally by Jekeli
%                    Gaussian weights of this radius (shLowLevel.shGaussianWeights)
%                    - soft edges reduce ringing/leakage of the truncated
%                    indicator (best practice: Taper ~ filter half-width)
%     R (6378136.3)  sphere radius [m] for the km conversions
%     OverSample (2) quadrature refinement factor: indicator masks are
%                    not band-limited; the default doubles the ring/
%                    meridian count to tame edge quantization (~18% cap-
%                    area error at factor 1 and Lmax 20, ~5% at factor 2)
%   Outputs
%     b     (P,1) double     kernel in IDX ordering
%     info  struct: areaFraction (quadrature area of the mask / 4pi),
%           nGridPoints, buffered/tapered flags
%
%   Note: with MinDegree > 0 in IDX the kernel misses the low degrees of
%   the indicator (incl. its mean); basinAverage's (B'B) normalization
%   accounts for that consistently.
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    idx (1,1) struct
    region
    opts.BufferKm (1,1) double = 0
    opts.TaperKm (1,1) double {mustBeNonnegative} = 0
    opts.R (1,1) double = 6378136.3
    opts.OverSample (1,1) double {mustBeInteger, mustBePositive} = 2
end
[mask, ~] = shLowLevel.evalMask(idx, region, BufferKm = opts.BufferKm, ...
    R = opts.R, OverSample = opts.OverSample);
[Y, w] = shLowLevel.synthesisMatrix(idx, ...
    NLat = opts.OverSample * (idx.Lmax + 1), ...
    NLon = opts.OverSample * (2 * idx.Lmax + 2));
b = Y' * (w .* mask);
if opts.TaperKm > 0
    Wn = shLowLevel.shGaussianWeights(idx.Lmax, opts.TaperKm);
    b = b .* Wn(idx.n + 1);
end
info = struct('areaFraction', w' * mask, 'nGridPoints', numel(mask), ...
    'buffered', opts.BufferKm ~= 0, 'tapered', opts.TaperKm > 0);
end
