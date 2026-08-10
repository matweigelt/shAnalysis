function spec = shOrderRMS(C, S, opts)
%SHORDERRMS Order-domain spectral diagnostics (striping direction).
%
%   SPEC = shLowLevel.shOrderRMS(C, S) sums coefficient power per ORDER m -
%   the complementary view to shLowLevel.shDegreeRMS and the natural axis for
%   GRACE striping (correlated errors live at high orders):
%
%       ordVariance(m) = sum_n ( C_nm^2 + S_nm^2 )
%       ordRMS         = sqrt(ordVariance)
%       ordAmplitude   = R (6378136.3) * ordRMS
%       cumVariance/cumRMS/cumAmplitude = running sums over m
%
%   Conventions mirror shDegreeRMS exactly (RMS = sqrt of the summed
%   power, amplitude = R * RMS; degree/order power partitions are two
%   marginals of the same total: sum(ordVariance) == sum(degVariance),
%   tested). Formal errors via sigmaC ([])/sigmaS ([]) give errVariance/errRMS/
%   errAmplitude per order.
%
%   Inputs
%     C, S        (n1,n1) double
%   Options
%     R (6378136.3) [m], n0 (0)  zero degrees below n0 first,
%     sigmaC, sigmaS ([])        formal errors
%   Outputs
%     spec  struct: order, ordVariance, ordRMS, ordAmplitude,
%           cumVariance, cumRMS, cumAmplitude, (err*), R,
%           domain = 'order'  (consumed by plotSHSpectrum)
%
%   Claude (Fable 5), 2026-08-07 (v2.4.1).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    C double
    S double
    opts.R (1,1) double = 6378136.3
    opts.n0 (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    opts.sigmaC double = []
    opts.sigmaS double = []
end
nmax = size(C, 1) - 1;
if opts.n0 > 0
    C(1:opts.n0, :) = 0; S(1:opts.n0, :) = 0;
end
ordVariance = (sum(C.^2, 1) + sum(S.^2, 1))';
spec.order = (0:nmax)';
spec.ordVariance = ordVariance;
spec.ordRMS = sqrt(ordVariance);
spec.ordAmplitude = opts.R .* spec.ordRMS;
spec.cumVariance = cumsum(ordVariance);
spec.cumRMS = sqrt(spec.cumVariance);
spec.cumAmplitude = opts.R .* spec.cumRMS;
spec.R = opts.R;
spec.domain = 'order';
if ~isempty(opts.sigmaC) && ~isempty(opts.sigmaS)
    sC = opts.sigmaC; sS = opts.sigmaS;
    sC(isnan(sC)) = 0; sS(isnan(sS)) = 0;
    if opts.n0 > 0
        sC(1:opts.n0, :) = 0; sS(1:opts.n0, :) = 0;
    end
    ev = (sum(sC.^2, 1) + sum(sS.^2, 1))';
    spec.errVariance = ev;
    spec.errRMS = sqrt(ev);
    spec.errAmplitude = opts.R .* spec.errRMS;
end
end
