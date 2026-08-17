function x = chi2Quantile(p, k)
%CHI2QUANTILE Chi-square quantile without toolboxes (Wilson-Hilferty).
%
%   X = shLowLevel.chi2Quantile(P, K) approximates the chi-square
%   quantile chi2inv(P, K) using the Wilson-Hilferty cube transform
%       x = k * (1 - 2/(9k) + z * sqrt(2/(9k)))^3,
%       z = sqrt(2) * erfinv(2p - 1)
%   in base MATLAB (erfinv is core; chi2inv is Statistics Toolbox -
%   same reason shLowLevel.pctile exists instead of prctile).
%
%   Accuracy, MEASURED against scipy.stats.chi2.ppf over the QC range
%   (tools/dev/validate_kalman_qc.py, check Q1): relative error
%   <= 3.5% for k <= 10, <= 0.2% at k = 50, <= 1e-4 for k >= 500. The
%   quality-control operating point of kalmanFilter is k = P (hundreds
%   to 1677 coefficients), where the approximation is essentially
%   exact; for k < 10 treat thresholds as approximate.
%
%   Inputs
%     p  (any size) double  probability level(s) in (0, 1)
%     k  (1 x 1) double     degrees of freedom, k >= 1
%
%   Outputs
%     x  (size of p) double  approximate chi-square quantile(s)
%
%   Example
%     thr = shLowLevel.chi2Quantile(1 - 1e-3, 1677);   % QC threshold
%
%   Reference: Wilson & Hilferty (1931), PNAS 17.
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-17, 23:55 UTC.

arguments
    p double {mustBeInRange(p, 0, 1, "exclusive")}
    k (1,1) double {mustBeInteger, mustBePositive}
end

z = sqrt(2) .* erfinv(2 .* p - 1);
x = k .* (1 - 2 ./ (9 * k) + z .* sqrt(2 ./ (9 * k))).^3;
x = max(x, 0);
end
