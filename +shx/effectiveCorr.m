function out = effectiveCorr(x, y)
%EFFECTIVECORR Correlation with AR(1)-corrected significance.
%
%   OUT = shx.effectiveCorr(X, Y) computes the Pearson correlation of
%   two time series together with a significance test that accounts for
%   serial correlation: monthly GRACE series are strongly autocorrelated,
%   so the effective sample size Neff = T*(1-r1*r2)/(1+r1*r2) (r1, r2 =
%   lag-1 autocorrelations) is well below T, and the naive t-test badly
%   over-states significance (Monte-Carlo: ~30% false positives at
%   phi = 0.7 vs ~6% corrected; validated in Python). The two-sided
%   p-value uses the incomplete beta function (base MATLAB betainc),
%   no Statistics Toolbox required.
%
%   Outputs
%     out        (1,1) struct  fields:
%                  .r     (1,1) double  Pearson correlation
%                  .neff  (1,1) double  effective sample size (clamped [4, T])
%                  .t     (1,1) double  t statistic on neff-2 dof
%                  .p     (1,1) double  two-sided p-value
%                  .n     (1,1) double  finite pairs used
%
%   Example
%     ec = shx.effectiveCorr(avg(1,:)', avg(2,:)');   % two basin series
%     fprintf("r = %.3f (p = %.3g, Neff = %.0f)\n", ec.r, ec.p, ec.neff)
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-08 (v2.6.0).
arguments
    x (:,1) double
    y (:,1) double
end
use = isfinite(x) & isfinite(y);
x = x(use); y = y(use);
T = numel(x);
if T < 8
    error('shx:effectiveCorr:tooShort', ...
        'Need at least 8 finite pairs (got %d).', T);
end
xc = x - mean(x); yc = y - mean(y);
r = (xc' * yc) / sqrt((xc' * xc) * (yc' * yc));
lag1 = @(v) (v(1:end-1)' * v(2:end)) / (v' * v);
r1 = lag1(xc); r2 = lag1(yc);
neff = min(max(T * (1 - r1 * r2) / (1 + r1 * r2), 4), T);
nu = neff - 2;
t = r * sqrt(nu / max(1 - r^2, eps));
p = betainc(nu / (nu + t^2), nu / 2, 0.5);
out = struct('r', r, 'neff', neff, 't', t, 'p', p, 'n', T);
end
