function v = pctile(x, p)
%PCTILE Percentile without the Statistics Toolbox (base MATLAB only).
%
%   V = shx.pctile(X, P) returns the P-th percentile(s) (0..100) of
%   the finite values of X(:), using the same linear interpolation
%   between order statistics as the Statistics Toolbox prctile: sample
%   k maps to percentile 100*(k-0.5)/n. P may be a vector (v2.5.1);
%   V has the same size as P. NaN/Inf are ignored; empty input yields
%   NaN(size(P)).
%
%   Introduced in v2.4.1 after an audit found prctile (Statistics
%   Toolbox) had slipped into plotSHMap/writeAnimation - the toolbox is
%   strictly base-MATLAB.
%
%   Inputs
%     x  double, any shape
%     p  (1,:) double in [0, 100], scalar or vector
%   Outputs
%     v  (size(p)) double  percentile value per requested p
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    x double
    p (1,:) double {mustBeInRange(p, 0, 100)}
end
x = sort(x(isfinite(x)));
n = numel(x);
if n == 0, v = NaN(size(p)); return, end
if n == 1, v = repmat(x, size(p)); return, end
q = 100 * ((1:n) - 0.5) / n;
v = interp1(q, x, min(max(p, q(1)), q(end)), 'linear');
end
