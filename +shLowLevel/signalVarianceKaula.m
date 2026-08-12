function [ab, info] = signalVarianceKaula(Cs, Ss, epochs, opts)
%SIGNALVARIANCEKAULA Cyclostationary Kaula signal variance (VDK input).
%
%   [AB, INFO] = shLowLevel.signalVarianceKaula(CS, SS, EPOCHS)
%   estimates the degree-dependent signal-variance approximation
%   sigma_M(l) = a * l^b (Horvath et al. 2018, Eq. 2) per CALENDAR
%   month from a coefficient series - the cyclostationary signal model
%   of the VDK/VADER decorrelation filter. Feed a PRE-FILTERED series
%   (e.g. DDK4, as in the paper): the estimator reads signal, so the
%   input must not be stripe-dominated.
%
%   Method, Python-prevalidated: per epoch, the degree variance
%   sigma_l^2 = sum_m (C_lm^2 + S_lm^2) is formed, the log of its
%   square root is corrected for the exact log-chi-square bias
%   c_l = 0.5 * (psi(k/2) - log(k/2)) with k = 2l + 1 (without it the
%   intercept a biases 6% low), and a, b come from a log-log
%   least-squares fit over LRange. The per-calendar-month result is
%   the median over all years (a bias-free to 3%, single-month scatter
%   up to ~10% - absorbed by the filter-strength alpha; VDK trend
%   results are robust against alpha, Horvath et al. 2018 Fig. 7).
%
%   Inputs
%     Cs, Ss  (n+1 x n+1 x T) coefficient stack, C(n+1, m+1) indexing
%     epochs  (T x 1) decimal years (calendar month via the fraction)
%
%   Options
%     LRange ([10, 60])  (1 x 2) degree range of the log-log fit; the
%                        lower bound excludes the low degrees where
%                        real signal is not Kaula-like and the
%                        chi-square dof are small
%
%   Outputs
%     ab   (12 x 2) double  [a, b] per calendar month (Jan..Dec);
%          evaluate as sigmaM = ab(mo, 1) * l.^ab(mo, 2)
%     info (1 x 1) struct  nPerMonth (12 x 1 epochs used), LRange,
%          abAll (T x 2 per-epoch fits, for inspection)
%
%   Example
%     [ab, ~] = shLowLevel.signalVarianceKaula(ts.Cs, ts.Ss, ts.epochs);
%     sigM = ab(4, 1) * (2:96)'.^ab(4, 2);   % April signal model
%
%   Error identifiers
%     shLowLevel:signalVarianceKaula:badRange  LRange exceeds the stack
%
%   Reference: Horvath, Murboeck, Pail, Horwath (2018), Geosciences 8,
%   323, doi:10.3390/geosciences8090323.
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.13.0).
arguments
    Cs (:,:,:) double
    Ss (:,:,:) double
    epochs (:,1) double
    opts.LRange (1,2) double = [10, 60]
end
nmax = size(Cs, 1) - 1;
if opts.LRange(2) > nmax || opts.LRange(1) < 1
    error('shLowLevel:signalVarianceKaula:badRange', ...
        'LRange [%g, %g] outside the stack degrees 1..%d.', ...
        opts.LRange, nmax);
end
T = numel(epochs);
ls = (opts.LRange(1):opts.LRange(2))';
k = 2 * ls + 1;
cl = 0.5 * (psi(k / 2) - log(k / 2));      % exact log-chi2 bias
A = [ones(numel(ls), 1), log(ls)];
abAll = nan(T, 2);
for t = 1:T
    dv = zeros(numel(ls), 1);
    for i = 1:numel(ls)
        l = ls(i);
        dv(i) = sum(Cs(l+1, 1:l+1, t).^2) + sum(Ss(l+1, 2:l+1, t).^2);
    end
    y = log(sqrt(dv)) - cl;
    co = A \ y;
    abAll(t, :) = [exp(co(1)), co(2)];
end
mo = max(1, min(12, 1 + floor(mod(epochs, 1) * 12)));
ab = nan(12, 2); nPerMonth = zeros(12, 1);
for m = 1:12
    sel = mo == m;
    nPerMonth(m) = nnz(sel);
    if any(sel)
        ab(m, :) = median(abAll(sel, :), 1);
    end
end
% calendar months absent from the series inherit the overall median -
% reported, not silent
if any(isnan(ab(:, 1)))
    ab(isnan(ab(:, 1)), 1) = median(abAll(:, 1));
    ab(isnan(ab(:, 2)), 2) = median(abAll(:, 2));
end
info = struct('nPerMonth', nPerMonth, 'LRange', opts.LRange, ...
    'abAll', abAll);
end
