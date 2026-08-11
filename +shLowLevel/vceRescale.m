function s = vceRescale(N, Xres, idx, opts)
%VCERESCALE Per-month variance factors for the noise covariance shape.
%   s = shLowLevel.vceRescale(N, Xres, idx) estimates one variance factor per
%   month such that N_t = s(t) * N. The factor is the robust (median-based,
%   chi-square(1) corrected) ratio between the observed squared residuals
%   and the formal variances, evaluated over noise-dominated coefficients
%   (n >= MinDegree, default round(2/3 * idx.Lmax)). This is the single-component
%   Helmert VCE specialization; it absorbs the month-to-month noise level
%   changes (orbit geometry, beta-prime angle, instrument state) that a
%   stationary filter ignores.
%
%   Options:
%     MinDegree (round(2/3 * idx.Lmax))  lower bound of the
%             noise-dominated band
%     Rows (true(0, 1))  extra logical row mask restricting the
%         variance-component estimation to a subset of the
%         coefficient vector; this is how VCEBands estimates one
%         factor per order band. Empty: no restriction
%   Outputs
%     s          (T x 1) or (nBands x T) double   monthly (band-wise) VCE noise variance factors
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    N   % (P,P) double, or block struct from buildNoiseCov(Assemble='blocks')
    Xres double
    idx (1,1) struct
    opts.MinDegree (1,1) double = round(2/3 * idx.Lmax)
    opts.Rows logical = true(0, 1)   % extra row mask (banded VCE, v2.2)
end

hi = idx.n >= opts.MinDegree;
if ~isempty(opts.Rows)
    hi = hi & opts.Rows(:);
end
assert(any(hi), 'No coefficients above MinDegree - lower it.');
if isstruct(N), d = N.diagN; else, d = diag(N); end
T  = size(Xres, 2);
s  = zeros(T, 1);
chi2med = 0.454936423119573;                 % median of chi-square(1)
for t = 1:T
    s(t) = median(Xres(hi, t).^2 ./ d(hi)) / chi2med;
end
s = max(s, eps);
end
