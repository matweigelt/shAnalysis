function spec = diffSpectrum(C1, S1, C2, S2)
%DIFFSPECTRUM Spectral comparison of two coefficient sets.
%
%   SPEC = shx.diffSpectrum(C1, S1, C2, S2) compares two solutions in
%   the spectral domain: the difference degree amplitude localizes the
%   disagreement in wavelength, the degree correlation separates "same
%   pattern, different amplitude" from genuine disagreement, and the
%   crossover of the difference against the signal spectrum gives the
%   effective resolution of agreement as a single scalar. Inputs must
%   share the same size; truncate beforehand for mixed degrees (the
%   shx.compareSolutions wrapper does this automatically).
%
%   Outputs
%     spec       (1,1) struct  fields:
%                  .n         (nmax+1 x 1) double  degree axis 0..nmax
%                  .diffAmp   (nmax+1 x 1) double  degree amplitude of (1)-(2)
%                  .amp1      (nmax+1 x 1) double  signal degree amplitude of (1)
%                  .amp2      (nmax+1 x 1) double  signal degree amplitude of (2)
%                  .degCorr   (nmax+1 x 1) double  per-degree correlation (NaN
%                                                  where a degree is all-zero)
%                  .ncross    (1,1) double  first degree where diffAmp exceeds
%                                           amp1 (NaN if never; agreement to nmax)
%
%   Example
%     gF = g.gaussian(350);
%     spec = shx.diffSpectrum(g.C, g.S, gF.C, gF.S);
%     semilogy(spec.n, [spec.amp1, spec.diffAmp]); grid on
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-08 (v2.6.0).
arguments
    C1 double
    S1 double
    C2 double
    S2 double
end
if ~isequal(size(C1), size(S1), size(C2), size(S2))
    error('shx:diffSpectrum:sizeMismatch', ...
        'All four coefficient matrices must share the same size.');
end
n1 = size(C1, 1);
n = (0:n1-1)';
diffAmp = zeros(n1, 1); amp1 = zeros(n1, 1); amp2 = zeros(n1, 1);
degCorr = nan(n1, 1);
for k = 1:n1
    c1 = C1(k, 1:k); s1 = S1(k, 1:k);
    c2 = C2(k, 1:k); s2 = S2(k, 1:k);
    amp1(k) = sqrt(sum(c1.^2 + s1.^2));
    amp2(k) = sqrt(sum(c2.^2 + s2.^2));
    diffAmp(k) = sqrt(sum((c1-c2).^2 + (s1-s2).^2));
    if amp1(k) > 0 && amp2(k) > 0
        degCorr(k) = (sum(c1.*c2) + sum(s1.*s2)) / (amp1(k) * amp2(k));
    end
end
ncross = NaN;
ix = find(diffAmp(3:end) > amp1(3:end), 1);      % start at degree 2
if ~isempty(ix), ncross = n(ix + 2); end
spec = struct('n', n, 'diffAmp', diffAmp, 'amp1', amp1, 'amp2', amp2, ...
    'degCorr', degCorr, 'ncross', ncross);
end
