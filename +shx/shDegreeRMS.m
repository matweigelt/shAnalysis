function spec = shDegreeRMS(C, S, varargin)
%SHDEGREERMS Spectral (degree-domain) analysis of spherical harmonic coefficients.
%
%   SPEC = SHDEGREERMS(C, S) computes, per degree n = 0..nmax:
%     .degVariance(n+1)  = sum_m (Cnm^2 + Snm^2)
%     .degRMS(n+1)       = sqrt(degVariance)
%     .degAmplitude(n+1) = R * degRMS            [geoid-height-equivalent, m]
%     .cumAmplitude(n+1) = R * sqrt(cumsum(degVariance))   [cumulative, m]
%
%   C, S are (nmax+1)x(nmax+1) coefficient matrices, C(n+1,m+1), fully
%   normalized (4-pi / ICGEM convention), lower-triangular (m<=n).
%
%   SPEC = SHDEGREERMS(C, S, 'sigmaC', sigC, 'sigmaS', sigS) additionally
%   computes the error degree amplitude .errAmplitude from formal errors,
%   using the same degree-RMS formula applied to sigC, sigS.
%
%   Name/value options:
%     'R'     reference radius [m], default 6378136.3
%     'n0'    lowest degree to include (coefficients below n0 zeroed out
%             before summation, e.g. n0=2 to drop degree 0/1), default 0
%
%   Claude (Sonnet 4.6), 2026-07-11; merged into +shx: Claude (Fable 5), 2026-08-07.
%   Outputs
%     spec       struct: degree (nmax+1 x 1), degVariance/degRMS/
%                degAmplitude, cumVariance/cumRMS/cumAmplitude, R,
%                domain; errVariance/errRMS/errAmplitude when sigmas
%                are supplied
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

p = inputParser;
addParameter(p, 'R', 6378136.3);
addParameter(p, 'n0', 0);
addParameter(p, 'sigmaC', []);
addParameter(p, 'sigmaS', []);
parse(p, varargin{:});
R = p.Results.R;
n0 = p.Results.n0;

nmax = size(C,1) - 1;
if n0 > 0
    C(1:n0,:) = 0;
    S(1:n0,:) = 0;
end

degVariance = zeros(nmax+1,1);
for n = 0:nmax
    degVariance(n+1) = sum(C(n+1,1:n+1).^2 + S(n+1,1:n+1).^2);
end
degRMS = sqrt(degVariance);

spec.degree = (0:nmax)';
spec.degVariance = degVariance;
spec.degRMS = degRMS;
spec.degAmplitude = R .* degRMS;
spec.cumVariance = cumsum(degVariance);              % v2.4.1
spec.cumRMS = sqrt(spec.cumVariance);                % v2.4.1
spec.cumAmplitude = R .* spec.cumRMS;
spec.R = R;
spec.domain = 'degree';                              % v2.4.1

if ~isempty(p.Results.sigmaC) && ~isempty(p.Results.sigmaS)
    sC = p.Results.sigmaC; sS = p.Results.sigmaS;
    sC(isnan(sC)) = 0; sS(isnan(sS)) = 0;
    if n0 > 0
        sC(1:n0,:) = 0; sS(1:n0,:) = 0;
    end
    errVariance = zeros(nmax+1,1);
    for n = 0:nmax
        errVariance(n+1) = sum(sC(n+1,1:n+1).^2 + sS(n+1,1:n+1).^2);
    end
    spec.errVariance = errVariance;                  % v2.4.1
    spec.errRMS = sqrt(errVariance);
    spec.errAmplitude = R .* spec.errRMS;
end

end
