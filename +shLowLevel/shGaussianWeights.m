function Wn = shGaussianWeights(nmax, radius_km, varargin)
%SHGAUSSIANWEIGHTS Degree-domain weights for isotropic Gaussian smoothing.
%
%   WN = SHGAUSSIANWEIGHTS(NMAX, RADIUS_KM) computes the Jekeli (1981)
%   Gaussian averaging function weights Wn, n = 0..NMAX, for an averaging
%   radius RADIUS_KM [km] (half-width at which the spatial weighting
%   function drops to 1/2). Standard GRACE/GRACE-FO mass-change smoothing.
%
%   WN = SHGAUSSIANWEIGHTS(..., 'R', R_km) sets the sphere radius used in
%   b = ln(2)/(1-cos(radius/R)), default R = 6378.1363 km.
%
%   Recursion (Jekeli 1981 / Wahr et al. 1998):
%     b   = ln(2) / (1 - cos(radius/R))
%     W0  = 1
%     W1  = (1+exp(-2b))/(1-exp(-2b)) - 1/b
%     W_{n+1} = -(2n+1)/b * Wn + W_{n-1}
%
%   Validated (see validate_filters.py) against direct numerical expansion
%   of the spatial-domain Gaussian kernel in Legendre polynomials:
%   max abs difference 2.4e-14 over n=0..60, radius=300 km.
%
%   Usage: [Cf, Sf] = shGaussianFilter(C, S, radius_km) applies these
%   weights to coefficient matrices directly.
%
%   Claude (Sonnet 4.6), 2026-07-11; merged into +shLowLevel: Claude (Fable 5), 2026-08-07.
%   Outputs
%     Wn         (nmax+1 x 1) double   Jekeli (1981) weights, W(1) = 1, monotone non-increasing
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

p = inputParser;
addParameter(p, 'R', 6378.1363);
parse(p, varargin{:});
R = p.Results.R;

b = log(2) / (1 - cos(radius_km / R));

Wn = zeros(nmax+1, 1);
Wn(1) = 1;
if nmax >= 1
    Wn(2) = (1+exp(-2*b))/(1-exp(-2*b)) - 1/b;
end
for n = 1:nmax-1
    Wn(n+2) = -(2*n+1)/b * Wn(n+1) + Wn(n);
    % The forward three-term recursion is unstable once the weights have
    % decayed to (relative) machine noise -- for large averaging radii or
    % very high degrees it then oscillates with exponentially growing
    % amplitude. The true kernel is positive and monotonically
    % non-increasing, so the first violation marks the breakdown point:
    % clamp everything from there on to 0 (standard practical safeguard,
    % cf. Wahr et al. 1998 usage of the Jekeli recursion).
    if Wn(n+2) < 0 || Wn(n+2) > Wn(n+1)
        Wn(n+2:end) = 0;
        break
    end
end

end
