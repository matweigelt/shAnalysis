function [Cf, Sf, Wn] = shGaussianFilter(C, S, radius_km, varargin)
%SHGAUSSIANFILTER Apply isotropic Gaussian smoothing to Stokes coefficients.
%
%   [CF, SF] = SHGAUSSIANFILTER(C, S, RADIUS_KM) multiplies each degree n
%   row of C, S by the Jekeli (1981) weight Wn(n), see SHGAUSSIANWEIGHTS.
%   Standard pre-processing step for GRACE/GRACE-FO mass-change maps
%   before spatial synthesis, used to suppress high-degree noise.
%
%   [CF, SF, WN] = ... also returns the weight vector used.
%
%   Name/value options: 'R' (sphere radius km, default 6378.1363),
%   passed through to SHGAUSSIANWEIGHTS.
%
%   Claude (Sonnet 4.6), 2026-07-11; merged into +shx: Claude (Fable 5), 2026-08-07.
%   Outputs
%     Cf         (nmax+1 x nmax+1) double   smoothed cosine coefficients
%     Sf         (nmax+1 x nmax+1) double   smoothed sine coefficients
%     Wn         (nmax+1 x 1) double   Jekeli degree weights applied
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

nmax = size(C,1) - 1;
Wn = shx.shGaussianWeights(nmax, radius_km, varargin{:});

Cf = C .* Wn;   % broadcast over rows (degree), same weight across all orders in a row
Sf = S .* Wn;

end
