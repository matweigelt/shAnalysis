function [Cf, Sf] = shFanFilter(C, S, rDegKm, rOrdKm, opts)
%SHFANFILTER Han fan filter: separable degree x order Gaussian smoothing.
%
%   [CF, SF] = shx.shFanFilter(C, S, RDEGKM, RORDKM) applies the fan
%   filter of Han et al. (2005): the product of a Jekeli Gaussian in
%   DEGREE and one in ORDER,
%       W_nm = W_n(rDeg) * W_m(rOrd),
%   damping high orders (the striping direction) harder than an
%   isotropic Gaussian of equal degree radius - the cheap anisotropic
%   benchmark between Jekeli and DDK/tvANS. Common choice:
%   rOrd = rDeg (original fan) or rOrd < rDeg for stronger destriping.
%
%   Inputs
%     C, S     (n1,n1) double
%     rDegKm   (1,1) double  degree half-response radius [km]
%     rOrdKm   (1,1) double  order half-response radius [km]
%     opts.R   (1,1) double = 6378136.3  [m]
%   Outputs
%     Cf, Sf   (n1,n1) double  filtered coefficients
%
%   Class conveniences: g.fan(r1, r2), ts.fan(r1, r2).
%
%   Claude (Fable 5), 2026-08-07 (v2.4).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    C double
    S double
    rDegKm (1,1) double {mustBePositive}
    rOrdKm (1,1) double {mustBePositive}
    opts.R (1,1) double = 6378136.3
end
nmax = size(C, 1) - 1;
Wn = shx.shGaussianWeights(nmax, rDegKm, R = opts.R / 1e3);  % expects km
Wm = shx.shGaussianWeights(nmax, rOrdKm, R = opts.R / 1e3);
Wmat = Wn(:) .* Wm(:)';                             % (n+1) x (m+1)
Cf = C .* Wmat;
Sf = S .* Wmat;
end
