function [Y, w, grid] = synthesisMatrix(idx, opts)
%SYNTHESISMATRIX Dense SH synthesis matrix on a Gauss-Legendre grid.
%
%   [Y, W, GRID] = shx.synthesisMatrix(IDX) returns Y (Ngrid x P) with
%   f = Y*x for coefficient vectors x in IDX ordering (shx.shIndex), and
%   quadrature weights W (Ngrid x 1), normalized such that the discrete
%   analysis A = Y' * diag(W) is EXACT for band-limited fields:
%   A*Y = eye(P) to machine precision (asserted in the test suite).
%
%   Built on the single toolbox Legendre engine shx.legendreALF
%   (vectorized over all Gauss-Legendre latitudes in one call): the GL
%   nodes t in cos(colatitude) map to geocentric latitude lat = asin(t).
%
%   Grid: NLat Gauss-Legendre latitude rings (default Lmax+1), NLon
%   equiangular longitudes (default 2*Lmax+2), latitude-major rows:
%   row = (i-1)*NLon + j.
%
%   Inputs
%     idx        struct   from shx.shIndex
%     opts.NLat  (1,1)    number of latitude rings   [default Lmax+1]
%     opts.NLon  (1,1)    number of longitudes       [default 2*Lmax+2]
%   Outputs
%     Y      (NLat*NLon x idx.P) double   synthesis matrix
%     w      (NLat*NLon x 1)     double   quadrature weights, sum(w) = 1
%     grid   struct: theta, lambda [rad], latDeg (NLat x 1), lonDeg (NLon x 1)
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    idx (1,1) struct
    opts.NLat (1,1) double {mustBeInteger, mustBePositive} = idx.Lmax + 1
    opts.NLon (1,1) double {mustBeInteger, mustBePositive} = 2*idx.Lmax + 2
end

[tGL, wGL] = shx.gaussLegendre(opts.NLat);
lat    = asin(tGL);                               % geocentric latitude [rad]
lambda = (0:opts.NLon-1)' * (2*pi/opts.NLon);

% all Legendre functions for all rings in ONE call to the unified engine
Pall = shx.legendreALF(idx.Lmax, lat);            % (L+1)x(L+1)xNLat
lin  = sub2ind([idx.Lmax+1, idx.Lmax+1], idx.n+1, idx.m+1);

% longitude trig block (identical for every ring)
Trig = zeros(opts.NLon, idx.P);
iC = idx.cs == 0;
Trig(:, iC)  = cos(lambda * idx.m(iC)');
Trig(:, ~iC) = sin(lambda * idx.m(~iC)');

Y = zeros(opts.NLat * opts.NLon, idx.P);
for i = 1:opts.NLat
    Pb = Pall(:, :, i);
    Y((i-1)*opts.NLon + (1:opts.NLon), :) = Trig .* Pb(lin)';
end

w = kron(wGL(:) / (2*opts.NLon), ones(opts.NLon, 1));   % sum(w) = 1

grid.theta  = pi/2 - lat;  grid.lambda = lambda;
grid.latDeg = rad2deg(lat);
grid.lonDeg = rad2deg(lambda);
end
