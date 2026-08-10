function [G, lam, info] = slepianBasis(idx, region, opts)
%SLEPIANBASIS Spatially concentrated (Slepian) basis for a region.
%
%   [G, LAM, INFO] = shLowLevel.slepianBasis(IDX, REGION) solves the spherical
%   concentration problem of Simons et al. (2006): find band-limited
%   functions g = Y*x maximizing the energy fraction inside REGION,
%       lambda = (x' K x) / (x' x),   K = Y' diag(w .* mask) Y,
%   by symmetric eigendecomposition of the localization kernel K
%   (exact Gauss-Legendre quadrature). Columns of G are the Slepian
%   taper coefficient vectors, orthonormal (G'G = I), sorted by
%   concentration LAM (in [0,1], descending). The Shannon number
%   N = sum(LAM) = areaFraction * P tells how many tapers are usefully
%   concentrated - regional analysis then estimates only ~N coefficients
%   instead of P, the well-posed alternative to Kaula-regularized
%   least squares on regional grids.
%
%   Typical use (regional estimate from global coefficients x):
%       a  = G(:,1:K)' * x;          % Slepian-domain coefficients
%       xR = G(:,1:K) * a;           % regional reconstruction
%
%   Inputs
%     idx     struct  shLowLevel.shIndex (P x P eig: fine to Lmax ~ 60-96)
%     region  handle | polygon | mask   (see shLowLevel.evalMask)
%   Options
%     NKeep (round(Shannon))  columns of G to return; [] = all P
%     BufferKm (0), R (6378136.3)  passed to shLowLevel.evalMask
%   Outputs
%     G     (P,NKeep) double  taper coefficient vectors, G'G = I
%     lam   (NKeep,1) double  concentrations, descending in [0,1]
%     info  struct: shannon (= areaFraction*P, exact by the addition
%           theorem), areaFraction, lamAll (P x 1)
%
%   Validated (Python, 30-deg polar cap, Lmax=12): lam in [0,1] to 4e-16,
%   G'G-I ~ 3e-15, trace(K) = quadrature-area * P exactly.
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    idx (1,1) struct
    region
    opts.NKeep double = []
    opts.BufferKm (1,1) double = 0
    opts.R (1,1) double = 6378136.3
    opts.OverSample (1,1) double {mustBeInteger, mustBePositive} = 2
end
[mask, ~] = shLowLevel.evalMask(idx, region, BufferKm = opts.BufferKm, ...
    R = opts.R, OverSample = opts.OverSample);
[Y, w] = shLowLevel.synthesisMatrix(idx, ...
    NLat = opts.OverSample * (idx.Lmax + 1), ...
    NLon = opts.OverSample * (2 * idx.Lmax + 2));
K = Y' * ((w .* mask) .* Y);
K = (K + K') / 2;
[G, D] = eig(K);
lamAll = min(max(real(diag(D)), 0), 1);
[lamAll, ord] = sort(lamAll, 'descend');
G = G(:, ord);
areaFrac = w' * mask;
shannon = areaFrac * idx.P;
nk = opts.NKeep;
if isempty(nk), nk = max(1, round(shannon)); end
nk = min(nk, idx.P);
info = struct('shannon', shannon, 'areaFraction', areaFrac, 'lamAll', lamAll);
lam = lamAll(1:nk);
G = G(:, 1:nk);
end
