function [S, grid, info] = seaLevelFingerprint(loadRegion, ocean, idx, opts)
%SEALEVELFINGERPRINT Elastic sea-level fingerprint (sea-level equation).
%
%   [S, GRID, INFO] = shLowLevel.seaLevelFingerprint(LOADREGION, OCEAN, IDX,
%       kn=kn, hn=hn) solves the elastic sea-level equation for a land
%   load: the ocean does NOT rise uniformly - it follows the perturbed
%   geoid minus the deformed crust under global mass conservation
%   (Farrell & Clark 1976, elastic limit):
%
%       S = O * ( N - U + kappa ),
%       N, U from the TOTAL load (land + ocean) via load Love numbers,
%       kappa from  integral( rho_w S ) = -integral( sigma_land )
%
%   iterated to convergence on the toolbox quadrature grid (exact
%   Gauss-Legendre analysis/synthesis pair, roundtrip 1.5e-13).
%   Python-validated: mass conserved to machine precision, ~11
%   iterations, and the classic pattern - sea level FALLS in the near
%   field of a melting mass (lost self-attraction + crustal rebound)
%   and rises above eustatic in the far field.
%
%   Inputs
%     loadRegion  handle f(latDeg,lonDeg) -> surface density CHANGE
%                 [kg/m^2] (negative = mass loss), polygon (Kx2 with
%                 opts.LoadValue), or Ngrid x 1 values on the solver grid
%     ocean       ocean region: handle | polygon | mask (shLowLevel.evalMask
%                 forms). The load is restricted to land (1 - ocean).
%     idx         shLowLevel.shIndex(Lmax, MinDegree = 0)  - MinDegree MUST be
%                 0 (mass conservation lives in degree 0/1)
%   Options
%     kn, hn (REQUIRED)    load Love numbers k'_n, h'_n (user-supplied)
%     LoadValue (1)        kg/m^2 for polygon loads
%     rho_water (1000), rho_ave (5517), R (6378136.3)
%     OverSample (2)       quadrature refinement (mask fidelity)
%     MaxIter (50), Tol (1e-8 * |eustatic|)  convergence on max|dS|
%   Outputs
%     S     (Ngrid,1) double  relative sea-level change over the ocean
%                             [m per unit load], 0 on land
%     grid  struct            latDeg/lonDeg rings of the solver grid
%     info  struct: iterations, dS, massResidual [kg/m^2 sphere mean],
%           eustatic [m], kappa, areaFractionOcean, S2D (nlat x nlon)
%
%   Limitations (documented): elastic only (no viscous GIA response),
%   no rotational feedback, fixed coastlines, band-limited load/ocean
%   masks (Gibbs at coasts - OverSample and Lmax control it).
%
%   Claude (Fable 5), 2026-08-07 (v2.4).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    loadRegion
    ocean
    idx (1,1) struct
    opts.kn double
    opts.hn double
    opts.LoadValue (1,1) double = 1
    opts.rho_water (1,1) double = 1000
    opts.rho_ave (1,1) double = 5517
    opts.R (1,1) double = 6378136.3
    opts.OverSample (1,1) double {mustBeInteger, mustBePositive} = 2
    opts.MaxIter (1,1) double {mustBeInteger, mustBePositive} = 50
    opts.Tol (1,1) double = 1e-8
end
assert(idx.minDegree == 0, 'shLowLevel:seaLevelFingerprint:badIndex', ...
    'idx must use MinDegree = 0 (mass conservation needs degrees 0-1).');
n1 = idx.Lmax + 1;
for fld = ["kn", "hn"]
    assert(numel(opts.(fld)) >= n1, 'shSynthesis:missingLoveNumbers', ...
        'Love numbers %s must cover degrees 0..%d.', fld, idx.Lmax);
end
[Y, w, grid] = shLowLevel.synthesisMatrix(idx, ...
    NLat = opts.OverSample * n1, NLon = opts.OverSample * (2 * idx.Lmax + 2));
[O, ~] = shLowLevel.evalMask(idx, ocean, OverSample = opts.OverSample);

latRow = repelem(grid.latDeg(:), numel(grid.lonDeg));
lonRow = repmat(grid.lonDeg(:), numel(grid.latDeg), 1);
if isa(loadRegion, 'function_handle')
    sigL = double(loadRegion(latRow, lonRow));
elseif isnumeric(loadRegion) && size(loadRegion, 2) == 2 && size(loadRegion, 1) >= 3
    sigL = opts.LoadValue * ...
        double(inpolygon(lonRow, latRow, loadRegion(:, 2), loadRegion(:, 1)));
elseif isnumeric(loadRegion) && numel(loadRegion) == numel(O)
    sigL = double(loadRegion(:));
else
    error('shLowLevel:seaLevelFingerprint:badLoad', ...
        'LOADREGION must be a function handle, polygon, or Ngrid x 1 values.');
end
sigL = sigL .* (1 - O);                            % land only

nn = (0:idx.Lmax)';
fSD = opts.R * opts.rho_ave / 3 * (2*nn + 1) ./ (1 + opts.kn(1:n1));
fN  = opts.R * ones(n1, 1);
fU  = opts.R * opts.hn(1:n1) ./ (1 + opts.kn(1:n1));
fSDr = fSD(idx.n + 1); fNr = fN(idx.n + 1); fUr = fU(idx.n + 1);

Aoc = w' * O;
Mland = w' * sigL;                                 % sphere-mean kg/m^2
eust = -Mland / (opts.rho_water * Aoc);
S = eust * O;
dS = Inf; kap = eust;
tolAbs = opts.Tol * max(abs(eust), realmin);
for it = 1:opts.MaxIter
    sig = sigL + opts.rho_water * S;
    cSig = Y' * (w .* sig);                        % quadrature analysis
    cSt = cSig ./ fSDr;                            % Stokes
    NU = Y * ((fNr - fUr) .* cSt);                 % geoid - uplift
    kap = (-Mland / opts.rho_water - w' * (O .* NU)) / Aoc;
    Snew = O .* (NU + kap);
    dS = max(abs(Snew - S));
    S = Snew;
    if dS < tolAbs, break; end
end
massRes = w' * (sigL + opts.rho_water * S);
S2D = reshape(S, numel(grid.lonDeg), numel(grid.latDeg))';
info = struct('iterations', it, 'dS', dS, 'massResidual', massRes, ...
    'eustatic', eust, 'kappa', kap, 'areaFractionOcean', Aoc, ...
    'S2D', S2D, 'converged', dS < tolAbs);
end
