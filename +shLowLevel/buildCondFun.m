function [f, info] = buildCondFun(idx, opts)
%BUILDCONDFUN Covariance conditioner in the EWH domain (Kvas 2019, 2.4).
%
%   [F, INFO] = shLowLevel.buildCondFun(IDX, kn=KN) builds the CondFun
%   handle for shLowLevel.estimateVAR that implements the covariance
%   conditioning of Kvas (2019), Sec. 2.4: the spectral covariance is
%   propagated to equivalent-water-height values on an exact
%   Gauss-Legendre grid, conditioned there, and propagated back:
%       Sigma~ = Fq * ( Z .* T .* (G Sigma G') ) * Fq'
%   with
%     G  = Y * D          spectral -> EWH grid (D the per-degree EWH
%                         kernel from shLowLevel.kernelFactors)
%     Fq = D^-1 * A       the EXACT inverse quadrature (A*Y = I is the
%                         synthesisMatrix contract), so Fq*G = I and
%                         with Z = 1, T = 1 the conditioner is the
%                         identity to machine precision (unit-tested)
%     Z  = region-block indicator, Z_ij = 1 iff points i, j share a
%                         region (his eq. 2.117: e.g. land/ocean -
%                         ocean-model errors are not correlated with
%                         continental hydrology)
%     T  = exp(-psi/psi0) distance-dependent correlation taper (his
%                         eq. 2.120), psi the spherical distance
%   Both Z (a Gram matrix of indicators) and T (an exponential kernel,
%   positive definite on the sphere) are PSD, so by the Schur product
%   theorem the conditioned covariance stays PSD; a singular empirical
%   Sigma(0) from a short series even becomes strictly positive
%   definite (the Hadamard product with a PD kernel is PD when the
%   diagonal is positive) - this is exactly why Kvas conditions: the
%   Yule-Walker solve stabilizes by orders of magnitude (his Fig. 2.5).
%
%   Physics: converting to EWH before windowing is deliberate - the
%   potential is global, but the GENERATING MASSES localize, so masking
%   and tapering are legitimate there (Kvas' argument). The EWH kernel
%   needs load Love numbers: kn is REQUIRED and user-supplied, never
%   hardcoded (fetchLoveNumbers / readLoveNumbers provide them).
%
%   Inputs
%     idx  (1 x 1) struct  from shLowLevel.shIndex - fixes Lmax and the
%          coefficient ordering the conditioner acts on
%
%   Options
%     kn (required)  (>= Lmax+1 x 1) double  load Love numbers k'_n,
%          kn(1) is degree 0 (kernelFactors contract)
%     Psi0Km (Inf)   (1 x 1) double  taper scale psi0 [km] of
%          exp(-psi/psi0); Inf disables the taper. Kvas' choice
%          corresponds to a half-width of ~1100 km
%     Regions ([])   region assignment for the block mask: a function
%          handle id = fun(latDeg, lonDeg) (array-valued, numeric ids;
%          e.g. @(la, lo) la >= 0) or an (M x 1) numeric vector on the
%          internal grid (INFO.grid, latitude-major rows as in
%          synthesisMatrix); [] disables masking
%     GM (3.986004415e14)  (1 x 1) double  [m^3/s^2] (overridable)
%     R (6378136.3)  (1 x 1) double  [m] reference radius; also
%          converts Psi0Km to spherical distance (geocentric sphere)
%     NLat ([])      (1 x 1) double  Gauss-Legendre rings (default
%          Lmax+1, the exactness minimum)
%     NLon ([])      (1 x 1) double  longitudes (default 2*Lmax+2)
%
%   Outputs
%     f     function handle  Sigma~ = f(Sigma) for (P x P) covariances
%           in IDX ordering - pass as estimateVAR(..., CondFun=f)
%     info  (1 x 1) struct  grid (synthesisMatrix grid struct), M
%           (1 x 1) number of grid points, region (M x 1 double or []),
%           W (M x M double) the elementwise spatial weight Z .* T,
%           memGB (1 x 1) storage of the precomputed operators
%
%   Example
%     kn  = shLowLevel.readLoveNumbers(lovePath);
%     idx = shLowLevel.shIndex(40);
%     land = @(la, lo) double(landMaskFun(la, lo));   % your land/ocean rule
%     cf  = shLowLevel.buildCondFun(idx, kn=kn, Psi0Km=1100, Regions=land);
%     model = shLowLevel.estimateVAR(X, Order=1, CondFun=cf);
%
%   Memory/cost: precomputes G, Fq (M x P each) and W (M x M); Lmax 40
%   gives M = 3362, ~180 MB total; every f(Sigma) costs two M x P x P
%   products. Numerics pre-validated in Python
%   (tools/dev/validate_condfun.py).
%
%   Reference: Kvas, TU Graz PhD thesis (2019), Sec. 2.4, eqs.
%   (2.117)-(2.120); Schur product theorem e.g. Horn & Johnson.
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-17, 21:30 UTC.

arguments
    idx (1,1) struct
    opts.kn double {mustBeNonempty}
    opts.Psi0Km (1,1) double {mustBePositive} = Inf
    opts.Regions = []
    opts.GM (1,1) double = 3.986004415e14
    opts.R (1,1) double = 6378136.3
    opts.NLat double = []
    opts.NLon double = []
end

% ---- exact synthesis/analysis pair on the Gauss-Legendre grid
nvY = {};
if ~isempty(opts.NLat), nvY = [nvY, {'NLat', opts.NLat}]; end
if ~isempty(opts.NLon), nvY = [nvY, {'NLon', opts.NLon}]; end
[Y, w, grd] = shLowLevel.synthesisMatrix(idx, nvY{:});
M = size(Y, 1);

% ---- per-degree EWH kernel (Love numbers user-supplied, never hardcoded)
kern = shLowLevel.kernelFactors('ewh', idx.Lmax, opts.GM, opts.R, kn = opts.kn);
d = kern(idx.n + 1);
if any(d == 0)
    error('shLowLevel:buildCondFun:zeroKernel', ...
        ['The EWH kernel vanishes for %d index entries (1 + kn = 0?): ' ...
         'the inverse quadrature is undefined there.'], nnz(d == 0));
end
G  = Y .* d.';                                 % spectral -> EWH grid
Fq = (Y .* (1 ./ d).').' .* w.';               % exact inverse: Fq * G = I

% ---- grid point coordinates (geocentric latitudes, toolbox convention;
%      latitude-major rows: row = (i-1)*NLon + j, as in synthesisMatrix)
la = reshape(deg2rad(repelem(grd.latDeg(:), numel(grd.lonDeg))), 1, []);
lo = reshape(deg2rad(repmat(grd.lonDeg(:), numel(grd.latDeg), 1)), 1, []);

% ---- elementwise spatial weight W = Z .* T
W = ones(M, 'like', Y);
if ~isempty(opts.Regions)
    if isa(opts.Regions, 'function_handle')
        region = double(opts.Regions(rad2deg(la(:)), rad2deg(lo(:))));
    else
        region = double(opts.Regions(:));
        if numel(region) ~= M
            error('shLowLevel:buildCondFun:badRegions', ...
                'Regions vector has %d entries, the grid has %d points.', ...
                numel(region), M);
        end
    end
    W = W .* double(region == region.');       % Gram of indicators: PSD
else
    region = [];
end
if isfinite(opts.Psi0Km)
    cosPsi = sin(la.') .* sin(la) + cos(la.') .* cos(la) .* cos(lo.' - lo);
    psiKm = acos(min(max(cosPsi, -1), 1)) * opts.R / 1000;
    W = W .* exp(-psiKm / opts.Psi0Km);        % PD kernel on the sphere
end

f = @(C) applyCond(C, G, Fq, W);
info = struct('grid', grd, 'M', M, 'region', region, 'W', W, ...
    'memGB', 8 * (2 * M * idx.P + M^2) / 1e9);
end

% ---------------------------------------------------------------- local
function Ct = applyCond(C, G, Fq, W)
S  = W .* (G * C * G.');
Ct = Fq * S * Fq.';
Ct = (Ct + Ct.') / 2;
end
