function [C, S, info] = shAnalysisGrid(grid, latVec, lonVec, nmax, opts)
%SHANALYSISGRID Spherical harmonic analysis: estimate Stokes coefficients
%   from gridded or scattered data (the inverse of shLowLevel.shSynthesis).
%
%   [C, S, INFO] = shLowLevel.shAnalysisGrid(GRID, LATVEC, LONVEC, NMAX, ...)
%
%   Two data layouts, auto-detected:
%     * ring grid:  GRID is numel(LATVEC) x numel(LONVEC), LONVEC uniform
%       covering the full circle. Longitude is diagonalized by FFT and one
%       small least-squares problem is solved per order over the latitude
%       rings ("rings" method) - fast and, for band-limited data with
%       enough rings (> NMAX-m per order) and nlon > 2*NMAX, EXACT
%       (round-trip validated to 8.7e-15 in Python). Any ring spacing
%       works: equiangular, Gauss-Legendre, irregular.
%     * scattered points: GRID, LATVEC, LONVEC equal-length vectors of
%       point values/positions - full least squares on (NMAX+1)^2
%       unknowns ("ls" method), normal equations accumulated in chunks.
%
%   Name/value options:
%     Method    "auto" (default) | "rings" | "ls"
%     quantity  'geoid' (default) | 'potential' | 'gravity_anomaly' |
%               'gravity_disturbance' | 'ewh'  - the grid's physical
%               quantity; the degree-dependent kernel (shLowLevel.kernelFactors)
%               is divided out so C,S are Stokes coefficients. Degrees
%               with zero kernel (e.g. n=0,1 for gravity_anomaly) are
%               unobservable and returned as 0 (listed in INFO).
%     GM, R     defaults 3.986004415e14, 6378136.3 (only used via kernel)
%     kn ([])        load Love numbers, required for quantity='ewh'
%     rho_ave, rho_water   EWH densities (5517, 1000)
%     Weights   "none" (default) | "coslat" - area weighting of rings/
%               points (recommended for equiangular ring grids: rows are
%               scaled by sqrt(cos lat), a diagonal approximation of the
%               area weight; for complete ring systems the LS solution is
%               weight-independent, weighting matters when redundant/noisy)
%     Kaula     0 (default: off). If > 0, Tikhonov regularization toward 0
%               with per-degree standard deviation Kaula*1e-5/n^2 (Kaula's
%               rule; Kaula=1 is the classic value). REQUIRED when the
%               sampling underdetermines the coefficients; without it,
%               rank deficiency raises shLowLevel:shAnalysisGrid:rankDeficient.
%     ChunkSize points per accumulation chunk for "ls" (default 2000)
%
%   Outputs
%     C     (nmax+1 x nmax+1) double  estimated cosine coefficients
%     S     (nmax+1 x nmax+1) double  estimated sine coefficients
%     info  (1,1) struct  fields: condest (1,1 double), nObs (1,1
%           double), kaula (1,1 double, NaN when unregularized)
%
%   Notes: analysis of an incomplete/regional grid is an ill-posed
%   problem; expect leakage across degrees unless Kaula (or restricting
%   NMAX to what the sampling supports) is used. S_n0 does not exist
%   (sin 0 = 0) and is always returned 0.
%
%   Claude (Fable 5), 2026-08-07.
%
%   Options
%     hn ([])  vertical-deformation Love numbers, degrees 0..nmax (user-supplied)
%     LatType ("geocentric")  "geocentric" | "geodetic": how the
%         latitude inputs are interpreted; geodetic values are
%         converted with Flattening before evaluation
%     Flattening (1/298.257223563)  flattening of the reference
%         ellipsoid for the geodetic/geocentric conversion (WGS84;
%         overridable, never silently assumed)
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    grid double
    latVec double
    lonVec double
    nmax (1,1) double {mustBeInteger, mustBeNonnegative}
    opts.Method (1,1) string {mustBeMember(opts.Method, ["auto","rings","ls"])} = "auto"
    opts.quantity {mustBeTextScalar} = 'geoid'
    opts.GM (1,1) double = 3.986004415e14
    opts.R (1,1) double = 6378136.3
    opts.kn double = []
    opts.hn double = []
    opts.rho_ave (1,1) double = 5517
    opts.rho_water (1,1) double = 1000
    opts.Weights (1,1) string {mustBeMember(opts.Weights, ["none","coslat"])} = "none"
    opts.Kaula (1,1) double {mustBeNonnegative} = 0
    opts.ChunkSize (1,1) double {mustBePositive} = 2000
    opts.LatType (1,1) string ...
        {mustBeMember(opts.LatType, ["geocentric","geodetic"])} = "geocentric"
    opts.Flattening (1,1) double = 1/298.257223563
end
if opts.LatType == "geodetic"
    % user grids (maps, mascons) are typically geodetic; the estimation
    % operates in geocentric latitude like all SH math here (v2.2)
    latVec = shLowLevel.geodetic2geocentric(latVec, Flattening = opts.Flattening);
end


lat = latVec(:); lon = lonVec(:);
kernel = shLowLevel.kernelFactors(opts.quantity, nmax, opts.GM, opts.R, ...
    kn = opts.kn, hn = opts.hn, rho_ave = opts.rho_ave, ...
    rho_water = opts.rho_water);
unobs = find(kernel == 0) - 1;                 % degrees carrying no signal

% ---- layout detection
isRing = ismatrix(grid) && isequal(size(grid), [numel(lat), numel(lon)]);
isScat = ~isRing && isvector(grid) && numel(grid) == numel(lat) ...
    && numel(lat) == numel(lon);
if ~isRing && ~isScat
    error('shLowLevel:shAnalysisGrid:badInput', ...
        ['GRID must be numel(lat) x numel(lon) (ring grid) or a vector ' ...
         'matching lat/lon point lists.']);
end
method = opts.Method;
if method == "auto"
    if isRing, method = "rings"; else, method = "ls"; end
end
if method == "rings" && ~isRing
    error('shLowLevel:shAnalysisGrid:badInput', 'Method="rings" needs a ring grid.');
end
if method == "ls" && isRing                    % flatten the ring grid
    [LO, LA] = meshgrid(lon, lat);
    lat = LA(:); lon = LO(:); grid = grid(:);
    isScat = true; %#ok<NASGU>
end

latRad = deg2rad(lat); lonRad = deg2rad(lon);
C = zeros(nmax+1); S = zeros(nmax+1);
info = struct('method', char(method), 'nPoints', numel(grid), ...
    'unobservedDegrees', unobs(:)', 'regularized', opts.Kaula > 0);

if method == "rings"
    nlat = numel(lat); nlon = numel(lon);
    dl = diff(lonRad);
    if nlon < 2 || any(abs(dl - dl(1)) > 1e-9) || dl(1) <= 0 ...
            || abs(nlon*dl(1) - 2*pi) > 1e-6
        error('shLowLevel:shAnalysisGrid:badGrid', ...
            'rings method requires a uniform full-circle longitude vector.');
    end
    if nlon <= 2*nmax
        error('shLowLevel:shAnalysisGrid:badGrid', ...
            'nlon=%d aliases orders for nmax=%d; need nlon > 2*nmax (or Method="ls").', ...
            nlon, nmax);
    end
    P = shLowLevel.legendreALF(nmax, latRad);         % (n+1)x(m+1)xnlat
    w = ones(nlat, 1);
    if opts.Weights == "coslat", w = sqrt(max(cos(latRad), 0)); end
    % FFT across longitude; account for the grid's first longitude
    G = fft(grid, [], 2) / nlon;               % nlat x nlon complex
    condMax = 0; res0 = 0; resM = 0;
    for m = 0:nmax
        ph = exp(-1i * m * lonRad(1));
        gm = G(:, m+1) * ph;
        if m == 0
            am = real(gm); bm = [];
        else
            am = 2*real(gm); bm = -2*imag(gm);
        end
        nn = m:nmax;
        A = zeros(nlat, numel(nn));
        for jj = 1:numel(nn)
            A(:, jj) = squeeze(P(nn(jj)+1, m+1, :)) * kernel(nn(jj)+1);
        end
        est = kernel(nn+1) ~= 0;
        Ae = A(:, est) .* w;
        [xc, cnd] = solveLS(Ae, w .* am, nn(est), opts.Kaula);
        C(nn(est)+1, m+1) = xc;
        condMax = max(condMax, cnd);
        if m == 0
            res0 = sum((am - A(:, est) * xc).^2);
        else
            resM = resM + sum((am - A(:, est) * xc).^2);
            xs = solveLS(Ae, w .* bm, nn(est), opts.Kaula);
            S(nn(est)+1, m+1) = xs;
            resM = resM + sum((bm - A(:, est) * xs).^2);
        end
    end
    info.condEst = condMax;
    % Parseval: sum_j f_j^2 = nlon*(a0^2 + sum_{m>0}(am^2+bm^2)/2)
    info.residRMS = sqrt((nlon*(res0 + resM/2)) / max(nlat*nlon, 1));
else
    % ---- scattered least squares via chunked normal equations
    if nmax < 1
        error('shLowLevel:shAnalysisGrid:badInput', ...
            'Method="ls" requires nmax >= 1 (use "rings" for nmax=0).');
    end
    idx = shLowLevel.shIndex(nmax, MinDegree = 0);
    est = kernel(idx.n + 1) ~= 0;
    Pn = nnz(est);
    ATA = zeros(Pn); ATb = zeros(Pn, 1);
    w2 = ones(numel(grid), 1);
    if opts.Weights == "coslat", w2 = max(cos(latRad), 0); end
    npts = numel(grid);
    for i0 = 1:opts.ChunkSize:npts
        ii = i0:min(i0 + opts.ChunkSize - 1, npts);
        Y = shLowLevel.ylm(latRad(ii), lonRad(ii), idx);       % npts_chunk x P
        Y = Y .* kernel(idx.n + 1)';
        Y = Y(:, est);
        ATA = ATA + Y' * (w2(ii) .* Y);
        ATb = ATb + Y' * (w2(ii) .* grid(ii));
    end
    if opts.Kaula > 0
        sigK = opts.Kaula * 1e-5 ./ max(idx.n(est), 1).^2;
        ATA = ATA + diag(1 ./ sigK.^2);
    end
    [R_, flag] = chol(ATA);
    if flag ~= 0
        error('shLowLevel:shAnalysisGrid:rankDeficient', ...
            ['Normal equations are singular: the sampling does not ' ...
             'determine all %d coefficients. Reduce nmax or set Kaula>0.'], Pn);
    end
    rc = rcond(ATA);
    if rc < 1e-13 && opts.Kaula == 0
        error('shLowLevel:shAnalysisGrid:rankDeficient', ...
            ['Normal equations numerically singular (rcond=%.1e). ' ...
             'Reduce nmax or set Kaula>0.'], rc);
    end
    x = R_ \ (R_' \ ATb);
    xf = zeros(idx.P, 1); xf(est) = x;
    [C, S] = shLowLevel.csFromVec(xf, idx);
    % residual RMS (second pass, chunked)
    residSS = 0;
    for i0 = 1:opts.ChunkSize:npts
        ii = i0:min(i0 + opts.ChunkSize - 1, npts);
        Y = shLowLevel.ylm(latRad(ii), lonRad(ii), idx) .* kernel(idx.n + 1)';
        residSS = residSS + sum((grid(ii) - Y * xf).^2);
    end
    info.residRMS = sqrt(residSS / npts);
    info.condEst = 1 / rc;
end
end

function [x, cnd] = solveLS(A, b, degs, kaula)
% per-order solve: plain LS, or Kaula-regularized normal equations
if isempty(A), x = zeros(0, 1); cnd = 0; return; end
if kaula > 0
    sigK = kaula * 1e-5 ./ max(degs(:), 1).^2;
    Nq = A'*A + diag(1 ./ sigK.^2);
    x = Nq \ (A' * b);
    cnd = cond(Nq);
else
    if rank(A) < size(A, 2)
        error('shLowLevel:shAnalysisGrid:rankDeficient', ...
            ['Rank-deficient per-order system (order block %d..%d): too few ' ...
             'latitude rings. Reduce nmax or set Kaula>0.'], degs(1), degs(end));
    end
    x = A \ b;
    cnd = cond(A);
end
end
