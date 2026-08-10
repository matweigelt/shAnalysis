function [sig, latDeg, lonDeg] = errorMap(M, idx, latDeg, lonDeg, opts)
%ERRORMAP Analytic formal-error maps from a full coefficient covariance.
%
%   [SIG, LAT, LON] = shLowLevel.errorMap(M, IDX, LAT, LON, quantity="ewh",
%       kn=kn) computes the pointwise 1-sigma of a synthesized field
%   from a full P x P covariance (e.g. shLowLevel.readSINEX(...,
%   Output="covariance", Index=IDX)):
%       sigma^2(lat,lon) = a' M a,   a = kernel- and Y-weighted row.
%   Exact and fast via Cholesky: sigma^2 = ||L' a||^2 with M = L L',
%   assembled per latitude for all longitudes at once - the analytic
%   counterpart (and cross-check) of shLowLevel.mcPropagate for maps.
%
%   Inputs
%     M       (P,P) double   covariance in IDX ordering (symmetric PSD)
%     idx     struct         shLowLevel.shIndex (analysis-style MinDegree
%                            matching M)
%     latDeg  (1,nlat), lonDeg (1,nlon) double [deg] geocentric
%   Options
%     quantity ("geoid"), kn, hn, GM (3.986004415e14), R (6378136.3),
%     Height (0)  - as in shLowLevel.kernelFactors
%   Outputs
%     sig  (nlat,nlon) double  1-sigma of the field [unit of quantity]
%
%   Cost: O(nlat * nlon * P^2) after one Cholesky; at n96 (P ~ 9400)
%   budget minutes for a 1-degree grid - truncate M/idx for quick looks.
%
%   Claude (Fable 5), 2026-08-07 (v2.4).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    M double
    idx (1,1) struct
    latDeg (1,:) double
    lonDeg (1,:) double
    opts.quantity (1,1) string = "geoid"
    opts.kn double = []
    opts.hn double = []
    opts.GM (1,1) double = 3.986004415e14
    opts.R (1,1) double = 6378136.3
    opts.Height (1,1) double = 0
end
P = idx.P;
assert(isequal(size(M), [P P]), 'shLowLevel:errorMap:badSize', ...
    'M must be P x P for the supplied idx (P = %d).', P);
Ms = (M + M') / 2;
[L, flag] = chol(Ms, 'lower');
if flag ~= 0
    % PSD-but-singular covariances: eigenvalue square root fallback
    [V, E] = eig(Ms);
    e = max(real(diag(E)), 0);
    L = V .* sqrt(e)';
end
nmax = idx.Lmax;
kernel = shLowLevel.kernelFactors(opts.quantity, nmax, opts.GM, opts.R, ...
    kn = opts.kn, hn = opts.hn, Height = opts.Height);
latRad = deg2rad(latDeg(:)');
lonRad = deg2rad(lonDeg(:)');
nlat = numel(latRad); nlon = numel(lonRad);
Pleg = shLowLevel.legendreALF(nmax, latRad);
m = (0:nmax)';
cosML = cos(m * lonRad); sinML = sin(m * lonRad);   % (nmax+1) x nlon
sig = zeros(nlat, nlon);
fk = kernel(idx.n + 1);                              % per-row factor
for k = 1:nlat
    Pk = Pleg(:, :, k);
    % row weights a(p) per lon: f * Pbar_{n,m} * trig_m
    base = fk .* Pk(sub2ind(size(Pk), idx.n + 1, idx.m + 1));  % P x 1
    isC = idx.cs == 0;
    A = zeros(P, nlon);
    A(isC, :) = base(isC) .* cosML(idx.m(isC) + 1, :);
    A(~isC, :) = base(~isC) .* sinML(idx.m(~isC) + 1, :);
    sig(k, :) = sqrt(sum((L' * A).^2, 1));
end
end
