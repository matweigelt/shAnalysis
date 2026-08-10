function [comps, info] = shSynthesisGradientTensor(C, S, GM, R, latDeg, lonDeg, opts)
%SHSYNTHESISGRADIENTTENSOR Full gravity gradient tensor in the local NEU frame.
%
%   [COMPS, INFO] = shLowLevel.shSynthesisGradientTensor(C, S, GM, R, LAT, LON,
%       Height=250e3) synthesizes all six independent components of the
%   disturbing gravity gradient tensor (GOCE-style diagnostics) in the
%   local north-east-up frame at radius r = R + Height:
%
%     G_uu = T_rr
%     G_nn = T_pp/r^2 + T_r/r
%     G_ee = T_ll/(r^2 cos^2 p) - tan(p) T_p/r^2 + T_r/r
%     G_un = T_rp/r - T_p/r^2
%     G_ue = T_rl/(r cos p) - T_l/(r^2 cos p)
%     G_ne = T_pl/(r^2 cos p) + tan(p) T_l/(r^2 cos p)
%
%   (p = latitude, l = longitude). Validation (Python, v2.4):
%   Laplace invariant trace(G) = 0 at 7e-16 of max|G|; G_uu identical to
%   the 'gravity_gradient_rr' kernel route (3e-16); the angular second
%   derivatives T_pp / T_pl match finite differences at the FD floor.
%   Second latitude derivatives of the Legendre functions come from
%   shLowLevel.legendreALFDeriv (frozen stencil applied twice).
%
%   Inputs
%     C, S    (n1,n1) double        Stokes coefficients (use RESIDUAL
%                                   fields for time-variable work)
%     GM, R   (1,1) double
%     latDeg  (1,nlat), lonDeg (1,nlon) double [deg], geocentric
%   Options
%     Height (250e3) [m]  evaluation altitude above R
%     nmin (2)            zero degrees below nmin
%   Outputs
%     comps  struct of (nlat,nlon) double [1/s^2]: uu, nn, ee, un, ue,
%            ne  (symmetric tensor; 1 Eotvos = 1e-9 1/s^2). Components
%            involving east are NaN at |lat| = 90.
%     info   struct: maxTraceResidual (max |trace|/max|G| over the grid,
%            excluding poles) - the built-in Laplace self-check
%
%   Claude (Fable 5), 2026-08-07 (v2.4).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    C double
    S double
    GM (1,1) double
    R (1,1) double
    latDeg (1,:) double
    lonDeg (1,:) double
    opts.Height (1,1) double {mustBeNonnegative} = 250e3
    opts.nmin (1,1) double {mustBeInteger, mustBeNonnegative} = 2
end
n1 = size(C, 1); nmax = n1 - 1;
r = R + opts.Height;
nn = (0:nmax)';
att = (GM/R) * (R/r).^(nn + 1);
att(1:min(opts.nmin, n1)) = 0;

fT   = att;
fTr  = -(nn + 1) / r .* att;
fTrr = (nn + 1) .* (nn + 2) / r^2 .* att;

latRad = deg2rad(latDeg(:)');
lonRad = deg2rad(lonDeg(:)');
nlat = numel(latRad); nlon = numel(lonRad);
[P, D, D2] = shLowLevel.legendreALFDeriv(nmax, latRad);
m = (0:nmax)';
cosML = cos(m * lonRad); sinML = sin(m * lonRad);

Z = zeros(nlat, nlon);
comps = struct('uu', Z, 'nn', Z, 'ee', Z, 'un', Z, 'ue', Z, 'ne', Z);
trRes = 0;
for k = 1:nlat
    Pk = P(:,:,k); Dk = D(:,:,k); D2k = D2(:,:,k);
    row = @(F, T_) deal( sum((F .* C) .* T_, 1), sum((F .* S) .* T_, 1) );
    [aT,  bT ] = row(fT,   Pk);
    [aTr, bTr] = row(fTr,  Pk);
    [aRR, bRR] = row(fTrr, Pk);
    [aP,  bP ] = row(fT,   Dk);
    [aPP, bPP] = row(fT,   D2k);
    [aRP, bRP] = row(fTr,  Dk);
    Tr  = aTr * cosML + bTr * sinML;
    Trr = aRR * cosML + bRR * sinML;
    Tp  = aP  * cosML + bP  * sinML;
    Tpp = aPP * cosML + bPP * sinML;
    Trp = aRP * cosML + bRP * sinML;
    Tl  = aT  * (-(m .* sinML)) + bT * (m .* cosML);
    Tll = aT  * (-(m.^2 .* cosML)) + bT * (-(m.^2 .* sinML));
    Trl = aTr * (-(m .* sinML)) + bTr * (m .* cosML);
    Tpl = aP  * (-(m .* sinML)) + bP  * (m .* cosML);
    cl = cos(latRad(k)); tl = tan(latRad(k));
    comps.uu(k, :) = Trr;
    comps.nn(k, :) = Tpp / r^2 + Tr / r;
    comps.un(k, :) = Trp / r - Tp / r^2;
    if abs(cl) < 1e-12
        comps.ee(k, :) = NaN; comps.ue(k, :) = NaN; comps.ne(k, :) = NaN;
    else
        comps.ee(k, :) = Tll / (r^2 * cl^2) - tl * Tp / r^2 + Tr / r;
        comps.ue(k, :) = Trl / (r * cl) - Tl / (r^2 * cl);
        comps.ne(k, :) = Tpl / (r^2 * cl) + tl * Tl / (r^2 * cl);
        tr = comps.uu(k, :) + comps.nn(k, :) + comps.ee(k, :);
        trRes = max(trRes, max(abs(tr)));
    end
end
gmax = max(abs([comps.uu(:); comps.nn(:); comps.ee(:)]), [], 'omitnan');
info = struct('maxTraceResidual', trRes / max(gmax, realmin), ...
    'r', r, 'Height', opts.Height);
end
