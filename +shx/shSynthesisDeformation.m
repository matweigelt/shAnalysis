function [up, north, east] = shSynthesisDeformation(C, S, R, latDeg, lonDeg, opts)
%SHSYNTHESISDEFORMATION Elastic load deformation: up / north / east.
%
%   [UP, NORTH, EAST] = shx.shSynthesisDeformation(C, S, R, LAT, LON,
%       kn=..., hn=..., ln=...) synthesizes the three components of the
%   elastic surface deformation caused by the load expressed in the
%   (residual) Stokes coefficients C, S (Wahr et al. 1998; Kusche &
%   Schrama 2005 - the GRACE <-> GNSS loading comparison):
%
%     up    = R * sum  hn/(1+kn) * dC * Ybar
%     north = R * sum  ln/(1+kn) * dC * dYbar/dphi / R * R  (= tangential)
%     east  = R * sum  ln/(1+kn) * dC * (1/cos phi) * dYbar/dlambda
%
%   Vertical needs only a degree factor (also available as
%   quantity="deformation_up" in shSynthesis); the horizontal components
%   need the exact spherical-harmonic gradient - dPbar/dphi comes from
%   shx.legendreALFDeriv (frozen validated identity), the longitude
%   derivative swaps and m-scales the trigonometric parts. Python-
%   validated point-wise against numerical gradients of the scalar
%   synthesis: north/east relative errors ~5e-10 (finite-difference
%   floor).
%
%   Inputs
%     C, S    (n1,n1) double        RESIDUAL Stokes coefficients (mean
%                                   field removed - deformation from the
%                                   static field is not meaningful)
%     R       (1,1) double [m]
%     latDeg  (1,nlat) double [deg] geocentric latitude
%     lonDeg  (1,nlon) double [deg]
%     opts.kn (>=n1,1) double       load Love numbers k'_n  (REQUIRED)
%     opts.hn (>=n1,1) double       load Love numbers h'_n  (REQUIRED)
%     opts.ln (>=n1,1) double       load Love numbers l'_n  (REQUIRED)
%     opts.nmin (1,1) double = 1    degrees below nmin are zeroed
%                                   (degree 0 carries no deformation;
%                                   include degree 1 ONLY if your
%                                   coefficients carry geocenter motion
%                                   in the CF frame consistently)
%     opts.Mode "grid"|"points" = "grid"
%                                   grid: nlat x nlon outer product;
%                                   points: equal-length LAT/LON vectors
%                                   evaluated pairwise (GNSS stations)
%   Outputs [m]
%     up, north, east   (nlat,nlon) double (grid) or (np,1) (points);
%                       EAST is NaN at |lat| = 90 (undefined direction)
%
%   All Love numbers are user-supplied - never hardcoded. Ensure kn, hn,
%   ln come from the SAME loading model (e.g. one PREM table) and the
%   same (CE/CF) frame convention for degree 1.
%
%   Claude (Fable 5), 2026-08-07.
%   Outputs
%     up         (nlat x nlon | npts) double   vertical elastic deformation [m]
%     north      same size   horizontal north component [m]
%     east       same size   horizontal east component [m]
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    C double
    S double
    R (1,1) double {mustBePositive}
    latDeg (1,:) double
    lonDeg (1,:) double
    opts.kn double
    opts.hn double
    opts.ln double
    opts.nmin (1,1) double {mustBeInteger, mustBeNonnegative} = 1
    opts.Mode (1,1) string {mustBeMember(opts.Mode, ["grid","points"])} = "grid"
end
n1 = size(C, 1);
nmax = n1 - 1;
assert(isequal(size(C), size(S)), 'shSynthesis:badInput', ...
    'C and S must be equally sized square matrices.');
for fld = ["kn", "hn", "ln"]
    v = opts.(fld);
    assert(numel(v) >= n1, 'shSynthesis:missingLoveNumbers', ...
        'Love numbers %s must cover degrees 0..%d (got %d values).', ...
        fld, nmax, numel(v));
end
pts = opts.Mode == "points";
if pts
    assert(numel(latDeg) == numel(lonDeg), 'shSynthesis:badInput', ...
        'Mode="points" needs equal-length LAT and LON.');
end

nn = (0:nmax)';
fU = R * opts.hn(1:n1) ./ (1 + opts.kn(1:n1));
fH = R * opts.ln(1:n1) ./ (1 + opts.kn(1:n1));
fU(1:min(opts.nmin, n1)) = 0;
fH(1:min(opts.nmin, n1)) = 0;
fU = fU(:); fH = fH(:);

KCu = fU .* C; KSu = fU .* S;
KCh = fH .* C; KSh = fH .* S;

latRad = deg2rad(latDeg(:)');
lonRad = deg2rad(lonDeg(:)');
nlat = numel(latRad); nlon = numel(lonRad);
[P, D] = shx.legendreALFDeriv(nmax, latRad);

m = (0:nmax)';
if pts
    up = zeros(nlat, 1); north = up; east = up;
    for k = 1:nlat
        Pk = P(:, :, k); Dk = D(:, :, k);
        cM = cos(m * lonRad(k)); sM = sin(m * lonRad(k));
        up(k)    = (sum(KCu .* Pk, 1) * cM) + (sum(KSu .* Pk, 1) * sM);
        north(k) = (sum(KCh .* Dk, 1) * cM) + (sum(KSh .* Dk, 1) * sM);
        cl = cos(latRad(k));
        if abs(cl) < 1e-12
            east(k) = NaN;
        else
            east(k) = ((sum(KCh .* Pk, 1) * (-m .* sM)) ...
                     + (sum(KSh .* Pk, 1) * ( m .* cM))) / cl;
        end
    end
else
    cosML = cos(m * lonRad);                        % (nmax+1) x nlon
    sinML = sin(m * lonRad);
    up = zeros(nlat, nlon); north = up; east = up;
    for k = 1:nlat
        Pk = P(:, :, k); Dk = D(:, :, k);
        Au = sum(KCu .* Pk, 1); Bu = sum(KSu .* Pk, 1);
        An = sum(KCh .* Dk, 1); Bn = sum(KSh .* Dk, 1);
        Ae = sum(KCh .* Pk, 1); Be = sum(KSh .* Pk, 1);
        up(k, :)    = Au * cosML + Bu * sinML;
        north(k, :) = An * cosML + Bn * sinML;
        cl = cos(latRad(k));
        if abs(cl) < 1e-12
            east(k, :) = NaN;
        else
            east(k, :) = (Ae * (-(m .* sinML)) + Be * (m .* cosML)) / cl;
        end
    end
end
end
