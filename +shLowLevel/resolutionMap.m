function res = resolutionMap(op, t, latDeg, lonDeg, opts)
%RESOLUTIONMAP Filter-kernel half-widths: per-point, per-azimuth resolution.
%
%   RES = shLowLevel.resolutionMap(OP, T, LATDEG, LONDEG) evaluates, for each
%   query point, the effective smoothing kernel of month T:
%       xf(P) = y_P' * W_t * x   =>   kernel coefficients k = W_t' * y_P
%   and finds the great-circle distance psi_1/2 at which the kernel drops
%   to half its peak, along NAz (8) azimuths. Replaces the single "equivalent
%   Gaussian radius" with a spatially and azimuthally resolved resolution
%   product (N-S vs E-W half-width ratio = striping fingerprint).
%   All transect evaluations per query point run through ONE vectorized
%   shLowLevel.ylm call.
%
%   Inputs
%     op       struct     from shLowLevel.tvANSFilter
%     t        (1,1)      month index
%     latDeg   (nPts,1)   query latitudes  [deg]
%     lonDeg   (nPts,1)   query longitudes [deg]
%     opts.NAz    (1,1)   azimuth count, default 8
%     opts.PsiMax (1,1)   radial search range [deg], default 20
%     opts.NPsi   (1,1)   radial samples, default 120
%   Outputs
%     res  (nlat x nlon) double  local resolution [km] of the filter
%
%   Claude (Fable 5), 2026-08-07.
%   Outputs
%     res        (nlat x nlon) double   local effective resolution [km] of the stored filter operator
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    op (1,1) struct
    t (1,1) double {mustBeInteger, mustBePositive}
    latDeg (:,1) double
    lonDeg (:,1) double
    opts.NAz (1,1) double = 8
    opts.PsiMax (1,1) double = 20
    opts.NPsi (1,1) double = 120
end

idx = op.idx;
R   = 6371;
az  = (0:opts.NAz-1) * (360/opts.NAz);
psi = linspace(0, opts.PsiMax, opts.NPsi);
nP  = numel(latDeg);

psiHalf = nan(nP, opts.NAz);
for p = 1:nP
    lat0 = deg2rad(latDeg(p));
    lon0 = deg2rad(lonDeg(p));
    yP = shLowLevel.ylm(lat0, lon0, idx)';                   % P x 1
    k  = shLowLevel.opApply(op, yP, t, 'transp');
    k0 = k' * yP;

    % all transect points (NAz x NPsi) in one shot
    [AZ, PS] = ndgrid(deg2rad(az), deg2rad(psi));
    lat2 = asin(sin(lat0)*cos(PS) + cos(lat0)*sin(PS).*cos(AZ));
    lon2 = lon0 + atan2(sin(AZ).*sin(PS)*cos(lat0), ...
                        cos(PS) - sin(lat0)*sin(lat2));
    V = reshape(shLowLevel.ylm(lat2(:), lon2(:), idx) * k, size(AZ));  % NAz x NPsi

    for a = 1:opts.NAz
        below = find(V(a, 2:end) <= k0/2, 1) + 1;
        if ~isempty(below)
            q = below;
            f = (V(a, q-1) - k0/2) / max(V(a, q-1) - V(a, q), eps);
            psiHalf(p, a) = psi(q-1) + f * (psi(q) - psi(q-1));
        end
    end
end

kmAll = deg2rad(psiHalf) * R;
res.psiHalfDeg = psiHalf;
res.azDeg  = az;
res.kmMean = mean(kmAll, 2, 'omitnan');
res.kmMin  = min(kmAll, [], 2, 'omitnan');
res.kmMax  = max(kmAll, [], 2, 'omitnan');
res.anisotropy = res.kmMax ./ res.kmMin;
end
