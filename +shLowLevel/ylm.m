function Y = ylm(latRad, lonRad, idx)
%YLM 4-pi normalized SH basis vectors at arbitrary points.
%
%   Y = shLowLevel.ylm(LATRAD, LONRAD, IDX) evaluates the basis at nPts points
%   (geocentric latitude/longitude in radians, equal-length vectors), so
%   that f(points) = Y * x for a coefficient vector x in IDX ordering.
%   Vectorized over points via one shLowLevel.legendreALF call.
%
%   Inputs
%     latRad  (nPts,1) double   geocentric latitudes  [rad]
%     lonRad  (nPts,1) double   longitudes            [rad]
%     idx     struct            from shLowLevel.shIndex
%   Outputs
%     Y       (nPts x idx.P) double
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    latRad (:,1) double
    lonRad (:,1) double
    idx (1,1) struct
end
assert(numel(latRad) == numel(lonRad), 'shLowLevel:ylm:sizeMismatch', ...
    'latRad and lonRad must have equal length.');

nPts = numel(latRad);
Pall = shLowLevel.legendreALF(idx.Lmax, latRad);         % (L+1)x(L+1)xnPts
lin  = sub2ind([idx.Lmax+1, idx.Lmax+1], idx.n+1, idx.m+1);

iC = idx.cs == 0;
Y = zeros(nPts, idx.P);
for k = 1:nPts
    Pb = Pall(:, :, k);
    pl = Pb(lin)';
    trig = zeros(1, idx.P);
    trig(iC)  = cos(idx.m(iC)'  * lonRad(k));
    trig(~iC) = sin(idx.m(~iC)' * lonRad(k));
    Y(k, :) = pl .* trig;
end
end
