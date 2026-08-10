function val = synthesisPoints(C, S, GM, R, latDeg, lonDeg, r, opts)
%SYNTHESISPOINTS Gravity field quantities at arbitrary 3-D points.
%
%   VAL = shx.synthesisPoints(C, S, GM, R, LATDEG, LONDEG, RGEO)
%   evaluates the spherical harmonic expansion POINTWISE at positions
%   given by equal-length vectors of geocentric latitude [deg],
%   longitude [deg] and geocentric radius [m] - satellite altitudes,
%   scattered stations, profiles - including the upward continuation
%   (R/r)^n that grid synthesis on the sphere cannot provide:
%     potential    T   = GM/r     sum_n (R/r)^n     sum_m Pnm (C cos + S sin)
%     disturbance  dg  = GM/r^2   sum_n (n+1)(R/r)^n ...   ( = -dT/dr )
%     anomaly      Dg  = GM/r^2   sum_n (n-1)(R/r)^n ...
%   Fully normalized coefficients, C(n+1,m+1) layout, geocentric
%   latitudes (house conventions). Python-validated: pinned reference
%   values at the surface and at 450 km altitude, and the radial
%   derivative identity dg = -dT/dr to 3e-10 (tools/validate_sh
%   pointwise block; pins in testCorrectness/testSynthesisPoints).
%
%   Options
%     Quantity ("potential")  "potential" [m^2/s^2] | "disturbance" |
%                             "anomaly" (both [m/s^2]; multiply by 1e5
%                             for mGal)
%     MinDegree (2)           first degree included (0 evaluates the
%                             full field including the central term)
%
%   Outputs
%     val        (P,1) double  quantity at each point
%
%   Example
%     dg = shx.synthesisPoints(g.C, g.S, g.GM, g.R, ...
%         [10; 12], [100; 101], (g.R + 450e3) * [1; 1], ...
%         Quantity = "disturbance") * 1e5;   % mGal at 450 km
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-10 (v3.0.0).
arguments
    C double
    S double
    GM (1,1) double {mustBePositive}
    R (1,1) double {mustBePositive}
    latDeg (:,1) double {mustBeGreaterThanOrEqual(latDeg, -90), ...
        mustBeLessThanOrEqual(latDeg, 90)}
    lonDeg (:,1) double
    r (:,1) double {mustBePositive}
    opts.Quantity (1,1) string {mustBeMember(opts.Quantity, ...
        ["potential", "disturbance", "anomaly"])} = "potential"
    opts.MinDegree (1,1) double {mustBeInteger, mustBeNonnegative} = 2
end
P = numel(latDeg);
if numel(lonDeg) ~= P || numel(r) ~= P
    error('shx:synthesisPoints:sizeMismatch', ...
        'latDeg, lonDeg and r must have equal length.');
end
if ~isequal(size(C), size(S))
    error('shx:synthesisPoints:sizeMismatch', 'C and S must match.');
end
nmax = size(C, 1) - 1;
n0 = min(opts.MinDegree, nmax);
val = zeros(P, 1);
mvec = 0:nmax;
for k = 1:P
    Pnm = shx.legendreALF(nmax, deg2rad(latDeg(k)));  % radians in!
    cm = cosd(mvec * lonDeg(k));
    sm = sind(mvec * lonDeg(k));
    q = (R / r(k)) .^ (0:nmax)';
    switch opts.Quantity
        case "potential"
            f = GM / r(k) .* q;
        case "disturbance"
            f = GM / r(k)^2 .* ((0:nmax)' + 1) .* q;
        case "anomaly"
            f = GM / r(k)^2 .* ((0:nmax)' - 1) .* q;
    end
    inner = sum(Pnm .* (C .* cm + S .* sm), 2);  % per-degree sums
    val(k) = f(n0+1:end)' * inner(n0+1:end);
end
end
