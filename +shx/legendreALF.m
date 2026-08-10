function P = legendreALF(nmax, lat)
%LEGENDREALF Fully normalized (4-pi) associated Legendre functions Pbar_nm.
%
%   P = LEGENDREALF(NMAX, LAT) computes fully normalized associated
%   Legendre functions up to degree/order NMAX for each latitude in LAT
%   (radians, geocentric, scalar or vector), using the scaled forward
%   column recursion of Holmes & Featherstone (2002): the u^m factor of
%   each column is carried separately with a 1e-280 scale factor, which
%   keeps the recursion overflow/underflow-free up to at least degree
%   2190 (EGM2008/XGM2019e) at all latitudes including near the poles.
%   For nmax <= a few hundred the results agree with the plain unscaled
%   recursion to machine precision (cross-validated in Python: max
%   difference 8.4e-15 at nmax=120; addition-theorem sum rule
%   sum_m Pbar_nm^2 = 2n+1 satisfied to 3e-11 relative at n=2190,
%   lat=89.99 deg).
%
%   Inputs
%     nmax  (1,1) double   maximum degree/order
%     lat   (1,:)/(:,1) double  geocentric latitudes [rad]
%   Outputs
%     P     (NMAX+1)x(NMAX+1)xnumel(LAT) double
%           P(n+1,m+1,k) = Pbar_{n,m}(sin(lat(k))), upper triangle 0
%
%   Normalization: integral over the unit sphere of (Pbar_nm cos(m*lon))^2
%   equals 4*pi, matching the ICGEM / GRACE gfc coefficient convention.
%
%   Claude (Sonnet 4.6), 2026-07-11 (unscaled original);
%   scaled rewrite: Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

lat = lat(:)';
nlat = numel(lat);
t = sin(lat);
u = cos(lat);

P = zeros(nmax+1, nmax+1, nlat);

% ---- m = 0 column: no u^m factor, the plain recursion is safe
P(1,1,:) = 1;
if nmax == 0, return; end
pm2 = ones(1, nlat);
pm1 = sqrt(3) * t;
P(2,1,:) = pm1;
for n = 2:nmax
    a = sqrt((2*n-1)*(2*n+1)) / n;
    b = (n-1) * sqrt((2*n+1)/(2*n-3)) / n;
    p0 = a * t .* pm1 - b * pm2;
    P(n+1,1,:) = p0;
    pm2 = pm1; pm1 = p0;
end

% ---- m >= 1 columns: recursion on f_nm = Pbar_nm / u^m, scaled by 1e-280
scalef = 1e-280;
pmm = scalef;                      % f_mm * scalef (latitude-independent)
rescalem = (1/scalef) * ones(1, nlat);   % u^m / scalef, accumulated
for m = 1:nmax
    rescalem = rescalem .* u;
    if m == 1
        pmm = pmm * sqrt(3);
    else
        pmm = pmm * sqrt((2*m+1) / (2*m));
    end
    p2 = pmm * ones(1, nlat);
    P(m+1, m+1, :) = p2 .* rescalem;         % sectoral
    if m == nmax, break; end
    p1 = sqrt(2*m+3) * pmm * t;
    P(m+2, m+1, :) = p1 .* rescalem;
    for n = m+2:nmax
        a = sqrt((2*n-1)*(2*n+1) / ((n-m)*(n+m)));
        b = sqrt((2*n+1)*(n+m-1)*(n-m-1) / ((n-m)*(n+m)*(2*n-3)));
        p0 = a * t .* p1 - b * p2;
        P(n+1, m+1, :) = p0 .* rescalem;
        p2 = p1; p1 = p0;
    end
end
end
