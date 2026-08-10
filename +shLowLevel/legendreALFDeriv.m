function [P, D, D2] = legendreALFDeriv(nmax, latRad)
%LEGENDREALFDERIV Fully normalized ALFs and their latitude derivatives.
%
%   [P, D] = shLowLevel.legendreALFDeriv(NMAX, LATRAD) returns the 4pi-fully-
%   normalized associated Legendre functions P (via shLowLevel.legendreALF,
%   scaled Holmes-Featherstone, stable to degree 2190) and their exact
%   first derivatives with respect to LATITUDE phi,
%
%     dPbar_nm/dphi = 1/2 * ( c1 * sqrt((n-m)(n+m+1)) * Pbar_{n,m+1}
%                           - c2 * sqrt((n+m)(n-m+1)) * Pbar_{n,m-1} )
%     c1 = sqrt(2) if m == 0 else 1     (Pbar_{n,1} normalization)
%     c2 = sqrt(2) if m == 1 else 1     (Pbar_{n,0} normalization)
%
%   with absent neighbors treated as zero. The coefficient pattern was
%   calibrated per (n,m) by least squares against numerical
%   differentiation and then frozen (residuals ~1e-9 = finite-difference
%   floor; Python validation, v2.3). Needed for horizontal load
%   deformation (shLowLevel.shSynthesisDeformation) and any future gradient/
%   vector synthesis.
%
%   Inputs
%     nmax    (1,1) double            maximum degree
%     latRad  (1,nlat) double [rad]   geocentric latitudes
%   Outputs
%     P  (nmax+1, nmax+1, nlat) double   Pbar_nm(sin lat)
%     D  (nmax+1, nmax+1, nlat) double   dPbar_nm/dphi [1/rad]
%     D2 (nmax+1, nmax+1, nlat) double   d2Pbar_nm/dphi2 [1/rad^2]:
%        the SAME column stencil applied to D - exact, because the
%        first-derivative identity holds degree-wise as a function of
%        phi and differentiates term by term (Python-validated vs
%        numerical second differences, rel ~2e-7 = FD floor; v2.4)
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    nmax (1,1) double {mustBeInteger, mustBeNonnegative}
    latRad (1,:) double
end
P = shLowLevel.legendreALF(nmax, latRad);
n1 = nmax + 1;
nn = (0:nmax)';
mm = 0:nmax;
% coefficient triangles (zero where neighbors do not exist)
f1 = 0.5 * sqrt(max((nn - mm) .* (nn + mm + 1), 0));   % -> P_{n,m+1}
f2 = 0.5 * sqrt(max((nn + mm) .* (nn - mm + 1), 0));   % -> P_{n,m-1}
f1(:, 1) = f1(:, 1) * sqrt(2);                         % c1, m = 0
if n1 >= 2
    f2(:, 2) = f2(:, 2) * sqrt(2);                     % c2, m = 1
end
tl = tril(true(n1));
f1(~tl) = 0; f2(~tl) = 0;
D = zeros(size(P));
for k = 1:numel(latRad)
    D(:, :, k) = stencil(P(:, :, k), f1, f2, n1);
end
if nargout >= 3
    D2 = zeros(size(P));
    for k = 1:numel(latRad)
        D2(:, :, k) = stencil(D(:, :, k), f1, f2, n1);
    end
end
end

function B = stencil(A, f1, f2, n1)
%STENCIL The frozen first-derivative column stencil (validated v2.3).
Ap = [A(:, 2:end), zeros(n1, 1)];                      % (n, m+1)
Am = [zeros(n1, 1), A(:, 1:end-1)];                    % (n, m-1)
B = f1 .* Ap - f2 .* Am;
end
