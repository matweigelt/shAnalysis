function [G, lam, info] = slepianBasis(idx, region, opts)
%SLEPIANBASIS Spatiospectral concentration (Slepian) basis for a region.
%
%   [G, LAM, INFO] = shLowLevel.slepianBasis(IDX, REGION) solves the
%   Slepian concentration eigenproblem D g = lambda g for the region:
%   D(a,b) = (1/4pi) * integral over the region of Ybar_a Ybar_b, built
%   by latitude-ring quadrature with the toolbox's 4-pi-normalized
%   Legendre functions (legendreALF, radians). The columns of G are the
%   Slepian coefficient vectors in IDX ordering, sorted by decreasing
%   concentration LAM (all in [0, 1]); the leading sum(LAM) ~ Shannon
%   number N = P * A/(4pi) columns span the signal the region can
%   resolve. Validated against a Python/mpmath-free reference: Shannon
%   number exact to 4 digits, concentration identity g'Dg/g'g = lambda
%   to 1e-10, polar-cap order decoupling to 1e-17 (tools/audit lineage).
%
%   Inputs
%     idx     (1 x 1) struct  index from shLowLevel.shIndex; use
%             MinDegree = 0 for a complete basis (P = (Lmax+1)^2)
%     region  (nLat x nLon logical | Q x 2 double) either a mask on the
%             grid implied by LatDeg/LonDeg, or a closed polygon
%             [lat lon] in degrees (lon in [0, 360))
%
%   Options
%     LatDeg ([])    (nLat x 1) mask-grid latitudes [deg]; [] uses the
%                    GridStep graticule
%     LonDeg ([])    (nLon x 1) mask-grid longitudes [deg]
%     GridStep (0.5) (1 x 1) quadrature step [deg] for polygon regions
%     K ([])         (1 x 1) return only the leading K columns; []
%                    returns all P
%     Quiet (true)   (1 x 1) suppress progress output
%
%   Outputs
%     G    (P x K) Slepian coefficient vectors, orthonormal columns
%     lam  (K x 1) concentration ratios, descending, in [0, 1]
%     info (1 x 1) struct with fields
%       shannon      (1 x 1) sum of ALL eigenvalues = P * areaFraction
%       areaFraction (1 x 1) region area / sphere area
%       P            (1 x 1) basis size
%
%   Example
%     idx = shLowLevel.shIndex(30, MinDegree = 0);
%     cap = [60 0; 60 90; 60 180; 60 270; 60 359.99];  % crude N cap
%     [G, lam, info] = shLowLevel.slepianBasis(idx, cap);
%     fprintf('Shannon %.1f, leading lambda %.4f\n', info.shannon, lam(1))
%
%   Error identifiers
%     shLowLevel:slepianBasis:badRegion   region neither mask nor Q x 2
%     shLowLevel:slepianBasis:emptyRegion mask selects no grid point
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.9.0).
arguments
    idx (1,1) struct
    region
    opts.LatDeg (:,1) double = []
    opts.LonDeg (:,1) double = []
    opts.GridStep (1,1) double {mustBePositive} = 0.5
    opts.K double {mustBeScalarOrEmpty} = []
    opts.Quiet (1,1) logical = true
end
% ---- region -> mask + grid
if islogical(region)
    mk = region;
    lat = opts.LatDeg; lon = opts.LonDeg;
    if isempty(lat) || isempty(lon)
        st = opts.GridStep;
        lat = (-90+st/2 : st : 90-st/2)'; lon = (st/2 : st : 360-st/2)';
    end
    if ~isequal(size(mk), [numel(lat), numel(lon)])
        error('shLowLevel:slepianBasis:badRegion', ...
            'mask is %d x %d but the grid is %d x %d.', ...
            size(mk,1), size(mk,2), numel(lat), numel(lon));
    end
elseif isnumeric(region) && size(region, 2) == 2
    st = opts.GridStep;
    lat = (-90+st/2 : st : 90-st/2)'; lon = (st/2 : st : 360-st/2)';
    [LO, LA] = meshgrid(lon, lat);
    lonW = mod(LO + 180, 360) - 180;
    pw = mod(region(:,2) + 180, 360) - 180;
    mk = inpolygon(lonW, LA, pw, region(:,1));
else
    error('shLowLevel:slepianBasis:badRegion', ...
        'region must be a logical mask or a Q x 2 [lat lon] polygon.');
end
if ~any(mk(:))
    error('shLowLevel:slepianBasis:emptyRegion', ...
        'the region selects no grid point at this GridStep.');
end
% ---- D by latitude-ring quadrature
P = idx.P;
dphi = deg2rad(abs(lat(2) - lat(1))); dlam = deg2rad(abs(lon(2) - lon(1)));
D = zeros(P, P);
areaFrac = 0;
mlon = deg2rad(lon(:)') .* idx.m(:);              % P x nLon phase table
isC = idx.cs(:) == 0;
for i = 1:numel(lat)
    cols = mk(i, :);
    if ~any(cols), continue; end
    Pl = shLowLevel.legendreALF(idx.Lmax, deg2rad(lat(i)));
    pv = Pl(sub2ind(size(Pl), idx.n(:)+1, idx.m(:)+1));   % P x 1
    T = zeros(P, nnz(cols));
    T(isC, :)  = cos(mlon(isC, cols));
    T(~isC, :) = sin(mlon(~isC, cols));
    Y = (pv .* T)';                                % nIn x P
    w = cosd(lat(i)) * dphi * dlam;
    D = D + w * (Y' * Y);
    areaFrac = areaFrac + w * nnz(cols);
end
D = D / (4*pi);
areaFrac = areaFrac / (4*pi);
% ---- eigendecomposition, descending
D = (D + D') / 2;
[V, L] = eig(D, 'vector');
[lamAll, si] = sort(L, 'descend');
V = V(:, si);
lamAll = min(max(lamAll, 0), 1);
K = opts.K; if isempty(K), K = P; end
K = min(K, P);
G = V(:, 1:K); lam = lamAll(1:K);
info = struct('shannon', sum(lamAll), 'areaFraction', areaFrac, 'P', P);
if ~opts.Quiet
    fprintf('slepianBasis: P=%d, Shannon %.2f (P*A/4pi = %.2f)\n', ...
        P, info.shannon, P*areaFrac);
end
end
