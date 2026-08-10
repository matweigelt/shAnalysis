function [grid, lat, lon, P] = shSynthesis(C, S, GM, R, latVec, lonVec, varargin)
%SHSYNTHESIS Spherical harmonic synthesis of a global/regional grid.
%
%   [GRID, LAT, LON] = SHSYNTHESIS(C, S, GM, R, LATVEC, LONVEC) synthesizes
%   the disturbing-potential-derived quantity (default: geoid height) on the
%   grid LATVEC x LONVEC (degrees) from fully normalized (ICGEM) Stokes
%   coefficients C, S at the reference sphere r = R.
%
%   [GRID, LAT, LON, P] = SHSYNTHESIS(...) also returns the precomputed
%   Legendre function array P (as from LEGENDREALF). For batch processing
%   a time series on the SAME grid and NMAX (e.g. monthly GRACE-FO
%   solutions), pass it back in on subsequent calls via 'P', P to skip
%   the recursion entirely:
%
%       [grid1, lat, lon, P] = shSynthesis(C1, S1, GM, R, latVec, lonVec);
%       for k = 2:nMonths
%           gridK = shSynthesis(Ck, Sk, GM, R, latVec, lonVec, 'P', P);
%       end
%
%   Name/value options:
%     'P'  precomputed Legendre array from a prior call/LEGENDREALF, for
%          the identical latVec and nmax -- skips recomputing the
%          recursion. No validation is done that latVec/nmax actually
%          match; passing a mismatched P silently gives wrong results.
%     'quantity'  one of:
%         'geoid'              N = R * sum Pbar_nm (Cnm cos m lon + Snm sin m lon)   [m]
%         'gravity_anomaly'    dg = GM/R^2 * sum (n-1) Pbar_nm (...)                 [m/s^2]
%         'gravity_disturbance' dg = GM/R^2 * sum (n+1) Pbar_nm (...)                [m/s^2]
%         'potential'          T  = GM/R * sum Pbar_nm (...)                        [m^2/s^2]
%         'ewh'                equivalent water height using love numbers 'kn'      [m]
%       default 'geoid'
%     'nmin','nmax'  degree range to include (default 0, size(C,1)-1)
%     'kn'           load Love number vector k_n (index 1 = degree 0), required
%                    for 'ewh'; NOT supplied by default -- see note below.
%     'rho_ave'      average Earth density, default 5517 kg/m^3 (used for 'ewh')
%     'rho_water'    water density, default 1000 kg/m^3 (used for 'ewh')
%
%   NOTE on 'ewh': accurate mass-change conversion requires the degree-
%   dependent load Love numbers k_n (e.g. from a PREM-based table, Wahr et
%   al. 1998 / Han & Wahr). This function does NOT hardcode a Love number
%   table to avoid silently using unvalidated constants -- pass 'kn'
%   explicitly. Without it, EWH synthesis will error.
%
%   Output GRID is numel(LATVEC) x numel(LONVEC).
%
%   'method'  'auto' (default) / 'direct' / 'fft': FFT-along-longitude
%             synthesis for uniform full-circle longitude grids (5-20x
%             faster at 1 degree); 'auto' picks it whenever applicable.
%             Symmetric latitude grids additionally reuse the Legendre
%             recursion for +/-lat pairs via the parity relation.
%
%   Claude (Sonnet 4.6), 2026-07-11; merged into +shx and FFT/parity
%   fast paths added: Claude (Fable 5), 2026-08-07.
%   Outputs
%     grid       (nlat x nlon) double   synthesized field in the requested quantity units
%     lat        (1 x nlat) double   geocentric latitudes used
%     lon        (1 x nlon) double   longitudes used
%     P          struct   reusable Legendre/cache handle for repeat calls
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

p = inputParser;
addParameter(p, 'quantity', 'geoid');
addParameter(p, 'nmin', 0);
addParameter(p, 'nmax', size(C,1)-1);
addParameter(p, 'kn', []);
addParameter(p, 'rho_ave', 5517);
addParameter(p, 'rho_water', 1000);
addParameter(p, 'P', []);
addParameter(p, 'method', 'auto');   % 'auto' | 'direct' | 'fft' (v2.1)
addParameter(p, 'MaxMemGB', 4);      % Legendre-block memory budget (v2.2)
addParameter(p, 'hn', []);           % load Love numbers h'_n (v2.3)
addParameter(p, 'Height', 0);        % evaluation height above R [m] (v2.3)
parse(p, varargin{:});
quantity = lower(p.Results.quantity);
nmin = p.Results.nmin;
nmax = p.Results.nmax;

if size(C,1) ~= size(C,2) || size(S,1) ~= size(S,2) || size(C,1) ~= size(S,1)
    error('shSynthesis:badInput', 'C and S must be square matrices of equal size.');
end
if nmax > size(C,1)-1
    error('shSynthesis:badInput', ...
        'Requested nmax (%d) exceeds available coefficients (max %d).', nmax, size(C,1)-1);
end

lat = latVec(:)';
lon = lonVec(:)';
nlat = numel(lat);
nlon = numel(lon);
latRad = deg2rad(lat);
lonRad = deg2rad(lon);

if isempty(p.Results.P)
    % symmetric-latitude parity trick: Pbar_nm(-lat) = (-1)^(n+m) Pbar_nm(lat)
    % -> compute the recursion only for unique |lat| and mirror (validated:
    % exact identity), halving the Legendre cost for symmetric grids.
    banded = [];                                 % Legendre in memory bands
else
    P = p.Results.P;
    if size(P,1) < nmax+1 || size(P,3) ~= nlat
        error('shSynthesis:badP', ...
            'Supplied P is incompatible with nmax/latVec (size %s).', mat2str(size(P)));
    end
    banded = false;
end

kernel = shx.kernelFactors(quantity, nmax, GM, R, kn = p.Results.kn, ...
    hn = p.Results.hn, Height = p.Results.Height, ...
    rho_ave = p.Results.rho_ave, rho_water = p.Results.rho_water);

% degree-factor-scaled coefficient copies, nmin truncation via zero rows
KC = kernel .* C(1:nmax+1, 1:nmax+1);
KS = kernel .* S(1:nmax+1, 1:nmax+1);
if nmin > 0
    KC(1:nmin, :) = 0;
    KS(1:nmin, :) = 0;
end

% per-order sums over degree: Am(m+1,k) = sum_n f_n C_nm Pbar_nm(lat_k)
% Memory model (v2.2): the (nmax+1)^2 x nlat Legendre stack is the memory
% driver - 7 GB at nmax=2190 / 181 latitudes. When it exceeds MaxMemGB,
% latitudes are processed in bands: compute the Legendre block, fold it
% into Am/Bm, discard. Only Am/Bm (2 x (nmax+1) x nlat) persist - the
% high-degree capability becomes real instead of nominal. Chunked ==
% monolithic exactly (Python-validated, 0.0); the parity trick still
% applies within each unique-|lat| chunk.
Am = zeros(nmax+1, nlat);
Bm = zeros(nmax+1, nlat);
if isempty(banded)
    needBytes = (nmax+1)^2 * nlat * 8;
    banded = needBytes > p.Results.MaxMemGB * 2^30;
end
if ~banded && isempty(p.Results.P)
    [uab, ~, iu] = unique(abs(latRad));
    if numel(uab) < nlat
        Pu = shx.legendreALF(nmax, uab);
        P = Pu(:, :, iu);
        parity = (-1).^((0:nmax)' + (0:nmax));      % (n+1) x (m+1)
        P(:, :, latRad < 0) = P(:, :, latRad < 0) .* parity;
    else
        P = shx.legendreALF(nmax, latRad);
    end
end
if banded
    P = [];   % Legendre blocks are streamed and discarded in banded mode;
              % request the 4th output only with a fitting memory budget
    [uab, ~, iu] = unique(abs(latRad));
    parity = (-1).^((0:nmax)' + (0:nmax));
    chunk = max(1, floor(p.Results.MaxMemGB * 2^30 / ((nmax+1)^2 * 8)));
    for c0 = 1:chunk:numel(uab)
        cc = c0:min(c0+chunk-1, numel(uab));
        Pu = shx.legendreALF(nmax, uab(cc));
        for j = find(ismember(iu', cc))
            Pk = Pu(:, :, iu(j) - c0 + 1);
            if latRad(j) < 0, Pk = Pk .* parity; end
            Am(:, j) = sum(KC .* Pk, 1)';
            Bm(:, j) = sum(KS .* Pk, 1)';
        end
    end
else
    for k = 1:nlat
        Pk = P(1:nmax+1, 1:nmax+1, k);
        Am(:, k) = sum(KC .* Pk, 1)';
        Bm(:, k) = sum(KS .* Pk, 1)';
    end
end

% method resolution: FFT along longitude for uniform full-circle grids
% (exact accumulation of aliased orders into mod(m, nlon) bins; validated
% against the direct sum to 3e-13 in Python), else direct trig product.
method = lower(p.Results.method);
if strcmp(method, 'auto')
    useFFT = false;
    if nlon >= 2
        dl = diff(lonRad);
        if all(abs(dl - dl(1)) < 1e-12) && dl(1) > 0 ...
                && abs(nlon*dl(1) - 2*pi) < 1e-9
            useFFT = true;
        end
    end
elseif strcmp(method, 'fft')
    dl = diff(lonRad);
    if nlon < 2 || any(abs(dl - dl(1)) >= 1e-12) || dl(1) <= 0 ...
            || abs(nlon*dl(1) - 2*pi) >= 1e-9
        error('shSynthesis:badMethod', ...
            'method=''fft'' requires a uniform, full-circle longitude vector.');
    end
    useFFT = true;
else
    useFFT = false;
end

if useFFT
    ph  = exp(1i * (0:nmax)' * lonRad(1));
    Csp = (Am - 1i*Bm) .* ph;                % (nmax+1) x nlat
    Cpad = zeros(nlon, nlat);
    bins = mod(0:nmax, nlon) + 1;
    for mm = 0:nmax
        Cpad(bins(mm+1), :) = Cpad(bins(mm+1), :) + Csp(mm+1, :);
    end
    grid = real(ifft(Cpad, [], 1)).' * nlon;
else
    m = (0:nmax)';
    cosMLon = cos(m * lonRad);               % (nmax+1) x nlon
    sinMLon = sin(m * lonRad);
    grid = Am.' * cosMLon + Bm.' * sinMLon;
end

end
