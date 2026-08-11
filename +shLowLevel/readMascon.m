function mas = readMascon(filename)
%READMASCON Read a GRACE/GRACE-FO mascon netCDF (JPL / CSR / GSFC style).
%
%   MAS = shLowLevel.readMascon(FILENAME) reads gridded mascon solutions for
%   comparison against the SH chain (base MATLAB: ncread/ncinfo, no
%   toolboxes). Variable names are auto-detected among the layouts used
%   by JPL RL06 ('lwe_thickness'), CSR RL06 ('lwe_thickness'), and
%   generic 'lwe'/'water_thickness'; time units "days since YYYY-MM-DD"
%   are converted to decimal years.
%
%   VERIFIED LAYOUTS. The dimension order differs between products, and
%   the reader permutes to (lat, lon, time) whatever it finds:
%     JPL RL06.3Mv04 CRI   lwe_thickness(time, lat, lon), 0.5 deg,
%                          lat ascending -89.75..89.75, lon 0.25..359.75,
%                          time "days since 2002-01-01", units cm
%   That file also carries uncertainty, land_mask, scale_factor,
%   mascon_ID and GAD, which this reader ignores - read them with ncread
%   if you need them.
%
%   LATITUDES ARE AS STORED. Mascon grids are conventionally GEODETIC
%   while this toolbox synthesises on GEOCENTRIC latitudes, and the two
%   differ by up to 0.19 deg (about 21 km, or 0.4 of a 0.5-degree
%   mascon cell) around 45 deg latitude. Convert with
%   shLowLevel.geodetic2geocentric before comparing against a
%   synthesised field, or the comparison carries a systematic
%   mid-latitude shift that looks like a real discrepancy.
%
%   Inputs
%     filename  char/string  .nc path
%   Outputs (struct)
%     lat    (nlat,1) double [deg]   (as stored - typically GEODETIC;
%            convert with shLowLevel.geodetic2geocentric before synthesizing
%            SH fields on this grid!)
%     lon    (nlon,1) double [deg, 0..360)
%     epoch  (T,1) double decimal years
%     ewh    (nlat,nlon,T) double [cm or m AS STORED - check mas.units]
%     units, name  strings
%
%   The step-by-step mascon-comparison workflow (epoch matching, EWH
%   synthesis on this grid, stats) is in the workflow guide PDF, sec.
%   "Comparing against mascons".
%
%   Errors loudly (shLowLevel:readMascon:*) rather than guessing field names.
%
%   Claude (Fable 5), 2026-08-07.
%   Outputs
%     mas        struct: lat (I x 1), lon (J x 1), t (T x 1 decimal years), ewh (I x J x T) [m], units/meta from the netCDF
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    filename {mustBeTextScalar}
end
filename = char(filename);
assert(isfile(filename), 'shLowLevel:readMascon:fileNotFound', 'File not found: %s', filename);
info = ncinfo(filename);
vars = string({info.Variables.Name});
pick = @(cands) cands(find(ismember(cands, vars), 1));
vE = pick(["lwe_thickness", "lwe", "water_thickness", "ewh"]);
vLa = pick(["lat", "latitude"]);
vLo = pick(["lon", "longitude"]);
vT = pick(["time", "t"]);
assert(~isempty(vE) && ~isempty(vLa) && ~isempty(vLo) && ~isempty(vT), ...
    'shLowLevel:readMascon:unknownLayout', ...
    'Could not find EWH/lat/lon/time variables among: %s', strjoin(vars, ', '));
mas.lat = double(ncread(filename, vLa)); mas.lat = mas.lat(:);
mas.lon = double(ncread(filename, vLo)); mas.lon = mod(mas.lon(:), 360);
tRaw = double(ncread(filename, vT)); tRaw = tRaw(:);
tu = '';
try, tu = ncreadatt(filename, vT, 'units'); catch, end
tok = regexp(char(tu), 'days\s+since\s+(\d{4})-(\d{2})-(\d{2})', 'tokens', 'once');
if ~isempty(tok)
    t0 = datetime(str2double(tok{1}), str2double(tok{2}), str2double(tok{3}));
    tt = t0 + days(tRaw);
    y0 = dateshift(tt, 'start', 'year');
    mas.epoch = year(tt) + days(tt - y0) ./ days(y0 + calyears(1) - y0);
else
    mas.epoch = tRaw;                            % assume decimal years
end
E = double(ncread(filename, vE));
% normalize to nlat x nlon x T
sz = size(E); nlat = numel(mas.lat); nlon = numel(mas.lon);
if isequal(sz(1:2), [nlon nlat])
    E = permute(E, [2 1 3]);
elseif ~isequal(sz(1:2), [nlat nlon])
    error('shLowLevel:readMascon:badShape', 'EWH variable shaped %s.', mat2str(sz));
end
mas.ewh = E;
try, mas.units = string(ncreadatt(filename, vE, 'units')); catch, mas.units = ""; end
mas.name = string(filename);
end
