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
%     CSR RL0603           lwe_thickness(lon, lat, time), 0.25 deg,
%                          time attribute named 'Units' (capital U) -
%                          attribute names are matched case-insensitively
%                          since netCDF producers differ here
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
%   Errors loudly (shLowLevel:readMascon:*) rather than guessing: unknown
%   variable layouts, missing or unrecognized time units, and ambiguous
%   grid orientation all raise identified errors instead of returning a
%   plausible-looking wrong result.
%
%   Claude (Fable 5), 2026-08-07; CSR support + loud time units: 2026-08-12 (v3.8.6).
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
% attribute names are NOT standardized in case: JPL writes 'units', CSR
% RL0603 writes 'Units'. ncreadatt is case-sensitive, so look the name up
% case-insensitively in the ncinfo metadata instead of guessing.
tu = localAtt(info, vT, 'units');
assert(~isempty(tu), 'shLowLevel:readMascon:noTimeUnits', ...
    ['Time variable ''%s'' carries no units attribute - cannot convert ' ...
     'to decimal years. Read the file with ncread and convert manually.'], vT);
tok = regexp(char(tu), 'days\s+since\s+(\d{4})-(\d{2})-(\d{2})', 'tokens', 'once');
if isempty(tok)
    error('shLowLevel:readMascon:unknownTimeUnits', ...
        ['Time units "%s" not recognized (expected "days since ' ...
         'YYYY-MM-DD..."). Refusing to guess: returning the raw values ' ...
         'as decimal years once shipped raw days silently.'], char(tu));
end
t0 = datetime(str2double(tok{1}), str2double(tok{2}), str2double(tok{3}));
tt = t0 + days(tRaw);
y0 = dateshift(tt, 'start', 'year');
mas.epoch = year(tt) + days(tt - y0) ./ days(y0 + calyears(1) - y0);
E = double(ncread(filename, vE));
% normalize to (lat, lon, time) by DIMENSION NAME, not by size: a size
% heuristic cannot disambiguate a square grid (nlat == nlon) and would
% silently transpose it. Names are authoritative in every real product.
iv = strcmp({info.Variables.Name}, vE);
dn = lower(string({info.Variables(iv).Dimensions.Name}));
pLat = find(dn == lower(string(vLa)) | dn == "lat" | dn == "latitude", 1);
pLon = find(dn == lower(string(vLo)) | dn == "lon" | dn == "longitude", 1);
pT   = find(dn == lower(string(vT))  | dn == "time" | dn == "t", 1);
if numel([pLat pLon pT]) == 3
    E = permute(E, [pLat pLon pT]);
else
    % fall back to the size heuristic only where it is unambiguous
    sz = size(E); nlat = numel(mas.lat); nlon = numel(mas.lon);
    assert(nlat ~= nlon, 'shLowLevel:readMascon:ambiguousShape', ...
        ['Cannot orient a square %dx%d grid without named dimensions ' ...
         '(EWH dims: %s).'], nlat, nlon, strjoin(cellstr(dn), ','));
    if isequal(sz(1:2), [nlon nlat])
        E = permute(E, [2 1 3]);
    elseif ~isequal(sz(1:2), [nlat nlon])
        error('shLowLevel:readMascon:badShape', 'EWH variable shaped %s.', mat2str(sz));
    end
end
mas.ewh = E;
mas.units = string(localAtt(info, vE, 'units'));
mas.name = string(filename);
end

function v = localAtt(info, vName, aName)
%LOCALATT Case-insensitive variable-attribute lookup from ncinfo metadata.
v = '';
iv = strcmp({info.Variables.Name}, vName);
A = info.Variables(iv).Attributes;
if isempty(A), return; end
ia = find(strcmpi({A.Name}, aName), 1);
if ~isempty(ia), v = A(ia).Value; end
end
