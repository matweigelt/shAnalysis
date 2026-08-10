function writeGrid(filename, grid, latDeg, lonDeg, opts)
%WRITEGRID Export gridded fields to CF-style netCDF (base MATLAB).
%
%   shLowLevel.writeGrid('ewh.nc', GRID, LAT, LON, Name ("field")="lwe_thickness",
%       Units ("")="m", Epoch=2010.29) writes a 2-D field or a 3-D stack
%   (nlat x nlon x T with Epochs ([])=) using nccreate/ncwrite - directly
%   readable back by shLowLevel.readMascon (roundtrip-tested), by GMT, panoply,
%   xarray, etc. Time is stored as "days since 2002-01-01" like the
%   mascon products.
%
%   Inputs
%     filename  char/string  .nc target (existing file is replaced;
%               missing parent folders are created, v2.4.1)
%     grid      (nlat,nlon) or (nlat,nlon,T) double
%     latDeg    (nlat,1) double [deg]  (document your latitude type -
%               geocentric unless you converted; attribute 'lat_type')
%     lonDeg    (nlon,1) double [deg, 0..360)
%   Options
%     Name ("field"), Units (""), Epoch/Epochs (decimal years, for 3-D),
%     LatType ("geocentric"), Description ("")
%     Sidecar (true)  write <file>.provenance.json alongside the output
%
%   Claude (Fable 5), 2026-08-07 (v2.4).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    filename {mustBeTextScalar}
    grid double
    latDeg (:,1) double
    lonDeg (:,1) double
    opts.Name (1,1) string = "field"
    opts.Units (1,1) string = ""
    opts.Epochs (:,1) double = []
    opts.LatType (1,1) string = "geocentric"
    opts.Description (1,1) string = ""
    opts.Sidecar (1,1) logical = true
end
filename = char(filename);
parent = fileparts(filename);
if ~isempty(parent) && ~isfolder(parent)
    mkdir(parent);                       % create target folder on demand
end
if isfile(filename), delete(filename); end
nlat = numel(latDeg); nlon = numel(lonDeg);
T = size(grid, 3);
assert(size(grid, 1) == nlat && size(grid, 2) == nlon, ...
    'shLowLevel:writeGrid:badSize', 'GRID must be nlat x nlon (x T).');
if T > 1
    assert(numel(opts.Epochs) == T, 'shLowLevel:writeGrid:badEpochs', ...
        '3-D stacks need Epochs (T = %d).', T);
end
vn = char(opts.Name);
nccreate(filename, 'lat', 'Dimensions', {'lat', nlat});
nccreate(filename, 'lon', 'Dimensions', {'lon', nlon});
ncwrite(filename, 'lat', latDeg);
ncwrite(filename, 'lon', mod(lonDeg, 360));
ncwriteatt(filename, 'lat', 'units', 'degrees_north');
ncwriteatt(filename, 'lat', 'standard_name', 'latitude');
ncwriteatt(filename, 'lat', 'axis', 'Y');
ncwriteatt(filename, 'lat', 'lat_type', char(opts.LatType));
ncwriteatt(filename, 'lon', 'units', 'degrees_east');
ncwriteatt(filename, 'lon', 'standard_name', 'longitude');
ncwriteatt(filename, 'lon', 'axis', 'X');
if T > 1
    nccreate(filename, 'time', 'Dimensions', {'time', T});
    t0 = datetime(2002, 1, 1);
    yy = floor(opts.Epochs);
    frac = opts.Epochs - yy;
    dts = datetime(yy, 1, 1) + years(0) + ...
        days(frac .* days(datetime(yy + 1, 1, 1) - datetime(yy, 1, 1)));
    ncwrite(filename, 'time', days(dts - t0));
    ncwriteatt(filename, 'time', 'units', 'days since 2002-01-01');
    ncwriteatt(filename, 'time', 'standard_name', 'time');
    ncwriteatt(filename, 'time', 'calendar', 'standard');
    ncwriteatt(filename, 'time', 'axis', 'T');
    nccreate(filename, vn, 'Dimensions', ...
        {'lon', nlon, 'lat', nlat, 'time', T});
    ncwrite(filename, vn, permute(grid, [2 1 3]));
else
    nccreate(filename, vn, 'Dimensions', {'lon', nlon, 'lat', nlat});
    ncwrite(filename, vn, grid.');
end
if strlength(opts.Units) > 0
    ncwriteatt(filename, vn, 'units', char(opts.Units));
end
ncwriteatt(filename, vn, 'long_name', ...
    ternary(strlength(opts.Description) > 0, char(opts.Description), vn));
ncwriteatt(filename, vn, 'coordinates', 'lat lon');
if strlength(opts.Description) > 0
    ncwriteatt(filename, '/', 'description', char(opts.Description));
end
v = shLowLevel.version();
ncwriteatt(filename, '/', 'Conventions', 'CF-1.8');
ncwriteatt(filename, '/', 'title', vn);
ncwriteatt(filename, '/', 'source', sprintf('%s v%s (%s)', ...
    v.Name, v.Version, v.Date));
ncwriteatt(filename, '/', 'history', sprintf('%s: written by %s v%s', ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
    v.Name, v.Version));
if opts.Sidecar
    writeSidecar(string(filename), struct('name', vn, ...
        'nlat', nlat, 'nlon', nlon, 'nEpochs', T, ...
        'units', char(opts.Units), 'latType', char(opts.LatType)));
end
end

function writeSidecar(target, extra)
% provenance JSON sidecar: <target>.provenance.json (v2.7.0)
v = shLowLevel.version();
p = struct('tool', sprintf('%s %s (%s)', v.Name, v.Version, v.Date), ...
    'provenance', v.Provenance, ...
    'created', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
    'matlab', sprintf('%s / %s', version, computer), ...
    'file', char(target));
fn = fieldnames(extra);
for k = 1:numel(fn), p.(fn{k}) = extra.(fn{k}); end
fid = fopen(char(string(target) + ".provenance.json"), 'w');
fprintf(fid, '%s\n', jsonencode(p, 'PrettyPrint', true));
fclose(fid);
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
