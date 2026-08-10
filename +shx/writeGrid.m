function writeGrid(filename, grid, latDeg, lonDeg, opts)
%WRITEGRID Export gridded fields to CF-style netCDF (base MATLAB).
%
%   shx.writeGrid('ewh.nc', GRID, LAT, LON, Name="lwe_thickness",
%       Units="m", Epoch=2010.29) writes a 2-D field or a 3-D stack
%   (nlat x nlon x T with Epochs=) using nccreate/ncwrite - directly
%   readable back by shx.readMascon (roundtrip-tested), by GMT, panoply,
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
    'shx:writeGrid:badSize', 'GRID must be nlat x nlon (x T).');
if T > 1
    assert(numel(opts.Epochs) == T, 'shx:writeGrid:badEpochs', ...
        '3-D stacks need Epochs (T = %d).', T);
end
vn = char(opts.Name);
nccreate(filename, 'lat', 'Dimensions', {'lat', nlat});
nccreate(filename, 'lon', 'Dimensions', {'lon', nlon});
ncwrite(filename, 'lat', latDeg);
ncwrite(filename, 'lon', mod(lonDeg, 360));
ncwriteatt(filename, 'lat', 'units', 'degrees_north');
ncwriteatt(filename, 'lat', 'lat_type', char(opts.LatType));
ncwriteatt(filename, 'lon', 'units', 'degrees_east');
if T > 1
    nccreate(filename, 'time', 'Dimensions', {'time', T});
    t0 = datetime(2002, 1, 1);
    yy = floor(opts.Epochs);
    frac = opts.Epochs - yy;
    dts = datetime(yy, 1, 1) + years(0) + ...
        days(frac .* days(datetime(yy + 1, 1, 1) - datetime(yy, 1, 1)));
    ncwrite(filename, 'time', days(dts - t0));
    ncwriteatt(filename, 'time', 'units', 'days since 2002-01-01');
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
if strlength(opts.Description) > 0
    ncwriteatt(filename, '/', 'description', char(opts.Description));
end
ncwriteatt(filename, '/', 'source', ...
    'shAnalysis v2.4 - Claude (Fable 5), 2026-08-07');
end
