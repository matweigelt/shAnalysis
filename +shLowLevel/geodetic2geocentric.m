function latGC = geodetic2geocentric(latGD, opts)
%GEODETIC2GEOCENTRIC Convert geodetic to geocentric latitude.
%
%   LATGC = shLowLevel.geodetic2geocentric(LATGD) converts geodetic (ellipsoidal)
%   latitudes to geocentric latitudes via
%       tan(latGC) = (1 - f)^2 * tan(latGD)
%   All spherical-harmonic routines in this toolbox (synthesis, analysis,
%   Legendre functions) expect GEOCENTRIC latitude; user grids from maps,
%   GIS or mascon products are almost always GEODETIC. Forgetting this
%   conversion biases mid-latitude values by up to ~0.19 deg of latitude
%   (~21 km) - a classic silent error source.
%
%   Inputs
%     latGD  double array   geodetic latitudes [deg]
%     opts.Flattening (1,1) double = 1/298.257223563  (WGS84; override
%            for other reference ellipsoids - a DEFAULT, not physics
%            baked in)
%   Outputs
%     latGC  double array   geocentric latitudes [deg], poles/equator exact
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    latGD double
    opts.Flattening (1,1) double {mustBePositive} = 1/298.257223563
end
latGC = atand((1 - opts.Flattening)^2 .* tand(latGD));
pole = abs(abs(latGD) - 90) < 1e-12;
latGC(pole) = latGD(pole);                       % tand(90) guard
end
