function latGD = geocentric2geodetic(latGC, opts)
%GEOCENTRIC2GEODETIC Convert geocentric to geodetic latitude.
%
%   LATGD = shLowLevel.geocentric2geodetic(LATGC) inverts
%   shLowLevel.geodetic2geocentric:  tan(latGD) = tan(latGC) / (1 - f)^2.
%
%   Inputs
%     latGC  double array   geocentric latitudes [deg]
%     opts.Flattening (1,1) double = 1/298.257223563 (WGS84 default)
%   Outputs
%     latGD  double array   geodetic latitudes [deg]
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    latGC double
    opts.Flattening (1,1) double {mustBePositive} = 1/298.257223563
end
latGD = atand(tand(latGC) ./ (1 - opts.Flattening)^2);
pole = abs(abs(latGC) - 90) < 1e-12;
latGD(pole) = latGC(pole);
end
