function [mk2, dmin] = erodeMask(mk, lat, lon, bufferKm)
%ERODEMASK Erode a region mask by a great-circle boundary buffer.
%   MK2 = shLowLevel.erodeMask(MK, LAT, LON, KM) keeps only mask pixels whose
%   great-circle distance to the nearest OUTSIDE pixel exceeds
%   BUFFERKM. Distances are computed against the boundary pixels of
%   the complement only (coast), so the cost is nCoast x nGrid.
%
%   Inputs
%     mk       (nLat x nLon) logical  region mask (true = keep)
%     lat, lon (nLat x 1), (nLon x 1) grid vectors [deg]
%     bufferKm (1 x 1) double  buffer distance [km]; 0 returns mk
%
%   Outputs
%     mk2   (nLat x nLon) logical  eroded mask
%     dmin  (nLat x nLon) double   distance to the nearest outside
%           coast pixel [km] (Inf where no outside pixel exists)
%
%   Example
%     oc = abs(LA) <= 66 & ~landBoxes;              % crude ocean mask
%     oc300 = shLowLevel.erodeMask(oc, lat, lon, 300);
%     nnz(oc) - nnz(oc300)                          % coast pixels shed
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-13 (v3.15.0).
Re = 6371;
if bufferKm <= 0, mk2 = mk; dmin = []; return, end
out = ~mk;
% coast = outside pixels with at least one inside neighbour (wrap in lon)
nb = circshift(mk, [0, 1]) | circshift(mk, [0, -1]) | ...
     [mk(2:end, :); false(1, size(mk, 2))] | ...
     [false(1, size(mk, 2)); mk(1:end-1, :)];
coast = out & nb;
[LO, LA] = meshgrid(deg2rad(lon), deg2rad(lat));
ci = find(coast);
dmin = inf(size(mk));
sinLA = sin(LA); cosLA = cos(LA);
for k = 1:numel(ci)
    i = ci(k);
    cp = sinLA(i) * sinLA + cosLA(i) * cosLA .* cos(LO - LO(i));
    d = Re * acos(min(1, max(-1, cp)));
    dmin = min(dmin, d);
end
mk2 = mk & dmin > bufferKm;
end
