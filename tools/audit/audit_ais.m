function audit_ais
%AUDIT_AIS GravIS Antarctica chain reproduction on the audited toolbox.
%   Mirrors audit_gravis.m (Greenland) on the SAME cached, corrected,
%   GIA-reduced trend grids (span <= 2023.099, chapter regime). The
%   reference is my independent fit to the 25 GravIS AIS basin series
%   (301..325, COST-G): total -146.9 Gt/yr over <=2023.10 (n=217);
%   full-record -135.4 (2002.29..2025.96, n=252).
%   Claude (Fable 5), 2026-08-12.

aux = 'E:\DATAPOOL\GravityField\GravIS';
S  = load(fullfile(aux, 'audit_chain_cache_g445.mat'));
Su = load(fullfile(aux, 'audit_chain_cache_unf.mat'));
lat = S.lat(:); lon = S.lon(:); Re = 6371e3;

% ---- AIS masks from the GravIS basin polygons themselves
J = jsondecode(fileread(fullfile(aux, 'basins_AIS.json')));
[AIS, cLon] = polys2mask(J, lat, lon);
% regional split by basin centroid longitude (lonW):
%   Peninsula: -75..-55 | WAIS: -140..-75 | EAIS: rest
PEN = false(size(AIS)); WAI = PEN; EAI = PEN;
for k = 1:numel(J.features)
    mk = polys2mask(subset(J, k), lat, lon);
    c = cLon(k);
    if c >= -75 && c <= -55, PEN = PEN | mk;
    elseif c >= -140 && c < -75, WAI = WAI | mk;
    else, EAI = EAI | mk;
    end
end
fprintf('AIS mask: %.2f%% of grid, area %.2fe6 km2 | PEN/WAIS/EAIS pixels %d/%d/%d\n', ...
    100*nnz(AIS)/numel(AIS), area_km2(lat, lon, AIS, Re)/1e6, nnz(PEN), nnz(WAI), nnz(EAI));

fprintf(['reference (my fits to GravIS AIS 301..325): total -146.9 (<=2023.10, n=217), ' ...
    '-135.4 full\n']);
fprintf('noise: sigTrend=%.3e m/yr (from cache; trend-matched)\n\n', S.sigTrend);

runs = struct('name',{},'Gt',{},'it',{},'stop',{});
    function do(name, Sx, mask, filt, noise, maxit)
        [m, info] = shLowLevel.leakageCorrect(Sx.trendGrid, lat, lon, ...
            Mask = mask, Filter = filt, NoiseLevel = noise, ...
            MaxIter = maxit, Quiet = true);
        Gt = mass_gt(m, lat, lon, mask, Re);
        runs(end+1) = struct('name',name,'Gt',Gt,'it',info.iterations, ...
            'stop',string(info.stoppedBy));
        fprintf('%-34s: %+8.1f Gt/yr  (%3d it, %s)\n', name, Gt, ...
            info.iterations, info.stoppedBy);
    end

% headline + the Greenland-lesson variants
do("AIS gauss445 discrepancy",        S,  AIS, "gauss445", S.sigTrend, 400);
do("AIS unfiltered none 300",         Su, AIS, "none",     S.sigTrend, 300);
do("AIS gauss445 fixed300",           S,  AIS, "gauss445", 0,          300);
do("AIS wrong noise (trend oceanRMS)",S,  AIS, "gauss445", S.trendOceanRMS, 400);
% GIA sensitivity: add the GIA trend field back (= uncorrected chain)
kn = readmatrix(fullfile('C:\Users\matth\Documents\MATLAB\shAnalysis', ...
    'tests','test_data','loadLoveNumbers_Gegout97.txt'), FileType='text', NumHeaderLines=2);
Gi = shLowLevel.readSHM(fullfile(aux,'GRAVIS-2B_COSTG_0200_GIA_ICE-6G_D_VM5a_0001.gz'), Nmax=90);
r445 = 4*111.195; w = shLowLevel.shGaussianWeights(90, r445); w = w(:);
giaTrend = shLowLevel.shSynthesis(Gi.C(1:91,1:91).*w, Gi.S(1:91,1:91).*w, ...
    3.986004415e14, 6378136.3, lat, lon, 'quantity','ewh','kn',kn,'nmin',0);
Sg = S; Sg.trendGrid = S.trendGrid + giaTrend;
do("AIS withOUT GIA correction",      Sg, AIS, "gauss445", S.sigTrend, 400);
fprintf('   (GIA field alone over AIS: %+7.1f Gt/yr apparent)\n', ...
    mass_gt(giaTrend, lat, lon, AIS, Re));
% regional split, headline settings
do("  Peninsula",                     S,  PEN, "gauss445", S.sigTrend, 400);
do("  WAIS",                          S,  WAI, "gauss445", S.sigTrend, 400);
do("  EAIS",                          S,  EAI, "gauss445", S.sigTrend, 400);

save(fullfile(aux, 'audit_ais_runs.mat'), 'runs');
end

% ------------------------------------------------------------------ local
function [mask, cLon] = polys2mask(J, lat, lon)
% union mask of all features; handles Polygon/MultiPolygon + dateline by
% recentring each part on its own mean longitude
[LO, LA] = meshgrid(lon, lat);
lonW = mod(LO + 180, 360) - 180;
mask = false(size(LA)); cLon = zeros(numel(J.features), 1);
for k = 1:numel(J.features)
    g = J.features(k).geometry;
    parts = {};
    if strcmp(g.type, 'Polygon')
        parts = ringcells(g.coordinates);
    else                                          % MultiPolygon
        if iscell(g.coordinates)
            for q = 1:numel(g.coordinates)
                parts = [parts, ringcells(g.coordinates{q})]; %#ok<AGROW>
            end
        else
            for q = 1:size(g.coordinates, 1)
                parts = [parts, ringcells(squeeze(g.coordinates(q,:,:,:)))]; %#ok<AGROW>
            end
        end
    end
    ml = false(size(LA)); allx = [];
    for q = 1:numel(parts)
        P = parts{q};                             % K x 2 [lon lat]
        c = atan2d(mean(sind(P(:,1))), mean(cosd(P(:,1))));
        px = mod(P(:,1) - c + 180, 360) - 180;
        gx = mod(lonW - c + 180, 360) - 180;
        ml = ml | inpolygon(gx, LA, px, P(:,2));
        allx = [allx; P(:,1)]; %#ok<AGROW>
    end
    mask = mask | ml;
    cLon(k) = atan2d(mean(sind(allx)), mean(cosd(allx)));
end
end

function parts = ringcells(coord)
% outer ring only (holes irrelevant for drainage basins)
if iscell(coord), C = coord{1}; else, C = squeeze(coord(1,:,:)); end
if size(C,2) ~= 2, C = squeeze(C); end
parts = {C};
end

function Jk = subset(J, k)
Jk = J; Jk.features = J.features(k);
end

function Gt = mass_gt(m, lat, lon, mask, Re)
dphi = deg2rad(abs(lat(2)-lat(1))); dlam = deg2rad(abs(lon(2)-lon(1)));
[~, LA] = meshgrid(lon, lat);
A = (Re^2) * cosd(LA) * dphi * dlam;
Gt = sum(m(mask) .* A(mask), 'omitnan') * 1000 / 1e12;
end

function a = area_km2(lat, lon, mask, Re)
dphi = deg2rad(abs(lat(2)-lat(1))); dlam = deg2rad(abs(lon(2)-lon(1)));
[~, LA] = meshgrid(lon, lat);
A = (Re^2) * cosd(LA) * dphi * dlam / 1e6;
a = sum(A(mask));
end
