function [gt, rep] = gravisRegionChain(region, folder, gravisFolder, opts)
%GRAVISREGIONCHAIN Shared engine of greenlandChain / antarcticaChain.
%
%   [GT, REP] = shLowLevel.gravisRegionChain(REGION, FOLDER,
%   GRAVISFOLDER, OPTS) is the common implementation behind the two
%   public ice chains; call those instead. REGION is "greenland" or
%   "antarctica"; OPTS is the (already argument-validated) option
%   struct of the wrappers. Stages: gravisL2B (with GIA) -> EWH
%   synthesis on a GridStep grid -> pixel-wise bias + trend + annual +
%   semiannual fit -> open-ocean residual RMS -> sigma_trend =
%   sigma_monthly / sqrt(Sxx) unless NoiseLevel is given -> regularised
%   leakage inversion over the basin mask (union with NeighbourBoxes) ->
%   mass integral over the basin mask [Gt/yr].
%
%   Inputs
%     region        (1,1) string  "greenland" | "antarctica"
%     folder        (1,1) string  monthly solution folder
%     gravisFolder  (1,1) string  GravIS aux folder
%     opts          (1,1) struct  see greenlandChain / antarcticaChain
%
%   Outputs
%     gt   (1,1) double  basin-mask mass trend [Gt/yr]
%     rep  (1,1) struct  trendGrid, lat, lon, sigTrend, residRMS, mask,
%          unionMask, m, iterations, stoppedBy, basins (antarctica),
%          steps, l2b, version, created
%
%   Example
%     % called for you by the wrappers; direct use mirrors them:
%     % [gt, rep] = shLowLevel.gravisRegionChain("greenland", ser, gd, opts)
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.8.8).
arguments
    region (1,1) string {mustBeMember(region, ["greenland","antarctica"])}
    folder (1,1) string
    gravisFolder (1,1) string
    opts (1,1) struct
end
if isempty(opts.kn)
    error('shLowLevel:gravisRegionChain:noKn', ...
        ['Load Love numbers are required (kn=). The chain will not ' ...
         'assume a reference frame for you - see fetchLoveNumbers.']);
end
steps = strings(1, 0);
say = @(s) fprintf('  %s\n', s);
if opts.Quiet, say = @(s) []; end
% ---- basin polygons ([lon lat] GeoJSON -> masks)
if strlength(gravisFolder) == 0
    gravisFolder = shLowLevel.gravisDataFolder();
end
bf = char(opts.BasinFile);
if ~isfile(bf), bf = fullfile(char(gravisFolder), bf); end
if ~isfile(bf)
    error('shLowLevel:gravisRegionChain:missingBasins', ...
        ['Basin GeoJSON not found: %s. Fetch it from ' ...
         'https://gravis.gfz.de/basins/GIS or /basins/AIS.'], opts.BasinFile);
end
J = jsondecode(fileread(bf));
% ---- trend grid: compute or reuse
if ~isempty(opts.TrendGrid)
    S = opts.TrendGrid;
    lat = S.lat(:); lon = S.lon(:);
    steps(end+1) = "reused a precomputed trend grid (TrendGrid=)";
else
    [ts, l2b] = shLowLevel.gravisL2B(folder, gravisFolder, ...
        GIAFile = opts.GIAFile, SpanEnd = opts.SpanEnd, Quiet = opts.Quiet);
    ep = ts.epochs(:); T = ts.nEpochs; nmax = ts.nmax;
    lat = (-90 + opts.GridStep/2 : opts.GridStep : 90)';
    lon = (opts.GridStep/2 : opts.GridStep : 360)';
    w = ones(nmax + 1, 1);
    if startsWith(opts.Filter, "gauss")
        rkm = double(extractAfter(opts.Filter, "gauss"));
        w = shLowLevel.shGaussianWeights(nmax, rkm); w = w(:);
    elseif opts.Filter ~= "none"
        error('shLowLevel:gravisRegionChain:badFilter', ...
            'Filter must be "none" or "gaussN" for the grid stage, got %s.', ...
            opts.Filter);
    end
    E = zeros(numel(lat), numel(lon), T); P = [];
    for k = 1:T
        Ck = ts.Cs(:,:,k) .* w; Sk = ts.Ss(:,:,k) .* w;
        if isempty(P)
            [E(:,:,k), ~, ~, P] = shLowLevel.shSynthesis(Ck, Sk, opts.GM, ...
                opts.R, lat, lon, 'quantity','ewh','kn',opts.kn,'nmin',0);
        else
            E(:,:,k) = shLowLevel.shSynthesis(Ck, Sk, opts.GM, opts.R, ...
                lat, lon, 'quantity','ewh','kn',opts.kn,'nmin',0,'P',P);
        end
    end
    steps(end+1) = sprintf("synthesized %d EWH grids (%s, %g deg)", ...
        T, opts.Filter, opts.GridStep);
    A = [ones(T,1), ep-mean(ep), cos(2*pi*ep), sin(2*pi*ep), ...
         cos(4*pi*ep), sin(4*pi*ep)];
    X = reshape(E, [], T)';
    coef = A \ X;
    res = X - A * coef;
    trendGrid = reshape(coef(2,:), numel(lat), numel(lon));
    residRMS = reshape(rms(res, 1), numel(lat), numel(lon));
    Sxx = sum((ep - mean(ep)).^2);
    % open-ocean sigma and the TREND-matched noise (guide V4c). The
    % ocean mask is USER-SUPPLIED (oceanRMS contract - no coastline is
    % assumed), so the sigma_trend policy requires OceanMask=.
    if isempty(opts.OceanMask) && isnan(opts.NoiseLevel)
        error('shLowLevel:gravisRegionChain:noOcean', ...
            ['The default noise policy derives sigma_trend from the ' ...
             'open-ocean residual RMS, and the ocean mask is a caller ' ...
             'decision (see oceanRMS). Give OceanMask= (logical grid or ' ...
             '@(lat,lon) handle) or a numeric NoiseLevel=.']);
    end
    sigMon = NaN;
    if ~isempty(opts.OceanMask)
        sigMon = shLowLevel.oceanRMS(residRMS, lat, lon, opts.OceanMask);
    end
    S = struct('trendGrid', trendGrid, 'lat', lat, 'lon', lon, ...
        'residRMS', residRMS, 'sigTrend', sigMon / sqrt(Sxx));
    steps(end+1) = sprintf("pixel fit: sigma_mon %.4g m -> sigma_trend %.3g m/yr", ...
        sigMon, S.sigTrend);
end
% ---- masks
[LO, LA] = meshgrid(lon, lat);
lonW = mod(LO + 180, 360) - 180;
mask = false(size(LA)); cLon = zeros(numel(J.features), 1);
basinMasks = cell(numel(J.features), 1);
for k = 1:numel(J.features)
    g = J.features(k).geometry;
    if strcmp(g.type, 'Polygon'), parts = {squeeze(g.coordinates(1,:,:))};
    else, parts = cellfun(@(c) squeeze(c(1,:,:)), g.coordinates, 'uni', 0);
    end
    mk = false(size(LA)); allx = [];
    for q = 1:numel(parts)
        Pq = parts{q};
        c = atan2d(mean(sind(Pq(:,1))), mean(cosd(Pq(:,1))));
        px = mod(Pq(:,1) - c + 180, 360) - 180;
        gx = mod(lonW - c + 180, 360) - 180;
        mk = mk | inpolygon(gx, LA, px, Pq(:,2));
        allx = [allx; Pq(:,1)]; %#ok<AGROW>
    end
    basinMasks{k} = mk; mask = mask | mk;
    cLon(k) = atan2d(mean(sind(allx)), mean(cosd(allx)));
end
unionMask = mask;
for b = 1:size(opts.NeighbourBoxes, 1)
    bx = opts.NeighbourBoxes(b, :);
    unionMask = unionMask | (LA >= bx(1) & LA <= bx(2) & ...
        LO >= bx(3) & LO <= bx(4));
end
steps(end+1) = sprintf("mask: %d basin pixels, union %d", nnz(mask), nnz(unionMask));
% ---- inversion + mass
noise = opts.NoiseLevel;
if isnan(noise), noise = S.sigTrend; end
[m, info] = shLowLevel.leakageCorrect(S.trendGrid, lat, lon, ...
    Mask = unionMask, Filter = opts.Filter, NoiseLevel = noise, ...
    MaxIter = opts.MaxIter, Quiet = true);
steps(end+1) = sprintf("leakage inversion: %d iterations, stopped by %s", ...
    info.iterations, info.stoppedBy);
gt = local_mass(m, lat, lon, mask);
say(sprintf('%s: %+.1f Gt/yr (%d it, %s)', region, gt, ...
    info.iterations, info.stoppedBy));
rep = struct('trendGrid', S.trendGrid, 'lat', lat, 'lon', lon, ...
    'sigTrend', S.sigTrend, 'mask', mask, 'unionMask', unionMask, ...
    'm', m, 'iterations', info.iterations, 'stoppedBy', string(info.stoppedBy), ...
    'steps', steps, 'version', shLowLevel.version(), ...
    'created', string(datetime('now')));
if isfield(S, 'residRMS'), rep.residRMS = S.residRMS; end
if exist('l2b', 'var'), rep.l2b = l2b; end
if region == "antarctica"
    K = numel(basinMasks);
    nm = strings(K,1); gtb = zeros(K,1);
    for k = 1:K
        nm(k) = string(J.features(k).properties.name);
        gtb(k) = local_mass(m, lat, lon, basinMasks{k});
    end
    rep.basins = table(nm, gtb, 'VariableNames', {'name', 'gt'});
end
end

function gt = local_mass(m, lat, lon, mask)
Re = 6371e3;
dphi = deg2rad(abs(lat(2)-lat(1))); dlam = deg2rad(abs(lon(2)-lon(1)));
[~, LA] = meshgrid(lon, lat);
A = (Re^2) * cosd(LA) * dphi * dlam;
gt = sum(m(mask) .* A(mask), 'omitnan') * 1000 / 1e12;
end

