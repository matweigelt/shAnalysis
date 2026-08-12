function audit_gravis(variant)
%AUDIT_GRAVIS Reproduce the guide's GravIS Greenland chain (audit script).
%   Claude (Fable 5) audit, 2026-08-12. Not part of the toolbox.
if nargin < 1, variant = "all"; end
R = 6378136.3; GM = 3.986004415e14; Re = 6371e3;      % Re: chapter's area radius
aux = 'E:\DATAPOOL\GravityField\GravIS\';
ser = 'E:\DATAPOOL\GravityField\icgem\series\02_COST-G__COST-G_Grace-Grace-FO_RL02.1';
kn  = readmatrix('E:\DATAPOOL\GravityField\loveNumbers\loadLoveNumbers_Gegout97.txt', ...
                 FileType='text', NumHeaderLines=2);
S = load_or_build(aux, ser, kn, R, GM, true, "g445");
lat = S.lat; lon = S.lon; GIS = S.GIS; UNION = S.UNION;
fprintf('grid %dx%d, months=%d, span %.3f..%.3f\n', numel(lat), numel(lon), ...
    numel(S.epochs), S.epochs(1), S.epochs(end));
fprintf('sigMon(oceanRMS,eroded)=%.5f m  -> sigTrend=%.3e m/yr (chapter: 0.0115 -> 1.25e-4)\n', ...
    S.sigMon, S.sigTrend);
runs = struct('name',{},'Gt',{},'it',{},'stop',{});
    function do(name, Sx, mask, filt, nl, mi)
        [m, info] = shLowLevel.leakageCorrect(Sx.trendGrid, lat, lon, ...
            Filter=filt, Mask=mask, Gain=2, MaxIter=mi, NoiseLevel=nl, Quiet=true);
        Gt = mass_gt(m, lat, lon, GIS, Re);
        runs(end+1) = struct('name',name,'Gt',Gt,'it',info.iterations,'stop',info.stoppedBy);
        fprintf('%-34s %+8.1f Gt/yr  it=%3d  stop=%s\n', name, Gt, info.iterations, info.stoppedBy);
    end
switch variant
    case "all"
        do("HEADLINE union/g445/sigTrend",  S, UNION, "gauss445", S.sigTrend, 400);
        do("chapter-noise control 1.25e-4", S, UNION, "gauss445", 1.25e-4,    400);
        do("wrong noise: trend oceanRMS",   S, UNION, "gauss445", S.trendOceanRMS, 400);
        do("fixed 300 iters (no criterion)",S, UNION, "gauss445", 0,          300);
        do("Greenland-only mask",           S, GIS,   "gauss445", S.sigTrend, 400);
        do("MIS-declared: filtered obs+none",S, UNION, "none",    S.sigTrend, 400);
    case "v4"   % the chapter's V4 regime: UNFILTERED obs, Filter=none, fixed 300
        U = load_or_build(aux, ser, kn, R, GM, false, "unf");
        do("V4 union/none/300 (exp -240.2)", U, UNION, "none", 0, 300);
        do("V4 GIS-only/none/300 (exp -259.8)", U, GIS, "none", 0, 300);
    case "kn"   % Earth-model sensitivity: ak135 (CM, deg0-1 pinned) vs Gegout97
        kn2 = readmatrix('E:\DATAPOOL\GravityField\loveNumbers\ak135_CM_pinned01.txt');
        A2 = load_or_build(aux, ser, kn2, R, GM, true, "ak135");
        do("HEADLINE with kn=CM_ak135", A2, UNION, "gauss445", A2.sigTrend, 400);
        fprintf('kn delta at n=2..6: %s\n', mat2str(round(kn2(3:7)-kn(3:7),4)'));
    otherwise
        error("unknown variant");
end
fprintf('reference (my basin fit, span<=2023.10): -231.15 | full-span ref: -220.82\n');
save('E:\DATAPOOL\GravityField\GravIS\audit_runs.mat', 'runs');
end

function Gt = mass_gt(m, lat, lon, mask, Re)
dphi = deg2rad(abs(lat(2)-lat(1))); dlam = deg2rad(abs(lon(2)-lon(1)));
[~, LA] = meshgrid(lon, lat);
A = (Re^2) * cosd(LA) * dphi * dlam;                  % m^2 per cell
Gt = sum(m(mask) .* A(mask), 'omitnan') * 1000 / 1e12; % m EWH * m^2 * rho / 1e12
end

function S = load_or_build(aux, ser, kn, R, GM, filtered, tag)
cache = fullfile(aux, "audit_chain_cache_" + tag + ".mat");
if isfile(cache), S = load(cache); return; end
maxTableEpoch = 2023.099;                              % chapter regime (derived)
ts = shSeries.fromFolder(ser);
ts = ts.select(ts.epochs <= maxTableEpoch + 0.05);
T = ts.nEpochs; ep = ts.epochs(:);
% ---- aux tables, matched on MJD-of-begin from the filenames
LD = read_table(fullfile(aux,'GRAVIS-2B_COSTG_0200_GRACE+SLR_LOW_DEGREES_0001.dat'), 14);
GC = read_table(fullfile(aux,'GRAVIS-2B_COSTG_0200_GEOCENTER_0001.dat'), 11);
Cs = ts.Cs; Ss = ts.Ss;                                % (n+1,n+1,T)
% begins from the filenames (fromFolder drops names -> rebuild via dir)
df = dir(fullfile(ser, 'GSM-2_*.gfc')); [~,si] = sort({df.name}); df = df(si);
tok = regexp({df.name}, 'GSM-2_(\d{4})(\d{3})-(\d{4})(\d{3})', 'tokens', 'once');
nf = numel(tok); begF = zeros(nf,1); midF = zeros(nf,1);
yd = @(y) 365 + double(mod(y,4)==0 & (mod(y,100)~=0 | mod(y,400)==0));
for k = 1:nf
    y1=str2double(tok{k}{1}); d1=str2double(tok{k}{2});
    y2=str2double(tok{k}{3}); d2=str2double(tok{k}{4});
    t1 = y1 + (d1-1)/yd(y1); t2 = y2 + d2/yd(y2);
    begF(k) = t1; midF(k) = (t1+t2)/2;
end
keep = midF <= maxTableEpoch + 0.05; begF = begF(keep); midF = midF(keep);
assert(numel(begF)==T && max(abs(sort(midF)-ep)) < 5e-3, ...
    'filename mids do not line up with series epochs');
[begF, si] = sort(begF); %#ok<TRSRT> % chronological, same order as ep
iL = zeros(T,1); iG = zeros(T,1);
for k = 1:T
    [dL, iL(k)] = min(abs(LD(:,2) - begF(k)));
    [dG, iG(k)] = min(abs(GC(:,2) - begF(k)));
    assert(dL < 2e-3 && dG < 2e-3, ...
        'no aux row at begin %.4f (dL=%.4f dG=%.4f)', begF(k), dL, dG);
end
assert(numel(unique(iL))==T && numel(unique(iG))==T, 'ambiguous aux matching');
for k = 1:T
    Cs(3,1,k) = LD(iL(k),3);  Cs(4,1,k) = LD(iL(k),6);
    Cs(3,2,k) = LD(iL(k),9);  Ss(3,2,k) = LD(iL(k),12);
    Cs(2,1,k) = GC(iG(k),3);  Cs(2,2,k) = GC(iG(k),6);  Ss(2,2,k) = GC(iG(k),9);
end
Mn = shLowLevel.readSHM(fullfile(aux,'GRAVIS-2B_COSTG_0200_MEAN_2002095-2020091_NFIL_0001.gz'), Nmax=90);
Gi = shLowLevel.readSHM(fullfile(aux,'GRAVIS-2B_COSTG_0200_GIA_ICE-6G_D_VM5a_0001.gz'), Nmax=90);
nmax = min([size(Cs,1)-1, 90]);
Cs = Cs(1:nmax+1,1:nmax+1,:); Ss = Ss(1:nmax+1,1:nmax+1,:);
MnC = Mn.C(1:nmax+1,1:nmax+1); MnS = Mn.S(1:nmax+1,1:nmax+1);
GiC = Gi.C(1:nmax+1,1:nmax+1); GiS = Gi.S(1:nmax+1,1:nmax+1);
for k = 1:T
    Cs(:,:,k) = Cs(:,:,k) - MnC - (ep(k)-2011)*GiC;
    Ss(:,:,k) = Ss(:,:,k) - MnS - (ep(k)-2011)*GiS;
end
% ---- EWH synthesis, gauss445 pre-applied, 1-deg grid
lat = (-89.5:89.5)'; lon = (0.5:359.5)';
if filtered
    r445 = 4*111.195; w = shLowLevel.shGaussianWeights(nmax, r445); w = w(:);
else
    w = ones(nmax+1, 1);
end
E = zeros(numel(lat), numel(lon), T); P = [];
for k = 1:T
    Ck = Cs(:,:,k) .* w; Sk = Ss(:,:,k) .* w;
    if isempty(P)
        [E(:,:,k), ~, ~, P] = shLowLevel.shSynthesis(Ck, Sk, GM, R, lat, lon, ...
            'quantity','ewh','kn',kn,'nmin',0);
    else
        E(:,:,k) = shLowLevel.shSynthesis(Ck, Sk, GM, R, lat, lon, ...
            'quantity','ewh','kn',kn,'nmin',0,'P',P);
    end
end
% ---- pixel-wise fit: bias+trend+annual+semi
A = [ones(T,1), ep-mean(ep), cos(2*pi*ep), sin(2*pi*ep), cos(4*pi*ep), sin(4*pi*ep)];
X = reshape(E, [], T)';                                % T x npix
coef = A \ X;
resid = X - A*coef;
trendGrid = reshape(coef(2,:), numel(lat), numel(lon));
residRMS  = reshape(sqrt(mean(resid.^2,1)), numel(lat), numel(lon));
% ---- masks
gj = jsondecode(fileread(fullfile(fileparts(mfilename('fullpath')), ...
    'tests','test_data','gravis_gis_basins.geojson')));
[LO, LA] = meshgrid(lon, lat);
GIS = false(size(LA)); lonW = mod(LO+180,360)-180;
for f = 1:numel(gj.features)
    G = gj.features(f).geometry;
    polys = G.coordinates; if strcmp(G.type,'Polygon'), polys = {polys}; end
    for p = 1:numel(polys)
        ring = polys{p}; if iscell(ring), ring = ring{1}; end
        ring = squeeze(ring); if size(ring,2)~=2, ring = ring'; end
        GIS = GIS | inpolygon(lonW, LA, ring(:,1), ring(:,2));
    end
end
box = @(la1,la2,lo1,lo2) LA>=la1 & LA<=la2 & lonW>=lo1 & lonW<=lo2;
UNION = GIS | box(60,84,-128,-60) | box(63,67,-25,-13) | box(76,81,10,34);
% ---- noise
ocean = @(la,lo) trueOcean(la,lo);
[sigMon, oi] = shLowLevel.oceanRMS(residRMS, lat, lon, ocean);
Sxx = sum((ep-mean(ep)).^2);
sigTrend = sigMon / sqrt(Sxx);
trendOceanRMS = shLowLevel.oceanRMS(trendGrid, lat, lon, ocean);
S = struct('lat',lat,'lon',lon,'epochs',ep,'trendGrid',trendGrid, ...
    'residRMS',residRMS,'GIS',GIS,'UNION',UNION,'sigMon',sigMon, ...
    'sigTrend',sigTrend,'trendOceanRMS',trendOceanRMS,'Sxx',Sxx, ...
    'oceanPix',oi.nPixels);
save(cache, '-struct', 'S');
end

function M = read_table(f, ncol)
L = readlines(f); M = [];
for i = 1:numel(L)
    v = sscanf(L(i), '%f');
    if numel(v) >= ncol && v(1) > 40000 && v(1) < 90000, M(end+1,1:ncol) = v(1:ncol)'; end %#ok<AGROW>
end
end

function tf = trueOcean(la, lo)
% coarse open-ocean: |lat|<=55, exclude crude continent boxes
lo = mod(lo+180,360)-180;
land = (la>-35 & la<37  & lo>-17 & lo<52)  | ... % Africa/Arabia
       (la>5   & la<75  & lo>25  & lo<180) | ... % Eurasia
       (la>15  & la<72  & lo>-168& lo<-52) | ... % N America
       (la>-56 & la<13  & lo>-82 & lo<-34) | ... % S America
       (la>-45 & la<-10 & lo>112 & lo<154) | ... % Australia
       (la<-60);                                  % Antarctica
tf = abs(la)<=55 & ~land;
end
