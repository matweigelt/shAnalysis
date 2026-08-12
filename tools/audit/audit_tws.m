function audit_tws
%AUDIT_TWS GravIS TWS chain comparison (full span, box-synthesis route).
%   Chain: COST-G RL02.1 -> GravIS-2B corrections (aux tables; uncovered
%   epochs dropped) -> minus NFIL mean -> filter in {none, gauss300,
%   DDK3(=1d12), DDK5(=1d11)} -> per-basin cos-weighted 1-degree mask
%   means in cm EWH, synthesized on per-basin bounding boxes with
%   Legendre recycling (0.001 s/epoch). Basin masks: GravIS rivbas
%   polygons; convention verified to reproduce the official gridded
%   product's basin means to 2 decimals.
%   Reference: GravIS TWS portal series (COST-G, VDK5+VDK3 blend with
%   earthquake steps removed - both declared differences).
%   Claude (Fable 5), 2026-08-12.

aux = 'E:\DATAPOOL\GravityField\GravIS';
ser = 'E:\DATAPOOL\GravityField\icgem\series\02_COST-G__COST-G_Grace-Grace-FO_RL02.1';
ddk = 'E:\DATAPOOL\GravityField\DDK';
kn  = readmatrix(fullfile('C:\Users\matth\Documents\MATLAB\shAnalysis', ...
    'tests','test_data','loadLoveNumbers_Gegout97.txt'), FileType='text', NumHeaderLines=2);
GM = 3.986004415e14; R = 6378136.3;

S = load_or_build(aux, ser);
ep = S.ep(:); T = numel(ep); nmax = size(S.Cs,1)-1;
fprintf('corrected series: T=%d, %.3f..%.3f, nmax=%d (dropped: %d)\n', ...
    T, ep(1), ep(end), nmax, S.nDropped);

% filtered coefficient sets
w300 = shLowLevel.shGaussianWeights(nmax, 300); w300 = w300(:);
sets = struct('none', [], 'gauss300', [], 'DDK3', [], 'DDK5', []);
sets.none = S;
G = S; for t = 1:T
    G.Cs(:,:,t) = S.Cs(:,:,t) .* w300; G.Ss(:,:,t) = S.Ss(:,:,t) .* w300;
end
sets.gauss300 = G;
for pair = ["DDK3","Wbd_2-120.a_1d12p_4"; "DDK5","Wbd_2-120.a_1d11p_4"]'
    W = shLowLevel.readDDK(fullfile(ddk, char(pair(2))));
    ts = shSeries(S.Cs, Ss = S.Ss, Epochs = ep);
    ts = ts.applyDDK(W);
    F = S; F.Cs = ts.Cs; F.Ss = ts.Ss;
    sets.(char(pair(1))) = F;
end

% basins
JR = jsondecode(fileread(fullfile(aux, 'basins_rivbas.json')));
allN = arrayfun(@(f) string(f.properties.name), JR.features);
pick = ["Amazonas","Congo","Danube","Ganges","Lena","Mississippi River", ...
        "Niger","Ob","Yangtze River (Chang Jiang)","Yenisei","Zambezi"];
refT = [ 0.01  0.58 -0.71 -1.22 -0.13  0.07  0.94  0.12  0.24  0.09 -0.12];
refA = [20.6   5.1   6.8  14.4   4.0   5.9   8.8   6.5   3.9   4.7  14.4];
K = numel(pick);

variants = ["none","gauss300","DDK3","DDK5"];
A6 = @(t)[ones(size(t)), t-mean(t), cos(2*pi*t), sin(2*pi*t), cos(4*pi*t), sin(4*pi*t)];
res = struct();
for v = variants, res.(char(v)) = struct('tr', zeros(K,1), 'am', zeros(K,1)); end
for k = 1:K
    g = JR.features(find(allN == pick(k), 1)).geometry;
    if strcmp(g.type, 'Polygon'), parts = {squeeze(g.coordinates(1,:,:))};
    else, parts = cellfun(@(c) squeeze(c(1,:,:)), g.coordinates, 'uni', 0); end
    Pall = vertcat(parts{:});
    la1 = floor(min(Pall(:,2)))-1; la2 = ceil(max(Pall(:,2)))+1;
    latB = (la1+0.5 : la2-0.5)';
    c0 = atan2d(mean(sind(Pall(:,1))), mean(cosd(Pall(:,1))));
    lw = mod(Pall(:,1) - c0 + 180, 360) - 180;
    lo1 = floor(min(lw))-1; lo2 = ceil(max(lw))+1;
    lonB = mod(c0 + (lo1+0.5 : lo2-0.5)', 360);      % absolute lon of box cols
    [LOb, LAb] = meshgrid(1:numel(lonB), latB);      % LOb: column index only
    lonWc = mod((lo1+0.5:lo2-0.5) , 1e9);            % centred offsets
    [LOc, ~] = meshgrid(lo1+0.5:lo2-0.5, latB);      % centred lon grid
    mkB = false(size(LAb));
    for q = 1:numel(parts)
        Pq = parts{q};
        pl = mod(Pq(:,1) - c0 + 180, 360) - 180;
        mkB = mkB | inpolygon(LOc, LAb, pl, Pq(:,2));
    end
    wB = cosd(LAb);
    for v = variants
        C3 = sets.(char(v)).Cs; S3 = sets.(char(v)).Ss;
        s = zeros(T, 1); Pl = [];
        for t = 1:T
            if isempty(Pl)
                [E,~,~,Pl] = shLowLevel.shSynthesis(C3(:,:,t), S3(:,:,t), GM, R, ...
                    latB, lonB, 'quantity','ewh','kn',kn,'nmin',0);
            else
                E = shLowLevel.shSynthesis(C3(:,:,t), S3(:,:,t), GM, R, ...
                    latB, lonB, 'quantity','ewh','kn',kn,'nmin',0,'P',Pl);
            end
            s(t) = sum(E(mkB) .* wB(mkB)) / sum(wB(mkB)) * 100;
        end
        x = A6(ep) \ s;
        res.(char(v)).tr(k) = x(2);
        res.(char(v)).am(k) = hypot(x(3), x(4));
    end
end

fprintf('\n%-27s | trend cm/yr:  none   g300   DDK3   DDK5    ref | amp cm:  none  g300  DDK3  DDK5   ref\n','basin');
for k = 1:K
    fprintf('%-27s | %+12.2f %+6.2f %+6.2f %+6.2f %+6.2f | %13.1f %5.1f %5.1f %5.1f %5.1f\n', ...
        pick(k), res.none.tr(k), res.gauss300.tr(k), res.DDK3.tr(k), res.DDK5.tr(k), refT(k), ...
        res.none.am(k), res.gauss300.am(k), res.DDK3.am(k), res.DDK5.am(k), refA(k));
end
fprintf('\n');
for v = variants
    r = res.(char(v));
    fprintf('%-9s: RMS dTrend %.3f cm/yr | amp ratio med %.3f (range %.2f..%.2f)\n', v, ...
        rms(r.tr - refT(:)), median(r.am ./ refA(:)), min(r.am./refA(:)), max(r.am./refA(:)));
end
save(fullfile(aux, 'audit_tws_runs.mat'), 'res', 'pick', 'refT', 'refA');
end

% ------------------------------------------------------------------ local
function S = load_or_build(aux, ser)
cache = fullfile(aux, 'audit_tws_cache_coeff.mat');
if isfile(cache), S = load(cache); return; end
ts = shSeries.fromFolder(ser);
ep = ts.epochs(:); T0 = ts.nEpochs;
LD = read_table(fullfile(aux, 'GRAVIS-2B_COSTG_0200_GRACE+SLR_LOW_DEGREES_0001.dat'), 14);
GC = read_table(fullfile(aux, 'GRAVIS-2B_COSTG_0200_GEOCENTER_0001.dat'), 11);
df = dir(fullfile(ser, 'GSM-2_*.gfc')); [~, si] = sort({df.name}); df = df(si);
tok = regexp({df.name}, 'GSM-2_(\d{4})(\d{3})-(\d{4})(\d{3})', 'tokens', 'once');
yd = @(y) 365 + double(mod(y,4)==0 & (mod(y,100)~=0 | mod(y,400)==0));
nf = numel(tok); begF = zeros(nf, 1); midF = zeros(nf, 1);
for k = 1:nf
    y1 = str2double(tok{k}{1}); d1 = str2double(tok{k}{2});
    y2 = str2double(tok{k}{3}); d2 = str2double(tok{k}{4});
    begF(k) = y1 + (d1-1)/yd(y1); midF(k) = (begF(k) + y2 + d2/yd(y2)) / 2;
end
[~, si2] = sort(midF); begF = begF(si2); midF = midF(si2);
assert(numel(midF) == T0 && max(abs(midF - ep)) < 5e-3, 'filename/epoch mismatch');
Cs = ts.Cs; Ss = ts.Ss;
keep = true(T0, 1); iL = zeros(T0, 1); iG = zeros(T0, 1);
for k = 1:T0
    [dL, iL(k)] = min(abs(LD(:, 2) - begF(k)));
    [dG, iG(k)] = min(abs(GC(:, 2) - begF(k)));
    if dL > 2e-3 || dG > 2e-3, keep(k) = false; end  % correction tables trail
end
for k = find(keep)'
    Cs(3,1,k) = LD(iL(k),3);  Cs(4,1,k) = LD(iL(k),6);
    Cs(3,2,k) = LD(iL(k),9);  Ss(3,2,k) = LD(iL(k),12);
    Cs(2,1,k) = GC(iG(k),3);  Cs(2,2,k) = GC(iG(k),6);  Ss(2,2,k) = GC(iG(k),9);
end
Cs = Cs(:,:,keep); Ss = Ss(:,:,keep); ep = ep(keep);
Mn = shLowLevel.readSHM(fullfile(aux, 'GRAVIS-2B_COSTG_0200_MEAN_2002095-2020091_NFIL_0001.gz'), Nmax = 90);
nmax = min(size(Cs,1)-1, 90);
Cs = Cs(1:nmax+1, 1:nmax+1, :); Ss = Ss(1:nmax+1, 1:nmax+1, :);
for k = 1:numel(ep)
    Cs(:,:,k) = Cs(:,:,k) - Mn.C(1:nmax+1, 1:nmax+1);
    Ss(:,:,k) = Ss(:,:,k) - Mn.S(1:nmax+1, 1:nmax+1);
end
S = struct('Cs', Cs, 'Ss', Ss, 'ep', ep, 'nDropped', T0 - numel(ep));
save(cache, '-struct', 'S', '-v7.3');
end

function M = read_table(fp, ncol)
M = readmatrix(fp, FileType = 'text', CommentStyle = '#');
M = M(:, 1:ncol);
end
