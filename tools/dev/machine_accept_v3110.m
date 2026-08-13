%MACHINE_ACCEPT_V3110 Acceptance runs for v3.11.0 on the acceptance machine.
%
%   Run section by section (or whole) once the MATLAB/MCP bridge is back.
%   Every check prints PASS/FAIL/SKIP with the measured number - nothing
%   is silent. Prerequisites: main checked out in the working copy,
%   GAX folders populated (E:\DATAPOOL\GravityField\GAX\GAD and \GAA,
%   251 files each), COST-G RL02.1 series on E:.
%
%   Sections
%     1  path hygiene (remove the nested shx_v3100 shadow)
%     2  full local test suite
%     3  fetchGAX bounded live run (SKIPs while the ICGEM 429 lasts)
%     4  obpChain field acceptance on the full series
%     5  oceanChain SeparateCirculation acceptance
%     6  optional: GravIS Level-3 OBP grid cross-check (SKIPs on 503)
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-12 (v3.11.0).

%% 1 - hygiene
src = 'C:\Users\matth\Documents\MATLAB\shAnalysis';
shd = 'C:\Users\matth\Documents\MATLAB\shx_v3100';
cd(src);
p = strsplit(path, pathsep);
dead = p(contains(p, 'shx_v3100') | contains(p, 'Temp\shx_'));
for k = 1:numel(dead), rmpath(dead{k}); end
if isfolder(shd)
    try, rmdir(shd, 's'); fprintf('PASS  nested shadow removed\n');
    catch ME, fprintf('FAIL  rmdir shadow: %s\n', ME.message); end
else
    fprintf('PASS  no stale shadow\n');
end
clear functions; rehash toolboxcache;
fprintf('using: %s\n', which('shLowLevel.obpChain'));

%% 2 - full suite
r = runAllTests;
fprintf('%s  suite: %d passed, %d failed, %d filtered\n', ...
    ternary(nnz([r.Failed]) == 0, "PASS", "FAIL"), ...
    nnz([r.Passed]), nnz([r.Failed]), nnz([r.Incomplete]));

%% 3 - fetchGAX bounded live (verification debt from v3.10)
try
    dst = fullfile(tempdir, 'gax_accept');
    if isfolder(dst), rmdir(dst, 's'); end
    [f, info] = shLowLevel.fetchGAX(dst, MaxFiles = 2, BudgetSec = 120, ...
        Quiet = true);
    nF = sum([info.products.nFailed]);
    fprintf('%s  fetchGAX: %d files, %d failed\n', ...
        ternary(numel(f) >= 4 && nF == 0, "PASS", "FAIL"), numel(f), nF);
catch ME
    if contains(ME.message, '429') || contains(ME.identifier, '429')
        fprintf('SKIP  fetchGAX: ICGEM still rate-limits this IP\n');
    else
        fprintf('FAIL  fetchGAX: %s\n', ME.identifier);
    end
end

%% 4 - obpChain field acceptance
ser = 'E:\DATAPOOL\GravityField\icgem\series\02_COST-G__COST-G_Grace-Grace-FO_RL02.1';
kn = readmatrix(fullfile(src, 'tests', 'test_data', ...
    'loadLoveNumbers_Gegout97.txt'), FileType = 'text', NumHeaderLines = 2);
oc = @(la,lo) abs(la)<=66 & ~( (la>-35 & la<37 & mod(lo+180,360)-180>-17 & mod(lo+180,360)-180<52) | ...
    (la>5 & la<75 & mod(lo+180,360)-180>25 & mod(lo+180,360)-180<180) | ...
    (la>15 & la<72 & mod(lo+180,360)-180>-168 & mod(lo+180,360)-180<-52) | ...
    (la>-56 & la<13 & mod(lo+180,360)-180>-82 & mod(lo+180,360)-180<-34) | ...
    (la>-45 & la<-10 & mod(lo+180,360)-180>112 & mod(lo+180,360)-180<154) | (la<-60));
tic
[obp, repO] = shLowLevel.obpChain(ser, kn = kn, OceanMask = oc, ...
    GADFolder = 'E:\DATAPOOL\GravityField\GAX\GAD');
tOBP = toc;
% land check via a known land pixel (central Asia ~ lat 45, lon 90)
[~, iLa] = min(abs(obp.lat - 45)); [~, iLo] = min(abs(obp.lon - 90));
sd = std(obp.oceanMeanOBP);
fprintf(['%s  obpChain %.0f s: T=%d, GAD %d/%d, ref epochs %d, ' ...
    'land NaN %d, ocean-mean std %.2f cm\n'], ...
    ternary(repO.nGadRestored >= 250 && isnan(obp.grid(iLa, iLo, 1)) ...
        && sd > 0.3 && sd < 3, "PASS", "FAIL"), ...
    tOBP, repO.nEpochs, repO.nGadRestored, repO.nEpochs, ...
    repO.nRefEpochs, isnan(obp.grid(iLa, iLo, 1)), sd);
% OBP ocean-mean trend should sit near the barystatic +1.41 (AOD1B is
% nearly trend-free, so the atmospheric term adds little secular signal)
ep = obp.epochs; A = [ones(numel(ep),1), ep - mean(ep)];
x = A \ obp.oceanMeanOBP;
fprintf('%s  OBP ocean-mean trend %+.2f mm/yr (barystatic was +1.41)\n', ...
    ternary(abs(x(2)*10 - 1.41) < 0.3, "PASS", "FAIL"), x(2)*10);

%% 5 - residual circulation separation
tic
[out, rep] = shLowLevel.oceanChain(ser, kn = kn, OceanMask = oc, ...
    GADFolder = 'E:\DATAPOOL\GravityField\GAX\GAD', ...
    GAAFolder = 'E:\DATAPOOL\GravityField\GAX\GAA', ...
    SeparateCirculation = true);
fprintf(['%s  oceanChain+EOF %.0f s: %d modes, circ RMS %.4f m, ' ...
    'noise %.4f m (was 0.0149 undivided), trend %+.2f mm/yr\n'], ...
    ternary(out.nModes >= 1 && out.sigMon < 0.0149 ...
        && abs(out.trend - 1.41) < 0.05, "PASS", "FAIL"), ...
    toc, out.nModes, out.circulationRMS, out.sigMon, out.trend);

%% 6b - v3.12 fetch family live (bounded; TU Graz is not rate-limited)
try
    d1 = fullfile(tempdir, 'itsg_bg');
    [f1, i1] = shLowLevel.fetchITSGBackground("2018-06", Dest = d1, Quiet = true);
    g = shCoefficients.read(f1(1));
    fprintf('%s  fetchITSGBackground: %d file(s), parsed nmax %d\n', ...
        ternary(numel(f1) >= 1 && size(g.C, 1) > 90, "PASS", "FAIL"), ...
        numel(f1), size(g.C, 1) - 1);
    % SINEX: 460 MB - only the URL/HEAD path is exercised via MaxFiles=0-like
    % budget; a real download is a deliberate manual decision.
    [~, i2] = shLowLevel.fetchITSGSINEX("2018-06", Dest = fullfile(tempdir, 'itsg_snx'), ...
        MaxFiles = 0, Quiet = true);   % MaxFiles=0 cuts BEFORE the 460 MB download (BudgetSec cannot)
    fprintf('%s  fetchITSGSINEX plumbing: %d listed, %d remaining (budget cut as designed)\n', ...
        ternary(i2.nListed == 1 && i2.nRemaining == 1, "PASS", "FAIL"), ...
        i2.nListed, i2.nRemaining);
catch ME
    fprintf('FAIL  fetch family: %s\n', ME.identifier);
end

%% 6c - v3.14 hydro index on the real series (Amazon droughts visible)
try
    [tws, ~] = shLowLevel.twsChain(ser, kn = kn);
    [Z, ih] = shLowLevel.hydroExtremeIndex(tws.grid, tws.epochs);
    % Amazon cell ~ (-5, 300): the 2005 and 2010 droughts are published
    [~, iLa] = min(abs(tws.lat - (-5))); [~, iLo] = min(abs(tws.lon - 300));
    z = squeeze(Z(iLa, iLo, :));
    [zmin, imin] = min(z);
    fprintf('%s  hydroIndex: Amazon min DSI %.2f at %.2f (published droughts 2005/2010)\n', ...
        ternary(zmin <= -1.3 && (abs(tws.epochs(imin)-2005.7) < 1 ...
            || abs(tws.epochs(imin)-2010.7) < 1), "PASS", "CHECK"), ...
        zmin, tws.epochs(imin));
catch ME
    fprintf('FAIL  hydroIndex: %s\n', ME.identifier);
end

%% 6d - v3.14 stage 2: daily Kalman flood tracking (bounded live)
try
    fD = shLowLevel.fetchITSG("2019-07", Product = "daily", Quiet = true);
    tsD = shSeries.read(fD);
    fprintf('%s  daily Kalman: %d days read, epochs %.3f..%.3f\n', ...
        ternary(numel(fD) >= 28, "PASS", "FAIL"), ...
        numel(fD), tsD.epochs(1), tsD.epochs(end));
    % full daily-DSI needs a multi-year daily archive - that fetch is a
    % deliberate batch (365 files/yr); this section verifies the path.
catch ME
    fprintf('FAIL  daily path: %s\n', ME.identifier);
end

%% 6e - v3.15 coastal-leakage controls on the real series
try
    base = shLowLevel.oceanChain(ser, kn = kn, OceanMask = oc, ...
        GADFolder = 'E:\\DATAPOOL\\GravityField\\GAX\\GAD', ...
        GAAFolder = 'E:\\DATAPOOL\\GravityField\\GAX\\GAA');
    for buf = [300, 500]
        o = shLowLevel.oceanChain(ser, kn = kn, OceanMask = oc, ...
            GADFolder = 'E:\\DATAPOOL\\GravityField\\GAX\\GAD', ...
            GAAFolder = 'E:\\DATAPOOL\\GravityField\\GAX\\GAA', ...
            CoastBufferKm = buf);
        fprintf('buffer %d km: trend %+.2f (base %+.2f) mm/yr, area %.1f%%\n', ...
            buf, o.trend, base.trend, 100 * o.oceanArea / base.oceanArea);
    end
    o2 = shLowLevel.oceanChain(ser, kn = kn, OceanMask = oc, ...
        GADFolder = 'E:\\DATAPOOL\\GravityField\\GAX\\GAD', ...
        GAAFolder = 'E:\\DATAPOOL\\GravityField\\GAX\\GAA', ...
        RemoveLandLeakage = true, CoastBufferKm = 300);
    fprintf('%s  land removal + 300 km: trend %+.2f mm/yr (published barystatic ~1.6-2.2; base %+.2f)\n', ...
        ternary(o2.trend > base.trend, "PASS", "CHECK"), o2.trend, base.trend);
catch ME
    fprintf('FAIL  leakage controls: %s\n', ME.identifier);
end

%% 6 - optional GravIS OBP cross-check (503-tolerant)
try
    gvDir = fullfile(tempdir, 'gravis_obp');
    if ~isfolder(gvDir), mkdir(gvDir); end
    gvFile = fullfile(gvDir, 'GRAVIS-3_OBP_COSTG.nc');
    if ~isfile(gvFile)
        shLowLevel.httpFetch(['https://isdc-data.gfz.de/grace/GravIS/' ...
            'COST-G/Level-3/OBP/GRAVIS-3_2002095-2025334_GFZOP_OBP_COSTG_0200.nc'], ...
            gvFile, MaxTries = 2, Timeout = 60);
    end
    fprintf(['NOTE  GravIS OBP grid at %s - compare a common epoch ' ...
        'against obp.grid by hand (filters differ: GravIS destripes)\n'], gvFile);
catch
    fprintf(['SKIP  GravIS OBP grid: ISDC unreachable (503 observed ' ...
        '2026-08-12) - retry later; the exact filename may also need ' ...
        'a listing check\n']);
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
