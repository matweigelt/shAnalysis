function reg = demo_shAnalysis(cases, opts)
%DEMO_SHANALYSIS Selectable demonstrations of every toolbox capability.
%
%   demo_shAnalysis                 runs the "core" set (D01 D02 D04 D05 D15)
%   demo_shAnalysis("all")          runs every case
%   demo_shAnalysis("list")        prints the case table, returns registry
%   demo_shAnalysis(["D04","D15"]) runs selected cases
%   demo_shAnalysis(..., Visible (true)=false)  hidden figures (CI/smoke tests)
%   demo_shAnalysis(..., OutDir (string(tempdir))=tempdir) target for D16 file exports
%   demo_shAnalysis(..., StopOnError (false)=true)  rethrow instead of the
%       default fail-and-continue (each case runs in try/catch; a
%       summary lists ok/failed cases at the end)
%
%   Cases (see the workflow guide's Demo Gallery for rendered examples):
%     D01  Read & spectral diagnostics   triangle, spectrum, Kaula, crossover
%     D02  Synthesis quantities & maps   geoid/EWH/anomaly, Hammer, Height
%     D03  Normal field -> geoid         subtractNormalField, +-100 m map
%     D04  Filter comparison             Gaussian vs fan vs destripe (vs DDK)
%     D05  Climatology & basin series    trend/annual maps, gap, sigma band
%     D06  tvANS pipeline                noise blocks, VCE, filtered series
%     D07  Analysis (inverse problem)    ring roundtrip, scattered + Kaula
%     D08  Basin tools                   kernel taper, deconvolve vs scaling
%     D09  Uncertainty                   errorMap vs mcPropagate, covariance
%     D10  Load deformation              station up/north/east series
%     D11  Gradient tensor               six NEU components at 250 km
%     D12  EOF analysis                  two-mode recovery, PCs
%     D13  Trend breakpoints             hinge fit, F-test map
%     D14  Multi-center combination      VCE weights, inter-center corr
%     D15  Sea-level fingerprint         S/eustatic map, conservation
%     D16  Export                        writeGrid netCDF, writeAnimation MP4
%
%   Data policy: D01-D04 use the REAL ITSG files shipped in
%   tests/test_data (GRACE 2008-04, GRACE-FO 2025-12, DDK3 Wbd) - D02/
%   D04 show the real 17.7-yr mass-change difference with real stripes.
%   D05/D06 use a real monthly series once shLowLevel.fetchITSG(years) has
%   populated tests/test_data/itsg_series (>= 24 months). D12/D13/D14
%   stay synthetic BY DESIGN: they demonstrate recovery of KNOWN truth
%   (modes, breaks, noise factors), which real data cannot provide.
%   Every case falls back to synthetic data and stays self-contained.
%   Love numbers in the demos are SYNTHETIC - real work requires a real
%   loading model (PREM table), supplied by you.
%
%   Claude (Fable 5), 2026-08-07 (v2.4).
%   Outputs
%     reg        (D x 3) table   demo registry: id, title, exercised API
%
%   Inputs
%     cases   (1 x n) string  demo IDs to run, e.g. "D01" or ["D01" "D13"];
%             "core" expands to the curated smoke set, "all" runs the
%             full registry in order
%
%   Options
%     Visible     (1 x 1) logical  show figures (false keeps them
%                 off-screen for batch/CI runs)
%     OutDir      (1 x 1) string   output folder for exported files;
%                 "" = a fresh folder under tempdir
%     StopOnError (1 x 1) logical  true stops at the first failing
%                 case; false (default) prints FAILED and continues -
%                 the fail-and-continue contract the suite exercises
%
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
arguments
    cases (1,:) string = "core"
    opts.Visible (1,1) logical = true
    opts.OutDir (1,1) string = string(tempdir)
    opts.StopOnError (1,1) logical = false
end
reg = registry();
if any(cases == "list")
    fprintf('\n%-5s %-34s %s\n', 'ID', 'Title', 'Exercises');
    for k = 1:numel(reg)
        fprintf('%-5s %-34s %s\n', reg(k).id, reg(k).title, reg(k).fns);
    end
    return
end
if any(cases == "core"), cases = ["D01" "D02" "D04" "D05" "D15"]; end
if any(cases == "all"),  cases = [reg.id]; end
okList = strings(1, 0); failList = strings(1, 0);
for c = cases
    k = find([reg.id] == c, 1);
    assert(~isempty(k), 'demo:unknownCase', 'Unknown case %s (try "list").', c);
    fprintf('\n=== %s: %s ===\n', reg(k).id, reg(k).title);
    try
        reg(k).run(opts.Visible, char(opts.OutDir));
        okList(end+1) = reg(k).id; %#ok<AGROW>
    catch err
        failList(end+1) = reg(k).id; %#ok<AGROW>
        fprintf(2, '  FAILED (%s): %s\n', err.identifier, err.message);
        if ~isempty(err.stack)
            fprintf(2, '    at %s (line %d)\n', ...
                err.stack(1).name, err.stack(1).line);
        end
        if opts.StopOnError, rethrow(err); end
        fprintf('  continuing with the next case (StopOnError=false)\n');
    end
end
fprintf('\n---- demo summary: %d ok', numel(okList));
if ~isempty(failList)
    fprintf(', %d FAILED (%s)', numel(failList), strjoin(failList, ', '));
end
fprintf(' ----\n');
v = shLowLevel.version();
fprintf('Provenance: %s v%s, Claude (Fable 5), %s.\n', ...
    v.Name, v.Version, v.Date);
if ~isempty(failList) && nargout == 0
    warning('demo:casesFailed', 'Failed cases: %s. Rerun individually or with StopOnError=true for debugging.', ...
        strjoin(failList, ', '));
end
end

% ---------------------------------------------------------------- registry
function reg = registry()
reg = [ ...
    mk("D01", "Read & spectral diagnostics", "shReadGFC, shDegreeRMS, plotSHCoeffTriangle, plotSHSpectrum", @d01)
    mk("D02", "Synthesis quantities & maps", "shSynthesis, kernelFactors, plotSHMap, Height", @d02)
    mk("D03", "Normal field -> geoid", "normalFieldCS, subtractNormalField, toReference", @d03)
    mk("D04", "Filter comparison", "shGaussianFilter, shFanFilter, shDestripe, readDDK, triangle RefC/RefS", @d04)
    mk("D05", "Climatology & basin series", "shSeries, climatology, basinKernel, plotBasinSeries", @d05)
    mk("D06", "tvANS pipeline", "buildNoiseCov, tvANSFilter, vceRescale, opApply", @d06)
    mk("D07", "Analysis (inverse problem)", "shAnalysisGrid, analysis, Kaula", @d07)
    mk("D08", "Basin tools", "basinKernel, basinDeconvolve, basinScaling", @d08)
    mk("D09", "Uncertainty", "errorMap, mcPropagate, plotCovariance", @d09)
    mk("D10", "Load deformation", "shSynthesisDeformation, deformation", @d10)
    mk("D11", "Gradient tensor", "shSynthesisGradientTensor, legendreALFDeriv", @d11)
    mk("D12", "EOF analysis", "eofAnalysis", @d12)
    mk("D13", "Trend breakpoints", "trendBreaks, fitDeterministicModel Breaks", @d13)
    mk("D14", "Multi-center combination", "combineCenters, buildNoiseCov", @d14)
    mk("D15", "Sea-level fingerprint", "seaLevelFingerprint, evalMask, synthesisMatrix", @d15)
    mk("D16", "Export", "writeGrid, writeAnimation", @d16)];
end
function s = mk(id, title, fns, run)
s = struct('id', id, 'title', title, 'fns', fns, 'run', run);
end
function f = newfig(vis, name)
f = figure('Visible', onoff(vis), 'Name', name, 'Position', [60 60 1000 460]);
end
function s = onoff(tf), if tf, s = 'on'; else, s = 'off'; end, end

function [kn, hn, ln] = synthLove(L)
% SYNTHETIC Love numbers - demo only; real work needs a PREM table.
n = (0:L)';
kn = -0.30 * n ./ (n + 3.0);
hn = -0.90 * n ./ (n + 2.0) - 0.1;
ln = -0.05 * n ./ (n + 4.0) - 0.01;
end

function g = demoField(L, seed)
% Kaula-decaying synthetic static-like field
rng(seed);
n = (0:L)';
scale = 1e-5 ./ max(n, 1).^2 .* sqrt(2*n + 1) / sqrt(2);
C = randn(L+1) .* scale .* tril(true(L+1));
S = randn(L+1) .* scale .* tril(true(L+1)); S(:, 1) = 0;
C(1, 1) = 1; C(2, :) = 0; S(2, :) = 0;
g = shCoefficients(C, S, Name = "synthetic demo field");
end

function ts = demoSeries(L, T, seed)
rng(seed);
n = (0:L)'; degScale = 1e-9 ./ (1 + n).^1.5;
mask = tril(true(L+1));
trendC = randn(L+1) .* degScale .* mask * 0.3;
annC = randn(L+1) .* degScale .* mask;
annS = randn(L+1) .* degScale .* mask; annS(:, 1) = 0;
ep = 2019 + (0:T-1)'/12;
Cs = zeros(L+1, L+1, T); Ss = Cs;
for t = 1:T
    tc = ep(t) - mean(ep);
    stripe = zeros(L+1); stripe(:, 9:2:end) = randn(L+1, numel(9:2:L+1)) * 2e-10;
    Cs(:,:,t) = trendC*tc + annC*cos(2*pi*tc) + (randn(L+1)*5e-11 + stripe) .* mask;
    Ss(:,:,t) = annS*sin(2*pi*tc) + randn(L+1)*5e-11 .* mask; Ss(:,1,t) = 0;
end
ts = shSeries(Cs, Ss = Ss, Epochs = ep, ProductType = "GSM");
end

function [g1, g2] = tryRealPair()
% the two shipped real months: GRACE 2008-04 and GRACE-FO 2025-12 (n60)
d = shLowLevel.testDataDir();
g1 = []; g2 = [];
f1 = dir(fullfile(d, 'ITSG-Grace2018_n60_*.gfc'));
f2 = dir(fullfile(d, 'ITSG-Grace_operational_n60_*.gfc'));
if isempty(f1) || isempty(f2), return, end
try
    g1 = shCoefficients.read(fullfile(f1(1).folder, f1(1).name));
    g2 = shCoefficients.read(fullfile(f2(1).folder, f2(1).name));
catch
    g1 = []; g2 = [];
end
end

function ts = tryRealSeries(Lcut)
% real monthly series if shLowLevel.fetchITSG has populated itsg_series/
ts = [];
cands = [string(fullfile(shLowLevel.dataFolder(), 'itsg_series')), ...
    string(fullfile(shLowLevel.testDataDir(), ...
    'itsg_series'))];                              % new + legacy location
for d = cands
    if isfolder(d) && numel(dir(fullfile(d, '*.gfc*'))) >= 24
        try
            ts = shSeries.fromFolder(d, Truncate = Lcut);
            return
        catch
            ts = [];
        end
    end
end
end

function g = tryRealITSG()
d = shLowLevel.testDataDir();
fl = [dir(fullfile(d, '*ITSG*2008*.gfc*')); dir(fullfile(d, '*ITSG*.gfc*'))];
if isempty(fl), g = []; return, end
try
    % shCoefficients.read (NOT the low-level shReadGFC, which returns a
    % plain struct without methods - root cause of the v2.4.1 D03 crash)
    g = shCoefficients.read(fullfile(fl(1).folder, fl(1).name));
catch
    g = [];
end
end

% ------------------------------------------------------------------- cases
function d01(vis, ~)
g = tryRealITSG();
if isempty(g), g = demoField(60, 11); fprintf('  (synthetic field)\n');
else, fprintf('  (real ITSG file)\n'); end
newfig(vis, 'D01 triangle & spectrum');
subplot(1, 2, 1);
shLowLevel.plotSHCoeffTriangle(g.C, g.S, 'ax', gca);
subplot(1, 2, 2);
if ~isempty(g.sigmaC)
    sC = g.sigmaC; sS = g.sigmaS;             % real formal errors
else
    sC = abs(g.C) * 0.02 + 1e-13; sS = sC;    % mock if none
end
spec = shLowLevel.shDegreeRMS(g.C, g.S, 'R', g.R, 'sigmaC', sC, 'sigmaS', sS);
shLowLevel.plotSHSpectrum(spec, 'ax', gca, 'Kaula', 1e-5, 'MarkCrossover', true);
end

function d02(vis, ~)
[g1, g2] = tryRealPair();
if isempty(g1)
    g = demoField(45, 12); gd = demoField(45, 121) * 1e-3;
    fprintf('  (synthetic fields)\n');
else
    g = g1.subtractNormalField();             % real geoid, +-100 m
    gd = g2 - g1;                             % real 2008->2025 mass change
    fprintf('  (real ITSG pair: %s-epoch difference spans %.1f yr)\n', ...
        'GRACE/GRACE-FO', g2.epoch - g1.epoch);
end
kn = synthLove(g.nmax);
lat = -89:2:89; lon = 0:2:358;
newfig(vis, 'D02 quantities');
subplot(1, 2, 1);
g.map(lat, lon, quantity = "geoid", Units = "m", Title = "geoid, plate");
subplot(1, 2, 2);
gd.gaussian(350).map(lat, lon, quantity = "ewh", kn = kn, ...
    Units = "m EWH", Projection = "hammer", ...
    Title = "mass change 2008-2025, Gaussian 350 (EWH; synthetic kn!)");
newfig(vis, 'D02 upward continuation');
g = g.truncate(min(g.nmax, 45));
ga0 = g.synthesis(lat, lon, quantity = "gravity_anomaly");
ga4 = g.synthesis(lat, lon, quantity = "gravity_anomaly", Height = 400e3);
subplot(1, 2, 1); shLowLevel.plotSHMap(ga0*1e5, lat, lon, Units="mGal", Title="anomaly, surface", ax=gca);
subplot(1, 2, 2); shLowLevel.plotSHMap(ga4*1e5, lat, lon, Units="mGal", Title="anomaly at 400 km", ax=gca);
fprintf('  attenuation at 400 km: rms ratio %.3f\n', sqrt(mean(ga4(:).^2))/sqrt(mean(ga0(:).^2)));
end

function d03(vis, ~)
g = tryRealITSG();
if isempty(g)
    % synthetic full field: rescaled normal field + Kaula residual
    gr = demoField(40, 13);
    [CnE, iN] = shLowLevel.normalFieldCS(40);
    Cell = zeros(41); Cell(:, 1) = CnE;
    CnR = shLowLevel.rescaleGMR(Cell, zeros(41), iN.GM, iN.a, gr.GM, gr.R);
    g = shCoefficients(gr.C*1e-1 + CnR, gr.S*1e-1, GM = gr.GM, R = gr.R);
    fprintf('  (synthetic full field)\n');
end
T = g.subtractNormalField();                  % WGS84, auto-rescaled
newfig(vis, 'D03 disturbing geoid');
T.map(-89:89, 0:359, quantity = "geoid", Units = "m", ...
    Title = "geoid undulation (normal field subtracted)");
fprintf('  C20 before %.3e -> after %.3e\n', g.C(3, 1), T.C(3, 1));
end

function d04(vis, ~)
[g1, g2] = tryRealPair();
if isempty(g1)
    ts = demoSeries(40, 24, 14);
    g = ts.at(12) - ts.mean();
    fprintf('  (synthetic month)\n');
else
    g = g2 - g1;                              % real field, real stripes
    g = g.truncate(40);
    fprintf('  (real GRACE-FO minus GRACE difference)\n');
end
gG = g.gaussian(350);
gF = g.fan(350, 200);
gD = g.destripe();
newfig(vis, 'D04 filter differences (removed signal)');
subplot(1, 3, 1);
shLowLevel.plotSHCoeffTriangle(gG.C, gG.S, 'RefC', g.C, 'RefS', g.S, 'ax', gca);
title('Gaussian 350 - raw');
subplot(1, 3, 2);
shLowLevel.plotSHCoeffTriangle(gF.C, gF.S, 'RefC', g.C, 'RefS', g.S, 'ax', gca);
title('fan 350/200 - raw');
subplot(1, 3, 3);
shLowLevel.plotSHCoeffTriangle(gD.C, gD.S, 'RefC', g.C, 'RefS', g.S, 'ax', gca);
title('destripe - raw');
try
    W = shLowLevel.readDDK("DDK3", Nmax = 40);            % name resolution (v2.4.1)
    haveDDK = true;
catch
    haveDDK = false;
end
if haveDDK
    gK = g.applyDDK(W);
    newfig(vis, 'D04 DDK3 (real Wbd file)');
    shLowLevel.plotSHCoeffTriangle(gK.C, gK.S, 'RefC', g.C, 'RefS', g.S, 'ax', gca);
    title('DDK3 - raw');
end
end

function d05(vis, ~)
ts = tryRealSeries(30);
synth = isempty(ts);
if synth
    ts = demoSeries(30, 60, 15);
    fprintf('  (synthetic series; shLowLevel.fetchITSG(2010:2016) enables real data)\n');
else
    fprintf('  (real ITSG series, %d months)\n', ts.nEpochs);
end
cl = ts.climatology();
lat = -89:3:88; lon = 0:3:357;
newfig(vis, 'D05 climatology maps');
subplot(1, 2, 1);
shLowLevel.plotSHMap(shLowLevel.shSynthesis(cl.trendC, cl.trendS, ts.GM, ts.R, lat, lon), ...
    lat, lon, Units = "geoid m/yr", Title = "trend", ax = gca);
subplot(1, 2, 2);
amp = sqrt(shLowLevel.shSynthesis(cl.cosAnnC, cl.cosAnnS, ts.GM, ts.R, lat, lon).^2 ...
    + shLowLevel.shSynthesis(cl.sinAnnC, cl.sinAnnS, ts.GM, ts.R, lat, lon).^2);
shLowLevel.plotSHMap(amp, lat, lon, Units = "geoid m", Title = "annual amplitude", ...
    CLim = [0, shLowLevel.pctile(amp, 98)], Colormap = "parula", ax = gca);
% basin series: synthetic gets an artificial gap; real series carry
% their own (2017-2018 mission gap, dropouts)
idx = shLowLevel.shIndex(30, MinDegree = 0);
b = shLowLevel.basinKernel(idx, [-5 90; -5 130; 25 130; 25 90]);
if synth
    keep = [1:30, 41:60];                      % artificial gap 31-40
else
    keep = 1:ts.nEpochs;
end
tsg = ts.select(keep);
c = zeros(numel(keep), 1);
for k = 1:numel(keep)
    gk = tsg.at(k);
    c(k) = b' * shLowLevel.vecFromCS(gk.C, gk.S, idx) / (b' * b);
end
newfig(vis, 'D05 basin series');
shLowLevel.plotBasinSeries(tsg.epochs, c, 0.1 * std(c) * ones(size(c)), ...
    Units = "geoid m", Label = "test basin");
end

function d06(vis, ~)
ts = tryRealSeries(30);
if isempty(ts)
    ts = demoSeries(30, 48, 16);
    fprintf('  (synthetic series; shLowLevel.fetchITSG(2010:2016) enables real data)\n');
else
    ts = ts - ts.mean();
    fprintf('  (real ITSG series, %d months)\n', ts.nEpochs);
end
[tsF, op] = ts.filter("tvANS", Blocks = "auto");
idx = op.idx;
b = shLowLevel.basinKernel(idx, [-5 90; -5 130; 25 130; 25 90]);
[avgHat, out] = tsF.basinAverage(b, Deconvolve = true, Op = op);
newfig(vis, 'D06 tvANS basin series with posterior sigma');
shLowLevel.plotBasinSeries(ts.epochs(:), avgHat(:), out.sigma(:), ...
    Units = "geoid m", Label = "tvANS + deconvolution");
fprintf('  VCE factors: min %.2f max %.2f\n', min(op.s), max(op.s));
end

function d07(vis, ~)
g = demoField(24, 17);
[xg, ~] = shLowLevel.gaussLegendre(25);
lat = asind(xg(:)'); lon = (0:49) * 360 / 50;
grid = g.synthesis(lat, lon, UseCache = false);
g2 = shCoefficients.analysis(grid, lat, lon, 24);
newfig(vis, 'D07 analysis roundtrip error');
shLowLevel.plotSHCoeffTriangle(g2.C, g2.S, 'RefC', g.C, 'RefS', g.S, 'ax', gca);
title(sprintf('recovered - true (max %.1e)', ...
    max(abs(g2.C(:) - g.C(:)))));
fprintf('  ring roundtrip max error %.2e\n', max(abs(g2.C(:) - g.C(:))));
end

function d08(vis, ~)
L = 24; idx = shLowLevel.shIndex(L, MinDegree = 0);
[b, infoB] = shLowLevel.basinKernel(idx, [-5 90; -5 130; 25 130; 25 90], ...
    TaperKm = 200);
newfig(vis, 'D08 basin kernel (tapered)');
lat = -60:2:60; lon = 40:2:180;
[Cb, Sb] = shLowLevel.csFromVec(b, idx);
shLowLevel.plotSHMap(shLowLevel.shSynthesis(Cb, Sb, 1, 1, lat, lon), lat, lon, ...
    Title = sprintf('basin kernel, taper 200 km (area %.4f)', infoB.areaFraction));
% scaling factor on a Gaussian operator
ts = demoSeries(L, 12, 18);
Wn = shLowLevel.shGaussianWeights(L, 400);
W = diag(Wn(idx.n + 1));
k = shLowLevel.basinScaling(W, b, ts, idx = idx, tYears = ts.epochs');
fprintf('  Gaussian-400 basin scaling factor k = %.3f\n', k);
end

function d09(vis, ~)
rng(19);
idx = shLowLevel.shIndex(8, MinDegree = 2);
A = randn(idx.P) * 1e-10; M = A * A' + 1e-22 * eye(idx.P);
newfig(vis, 'D09 covariance & error map');
subplot(1, 2, 1);
shLowLevel.plotCovariance(M, idx, ax = gca);
subplot(1, 2, 2);
lat = -85:5:85; lon = 0:5:355;
sig = shLowLevel.errorMap(M, idx, lat, lon, quantity = "geoid");
shLowLevel.plotSHMap(sig, lat, lon, Units = "m", Title = "formal sigma (geoid)", ...
    CLim = [0, shLowLevel.pctile(sig, 98)], Colormap = "parula", ax = gca);
end

function d10(vis, ~)
L = 30; [kn, hn, ln] = synthLove(L);
ts = demoSeries(L, 36, 20);
tsr = ts - ts.mean();
sta = [47.1 8.6; -16.5 291.9; 64.1 338.0];     % demo stations
T = tsr.nEpochs;
up = zeros(3, T); no = up; ea = up;
for t = 1:T
    gk = tsr.at(t);
    [up(:, t), no(:, t), ea(:, t)] = gk.deformation(sta(:, 1)', sta(:, 2)', ...
        kn = kn, hn = hn, ln = ln, Mode = "points");
end
newfig(vis, 'D10 station deformation series');
plot(tsr.epochs, up * 1e3, '.-'); grid on;
xlabel('epoch [yr]'); ylabel('up [mm]');
legend("47N 9E", "17S 68W", "64N 22W");
title('elastic vertical deformation (synthetic Love numbers!)');
fprintf('  horizontal/vertical rms ratio: %.2f\n', ...
    sqrt(mean([no(:); ea(~isnan(ea))].^2)) / sqrt(mean(up(:).^2)));
end

function d11(vis, ~)
g = demoField(30, 21);
lat = -85:5:85; lon = 0:5:355;
[G, info] = shLowLevel.shSynthesisGradientTensor(g.C, g.S, g.GM, g.R, lat, lon, ...
    Height = 250e3);
newfig(vis, 'D11 gradient tensor at 250 km (Eotvos)');
f = ["uu" "nn" "ee" "un" "ue" "ne"];
for k = 1:6
    subplot(2, 3, k);
    shLowLevel.plotSHMap(G.(f(k)) * 1e9, lat, lon, Title = "G_{" + f(k) + "}", ...
        Coast = false, ax = gca);
end
fprintf('  Laplace trace residual: %.1e (built-in self-check)\n', ...
    info.maxTraceResidual);
end

function d12(vis, ~)
% synthetic BY DESIGN: EOF recovery is only checkable with known modes
rng(22);
L = 16; n1 = L + 1; T = 48;
P1 = tril(randn(n1)) * 1e-9; P2 = tril(randn(n1)) * 1e-9;
P2 = P2 - P1 * (P1(:)' * P2(:)) / (P1(:)' * P1(:));
a1 = 3 * sin(2*pi*(1:T)/16)'; a2 = randn(T, 1);
Cs = zeros(n1, n1, T);
for t = 1:T, Cs(:, :, t) = a1(t) * P1 + a2(t) * P2; end
ts = shSeries(Cs, Ss = zeros(n1, n1, T), Epochs = 2015 + (1:T)'/12);
[modes, pcs, ve] = shLowLevel.eofAnalysis(ts, NModes = 2);
newfig(vis, 'D12 EOF modes');
lat = -87:3:87; lon = 0:3:357;
subplot(2, 2, 1); modes{1}.map(lat, lon, Title = sprintf("mode 1 (%.0f%%)", 100*ve(1)));
subplot(2, 2, 2); modes{2}.map(lat, lon, Title = sprintf("mode 2 (%.0f%%)", 100*ve(2)));
subplot(2, 1, 2); plot(ts.epochs, pcs, '.-'); grid on;
legend('PC1', 'PC2'); xlabel('epoch [yr]'); title('principal components');
end

function d13(vis, ~)
% synthetic BY DESIGN: the F-test needs a known break as ground truth
rng(23);
L = 8; n1 = L + 1; T = 120;
ep = 2008 + (0:T-1)'/12; tb = 2013;
Cs = randn(n1, n1, T) * 2e-11;
hinge = max(ep - tb, 0);
for t = 1:T
    Cs(4, 2, t) = Cs(4, 2, t) + 1e-10 * hinge(t);    % C31 break
end
ts = shSeries(Cs, Ss = zeros(n1, n1, T), Epochs = ep);
out = ts.trendBreaks(Breaks = tb);
newfig(vis, 'D13 breakpoint F-test');
subplot(1, 2, 1);
imagesc(log10(max(out.F, 1e-2))); colorbar; axis image;
title('log_{10} F (C part)'); xlabel('order+1'); ylabel('degree+1');
subplot(1, 2, 2);
y = squeeze(Cs(4, 2, :));
plot(ep, y, '.'); hold on;
yfit = out.trend1C(4,2)*(ep - out.t0) + out.hingeC(4,2)*hinge + ...
    mean(y) - mean(out.trend1C(4,2)*(ep - out.t0) + out.hingeC(4,2)*hinge);
plot(ep, yfit, 'r-', 'LineWidth', 1.2); grid on;
title(sprintf('C_{31}: F = %.0f, p = %.1e', out.F(4,2), out.pValue(4,2)));
fprintf('  break at %.1f: F(C31) = %.1f, p = %.2e\n', tb, out.F(4, 2), out.pValue(4, 2));
end

function d14(vis, ~)
% synthetic BY DESIGN: VCE factor recovery needs known noise ratios
rng(24);
L = 10; n1 = L + 1; T = 18;
base = demoSeries(L, T, 24);
tsC = cell(1, 3); fac = [1, 2, 1.5];
for c = 1:3
    Cs = base.Cs; Ss = base.Ss;
    for t = 1:T
        Cs(:,:,t) = Cs(:,:,t) + tril(randn(n1)) * 1e-10 * fac(c);
        Ns = tril(randn(n1), -1) * 1e-10 * fac(c); Ns(:, 1) = 0;
        Ss(:,:,t) = Ss(:,:,t) + Ns;
    end
    tsC{c} = shSeries(Cs, Ss = Ss, Epochs = base.epochs);
end
[~, info] = shLowLevel.combineCenters(tsC);
newfig(vis, 'D14 multi-center VCE');
subplot(1, 2, 1);
plot(info.epochs, info.weights', '.-'); grid on;
legend('center 1', 'center 2', 'center 3');
xlabel('epoch [yr]'); ylabel('relative weight');
title('per-(center,month) VCE weights');
subplot(1, 2, 2);
imagesc(info.interCenterCorr); colorbar; clim([-1 1]); axis image;
title('inter-center residual correlation (honesty diagnostic)');
fprintf('  median factors: %s (true noise ratios 1 : 4 : 2.25)\n', ...
    mat2str(round(median(info.s2, 2)' / median(info.s2(1, :)), 2)));
end

function d15(vis, ~)
L = 32; [kn, hn] = synthLove(L);
idx = shLowLevel.shIndex(L, MinDegree = 0);
ocean = @(la, lo) ~(la > 55 & (lo < 60 | lo > 300));
loadF = @(la, lo) -100 * double(la > 65 & lo < 40);
[~, grid, info] = shLowLevel.seaLevelFingerprint(loadF, ocean, idx, ...
    kn = kn, hn = hn);
newfig(vis, 'D15 sea-level fingerprint');
Sn = info.S2D / info.eustatic; Sn(info.S2D == 0) = NaN;
shLowLevel.plotSHMap(Sn, grid.latDeg, grid.lonDeg, Units = "S / eustatic", ...
    Title = sprintf('fingerprint (%d iters, mass residual %.0e)', ...
    info.iterations, abs(info.massResidual)), CLim = [-1.5 1.5]);
fprintf('  eustatic %.4g m | conservation residual %.1e | converged: %d\n', ...
    info.eustatic, info.massResidual, info.converged);
end

function d16(~, outdir)
ts = demoSeries(16, 3, 26);
lat = (-88:4:88)'; lon = (0:4:356)';
G = zeros(numel(lat), numel(lon), 3);
for t = 1:3
    gk = ts.at(t);
    G(:, :, t) = shLowLevel.shSynthesis(gk.C, gk.S, gk.GM, gk.R, lat', lon');
end
ncf = fullfile(outdir, 'shx_demo_series.nc');
shLowLevel.writeGrid(ncf, G, lat, lon, Name = "geoid", Units = "m", ...
    Epochs = ts.epochs);
fprintf('  wrote %s\n', ncf);
mp4 = fullfile(outdir, 'shx_demo_anim.mp4');
shLowLevel.writeAnimation(ts, mp4, quantity = "geoid", lat = -85:5:85, ...
    lon = 0:5:355, FrameRate = 2, Units = "m");
fprintf('  wrote %s\n', mp4);
end
