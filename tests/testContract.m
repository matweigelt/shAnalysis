function tests = testContract
%TESTCONTRACT API contract tests for the shAnalysis v2 toolbox.
%
%   Verifies the *interface*, independent of numerics: every documented
%   error path throws its documented identifier, value objects are
%   immutable (methods return new objects, originals untouched), output
%   dimensions match the documented contracts, and the compat/ wrappers
%   remain signature-identical to shLowLevel internals.
%
%   Run:  results = runtests('testContract');
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

tests = functiontests(localfunctions);
end

% ---------------------------------------------------------------- fixture
function setupOnce(tc)
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);
tc.TestData.root = root;
tc.TestData.dataDir = fullfile(here, 'test_data');
shLowLevel.legendreCached('clear');
end

function g = randomField(nmax, epoch)
rng(42);
C = tril(randn(nmax+1)) * 1e-9; S = tril(randn(nmax+1), -1) * 1e-9;
C(1,1) = 1;
g = shCoefficients(C, S, GM = 3.986004415e14, R = 6.3781363e6, ...
    Epoch = epoch, ProductType = "GSM");
end

function ts = randomSeries(nmax, T)
arr = shCoefficients.empty(0, 1);
for k = 1:T
    arr(k,1) = randomField(nmax, 2020 + (k-1)/12);
end
ts = shSeries(arr);
end

% ----------------------------------------------------- shCoefficients IDs
function testConstructorBadInput(tc)
verifyError(tc, @() shCoefficients(ones(3,4), ones(3,4)), ...
    'shCoefficients:badInput');                      % non-square
verifyError(tc, @() shCoefficients(ones(3), ones(4)), ...
    'shCoefficients:badInput');                      % C/S size differ
verifyError(tc, @() shCoefficients(ones(3), ones(3), SigmaC = ones(4)), ...
    'shCoefficients:badInput');
verifyError(tc, @() shCoefficients(ones(3), ones(3), SigmaS = ones(4)), ...
    'shCoefficients:badInput');
end

function testTruncationAndIndexIDs(tc)
g = randomField(10, 2020);
verifyError(tc, @() g.truncate(11), 'shCoefficients:badTruncation');
verifyError(tc, @() g.setCoefficient(11, 0, 1, NaN), 'shCoefficients:badIndex');
verifyError(tc, @() g.setCoefficient(2, 3, 1, 1), 'shCoefficients:badIndex');
end

function testPairCheckIDs(tc)
a = randomField(10, 2020.0);
b = randomField(12, 2020.0);
verifyError(tc, @() a - b, 'shCoefficients:sizeMismatch');
c = shCoefficients(a.C, a.S, GM = a.GM * 1.001, R = a.R, Epoch = 2020);
verifyError(tc, @() a - c, 'shCoefficients:constantsMismatch');
d = randomField(10, 2021.0);        % plus enforces epoch agreement
verifyError(tc, @() a + d, 'shCoefficients:epochMismatch');
end

function testTN14IDs(tc)
tnFile = writeSyntheticTN14(tc);
g0 = randomField(10, NaN);
verifyError(tc, @() g0.applyTN14(tnFile), 'shCoefficients:noEpoch');
g1 = randomField(10, 2035.0);       % far outside table
verifyError(tc, @() g1.applyTN14(tnFile), 'shCoefficients:epochNotInTable');
g2 = randomField(10, 2020.13);      % nearest row has NaN C30
verifyError(tc, @() g2.applyTN14(tnFile, ReplaceC30 = "always"), ...
    'shCoefficients:noC30');
verifyError(tc, @() shLowLevel.readTN14('no_such_file.txt'), ...
    'shLowLevel:readTN14:fileNotFound');
end

function testCrossoverNoSigmasID(tc)
g = randomField(10, 2020);          % constructed without sigmas
verifyError(tc, @() g.crossover(), 'shCoefficients:noSigmas');
end

function testSynthesisIDs(tc)
g = randomField(10, 2020);
verifyError(tc, ...
    @() g.synthesis(-80:10:80, 0:10:350, quantity = "ewh"), ...
    'shSynthesis:missingLoveNumbers');
verifyError(tc, ...
    @() shLowLevel.shSynthesis(g.C, g.S, g.GM, g.R, 0, 0, 'quantity', 'nonsense'), ...
    'shSynthesis:badQuantity');
end

function testReadIDs(tc)
verifyError(tc, @() shCoefficients.read('no_such_file.gfc'), ...
    'shReadGFC:fileNotFound');
end

% ---------------------------------------------------------- shSeries IDs
function testSeriesConstructionIDs(tc)
verifyError(tc, @() shSeries(rand(3,3,4)), 'shSeries:badInput');   % no Epochs
verifyError(tc, @() shSeries(rand(3,3,4), Ss = rand(3,3,4), ...
    Epochs = [1 2 3]), 'shSeries:badInput');                       % T mismatch
a = randomField(8, 2020.0); b = randomField(10, 2020.1);
verifyError(tc, @() shSeries([a; b]), 'shSeries:sizeMismatch');
c = shCoefficients(a.C, a.S, GM = a.GM, R = a.R * 1.01, Epoch = 2020.1);
verifyError(tc, @() shSeries([a; c]), 'shSeries:constantsMismatch');
end

function testSeriesMethodIDs(tc)
ts = randomSeries(8, 8);
verifyError(tc, @() ts.at(9), 'shSeries:badIndex');
verifyError(tc, @() ts.basinAverage(ones(3, 1)), 'shSeries:badBasis');
verifyError(tc, @() ts.basinAverage(ones(3, 1), Deconvolve = true), ...
    'shSeries:missingOp');
verifyError(tc, @() shSeries.read('zz_no_such_*.gfc'), 'shSeries:noFiles');

short = randomSeries(8, 3);
verifyError(tc, @() short.climatology(), 'shSeries:tooFewEpochs');

nanEp = shSeries(ts.Cs, Ss = ts.Ss, Epochs = [NaN, 2:8], ...
    GM = ts.GM, R = ts.R);
verifyError(tc, @() nanEp.climatology(), 'shSeries:noEpoch');
verifyError(tc, @() nanEp.filter("tvANS"), 'shSeries:noEpoch');

Cs = ts.Cs; Cs(2,1,3) = NaN;
nanTs = shSeries(Cs, Ss = ts.Ss, Epochs = ts.epochs, GM = ts.GM, R = ts.R);
verifyError(tc, @() nanTs.destripe(), 'shSeries:nanInSeries');
end

function testRestoreIDs(tc)
gsm = randomSeries(8, 6);
gad = randomSeries(10, 6);
verifyError(tc, @() gsm.restore(gad), 'shSeries:sizeMismatch');
gad8 = gsm.truncate(8);             % same nmax, but shift epochs far away
shifted = shSeries(gad8.Cs, Ss = gad8.Ss, Epochs = gad8.epochs + 5, ...
    GM = gad8.GM, R = gad8.R, ProductType = "GAD");
verifyError(tc, @() gsm.restore(shifted), 'shSeries:epochMismatch');
end

% ----------------------------------------------------------- immutability
function testImmutability(tc)
g = randomField(20, 2020);
C0 = g.C; S0 = g.S; h0 = g.history;
g2 = g.destripe(minOrder = 4); %#ok<NASGU>
g3 = g.gaussian(300);          %#ok<NASGU>
g4 = g.truncate(10);           %#ok<NASGU>
verifyEqual(tc, g.C, C0); verifyEqual(tc, g.S, S0);
verifyEqual(tc, g.history, h0);

ts = randomSeries(8, 6);
Cs0 = ts.Cs;
ts2 = ts.gaussian(500); %#ok<NASGU>
verifyEqual(tc, ts.Cs, Cs0);
end

function testFluentHistory(tc)
g = randomField(20, 2020);
g2 = g.destripe().gaussian(300).truncate(15);
verifyClass(tc, g2, 'shCoefficients');
verifyEqual(tc, numel(g2.history), numel(g.history) + 3);
verifySubstring(tc, char(g2.history(end)), 'truncate');
end

% ------------------------------------------------------------- dimensions
function testOutputDimensions(tc)
g = randomField(15, 2020);
lat = -60:5:60; lon = 0:5:355;
[grid, latO, lonO] = g.synthesis(lat, lon);
verifySize(tc, grid, [numel(lat), numel(lon)]);
verifyEqual(tc, latO(:)', lat); verifyEqual(tc, lonO(:)', lon);

spec = g.degreeRMS;
n = (spec.degree(1):spec.degree(end))';
verifyEqual(tc, numel(spec.degRMS), numel(n));
verifyEqual(tc, numel(spec.degAmplitude), numel(n));

idx = shLowLevel.shIndex(g.nmax, MinDegree = 2);
x = g.vec(idx);
verifySize(tc, x, [idx.P, 1]);
g5 = shCoefficients.fromVec(x, idx, g);
verifyEqual(tc, g5.C(3:end, :), g.C(3:end, :), AbsTol = 0);

ts = randomSeries(8, 10);
[clim, resid] = ts.climatology();
verifyClass(tc, clim, 'shClimatology');
verifyEqual(tc, resid.nEpochs, 10);
m = ts.mean;
verifyClass(tc, m, 'shCoefficients');
verifySize(tc, m.C, size(ts.Cs(:, :, 1)));
end

% ------------------------------------------------- compat wrapper contract
% ----------------------------------------------------------------- helper
function f = writeSyntheticTN14(tc)
% Minimal TN-14-like table: 2 months; second month has no C30 (NaN).
f = fullfile(tempdir, 'tn14_contract.txt');
fid = fopen(f, 'w');
fprintf(fid, 'Product: synthetic TN-14 for contract tests\n');
fprintf(fid, 'end of header ===========================\n');
% MJD start, frac yr start, C20, dC20*1e10, sigC20*1e10, C30, dC30*1e10,
% sigC30*1e10, MJD end, frac yr end   (10-column TN-14 layout)
% -> mid-epochs 2020.0417 and 2020.1250
fprintf(fid, ['58849.0 2020.0000 -4.841694e-04 0.5 0.3 ' ...
    '-9.571e-07 0.4 0.2 58880.0 2020.0833\n']);
fprintf(fid, ['58880.0 2020.0833 -4.841700e-04 0.6 0.3 ' ...
    'NaN NaN NaN 58909.0 2020.1666\n']);
fclose(fid);
tc.addTeardown(@() delete(f));
end

% =================================================================== v2.1
function testAnalysisBadInputID(testCase)
verifyError(testCase, @() shLowLevel.shAnalysisGrid(ones(3,4), 1:5, 1:4, 2), ...
    'shLowLevel:shAnalysisGrid:badInput');
end

function testAnalysisBadGridID(testCase)
% non-uniform longitude must be rejected by the rings method
g = ones(8, 10);
lat = linspace(-60, 60, 8);
lon = [0 30 70 100 140 180 220 260 300 340];
verifyError(testCase, ...
    @() shLowLevel.shAnalysisGrid(g, lat, lon, 2, Method = "rings"), ...
    'shLowLevel:shAnalysisGrid:badGrid');
end

function testAnalysisAliasedGridID(testCase)
% nlon <= 2*nmax aliases orders
g = ones(20, 10); lat = linspace(-80, 80, 20); lon = (0:9)*36;
verifyError(testCase, @() shLowLevel.shAnalysisGrid(g, lat, lon, 6, Method = "rings"), ...
    'shLowLevel:shAnalysisGrid:badGrid');
end

function testSynthesisBadMethodID(testCase)
C = zeros(3); S = zeros(3);
verifyError(testCase, @() shLowLevel.shSynthesis(C, S, 1, 1, 0, [0 10 50], ...
    'method', 'fft'), 'shSynthesis:badMethod');
end

function testReadTN13ErrorIDs(testCase)
verifyError(testCase, @() shLowLevel.readTN13('/nonexistent/tn13.txt'), ...
    'shLowLevel:readTN13:fileNotFound');
f = [tempname '.txt'];
fid = fopen(f, 'w'); fprintf(fid, 'header only\nno data here\n'); fclose(fid);
cleanup = onCleanup(@() delete(f)); %#ok<NASGU>
verifyError(testCase, @() shLowLevel.readTN13(f), 'shLowLevel:readTN13:noData');
end

function testReadSINEXErrorIDs(testCase)
verifyError(testCase, @() shLowLevel.readSINEX('/nonexistent/x.snx'), ...
    'shLowLevel:readSINEX:fileNotFound');
f = [tempname '.snx'];
fid = fopen(f, 'w'); fprintf(fid, '%%=SNX 2.02\n+SOME/BLOCK\n-SOME/BLOCK\n%%ENDSNX\n'); fclose(fid);
cleanup = onCleanup(@() delete(f)); %#ok<NASGU>
verifyError(testCase, @() shLowLevel.readSINEX(f), 'shLowLevel:readSINEX:noData');
end

function testEvalGFCTEpochOutsideID(testCase)
f = [tempname '.gfc'];
fid = fopen(f, 'w');
fprintf(fid, ['product_type gravity_field\nmodelname T2\n' ...
    'earth_gravity_constant 3.986004415E+14\nradius 6378136.3\n' ...
    'max_degree 2\nformat icgem2.0\nend_of_head ====\n' ...
    'gfct 2 0 -4.84e-4 0.0 0 0 20100101 20150101\n']);
fclose(fid);
cleanup = onCleanup(@() delete(f)); %#ok<NASGU>
model = shLowLevel.shReadGFC(f);
verifyError(testCase, @() shLowLevel.shEvalGFCT(model, 2020.0), ...
    'shEvalGFCT:epochOutside');
end

function testTvANSBlocksUnavailableID(testCase)
idx = shLowLevel.shIndex(4);
X = randn(idx.P, 12); t = 2002 + (0:11)'/12;
verifyError(testCase, @() shLowLevel.tvANSFilter(X, t, idx, Blocks = 'on', ...
    Constraints = ones(idx.P, 1)), 'shLowLevel:tvANSFilter:blocksUnavailable');
end

function testNoiseCovNotBlockDiagonalID(testCase)
% external NoiseCov on Blocks='on' must be block-diagonal in the
% (order, C/S, parity) partition; Blocks='auto' falls back quietly
rng(30);
idx = shLowLevel.shIndex(4);
X = randn(idx.P, 12); t = 2002 + (0:11)'/12;
Nfull = eye(idx.P) + 0.3 * ones(idx.P);          % dense: NOT block-diag
verifyError(testCase, @() shLowLevel.tvANSFilter(X, t, idx, Blocks = 'on', ...
    NoiseCov = Nfull), 'shLowLevel:buildNoiseCov:notBlockDiagonal');
[Xf, op] = shLowLevel.tvANSFilter(X, t, idx, Blocks = 'auto', NoiseCov = Nfull);
verifyTrue(testCase, all(isfinite(Xf(:))));
verifyEqual(testCase, op.layout, 'full');        % fell back to full path
% a DIAGONAL external N is block-diagonal under any partition: block
% path engages
[~, op2] = shLowLevel.tvANSFilter(X, t, idx, Blocks = 'on', NoiseCov = eye(idx.P));
verifyEqual(testCase, op2.layout, 'blocks');
end

function testSelectBadInputID(testCase)
n1 = 3; T = 5;
ts = shSeries(zeros(n1,n1,T), Ss = zeros(n1,n1,T), Epochs = 2002+(0:T-1)');
verifyError(testCase, @() ts.select([0 7 9]), 'shSeries:badInput');
end

function testClimatologyPeriodicBadIndexID(testCase)
T = 24; n1 = 3;
ts = shSeries(1e-9*randn(n1,n1,T), Ss = zeros(n1,n1,T), ...
    Epochs = 2002 + (0:T-1)'/12);
clim = ts.climatology();
verifyError(testCase, @() clim.periodic(1), 'shClimatology:badIndex');
end

function testAddDegree1NoEpochID(testCase)
g = shCoefficients(zeros(3), zeros(3));
tn = struct('epoch', 2008.29, 'C10', 1e-10, 'C11', 2e-10, 'S11', 3e-10, ...
    'sigC10', 1e-12, 'sigC11', 1e-12, 'sigS11', 1e-12);
verifyError(testCase, @() g.addDegree1(tn), 'shCoefficients:noEpoch');
end

% =============================================================== v2.4
function testTriangleDiffModeFullSineWing(testCase)
% Regression (v2.4.2): the RefC/RefS difference branch masked S with
% tril(...,-1), blanking the sine SECTORALS S_nn and rendering the left
% wing one column narrower than the C wing in every diff triangle
% (demo D04/D07 figures). Pin the rendered CData: every valid sine cell
% (n >= m >= 1, sectorals included) must be finite; m=0 must be NaN on
% the sine side; the C wing must span m = 0..n.
fig = figure('Visible', 'off');
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
L = 8; n1 = L + 1;
C = tril(ones(n1));                              % all valid cells nonzero
S = tril(ones(n1)); S(:, 1) = 0;                 % incl. sectorals S_nn
ax = axes('Parent', fig);
shLowLevel.plotSHCoeffTriangle(C, S, 'RefC', 0*C, 'RefS', 0*S, 'ax', ax);
im = findobj(ax, 'Type', 'image');
verifyNotEmpty(testCase, im, 'no image object rendered');
img = im(1).CData;                               % (L+1) x (2L+1)
verifySize(testCase, img, [n1, 2*L + 1]);
for n = 0:L
    for m = 1:n                                  % sine wing at x = -m
        verifyTrue(testCase, isfinite(img(n+1, L+1-m)), sprintf( ...
            'sine cell S(%d,%d) missing from diff triangle', n, m));
    end
    for m = 0:n                                  % cosine wing at x = +m
        verifyTrue(testCase, isfinite(img(n+1, L+1+m)), sprintf( ...
            'cosine cell C(%d,%d) missing from diff triangle', n, m));
    end
    if n < L                                     % outside the triangle
        verifyTrue(testCase, all(isnan(img(n+1, [1:L-n, L+2+n:end]))), ...
            sprintf('degree %d shows cells beyond |m| = n', n));
    end
end
end

function testPlotSuiteRuns(testCase)
rng(81);
fig = figure('Visible', 'off');
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
L = 12; n1 = L + 1;
C = tril(randn(n1)) * 1e-8; S = tril(randn(n1), -1) * 1e-8; S(:, 1) = 0;
grid = shLowLevel.shSynthesis(C, S, 3.986004415e14, 6378136.3, ...
    -85:5:85, 0:5:355, 'quantity', 'geoid');
ax1 = axes('Parent', fig);
h = shLowLevel.plotSHMap(grid, -85:5:85, 0:5:355, ax = ax1, Units = "m");
verifyClass(testCase, h, 'matlab.graphics.axis.Axes');
cla(ax1);
h2 = shLowLevel.plotSHMap(grid, -85:5:85, 0:5:355, ax = ax1, ...
    Projection = "hammer");
verifyClass(testCase, h2, 'matlab.graphics.axis.Axes');
% basin series with gap + band
t = [2010 + (0:40)/12, 2015 + (0:30)/12];
c = sin(2*pi*t(:)) + 0.01*(t(:) - 2012);
cla(ax1);
h3 = shLowLevel.plotBasinSeries(t(:), c, 0.2 + 0*c, ax = ax1, Units = "cm");
verifyClass(testCase, h3, 'matlab.graphics.axis.Axes');
% covariance plot
idx = shLowLevel.shIndex(6, MinDegree = 2);
A = randn(idx.P); M = A*A';
cla(ax1);
h4 = shLowLevel.plotCovariance(M, idx, ax = ax1);
verifyClass(testCase, h4, 'matlab.graphics.axis.Axes');
% triangle diff mode
cla(ax1);
h5 = shLowLevel.plotSHCoeffTriangle(C, S, 'RefC', 0.9*C, 'RefS', 0.9*S, 'ax', ax1);
verifyClass(testCase, h5, 'matlab.graphics.axis.Axes');
verifyError(testCase, ...
    @() shLowLevel.plotSHCoeffTriangle(C, S, 'RefC', C, 'ax', axes('Parent', fig)), ...
    'shLowLevel:plotSHCoeffTriangle:needRefS');
end

function testSpectrumOverlays(testCase)
rng(82);
fig = figure('Visible', 'off');
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
L = 30; n1 = L + 1;
C = tril(randn(n1)) * 1e-9 ./ max(1, (0:L)').^2;
S = tril(randn(n1), -1) * 1e-9 ./ max(1, (0:L)').^2; S(:, 1) = 0;
sC = 1e-11 * max(1, (0:L)') / 10 .* ones(n1) .* tril(ones(n1));
spec = shLowLevel.shDegreeRMS(C, S, 'R', 6378136.3, ...
    'sigmaC', sC, 'sigmaS', sC);
ax = axes('Parent', fig);
h = shLowLevel.plotSHSpectrum(spec, 'ax', ax, 'Kaula', 1e-5, 'MarkCrossover', true);
verifyClass(testCase, h(1), 'matlab.graphics.chart.primitive.Line');
end

function testWriteAnimationSmoke(testCase)
rng(83);
L = 8; n1 = L + 1; T = 2;
Cs = randn(n1, n1, T) * 1e-9; Ss = randn(n1, n1, T) * 1e-9;
for t = 1:T
    Cs(:, :, t) = tril(Cs(:, :, t));
    Ss(:, :, t) = tril(Ss(:, :, t), -1); Ss(:, 1, t) = 0;
end
ts = shSeries(Cs, Ss = Ss, Epochs = [2010.1 2010.2]);
tmp = fullfile(tempdir, sprintf('shx_anim_%d.avi', randi(1e9)));  % CI-portable
cleanup = onCleanup(@() deleteIfThere(tmp)); %#ok<NASGU>
shLowLevel.writeAnimation(ts, tmp, quantity = "geoid", ...
    lat = -80:10:80, lon = 0:10:350, FrameRate = 2);
verifyTrue(testCase, isfile(tmp));
d = dir(tmp);
verifyGreaterThan(testCase, d.bytes, 1000);
verifyTrue(testCase, isfile([tmp '.provenance.json']));   % v2.7.0
delete([tmp '.provenance.json']);
end

function deleteIfThere(f)
if isfile(f), delete(f); end
end

function testBasinScalingRecoversFactor(testCase)
rng(84);
% operator: pure Gaussian attenuation as a tvANS-like op via matFilter
L = 20; idx = shLowLevel.shIndex(L, MinDegree = 0);
Wn = shLowLevel.shGaussianWeights(L, 500);
W = diag(Wn(idx.n + 1));
T = 10; ep = 2010 + (0:T-1)/12;
% static matrix operator (Gaussian in idx ordering)
% model series: fixed pattern, varying amplitude
n1 = L + 1;
pat = tril(randn(n1)) * 1e-9; patS = tril(randn(n1), -1) * 1e-9;
patS(:, 1) = 0;
amp = 1 + 0.5 * sin(2*pi*(1:T)/T);
Cs = zeros(n1, n1, T); Ss = Cs;
for t = 1:T
    Cs(:, :, t) = amp(t) * pat; Ss(:, :, t) = amp(t) * patS;
end
tsM = shSeries(Cs, Ss = Ss, Epochs = ep);
b = shLowLevel.basinKernel(idx, [ -10 100; -10 140; 25 140; 25 100 ]);
[k, info] = shLowLevel.basinScaling(W, b, tsM, idx = idx, tYears = ep);
% k must invert the operator's basin-average attenuation: with a FIXED
% pattern, aTrue = k * aFilt exactly, so sigmaK ~ 0
x1 = shLowLevel.vecFromCS(pat, patS, idx);
kExact = (b' * x1) / (b' * (W * x1));
verifyEqual(testCase, k, kExact, 'RelTol', 1e-10);
verifyLessThan(testCase, info.sigmaK, 1e-8 * abs(k));
verifyEqual(testCase, info.nMatched, T);
verifyError(testCase, ...
    @() shLowLevel.basinScaling(W, b(1:3), tsM, idx = idx, tYears = ep), ...
    'shLowLevel:basinScaling:badKernel');
end

function testDemoRegistryAndSmoke(testCase)
reg = demo_shAnalysis("list");
verifyEqual(testCase, numel(reg), 16);
verifyEqual(testCase, numel(unique([reg.id])), 16);
verifyTrue(testCase, all(arrayfun(@(r) isa(r.run, 'function_handle'), reg)));
% two cheap cases run headless without touching the screen; a failing
% demo raises demo:casesFailed, which is a hard test failure here
% (D01 failing via an unqualified compat-era plot call slipped through
% this smoke test once - never again)
verifyWarningFree(testCase, ...
    @() demo_shAnalysis(["D01", "D13"], Visible = false));
close all hidden;
verifyError(testCase, @() demo_shAnalysis("D99"), 'demo:unknownCase');
% fail-and-continue contract: D16 into an impossible target (a FILE
% posing as OutDir) must not abort the run (warning instead);
% StopOnError=true must rethrow. writeGrid itself creates missing
% folders since v2.4.1, so a plain nonexistent dir no longer fails.
badDir = [tempname, '.txt'];
fid = fopen(badDir, 'w'); fclose(fid);
cleanupBad = onCleanup(@() delete(badDir)); %#ok<NASGU>
fprintf(['\n  NOTE: the following D16 "FAILED" line is INTENTIONAL - ' ...
    'it exercises the fail-and-continue contract.\n']);
verifyWarning(testCase, ...
    @() demo_shAnalysis("D16", Visible = false, OutDir = badDir), ...
    'demo:casesFailed');
close all hidden;
verifyError(testCase, ...
    @() demo_shAnalysis("D16", Visible = false, OutDir = badDir, ...
    StopOnError = true), ?MException);
end

function testAPIExamplesRun(testCase)
% v2.5: every runnable example of the API reference (docs/
% apiExamples.json - the SAME file the workflow-guide Part IV is built
% from) must execute against the canonical fixture workspace. This is
% the executable link between documentation and contract: an API change
% that breaks a documented example fails runAllTests.
root = fileparts(fileparts(mfilename('fullpath')));
J = jsondecode(fileread(fullfile(root, 'docs', 'apiExamples.json')));
d = testCase.TestData.dataDir;
tmp = tempname; mkdir(tmp);
old = cd(tmp);
cleanup = onCleanup(@() cd(old));
cleanup2 = onCleanup(@() rmIfFolder(tmp)); %#ok<NASGU>
% bare-filename fixtures into the cwd
for f = ["ITSG-Grace2018_n60_2008-04.gfc", ...
         "ITSG-Grace_operational_n60_2025-12.gfc", ...
         "TN-13_GEOC_CSR_RL06.3.txt", "TN-14_C30_C20_SLR_GSFC.txt", ...
         "test_variable.gfct", ...
         "ITSG-Grace2018_n96_2008-04_head12.snx"]
    copyfile(fullfile(d, f), tmp);
end
ws = buildAPIWorkspace(tmp);
set(groot, 'defaultFigureVisible', 'off');
cleanup3 = onCleanup(@() set(groot, 'defaultFigureVisible', 'on')); %#ok<NASGU>
failures = strings(1, 0);
for k = 1:numel(J)
    if ~J(k).run, continue; end
    try
        runAPISnippet(J(k).example, ws);
    catch err
        failures(end+1) = J(k).name + ": " + err.message; %#ok<AGROW>
    end
    close all
end
verifyEmpty(testCase, failures, "API examples failed:" + newline + ...
    strjoin(failures, newline));
end

function ws = buildAPIWorkspace(tmp) %#ok<INUSD>
% canonical fixture workspace for the API examples (real data chain)
ws = struct();
ws.g = shCoefficients.read("ITSG-Grace2018_n60_2008-04.gfc", ...
    Epoch = 2008 + 3.5/12);
ws.gF = ws.g.gaussian(350);
ws.g10 = ws.g.truncate(10);
ws.idx = shLowLevel.shIndex(10);
ws.x = shLowLevel.vecFromCS(ws.g10.C, ws.g10.S, ws.idx);
ws.GM = ws.g.GM; ws.R = ws.g.R;
ws.GMref = 3.986004415e14; ws.Rref = 6.3781363e6;
nn = (0:120)';
ws.kn = -0.3 * nn ./ (nn + 6);
% Love-number table via the reader (real-life G7 path)
lnRows = compose("%d %.6f %.6f %.6f", nn, -0.9 * nn ./ (nn + 4), ...
    0.1 * nn ./ (nn + 8), -0.3 * nn ./ (nn + 6));
writelines(["# n h l k"; lnRows], "ln_table.txt");
ws.LN = shLowLevel.readLoveNumbers("ln_table.txt", MaxDegree = 120);
copyfile("ln_table.txt", "prem_load.txt");   % readLoveNumbers help example
ws.lat = -88:8:88; ws.lon = 0:12:348;    % = the grid axes below
ws.latPts = [10 20 30]; ws.lonPts = [40 50 60];
ws.latGc = [10 20 30]; ws.latGd = [10.06 20.10 30.12];
ws.theta = deg2rad(30:15:60);
ws.grid = ws.g10.synthesis(-88:8:88, 0:12:348);
ws.tn13 = shLowLevel.readTN13("TN-13_GEOC_CSR_RL06.3.txt");
ws.tn14 = shLowLevel.readTN14("TN-14_C30_C20_SLR_GSFC.txt");
ws.model = shLowLevel.shReadGFC("test_variable.gfct");   % struct incl. variableTerms
% synthetic monthly series + tvANS chain
rng(31);
L = 10; T = 48; n1 = L + 1;              % same L as ws.idx = shIndex(10)
tY = 2019 + (0:T-1)'/12;                 % continuous TN-13/14 coverage
mL = tril(true(n1)); mL1 = mL; mL1(:, 1) = false;
Cs = zeros(n1, n1, T); Ss = Cs;
for t = 1:T
    Cs(:,:,t) = (mL .* randn(n1)) * 1e-9 + (mL * 2e-9) * cos(2*pi*tY(t));
    Ss(:,:,t) = (mL1 .* randn(n1)) * 1e-9;
end
ws.ts = shSeries(Cs, Ss = Ss, Epochs = tY);
ws.Cs = Cs; ws.Ss = Ss;                  % shSeries constructor example
ws.C = ws.g10.C; ws.S = ws.g10.S;        % shCoefficients constructor example
ws.gad = ws.g;                           % plus: epoch-matched background
ws.gMean = ws.g;                         % minus: any compatible field
ws.w = ones(size(ws.g.C));               % times: per-coefficient weights
ws.c20slr = ws.tn14.C20(1);              % setCoefficient
ws.tsGAD = ws.ts;                        % restore: epoch-matched series
[ws.clim, ws.resid] = ws.ts.climatology(Periods = 161/365.25);
idxF = shLowLevel.shIndex(L);
ws.Ac = exp(-idxF.n / 3) .* randn(idxF.P, 1);
[ws.tsF, ws.op, ws.info] = ws.ts.filter("tvANS", Blocks = "off");
ws.B = randn(idxF.P, 2);
ws.b = ws.B(:, 1);
[ws.avg, ws.out] = shLowLevel.basinDeconvolve(ws.B, ws.op);
ws.Y3 = [ws.avg(1, :)', ws.avg(1, :)' + 2e-4 * randn(T, 1), ...
    ws.avg(1, :)' + 5e-4 * randn(T, 1)];   % 3-center stack for TCH
ws.tYears = tY;
XV = zeros(ws.idx.P, T);                 % idx-ordered raw stack (L = 10)
for t = 1:T
    XV(:, t) = shLowLevel.vecFromCS(Cs(:,:,t), Ss(:,:,t), ws.idx);
end
ws.X = XV;
[~, ws.Xres] = shLowLevel.fitDeterministicModel(XV, tY);   % idx-ordered residuals
ws.N = shLowLevel.buildNoiseCov(ws.Xres, ws.idx);
A0 = randn(ws.idx.P);
ws.M = (A0 * A0' / ws.idx.P) * 1e-20;   % SPD "covariance" for errorMap
% three 'centers' for combineCenters: identical series suffice for
% an example run (shSeries has no plus by design - restore covers it)
ws.tsCSR = ws.ts; ws.tsGFZ = ws.ts; ws.tsJPL = ws.ts;
ws.giaModel = 1e-10 * ws.g10;
ws.gLoad = ws.g10; ws.oceanMask = [];   % fingerprint example not run
ws.res = randn(200, 1);
ws.coef = []; ws.t0 = 2016.0;           % fromCoef example not run
ws.gT10 = ws.g10;
end

function runAPISnippet(code, ws)
% evaluate one documented example inside the fixture workspace
fn = fieldnames(ws);
for kk = 1:numel(fn)
    eval(sprintf('%s = ws.%s;', fn{kk}, fn{kk}));  %#ok<EVLDOT>
end
eval(code);  %#ok<EVLC>
end

function testSetupDryRunPlan(testCase)
% setup_shAnalysis(DryRun=true) must have ZERO side effects and return
% the full cumulative plan for the requested level.
p0 = path;
restore = onCleanup(@() path(p0)); %#ok<NASGU>
s = setup_shAnalysis(DryRun = true, Download = "starter", ...
    DDK = [2 5], Permanent = true, Docs = true, Quiet = true);
verifyEqual(testCase, path, p0);                 % path untouched
verifyEqual(testCase, s.pathAction, "planned");
txt = strjoin(s.plan, newline);
verifyTrue(testCase, contains(txt, "TN-14"));
verifyTrue(testCase, contains(txt, "TN-13"));
verifyTrue(testCase, contains(txt, "DDK2, DDK5"));
verifyTrue(testCase, contains(txt, "ITSG n96"));
verifyTrue(testCase, contains(txt, "2008-04, 2025-12"));
verifyTrue(testCase, contains(txt, "builddocsearchdb"));
verifyTrue(testCase, s.ok && isempty(s.fetched) && isempty(s.failed));
% cumulative levels: "core" plans TN but no DDK/ITSG lines
s2 = setup_shAnalysis(DryRun = true, Download = "core", Quiet = true);
t2 = strjoin(s2.plan, newline);
verifyTrue(testCase, contains(t2, "TN-13") && ~contains(t2, "DDK") ...
    && ~contains(t2, "ITSG"));
% default: path line only
s3 = setup_shAnalysis(DryRun = true, Quiet = true);
verifyEqual(testCase, numel(s3.plan), 1);
verifyTrue(testCase, startsWith(s3.plan(1), "path:"));
% Update=true is announced in the plan lines
s4 = setup_shAnalysis(DryRun = true, Download = "core", ...
    Update = true, Quiet = true);
verifyTrue(testCase, contains(strjoin(s4.plan, newline), ...
    "(update existing)"));
end

function testComparePlotsSmoke(testCase)
% v2.6.0 visualization smoke: taylorDiagram and the 4-panel figures of
% both compare aggregators build and close without error
old = get(0, 'DefaultFigureVisible');
cl = onCleanup(@() set(0, 'DefaultFigureVisible', old)); %#ok<NASGU>
set(0, 'DefaultFigureVisible', 'off');
figure;
h = shLowLevel.taylorDiagram(1.0, [0.8, 1.1], [0.9, 0.95], ...
    Labels = ["a", "b"], Normalize = true);
verifyTrue(testCase, isgraphics(h));
d = fullfile(fileparts(mfilename('fullpath')), 'test_data');
fG = fullfile(d, 'ITSG-Grace2018_n60_2008-04.gfc');
assumeTrue(testCase, isfile(fG));
g = shCoefficients.read(fG, Epoch = 2008.29);
[~, h1] = shLowLevel.compareSolutions(g, g.gaussian(500), Plot = true, ...
    LatDeg = -85:5:85, LonDeg = 0:9:351);
verifyTrue(testCase, isgraphics(h1));
rng(11); n1 = 7; T = 24;
Cs = 1e-9 * randn(n1, n1, T); Ss = 1e-9 * randn(n1, n1, T);
mk = @(s) shSeries(Cs + s * 1e-11 * randn(n1, n1, T), ...
    Ss = Ss + s * 1e-11 * randn(n1, n1, T), Epochs = 2019 + (0:T-1)'/12);
[~, h2] = shLowLevel.compareSeries({mk(1), mk(2), mk(3)}, Plot = true, ...
    LatDeg = -80:20:80, LonDeg = 0:30:330);
verifyTrue(testCase, isgraphics(h2));
close all
end

function testProvenanceSidecars(testCase)
% v2.7.0: writers emit <file>.provenance.json (Sidecar=false disables)
% and writeGrid netCDF is CF-1.8 complete with a dynamic source stamp
rng(22); L = 6;
C = 1e-9 * randn(L+1); S = 1e-9 * randn(L+1); S(:, 1) = 0;
tmp = tempname; mkdir(tmp);
cl = onCleanup(@() rmIfFolder(tmp)); %#ok<NASGU>
fg = fullfile(tmp, 'out.gfc');
shLowLevel.writeGFC(fg, C, S, 3.986004415e14, 6378136.3);
sj = [fg '.provenance.json'];
verifyTrue(testCase, isfile(sj));
p = jsondecode(fileread(sj));
verifyTrue(testCase, contains(p.tool, "shAnalysis"));
verifyEqual(testCase, p.nmax, L);
verifyTrue(testCase, isfield(p, 'matlab') && isfield(p, 'created'));
% opt-out
fg2 = fullfile(tmp, 'out2.gfc');
shLowLevel.writeGFC(fg2, C, S, 3.986004415e14, 6378136.3, Sidecar = false);
verifyFalse(testCase, isfile([fg2 '.provenance.json']));
% netCDF: CF-1.8 attributes + sidecar
lat = -80:20:80; lon = 0:30:330;
G = randn(numel(lat), numel(lon));
fn = fullfile(tmp, 'grid.nc');
shLowLevel.writeGrid(fn, G, lat, lon, Name = "ewh", Units = "m");
verifyEqual(testCase, ncreadatt(fn, '/', 'Conventions'), 'CF-1.8');
verifyEqual(testCase, ncreadatt(fn, 'lat', 'standard_name'), 'latitude');
verifyEqual(testCase, ncreadatt(fn, 'lon', 'axis'), 'X');
verifyEqual(testCase, ncreadatt(fn, 'ewh', 'coordinates'), 'lat lon');
src = ncreadatt(fn, '/', 'source');
v = shLowLevel.version();
verifyTrue(testCase, contains(src, char(v.Version)));   % no stale stamps
verifyTrue(testCase, isfile([fn '.provenance.json']));
end

function testStandardChain(testCase)
% v3.1.0: the canonical pipeline entry point on fixtures - order,
% report, GIA subtraction, all three filter forms, error contracts
d = fullfile(fileparts(mfilename('fullpath')), 'test_data');
src = fullfile(d, 'ITSG-Grace2018_n60_2008-04.gfc');
tn14 = fullfile(d, 'TN-14_C30_C20_SLR_GSFC.txt');
% the shipped GFZ fixture spells the release with an UNDERSCORE
% (RL06_3), unlike CSR and JPL. A dot here made assumeTrue filter this
% whole test out silently for four releases - a skipped test looks
% exactly like a green one in the summary, so verify the fixtures exist
% rather than assuming them away.
tn13 = fullfile(d, 'TN-13_GEOC_GFZ_RL06_3.txt');
verifyTrue(testCase, isfile(src), 'missing fixture: ' + string(src));
verifyTrue(testCase, isfile(tn14), 'missing fixture: ' + string(tn14));
verifyTrue(testCase, isfile(tn13), 'missing fixture: ' + string(tn13));
fol = tempname; mkdir(fol);
cl = onCleanup(@() rmIfFolder(fol)); %#ok<NASGU>
copyfile(src, fullfile(fol, 'ITSG-Grace2018_n60_2008-04.gfc'));
copyfile(src, fullfile(fol, 'ITSG-Grace2018_n60_2008-05.gfc'));
[ts, rep] = shLowLevel.standardChain(fol, TN14File = tn14, ...
    Degree1File = tn13, Filter = "gauss300", Quiet = true);
verifyEqual(testCase, ts.nEpochs, 2);
verifyEqual(testCase, numel(rep.steps), 4);      % read/TN14/deg1/filter
verifyTrue(testCase, contains(rep.version, "shAnalysis"));
verifyTrue(testCase, any(contains(string(ts.history), "gaussian", ...
    'IgnoreCase', true)) || numel(ts.history) >= 3);
% GIA: with t0 = first epoch, only the second epoch shifts
gT = shCoefficients(zeros(61), zeros(61));
gT = gT.setCoefficient(2, 0, 1e-9, NaN);
[t2, ~] = shLowLevel.standardChain(fol, TN14 = false, Degree1 = "none", ...
    GIA = gT, GIAEpoch = NaN, Filter = "none", Quiet = true);
[t0, ~] = shLowLevel.standardChain(fol, TN14 = false, Degree1 = "none", ...
    Filter = "none", Quiet = true);
dt = t2.epochs(2) - t2.epochs(1);
g2 = t2.at(2); g0 = t0.at(2);
verifyEqual(testCase, g0.C(3, 1) - g2.C(3, 1), (dt / 2) * 1e-9, ...
    'RelTol', 1e-10);                            % t0 = mean epoch
% custom W path (designFilter output drops in)
g1 = t0.at(1);
sc = 3e-11 * ones(61);
W = shLowLevel.designFilter(sc, sc, Kaula = 1e-6);
[tw, repW] = shLowLevel.standardChain(fol, TN14 = false, ...
    Degree1 = "none", Filter = W, Quiet = true);
verifyEqual(testCase, tw.nEpochs, 2);
verifyTrue(testCase, any(contains(repW.steps, "custom W")));
verifyError(testCase, @() shLowLevel.standardChain(fol, ...
    Filter = "bogus", TN14 = false, Degree1 = "none", Quiet = true), ...
    'shLowLevel:standardChain:badFilter');
verifyError(testCase, @() shLowLevel.standardChain(fol, ...
    TN14File = "no/such/file.txt", Degree1 = "none", Quiet = true), ...
    'shLowLevel:standardChain:noTN14');
end

function testVersionMetadata(testCase)
% shLowLevel.version reports toolbox metadata parsed from Contents.m - the
% single source of truth also honoured by MATLAB's ver().
v = shLowLevel.version();
verifyEqual(testCase, v.Name, "shAnalysis");
verifyTrue(testCase, ~isempty(regexp(v.Version, '^\d+\.\d+', 'once')));
verifyTrue(testCase, isfolder(v.Root));
txt = fileread(fullfile(v.Root, 'Contents.m'));
tok = regexp(txt, '^%\s*Version\s+(\S+).*?(\d{2}-\w{3}-\d{4})\s*$', ...
    'tokens', 'once', 'lineanchors');
verifyFalse(testCase, isempty(tok));
verifyEqual(testCase, v.Version, string(tok{1}));
verifyEqual(testCase, v.Date, string(tok{2}));
verifyTrue(testCase, contains(v.Provenance, "Claude"));
end

function testFetchITSGInputIDs(testCase)
% input validation only - no network in unit tests
verifyError(testCase, @() shLowLevel.fetchITSG("2010-4"), 'shLowLevel:fetchITSG:badMonth');
verifyError(testCase, @() shLowLevel.fetchITSG("April 2010"), 'shLowLevel:fetchITSG:badMonth');
verifyError(testCase, @() shLowLevel.fetchITSG(1990), 'shLowLevel:fetchITSG:badMonth');
% daily product (v2.4.1): n40 only; monthly rejects n40
verifyError(testCase, ...
    @() shLowLevel.fetchITSG("2010-04", Product = "daily", Nmax = 96), ...
    'shLowLevel:fetchITSG:badNmax');
verifyError(testCase, @() shLowLevel.fetchITSG("2010-04", Nmax = 40), ...
    'shLowLevel:fetchITSG:badNmax');
end

function testSpectrumPlotOptionsAndAxes(testCase)
rng(102);
fig = figure('Visible', 'off');
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
L = 15; n1 = L + 1;
C = tril(randn(n1)) * 1e-9; S = tril(randn(n1), -1) * 1e-9; S(:, 1) = 0;
sd = shLowLevel.shDegreeRMS(C, S, 'sigmaC', abs(C)*0.1, 'sigmaS', abs(S)*0.1);
so = shLowLevel.shOrderRMS(C, S);
ax = axes('Parent', fig);
% linear x-axis contract (v2.4.1)
shLowLevel.plotSHSpectrum(sd, 'ax', ax);
verifyEqual(testCase, ax.XScale, 'linear');
verifyEqual(testCase, ax.YScale, 'log');
% every quantity renders, degree and order domain
for q = ["amplitude", "rms", "variance", "cumamplitude", "cumrms", "cumvariance"]
    cla(ax); shLowLevel.plotSHSpectrum(sd, 'ax', ax, 'Quantity', char(q));
    cla(ax); shLowLevel.plotSHSpectrum(so, 'ax', ax, 'Quantity', char(q));
end
% contracts
verifyError(testCase, ...
    @() shLowLevel.plotSHSpectrum({sd, so}, 'ax', axes('Parent', fig)), ...
    'shLowLevel:plotSHSpectrum:mixedDomains');
verifyError(testCase, ...
    @() shLowLevel.plotSHSpectrum(sd, 'ax', axes('Parent', fig), ...
    'Quantity', 'variance', 'Kaula', 1e-5), ...
    'shLowLevel:plotSHSpectrum:kaulaDomain');
verifyError(testCase, ...
    @() shLowLevel.plotSHSpectrum(sd, 'ax', axes('Parent', fig), 'Quantity', 'xxx'), ...
    'shLowLevel:plotSHSpectrum:badQuantity');
% triangle: no center line anymore (v2.4.1)
ax2 = axes('Parent', fig);
shLowLevel.plotSHCoeffTriangle(C, S, 'ax', ax2);
verifyEmpty(testCase, findobj(ax2, 'Type', 'constantline'));
end

function testDataFolderAndDDKNames(testCase)
% preserve any user preference across the test
had = ispref('shAnalysis', 'dataFolder');
if had
    old = getpref('shAnalysis', 'dataFolder');
    restore = onCleanup(@() setpref('shAnalysis', 'dataFolder', old)); %#ok<NASGU>
else
    restore = onCleanup(@() shLowLevel.dataFolder("reset")); %#ok<NASGU>
end
tmp = fullfile(tempdir, sprintf('shx_data_%d', randi(1e9)));
f = shLowLevel.dataFolder(tmp);
verifyEqual(testCase, char(f), char(string(tmp)));
verifyTrue(testCase, isfolder(f));
verifyEqual(testCase, char(shLowLevel.dataFolder()), char(string(tmp)));
f2 = shLowLevel.dataFolder("reset");
verifyTrue(testCase, endsWith(f2, "data"));
% DDK mapping: 8 unique released names, DDK3 matches the shipped file
nm = shLowLevel.ddkNames();
verifyEqual(testCase, numel(nm), 8);
verifyEqual(testCase, numel(unique(nm)), 8);
verifyEqual(testCase, nm(3), "Wbd_2-120.a_1d12p_4");
% name resolution: DDK3 loads from shipped test data even with a fresh
% (empty) data folder; an unfetched filter errors with the fetch hint
shLowLevel.dataFolder(tmp);
W = shLowLevel.readDDK("DDK3", Nmax = 30);
verifyEqual(testCase, W.nmax, 30);
verifyError(testCase, @() shLowLevel.readDDK("DDK7"), 'shLowLevel:readDDK:notFetched');
verifyError(testCase, @() shLowLevel.fetchDDK(0), 'MATLAB:validators:mustBeInRange');
end


% ------------------------------------------------- tvANS option forwarding
function testFilterForwardsTvANSOptions(testCase)
%TESTFILTERFORWARDSTVANSOPTIONS ts.filter is the full single point of access.
%   The workflow guide has advertised
%   ts.filter("tvANS", Blocks="auto", VCEBands=[...]) since Edition 2 while
%   the method rejected the option (MATLAB:TooManyInputs). Pin the three
%   forwarded tuning options AND the "[] does not override" contract, so
%   the defaults keep exactly one home (shLowLevel.tvANSFilter).
rng(11);
L = 6; T = 30; n1 = L + 1;
tY = 2019 + (0:T-1)'/12;
mL = tril(true(n1)); mL1 = mL; mL1(:, 1) = false;
Cs = zeros(n1, n1, T); Ss = Cs;
for t = 1:T
    Cs(:,:,t) = (mL .* randn(n1)) * 1e-9 + (mL * 2e-9) * cos(2*pi*tY(t));
    Ss(:,:,t) = (mL1 .* randn(n1)) * 1e-9;
end
ts = shSeries(Cs, Ss = Ss, Epochs = tY);

% 1. the three options are accepted (this is what used to throw)
[tsB, ~, infoB] = ts.filter("tvANS", Blocks = "on", VCEBands = [0 3 7]);
verifyClass(testCase, tsB, 'shSeries');
verifySize(testCase, tsB.Cs, [n1 n1 T]);

% 2. banding actually reaches tvANSFilter: the banded run must differ from
%    the unbanded one (a silently dropped option would give equality)
[tsU, ~, infoU] = ts.filter("tvANS", Blocks = "on");
verifyTrue(testCase, max(abs(tsB.Cs(:) - tsU.Cs(:))) > 0);
verifyTrue(testCase, isfield(infoB, 'sigmaXfres') && ...
    isfield(infoU, 'sigmaXfres'));

% 3. [] means "do not override": explicit empties reproduce the defaults
%    bit for bit, so no default is duplicated in shSeries
tsE = ts.filter("tvANS", Blocks = "on", Shrinkage = [], ...
    VCEMinDegree = [], VCEBands = []);
verifyEqual(testCase, tsE.Cs, tsU.Cs);
verifyEqual(testCase, tsE.Ss, tsU.Ss);

% 4. Shrinkage and VCEMinDegree are forwarded too
tsS = ts.filter("tvANS", Blocks = "on", Shrinkage = 0.5);
verifyTrue(testCase, max(abs(tsS.Cs(:) - tsU.Cs(:))) > 0);
tsV = ts.filter("tvANS", Blocks = "on", VCEMinDegree = 2);
verifyTrue(testCase, max(abs(tsV.Cs(:) - tsU.Cs(:))) > 0);

% 5. scalar validators still bite
verifyError(testCase, @() ts.filter("tvANS", Shrinkage = [0.1 0.2]), ...
    'MATLAB:validators:mustBeScalarOrEmpty');

% 6. the low-level entry point rejects banding off the block path -
%    the forwarded option must not swallow that guard
verifyError(testCase, ...
    @() ts.filter("tvANS", Blocks = "off", VCEBands = [0 3 7]), ...
    'shLowLevel:tvANSFilter:bandsNeedBlocks');
end

% --------------------------------------------- documentation-metadata sync
function testVersionMetadataIsConsistent(testCase)
%TESTVERSIONMETADATAISCONSISTENT One version, parsed from Contents.m.
%   Contents.m is the single source of truth (five bugs died from
%   hardcoded version strings). CITATION.cff, the CHANGELOG top section
%   and the generated API reference must agree with it, and ver() must
%   report a short product name rather than a truncated sentence.
root = fileparts(fileparts(mfilename('fullpath')));
v = shLowLevel.version();
verifyMatches(testCase, v.Version, '^\d+\.\d+\.\d+$');

cff = string(fileread(fullfile(root, 'CITATION.cff')));
tok = regexp(cff, 'version:\s*"([^"]+)"', 'tokens', 'once');
verifyEqual(testCase, string(tok{1}), v.Version, ...
    'CITATION.cff version must match Contents.m');

chg = string(fileread(fullfile(root, 'CHANGELOG.md')));
tok = regexp(chg, '##\s*\[([^\]]+)\]', 'tokens', 'once');
verifyEqual(testCase, string(tok{1}), v.Version, ...
    'the top CHANGELOG section must be the current version');
openSec = regexp(chg, '^##\s*\[[^\]]+\]\s*-\s*Unreleased', ...
    'match', 'lineanchors');
verifyEmpty(testCase, openSec, ...
    'a tagged release must not leave "Unreleased" sections behind');

api = string(fileread(fullfile(root, 'html', 'apiReference.html')));
verifyTrue(testCase, contains(api, "API reference (v" + v.Version + ")"), ...
    'html/apiReference.html is stale - regenerate it');

% ver() name: short, no version number, no dangling sentence
mv = ver('shAnalysis');
assumeNotEmpty(testCase, mv, ...
    'ver() needs the toolbox folder to be named shAnalysis');
verifyEqual(testCase, numel(mv), 1);
verifyEqual(testCase, string(mv.Version), v.Version);
verifyLessThan(testCase, strlength(string(mv.Name)), 60);
verifyFalse(testCase, contains(string(mv.Name), digitsPattern));
end

function testHelpBrowserTocIsComplete(testCase)
%TESTHELPBROWSERTOCISCOMPLETE Every help page is reachable from helptoc.xml.
%   Eleven pages were orphaned (invisible in the Help browser) while the
%   files sat in html/. Both directions are pinned: no orphan pages, no
%   dead targets.
root = fileparts(fileparts(mfilename('fullpath')));
hdir = fullfile(root, 'html');
toc = string(fileread(fullfile(hdir, 'helptoc.xml')));
tk = regexp(toc, 'target="([^"]+)"', 'tokens');
targets = unique(string(cellfun(@(c) c{1}, tk, 'UniformOutput', false)));
d = dir(fullfile(hdir, '*.html'));
files = string({d.name}');
verifyEmpty(testCase, setdiff(files, targets), ...
    'help pages missing from helptoc.xml');
verifyEmpty(testCase, setdiff(targets, files), ...
    'helptoc.xml points at files that do not exist');
% the removed compat/ folder must not be advertised anywhere as current
lp = string(fileread(fullfile(hdir, 'shAnalysis.html')));
verifyFalse(testCase, contains(lp, "remains available"), ...
    'the landing page still advertises the removed compat/ folder');
cm = string(fileread(fullfile(root, 'Contents.m')));
verifyFalse(testCase, contains(cm, "folder compat/, add to path"), ...
    'Contents.m still advertises the removed compat/ folder');
end

% ------------------------------------------------ in-file help completeness
function testHelpHasNoPlaceholders(testCase)
%TESTHELPHASNOPLACEHOLDERS No option is documented as "see arguments block".
%   58 name-value options across 23 entities used to be "documented" with
%   a pointer to the code, so `help` said nothing while the generated
%   apiReference.html showed size, type and default. Pin that to zero:
%   the in-file help is the primary source, not a redirect.
root = fileparts(fileparts(mfilename('fullpath')));
files = [ ...
    string(fullfile(root, {'shCoefficients.m', 'shSeries.m', ...
                           'shClimatology.m', 'setup_shAnalysis.m'})), ...
    string(fullfile(root, '+shLowLevel', {dir(fullfile(root, ...
        '+shLowLevel', '*.m')).name}))];
bad = strings(0, 1);
for f = files
    if ~isfile(f), continue; end
    txt = string(fileread(f));
    if contains(txt, "see arguments block")
        bad(end+1) = f; %#ok<AGROW>
    end
end
verifyEmpty(testCase, bad, sprintf( ...
    'placeholder help ("see arguments block") in: %s', ...
    strjoin(bad, ', ')));
end

% -------------------------------------------------- documentation sync gate
function testDocSyncAuditIsWired(testCase)
%TESTDOCSYNCAUDITISWIRED The sixth gate exists and runs in CI.
%   The five older gates were all green while the API reference was
%   stale, eleven help pages were unreachable and the guide advertised a
%   call that threw. tools/doc_sync_audit.py closes that hole, so pin
%   that it is present and actually wired into the required CI job -
%   a gate nobody runs is not a gate.
root = fileparts(fileparts(mfilename('fullpath')));
gate = fullfile(root, 'tools', 'doc_sync_audit.py');
verifyTrue(testCase, isfile(gate), 'tools/doc_sync_audit.py is missing');

wf = fullfile(root, '.github', 'workflows', 'ci.yml');
assumeTrue(testCase, isfile(wf), 'no workflow file in this checkout');
txt = string(fileread(wf));
verifyTrue(testCase, contains(txt, "tools/doc_sync_audit.py"), ...
    'doc_sync_audit.py is not run by the CI workflow');
% it must run BEFORE MATLAB: a python gate that only fires after a
% 2-minute toolbox install wastes the fast feedback it exists for
iGate = strfind(txt, "tools/doc_sync_audit.py");
iML = strfind(txt, "matlab-actions/setup-matlab");
verifyTrue(testCase, ~isempty(iML) && iGate(1) < iML(1), ...
    'doc_sync_audit.py must run before the MATLAB setup step');
end

% ------------------------------------------------------------- safeMove
function testSafeMoveContract(testCase)
%TESTSAFEMOVECONTRACT The fetchers' final swap: verify, retry, report.
%   safeMove replaces a bare movefile in all six fetchers. It must move
%   the file, VERIFY the outcome rather than trust the status flag,
%   retry a configurable number of times, and fail with an actionable
%   identifier rather than a bare OS message.
d = fullfile(tempdir, sprintf('shx_move_%d', randi(1e9)));
mkdir(d);
cl = onCleanup(@() rmIfFolder(d)); %#ok<NASGU>
src = fullfile(d, 'a.bin');
dst = fullfile(d, 'b.bin');

% happy path: one attempt, file actually moved
writematrix(magic(4), src, 'FileType', 'text');
n = shLowLevel.safeMove(src, dst);
verifyEqual(testCase, n, 1);
verifyTrue(testCase, isfile(dst));
verifyFalse(testCase, isfile(src));

% overwriting an existing destination is the normal fetcher case
writematrix(magic(3), src, 'FileType', 'text');
verifyEqual(testCase, shLowLevel.safeMove(src, dst), 1);
verifyEqual(testCase, readmatrix(dst, 'FileType', 'text'), magic(3));

% a missing source is a caller bug, not a lock: say so immediately
verifyError(testCase, @() shLowLevel.safeMove(fullfile(d, 'nope'), dst), ...
    'shLowLevel:safeMove:noSource');

% a missing destination folder cannot be cured by waiting: report it
% at once, and do not blame a scanner for a caller bug
writematrix(1, src, 'FileType', 'text');
bad = fullfile(d, 'no_such_dir', 'c.bin');
t0 = tic;
verifyError(testCase, @() shLowLevel.safeMove(src, bad), ...
    'shLowLevel:safeMove:noDestFolder');
verifyLessThan(testCase, toc(t0), 1, ...
    'an impossible move must not spend the retry budget');
verifyTrue(testCase, isfile(src), ...
    'a failed move must leave the source in place');
end


function testSafeMoveSurvivesALockedDestination(testCase)
%TESTSAFEMOVESURVIVESALOCKEDDESTINATION The case this function exists for.
%   On Windows an antivirus or sync client holds a freshly written file
%   for a moment and movefile fails with a sharing violation. Reproduced
%   here by keeping the destination open and releasing it from a timer
%   while safeMove retries.
%
%   This is the ONLY test of the retry loop, deliberately. POSIX renames
%   over open files, and permission tricks do not substitute: the CI
%   runner is root, so a read-only destination folder does not block it
%   either (tried; the move succeeded). Rather than fake a failure on
%   Linux, the loop is tested on the platform whose behaviour it exists
%   for - which is also the acceptance machine.
assumeTrue(testCase, ispc, ...
    'file locking during a move is Windows behaviour');
d = fullfile(tempdir, sprintf('shx_lock_%d', randi(1e9)));
mkdir(d);
cl = onCleanup(@() rmIfFolder(d)); %#ok<NASGU>
src = fullfile(d, 'a.bin');
dst = fullfile(d, 'b.bin');
writematrix(magic(4), src, 'FileType', 'text');
fid = fopen(dst, 'w');
fwrite(fid, 'held by a scanner');
verifyFalse(testCase, movefile(src, dst, 'f'), ...
    'the destination is not actually locked - test cannot conclude');

% Retries = 0 is a plain movefile: it must fail against the held file,
% spend no time waiting, and leave the source in place
t0 = tic;
verifyError(testCase, @() shLowLevel.safeMove(src, dst, Retries = 0), ...
    'shLowLevel:safeMove:locked');
verifyLessThan(testCase, toc(t0), 1);
verifyTrue(testCase, isfile(src), ...
    'a failed move must leave the source in place');

% with retries, the lock is outlived: released at 0.4 s from a timer
rel = timer('StartDelay', 0.4, 'TimerFcn', @(~, ~) fclose(fid));
clT = onCleanup(@() delete(rel)); %#ok<NASGU>
start(rel);
t0 = tic;
n = shLowLevel.safeMove(src, dst, Pause = 0.15, Retries = 6);
verifyGreaterThan(testCase, n, 1, ...
    'the first attempt must have failed against the held file');
verifyGreaterThanOrEqual(testCase, toc(t0), 0.15);   % backoff was spent
verifyTrue(testCase, isfile(dst));
verifyFalse(testCase, isfile(src));
end

% ------------------------------------- correction tables trail the data
function testStandardChainHandlesTrailingCorrectionTables(testCase)
%TESTSTANDARDCHAINHANDLESTRAILINGCORRECTIONTABLES The newest months.
%   TN-13 and TN-14 are published weeks after the monthly solutions, so
%   the newest one or two epochs of a fresh series are ALWAYS uncovered.
%   The chain used to stop with "No TN-13 entry within 0.050 yr", which
%   made routine use of an up-to-date series impossible. Uncovered
%   epochs are now dropped and the fact recorded in the report.
dd = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'tests', 'test_data');
fol = fullfile(tempdir, sprintf('shx_chain_%d', randi(1e9)));
mkdir(fol);
cl = onCleanup(@() rmIfFolder(fol)); %#ok<NASGU>

% one solution the tables cover, and one far in the future that no
% table can reach - the situation a fresh download produces
% epoch is read-only on the value class; a series takes its epochs from
% the ITSG-style filename, which is also how a real folder works
g = shCoefficients.read(fullfile(dd, 'ITSG-Grace2018_n60_2008-04.gfc'));
g.write(fullfile(fol, 'ITSG-Grace2018_n60_2008-04.gfc'), Sidecar = false);
g.write(fullfile(fol, 'ITSG-Grace2018_n60_2099-07.gfc'), Sidecar = false);

tn14 = fullfile(dd, 'TN-14_C30_C20_SLR_GSFC.txt');
tn13 = fullfile(dd, 'TN-13_GEOC_CSR_RL06.3.txt');

[ts, rep] = shLowLevel.standardChain(fol, TN14File = tn14, ...
    Degree1File = tn13, Filter = "none", Quiet = true);
verifyEqual(testCase, ts.nEpochs, 1, ...
    'the uncovered epoch must be dropped, not carried uncorrected');
verifyLessThan(testCase, ts.epochs, 2099);
verifyTrue(testCase, any(contains(rep.steps, "dropped")), ...
    'a dropped epoch must be visible in the provenance report');

% the old behaviour is still available for anyone who wants it
verifyError(testCase, @() shLowLevel.standardChain(fol, ...
    TN14File = tn14, Degree1File = tn13, Filter = "none", ...
    OnMissing = "error", Quiet = true), ...
    'shLowLevel:standardChain:uncoveredEpochs');

% a series with no overlap at all is a different mistake and says so
folF = fullfile(tempdir, sprintf('shx_chain_f_%d', randi(1e9)));
mkdir(folF);
clF = onCleanup(@() rmIfFolder(folF)); %#ok<NASGU>
g.write(fullfile(folF, 'ITSG-Grace2018_n60_2099-07.gfc'), Sidecar = false);
verifyError(testCase, @() shLowLevel.standardChain(folF, ...
    TN14File = tn14, Degree1File = tn13, Filter = "none", ...
    Quiet = true), 'shLowLevel:standardChain:noCoveredEpochs');
end

function testProvenanceSidecarsDoNotBreakTheRoundTrip(testCase)
%TESTPROVENANCESIDECARSDONOTBREAKTHEROUNDTRIP Write a folder, read it back.
%   writeGFC drops "<file>.gfc.provenance.json" beside every export, and
%   the read pattern has to be loose enough for ".gfc.gz" - so "*.gfc*"
%   matched the sidecars and a folder written BY the toolbox could not be
%   read back BY the toolbox. Also pins that the class API can switch the
%   sidecar off, which it previously could not.
dd = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'tests', 'test_data');
fol = fullfile(tempdir, sprintf('shx_sidecar_%d', randi(1e9)));
mkdir(fol);
cl = onCleanup(@() rmIfFolder(fol)); %#ok<NASGU>
g = shCoefficients.read(fullfile(dd, 'ITSG-Grace2018_n60_2008-04.gfc'));

g.write(fullfile(fol, 'ITSG-Grace2018_n60_2008-04.gfc'));   % sidecar on
g.write(fullfile(fol, 'ITSG-Grace2018_n60_2009-04.gfc'));
verifyTrue(testCase, isfile(fullfile(fol, ...
    'ITSG-Grace2018_n60_2008-04.gfc.provenance.json')), ...
    'the sidecar is expected - this test is about reading past it');

ts = shSeries.fromFolder(fol);
verifyEqual(testCase, ts.nEpochs, 2, ...
    'sidecars must not be read as solutions');
verifyEqual(testCase, sort(ts.epochs(:)), [2008.292; 2009.292], ...
    'AbsTol', 0.01);

% and the class API can suppress it
g.write(fullfile(fol, 'plain.gfc'), Sidecar = false);
verifyFalse(testCase, isfile(fullfile(fol, 'plain.gfc.provenance.json')));
end

% ------------------------------------------------- generated tutorials
function testMakeTutorialsFollowsTheRegistry(testCase)
%TESTMAKETUTORIALSFOLLOWSTHEREGISTRY Tutorials are generated, not curated.
%   demo_shAnalysis is the single source of truth for what the toolbox
%   demonstrates. A parallel set of hand-written tutorials would drift
%   from it within a release - and .mlx is a binary zip, so the drift
%   would not even show in a diff. Pin that every case gets a tutorial
%   and that the content comes from the registry.
d = fullfile(tempdir, sprintf('shx_tut_%d', randi(1e9)));
cl = onCleanup(@() rmIfFolder(d)); %#ok<NASGU>

out = shLowLevel.makeTutorials(Dest = d, Cases = ["D01" "D05"], ...
    Convert = false, Quiet = true);
verifyEqual(testCase, numel(out.files), 2);
verifyEmpty(testCase, out.mlx);
verifyTrue(testCase, all(isfile(out.files)));

reg = demo_shAnalysis("list");
r1 = reg(string({reg.id}) == "D01");
txt = string(fileread(out.files(1)));
verifyTrue(testCase, contains(txt, r1.title), ...
    'the title must come from the registry');
verifyTrue(testCase, contains(txt, r1.fns), ...
    'the function list must come from the registry');
verifyTrue(testCase, contains(txt, "demo_shAnalysis(""D01"""), ...
    'the tutorial must actually run its case');
verifyTrue(testCase, startsWith(txt, "%%"), ...
    'Live Editor section markers are what make this a live script');
verifyTrue(testCase, contains(txt, "overwritten on the next run"), ...
    'a generated file must say it is generated');

% "all" covers every registered case, one file each
outAll = shLowLevel.makeTutorials(Dest = d, Convert = false, Quiet = true);
verifyEqual(testCase, numel(outAll.files), numel(reg));

verifyError(testCase, @() shLowLevel.makeTutorials(Dest = d, ...
    Cases = "D99", Convert = false, Quiet = true), ...
    'shLowLevel:makeTutorials:unknownCase');
end
