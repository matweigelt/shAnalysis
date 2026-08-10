function tests = testContract
%TESTCONTRACT API contract tests for the shAnalysis v2 toolbox.
%
%   Verifies the *interface*, independent of numerics: every documented
%   error path throws its documented identifier, value objects are
%   immutable (methods return new objects, originals untouched), output
%   dimensions match the documented contracts, and the compat/ wrappers
%   remain signature-identical to shx internals.
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
cdir = fullfile(root, 'compat');            % private, optional
if isfolder(cdir), addpath(cdir); end
tc.TestData.root = root;
tc.TestData.dataDir = fullfile(here, 'test_data');
shx.legendreCached('clear');
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
verifyError(tc, @() shx.readTN14('no_such_file.txt'), ...
    'shx:readTN14:fileNotFound');
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
    @() shx.shSynthesis(g.C, g.S, g.GM, g.R, 0, 0, 'quantity', 'nonsense'), ...
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

idx = shx.shIndex(g.nmax, MinDegree = 2);
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
function testCompatWrappersDelegate(tc)
% cross-validation against the PRIVATE v1 reference implementations
% (compat/ is not published); runs locally, filtered on CI
assumeTrue(tc, isfolder(fullfile(fileparts( ...
    fileparts(mfilename('fullpath'))), 'compat')), ...
    'compat reference implementations not present (unpublished)');
g = randomField(12, 2020);
[c1, s1] = shDestripe(g.C, g.S, 'minOrder', 4);        % compat
[c2, s2] = shx.shDestripe(g.C, g.S, 'minOrder', 4);    % internal
verifyEqual(tc, c1, c2, AbsTol = 0); verifyEqual(tc, s1, s2, AbsTol = 0);
Wn1 = shGaussianWeights(20, 300); Wn2 = shx.shGaussianWeights(20, 300);
verifyEqual(tc, Wn1, Wn2, AbsTol = 0);
P1 = legendreALF(10, [-0.5, 0, 0.5]); P2 = shx.legendreALF(10, [-0.5, 0, 0.5]);
verifyEqual(tc, P1, P2, AbsTol = 0);
end

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
verifyError(testCase, @() shx.shAnalysisGrid(ones(3,4), 1:5, 1:4, 2), ...
    'shx:shAnalysisGrid:badInput');
end

function testAnalysisBadGridID(testCase)
% non-uniform longitude must be rejected by the rings method
g = ones(8, 10);
lat = linspace(-60, 60, 8);
lon = [0 30 70 100 140 180 220 260 300 340];
verifyError(testCase, ...
    @() shx.shAnalysisGrid(g, lat, lon, 2, Method = "rings"), ...
    'shx:shAnalysisGrid:badGrid');
end

function testAnalysisAliasedGridID(testCase)
% nlon <= 2*nmax aliases orders
g = ones(20, 10); lat = linspace(-80, 80, 20); lon = (0:9)*36;
verifyError(testCase, @() shx.shAnalysisGrid(g, lat, lon, 6, Method = "rings"), ...
    'shx:shAnalysisGrid:badGrid');
end

function testSynthesisBadMethodID(testCase)
C = zeros(3); S = zeros(3);
verifyError(testCase, @() shx.shSynthesis(C, S, 1, 1, 0, [0 10 50], ...
    'method', 'fft'), 'shSynthesis:badMethod');
end

function testReadTN13ErrorIDs(testCase)
verifyError(testCase, @() shx.readTN13('/nonexistent/tn13.txt'), ...
    'shx:readTN13:fileNotFound');
f = [tempname '.txt'];
fid = fopen(f, 'w'); fprintf(fid, 'header only\nno data here\n'); fclose(fid);
cleanup = onCleanup(@() delete(f)); %#ok<NASGU>
verifyError(testCase, @() shx.readTN13(f), 'shx:readTN13:noData');
end

function testReadSINEXErrorIDs(testCase)
verifyError(testCase, @() shx.readSINEX('/nonexistent/x.snx'), ...
    'shx:readSINEX:fileNotFound');
f = [tempname '.snx'];
fid = fopen(f, 'w'); fprintf(fid, '%%=SNX 2.02\n+SOME/BLOCK\n-SOME/BLOCK\n%%ENDSNX\n'); fclose(fid);
cleanup = onCleanup(@() delete(f)); %#ok<NASGU>
verifyError(testCase, @() shx.readSINEX(f), 'shx:readSINEX:noData');
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
model = shx.shReadGFC(f);
verifyError(testCase, @() shx.shEvalGFCT(model, 2020.0), ...
    'shEvalGFCT:epochOutside');
end

function testTvANSBlocksUnavailableID(testCase)
idx = shx.shIndex(4);
X = randn(idx.P, 12); t = 2002 + (0:11)'/12;
verifyError(testCase, @() shx.tvANSFilter(X, t, idx, Blocks = 'on', ...
    Constraints = ones(idx.P, 1)), 'shx:tvANSFilter:blocksUnavailable');
end

function testNoiseCovNotBlockDiagonalID(testCase)
% external NoiseCov on Blocks='on' must be block-diagonal in the
% (order, C/S, parity) partition; Blocks='auto' falls back quietly
rng(30);
idx = shx.shIndex(4);
X = randn(idx.P, 12); t = 2002 + (0:11)'/12;
Nfull = eye(idx.P) + 0.3 * ones(idx.P);          % dense: NOT block-diag
verifyError(testCase, @() shx.tvANSFilter(X, t, idx, Blocks = 'on', ...
    NoiseCov = Nfull), 'shx:buildNoiseCov:notBlockDiagonal');
[Xf, op] = shx.tvANSFilter(X, t, idx, Blocks = 'auto', NoiseCov = Nfull);
verifyTrue(testCase, all(isfinite(Xf(:))));
verifyEqual(testCase, op.layout, 'full');        % fell back to full path
% a DIAGONAL external N is block-diagonal under any partition: block
% path engages
[~, op2] = shx.tvANSFilter(X, t, idx, Blocks = 'on', NoiseCov = eye(idx.P));
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
plotSHCoeffTriangle(C, S, 'RefC', 0*C, 'RefS', 0*S, 'ax', ax);
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
grid = shx.shSynthesis(C, S, 3.986004415e14, 6378136.3, ...
    -85:5:85, 0:5:355, 'quantity', 'geoid');
ax1 = axes('Parent', fig);
h = shx.plotSHMap(grid, -85:5:85, 0:5:355, ax = ax1, Units = "m");
verifyClass(testCase, h, 'matlab.graphics.axis.Axes');
cla(ax1);
h2 = shx.plotSHMap(grid, -85:5:85, 0:5:355, ax = ax1, ...
    Projection = "hammer");
verifyClass(testCase, h2, 'matlab.graphics.axis.Axes');
% basin series with gap + band
t = [2010 + (0:40)/12, 2015 + (0:30)/12];
c = sin(2*pi*t(:)) + 0.01*(t(:) - 2012);
cla(ax1);
h3 = shx.plotBasinSeries(t(:), c, 0.2 + 0*c, ax = ax1, Units = "cm");
verifyClass(testCase, h3, 'matlab.graphics.axis.Axes');
% covariance plot
idx = shx.shIndex(6, MinDegree = 2);
A = randn(idx.P); M = A*A';
cla(ax1);
h4 = shx.plotCovariance(M, idx, ax = ax1);
verifyClass(testCase, h4, 'matlab.graphics.axis.Axes');
% triangle diff mode
cla(ax1);
h5 = plotSHCoeffTriangle(C, S, 'RefC', 0.9*C, 'RefS', 0.9*S, 'ax', ax1);
verifyClass(testCase, h5, 'matlab.graphics.axis.Axes');
verifyError(testCase, ...
    @() plotSHCoeffTriangle(C, S, 'RefC', C, 'ax', axes('Parent', fig)), ...
    'shx:plotSHCoeffTriangle:needRefS');
end

function testSpectrumOverlays(testCase)
rng(82);
fig = figure('Visible', 'off');
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
L = 30; n1 = L + 1;
C = tril(randn(n1)) * 1e-9 ./ max(1, (0:L)').^2;
S = tril(randn(n1), -1) * 1e-9 ./ max(1, (0:L)').^2; S(:, 1) = 0;
sC = 1e-11 * max(1, (0:L)') / 10 .* ones(n1) .* tril(ones(n1));
spec = shx.shDegreeRMS(C, S, 'R', 6378136.3, ...
    'sigmaC', sC, 'sigmaS', sC);
ax = axes('Parent', fig);
h = plotSHSpectrum(spec, 'ax', ax, 'Kaula', 1e-5, 'MarkCrossover', true);
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
shx.writeAnimation(ts, tmp, quantity = "geoid", ...
    lat = -80:10:80, lon = 0:10:350, FrameRate = 2);
verifyTrue(testCase, isfile(tmp));
d = dir(tmp);
verifyGreaterThan(testCase, d.bytes, 1000);
end

function deleteIfThere(f)
if isfile(f), delete(f); end
end

function testBasinScalingRecoversFactor(testCase)
rng(84);
% operator: pure Gaussian attenuation as a tvANS-like op via matFilter
L = 20; idx = shx.shIndex(L, MinDegree = 0);
Wn = shx.shGaussianWeights(L, 500);
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
b = shx.basinKernel(idx, [ -10 100; -10 140; 25 140; 25 100 ]);
[k, info] = shx.basinScaling(W, b, tsM, idx = idx, tYears = ep);
% k must invert the operator's basin-average attenuation: with a FIXED
% pattern, aTrue = k * aFilt exactly, so sigmaK ~ 0
x1 = shx.vecFromCS(pat, patS, idx);
kExact = (b' * x1) / (b' * (W * x1));
verifyEqual(testCase, k, kExact, 'RelTol', 1e-10);
verifyLessThan(testCase, info.sigmaK, 1e-8 * abs(k));
verifyEqual(testCase, info.nMatched, T);
verifyError(testCase, ...
    @() shx.basinScaling(W, b(1:3), tsM, idx = idx, tYears = ep), ...
    'shx:basinScaling:badKernel');
end

function testDemoRegistryAndSmoke(testCase)
reg = demo_shAnalysis("list");
verifyEqual(testCase, numel(reg), 16);
verifyEqual(testCase, numel(unique([reg.id])), 16);
verifyTrue(testCase, all(arrayfun(@(r) isa(r.run, 'function_handle'), reg)));
% two cheap cases run headless without touching the screen
demo_shAnalysis(["D01", "D13"], Visible = false);
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
cleanup2 = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>
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
ws.idx = shx.shIndex(10);
ws.x = shx.vecFromCS(ws.g10.C, ws.g10.S, ws.idx);
ws.GM = ws.g.GM; ws.R = ws.g.R;
ws.GMref = 3.986004415e14; ws.Rref = 6.3781363e6;
nn = (0:120)';
ws.kn = -0.3 * nn ./ (nn + 6);
% Love-number table via the reader (real-life G7 path)
lnRows = compose("%d %.6f %.6f %.6f", nn, -0.9 * nn ./ (nn + 4), ...
    0.1 * nn ./ (nn + 8), -0.3 * nn ./ (nn + 6));
writelines(["# n h l k"; lnRows], "ln_table.txt");
ws.LN = shx.readLoveNumbers("ln_table.txt", MaxDegree = 120);
copyfile("ln_table.txt", "prem_load.txt");   % readLoveNumbers help example
ws.lat = -88:8:88; ws.lon = 0:12:348;    % = the grid axes below
ws.latPts = [10 20 30]; ws.lonPts = [40 50 60];
ws.latGc = [10 20 30]; ws.latGd = [10.06 20.10 30.12];
ws.theta = deg2rad(30:15:60);
ws.grid = ws.g10.synthesis(-88:8:88, 0:12:348);
ws.tn13 = shx.readTN13("TN-13_GEOC_CSR_RL06.3.txt");
ws.tn14 = shx.readTN14("TN-14_C30_C20_SLR_GSFC.txt");
ws.model = shx.shReadGFC("test_variable.gfct");   % struct incl. variableTerms
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
idxF = shx.shIndex(L);
ws.Ac = exp(-idxF.n / 3) .* randn(idxF.P, 1);
[ws.tsF, ws.op, ws.info] = ws.ts.filter("tvANS", Blocks = "off");
ws.B = randn(idxF.P, 2);
ws.b = ws.B(:, 1);
[ws.avg, ws.out] = shx.basinDeconvolve(ws.B, ws.op);
ws.Y3 = [ws.avg(1, :)', ws.avg(1, :)' + 2e-4 * randn(T, 1), ...
    ws.avg(1, :)' + 5e-4 * randn(T, 1)];   % 3-center stack for TCH
ws.tYears = tY;
XV = zeros(ws.idx.P, T);                 % idx-ordered raw stack (L = 10)
for t = 1:T
    XV(:, t) = shx.vecFromCS(Cs(:,:,t), Ss(:,:,t), ws.idx);
end
ws.X = XV;
[~, ws.Xres] = shx.fitDeterministicModel(XV, tY);   % idx-ordered residuals
ws.N = shx.buildNoiseCov(ws.Xres, ws.idx);
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
h = shx.taylorDiagram(1.0, [0.8, 1.1], [0.9, 0.95], ...
    Labels = ["a", "b"], Normalize = true);
verifyTrue(testCase, isgraphics(h));
d = fullfile(fileparts(mfilename('fullpath')), 'test_data');
fG = fullfile(d, 'ITSG-Grace2018_n60_2008-04.gfc');
assumeTrue(testCase, isfile(fG));
g = shCoefficients.read(fG, Epoch = 2008.29);
[~, h1] = shx.compareSolutions(g, g.gaussian(500), Plot = true, ...
    LatDeg = -85:5:85, LonDeg = 0:9:351);
verifyTrue(testCase, isgraphics(h1));
rng(11); n1 = 7; T = 24;
Cs = 1e-9 * randn(n1, n1, T); Ss = 1e-9 * randn(n1, n1, T);
mk = @(s) shSeries(Cs + s * 1e-11 * randn(n1, n1, T), ...
    Ss = Ss + s * 1e-11 * randn(n1, n1, T), Epochs = 2019 + (0:T-1)'/12);
[~, h2] = shx.compareSeries({mk(1), mk(2), mk(3)}, Plot = true, ...
    LatDeg = -80:20:80, LonDeg = 0:30:330);
verifyTrue(testCase, isgraphics(h2));
close all
end

function testVersionMetadata(testCase)
% shx.version reports toolbox metadata parsed from Contents.m - the
% single source of truth also honoured by MATLAB's ver().
v = shx.version();
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
verifyError(testCase, @() shx.fetchITSG("2010-4"), 'shx:fetchITSG:badMonth');
verifyError(testCase, @() shx.fetchITSG("April 2010"), 'shx:fetchITSG:badMonth');
verifyError(testCase, @() shx.fetchITSG(1990), 'shx:fetchITSG:badMonth');
% daily product (v2.4.1): n40 only; monthly rejects n40
verifyError(testCase, ...
    @() shx.fetchITSG("2010-04", Product = "daily", Nmax = 96), ...
    'shx:fetchITSG:badNmax');
verifyError(testCase, @() shx.fetchITSG("2010-04", Nmax = 40), ...
    'shx:fetchITSG:badNmax');
end

function testSpectrumPlotOptionsAndAxes(testCase)
rng(102);
fig = figure('Visible', 'off');
cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
L = 15; n1 = L + 1;
C = tril(randn(n1)) * 1e-9; S = tril(randn(n1), -1) * 1e-9; S(:, 1) = 0;
sd = shx.shDegreeRMS(C, S, 'sigmaC', abs(C)*0.1, 'sigmaS', abs(S)*0.1);
so = shx.shOrderRMS(C, S);
ax = axes('Parent', fig);
% linear x-axis contract (v2.4.1)
plotSHSpectrum(sd, 'ax', ax);
verifyEqual(testCase, ax.XScale, 'linear');
verifyEqual(testCase, ax.YScale, 'log');
% every quantity renders, degree and order domain
for q = ["amplitude", "rms", "variance", "cumamplitude", "cumrms", "cumvariance"]
    cla(ax); plotSHSpectrum(sd, 'ax', ax, 'Quantity', char(q));
    cla(ax); plotSHSpectrum(so, 'ax', ax, 'Quantity', char(q));
end
% contracts
verifyError(testCase, ...
    @() plotSHSpectrum({sd, so}, 'ax', axes('Parent', fig)), ...
    'shx:plotSHSpectrum:mixedDomains');
verifyError(testCase, ...
    @() plotSHSpectrum(sd, 'ax', axes('Parent', fig), ...
    'Quantity', 'variance', 'Kaula', 1e-5), ...
    'shx:plotSHSpectrum:kaulaDomain');
verifyError(testCase, ...
    @() plotSHSpectrum(sd, 'ax', axes('Parent', fig), 'Quantity', 'xxx'), ...
    'shx:plotSHSpectrum:badQuantity');
% triangle: no center line anymore (v2.4.1)
ax2 = axes('Parent', fig);
plotSHCoeffTriangle(C, S, 'ax', ax2);
verifyEmpty(testCase, findobj(ax2, 'Type', 'constantline'));
end

function testDataFolderAndDDKNames(testCase)
% preserve any user preference across the test
had = ispref('shAnalysis', 'dataFolder');
if had
    old = getpref('shAnalysis', 'dataFolder');
    restore = onCleanup(@() setpref('shAnalysis', 'dataFolder', old)); %#ok<NASGU>
else
    restore = onCleanup(@() shx.dataFolder("reset")); %#ok<NASGU>
end
tmp = fullfile(tempdir, sprintf('shx_data_%d', randi(1e9)));
f = shx.dataFolder(tmp);
verifyEqual(testCase, char(f), char(string(tmp)));
verifyTrue(testCase, isfolder(f));
verifyEqual(testCase, char(shx.dataFolder()), char(string(tmp)));
f2 = shx.dataFolder("reset");
verifyTrue(testCase, endsWith(f2, "data"));
% DDK mapping: 8 unique released names, DDK3 matches the shipped file
nm = shx.ddkNames();
verifyEqual(testCase, numel(nm), 8);
verifyEqual(testCase, numel(unique(nm)), 8);
verifyEqual(testCase, nm(3), "Wbd_2-120.a_1d12p_4");
% name resolution: DDK3 loads from shipped test data even with a fresh
% (empty) data folder; an unfetched filter errors with the fetch hint
shx.dataFolder(tmp);
W = shx.readDDK("DDK3", Nmax = 30);
verifyEqual(testCase, W.nmax, 30);
verifyError(testCase, @() shx.readDDK("DDK7"), 'shx:readDDK:notFetched');
verifyError(testCase, @() shx.fetchDDK(0), 'MATLAB:validators:mustBeInRange');
end

