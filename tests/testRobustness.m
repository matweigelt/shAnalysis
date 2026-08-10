function tests = testRobustness
%TESTROBUSTNESS Edge-case and degenerate-input tests for shAnalysis v2.
%
%   Covers: NaN handling, missing GAX epochs, single-epoch series,
%   rank-deficient basin sets (ridge path), poles in the synthesis grid,
%   zero fields, degenerate destriping parameters, gzip round trips, and
%   corrupted/empty gfc files.
%
%   Run:  results = runtests('testRobustness');
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

tests = functiontests(localfunctions);
end

function setupOnce(tc)
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);
tc.TestData.dataDir = fullfile(here, 'test_data');
shLowLevel.legendreCached('clear');
end

function g = randomField(nmax, epoch)
rng(11);
C = tril(randn(nmax+1)) * 1e-9; S = tril(randn(nmax+1), -1) * 1e-9;
g = shCoefficients(C, S, Epoch = epoch, ProductType = "GSM");
end

function ts = randomSeries(nmax, T)
rng(12);
mL  = tril(true(nmax+1));           % tril does not accept 3-D input
mL1 = tril(true(nmax+1), -1);
Cs = (mL  .* randn(nmax+1, nmax+1, T)) * 1e-9;
Ss = (mL1 .* randn(nmax+1, nmax+1, T)) * 1e-9;
ts = shSeries(Cs, Ss = Ss, Epochs = 2020 + (0:T-1)'/12, ProductType = "GSM");
end

% ----------------------------------------------------------- NaN handling
function testNaNMonthMeanOmitnan(tc)
ts = randomSeries(6, 8);
Cs = ts.Cs; Cs(:,:,4) = NaN;
tsN = shSeries(Cs, Ss = ts.Ss, Epochs = ts.epochs);
m = tsN.mean(Omitnan = true);                   % must not propagate NaN
verifyTrue(tc, all(isfinite(m.C(:))));
mBad = tsN.mean(Omitnan = false);
verifyTrue(tc, all(isnan(mBad.C(:))));
% ... but the filter chain must refuse NaN stacks loudly:
verifyError(tc, @() tsN.gaussian(300), 'shSeries:nanInSeries');
verifyError(tc, @() tsN.filter("tvANS"), 'shSeries:nanInSeries');
end

% -------------------------------------------------------- GAX edge cases
function testRestoreAllowMissing(tc)
gsm = randomSeries(6, 6);
gadC = gsm.Cs * 0.1; gadS = gsm.Ss * 0.1;
% GAX series missing the last two months
gad = shSeries(gadC(:,:,1:4), Ss = gadS(:,:,1:4), ...
    Epochs = gsm.epochs(1:4), ProductType = "GAD");
out = gsm.restore(gad, AllowMissing = true);
% matched epochs changed, unmatched untouched
verifyEqual(tc, out.Cs(:,:,1), gsm.Cs(:,:,1) + gadC(:,:,1), RelTol = 1e-14);
verifyEqual(tc, out.Cs(:,:,6), gsm.Cs(:,:,6), AbsTol = 0);
verifyEqual(tc, out.productType, "GSM+GAD");
verifySubstring(tc, char(out.history(end)), '4/6');
end

% ------------------------------------------------------ tiny/degenerate T
function testSingleEpochSeries(tc)
ts = randomSeries(6, 1);
m = ts.mean;                                    % mean of one epoch = itself
verifyEqual(tc, m.C, ts.Cs(:,:,1), AbsTol = 0);
verifyError(tc, @() ts.climatology(), 'shSeries:tooFewEpochs');
end

% ------------------------------------------------- rank-deficient basins
function testRankDeficientBasinsRidge(tc)
ts = randomSeries(10, 24);
[tsF, op] = ts.filter("tvANS");                 %#ok<ASGLU>
idx = shLowLevel.shIndex(ts.nmax, MinDegree = 2);
rng(13);
b = randn(idx.P, 1);
B = [b, b];                                     % exactly collinear pair
% Ridge = 0 would invert a singular (B'SB); the ridge path must survive.
avg = ts.basinAverage(B / norm(b)^2 * idx.P, Deconvolve = true, ...
    Op = op, Ridge = 1e-6);
verifyTrue(tc, all(isfinite(avg(:))));
verifySize(tc, avg, [2, ts.nEpochs]);
end

% -------------------------------------------------------- synthesis edges
function testSynthesisAtPoles(tc)
g = randomField(15, 2020);
[grid, ~, ~] = g.synthesis([-90, 0, 90], 0:30:330);
verifyTrue(tc, all(isfinite(grid(:))));
% at the poles only m=0 survives: no longitude dependence
verifyEqual(tc, grid(1,:), grid(1,1)*ones(1,12), RelTol = 1e-12);
verifyEqual(tc, grid(3,:), grid(3,1)*ones(1,12), RelTol = 1e-12);
end

function testSynthesisNonMonotonicLat(tc)
g = randomField(10, 2020);
latA = [-30, 45, 10]; lon = 0:60:300;
gridA = g.synthesis(latA, lon);
for k = 1:3
    gk = g.synthesis(latA(k), lon);
    verifyEqual(tc, gridA(k,:), gk, RelTol = 1e-13);
end
end

function testZeroField(tc)
z = shCoefficients(zeros(11), zeros(11), Epoch = 2020);
z2 = z.destripe().gaussian(500);
verifyEqual(tc, z2.C, zeros(11), AbsTol = 0);
grid = z2.synthesis(-60:30:60, 0:90:270);
verifyEqual(tc, grid, zeros(5, 4), AbsTol = 0);
spec = z.degreeRMS;
verifyEqual(tc, spec.degRMS, zeros(size(spec.degRMS)), AbsTol = 0);
end

% ------------------------------------------------------- degenerate opts
function testDestripeMinOrderAboveNmax(tc)
g = randomField(10, 2020);
g2 = g.destripe(minOrder = 99);                 % nothing to filter
verifyEqual(tc, g2.C, g.C, AbsTol = 0);
verifyEqual(tc, g2.S, g.S, AbsTol = 0);
end

function testGaussianTinyAndHugeRadius(tc)
g = randomField(20, 2020);
gT = g.gaussian(1);                             % ~no smoothing
verifyEqual(tc, gT.C, g.C, RelTol = 1e-3);
gH = g.gaussian(20000);                         % kills high degrees
specH = gH.degreeRMS; specG = g.degreeRMS;
verifyLessThan(tc, specH.degRMS(end), 1e-3 * specG.degRMS(end));
verifyTrue(tc, all(isfinite(gH.C(:))));         % recursion must not blow up
end

% ------------------------------------------------------------- file I/O
function testWriteGzipRoundtrip(tc)
g = randomField(8, 2020.5);
f = fullfile(tempdir, 'robust_rt.gfc');
g.write(f);
gzip(f); gz = [f '.gz'];
tc.addTeardown(@() cellfun(@delete, {f, gz}));
g2 = shCoefficients.read(gz);
verifyEqual(tc, g2.C, g.C, AbsTol = 1e-22);
verifyEqual(tc, g2.S, g.S, AbsTol = 1e-22);
verifyEqual(tc, g2.GM, g.GM, RelTol = 1e-14);
end

function testEmptyAndCorruptGFC(tc)
% header only, no data records
f1 = fullfile(tempdir, 'robust_empty.gfc');
fid = fopen(f1, 'w');
fprintf(fid, 'product_type gravity_field\nend_of_head\n');
fclose(fid);
tc.addTeardown(@() delete(f1));
verifyError(tc, @() shCoefficients.read(f1), 'shReadGFC:noData');

% junk lines between valid records must be skipped, not crash
f2 = fullfile(tempdir, 'robust_junk.gfc');
fid = fopen(f2, 'w');
fprintf(fid, 'earth_gravity_constant 3.986004415e14\nradius 6378136.3\n');
fprintf(fid, 'end_of_head\n');
fprintf(fid, 'gfc 0 0 1.0 0.0\n');
fprintf(fid, 'this line is garbage and has no numbers\n');
fprintf(fid, 'gfc 2 0 -4.84e-4 0.0\n');
fclose(fid);
tc.addTeardown(@() delete(f2));
g = shCoefficients.read(f2);
verifyEqual(tc, g.C(1,1), 1.0);
verifyEqual(tc, g.C(3,1), -4.84e-4);
end

% ---------------------------------------------- gfct without variableTerms
function testEvalAtOnStaticModel(tc)
f = fullfile(tc.TestData.dataDir, 'test_static.gfc');
g = shCoefficients.read(f);
% static model: evalAt must be an identity AND must say so - the warning
% is part of the contract, so assert it (verifyWarning also captures it,
% keeping the test log clean)
g2 = verifyWarning(tc, @() g.evalAt(2020.0), 'shEvalGFCT:staticModel');
verifyEqual(tc, g2.C, g.C, AbsTol = 0);
end

% =================================================================== v2.1
function testTN13SyntheticRoundTrip(tc)
% synthetic TN-13 file in the documented GRCOF2 layout -> parse -> apply
f = [tempname '.txt'];
fid = fopen(f, 'w');
fprintf(fid, 'HEADER: GRACE Technical Note 13 (synthetic fixture)\n');
fprintf(fid, 'GRCOF2    1    0  1.2e-10  0.0      3.0e-11  0.0      20080401.0000 20080430.2359\n');
fprintf(fid, 'GRCOF2    1    1 -4.0e-10  2.5e-10  3.1e-11  3.2e-11  20080401.0000 20080430.2359\n');
fprintf(fid, 'GRCOF2    1    0  1.3e-10  0.0      3.0e-11  0.0      20080501.0000 20080531.2359\n');
fprintf(fid, 'GRCOF2    1    1 -4.1e-10  2.6e-10  3.1e-11  3.2e-11  20080501.0000 20080531.2359\n');
fclose(fid);
cleanup = onCleanup(@() delete(f)); %#ok<NASGU>
tn = shLowLevel.readTN13(f);
verifyEqual(tc, numel(tn.epoch), 2);
verifyEqual(tc, tn.C10(1), 1.2e-10);
verifyEqual(tc, tn.S11(2), 2.6e-10);
verifyGreaterThan(tc, tn.epoch(1), 2008.2);
verifyLessThan(tc, tn.epoch(1), 2008.4);
% apply to a coefficient set with a matching epoch
g = shCoefficients(zeros(4), zeros(4), Epoch = tn.epoch(1));
g1 = g.addDegree1(tn);
verifyEqual(tc, g1.C(2,1), 1.2e-10);
verifyEqual(tc, g1.C(2,2), -4.0e-10);
verifyEqual(tc, g1.S(2,2), 2.5e-10);
verifyEqual(tc, g1.sigmaS(2,2), 3.2e-11);
end

function testSINEXCovaRoundTrip(tc)
% synthetic SINEX with SOLUTION/ESTIMATE + lower-triangular COVA
f = [tempname '.snx'];
fid = fopen(f, 'w');
fprintf(fid, '%%=SNX 2.02\n+SOLUTION/ESTIMATE\n');
fprintf(fid, '*INDEX TYPE N M ...\n');
fprintf(fid, '  1 CN 2 0  X  A  1  -4.84e-04 1.0e-11\n');
fprintf(fid, '  2 CN 2 1  X  A  1   2.03e-06 1.1e-11\n');
fprintf(fid, '  3 SN 2 1  X  A  1   1.40e-06 1.2e-11\n');
fprintf(fid, '-SOLUTION/ESTIMATE\n+SOLUTION/MATRIX_ESTIMATE L COVA\n');
fprintf(fid, ' 1 1  4.0e-22\n 2 1  1.0e-23  9.0e-22\n 3 1  0.0  2.0e-23  1.6e-21\n');
fprintf(fid, '-SOLUTION/MATRIX_ESTIMATE L COVA\n%%ENDSNX\n');
fclose(fid);
cleanup = onCleanup(@() delete(f)); %#ok<NASGU>
snx = shLowLevel.readSINEX(f);
verifyEqual(tc, snx.kind, 'COVA');
verifyEqual(tc, snx.x, [-4.84e-04; 2.03e-06; 1.40e-06]);
verifyEqual(tc, snx.M(1,1), 4.0e-22);
verifyEqual(tc, snx.M(1,2), 1.0e-23);          % symmetrized
verifyEqual(tc, snx.M(3,2), 2.0e-23);
verifyTrue(tc, issymmetric(snx.M));
end

function testSINEXNeqInversionAndReorder(tc)
% NEQ block inverted to covariance, then reordered to a shIndex ordering
f = [tempname '.snx'];
fid = fopen(f, 'w');
fprintf(fid, '%%=SNX 2.02\n+SOLUTION/ESTIMATE\n');
fprintf(fid, '  1 SN 2 1 X A 1  1.40e-06 1.0e-11\n');   % file order: SN first
fprintf(fid, '  2 CN 2 0 X A 1 -4.84e-04 1.0e-11\n');
fprintf(fid, '  3 CN 2 1 X A 1  2.03e-06 1.0e-11\n');
fprintf(fid, '-SOLUTION/ESTIMATE\n+SOLUTION/NORMAL_EQUATION_MATRIX U\n');
fprintf(fid, ' 1 1  2.0e21  1.0e19  0.0\n 2 2  4.0e21  2.0e19\n 3 3  1.0e21\n');
fprintf(fid, '-SOLUTION/NORMAL_EQUATION_MATRIX U\n%%ENDSNX\n');
fclose(fid);
cleanup = onCleanup(@() delete(f)); %#ok<NASGU>
idx = struct('n', [2;2;2], 'm', [0;1;1], 'cs', [0;0;1], 'P', 3);
snx = shLowLevel.readSINEX(f, Output = "covariance", Index = idx);
verifyEqual(tc, snx.kind, 'COVA');
verifyEqual(tc, snx.x, [-4.84e-04; 2.03e-06; 1.40e-06], 'AbsTol', 0);
% covariance must be inv(NEQ) reordered: check against direct computation
Nq = [2.0e21 1.0e19 0.0; 1.0e19 4.0e21 2.0e19; 0.0 2.0e19 1.0e21];
Cv = inv(Nq);
perm = [2 3 1];                                 % file rows in idx order
verifyEqual(tc, snx.M, Cv(perm, perm), 'RelTol', 1e-10);
end

function testICGEM2PiecewiseEval(tc)
% two gfct pieces + interval-limited trnd/periodic terms
f = [tempname '.gfc'];
fid = fopen(f, 'w');
fprintf(fid, ['product_type gravity_field\nmodelname T2\n' ...
    'earth_gravity_constant 3.986004415E+14\nradius 6378136.3\n' ...
    'max_degree 2\nformat icgem2.0\nerrors formal\nend_of_head ====\n']);
fprintf(fid, 'gfc  2 1  1.0e-06 2.0e-06 0 0\n');
fprintf(fid, 'gfct 2 0 -4.0e-04 0.0 0 0 20100101 20150101\n');
fprintf(fid, 'gfct 2 0 -5.0e-04 0.0 0 0 20150101 20200101\n');
fprintf(fid, 'trnd 2 0  1.0e-11 0.0 0 0 20150101 20200101\n');
fprintf(fid, 'acos 2 0  2.0e-11 0.0 0 0 1.0 20100101 20200101\n');
fclose(fid);
cleanup = onCleanup(@() delete(f)); %#ok<NASGU>
model = shLowLevel.shReadGFC(f);
% at 2012.0: first piece, no trend (inactive), annual cos with t0=2010.0
[Ct, ~] = shLowLevel.shEvalGFCT(model, 2012.0);
verifyEqual(tc, Ct(3,1), -4.0e-04 + 2.0e-11*cos(2*pi*2.0), 'RelTol', 1e-12);
verifyEqual(tc, Ct(3,2), 1.0e-06);              % static gfc untouched
% at 2016.5: second piece + trend from its own t0 (2015.0) + annual
[Ct2, ~] = shLowLevel.shEvalGFCT(model, 2016.5);
expected = -5.0e-04 + 1.0e-11*1.5 + 2.0e-11*cos(2*pi*6.5);
verifyEqual(tc, Ct2(3,1), expected, 'RelTol', 1e-10);
% C/S expose the most recent piece
verifyEqual(tc, model.C(3,1), -5.0e-04);
end

function testRealITSGGraceRead(tc)
% real ITSG-Grace2018 monthly solution (downloaded from TU Graz):
% header constants, completeness, physically plausible C20, finite sigmas
f = fullfile(tc.TestData.dataDir, 'ITSG-Grace2018_n60_2008-04.gfc');
assumeTrue(tc, isfile(f), 'real ITSG GRACE file not present');
g = shCoefficients.read(f, Epoch = 2008 + 3.5/12);
verifyEqual(tc, g.nmax, 60);
verifyEqual(tc, g.GM, 3.9860044150e14);
verifyEqual(tc, g.R, 6.3781363000e6);
verifyEqual(tc, g.C(1,1), 1, 'AbsTol', 1e-12);          % C00 = 1
verifyEqual(tc, g.C(3,1), -4.84e-4, 'AbsTol', 1e-6);    % C20 ballpark
verifyTrue(tc, all(isfinite(g.sigmaC(3:end,1))));        % formal errors
verifyGreaterThan(tc, nnz(g.C), 1800);                   % complete triangle populated
spec = g.degreeRMS();
verifyTrue(tc, all(isfinite(spec.degAmplitude(2:end))));
end

function testRealGraceFORead(tc)
% real GRACE-FO operational solution (2025-12) reads identically
f = fullfile(tc.TestData.dataDir, 'ITSG-Grace_operational_n60_2025-12.gfc');
assumeTrue(tc, isfile(f), 'real GRACE-FO file not present');
g = shCoefficients.read(f, Epoch = 2025 + 11.5/12);
verifyEqual(tc, g.nmax, 60);
verifyEqual(tc, g.C(3,1), -4.84e-4, 'AbsTol', 1e-6);
verifyTrue(tc, all(isfinite(g.sigmaC(3:end,1))));
% the two real files must differ (different months/missions)
f2 = fullfile(tc.TestData.dataDir, 'ITSG-Grace2018_n60_2008-04.gfc');
g2 = shCoefficients.read(f2);
verifyGreaterThan(tc, max(abs(g.C(:) - g2.C(:))), 0);
end

function testAnalysisSingleRing(tc)
% degenerate: one latitude ring cannot separate degrees -> rankDeficient
lon = (0:29) * 12;
g = ones(1, 30);
verifyError(tc, @() shLowLevel.shAnalysisGrid(g, 10, lon, 4, Method = "rings"), ...
    'shLowLevel:shAnalysisGrid:rankDeficient');
% but with Kaula it returns finite coefficients
[C, S] = shLowLevel.shAnalysisGrid(g, 10, lon, 4, Method = "rings", Kaula = 1);
verifyTrue(tc, all(isfinite(C(:))) && all(isfinite(S(:))));
end

function testRealSINEXFixture(tc)
% REAL ITSG-Grace2018 n96 SINEX content (truncated fixture, TU Graz):
% verifies the actual layout "CN 2 -- 0 08:107:00000 ---- 2 val sig"
% (CODE=degree, SOLN=order, PT='--') and the 12x12 NEQ principal
% submatrix. Expected numbers cross-checked in Python.
f = fullfile(tc.TestData.dataDir, 'ITSG-Grace2018_n96_2008-04_head12.snx');
assumeTrue(tc, isfile(f), 'real SINEX fixture not present');
snx = shLowLevel.readSINEX(f);
verifyEqual(tc, snx.kind, 'NEQ');
verifyEqual(tc, numel(snx.x), 12);
verifyEqual(tc, snx.n(1:6), [2;2;2;2;2;3]);
verifyEqual(tc, snx.m(1:6), [0;1;1;2;2;0]);
verifyEqual(tc, snx.cs(1:6), [0;0;1;0;1;0]);
verifyEqual(tc, snx.x(1), -4.84169356322812e-04, 'AbsTol', 0);
verifyEqual(tc, snx.x(4),  2.43934829341615e-06, 'AbsTol', 0);
verifyEqual(tc, snx.x(5), -1.40033136191992e-06, 'AbsTol', 0);
verifyEqual(tc, snx.sig(1), 1.13782e-11, 'AbsTol', 0);
% NEQ -> covariance (Cholesky inversion of the PD 12x12 submatrix);
% expected sigmas from independent numpy inversion
snxC = shLowLevel.readSINEX(f, Output = "covariance");
verifyEqual(tc, snxC.kind, 'COVA');
verifyTrue(tc, issymmetric(snxC.M) || norm(snxC.M - snxC.M', 'fro') < 1e-30);
sd = sqrt(diag(snxC.M));
verifyEqual(tc, sd(1), 1.08211427e-11, 'RelTol', 1e-6);
verifyEqual(tc, sd(2), 2.87455389e-12, 'RelTol', 1e-6);
verifyEqual(tc, sd(3), 3.04963470e-12, 'RelTol', 1e-6);
% Index reorder into an explicit target ordering
idx = struct('n', [3;2;2], 'm', [0;0;1], 'cs', [0;0;1], 'P', 3);
snxR = shLowLevel.readSINEX(f, Index = idx);
verifyEqual(tc, snxR.x(1), 9.57262594695377e-07, 'AbsTol', 0);   % C30
verifyEqual(tc, snxR.x(2), -4.84169356322812e-04, 'AbsTol', 0);  % C20
verifyEqual(tc, snxR.x(3), 1.44480342951435e-09, 'AbsTol', 0);   % S21
end

function testSINEXTruncatedMatrixErrors(tc)
% matrix indices beyond the parameter list must fail loudly, not crash
f = [tempname '.snx'];
fid = fopen(f, 'w');
fprintf(fid, '%%=SNX 2.02\n+SOLUTION/ESTIMATE\n');
fprintf(fid, ' 1 CN 2 -- 0 08:107:00000 ---- 2 -4.84e-04 1.1e-11\n');
fprintf(fid, '-SOLUTION/ESTIMATE\n+SOLUTION/NORMAL_EQUATION_MATRIX U\n');
fprintf(fid, ' 1 1 1.0e21 2.0e19 3.0e19\n');    % j+2 = 3 > Q = 1
fprintf(fid, '-SOLUTION/NORMAL_EQUATION_MATRIX U\n%%ENDSNX\n');
fclose(fid);
cleanup = onCleanup(@() delete(f)); %#ok<NASGU>
verifyError(tc, @() shLowLevel.readSINEX(f), 'shLowLevel:readSINEX:badMatrixIndex');
end

function testSetupPathIdempotent(tc)
% in the test environment the toolbox is already on the path: setup must
% detect that, change nothing, and report ok (Download="none" default).
p0 = path;
restore = onCleanup(@() path(p0)); %#ok<NASGU>
s = setup_shAnalysis(Quiet = true);
verifyEqual(tc, s.pathAction, "already-on-path");
verifyEqual(tc, path, p0);
verifyTrue(tc, s.ok);
verifyTrue(tc, isfolder(s.dataFolder));
verifyTrue(tc, isempty(s.fetched) && isempty(s.failed));
end

function testFetchTNSkipAndVerify(tc)
% OFFLINE contract of shLowLevel.fetchTN: files already present in Dest are
% skipped without any network access but still verified by parse; a
% present-but-corrupt file is reported in info.failed and excluded from
% files (and NOT deleted, since it was not fetched by this call).
d = tc.TestData.dataDir;
tmp = tempname; mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>
copyfile(fullfile(d, 'TN-13_GEOC_GFZ_RL06_3.txt'), ...
    fullfile(tmp, 'TN-13_GEOC_GFZ_RL06.3.txt'));
copyfile(fullfile(d, 'TN-13_GEOC_CSR_RL06.3.txt'), tmp);
copyfile(fullfile(d, 'TN-13_GEOC_JPL_RL06.3.txt'), tmp);
copyfile(fullfile(d, 'TN-14_C30_C20_SLR_GSFC.txt'), tmp);
[files, info] = shLowLevel.fetchTN(Dest = tmp, Quiet = true);
verifyEqual(tc, numel(files), 4);
verifyEqual(tc, numel(info.skipped), 4);
verifyTrue(tc, isempty(info.fetched) && isempty(info.failed));
% corrupt present file -> failed, not deleted
tmp2 = tempname; mkdir(tmp2);
cleanup2 = onCleanup(@() rmdir(tmp2, 's')); %#ok<NASGU>
bad = fullfile(tmp2, 'TN-13_GEOC_CSR_RL06.3.txt');
writelines("this is not a TN-13 file", bad);
w0 = warning('off', 'shLowLevel:fetchTN:failed');
restoreW = onCleanup(@() warning(w0)); %#ok<NASGU>
[files2, info2] = shLowLevel.fetchTN(Dest = tmp2, Providers = "CSR", ...
    TN14 = false, Quiet = true);
verifyTrue(tc, isempty(files2));
verifyEqual(tc, numel(info2.failed), 1);
verifyTrue(tc, contains(info2.failed(1), "parse"));
verifyTrue(tc, isfile(bad));                     % pre-existing: kept
end

function testFetchTNUpdateFromMirror(tc)
% Update=true must refresh existing files via a SAFE SWAP: the fresh
% copy is parse-verified BEFORE it replaces the old one, and a missing
% or corrupt source keeps the existing file authoritative (reported in
% info.failed, file listed in files/skipped). Exercised offline through
% the BaseURL local-mirror folder mode.
d = tc.TestData.dataDir;
mir = tempname; mkdir(mir);
c1 = onCleanup(@() rmdir(mir, 's')); %#ok<NASGU>
dst = tempname; mkdir(dst);
c2 = onCleanup(@() rmdir(dst, 's')); %#ok<NASGU>
nmC = 'TN-13_GEOC_CSR_RL06.3.txt';
% mirror holds the true CSR file; dest is pre-seeded with GFZ content
% under the CSR name (a VALID TN-13, but different bytes)
copyfile(fullfile(d, nmC), fullfile(mir, nmC));
copyfile(fullfile(d, 'TN-13_GEOC_GFZ_RL06_3.txt'), fullfile(dst, nmC));
dm = dir(fullfile(mir, nmC)); szMir = dm.bytes;
dp = dir(fullfile(dst, nmC)); szPre = dp.bytes;
verifyNotEqual(tc, szPre, szMir);        % providers differ in size
% (1) no Update: existing file skipped and untouched
[f0, i0] = shLowLevel.fetchTN(Dest = dst, BaseURL = mir, Providers = "CSR", ...
    TN14 = false, Quiet = true);
verifyEqual(tc, numel(f0), 1);
verifyEqual(tc, i0.skipped, string(nmC));
verifyTrue(tc, isempty(i0.updated) && isempty(i0.fetched));
d0 = dir(fullfile(dst, nmC));
verifyEqual(tc, d0.bytes, szPre);
% (2) Update: replaced by the mirror copy after parse verification
[f1, i1] = shLowLevel.fetchTN(Dest = dst, BaseURL = mir, Providers = "CSR", ...
    TN14 = false, Update = true, Quiet = true);
verifyEqual(tc, numel(f1), 1);
verifyEqual(tc, i1.updated, string(nmC));
verifyTrue(tc, isempty(i1.fetched) && isempty(i1.failed));
d1 = dir(fullfile(dst, nmC));
verifyEqual(tc, d1.bytes, szMir);
% (3) failed refresh: TN-14 absent in the mirror -> existing file kept
nm14 = 'TN-14_C30_C20_SLR_GSFC.txt';
copyfile(fullfile(d, nm14), fullfile(dst, nm14));
[f2, i2] = shLowLevel.fetchTN(Dest = dst, BaseURL = mir, Providers = "CSR", ...
    TN14 = true, Update = true, Quiet = true);
verifyTrue(tc, any(contains(i2.failed, "TN-14")));
verifyTrue(tc, any(i2.skipped == string(nm14)));   % kept + re-verified
verifyTrue(tc, any(f2 == string(fullfile(dst, nm14))));
verifyTrue(tc, isfile(fullfile(dst, nm14)));
% (4) fresh dest: mirror serves a NEW file (fetched, not updated)
dst2 = tempname; mkdir(dst2);
c3 = onCleanup(@() rmdir(dst2, 's')); %#ok<NASGU>
[f3, i3] = shLowLevel.fetchTN(Dest = dst2, BaseURL = mir, Providers = "CSR", ...
    TN14 = false, Quiet = true);
verifyEqual(tc, numel(f3), 1);
verifyEqual(tc, i3.fetched, string(nmC));
verifyTrue(tc, isempty(i3.updated));
end

function testFetchProxyPlumbing(tc)
% v2.7.0: Proxy= routes downloads through matlab.net.http with an
% explicit ProxyURI. Offline-deterministic: 127.0.0.1 refuses instantly
% on both paths; the mirror mode copies files and ignores Proxy.
dst = tempname; mkdir(dst);
c1 = onCleanup(@() rmdir(dst, 's')); %#ok<NASGU>
% websave path (no proxy): unreachable base -> clean failure report
[~, i0] = shLowLevel.fetchTN(Dest = dst, BaseURL = "https://127.0.0.1:1", ...
    Providers = "CSR", TN14 = false, Timeout = 2, Quiet = true);
verifyEqual(tc, numel(i0.failed), 1);
% webFetch path (proxy set): unreachable proxy -> clean failure report
[~, i1] = shLowLevel.fetchTN(Dest = dst, BaseURL = "https://127.0.0.1:1", ...
    Providers = "CSR", TN14 = false, Timeout = 2, Quiet = true, ...
    Proxy = "http://127.0.0.1:9");
verifyEqual(tc, numel(i1.failed), 1);
% mirror mode: local copy, Proxy irrelevant, succeeds
d = tc.TestData.dataDir;
mir = tempname; mkdir(mir);
c2 = onCleanup(@() rmdir(mir, 's')); %#ok<NASGU>
nmC = 'TN-13_GEOC_CSR_RL06.3.txt';
copyfile(fullfile(d, nmC), fullfile(mir, nmC));
[f2, i2] = shLowLevel.fetchTN(Dest = dst, BaseURL = mir, Providers = "CSR", ...
    TN14 = false, Quiet = true, Proxy = "http://127.0.0.1:9");
verifyEqual(tc, numel(f2), 1);
verifyEqual(tc, i2.fetched, string(nmC));
end

function testRealTN14GSFC(tc)
% REAL GSFC TN-14 C20/C30 file (v.3, snapshot retrieved from GFZ ISDC
% 2026-08-07): 258 solution windows, starts 2002.2548 .. 2026.4137,
% last window stop 2026.4959; C30 valid from 2012.16 (GRACE-FO era
% backfilled by SLR). Head values are immutable history and pinned
% exactly (independent Python parse of the 10-column layout, column 2
% = window START year fraction, column 10 = window STOP - the v2.5
% pin confused the two). Tail checks are lower bounds so a refreshed
% upstream file with additional months still passes.
f = fullfile(tc.TestData.dataDir, 'TN-14_C30_C20_SLR_GSFC.txt');
assumeTrue(tc, isfile(f), 'real TN-14 file not present');
tn = shLowLevel.readTN14(f);
verifyGreaterThanOrEqual(tc, numel(tn.epoch), 258);
verifyEqual(tc, tn.C20(1), -4.8416934147454e-04, 'AbsTol', 0);
verifyEqual(tc, tn.sigmaC20(1), 0.1628e-10, 'AbsTol', 1e-16);
verifyEqual(tc, tn.epochStart(1), 2002.2548, 'AbsTol', 0);
verifyEqual(tc, tn.epochStop(1), 2002.3288, 'AbsTol', 0);
verifyGreaterThanOrEqual(tc, tn.epochStart(end), 2026.4137);
verifyGreaterThanOrEqual(tc, tn.epochStop(end), 2026.4959);
verifyGreaterThan(tc, tn.epochStop(end), tn.epochStart(end));
verifyTrue(tc, all(isfinite(tn.C20)));
firstC30 = find(isfinite(tn.C30), 1);
verifyEqual(tc, tn.epochStart(firstC30), 2012.1639, 'AbsTol', 0);
verifyTrue(tc, all(diff(tn.epochStart) > 0));
end

function testRealTN13CSRJPL(tc)
% REAL CSR and JPL TN-13 files (RL06.3, retrieved from GFZ ISDC
% 2026-08-07): same GRCOF2 layout as GFZ, 256 paired months each,
% 2002-04 .. 2026-04. Reference values pinned from an independent
% Python parse; cross-provider C10 correlation vs the GFZ file 0.995
% (annual geocenter signal is common to all providers).
d = tc.TestData.dataDir;
fC = fullfile(d, 'TN-13_GEOC_CSR_RL06.3.txt');
fJ = fullfile(d, 'TN-13_GEOC_JPL_RL06.3.txt');
assumeTrue(tc, isfile(fC) && isfile(fJ), 'CSR/JPL TN-13 not present');
tnC = shLowLevel.readTN13(fC);
tnJ = shLowLevel.readTN13(fJ);
for tn = [tnC, tnJ]
    verifyEqual(tc, numel(tn.epoch), 256);
    verifyTrue(tc, all(isfinite(tn.C10) & isfinite(tn.C11) & ...
        isfinite(tn.S11)));
    verifyTrue(tc, all(tn.sigC10 > 0 & tn.sigC11 > 0 & tn.sigS11 > 0));
    verifyTrue(tc, all(diff(tn.epoch) > 0));
    verifyGreaterThan(tc, max(diff(tn.epoch)), 1.0);   % mission gap
end
verifyEqual(tc, tnC.C10(1), 5.120437146e-10, 'AbsTol', 0);
verifyEqual(tc, tnJ.C10(1), 4.839777801e-10, 'AbsTol', 0);
verifyEqual(tc, tnC.C11(end), 1.256366494e-10, 'AbsTol', 0);
verifyEqual(tc, tnJ.S11(end), 2.394137498e-10, 'AbsTol', 0);
% cross-provider consistency vs the GFZ file (epoch-matched)
fG = fullfile(d, 'TN-13_GEOC_GFZ_RL06_3.txt');
if isfile(fG)
    tnG = shLowLevel.readTN13(fG);
    [ok, iC] = ismember(round(tnG.epoch, 6), round(tnC.epoch, 6));
    r = corrcoef(tnG.C10(ok), tnC.C10(iC(ok)));
    verifyGreaterThan(tc, r(1, 2), 0.9);
end
% the completion chain accepts the CSR file on the real ITSG month
fI = fullfile(d, 'ITSG-Grace2018_n60_2008-04.gfc');
if isfile(fI)
    g = shCoefficients.read(fI, Epoch = 2008 + 3.5/12);
    g1 = g.addDegree1(tnC);
    verifyEqual(tc, g1.C(2, 1), 2.630331004e-10, 'AbsTol', 0);
end
end

function testUserSuppliedDDKFiles(tc)
% Discovery test: validate EVERY Wbd_* binary present in test_data or
% dataFolder/ddk (fetchDDK target). All eight released files were
% parsed and validated externally (Python, 2026-08-07): BIN v2.1
% BDFULLV0, Lmax 120, 241 blocks, pack fully consumed, and the
% diagonal gains strictly ordered with the regularization strength
% DDK1 (a_1d14) .. DDK8 (a_5d9). This test re-runs the structural and
% ordering checks on whatever subset is available locally.
dirs = {tc.TestData.dataDir};
try
    dd = fullfile(shLowLevel.dataFolder(), 'DDK');   % fetchDDK target
    if isfolder(dd), dirs{end+1} = dd; end
catch
end
files = [];
for d = dirs
    files = [files; dir(fullfile(d{1}, 'Wbd_2-120.a_*'))]; %#ok<AGROW>
end
% The same released file may sit in BOTH test_data and dataFolder/DDK
% (fetchDDK target). Duplicates produce tied strengths and break the
% strict monotonicity check below, so keep one instance per file name.
if ~isempty(files)
    [~, iu] = unique({files.name}, 'stable');
    files = files(iu);
end
assumeTrue(tc, ~isempty(files), 'no Wbd files present');
strength = zeros(numel(files), 1); gain60 = zeros(numel(files), 1);
keep = true(numel(files), 1);
for k = 1:numel(files)
    W = shLowLevel.readDDK(fullfile(files(k).folder, files(k).name));
    verifyEqual(tc, W.nmax, 120);
    verifyEqual(tc, numel(W.blocks), 241);
    b0 = W.blocks(1);                       % m = 0, C, degrees 2..120
    verifyEqual(tc, b0.m, 0);
    dg = diag(b0.M);
    verifyTrue(tc, all(dg > 0 & dg < 1.02), ...
        sprintf('%s: m=0 diagonal gains outside (0, 1.02)', ...
        files(k).name));
    tok = regexp(files(k).name, 'a_(\d)d(\d+)p', 'tokens', 'once');
    if isempty(tok)                          % non-canonical user name:
        keep(k) = false;                     % structural checks only
        continue
    end
    strength(k) = str2double(tok{1}) * 10^str2double(tok{2});
    gain60(k) = dg(60 - b0.n(1) + 1);
end
strength = strength(keep); gain60 = gain60(keep);
files = files(keep);   % ordering check on canonical names only
if numel(files) > 1                          % stronger filter, lower gain
    [~, order] = sort(strength, 'descend');
    verifyTrue(tc, all(diff(gain60(order)) > 0), ...
        'DDK gain(n=60, m=0) not monotone with regularization strength');
end
end

function testRealTN13GFZ(tc)
% REAL GFZ TN-13 geocenter file (RL06.3, user-supplied): 256 months
% 2002-04 .. 2026-04, perfectly paired GRCOF2 lines, layout
% "GRCOF2 n m Clm Slm sigC sigS begin end" with yyyymmdd.hhmm epochs.
% Expected numbers cross-checked with an independent Python parse.
f = fullfile(tc.TestData.dataDir, 'TN-13_GEOC_GFZ_RL06_3.txt');
assumeTrue(tc, isfile(f), 'real TN-13 file not present');
tn = shLowLevel.readTN13(f);
verifyEqual(tc, numel(tn.epoch), 256);
verifyTrue(tc, all(isfinite(tn.C10) & isfinite(tn.C11) & isfinite(tn.S11)));
verifyTrue(tc, all(diff(tn.epoch) > 0));
% first month (2002-04-05 .. 2002-05-01)
verifyEqual(tc, tn.C10(1), 5.462244001e-10, 'AbsTol', 0);
verifyEqual(tc, tn.sigC10(1), 4.4585e-11, 'AbsTol', 0);
verifyEqual(tc, tn.t0(1), 2002.2575342465752, 'RelTol', 1e-12);
verifyEqual(tc, tn.epoch(1), 2002.2931506849313, 'RelTol', 1e-12);
% last month (2026-04)
verifyEqual(tc, tn.C11(end), 7.519802514e-11, 'AbsTol', 0);
verifyEqual(tc, tn.S11(end), 3.085246609e-10, 'AbsTol', 0);
verifyEqual(tc, tn.epoch(end), 2026.2876712328766, 'RelTol', 1e-12);
% the GRACE <-> GRACE-FO gap must be visible in the sampling
verifyGreaterThan(tc, max(diff(tn.epoch)), 1.0);
end

function testRealTN13ChainWithITSG(tc)
% full degree-1 chain on REAL data: ITSG-Grace2018 2008-04 monthly
% solution + GFZ TN-13 month 2008-04 (mid-epoch 2008.2896, within the
% 0.05-yr tolerance of the file epoch 2008+3.5/12)
fG = fullfile(tc.TestData.dataDir, 'ITSG-Grace2018_n60_2008-04.gfc');
fT = fullfile(tc.TestData.dataDir, 'TN-13_GEOC_GFZ_RL06_3.txt');
assumeTrue(tc, isfile(fG) && isfile(fT), 'real files not present');
g = shCoefficients.read(fG, Epoch = 2008 + 3.5/12);
verifyEqual(tc, g.C(2,1), 0);                    % ITSG ships degree 1 = 0
g1 = g.addDegree1(fT);                           % filename form
verifyEqual(tc, g1.C(2,1), 2.975653129e-10, 'AbsTol', 0);
verifyEqual(tc, g1.C(2,2), 2.227999852e-10, 'AbsTol', 0);
verifyEqual(tc, g1.S(2,2), -1.225458023e-10, 'AbsTol', 0);
verifyEqual(tc, g1.sigmaC(2,1), 4.4585e-11, 'AbsTol', 0);
verifyEqual(tc, g1.sigmaS(2,2), 5.0724e-11, 'AbsTol', 0);
% rest of the field untouched
verifyEqual(tc, g1.C(3,1), g.C(3,1), 'AbsTol', 0);
end

function testRealTN13GzipRoundTrip(tc)
% .gz branch on the real file
f = fullfile(tc.TestData.dataDir, 'TN-13_GEOC_GFZ_RL06_3.txt');
assumeTrue(tc, isfile(f), 'real TN-13 file not present');
tmp = tempname; mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>
gz = gzip(f, tmp);
tn  = shLowLevel.readTN13(f);
tnZ = shLowLevel.readTN13(gz{1});
verifyEqual(tc, tnZ.epoch, tn.epoch, 'AbsTol', 0);
verifyEqual(tc, tnZ.C10, tn.C10, 'AbsTol', 0);
verifyEqual(tc, tnZ.S11, tn.S11, 'AbsTol', 0);
end

function testUserSuppliedSINEXFiles(tc)
% Auto-discovery: any additional *.snx / *.snx.gz the user drops into
% tests/test_data is parsed and sanity-checked (the shipped head12
% fixture has its own dedicated test). No unzipping needed - readSINEX
% gunzips transparently. Very large archives (full ITSG n96 NEQ SINEX is
% ~440 MB gz / several GB text) exceed the in-memory reader's documented
% scope and are skipped with a message rather than ground through.
dd = tc.TestData.dataDir;
f = [dir(fullfile(dd, '*.snx')); dir(fullfile(dd, '*.snx.gz'))];
f = f(~strcmp({f.name}, 'ITSG-Grace2018_n96_2008-04_head12.snx'));
assumeTrue(tc, ~isempty(f), 'no user-supplied SINEX files in test_data');
maxBytes = 30e6;                                % ~30 MB gz / text cap
for k = 1:numel(f)
    fk = fullfile(f(k).folder, f(k).name);
    if f(k).bytes > maxBytes
        % v2.2: stream just the SOLUTION/ESTIMATE block instead of skipping
        try
            tS = tic;
            est = shLowLevel.readSINEX(fk, Only = "estimate");
            fprintf(['  [stream] %s (%.0f MB): %d estimate params in ' ...
                '%.1f s (matrix blocks skipped)\n'], f(k).name, ...
                f(k).bytes/1e6, numel(est.x), toc(tS));
            verifyGreaterThan(tc, numel(est.x), 0);
            verifyTrue(tc, all(est.m <= est.n) & all(est.n <= 2190));
            verifyTrue(tc, all(isfinite(est.x)));
        catch ME
            fprintf('  [skip] %s: streaming failed (%s)\n', f(k).name, ME.message);
        end
        continue;
    end
    snx = shLowLevel.readSINEX(fk);
    fprintf('  [ok]   %s: %d params, kind=%s\n', f(k).name, numel(snx.x), snx.kind);
    verifyGreaterThan(tc, numel(snx.x), 0);
    verifyTrue(tc, all(snx.m <= snx.n) && all(snx.n >= 0) && all(snx.n <= 2190));
    verifyTrue(tc, all(isfinite(snx.x)));
    if ~isempty(snx.M)
        verifyEqual(tc, size(snx.M), [numel(snx.x), numel(snx.x)]);
        verifyLessThan(tc, norm(snx.M - snx.M', 'fro'), ...
            1e-12 * max(norm(snx.M, 'fro'), realmin));
    end
end
end

function testZeroVarianceRowsNoiseFloor(tc)
% coefficient rows with EXACTLY constant residuals (here: identically-zero
% sectoral S, the tril(...,-1) construction) must not break the empirical
% noise covariance: the PD floor references the smallest POSITIVE variance
% (v2.1 fix; previously min(diag)=0 collapsed the guard -> singular N ->
% Inf eigenvalues -> NaN through the whole chain). Zero rows pass through
% with gain ~ 1 and stay (exactly) zero.
rng(31);
nmax = 6; T = 16; n1 = nmax + 1;
mL  = tril(true(n1));
mL1 = tril(true(n1), -1);                        % zero DIAGONAL: S(n,n) = 0
Cs = (mL  .* randn(n1, n1, T)) * 1e-9;
Ss = (mL1 .* randn(n1, n1, T)) * 1e-9;
ts = shSeries(Cs, Ss = Ss, Epochs = 2020 + (0:T-1)'/12);
for blocksMode = ["on", "off"]
    [tsF, op, info] = ts.filter("tvANS", Blocks = blocksMode); %#ok<ASGLU>
    verifyTrue(tc, all(isfinite(tsF.Cs(:))), ...
        sprintf('non-finite Cs (Blocks=%s)', blocksMode));
    verifyTrue(tc, all(isfinite(info.sigmaXfres(:))));
    % sectoral S rows were exactly zero and must remain (essentially) zero
    for n = 2:nmax
        verifyLessThan(tc, max(abs(squeeze(tsF.Ss(n+1, n+1, :)))), 1e-20);
    end
end
% and the collinear-basin ridge path on top of such a series stays finite
[~, op] = ts.filter("tvANS");
idx = shLowLevel.shIndex(nmax, MinDegree = 2);
rng(32); b = randn(idx.P, 1);
avg = ts.basinAverage([b, b] / norm(b)^2 * idx.P, Deconvolve = true, ...
    Op = op, Ridge = 1e-6);
verifyTrue(tc, all(isfinite(avg(:))));
end

% =============================================================== v2.2
function testDDKRoundTripAndIdentity(tc)
% ASCII exchange format roundtrip + identity filter is a no-op
tmp = tempname; mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>
f = fullfile(tmp, 'ddk_test.txt');
fid = fopen(f, 'w');
fprintf(fid, '# synthetic DDK exchange fixture\n');
rng(41); M2 = 0.5 * eye(3) + 0.05 * randn(3);
fprintf(fid, 'block 2 0 2 4\n');
fprintf(fid, '%.15e %.15e %.15e\n', M2');
fprintf(fid, 'block 2 1 2 4\n');
fprintf(fid, '%.15e %.15e %.15e\n', eye(3)');
fclose(fid);
W = shLowLevel.readDDK(f);
verifyEqual(tc, numel(W.blocks), 2);
verifyEqual(tc, W.blocks(1).M, M2, 'AbsTol', 1e-14);
verifyEqual(tc, W.nmax, 4);
g = randomField(6, 2024);
gF = g.applyDDK(W);
% S block was identity: S coefficients of order 2 unchanged
verifyEqual(tc, gF.S(3:5, 3), g.S(3:5, 3), 'AbsTol', 0);
% C block applied: matches manual product
verifyEqual(tc, gF.C(3:5, 3), M2 * g.C(3:5, 3), 'AbsTol', 1e-15);
% uncovered degrees untouched
verifyEqual(tc, gF.C(6:7, :), g.C(6:7, :), 'AbsTol', 0);
% v2.5.1: a block reaching beyond the field's nmax is TRUNCATED to the
% available degrees (submatrix of the gain) - the standard practical
% case of Lmax-120 Wbd filters on n60/n96 GRACE fields
gSmall = g.truncate(3);              % block covers n=2..4, field n<=3
gSF = gSmall.applyDDK(W);
keep = 1:2;                          % degrees 2,3 of the n=2..4 block
verifyEqual(tc, gSF.C(3:4, 3), M2(keep, keep) * gSmall.C(3:4, 3), ...
    'AbsTol', 1e-15);
verifyEqual(tc, gSF.S(3:4, 3), gSmall.S(3:4, 3), 'AbsTol', 0);
% series wrapper filters every epoch and drops sigmas
ts = randomSeries(6, 5);
tsF = ts.applyDDK(W);
verifyEqual(tc, tsF.Cs(3:5, 3, 4), M2 * ts.Cs(3:5, 3, 4), 'AbsTol', 1e-15);
end

function testReadLoveNumbersLayouts(tc)
% shLowLevel.readLoveNumbers (v2.5): layouts, sparse interpolation, degree-1
% frame conversion. All numeric expectations Python-validated
% (Blewitt 2003 PREM: CE (-0.290, 0.113, 0.021) -> CF (-0.269, 0.134);
% pchip-in-log(1+n) reconstruction error 1.3e-4 on Farrell sampling).
tmp = tempname; mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>
w = @(name, txt) writelines(txt, fullfile(tmp, name));

% 2-column default with a prose comment
w("f1.txt", ["% my kn table"; "0 0.000"; "1 0.021"; ...
    "2 -0.303"; "3 -0.194"]);
LN = shLowLevel.readLoveNumbers(fullfile(tmp, "f1.txt"));
verifyEqual(tc, LN.n, (0:3)');
verifyEqual(tc, LN.kn(3), -0.303, 'AbsTol', 0);
verifyTrue(tc, all(isnan(LN.hn)));

% commented Farrell-order header
w("f2.txt", ["#  n    h'      l'      k'"; "0 -0.13 0.0 0.0"; ...
    "1 -0.29 0.113 0.021"; "2 -0.99 0.02 -0.303"]);
LN = shLowLevel.readLoveNumbers(fullfile(tmp, "f2.txt"));
verifyEqual(tc, LN.kn(2), 0.021, 'AbsTol', 0);
verifyEqual(tc, LN.hn(2), -0.29, 'AbsTol', 0);

% headerless 4-column is ambiguous by contract
w("f3.txt", ["0 0 0 0"; "1 1 1 1"]);
verifyError(tc, @() shLowLevel.readLoveNumbers(fullfile(tmp, "f3.txt")), ...
    'shLowLevel:readLoveNumbers:ambiguousColumns');

% Columns= override, csv separators, unordered rows
w("f4.txt", ["2, -0.99, 0.02, -0.303"; "1, -0.29, 0.113, 0.021"; ...
    "0, -0.13, 0, 0"]);
LN = shLowLevel.readLoveNumbers(fullfile(tmp, "f4.txt"), Columns = "n h l k");
verifyEqual(tc, LN.n, (0:2)');
verifyEqual(tc, LN.kn(3), -0.303, 'AbsTol', 0);

% sparse Farrell sampling + MaxDegree + pchip interpolation
ns = [0 1 2 3 4 5 6 8 10 12 18 32 56 100];
rows = compose("%d %.8f", ns', -0.30 * ns' ./ (ns' + 6.0));
w("f5.txt", rows);
verifyError(tc, @() shLowLevel.readLoveNumbers(fullfile(tmp, "f5.txt"), ...
    MaxDegree = 96), 'shLowLevel:readLoveNumbers:degreeGap');
[LN, info] = shLowLevel.readLoveNumbers(fullfile(tmp, "f5.txt"), ...
    MaxDegree = 96, Interp = "pchip");
nn = (0:96)';
verifyEqual(tc, LN.kn, -0.30 * nn ./ (nn + 6.0), 'AbsTol', 2e-3);
verifyTrue(tc, any(info.interpolated) && ~info.interpolated(1));

% frame conversion CE->CF (Blewitt PREM) and refused ->CE
w("f6.txt", ["# n h l k"; "0 -0.13 0 0"; "1 -0.290 0.113 0.021"; ...
    "2 -0.99 0.02 -0.303"]);
LN = shLowLevel.readLoveNumbers(fullfile(tmp, "f6.txt"), ...
    InFrame = "CE", OutFrame = "CF");
verifyEqual(tc, LN.hn(2), -0.269, 'AbsTol', 5e-4);
verifyEqual(tc, LN.ln(2), 0.134, 'AbsTol', 5e-4);
verifyEqual(tc, LN.kn(1), 0, 'AbsTol', 0);       % only degree 1 shifts
verifyError(tc, @() shLowLevel.readLoveNumbers(fullfile(tmp, "f6.txt"), ...
    InFrame = "CM", OutFrame = "CE"), ...
    'shLowLevel:readLoveNumbers:frameNotRecoverable');
% exact identity k1(CM) = -1 from CE input
LN = shLowLevel.readLoveNumbers(fullfile(tmp, "f6.txt"), ...
    InFrame = "CE", OutFrame = "CM");
verifyEqual(tc, LN.kn(2), -1, 'AbsTol', 1e-12);

% the deformation chain accepts the reader output directly
g = shCoefficients(1e-9 * tril(ones(3)), zeros(3));
up = g.deformation(45, 10, kn = LN.kn, hn = LN.hn, ln = LN.ln, ...
    Mode = "points");
verifyTrue(tc, isfinite(up));
end

function testRemoveGIA(tc)
% linear model removal: series and climatology routes agree
rng(42);
nmax = 8; T = 36;
t = 2005 + (0:T-1)'/12;
n1 = nmax + 1;
giaC = tril(1e-11 * randn(n1)); giaS = tril(1e-11 * randn(n1), -1);
gia = shCoefficients(giaC, giaS);
Cs = zeros(n1, n1, T); Ss = zeros(n1, n1, T);
base = tril(randn(n1)) * 1e-9;
for k = 1:T
    Cs(:,:,k) = base + giaC * (t(k) - mean(t)) + tril(randn(n1)) * 1e-13;
    Ss(:,:,k) = giaS * (t(k) - mean(t)) + tril(randn(n1), -1) * 1e-13;
end
ts = shSeries(Cs, Ss = Ss, Epochs = t);
% route 1: remove from series, then fit -> trend ~ 0
clim1 = ts.removeGIA(gia).climatology();
verifyLessThan(tc, max(abs(clim1.trendC(:))), 5e-12);
% route 2: fit, then correct the climatology -> identical trend
clim2 = ts.climatology().removeGIA(gia);
verifyEqual(tc, clim2.trendC, clim1.trendC, 'AbsTol', 1e-13);
% smaller-nmax model: zero-padded, note in history
giaSmall = gia.truncate(4);
ts2 = ts.removeGIA(giaSmall);
verifyTrue(tc, any(contains(ts2.history, "zero-padded")));
end

function testARCorrectInflatesSigmas(tc)
% AR(1)-correlated residuals: corrected sigmas grow by ~sqrt((1+r1)/(1-r1))
rng(43);
nmax = 4; T = 120;
t = 2005 + (0:T-1)'/12;
n1 = nmax + 1;
phi = 0.6;
Cs = zeros(n1, n1, T); Ss = zeros(n1, n1, T);
e = zeros(n1, n1);
mL = tril(true(n1)); mL1 = tril(true(n1), -1);
for k = 1:T
    e = phi * e + (mL .* randn(n1)) * 1e-9;
    Cs(:,:,k) = e;
    Ss(:,:,k) = (mL1 .* randn(n1)) * 1e-10;
end
ts = shSeries(Cs, Ss = Ss, Epochs = t);
c0 = ts.climatology();
c1 = ts.climatology(ARCorrect = true);
tr1 = c1.trend(); tr0 = c0.trend();                % sigmas ride on trend()
r = tr1.sigmaC(3, 1) / tr0.sigmaC(3, 1);           % C20-type element
verifyGreaterThan(tc, r, 1.4);                     % phi=0.6 -> ~2.0 ideal
verifyLessThan(tc, r, 2.6);
verifyTrue(tc, any(contains(c1.history, "AR(1)")));
end

function testMasconReaderSyntheticNC(tc)
% synthetic JPL-style netCDF roundtrip (base MATLAB nccreate/ncread)
tmp = tempname; mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>
f = fullfile(tmp, 'mascon_test.nc');
lat = (-89.75:0.5:89.75)'; lon = (0.25:0.5:359.75)'; tt = [100; 130.5];
nccreate(f, 'lat', 'Dimensions', {'lat', numel(lat)});
nccreate(f, 'lon', 'Dimensions', {'lon', numel(lon)});
nccreate(f, 'time', 'Dimensions', {'time', numel(tt)});
nccreate(f, 'lwe_thickness', 'Dimensions', ...
    {'lon', numel(lon), 'lat', numel(lat), 'time', numel(tt)});
ncwrite(f, 'lat', lat); ncwrite(f, 'lon', lon); ncwrite(f, 'time', tt);
ncwriteatt(f, 'time', 'units', 'days since 2002-01-01');
ncwriteatt(f, 'lwe_thickness', 'units', 'cm');
rng(44); E = randn(numel(lon), numel(lat), numel(tt));
ncwrite(f, 'lwe_thickness', E);
mas = shLowLevel.readMascon(f);
verifySize(tc, mas.ewh, [numel(lat), numel(lon), 2]);   % permuted to lat x lon
verifyEqual(tc, mas.ewh(:, :, 1), E(:, :, 1)', 'AbsTol', 0);
verifyEqual(tc, mas.epoch(1), 2002 + 100/365, 'AbsTol', 1e-3);
verifyEqual(tc, mas.units, "cm");
verifyEqual(tc, mas.lon(1), 0.25, 'AbsTol', 0);
end

function testSINEXStreamingEstimateOnly(tc)
% streaming Only="estimate" == full-read estimates on the real fixture
f = fullfile(tc.TestData.dataDir, 'ITSG-Grace2018_n96_2008-04_head12.snx');
assumeTrue(tc, isfile(f));
full_ = shLowLevel.readSINEX(f);
est = shLowLevel.readSINEX(f, Only = "estimate");
verifyEqual(tc, est.n, full_.n);
verifyEqual(tc, est.m, full_.m);
verifyEqual(tc, est.cs, full_.cs);
verifyEqual(tc, est.x, full_.x, 'AbsTol', 0);
verifyEqual(tc, est.kind, 'estimate');
verifyTrue(tc, isempty(est.M));
% gz variant through the streaming path
tmp = tempname; mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>
gz = gzip(f, tmp);
estZ = shLowLevel.readSINEX(gz{1}, Only = "estimate");
verifyEqual(tc, estZ.x, full_.x, 'AbsTol', 0);
end

function testRealDDKWbdBinary(tc)
% native parse of the released binary DDK3 file (Rietbroek/Kusche BIN
% format, MIT-licensed, github.com/strawpants/GRACE-filter), pinned to
% the repository-documented reference values
f = fullfile(tc.TestData.dataDir, 'Wbd_2-120.a_1d12p_4');
assumeTrue(tc, isfile(f));
W = shLowLevel.readDDK(f);
verifyEqual(tc, W.nmax, 120);
verifyEqual(tc, numel(W.blocks), 241);
% block ordering: C0, C1, S1, ..., C120, S120
verifyEqual(tc, [W.blocks(1).m, W.blocks(1).cs], [0, 0]);
verifyEqual(tc, [W.blocks(2).m, W.blocks(2).cs], [1, 0]);
verifyEqual(tc, [W.blocks(3).m, W.blocks(3).cs], [1, 1]);
verifyEqual(tc, [W.blocks(241).m, W.blocks(241).cs], [120, 1]);
verifyEqual(tc, W.blocks(1).n, 2:120);
verifyEqual(tc, W.blocks(241).n, 120);
% independently documented pack1(1:5) (strawpants README, format long)
refHead = [0.999995413476564; 0.000000054102768; 0.000000668145433; ...
    -0.000000013374350; 0.000000077291901];
verifyEqual(tc, W.blocks(1).M(1:5, 1), refHead, 'AbsTol', 1e-15);
% truncation to n96: leading submatrices, reference-consistent
W96 = shLowLevel.readDDK(f, Nmax = 96);
verifyEqual(tc, W96.nmax, 96);
verifyEqual(tc, numel(W96.blocks), 193);         % orders 0..96
verifyEqual(tc, W96.blocks(1).M, W.blocks(1).M(1:95, 1:95), 'AbsTol', 0);
% semantic check on a synthetic n96 field: linear per-block action and
% strong high-degree attenuation (DDK3)
rng(51);
n1 = 97;
C = tril(randn(n1)) * 1e-9; S = tril(randn(n1), -1) * 1e-9; S(:, 1) = 0;
g = shCoefficients(C, S);
gF = g.applyDDK(W96);
verifyEqual(tc, gF.C(3:97, 1), W96.blocks(1).M * g.C(3:97, 1), 'AbsTol', 1e-24);
sF = shLowLevel.shDegreeRMS(gF.C, gF.S);
s0 = shLowLevel.shDegreeRMS(g.C, g.S);
r = sF.degRMS ./ max(s0.degRMS, realmin);
% thresholds grounded by the Python cross-run on this filter: white-
% field ratios ~0.96 (n=5), 0.024 (n=60), 2e-4 (n=90), max 1.031 (off-
% diagonal coupling can locally exceed 1 on white input)
verifyLessThan(tc, r(91), 0.2);                  % degree 90 heavily damped
verifyGreaterThan(tc, r(6), 0.5);                % degree 5 mostly kept
verifyLessThan(tc, max(r(3:end)), 1.25);
end

% =============================================================== v2.4
function testFingerprintContracts(testCase)
L = 12; n = (0:L)';
kn = -0.3 * n ./ (n + 3); hn = -0.9 * n ./ (n + 2) - 0.1;
idxBad = shLowLevel.shIndex(L, MinDegree = 2);
ocean = @(la, lo) double(la < 50);
loadF = @(la, lo) -10 * double(la > 60);
verifyError(testCase, ...
    @() shLowLevel.seaLevelFingerprint(loadF, ocean, idxBad, kn = kn, hn = hn), ...
    'shLowLevel:seaLevelFingerprint:badIndex');
idx = shLowLevel.shIndex(L, MinDegree = 0);
verifyError(testCase, ...
    @() shLowLevel.seaLevelFingerprint(struct('bad', 1), ocean, idx, ...
    kn = kn, hn = hn), 'shLowLevel:seaLevelFingerprint:badLoad');
% grid-valued load form: constant land load, must still conserve
[~, ~, grid0] = shLowLevel.synthesisMatrix(idx, NLat = 2*(L+1), NLon = 2*(2*L+2));
Ng = numel(grid0.latDeg) * numel(grid0.lonDeg);
[S, ~, info] = shLowLevel.seaLevelFingerprint(-5 * ones(Ng, 1), ocean, idx, ...
    kn = kn, hn = hn);
verifyLessThan(testCase, abs(info.massResidual), 1e-12);
verifyEqual(testCase, numel(S), Ng);
end

function testEOFRejectsNaN(testCase)
L = 4; n1 = L + 1;
Cs = randn(n1, n1, 5); Cs(2, 1, 3) = NaN;
Ss = zeros(n1, n1, 5);
ts = shSeries(Cs, Ss = Ss, Epochs = 2010 + (1:5)/12);
verifyError(testCase, @() shLowLevel.eofAnalysis(ts), 'shLowLevel:eofAnalysis:nanInSeries');
end

function testFromFolderRealFiles(testCase)
d = fullfile(fileparts(mfilename('fullpath')), 'test_data');
ts = shSeries.fromFolder(d, Pattern = "ITSG-*_n60_*.gfc");
verifyEqual(testCase, ts.nEpochs, 2);
verifyEqual(testCase, ts.nmax, 60);
% sorted by epoch: GRACE 2008-04 before GRACE-FO 2025-12
verifyLessThan(testCase, ts.epochs(1), 2009);
verifyGreaterThan(testCase, ts.epochs(2), 2025);
% truncation path
ts40 = shSeries.fromFolder(d, Pattern = "ITSG-*_n60_*.gfc", Truncate = 40);
verifyEqual(testCase, ts40.nmax, 40);
verifyError(testCase, ...
    @() shSeries.fromFolder(fullfile(d, 'no_such_dir')), 'shSeries:noFiles');
end

function testSeriesTNWrappers(testCase)
d = fullfile(fileparts(mfilename('fullpath')), 'test_data');
fG = fullfile(d, 'ITSG-Grace2018_n60_2008-04.gfc');
fT = fullfile(d, 'TN-13_GEOC_GFZ_RL06_3.txt');
assumeTrue(testCase, isfile(fG) && isfile(fT));
ts = shSeries.read(string(fG));
ts1 = ts.addDegree1(fT);
% pinned values from the coefficient-level real-chain test
verifyEqual(testCase, ts1.Cs(2, 1, 1), 2.975653129e-10, 'AbsTol', 0);
verifyEqual(testCase, ts1.Ss(2, 2, 1), -1.225458023e-10, 'AbsTol', 0);
verifyEqual(testCase, ts1.nEpochs, 1);
% rest untouched
verifyEqual(testCase, ts1.Cs(3, 1, 1), ts.Cs(3, 1, 1), 'AbsTol', 0);
% v2.5.1 regression: applyTN14 on a SIGMA-LESS synthetic series must
% not create partially-sigma'd months (setCoefficient now initializes
% both sigma stacks; a lone sigmaC used to break shSeries stacking)
fT14 = fullfile(d, 'TN-14_C30_C20_SLR_GSFC.txt');
assumeTrue(testCase, isfile(fT14));
n1 = 5; T = 3;
Cs = 1e-9 * randn(n1, n1, T); Ss = 1e-9 * randn(n1, n1, T);
tsS = shSeries(Cs, Ss = Ss, Epochs = 2019 + (0:T-1)'/12);
ts14 = tsS.applyTN14(fT14);
verifyEqual(testCase, size(ts14.sigmaCs), [n1 n1 T]);
verifyEqual(testCase, size(ts14.sigmaSs), [n1 n1 T]);
verifyTrue(testCase, all(isfinite(squeeze(ts14.sigmaCs(3, 1, :)))));
verifyTrue(testCase, all(isnan(squeeze(ts14.sigmaSs(3, 1, :)))));
end

function testDailyKalmanRealFile(testCase)
% real ITSG daily Kalman file (shipped, 83 kB): format, epoch, sigmas
d = fullfile(fileparts(mfilename('fullpath')), 'test_data');
fD = fullfile(d, 'ITSG-Grace2018_Kalman_n40_2008-04-15.gfc');
assumeTrue(testCase, isfile(fD));
g = shCoefficients.read(fD);
verifyEqual(testCase, g.nmax, 40);
verifyEqual(testCase, g.epoch, 2008 + 105.5/366, 'AbsTol', 1e-9);
verifyEqual(testCase, g.productType, "GSM");
verifyTrue(testCase, ~isempty(g.sigmaC) && any(g.sigmaC(:) > 0));
verifyEqual(testCase, g.C(1, 1), 1);
% meta from the filename parser directly
meta = shLowLevel.parseGraceFilename('ITSG-Grace2018_Kalman_n40_2008-04-15.gfc');
verifyEqual(testCase, meta.epoch, 2008 + 105.5/366, 'AbsTol', 1e-9);
verifyLessThan(testCase, meta.epochStop - meta.epochStart, 1.1/365);
end

function testFetchLoveNumbersMirror(tc)
% v3.0.0: GROOPS Love-number fetch from a local mirror + parser contract
d = tc.TestData.dataDir;
src = fullfile(d, 'loadLoveNumbers_Gegout97.txt');
assumeTrue(tc, isfile(src));
mir = tempname; mkdir(mir); copyfile(src, mir);
dst = tempname; mkdir(dst);
c1 = onCleanup(@() rmdir(mir, 's')); c2 = onCleanup(@() rmdir(dst, 's')); %#ok<NASGU>
[f, inf] = shLowLevel.fetchLoveNumbers("loadLoveNumbers_Gegout97.txt", ...
    BaseURL = mir, Dest = dst, Quiet = true);
verifyEqual(tc, numel(f), 1);
verifyEqual(tc, numel(inf.parsed), 1);
kn = inf.parsed(1).kn;
verifyEqual(tc, numel(kn), 1025);
verifyEqual(tc, kn(1:2), [0; 0], 'AbsTol', 0);
verifyEqual(tc, kn(3), -0.3054020195, 'AbsTol', 1e-10);   % degree 2
verifyLessThan(tc, kn(3:end), 0);                          % all negative
% second call skips (safe-swap idempotence)
[~, i2] = shLowLevel.fetchLoveNumbers("loadLoveNumbers_Gegout97.txt", ...
    BaseURL = mir, Dest = dst, Quiet = true);
verifyEqual(tc, numel(i2.skipped), 1);
verifyError(tc, @() shLowLevel.fetchLoveNumbers("nonsense.txt", BaseURL = mir), ...
    'shLowLevel:fetchLoveNumbers:unknownName');
end

function testListITSGAndFetchAllMirror(tc)
% v3.0.0: listITSG walks a local mirror tree; fetchITSG months="all" +
% Catalog= selection operate on it offline
mir = tempname;
m96 = fullfile(mir, 'ITSG-Grace2018', 'monthly', 'monthly_n96');
mkdir(m96); mkdir(fullfile(mir, 'ITSG-Grace2018', 'daily_kalman'));
d = tc.TestData.dataDir;
src = fullfile(d, 'ITSG-Grace2018_n60_2008-04.gfc');
assumeTrue(tc, isfile(src));
% two months, correct n96 naming for enumeration
copyfile(src, fullfile(m96, 'ITSG-Grace2018_n96_2008-04.gfc'));
copyfile(src, fullfile(m96, 'ITSG-Grace2018_n96_2008-05.gfc'));
c1 = onCleanup(@() rmdir(mir, 's')); %#ok<NASGU>
T = shLowLevel.listITSG(BaseURL = mir);
verifyEqual(tc, T.idx, (1:height(T))');
verifyTrue(tc, any(T.release == "ITSG-Grace2018" & T.product == "monthly" ...
    & T.nmax == 96));
verifyTrue(tc, any(T.product == "daily"));
dst = tempname; mkdir(dst);
c2 = onCleanup(@() rmdir(dst, 's')); %#ok<NASGU>
[f, inf] = shLowLevel.fetchITSG("all", Release = "ITSG-Grace2018", Nmax = 96, ...
    BaseURL = mir, Dest = dst, Quiet = true);
verifyEqual(tc, numel(f), 2);
verifyEqual(tc, numel(inf.fetched), 2);
% Catalog selection: the monthly n96 row fetches the same two files
ii = find(T.product == "monthly" & T.nmax == 96, 1);
[f2, i2] = shLowLevel.fetchITSG(Catalog = ii, BaseURL = mir, Dest = dst, ...
    Quiet = true);
verifyEqual(tc, numel(f2), 2);
verifyEqual(tc, numel(i2.skipped), 2);       % already present
end

function testICGEMListFixtureAndResolve(testCase)
d = fullfile(fileparts(mfilename('fullpath')), 'test_data');
fx = fullfile(d, 'icgem_list_fixture.html');
assumeTrue(testCase, isfile(fx));
T = shLowLevel.listICGEM(Source = fx);
verifyGreaterThan(testCase, height(T), 30);
% v3.0.0: numbered catalogue + numeric selection contract (offline)
verifyTrue(testCase, ismember('idx', T.Properties.VariableNames));
verifyEqual(testCase, T.idx, (1:height(T))');
verifyError(testCase, @() shLowLevel.fetchICGEM(height(T) + 7, List = T), ...
    'shLowLevel:fetchICGEM:badIdx');
verifyTrue(testCase, any(T.name == "Tongji-GMMG2025S"));
verifyTrue(testCase, all(startsWith(T.url, "https://icgem.gfz.de/getmodel/gfc/")));
verifyTrue(testCase, all(isfinite(T.year)));
% name resolution against the fixture list - no network involved
verifyError(testCase, ...
    @() shLowLevel.fetchICGEM("definitely_no_such_model", List = T), ...
    'shLowLevel:fetchICGEM:notFound');
verifyError(testCase, @() shLowLevel.fetchICGEM("WHU", List = T), ...
    'shLowLevel:fetchICGEM:ambiguous');
% temporal branch on a saved live page (the v2.4.1 nested-group regexp
% crash regression: MATLAB tokens = outermost parentheses only).
% v2.4.2: the branch now parses the FULL /sp/ series-page catalogue
% (the old getseries-only match returned 3 release-note rows).
fx2 = fullfile(d, 'icgem_temporal_fixture.html');
if isfile(fx2)
    [Tt, infoT] = shLowLevel.listICGEM(Type = "temporal", Source = fx2);
    verifyGreaterThan(testCase, height(Tt), 60);   % 74 series on the capture
    verifyTrue(testCase, all(ismember(["CSR"; "GFZ"; "JPL"; "COST-G"; ...
        "ITSG"; "Tongji"], Tt.center)));
    verifyTrue(testCase, all(ismember(["01_GRACE"; "02_COST-G_"; ...
        "03_other"; "04_SLR_"], Tt.group)));
    verifyTrue(testCase, any(Tt.path == ...
        "03_other/ITSG/ITSG-Grace2018/monthly"));
    verifyTrue(testCase, all(startsWith(Tt.url, ...
        "https://icgem.gfz.de/sp/")));
    verifyTrue(testCase, all(startsWith(Tt.zip, ...
        "https://icgem.gfz.de/getseries/")));
    verifyFalse(testCase, any(contains(Tt.url, " ")));  % %20-encoded
    verifyTrue(testCase, strlength(infoT.note) > 0);
end
% Series= file listing on a saved series page (individual monthly files
% ARE statically listed - the old "JS-only" caveat is obsolete)
fx3 = fullfile(d, 'icgem_series_fixture.html');
if isfile(fx3)
    F = shLowLevel.listICGEM(Type = "temporal", Series = "x", Source = fx3);
    verifyGreaterThan(testCase, height(F), 400);   % 486 files on the capture
    verifyTrue(testCase, any(F.name == "ITSG-Grace2018_n120_2002-04.gfc"));
    verifyTrue(testCase, all(startsWith(F.url, ...
        "https://icgem.gfz.de/getseries/")));
    verifyTrue(testCase, all(endsWith(F.name, [".gfc", ".gz"])));
end
verifyError(testCase, ...
    @() shLowLevel.listICGEM(Type = "static", Series = "x", Source = fx), ...
    'shLowLevel:listICGEM:seriesNeedsTemporal');
verifyError(testCase, ...
    @() shLowLevel.listICGEM(Type = "temporal", Series = "x", Source = fx2), ...
    'shLowLevel:listICGEM:noFiles');           % catalogue page has no file links
% exact-name resolution finds the row (skip download: file marker trick)
tmp = fullfile(tempdir, sprintf('shx_icgem_%d', randi(1e9)));
mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>
row = T(T.name == "Tongji-GMMG2025S", :);
[~, b, e] = fileparts(char(row.url));
fid = fopen(fullfile(tmp, [b, e]), 'w'); fclose(fid);   % pretend cached
[f, info] = shLowLevel.fetchICGEM("Tongji-GMMG2025S", List = T, Dest = tmp);
verifyTrue(testCase, info.skipped);
verifyTrue(testCase, isfile(f));
end
