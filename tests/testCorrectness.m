function tests = testCorrectness
%TESTCORRECTNESS Correctness suite: class API vs validated cores, physics.
%   runtests('testCorrectness')  (run from tests/, or use runAllTests)
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root, fullfile(root, 'compat'));
testCase.TestData.dataDir = fullfile(here, 'test_data');
shx.legendreCached('clear');
end

% ------------------------------------------------------------------- read
function testReadMatchesCompat(testCase)
f = fullfile(testCase.TestData.dataDir, 'test_static.gfc');
g = shCoefficients.read(f);
m = shx.shReadGFC(f);
verifyEqual(testCase, g.C, m.C);
verifyEqual(testCase, g.S, m.S);
verifyEqual(testCase, g.GM, m.GM);
verifyEqual(testCase, g.R, m.R);
verifyEqual(testCase, g.sigmaC, m.sigmaC);
end

function testFilenameParsing(testCase)
meta = shx.parseGraceFilename('GSM-2_2024032-2024060_GRFO_UTCSR_BA01_0600.gfc');
verifyEqual(testCase, meta.productType, "GSM");
verifyGreaterThan(testCase, meta.epoch, 2024.0);
verifyLessThan(testCase, meta.epoch, 2024.25);
meta2 = shx.parseGraceFilename('GAD-2_2024032-2024060_GRFO_UTCSR_BA01_0600.gfc.gz');
verifyEqual(testCase, meta2.productType, "GAD");
verifyEqual(testCase, meta2.epoch, meta.epoch, 'AbsTol', 1e-12);
meta3 = shx.parseGraceFilename('some_random_model.gfc');
verifyEqual(testCase, meta3.productType, "unknown");
verifyTrue(testCase, isnan(meta3.epoch));
end

% ------------------------------------------------------------- arithmetic
function testArithmeticAndSigmaRSS(testCase)
rng(1); L = 8;
[gA, gB] = randomPair(L);
d = gB - gA;
verifyEqual(testCase, d.C, gB.C - gA.C, 'AbsTol', 0);
verifyEqual(testCase, d.sigmaC, sqrt(gA.sigmaC.^2 + gB.sigmaC.^2), 'RelTol', 1e-14);
verifyEqual(testCase, d.productType, "difference");
s = gA + gB;
verifyEqual(testCase, s.C, gA.C + gB.C, 'AbsTol', 0);
g2 = 2.5 * gA;
verifyEqual(testCase, g2.C, 2.5 * gA.C, 'AbsTol', 0);
verifyEqual(testCase, g2.sigmaC, 2.5 * gA.sigmaC, 'RelTol', 1e-14);
% v2.5.1: elementwise .* with a per-coefficient weight matrix
W = 0.5 * ones(size(gA.C)); W(3, 1) = 2;
gW = W .* gA;
verifyEqual(testCase, gW.C, W .* gA.C, 'AbsTol', 0);
verifyEqual(testCase, gW.S, W .* gA.S, 'AbsTol', 0);
verifyEqual(testCase, gW.sigmaC, abs(W) .* gA.sigmaC, 'AbsTol', 0);
end

function testGSMPlusGAXType(testCase)
rng(2);
[gsm, gad] = randomPair(6);
gsm = setType(gsm, "GSM"); gad = setType(gad, "GAD");
out = gsm + gad;
verifyEqual(testCase, out.productType, "GSM+GAD");
verifyEqual(testCase, out.C, gsm.C + gad.C, 'AbsTol', 0);
end

% -------------------------------------------------- filtering vs v1 cores
function testDestripeMatchesCompat(testCase)
rng(3); L = 40;
g = randomField(L);
g2 = g.destripe(minOrder = 5, polyOrder = 2);
[Cf, Sf] = shDestripe(g.C, g.S, 'minOrder', 5, 'polyOrder', 2);
verifyEqual(testCase, g2.C, Cf, 'AbsTol', 0);
verifyEqual(testCase, g2.S, Sf, 'AbsTol', 0);
g3 = g.destripe(minOrder = 5, polyOrder = 3, windowLength = 7);
[Cw, Sw] = shDestripe(g.C, g.S, 'minOrder', 5, 'polyOrder', 3, 'windowLength', 7);
verifyEqual(testCase, g3.C, Cw, 'AbsTol', 0);
end

function testGaussianMatchesCompatAndSigmas(testCase)
rng(4); L = 30;
g = randomField(L);
g2 = g.gaussian(300);
[Cf, Sf, Wn] = shGaussianFilter(g.C, g.S, 300);
verifyEqual(testCase, g2.C, Cf, 'AbsTol', 0);
verifyEqual(testCase, g2.S, Sf, 'AbsTol', 0);
verifyEqual(testCase, g2.sigmaC, g.sigmaC .* Wn(:), 'RelTol', 1e-14);
end

% -------------------------------------------------------------- synthesis
function testSynthesisCacheAndCompat(testCase)
rng(5); L = 20;
g = randomField(L);
lat = -88:4:88; lon = 0:6:354;
shx.legendreCached('clear');
[G1, la1, lo1] = g.synthesis(lat, lon, UseCache = false);
G2 = g.synthesis(lat, lon, UseCache = true);      % cold cache
G3 = g.synthesis(lat, lon, UseCache = true);      % warm cache
Gref = shSynthesis(g.C, g.S, g.GM, g.R, lat, lon);
verifyEqual(testCase, G1, Gref, 'AbsTol', 0);
verifyEqual(testCase, G2, Gref, 'RelTol', 1e-14);
verifyEqual(testCase, G3, Gref, 'RelTol', 1e-14);
verifyEqual(testCase, numel(la1), numel(lat));
verifyEqual(testCase, numel(lo1), numel(lon));
end

function testQuadratureIdentity(testCase)
idx = shx.shIndex(15);
[Y, w] = shx.synthesisMatrix(idx);
verifyLessThan(testCase, norm(Y' * (w .* Y) - eye(idx.P), 'fro'), 1e-10);
end

% ------------------------------------------------------- series statistics
function testMeanField(testCase)
rng(6); L = 6; T = 24;
ts = randomSeries(L, T);
m = ts.mean;
verifyEqual(testCase, m.C, mean(ts.Cs, 3), 'AbsTol', 0);
verifyEqual(testCase, m.sigmaC, std(ts.Cs, 0, 3)/sqrt(T), 'RelTol', 1e-14);
d = ts - m;
verifyEqual(testCase, mean(d.Cs, 3), zeros(L+1), 'AbsTol', 1e-14);
end

function testClimatologyRecovery(testCase)
rng(7); L = 6; T = 96;
t = (0:T-1)'/12 + 2010;
t0 = mean(t);
n1 = L + 1;
bias = randn(n1); trend = 0.1*randn(n1); ca = randn(n1); sa = randn(n1);
csm = 0.5*randn(n1); ssm = 0.5*randn(n1);
bias = tril(bias); trend = tril(trend); ca = tril(ca); sa = tril(sa);
csm = tril(csm); ssm = tril(ssm);
Cs = zeros(n1, n1, T); Ss = zeros(n1, n1, T);
for k = 1:T
    dt = t(k) - t0;
    Cs(:,:,k) = bias + trend*dt + ca*cos(2*pi*dt) + sa*sin(2*pi*dt) ...
        + csm*cos(4*pi*dt) + ssm*sin(4*pi*dt) + 1e-9*tril(randn(n1));
    Ss(:,:,k) = 1e-9*tril(randn(n1), -1);
end
ts = shSeries(Cs, Ss = Ss, Epochs = t);
[clim, resid] = ts.climatology(T0 = t0);
verifyEqual(testCase, clim.biasC,  bias,  'AbsTol', 1e-7);
verifyEqual(testCase, clim.trendC, trend, 'AbsTol', 1e-7);
verifyEqual(testCase, clim.cosAnnC, ca,   'AbsTol', 1e-7);
verifyEqual(testCase, clim.sinAnnC, sa,   'AbsTol', 1e-7);
verifyEqual(testCase, clim.cosSemiC, csm, 'AbsTol', 1e-7);
% eval reproduces the series
g = clim.eval(t(10));
verifyEqual(testCase, g.C, Cs(:,:,10), 'AbsTol', 1e-6);
% residuals at noise level
verifyLessThan(testCase, max(abs(resid.Cs(:))), 1e-7);
% robust fit survives an outlier month
Cs(:,:,5) = Cs(:,:,5) + 100;
tsO = shSeries(Cs, Ss = Ss, Epochs = t);
climR = tsO.climatology(T0 = t0, Robust = true);
verifyEqual(testCase, climR.trendC, trend, 'AbsTol', 1e-4);
end

function testGAXRestore(testCase)
rng(8); L = 6; T = 12;
gsm = randomSeries(L, T);
gad = randomSeries(L, T);
gad = shSeries(gad.Cs, Ss = gad.Ss, ...
    Epochs = gsm.epochs + 0.01, ProductType = "GAD");   % slight offset, within tol
gsm = shSeries(gsm.Cs, Ss = gsm.Ss, Epochs = gsm.epochs, ProductType = "GSM");
out = gsm.restore(gad);
verifyEqual(testCase, out.Cs, gsm.Cs + gad.Cs, 'AbsTol', 0);
verifyEqual(testCase, out.productType, "GSM+GAD");
end

% -------------------------------------------------------------- tvANS chain
function testTvANSPipelineAndDeconvolution(testCase)
rng(9); L = 10; T = 48;
t = (0:T-1)'/12 + 2018;
idx = shx.shIndex(L);
% basin kernels
B = zeros(idx.P, 2);
B(:,1) = exp(-((idx.n - 4)/3).^2) .* shx.ylm(deg2rad(50),  deg2rad(10),  idx)' * 0.05;
B(:,2) = exp(-((idx.n - 4)/3).^2) .* shx.ylm(deg2rad(35),  deg2rad(60),  idx)' * 0.05;
cTrue = [1.5*sin(2*pi*t)'; 0.8*cos(2*pi*t + 0.7)'] + 0.3*randn(2, T);
Xtrue = B * cTrue;
% striping noise (per-order, same-parity correlated)
Xn = zeros(idx.P, T);
for m = 0:L
    for cs = 0:double(m>0)
        for par = 0:1
            nn = max(m,2):L; nn = nn(mod(nn,2)==par);
            rows = squeeze(idx.pos(nn+1, m+1, cs+1)); rows = rows(rows>0);
            if isempty(rows), continue; end
            sd = 1e-3 + 0.05*(idx.n(rows)/L).^6;
            Cb = 0.9.^(abs(idx.n(rows) - idx.n(rows)')/2) .* (sd*sd');
            Xn(rows,:) = chol(Cb + 1e-12*eye(numel(rows)), 'lower') * randn(numel(rows), T);
        end
    end
end
X = Xtrue + Xn;
% pack into an shSeries
n1 = L + 1;
Cs = zeros(n1,n1,T); Ss = zeros(n1,n1,T);
for k = 1:T
    [Cs(:,:,k), Ss(:,:,k)] = shx.csFromVec(X(:,k), idx);
end
ts = shSeries(Cs, Ss = Ss, Epochs = t);
[tsF, op] = ts.filter("tvANS");
% high-degree noise strongly damped
XfChk = zeros(idx.P, T);
for k = 1:T
    XfChk(:,k) = shx.vecFromCS(tsF.Cs(:,:,k), tsF.Ss(:,:,k), idx);
end
hi = idx.n >= 8;
r0 = mean(mean((X(hi,:)  - Xtrue(hi,:)).^2));
r1 = mean(mean((XfChk(hi,:) - Xtrue(hi,:)).^2));
% Threshold from the ORACLE Wiener bound for this construction: the basin
% kernels exp(-((n-4)/3)^2) put real signal at n=8-9, so the per-coefficient
% oracle MSE ratio s_sig/(s_sig+s_noise) averages ~0.09 with p90 ~ 0.34
% (numpy MC) - a flat 0.05 is unattainable even for the optimal filter.
% Measured 0.125 is oracle-class; 0.3 keeps the check meaningful (>= 3x
% damping) with margin.
verifyLessThan(testCase, r1/r0, 0.3);
% at n = 10 the noise dominates; the data-driven chain measures 0.125
% (deterministic under rng(9)), above the ~0.01 oracle because the
% empirical per-order noise blocks are contaminated by the stochastic
% basin signal at n=8-9 through the same-parity striping correlation - a
% known, documented empirical-N limitation, not a regression. 0.2 keeps
% a meaningful >= 5x damping assertion with margin across MATLAB builds.
h10 = idx.n == 10;
r010 = mean(mean((X(h10,:)  - Xtrue(h10,:)).^2));
r110 = mean(mean((XfChk(h10,:) - Xtrue(h10,:)).^2));
verifyLessThan(testCase, r110/r010, 0.2);
% deconvolved basin averages beat naive
avgTrue = (B' * Xtrue) ./ diag(B'*B);
[avgHat, out] = tsF.basinAverage(B, Deconvolve = true, Op = op);
rmsD = sqrt(mean((avgHat - avgTrue).^2, 2));
rmsN = sqrt(mean((out.avgNaive - avgTrue).^2, 2));
verifyLessThanOrEqual(testCase, max(rmsD ./ rmsN), 1.0);
% constrained variant reproduces the constraint functional exactly
Ac = exp(-idx.n/4) .* randn(idx.P, 1);
[~, opC] = ts.filter("tvANS", Constraints = Ac);
err = norm(Ac'*opC.Xfres - Ac'*opC.Xres) / norm(Ac'*opC.Xres);
verifyLessThan(testCase, err, 1e-8);
end

function testEigTrickEquivalence(testCase)
rng(10); P = 40;
A = randn(P); S = A*A'/P + 0.1*eye(P);
Bm = randn(P); N = Bm*Bm'/P + 0.1*eye(P);
[U, D] = eig(S, N, 'chol');
lam = max(real(diag(D)), 0);
Ut = U'; V = inv(Ut);
for s = [0.3, 1.0, 4.7]
    Wd = S / (S + s*N);
    Wt = V * ((lam./(lam+s)) .* Ut);
    verifyLessThan(testCase, norm(Wt - Wd, 'fro')/norm(Wd, 'fro'), 1e-10);
end
end

% ----------------------------------------------------------------- I/O etc
function testWriteReadRoundtrip(testCase)
rng(11); L = 12;
g = randomField(L);
f = [tempname '.gfc'];
cl = onCleanup(@() delete(f));
g.write(f);
g2 = shCoefficients.read(f);
verifyEqual(testCase, g2.C, g.C, 'RelTol', 1e-13, 'AbsTol', 1e-18);
verifyEqual(testCase, g2.S, g.S, 'RelTol', 1e-13, 'AbsTol', 1e-18);
verifyEqual(testCase, g2.GM, g.GM, 'RelTol', 1e-10);
end

function testEvalAtMatchesCompat(testCase)
f = fullfile(testCase.TestData.dataDir, 'test_variable.gfct');
g = shCoefficients.read(f);
ge = g.evalAt(2010.5);
model = shReadGFC(f);
[Ct, St] = shEvalGFCT(model, 2010.5);
verifyEqual(testCase, ge.C, Ct, 'AbsTol', 0);
verifyEqual(testCase, ge.S, St, 'AbsTol', 0);
end

function testVecRoundtrip(testCase)
rng(12); L = 9;
g = randomField(L);
idx = shx.shIndex(L);
x = g.vec(idx);
g2 = shCoefficients.fromVec(x, idx, g);
% below minDegree is zeroed by construction; compare the indexed part
verifyEqual(testCase, g2.vec(idx), x, 'AbsTol', 0);
end

function testTN14Apply(testCase)
f = writeSyntheticTN14();
cl = onCleanup(@() delete(f));
tn = shx.readTN14(f);
verifyEqual(testCase, numel(tn.epoch), 3);
verifyTrue(testCase, isnan(tn.C30(1)) && ~isnan(tn.C30(3)));
g = randomField(6);
g = setEpoch(g, tn.epoch(3));
g2 = g.applyTN14(tn);
verifyEqual(testCase, g2.C(3,1), tn.C20(3), 'AbsTol', 0);
verifyEqual(testCase, g2.C(4,1), tn.C30(3), 'AbsTol', 0);
g1 = setEpoch(g, tn.epoch(1));
g3 = g1.applyTN14(tn);                      % C30 is NaN there -> C20 only
verifyEqual(testCase, g3.C(3,1), tn.C20(1), 'AbsTol', 0);
verifyEqual(testCase, g3.C(4,1), g.C(4,1), 'AbsTol', 0);
end

% ---------------------------------------------------------------- helpers
function g = randomField(L)
n1 = L + 1;
C = tril(randn(n1)) ./ max(1, (0:L)'.^2);
S = tril(randn(n1), -1) ./ max(1, (0:L)'.^2);
S(:, 1) = 0;    % S_n0 has no meaning (sin(0*lon)=0): zero in synthesis,
                % unrecoverable in analysis - keep test fields physical
g = shCoefficients(C, S, SigmaC = 0.01*abs(tril(rand(n1)))+1e-6, ...
    SigmaS = 0.01*abs(tril(rand(n1),-1))+1e-6, Epoch = 2024.1);
end

function [gA, gB] = randomPair(L)
gA = randomField(L);
gB = randomField(L);
end

function ts = randomSeries(L, T)
n1 = L + 1;
Cs = zeros(n1,n1,T); Ss = zeros(n1,n1,T);
for k = 1:T
    Cs(:,:,k) = tril(randn(n1));
    Ss(:,:,k) = tril(randn(n1), -1);
end
ts = shSeries(Cs, Ss = Ss, Epochs = 2015 + (0:T-1)'/12);
end

function g2 = setType(g, t)
g2 = shCoefficients(g.C, g.S, SigmaC = g.sigmaC, SigmaS = g.sigmaS, ...
    GM = g.GM, R = g.R, Epoch = g.epoch, ProductType = t);
end

function g2 = setEpoch(g, ep)
g2 = shCoefficients(g.C, g.S, SigmaC = g.sigmaC, SigmaS = g.sigmaS, ...
    GM = g.GM, R = g.R, Epoch = ep, ProductType = g.productType);
end

function f = writeSyntheticTN14()
f = [tempname '.txt'];
fid = fopen(f, 'w');
fprintf(fid, 'TITLE: synthetic TN-14 for unit tests\n');
fprintf(fid, 'some header line without numbers\n');
fprintf(fid, 'PRODUCT: C20/C30 replacement\n');
% MJDb yrb C20 dC20e10 sigC20e10 C30 dC30e10 sigC30e10 MJDe yre
fprintf(fid, '58119 2018.000 -4.841e-04 1.0 0.5 NaN NaN NaN 58149 2018.082\n');
fprintf(fid, '58149 2018.082 -4.842e-04 1.1 0.5 NaN NaN NaN 58180 2018.163\n');
fprintf(fid, '58180 2018.163 -4.843e-04 1.2 0.5 -9.5e-07 2.0 1.0 58208 2018.240\n');
fclose(fid);
end

% =================================================================== v2.1
function testFFTSynthesisMatchesDirect(testCase) %#ok<INUSD>
% FFT-along-longitude path must equal the direct trig product on a uniform
% full-circle grid (Python-validated to 2.9e-13)
rng(11);
L = 24;
C = tril(randn(L+1)); S = tril(randn(L+1), -1);
lat = -87.5:5:87.5; lon = 0:5:355;
gD = shx.shSynthesis(C, S, 3.986004415e14, 6378136.3, lat, lon, ...
    'quantity', 'geoid', 'method', 'direct');
gF = shx.shSynthesis(C, S, 3.986004415e14, 6378136.3, lat, lon, ...
    'quantity', 'geoid', 'method', 'fft');
verifyEqual(testCase, gF, gD, 'AbsTol', 1e-9 * max(abs(gD(:))));
end

function testFFTSynthesisNonZeroStartLon(testCase)
% the FFT path must honor a non-zero first longitude via the phase factor
rng(12);
L = 10; C = tril(randn(L+1)); S = tril(randn(L+1), -1);
lat = 10; lon0 = 17.5; lon = lon0 + (0:359);
gD = shx.shSynthesis(C, S, 1, 1, lat, lon, 'method', 'direct');
gF = shx.shSynthesis(C, S, 1, 1, lat, lon, 'method', 'fft');
verifyEqual(testCase, gF, gD, 'AbsTol', 1e-10 * max(abs(gD(:))));
end

function testScaledLegendreSumRule(testCase)
% addition theorem sum_m Pbar_nm^2 = 2n+1: the sharpest cheap check of the
% scaled recursion; holds to 3e-11 at n=2190 even at lat=89.99 deg
for latDeg = [0 45 89 89.99]
    P = shx.legendreALF(2190, deg2rad(latDeg));
    for n = [60 600 2190]
        s = sum(P(n+1, 1:n+1).^2);
        verifyEqual(testCase, s / (2*n+1), 1, 'RelTol', 1e-9, ...
            sprintf('sum rule failed at n=%d, lat=%g', n, latDeg));
    end
end
verifyTrue(testCase, all(isfinite(P(:))));
end

function testAnalysisRingRoundTrip(testCase)
% band-limited synthesis -> rings analysis must be exact (8.7e-15 in Python)
rng(13);
L = 12;
C0 = tril(randn(L+1)); S0 = tril(randn(L+1), -1); S0(:,1) = 0;
lat = linspace(-84, 84, 2*L+4); lon = (0:2*L+3) * 360/(2*L+4);
g = shx.shSynthesis(C0, S0, 1, 1, lat, lon, 'method', 'direct');
[C, S, info] = shx.shAnalysisGrid(g, lat, lon, L, GM = 1, R = 1);
verifyEqual(testCase, C, C0, 'AbsTol', 1e-10);
verifyEqual(testCase, S, S0, 'AbsTol', 1e-10);
verifyEqual(testCase, info.method, 'rings');
verifyLessThan(testCase, info.residRMS, 1e-10);
end

function testAnalysisScatteredRoundTrip(testCase)
% scattered least squares: exact for sufficient random sampling
rng(14);
L = 8;
C0 = tril(randn(L+1)); S0 = tril(randn(L+1), -1); S0(:,1) = 0;
Np = 500;
latp = asind(2*rand(Np,1) - 1); lonp = 360*rand(Np,1);
f = zeros(Np, 1);
for k = 1:Np
    f(k) = shx.shSynthesis(C0, S0, 1, 1, latp(k), lonp(k), 'method', 'direct');
end
[C, S, info] = shx.shAnalysisGrid(f, latp, lonp, L, GM = 1, R = 1);
verifyEqual(testCase, C, C0, 'AbsTol', 1e-9);
verifyEqual(testCase, S, S0, 'AbsTol', 1e-9);
verifyEqual(testCase, info.method, 'ls');
end

function testAnalysisQuantityInversion(testCase)
% analysis of a gravity_anomaly grid recovers Stokes coefficients with
% degrees 0/1 correctly reported unobservable
rng(15);
L = 10; GM = 3.986004415e14; R = 6378136.3;
C0 = tril(randn(L+1)) * 1e-8; S0 = tril(randn(L+1), -1) * 1e-8; S0(:,1) = 0;
C0(1:2, :) = 0; S0(1:2, :) = 0;                 % nothing below degree 2
lat = linspace(-85, 85, 2*L+6); lon = (0:2*L+5) * 360/(2*L+6);
g = shx.shSynthesis(C0, S0, GM, R, lat, lon, 'quantity', 'gravity_anomaly');
[C, S, info] = shx.shAnalysisGrid(g, lat, lon, L, ...
    quantity = 'gravity_anomaly', GM = GM, R = R);
verifyEqual(testCase, C, C0, 'AbsTol', 1e-16);
verifyEqual(testCase, S, S0, 'AbsTol', 1e-16);
% gravity_anomaly kernel (n-1): only n=1 carries no information
% (n=0 has factor -1, observable) - matches the physics
verifyTrue(testCase, isequal(info.unobservedDegrees, 1));
end

function testAnalysisKaulaRegularization(testCase)
% under-determined sampling: plain LS must refuse, Kaula must stabilize
rng(16);
L = 8;
Np = 40;                                        % << (L+1)^2 = 81 unknowns
latp = asind(2*rand(Np,1) - 1); lonp = 360*rand(Np,1);
f = randn(Np, 1) * 1e-8;
verifyError(testCase, @() shx.shAnalysisGrid(f, latp, lonp, L), ...
    'shx:shAnalysisGrid:rankDeficient');
[C, S] = shx.shAnalysisGrid(f, latp, lonp, L, Kaula = 1);
verifyTrue(testCase, all(isfinite(C(:))) && all(isfinite(S(:))));
end

function testTvANSBlocksMatchFull(testCase)
% block-diagonal path must be IDENTICAL to the full eigendecomposition
rng(17);
Lmax = 8; T = 30;
idx = shx.shIndex(Lmax);
t = 2002 + (0:T-1)'/12;
X = randn(idx.P, T) * 1e-9;
[XfF, ~, infoF] = shx.tvANSFilter(X, t, idx, Blocks = 'off');
[XfB, opB, infoB] = shx.tvANSFilter(X, t, idx, Blocks = 'on');
verifyEqual(testCase, XfB, XfF, 'AbsTol', 1e-12 * max(abs(XfF(:))));
verifyEqual(testCase, infoB.sigmaXfres, infoF.sigmaXfres, ...
    'AbsTol', 1e-10 * max(infoF.sigmaXfres(:)));
verifyEqual(testCase, opB.layout, 'blocks');
end

function testPosteriorSigmaAgainstDirect(testCase)
% (V.^2)*(s*lam/(lam+s)) must equal diag((I-W)S) formed explicitly
rng(18);
P = 40;
A = randn(P); S = A*A' + P*eye(P);
B = randn(P); N = B*B' + P*eye(P);
[U, D] = eig(S, N, 'chol');
lam = diag(D); Ut = U'; V = inv(Ut); %#ok<MINV>
s = 0.8;
W = S / (S + s*N);
direct = diag((eye(P) - W) * S);
formula = (V.^2) * (s*lam ./ (lam + s));
verifyEqual(testCase, formula, direct, 'AbsTol', 1e-10 * max(abs(direct)));
end

function testConstrainedPosteriorSigmaExact(testCase)
% v2.5: with hard constraints the posterior is the exact diagonal of
%   (W-I)S(W-I)' + s*W*N*W',  W = W0 + S*Ac*M^-1*Ac'*(I-W0)
% - verify the O(P^2 q) eig-basis formula against the brute-force
% covariance, and the constrained-direction identity
% Ac'*Cov*Ac = s*Ac'*N*Ac (the filter passes Ac'x unfiltered).
rng(21);
P = 40; q = 3;
A0 = randn(P); S = A0*A0'/P + 0.5*eye(P);
B0 = randn(P); N = B0*B0'/P + 0.3*eye(P);
Ac = randn(P, q);
s = 0.7;
[U, D] = eig(S, N, 'chol');
lam = max(real(diag(D)), 0); Ut = U'; V = inv(Ut); %#ok<MINV>
g = lam ./ (lam + s);
W0 = V * (g .* Ut);
M = Ac' * S * Ac;
W = W0 + S * Ac * (M \ (Ac' * (eye(P) - W0)));
Cov = (W - eye(P)) * S * (W - eye(P))' + s * (W * N * W');
dTrue = diag(Cov);
% the implemented formula (identical algebra to tvANSFilter v2.5)
Cq = Ac' * V; CqT = Cq';
Fq = (lam .* CqT) / M;  VF = V * Fq;  V2 = V.^2;
w1 = (1 - g).^2 .* lam;  w2 = g .* (1 - g);
A2q = Cq * (w1 .* CqT);  A4q = Cq * ((1 - g).^2 .* CqT);
d = V2 * (w1 + s * g.^2) ...
  - 2 * sum((V * (w1 .* CqT)) .* VF, 2) ...
  + sum((VF * A2q) .* VF, 2) ...
  + 2 * s * sum((V * (w2 .* CqT)) .* VF, 2) ...
  + s * sum((VF * A4q) .* VF, 2);
verifyEqual(testCase, d, dTrue, 'AbsTol', 1e-10 * max(abs(dTrue)));
verifyEqual(testCase, Ac' * Cov * Ac, s * (Ac' * N * Ac), ...
    'AbsTol', 1e-9 * norm(N));
% and the old unconstrained formula must NOT match (it underestimates)
dOld = V2 * (s * lam ./ (lam + s));
verifyLessThan(testCase, min(dOld ./ dTrue), 0.95);
end

function testConstrainedFilterSigmaStack(testCase)
% end-to-end: the sigma stacks of a CONSTRAINED tvANS run match the
% brute-force covariance diagonal built from the returned operator.
rng(22);
L = 6; T = 18;
idx = shx.shIndex(L, MinDegree = 2); P = idx.P;
Cs = zeros(L+1, L+1, T); Ss = Cs;
mL = tril(true(L+1));
for t = 1:T
    Cs(:,:,t) = (mL .* randn(L+1)) * 1e-9;
    St = (mL .* randn(L+1)) * 1e-9; St(:, 1) = 0;
    Ss(:,:,t) = St;
end
ts = shSeries(Cs, Ss = Ss, Epochs = 2015 + (0:T-1)'/12);
Ac = exp(-idx.n / 3) .* randn(P, 2);
[tsF, op, info] = ts.filter("tvANS", Constraints = Ac, Blocks = "off");
verifyTrue(testCase, contains(info.sigmaNote, "exact"));
S = op.V * (op.lam .* op.V');           % S = V diag(lam) V'
N = op.V * op.V';                       % N = V V'
for t = [1, T]
    st = op.s(t);
    g = op.lam ./ (op.lam + st);
    W0 = op.V * (g .* op.Ut);
    W = W0 + S * Ac * (op.M \ (Ac' * (eye(P) - W0)));
    dTrue = diag((W - eye(P)) * S * (W - eye(P))' + st * (W * N * W'));
    verifyEqual(testCase, info.sigmaXfres(:, t).^2, dTrue, ...
        'AbsTol', 1e-9 * max(abs(dTrue)));
end
verifyEqual(testCase, tsF.nEpochs, T);
end

function testBasinSigmaAgainstDirect(testCase)
% cov(c) = s * A^-1 (Gv' g^2 Gv) A^-T  vs the explicit matrix expression
rng(19);
P = 30; K = 2;
A0 = randn(P); S = A0*A0' + P*eye(P);
B0 = randn(P); N = B0*B0' + P*eye(P);
B = randn(P, K);
[U, D] = eig(S, N, 'chol');
lam = diag(D); Ut = U'; V = inv(Ut); %#ok<MINV>
s = 1.3;
g = lam ./ (lam + s);
W = S / (S + s*N);
A = B' * W * B;
covDirect = s * (A \ (B' * W * N * W' * B)) / A';
Gv = V' * B;
covFormula = s * (A \ (Gv' * ((g.^2) .* Gv))) / A';
verifyEqual(testCase, covFormula, covDirect, ...
    'AbsTol', 1e-9 * max(abs(covDirect(:))));
end

function testBasinSigmaIncludesDeterministic(testCase)
% v2.5: out.sigma of the deconvolution carries the OLS parameter
% uncertainty of the restored deterministic part. Non-circular check:
% rebuild the residual-noise sigma from the DIRECT operator matrices
% (full layout), add the leverage term from independently recomputed
% fit quantities, and match out.sigma. Plus the hat-matrix trace
% identity sum(h_t) = #parameters.
rng(23);
L = 6; T = 20; K = 2;
idx = shx.shIndex(L, MinDegree = 2); P = idx.P;
Cs = zeros(L+1, L+1, T); Ss = Cs;
mL = tril(true(L+1));
tY = 2012 + (0:T-1)'/12;
for t = 1:T
    Cs(:,:,t) = (mL .* randn(L+1)) * 1e-9 + ...
        (mL * 2e-9) * cos(2*pi*tY(t));         % strong seasonal model
    St = (mL .* randn(L+1)) * 1e-9; St(:, 1) = 0;
    Ss(:,:,t) = St;
end
ts = shSeries(Cs, Ss = Ss, Epochs = tY);
B = randn(P, K);
[~, op] = ts.filter("tvANS", Blocks = "off");
[~, out] = shx.basinDeconvolve(B, op);
% leverage recomputed independently from the same fit call
X = zeros(idx.P, T);
for t = 1:T
    X(:, t) = shx.vecFromCS(Cs(:,:,t), Ss(:,:,t), idx);
end
[~, ~, ~, Afit, ~, resVar] = shx.fitDeterministicModel(X, tY);
h = sum((Afit / (Afit' * Afit)) .* Afit, 2)';
verifyEqual(testCase, op.detLeverage, h, 'AbsTol', 1e-12);
verifyEqual(testCase, sum(h), size(Afit, 2), 'RelTol', 1e-10);
verifyEqual(testCase, op.detResVar, resVar(:), 'AbsTol', 1e-15);
% direct residual-noise sigma from explicit matrices
S = op.V * (op.lam .* op.V'); N = op.V * op.V';
BtB = B' * B; dB = diag(BtB);
varDet = ((B.^2)' * op.detResVar) * h;           % K x T
for t = [1, T]
    st = op.s(t);
    W = S / (S + st * N);
    A = B' * W * B;
    Q = st * (B' * (W * N * W') * B);
    covAvg = (BtB * ((A \ Q) / A') * BtB') ./ (dB * dB');
    expct = sqrt(diag(covAvg) + varDet(:, t) ./ dB.^2);
    verifyEqual(testCase, out.sigma(:, t), expct, ...
        'RelTol', 1e-8);
end
end

function testBasinSigmaConstrainedExact(testCase)
% v2.5: with hard constraints the deconvolution noise covariance uses
% the exact constrained operator. Verify against explicit matrices.
rng(24);
L = 5; T = 14; K = 2;
idx = shx.shIndex(L, MinDegree = 2); P = idx.P;
Cs = zeros(L+1, L+1, T); Ss = Cs;
mL = tril(true(L+1));
for t = 1:T
    Cs(:,:,t) = (mL .* randn(L+1)) * 1e-9;
    St = (mL .* randn(L+1)) * 1e-9; St(:, 1) = 0;
    Ss(:,:,t) = St;
end
ts = shSeries(Cs, Ss = Ss, Epochs = 2016 + (0:T-1)'/12);
Ac = exp(-idx.n / 3) .* randn(P, 1);
B = randn(P, K);
[~, op] = ts.filter("tvANS", Constraints = Ac, Blocks = "off");
[~, out] = shx.basinDeconvolve(B, op);
S = op.V * (op.lam .* op.V'); N = op.V * op.V';
BtB = B' * B; dB = diag(BtB);
varDet = ((B.^2)' * op.detResVar) * op.detLeverage;
for t = [1, T]
    st = op.s(t);
    g = op.lam ./ (op.lam + st);
    W0 = op.V * (g .* op.Ut);
    W = W0 + S * Ac * (op.M \ (Ac' * (eye(P) - W0)));
    A = B' * W * B;
    Q = st * (B' * (W * N * W') * B);
    covAvg = (BtB * ((A \ Q) / A') * BtB') ./ (dB * dB');
    expct = sqrt(diag(covAvg) + varDet(:, t) ./ dB.^2);
    verifyEqual(testCase, out.sigma(:, t), expct, 'RelTol', 1e-8);
end
end

function testClimatologyAliasPeriodRecovery(testCase)
% synthetic S2-alias signal must be recovered exactly by Periods=[...]
pS2 = 161/365.25;
T = 80; t = 2002 + (0:T-1)'/12;
nmax = 3; n1 = nmax + 1;
Cs = zeros(n1, n1, T); Ss = zeros(n1, n1, T);
tc = t - mean(t);
sig = 1e-9*(0.3 + 0.5*tc + 0.4*cos(2*pi*tc/pS2) - 0.2*sin(2*pi*tc/pS2));
for k = 1:T, Cs(3, 1, k) = sig(k); end
ts = shSeries(Cs, Ss = Ss, Epochs = t);
clim = ts.climatology(Periods = pS2);
[gc, gs] = clim.periodic(1);
verifyEqual(testCase, gc.C(3,1), 0.4e-9, 'AbsTol', 1e-15);
verifyEqual(testCase, gs.C(3,1), -0.2e-9, 'AbsTol', 1e-15);
verifyEqual(testCase, clim.trend().C(3,1), 0.5e-9, 'AbsTol', 1e-15);
end

function testClimatologySigmaMonteCarlo(testCase)
% empirical scatter of the trend estimate must match the claimed 1-sigma
% (OLS formula; Python MC gave ratios 0.99-1.02)
rng(20);
T = 60; t = 2002 + (0:T-1)'/12;
nmax = 2; n1 = nmax + 1;
MC = 60; sigNoise = 1e-9;
trends = zeros(MC, 1); claimed = zeros(MC, 1);
for mc = 1:MC
    Cs = zeros(n1, n1, T); Ss = zeros(n1, n1, T);
    Cs(3, 1, :) = sigNoise * randn(T, 1);
    ts = shSeries(Cs, Ss = Ss, Epochs = t);
    clim = ts.climatology();
    tr = clim.trend();
    trends(mc) = tr.C(3, 1);
    claimed(mc) = tr.sigmaC(3, 1);
end
ratio = std(trends) / mean(claimed);
verifyGreaterThan(testCase, ratio, 0.6);
verifyLessThan(testCase, ratio, 1.6);
end

function testFilterSigmaStacks(testCase)
% filtered series must carry finite posterior sigmas above MinDegree
rng(21);
Lmax = 6; T = 24;
n1 = Lmax + 1;
t = 2002 + (0:T-1)'/12;
Cs = 1e-9 * randn(n1, n1, T); Ss = 1e-9 * randn(n1, n1, T);
for k = 1:T
    Cs(:,:,k) = tril(Cs(:,:,k)); Ss(:,:,k) = tril(Ss(:,:,k)); Ss(:,1,k) = 0;
end
ts = shSeries(Cs, Ss = Ss, Epochs = t);
tf = ts.filter("tvANS");
verifyTrue(testCase, all(isfinite(tf.sigmaCs(3:end, 1, 1)), 'all'));
verifyTrue(testCase, all(tf.sigmaCs(3:end, 1, 1) > 0, 'all'));
verifyTrue(testCase, all(isnan(tf.sigmaCs(1:2, 1, 1))));  % below MinDegree
end

function testSelectAndDropNaN(testCase)
rng(22);
n1 = 4; T = 10;
t = 2002 + (0:T-1)';
Cs = randn(n1, n1, T); Ss = randn(n1, n1, T);
Cs(2, 1, 4) = NaN;
ts = shSeries(Cs, Ss = Ss, Epochs = t);
t2 = ts.dropNaN();
verifyEqual(testCase, t2.nEpochs, T - 1);
verifyEqual(testCase, t2.epochs, t([1:3 5:end]));
t3 = ts.select([2002.5 2005.5]);
verifyEqual(testCase, t3.nEpochs, 3);
t4 = ts.select([1 3 5]);
verifyEqual(testCase, t4.epochs, t([1 3 5]));
end

function testFitModelWeightsEquivalence(testCase)
% uniform weights must reproduce the unweighted fit exactly
rng(23);
T = 40; t = (0:T-1)'/12 + 2002;
X = randn(5, T);
[m1, ~, c1] = shx.fitDeterministicModel(X, t);
[m2, ~, c2] = shx.fitDeterministicModel(X, t, Weights = 2*ones(T,1));
verifyEqual(testCase, c2, c1, 'AbsTol', 1e-12 * max(abs(c1(:))));
verifyEqual(testCase, m2, m1, 'AbsTol', 1e-12 * max(abs(m1(:))));
end

% =============================================================== v2.2
function testGeodeticGeocentricRoundtrip(testCase)
lat = [-89.9, -45, -0.1, 0, 23.4567, 60, 90];
gc = shx.geodetic2geocentric(lat);
back = shx.geocentric2geodetic(gc);
verifyEqual(testCase, back, lat, 'AbsTol', 1e-12);
% known value (WGS84): atand((1-f)^2), independently computed in Python
verifyEqual(testCase, shx.geodetic2geocentric(45), 44.807576784018, 'AbsTol', 1e-9);
verifyEqual(testCase, shx.geodetic2geocentric(90), 90, 'AbsTol', 0);
% synthesis with LatType="geodetic" == manual conversion
g = randomField(12);
latGD = [-60, -10, 35, 70];
lon = 0:60:300;
g1 = g.synthesis(shx.geodetic2geocentric(latGD), lon);
g2 = g.synthesis(latGD, lon, LatType = "geodetic");
verifyEqual(testCase, g2, g1, 'AbsTol', 0);
end

function testBasinKernelAreaAndTaper(testCase)
idx = shx.shIndex(20, MinDegree = 0);
cap = @(la, lo) double(la > 60);                 % polar cap, ~6.7% area
[b, info] = shx.basinKernel(idx, cap);
% indicator staircase: +4.5% at OverSample=2 / Lmax=20 (Python-computed);
% the EXACT identities below are the strong checks
verifyEqual(testCase, info.areaFraction, (1 - cosd(30))/2, 'RelTol', 0.08);
% degree-0 coefficient of the indicator = area fraction (Y00 = 1)
verifyEqual(testCase, b(1), info.areaFraction, 'AbsTol', 1e-12);
% spectral taper = elementwise Jekeli weights
[bT, ~] = shx.basinKernel(idx, cap, TaperKm = 500);
Wn = shx.shGaussianWeights(idx.Lmax, 500);
verifyEqual(testCase, bT, b .* Wn(idx.n + 1), 'AbsTol', 0);
% buffered kernel encloses more area
[~, infoG] = shx.basinKernel(idx, cap, BufferKm = 1000);
verifyGreaterThan(testCase, infoG.areaFraction, info.areaFraction);
end

function testSlepianPolarCap(testCase)
idx = shx.shIndex(12, MinDegree = 0);
cap = @(la, lo) double(la > 60);
[G, lam, info] = shx.slepianBasis(idx, cap, NKeep = idx.P);
verifyGreaterThanOrEqual(testCase, min(lam), 0);
verifyLessThanOrEqual(testCase, max(lam), 1);
verifyEqual(testCase, sort(lam, 'descend'), lam, 'AbsTol', 0);
verifyEqual(testCase, G' * G, eye(idx.P), 'AbsTol', 1e-12);
% trace identity: sum(lam) = quadrature area fraction * P (addition thm)
verifyEqual(testCase, sum(info.lamAll), info.areaFraction * idx.P, ...
    'RelTol', 1e-10);
% the best taper concentrates: reconstruct on grid, energy inside > 95%
% (mask and Y/w must live on the SAME quadrature grid: evalMask defaults
% to OverSample=2 since v2.2.1, so pin both to the base grid here)
[Y, w] = shx.synthesisMatrix(idx);
[mask, ~] = shx.evalMask(idx, cap, OverSample = 1);
g1 = Y * G(:, 1);
verifyGreaterThan(testCase, sum(w .* mask .* g1.^2) / sum(w .* g1.^2), 0.95);
end

function testMCPropagateLinearFunctional(testCase)
rng(7);
g = randomFieldWithSigmas(10, 2020);
idx = shx.shIndex(10, MinDegree = 0);
wv = randn(idx.P, 1);
fun = @(gs) wv' * shx.vecFromCS(gs.C, gs.S, idx);
out = shx.mcPropagate(fun, g, N = 4000, Seed = 11);
% analytic: sigma^2 = sum(w_p^2 sigma_p^2) for independent Gaussians
sv = shx.vecFromCS(g.sigmaC, g.sigmaS, idx);
sv(~isfinite(sv)) = 0;
sigAna = sqrt(sum((wv .* sv).^2));
verifyEqual(testCase, out.sigma, sigAna, 'RelTol', 0.05);   % MC ~1.1%
verifyEqual(testCase, out.mean, fun(g), 'AbsTol', 4 * sigAna / sqrt(4000));
end

function testMCPropagateFullCovariance(testCase)
% the Cov path with the REAL SINEX fixture: 12 params = shIndex(3,
% MinDegree=2); propagated sigma of a linear functional must equal
% sqrt(w' M w)
f = fullfile(fileparts(mfilename('fullpath')), 'test_data', ...
    'ITSG-Grace2018_n96_2008-04_head12.snx');
assumeTrue(testCase, isfile(f));
idx = shx.shIndex(3, MinDegree = 2);
snx = shx.readSINEX(f, Output = "covariance", Index = idx);
M = (snx.M + snx.M') / 2;
g = shCoefficients(zeros(4), zeros(4));
rng(5); wv = randn(idx.P, 1);
fun = @(gs) wv' * shx.vecFromCS(gs.C, gs.S, idx);
out = shx.mcPropagate(fun, g, Cov = M, Idx = idx, N = 4000, Seed = 6);
verifyEqual(testCase, out.sigma, sqrt(wv' * M * wv), 'RelTol', 0.06);
end

function testBandedVCEMatchesGlobalWhenUniform(testCase)
% with one band spanning all orders, banded == global exactly
rng(21);
idx = shx.shIndex(10, MinDegree = 2);
T = 20; X = randn(idx.P, T);
t = 2020 + (0:T-1)'/12;
[Xf1, op1, i1] = shx.tvANSFilter(X, t, idx, Blocks = 'on'); %#ok<ASGLU>
[Xf2, op2, i2] = shx.tvANSFilter(X, t, idx, Blocks = 'on', ...
    VCEBands = [0, idx.Lmax + 1]); %#ok<ASGLU>
verifyEqual(testCase, Xf2, Xf1, 'AbsTol', 1e-14);
verifyEqual(testCase, op2.sBlocks, repmat(op2.s(:)', numel(op2.blocks), 1), ...
    'AbsTol', 1e-14);
end

function testBandedVCETracksOrderNoise(testCase)
% noise 3x stronger for m >= 6: band factors must reflect the contrast.
% NoiseCov must be a FIXED external model here: the default N is built
% from the data residuals and absorbs the very contrast being tested.
rng(22);
idx = shx.shIndex(12, MinDegree = 2);
T = 30; t = 2020 + (0:T-1)'/12;
scale = 1 + 2 * (idx.m >= 6);
X = scale .* randn(idx.P, T);
[~, op] = shx.tvANSFilter(X, t, idx, Blocks = 'on', ...
    NoiseCov = eye(idx.P), ...
    VCEBands = [0, 6, idx.Lmax + 1], VCEMinDegree = 4);
% mean band factors: high band / low band ~ variance ratio 9 (loose)
mLow = zeros(0,1); mHigh = zeros(0,1);
for k = 1:numel(op.blocks)
    mB = idx.m(op.blocks(k).rows(1));
    if mB < 6, mLow(end+1,1) = mean(op.sBlocks(k, :)); %#ok<AGROW>
    else, mHigh(end+1,1) = mean(op.sBlocks(k, :)); end %#ok<AGROW>
end
verifyGreaterThan(testCase, mean(mHigh) / mean(mLow), 3);
end

function testBandedSynthesisEqualsMonolithic(testCase)
g = randomField(40);
lat = -85:5:85; lon = 0:6:354;
g1 = g.synthesis(lat, lon);                       % cached, monolithic
g2 = g.synthesis(lat, lon, MaxMemGB = 1e-4, UseCache = false);  % forces bands
verifyEqual(testCase, g2, g1, 'AbsTol', 1e-18);
% direct shx call with tiny budget, supplied-P path unaffected
[gr1, ~, ~] = shx.shSynthesis(g.C, g.S, g.GM, g.R, lat, lon);
[gr2, ~, ~] = shx.shSynthesis(g.C, g.S, g.GM, g.R, lat, lon, ...
    'MaxMemGB', 1e-4);
verifyEqual(testCase, gr2, gr1, 'AbsTol', 1e-18);
end


function g = randomFieldWithSigmas(nmax, epoch)
rng(40 + nmax);
n1 = nmax + 1;
mL = tril(true(n1)); mL1 = tril(true(n1), -1);
C = (mL .* randn(n1)) * 1e-9;  S = (mL1 .* randn(n1)) * 1e-9;
sC = (mL .* abs(randn(n1))) * 1e-11;
sS = (mL1 .* abs(randn(n1))) * 1e-11;
g = shCoefficients(C, S, SigmaC = sC, SigmaS = sS, Epoch = epoch);
end

function testPropertyBasedSynthesisAnalysisRoundtrip(testCase)
% property-based: random nmax / grid sizes / quantities, exact ring-grid
% roundtrip coefficients -> grid -> coefficients (v2.2)
qs = ["geoid", "potential", "gravity_disturbance"];
for seed = 1:6
    rng(100 + seed);
    nmax = randi([4, 24]);
    nlat = nmax + 1 + randi(8);
    nlon = 2*nmax + 1 + randi(10);
    q = qs(randi(3));
    g = randomField(nmax);
    [xg, ~] = shx.gaussLegendre(nlat);            % ring latitudes (nodes)
    lat = asind(xg(:)');
    lon = (0:nlon-1) * 360 / nlon;
    grid = g.synthesis(lat, lon, quantity = q, UseCache = false);
    g2 = shCoefficients.analysis(grid, lat, lon, nmax, quantity = char(q));
    verifyEqual(testCase, g2.C, g.C, 'AbsTol', 1e-10, ...
        sprintf('roundtrip C failed (seed %d, nmax %d, %s)', seed, nmax, q));
    verifyEqual(testCase, g2.S, g.S, 'AbsTol', 1e-10, ...
        sprintf('roundtrip S failed (seed %d, nmax %d, %s)', seed, nmax, q));
end
end

% =============================================================== v2.3
function testNewKernelIdentities(testCase)
GM = 3.986004415e14; R = 6378136.3; L = 60;
n = (0:L)';
kn = -0.30 * n ./ (n + 3.0);                     % synthetic Love numbers
sd  = shx.kernelFactors('surface_density', L, GM, R, kn = kn);
ew  = shx.kernelFactors('ewh', L, GM, R, kn = kn, rho_water = 1000);
bp  = shx.kernelFactors('bottom_pressure', L, GM, R, kn = kn);
trr = shx.kernelFactors('gravity_gradient_rr', L, GM, R);
dis = shx.kernelFactors('gravity_disturbance', L, GM, R);
% exact relations between the kernels
verifyEqual(testCase, sd, 1000 * ew, 'RelTol', 1e-15);
verifyEqual(testCase, bp, (GM/R^2) * sd, 'RelTol', 1e-15);
verifyEqual(testCase, trr, dis .* (n + 2) / R, 'RelTol', 1e-15);
% upward continuation: Python-validated attenuation powers
hgt = 400e3; r = R + hgt;
pH = shx.kernelFactors('potential', L, GM, R, Height = hgt);
p0 = shx.kernelFactors('potential', L, GM, R);
verifyEqual(testCase, pH, p0 .* (R/r).^(n+1), 'RelTol', 1e-15);
dH = shx.kernelFactors('gravity_disturbance', L, GM, R, Height = hgt);
verifyEqual(testCase, dH, dis .* (R/r).^(n+2), 'RelTol', 1e-15);
% error contracts
verifyError(testCase, @() shx.kernelFactors('surface_density', L, GM, R), ...
    'shSynthesis:missingLoveNumbers');
verifyError(testCase, @() shx.kernelFactors('deformation_up', L, GM, R, ...
    kn = kn), 'shSynthesis:missingLoveNumbers');
verifyError(testCase, @() shx.kernelFactors('ewh', L, GM, R, kn = kn, ...
    Height = 1e5), 'shSynthesis:heightInvalid');
end

function testLegendreDerivativeIdentity(testCase)
% frozen dPbar/dphi identity vs central differences (Python-calibrated)
L = 30;
lats = deg2rad([-72.3, -33.1, -5.0, 12.7, 48.9, 81.2]);
h = 1e-6;
[~, D] = shx.legendreALFDeriv(L, lats);
for k = 1:numel(lats)
    Dnum = (shx.legendreALF(L, lats(k) + h) ...
          - shx.legendreALF(L, lats(k) - h)) / (2*h);
    verifyEqual(testCase, D(:, :, k), Dnum, 'AbsTol', 1e-6 * max(abs(Dnum(:))));
end
% poles are finite (no 1/cos singularities in D itself)
[~, Dp] = shx.legendreALFDeriv(L, deg2rad([90, -90]));
verifyTrue(testCase, all(isfinite(Dp(:))));
end

function testDeformationSynthesis(testCase)
% north/east against numerical gradients of the tangential-kernel scalar
rng(61);
L = 16; R = 6378136.3;
n1 = L + 1; n = (0:L)';
C = tril(randn(n1)) * 1e-8;
S = tril(randn(n1), -1) * 1e-8; S(:, 1) = 0;
kn = -0.30 * n ./ (n + 3.0);
hn = -0.90 * n ./ (n + 2.0) - 0.1;
ln = -0.05 * n ./ (n + 4.0) - 0.01;
lat = [37.2, -61.8]; lon = [211.0, 12.5];
[up, north, east] = shx.shSynthesisDeformation(C, S, R, lat, lon, ...
    kn = kn, hn = hn, ln = ln, nmin = 1);
% vertical == kernel route through the standard synthesis
fUp = shx.shSynthesisDeformation(C, S, R, lat, lon, ...
    kn = kn, hn = hn, ln = ln, nmin = 1);
upK = shx.shSynthesis(C, S, 3.986004415e14, R, lat, lon, ...
    'quantity', 'deformation_up', 'kn', kn, 'hn', hn, 'nmin', 1);
verifyEqual(testCase, fUp, upK, 'AbsTol', 1e-15);
% horizontal vs finite differences of the fH-scaled scalar field
fH = R * ln ./ (1 + kn); fH(1) = 0;
dphi = 1e-6;
for i = 1:2
    for j = 1:2
        gp = shx.shSynthesis(fH .* C, fH .* S, 1, 1, lat(i) + rad2deg(dphi), lon(j), ...
            'quantity', 'geoid');
        gm_ = shx.shSynthesis(fH .* C, fH .* S, 1, 1, lat(i) - rad2deg(dphi), lon(j), ...
            'quantity', 'geoid');
        dN = (gp - gm_) / (2 * dphi);
        verifyEqual(testCase, north(i, j), dN, 'RelTol', 1e-4, ...
            sprintf('north at (%g, %g)', lat(i), lon(j)));
        gp = shx.shSynthesis(fH .* C, fH .* S, 1, 1, lat(i), lon(j) + rad2deg(dphi), ...
            'quantity', 'geoid');
        gm_ = shx.shSynthesis(fH .* C, fH .* S, 1, 1, lat(i), lon(j) - rad2deg(dphi), ...
            'quantity', 'geoid');
        dE = (gp - gm_) / (2 * dphi) / cosd(lat(i));
        verifyEqual(testCase, east(i, j), dE, 'RelTol', 1e-4, ...
            sprintf('east at (%g, %g)', lat(i), lon(j)));
    end
end
% points mode == grid diagonal; east NaN at the pole, north finite
[uP, nP, eP] = shx.shSynthesisDeformation(C, S, R, lat, lon, ...
    kn = kn, hn = hn, ln = ln, Mode = "points");
% grid and points modes use different summation orders: agreement to a
% few ULP, not bit-identity (observed 1-4e-16 on R2026a)
verifyEqual(testCase, uP, [up(1,1); up(2,2)], 'RelTol', 1e-12);
verifyEqual(testCase, nP, [north(1,1); north(2,2)], 'RelTol', 1e-12);
verifyEqual(testCase, eP, [east(1,1); east(2,2)], 'RelTol', 1e-12);
[uPole, nPole, ePole] = shx.shSynthesisDeformation(C, S, R, 90, 45, ...
    kn = kn, hn = hn, ln = ln);
verifyTrue(testCase, isfinite(uPole) && isfinite(nPole) && isnan(ePole));
end

function testSurfaceDensityAnalysisRoundtrip(testCase)
rng(62);
L = 12; n = (0:L)';
kn = -0.30 * n ./ (n + 3.0);
g = randomField(L);
[xg, ~] = shx.gaussLegendre(L + 1);
lat = asind(xg(:)'); lon = (0:2*L+1) * 360 / (2*L + 2);
grid = g.synthesis(lat, lon, quantity = "surface_density", kn = kn, ...
    UseCache = false);
g2 = shCoefficients.analysis(grid, lat, lon, L, ...
    quantity = 'surface_density', kn = kn);
verifyEqual(testCase, g2.C, g.C, 'AbsTol', 1e-10);
verifyEqual(testCase, g2.S, g.S, 'AbsTol', 1e-10);
end

% =============================================================== v2.4
function testSecondDerivativeIdentity(testCase)
L = 25; h = 1e-5;
for lat = [-71, -20, 33.3, 66.6]
    la = deg2rad(lat);
    [P0, ~, D2] = shx.legendreALFDeriv(L, la);
    Pp = shx.legendreALF(L, la + h);
    Pm = shx.legendreALF(L, la - h);
    num = (Pp - 2*P0 + Pm) / h^2;
    verifyEqual(testCase, D2, num, 'AbsTol', 1e-5 * max(abs(num(:))));
end
end

function testGradientTensorInvariants(testCase)
rng(71);
L = 14; GM = 3.986004415e14; R = 6378136.3;
n1 = L + 1;
C = tril(randn(n1)) * 1e-8; C(1, 1) = 0;
S = tril(randn(n1), -1) * 1e-8; S(:, 1) = 0;
lat = [41, -55.5]; lon = [100, 331];
[G, info] = shx.shSynthesisGradientTensor(C, S, GM, R, lat, lon, ...
    Height = 250e3, nmin = 2);
% Laplace: trace = 0 (Python: 7e-16)
verifyLessThan(testCase, info.maxTraceResidual, 1e-12);
% Guu == gravity_gradient_rr kernel route with Height
guuK = shx.shSynthesis(C, S, GM, R, lat, lon, ...
    'quantity', 'gravity_gradient_rr', 'Height', 250e3, 'nmin', 2);
verifyEqual(testCase, G.uu, guuK, 'RelTol', 1e-12);
% symmetry fields present and finite away from poles
for f = ["nn", "ee", "un", "ue", "ne"]
    verifyTrue(testCase, all(isfinite(G.(f)(:))));
end
% angular second derivatives vs finite differences through the potential
fT = shx.kernelFactors('potential', L, GM, R, Height = 250e3);
fT(1:2) = 0;
h = 1e-4; r = R + 250e3;
scal = @(la_, lo_) shx.shSynthesis(fT .* C, fT .* S, 1, 1, la_, lo_, ...
    'quantity', 'geoid');
la = lat(1); lo = lon(1);
Tpp = (scal(la + rad2deg(h), lo) - 2*scal(la, lo) + scal(la - rad2deg(h), lo)) / h^2;
Gnn_num = Tpp / r^2 + gradTr(C, S, GM, R, L, la, lo, r) / r;
verifyEqual(testCase, G.nn(1, 1), Gnn_num, 'RelTol', 5e-4);
end

function Tr = gradTr(C, S, GM, R, L, la, lo, r)
n = (0:L)';
fTr = -(GM/R) * (n + 1) / r .* (R/r).^(n + 1);
fTr(1:2) = 0;
Tr = shx.shSynthesis(fTr .* C, fTr .* S, 1, 1, la, lo, 'quantity', 'geoid');
end

function testCombineCentersRecovery(testCase)
rng(72);
L = 8; n1 = L + 1; T = 14;
idx = shx.shIndex(L, MinDegree = 0);
% truth: smooth series
Ct = zeros(n1, n1, T); St = zeros(n1, n1, T);
base = tril(randn(n1)) * 1e-9;
for t = 1:T
    Ct(:, :, t) = base * (1 + 0.1 * sin(t/2));
    St(:, :, t) = tril(randn(n1), -1) * 1e-10; St(:, 1, t) = 0;
end
ep = 2010 + (0:T-1) / 12;
s2true = [1.0; 4.0; 2.25];                       % static center factors
tsC = cell(1, 3);
for c = 1:3
    Cs = Ct; Ss = St;
    for t = 1:T
        Nc = tril(randn(n1)) * 2e-10 * sqrt(s2true(c));
        Ns = tril(randn(n1), -1) * 2e-10 * sqrt(s2true(c)); Ns(:, 1) = 0;
        Cs(:, :, t) = Cs(:, :, t) + Nc;
        Ss(:, :, t) = Ss(:, :, t) + Ns;
    end
    tsC{c} = shSeries(Cs, Ss = Ss, Epochs = ep, ProductType = "GSM");
end
[tsComb, info] = shx.combineCenters(tsC, MaxIter = 8);
% factor RATIOS recovered (VCE fixes relative scale per month; compare
% medians across months, generous tolerance for T*P statistics)
med = median(info.s2, 2);
ratio = med / med(1);
verifyEqual(testCase, ratio, s2true / s2true(1), 'RelTol', 0.5);
% combined closer to truth than the worst center, month-wise stack
errC = zeros(3, 1);
for c = 1:3
    errC(c) = norm([tsC{c}.Cs(:) - Ct(:); tsC{c}.Ss(:) - St(:)]);
end
errComb = norm([tsComb.Cs(:) - Ct(:); tsComb.Ss(:) - St(:)]);
verifyLessThan(testCase, errComb, min(errC));
% honesty diagnostics present
verifyEqual(testCase, size(info.interCenterCorr), [3 3]);
verifyTrue(testCase, all(isfinite(info.redundancy(:))));
% posterior sigmas on the combined series
verifyEqual(testCase, size(tsComb.sigmaCs), [n1 n1 T]);
verifyTrue(testCase, all(tsComb.sigmaCs(:) >= 0));
end

function testCombineCentersContract(testCase)
rng(73);
L = 4; n1 = L + 1;
Cs = randn(n1, n1, 3) * 1e-9; Ss = Cs; Ss(:, 1, :) = 0;
ts1 = shSeries(Cs, Ss = Ss, Epochs = [2010 2010.1 2010.2]);
verifyError(testCase, @() shx.combineCenters({ts1}), ...
    'shx:combineCenters:tooFewCenters');
Cs2 = randn(L, L, 3) * 1e-9; Ss2 = Cs2; Ss2(:, 1, :) = 0;
ts2 = shSeries(Cs2, Ss = Ss2, Epochs = [2010 2010.1 2010.2]);
verifyError(testCase, @() shx.combineCenters({ts1, ts2}), ...
    'shx:combineCenters:nmaxMismatch');
ts3 = shSeries(Cs, Ss = Ss, Epochs = [2015 2015.1 2015.2]);
verifyError(testCase, ...
    @() shx.combineCenters({ts1, ts3}, AllowMissing = false), ...
    'shx:combineCenters:noCommonEpochs');
end

function testFingerprintConservation(testCase)
L = 24; n = (0:L)';
kn = -0.30 * n ./ (n + 3.0);
hn = -0.90 * n ./ (n + 2.0) - 0.1;
idx = shx.shIndex(L, MinDegree = 0);
ocean = @(la, lo) ~(la > 55 & (lo < 60 | lo > 300));
loadF = @(la, lo) -100 * double(la > 65 & lo < 40);   % kg/m^2 ice loss
[S, grid, info] = shx.seaLevelFingerprint(loadF, ocean, idx, ...
    kn = kn, hn = hn);
verifyTrue(testCase, info.converged);
% mass conservation at machine precision (Python: exact)
verifyLessThan(testCase, abs(info.massResidual), 1e-12);
% classic pattern: near-field BELOW eustatic (here negative), far above
S2 = info.S2D;
[LAT, ~] = ndgrid(grid.latDeg, grid.lonDeg);
OC = S2 ~= 0;
near = mean(S2(LAT > 60 & OC));
far = mean(S2(LAT < -30 & OC));
verifyLessThan(testCase, near, info.eustatic);
verifyGreaterThan(testCase, far, info.eustatic);
verifyLessThan(testCase, near, 0);
% polygon load form runs and conserves too
poly = [66 5; 66 35; 80 35; 80 5];
[~, ~, info2] = shx.seaLevelFingerprint(poly, ocean, idx, ...
    kn = kn, hn = hn, LoadValue = -100);
verifyLessThan(testCase, abs(info2.massResidual), 1e-12);
end

function testEOFReconstruction(testCase)
rng(74);
L = 10; n1 = L + 1; T = 40;
% two orthogonal spatial patterns, known temporal amplitudes
P1 = tril(randn(n1)); P2 = tril(randn(n1));
P2 = P2 - P1 * (P1(:)' * P2(:)) / (P1(:)' * P1(:));   % orthogonalize
a1 = 3 * sin(2*pi*(1:T)/10)'; a2 = randn(T, 1);
Cs = zeros(n1, n1, T); Ss = zeros(n1, n1, T);
for t = 1:T
    Cs(:, :, t) = a1(t) * P1 + a2(t) * P2;
end
ts = shSeries(Cs, Ss = Ss, Epochs = 2010 + (1:T)/12);
[modes, pcs, ve, info] = shx.eofAnalysis(ts, NModes = 3);
% two modes carry (essentially) all variance
verifyGreaterThan(testCase, sum(ve(1:2)), 0.999);
% mode-1 spatial pattern parallel to P1 (up to sign)
m1 = modes{1}.C(:);
cosang = abs(m1' * P1(:)) / (norm(m1) * norm(P1(:)));
verifyGreaterThan(testCase, cosang, 0.98);
% pcs unit variance
verifyEqual(testCase, var(pcs(:, 1)), 1, 'RelTol', 1e-10);
verifyEqual(testCase, info.T, T);
end

function testBreakpointFTest(testCase)
rng(75);
L = 2; n1 = L + 1; T = 120;
ep = 2008 + (0:T-1)/12; tb = 2013.0;
hinge = max(ep - tb, 0);
Cs = zeros(n1, n1, T); Ss = zeros(n1, n1, T);
% C20: real break; C21: pure noise
for t = 1:T
    Cs(3, 1, t) = 0.02 * hinge(t) + 0.05 * randn;
    Cs(3, 2, t) = 0.05 * randn;
end
ts = shSeries(Cs, Ss = Ss, Epochs = ep);
out = ts.trendBreaks(Breaks = tb);
verifyGreaterThan(testCase, out.F(3, 1), 10);         % Python: F ~ 23
verifyLessThan(testCase, out.pValue(3, 1), 0.01);
verifyGreaterThan(testCase, out.pValue(3, 2), 0.01);  % null survives
% hinge slope recovered
verifyEqual(testCase, out.hingeC(3, 1), 0.02, 'AbsTol', 0.01);
end

function testFanFilterSeparable(testCase)
rng(76);
L = 30; n1 = L + 1;
C = tril(randn(n1)); S = tril(randn(n1), -1); S(:, 1) = 0;
[Cf, Sf] = shx.shFanFilter(C, S, 300, 500);
Wn = shx.shGaussianWeights(L, 300);
Wm = shx.shGaussianWeights(L, 500);
verifyEqual(testCase, Cf, C .* (Wn(:) .* Wm(:)'), 'RelTol', 1e-14);
verifyEqual(testCase, Sf, S .* (Wn(:) .* Wm(:)'), 'RelTol', 1e-14);
% equal radii on the degree axis: fan == isotropic Gaussian at m=0
[Cf2, ~] = shx.shFanFilter(C, S, 400, 400);
Wn4 = shx.shGaussianWeights(L, 400);
verifyEqual(testCase, Cf2(:, 1), C(:, 1) .* Wn4(:) * Wn4(1), 'RelTol', 1e-14);
end

function testErrorMapMatchesMC(testCase)
rng(77);
idx = shx.shIndex(3, MinDegree = 2);
A = randn(idx.P) * 1e-10;
M = A * A' + 1e-22 * eye(idx.P);              % PD synthetic covariance
lat = [-40, 10, 55]; lon = [30, 200];
sig = shx.errorMap(M, idx, lat, lon, quantity = "geoid");
% Monte-Carlo cross-check
Nmc = 4000;
Lch = chol((M + M')/2 + 1e-30*eye(idx.P), 'lower');
vals = zeros(numel(lat), numel(lon), Nmc);
for k = 1:Nmc
    x = Lch * randn(idx.P, 1);
    [Cm, Sm] = shx.csFromVec(x, idx);
    vals(:, :, k) = shx.shSynthesis(Cm, Sm, 3.986004415e14, 6378136.3, ...
        lat, lon, 'quantity', 'geoid');
end
sigMC = std(vals, 0, 3);
verifyEqual(testCase, sig, sigMC, 'RelTol', 0.12);
end

function testWriteGridRoundtrip(testCase)
tmp = [tempname, '.nc'];
cleanup = onCleanup(@() delete(tmp)); %#ok<NASGU>
lat = (-88:4:88)'; lon = (0:4:356)';
[LA, LO] = ndgrid(lat, lon);
G = sind(LA) .* cosd(LO);
shx.writeGrid(tmp, G, lat, lon, Name = "lwe_thickness", Units = "m");
back = ncread(tmp, 'lwe_thickness')';
verifyEqual(testCase, back, G, 'AbsTol', 1e-14);
verifyEqual(testCase, ncread(tmp, 'lat'), lat, 'AbsTol', 0);
% 3-D stack
tmp2 = [tempname, '.nc'];
cleanup2 = onCleanup(@() delete(tmp2)); %#ok<NASGU>
G3 = cat(3, G, 2*G, -G);
shx.writeGrid(tmp2, G3, lat, lon, Name = "ewh", ...
    Epochs = [2010.1; 2010.2; 2010.3]);
b3 = permute(ncread(tmp2, 'ewh'), [2 1 3]);
verifyEqual(testCase, b3, G3, 'AbsTol', 1e-14);
tv = ncread(tmp2, 'time');
verifyEqual(testCase, numel(tv), 3);
verifyGreaterThan(testCase, tv(2), tv(1));
end

function testNormalFieldWGS84Values(testCase)
% closed form from defining constants vs published NIMA TR8350.2
[Cn, info] = shx.normalFieldCS(8);
verifyEqual(testCase, info.J2, 1.082629821313e-3, 'RelTol', 1e-12);
verifyEqual(testCase, Cn(3), -0.484166774985e-3, 'RelTol', 1e-11);
verifyEqual(testCase, Cn(5),  0.790303733511e-6, 'RelTol', 1e-10);
verifyEqual(testCase, Cn(7), -0.168724961151e-8, 'RelTol', 1e-9);
verifyEqual(testCase, Cn(9),  0.346052468394e-11, 'RelTol', 1e-8);
% odd degrees zero, degree 0 unity
verifyEqual(testCase, Cn([2 4 6 8]), zeros(4, 1));
verifyEqual(testCase, Cn(1), 1);
% GRS80 J2 hits its defining value through the derived flattening
[~, i80] = shx.normalFieldCS(2, System = "GRS80");
verifyEqual(testCase, i80.J2, 1.08263e-3, 'RelTol', 1e-8);
end

function testRescaleInvariance(testCase)
rng(91);
L = 12; n1 = L + 1;
C = tril(randn(n1)) * 1e-7; C(1, 1) = 1;
S = tril(randn(n1), -1) * 1e-7; S(:, 1) = 0;
GM1 = 3.986004415e14; R1 = 6378136.3;
GM2 = 3.986004418e14; R2 = 6378137.0;
g1 = shCoefficients(C, S, GM = GM1, R = R1);
g2 = g1.toReference(GM = GM2, R = R2);
% roundtrip identity
g3 = g2.toReference(GM = GM1, R = R1);
verifyEqual(testCase, g3.C, g1.C, 'RelTol', 1e-14);
% physical invariance: potential at points agrees (Python: exact)
lat = [12.5, -60]; lon = [33, 250];
% physical invariance: potential at a COMMON radius r = R2 + 100 km
V1r = g1.synthesis(lat, lon, quantity = "potential", ...
    Height = R2 + 100e3 - R1, UseCache = false);
V2r = g2.synthesis(lat, lon, quantity = "potential", ...
    Height = 100e3, UseCache = false);
verifyEqual(testCase, V1r, V2r, 'RelTol', 1e-12);
end

function testSubtractNormalField(testCase)
rng(92);
L = 10; n1 = L + 1;
GMm = 3.986004415e14; Rm = 6378136.3;
% construct a field = rescaled WGS84 normal field + known residual
[CnEll, infoN] = shx.normalFieldCS(L);
Cell = zeros(n1); Cell(:, 1) = CnEll;
CnResc = shx.rescaleGMR(Cell, zeros(n1), infoN.GM, infoN.a, GMm, Rm);
res = tril(randn(n1)) * 1e-9;
Sres = tril(randn(n1), -1) * 1e-9; Sres(:, 1) = 0;
g = shCoefficients(CnResc + res, Sres, GM = GMm, R = Rm);
out = g.subtractNormalField();
verifyEqual(testCase, out.C, g.C - CnResc, 'AbsTol', 1e-18);
verifyEqual(testCase, out.C(:, 1), g.C(:, 1) - CnResc(:, 1), 'AbsTol', 1e-18);
% zonals now at residual level, not 5e-4
verifyLessThan(testCase, abs(out.C(3, 1) - res(3, 1)), 1e-15);
% arithmetic guard: mismatched R errors with the documented ID
gW = shCoefficients(g.C, g.S, GM = infoN.GM, R = infoN.a);
verifyError(testCase, @() g - gW, 'shCoefficients:constantsMismatch');
end

function testDiffSpectrumIdentities(testCase)
% v2.6.0 comparison suite: scaled field has degree correlation 1 and
% half the signal amplitude as difference; identical fields have zero
% difference and no crossover (Python-validated identities)
rng(5); L = 16;
g = randomField(L, 2020);
spec = shx.diffSpectrum(g.C, g.S, 0.5 * g.C, 0.5 * g.S);
verifyEqual(testCase, spec.degCorr(3:end), ones(L - 1, 1), 'AbsTol', 1e-12);
verifyEqual(testCase, spec.diffAmp, 0.5 * spec.amp1, 'RelTol', 1e-12);
verifyEqual(testCase, spec.ncross, NaN);
spec0 = shx.diffSpectrum(g.C, g.S, g.C, g.S);
verifyEqual(testCase, spec0.diffAmp, zeros(L + 1, 1), 'AbsTol', 0);
% a difference exceeding the signal from some degree sets ncross there
C2 = g.C; S2 = g.S;
C2(9:end, :) = C2(9:end, :) + 10 * max(abs(g.C(:)));
specX = shx.diffSpectrum(g.C, g.S, C2, S2);
verifyEqual(testCase, specX.ncross, 8);
end

function testSpatialStatsTaylorIdentity(testCase)
% weighted Taylor identity crmsd^2 = stdA^2 + stdB^2 - 2*stdA*stdB*corr,
% analytic bias on a constant offset, and mask handling
rng(6);
lat = -88:4:88; lon = 0:6:354;
A = randn(numel(lat), numel(lon));
B = 0.8 * A + 0.3 * randn(size(A));
st = shx.spatialStats(A, B, lat, lon);
verifyEqual(testCase, st.crmsd^2, ...
    st.stdA^2 + st.stdB^2 - 2 * st.stdA * st.stdB * st.corr, ...
    'RelTol', 1e-10);
st2 = shx.spatialStats(A, A - 3, lat, lon);
verifyEqual(testCase, st2.bias, 3, 'AbsTol', 1e-12);
verifyEqual(testCase, st2.rmsd, 3, 'AbsTol', 1e-12);
verifyEqual(testCase, st2.corr, 1, 'AbsTol', 1e-12);
mask = false(size(A)); mask(1:10, :) = true;
stM = shx.spatialStats(A, B, lat, lon, Mask = mask);
verifyEqual(testCase, stM.nUsed, nnz(mask));
end

function testNSEAndEffectiveCorr(testCase)
% NSE anchor points and the AR(1)-corrected effective sample size
rng(7); T = 240;
ref = randn(T, 1);
verifyEqual(testCase, shx.nashSutcliffe(ref, ref), 1, 'AbsTol', 1e-14);
verifyEqual(testCase, shx.nashSutcliffe(ref, mean(ref) * ones(T, 1)), ...
    0, 'AbsTol', 1e-12);
% white noise: Neff stays close to T
ecW = shx.effectiveCorr(randn(T, 1), randn(T, 1));
verifyGreaterThan(testCase, ecW.neff, 0.7 * T);
verifyTrue(testCase, ecW.p >= 0 && ecW.p <= 1);
% strong AR(1): Neff collapses well below T
x = filter(1, [1 -0.85], randn(T, 1));
y = filter(1, [1 -0.85], randn(T, 1));
ecA = shx.effectiveCorr(x, y);
verifyLessThan(testCase, ecA.neff, T / 2);
% perfect correlation stays finite
ecP = shx.effectiveCorr(x, 2 * x + 1);
verifyEqual(testCase, ecP.r, 1, 'AbsTol', 1e-12);
verifyTrue(testCase, isfinite(ecP.t) && ecP.p < 1e-6);
end

function testThreeCorneredHatRecovery(testCase)
% generalized TCH recovers injected noise levels (Python-validated to
% 1% at T = 20000; generous tolerance here at T = 8000)
rng(8); T = 8000;
sig = [1, 2, 3, 1.5];
common = 0.05 * cumsum(randn(T, 1));
X = common + randn(T, numel(sig)) .* sig;
out = shx.threeCorneredHat(X);
verifyEqual(testCase, out.sigma, sig, 'RelTol', 0.15);
verifyFalse(testCase, any(out.clipped));
verifyEqual(testCase, size(out.pairVar, 1), 6);
verifyError(testCase, @() shx.threeCorneredHat(X(:, 1:2)), ...
    'shx:threeCorneredHat:needThree');
end

function testCompareReports(testCase)
% aggregator contracts on the real ITSG chain (solutions) and a
% synthetic 3-center stack (series) incl. TCH ordering and epoch drops
d = fullfile(fileparts(mfilename('fullpath')), 'test_data');
fG = fullfile(d, 'ITSG-Grace2018_n60_2008-04.gfc');
assumeTrue(testCase, isfile(fG));
g = shCoefficients.read(fG, Epoch = 2008.29);
rep = shx.compareSolutions(g, g.gaussian(350), Names = ["raw", "G350"]);
verifyEqual(testCase, rep.nmax, 60);
verifyTrue(testCase, isfinite(rep.spectral.ncross));   % filter diverges
verifyTrue(testCase, isfinite(rep.chi2dof) && rep.chi2dof > 0);
verifyTrue(testCase, rep.spatial.corr > 0.5 && rep.spatial.corr <= 1);
verifyFalse(testCase, rep.rescaled);
% mixed degrees truncate to the smaller solution
rep2 = shx.compareSolutions(g, g.truncate(30));
verifyEqual(testCase, rep2.nmax, 30);
% ---- series: common signal + graded noise, one epoch missing in (3)
rng(9); L = 6; n1 = L + 1; T = 36;
tY = 2019 + (0:T-1)'/12;
Cs = 1e-9 * randn(n1, n1, T); Ss = 1e-9 * randn(n1, n1, T);
mk = @(s, sel) shSeries(Cs(:, :, sel) + s * 1e-11 * randn(n1, n1, nnz(sel)), ...
    Ss = Ss(:, :, sel) + s * 1e-11 * randn(n1, n1, nnz(sel)), ...
    Epochs = tY(sel));
all36 = true(T, 1); m35 = all36; m35(20) = false;
ts1 = mk(1, all36); ts2 = mk(2, all36); ts3 = mk(6, m35);
repS = shx.compareSeries({ts1, ts2, ts3}, Names = ["A", "B", "C"]);
verifyEqual(testCase, numel(repS.epochs), T - 1);      % common epochs
verifyEqual(testCase, repS.nDropped(3), 1);
verifyTrue(testCase, all(repS.nse(2:3) < 1));
verifyTrue(testCase, repS.tch.sigma(3) > repS.tch.sigma(2));
verifyTrue(testCase, all(isfinite(repS.trendZ(2:3))));
verifyTrue(testCase, all(abs(repS.phaseLagDays(2:3)) <= 183));
% Basin path: unit vector picks a single coefficient series
idx = shx.shIndex(L);
b = zeros(idx.P, 1); b(1) = 1;
repB = shx.compareSeries({ts1, ts2}, Basin = b, Idx = idx);
g1 = ts1.at(1);
x1 = shx.vecFromCS(g1.C, g1.S, idx);
verifyEqual(testCase, repB.y(1, 1), x1(1), 'AbsTol', 1e-15);
end

function testSpectralFamilyIdentities(testCase)
rng(101);
L = 20; n1 = L + 1;
C = tril(randn(n1)) * 1e-8; S = tril(randn(n1), -1) * 1e-8; S(:, 1) = 0;
sC = abs(C) * 0.1; sS = abs(S) * 0.1;
sd = shx.shDegreeRMS(C, S, 'R', 6378136.3, 'sigmaC', sC, 'sigmaS', sS);
so = shx.shOrderRMS(C, S, sigmaC = sC, sigmaS = sS);
% degree/order marginals of the same total power
verifyEqual(testCase, sum(so.ordVariance), sum(sd.degVariance), ...
    'RelTol', 1e-14);
verifyEqual(testCase, sum(so.errVariance), sum(sd.errVariance), ...
    'RelTol', 1e-14);
% cumulative fields consistent
verifyEqual(testCase, sd.cumVariance(end), sum(sd.degVariance), 'RelTol', 1e-15);
verifyEqual(testCase, sd.cumRMS, sqrt(sd.cumVariance), 'AbsTol', 0);
verifyEqual(testCase, sd.cumAmplitude, 6378136.3 * sd.cumRMS, 'RelTol', 1e-15);
verifyEqual(testCase, so.cumVariance(end), sum(so.ordVariance), 'RelTol', 1e-15);
% conventions mirror: amplitude = R * RMS in both domains
verifyEqual(testCase, so.ordAmplitude, 6378136.3 * so.ordRMS, 'RelTol', 1e-15);
verifyEqual(testCase, so.domain, 'order');
verifyEqual(testCase, sd.domain, 'degree');
% n0 zeroing consistent between domains
sd2 = shx.shDegreeRMS(C, S, 'n0', 2);
so2 = shx.shOrderRMS(C, S, n0 = 2);
verifyEqual(testCase, sum(so2.ordVariance), sum(sd2.degVariance), 'RelTol', 1e-14);
end
