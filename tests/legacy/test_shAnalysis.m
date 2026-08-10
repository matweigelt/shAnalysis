classdef test_shAnalysis < matlab.unittest.TestCase
%TEST_SHANALYSIS Unit tests for legendreALF, shDegreeRMS, shSynthesis.
%   Run with: results = runtests('test_shAnalysis')
%
%   Claude (Sonnet 4.6), 2026-07-11.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

methods (Test)

    function testP00P10P11(testCase)
        lat = deg2rad([0 30 90]);
        P = legendreALF(1, lat);
        testCase.verifyEqual(squeeze(P(1,1,:)), ones(3,1), 'AbsTol', 1e-12);
        testCase.verifyEqual(squeeze(P(2,1,:)), sqrt(3)*sin(lat)', 'AbsTol', 1e-12);
        testCase.verifyEqual(squeeze(P(2,2,:)), sqrt(3)*cos(lat)', 'AbsTol', 1e-12);
    end

    function testOrthonormality(testCase)
        % Gauss-Legendre quadrature check of 4-pi normalization, analogous
        % to the independent Python cross-validation.
        nmax = 20;
        n_nodes = nmax + 5;
        beta = (1:n_nodes-1) ./ sqrt(4*(1:n_nodes-1).^2 - 1);
        T = diag(beta,1) + diag(beta,-1);
        [V,D] = eig(T);
        nodes = diag(D);
        weights = 2*(V(1,:)').^2;
        lat = asin(nodes);
        P = legendreALF(nmax, lat);
        for n = [0 5 10 20]
            for m = [0 min(n,3) n]
                lambdaFactor = 2*pi;
                if m ~= 0, lambdaFactor = pi; end
                Pnm = squeeze(P(n+1,m+1,:));
                integral = sum(weights .* Pnm.^2) * lambdaFactor;
                testCase.verifyEqual(integral, 4*pi, 'RelTol', 1e-8);
            end
        end
    end

    function testDegreeRMSFlatSpectrum(testCase)
        nmax = 5;
        C = zeros(nmax+1); S = zeros(nmax+1);
        C(3,1) = 1; % n=2,m=0
        S(4,3) = 2; % n=3,m=2
        spec = shDegreeRMS(C, S, 'R', 1);
        testCase.verifyEqual(spec.degVariance(3), 1, 'AbsTol', 1e-12);
        testCase.verifyEqual(spec.degVariance(4), 4, 'AbsTol', 1e-12);
        testCase.verifyEqual(spec.degVariance([1 2 5 6]), zeros(4,1), 'AbsTol', 1e-12);
    end

    function testSynthesisConstantTerm(testCase)
        % C(0,0) alone should produce a spatially constant geoid of value R*C00
        nmax = 3;
        C = zeros(nmax+1); S = zeros(nmax+1);
        C(1,1) = 2.5e-3;
        R = 6378136.3; GM = 3.986004415e14;
        [grid,~,~] = shSynthesis(C, S, GM, R, -80:20:80, 0:60:300, 'quantity', 'geoid');
        expected = R * C(1,1);
        testCase.verifyEqual(grid, expected*ones(size(grid)), 'RelTol', 1e-10);
    end

    function testSynthesisSectoralPattern(testCase)
        % pure C(2,2) term must vanish at the poles (Pbar_22(sin(+-90deg))=0)
        nmax = 2;
        C = zeros(nmax+1); S = zeros(nmax+1);
        C(3,3) = 1;
        R = 6378136.3; GM = 3.986004415e14;
        [grid,~,~] = shSynthesis(C, S, GM, R, [-90 90], [0 90], 'quantity', 'geoid');
        testCase.verifyEqual(grid, zeros(size(grid)), 'AbsTol', 1e-9);
    end

    function testSynthesisPrecomputedPMatchesFresh(testCase)
        nmax = 10;
        C = 1e-8*randn(nmax+1); S = 1e-8*randn(nmax+1);
        C = tril(C); S = tril(S);
        R = 6378136.3; GM = 3.986004415e14;
        lat = -80:10:80; lon = 0:20:340;
        [g1,~,~,P] = shSynthesis(C, S, GM, R, lat, lon, 'quantity', 'geoid');
        g2 = shSynthesis(C, S, GM, R, lat, lon, 'quantity', 'geoid', 'P', P);
        testCase.verifyEqual(g1, g2, 'AbsTol', 1e-15);
    end

    function testGaussianWeightsMonotonicAndUnitAtZero(testCase)
        Wn = shGaussianWeights(60, 300);
        testCase.verifyEqual(Wn(1), 1, 'AbsTol', 1e-12);
        testCase.verifyTrue(all(diff(Wn) <= 1e-10)); % monotonically non-increasing
        Wn_narrow = shGaussianWeights(60, 100);
        Wn_wide   = shGaussianWeights(60, 500);
        % a wider averaging radius must suppress high degrees more strongly
        testCase.verifyTrue(Wn_wide(end) < Wn_narrow(end));
    end

    function testGaussianFilterAppliesWeights(testCase)
        nmax = 5;
        C = ones(nmax+1); S = ones(nmax+1);
        [Cf, Sf, Wn] = shGaussianFilter(C, S, 300);
        for n = 0:nmax
            testCase.verifyEqual(Cf(n+1,:), Wn(n+1)*C(n+1,:), 'AbsTol', 1e-14);
            testCase.verifyEqual(Sf(n+1,:), Wn(n+1)*S(n+1,:), 'AbsTol', 1e-14);
        end
    end

    function testDestripeRemovesSmoothTrend(testCase)
        % a purely smooth (quadratic-in-n) coefficient trend at one order,
        % split by parity, should be almost entirely removed
        nmax = 40; m = 10;
        C = zeros(nmax+1); S = zeros(nmax+1);
        n = (m:nmax)';
        C(n+1, m+1) = 0.01*(n-20).^2 * 1e-9; % smooth trend, no real high-freq content
        [Cf, ~] = shDestripe(C, S, 'minOrder', m, 'polyOrder', 2);
        testCase.verifyLessThan(max(abs(Cf(n+1,m+1))), 1e-13);
    end

    function testDestripePreservesLowOrderUntouched(testCase)
        nmax = 20;
        C = 1e-9*randn(nmax+1); S = 1e-9*randn(nmax+1);
        C = tril(C); S = tril(S);
        [Cf, Sf] = shDestripe(C, S, 'minOrder', 8);
        testCase.verifyEqual(Cf(:,1:8), C(:,1:8), 'AbsTol', 1e-20);
        testCase.verifyEqual(Sf(:,1:8), S(:,1:8), 'AbsTol', 1e-20);
    end

    function testReadGFCStatic(testCase)
        thisDir = fileparts(mfilename('fullpath'));
        model = shReadGFC(fullfile(thisDir, '..', 'test_data', 'test_static.gfc'));
        testCase.verifyEqual(model.nmax, 3);
        testCase.verifyEqual(model.C(1,1), 1.0, 'AbsTol', 1e-15);
        testCase.verifyEqual(model.C(3,1), -4.84165143790815e-04, 'AbsTol', 1e-18);
        testCase.verifyEqual(model.GM, 3.986004415e14, 'RelTol', 1e-12);
        testCase.verifyEqual(model.R, 6378136.30, 'RelTol', 1e-12);
        testCase.verifyEmpty(model.variableTerms);
    end

    function testReadGFCGzip(testCase)
        thisDir = fileparts(mfilename('fullpath'));
        m1 = shReadGFC(fullfile(thisDir, '..', 'test_data', 'test_static.gfc'));
        m2 = shReadGFC(fullfile(thisDir, '..', 'test_data', 'test_static.gfc.gz'));
        testCase.verifyEqual(m1.C, m2.C, 'AbsTol', 1e-20);
        testCase.verifyEqual(m1.S, m2.S, 'AbsTol', 1e-20);
    end

    function testReadGFCTVariableTerms(testCase)
        thisDir = fileparts(mfilename('fullpath'));
        model = shReadGFC(fullfile(thisDir, '..', 'test_data', 'test_variable.gfct'));
        testCase.verifyEqual(numel(model.variableTerms), 3);
        testCase.verifyEqual(model.variableTerms(1).type, 'trnd');
        testCase.verifyEqual(model.variableTerms(1).t0, 2010.0);
        testCase.verifyEqual(model.variableTerms(1).t1, 2020.0);
    end

    function testEvalGFCT(testCase)
        thisDir = fileparts(mfilename('fullpath'));
        model = shReadGFC(fullfile(thisDir, '..', 'test_data', 'test_variable.gfct'));
        [Ct, ~] = shEvalGFCT(model, 2015.0);
        expected = -4.84165098790815e-04; % cross-validated in validate_destripe-style Python check
        testCase.verifyEqual(Ct(3,1), expected, 'RelTol', 1e-10);
    end

    function testSpectralCrossoverFindsKnownCrossing(testCase)
        spec.degree = (0:50)';
        spec.degAmplitude = 100 * exp(-0.05*spec.degree); % decaying signal
        spec.errAmplitude = 0.5 * exp(0.08*spec.degree);  % growing error
        [nCrossover, degInterp] = shSpectralCrossover(spec);
        testCase.verifyGreaterThan(nCrossover, 0);
        testCase.verifyLessThan(abs(degInterp - nCrossover), 1.5);
    end

    function testSpectralCrossoverNoCrossing(testCase)
        spec.degree = (0:20)';
        spec.degAmplitude = 100*ones(21,1);
        spec.errAmplitude = 1e-6*ones(21,1);
        [nCrossover, degInterp] = shSpectralCrossover(spec);
        testCase.verifyTrue(isnan(nCrossover));
        testCase.verifyTrue(isnan(degInterp));
    end

    function testDegreeRMSErrorSpectrum(testCase)
        nmax = 4;
        C = zeros(nmax+1); S = zeros(nmax+1);
        sigmaC = zeros(nmax+1); sigmaS = zeros(nmax+1);
        sigmaC(4,2) = 3; sigmaS(4,2) = 4; % n=3,m=1 -> sqrt(9+16)=5
        spec = shDegreeRMS(C, S, 'R', 2, 'sigmaC', sigmaC, 'sigmaS', sigmaS);
        testCase.verifyEqual(spec.errRMS(4), 5, 'AbsTol', 1e-12);
        testCase.verifyEqual(spec.errAmplitude(4), 10, 'AbsTol', 1e-12); % R=2
        testCase.verifyEqual(spec.errRMS([1 2 3 5]), zeros(4,1), 'AbsTol', 1e-12);
    end

    function testSynthesisKernelRatiosMatchAnalyticFormula(testCase)
        % pure single-degree zonal term: ratio of quantities at any grid
        % point must equal the ratio of their analytic degree kernels
        nmax = 4; n = 3;
        C = zeros(nmax+1); S = zeros(nmax+1);
        C(n+1,1) = 1e-8;
        R = 6378136.3; GM = 3.986004415e14;
        lat = -60:30:60; lon = 0:90:270;
        gGeoid = shSynthesis(C, S, GM, R, lat, lon, 'quantity', 'geoid');
        gPot   = shSynthesis(C, S, GM, R, lat, lon, 'quantity', 'potential');
        gAnom  = shSynthesis(C, S, GM, R, lat, lon, 'quantity', 'gravity_anomaly');
        gDist  = shSynthesis(C, S, GM, R, lat, lon, 'quantity', 'gravity_disturbance');

        % v2.1 patch (documented): the scaled Legendre engine returns EXACT
        % zeros at symmetry points (Pbar_30(sin 0) = 0 in floating point),
        % so the equator samples of this zonal term are 0/0 - undefined,
        % not wrong. The v1 engine's colatitude rounding fuzz (~6e-17)
        % masked this. Ratios are checked where the field is nonzero, and
        % the zero rows are verified to be consistently zero everywhere.
        nz = abs(gGeoid(:)) > 0;
        testCase.verifyTrue(any(nz));
        testCase.verifyEqual(gPot(~nz),  zeros(nnz(~nz), 1), 'AbsTol', 0);
        testCase.verifyEqual(gAnom(~nz), zeros(nnz(~nz), 1), 'AbsTol', 0);
        testCase.verifyEqual(gDist(~nz), zeros(nnz(~nz), 1), 'AbsTol', 0);

        ratioPotGeoid = gPot(nz) ./ gGeoid(nz);   % expect (GM/R)/R = GM/R^2
        ratioAnomGeoid = gAnom(nz) ./ gGeoid(nz); % expect (GM/R^2)*(n-1)/R
        ratioDistGeoid = gDist(nz) ./ gGeoid(nz); % expect (GM/R^2)*(n+1)/R

        testCase.verifyEqual(ratioPotGeoid(:), (GM/R)/R*ones(nnz(nz),1), 'RelTol', 1e-10);
        testCase.verifyEqual(ratioAnomGeoid(:), (GM/R^2)*(n-1)/R*ones(nnz(nz),1), 'RelTol', 1e-10);
        testCase.verifyEqual(ratioDistGeoid(:), (GM/R^2)*(n+1)/R*ones(nnz(nz),1), 'RelTol', 1e-10);
    end

    function testSynthesisEWHRequiresLoveNumbers(testCase)
        nmax = 2;
        C = zeros(nmax+1); S = zeros(nmax+1); C(3,1) = 1e-9;
        R = 6378136.3; GM = 3.986004415e14;
        testCase.verifyError(@() shSynthesis(C, S, GM, R, 0, 0, 'quantity', 'ewh'), ...
            'shSynthesis:missingLoveNumbers');
    end

    function testSynthesisNminTruncation(testCase)
        nmax = 3;
        C = zeros(nmax+1); S = zeros(nmax+1);
        C(1,1) = 5; C(3,1) = 1e-8; % degree 0 and degree 2
        R = 6378136.3; GM = 3.986004415e14;
        gFull = shSynthesis(C, S, GM, R, 0, 0, 'quantity', 'geoid', 'nmin', 0);
        gNoDeg0 = shSynthesis(C, S, GM, R, 0, 0, 'quantity', 'geoid', 'nmin', 1);
        testCase.verifyEqual(gFull - gNoDeg0, R*C(1,1), 'AbsTol', 1e-9);
    end

    function testSynthesisBadInputValidation(testCase)
        C = zeros(3,3); S = zeros(3,3);
        R = 6378136.3; GM = 3.986004415e14;
        testCase.verifyError(@() shSynthesis(C, S(1:2,1:2), GM, R, 0, 0), 'shSynthesis:badInput');
        testCase.verifyError(@() shSynthesis(C, S, GM, R, 0, 0, 'nmax', 10), 'shSynthesis:badInput');
    end

    function testDestripeMovingWindowReducesVariance(testCase)
        nmax = 40; m = 10;
        C = zeros(nmax+1); S = zeros(nmax+1);
        n = (m:nmax)';
        trend = 0.01*(n-20).^2 * 1e-9;
        noise = 1e-11*randn(size(n));
        C(n+1, m+1) = trend + noise;
        [Cf, ~] = shDestripe(C, S, 'minOrder', m, 'polyOrder', 3, 'windowLength', 7);
        varBefore = var(C(n+1,m+1));
        varAfter = var(Cf(n+1,m+1));
        testCase.verifyLessThan(varAfter, varBefore);
    end

    function testReadGFCMissingFileErrors(testCase)
        testCase.verifyError(@() shReadGFC('this_file_does_not_exist_12345.gfc'), ...
            'shReadGFC:fileNotFound');
    end

    function testPlotSHSpectrumRuns(testCase)
        fig = figure('Visible', 'off');
        cleanupObj = onCleanup(@() close(fig)); %#ok<NASGU>
        spec.degree = (0:10)';
        spec.degAmplitude = linspace(1, 0.1, 11)';
        h = plotSHSpectrum(spec, 'ax', axes('Parent', fig));
        testCase.verifyEqual(numel(h), 1);
        testCase.verifyEqual(get(h,'YData'), (spec.degAmplitude(2:end)*1000)', 'AbsTol', 1e-9);
    end

    function testPlotSHCoeffTriangleRuns(testCase)
        fig = figure('Visible', 'off');
        cleanupObj = onCleanup(@() close(fig)); %#ok<NASGU>
        nmax = 5;
        C = magic(nmax+1)*1e-9; S = magic(nmax+1)'*1e-9;
        ax = plotSHCoeffTriangle(C, S, 'nmax', nmax, 'ax', axes('Parent', fig));
        testCase.verifyClass(ax, 'matlab.graphics.axis.Axes');
        % v2.3.1 orientation contract: apex (n=0) on top, square cells
        testCase.verifyEqual(ax.YDir, 'reverse');
        testCase.verifyEqual(daspect(ax), [1 1 1]);
    end

end
end
