% shAnalysis v3.0.0 - Spherical harmonic (Stokes coefficient) analysis for
% Version 3.0.0 (R2026a-compatible) 10-Aug-2026
% GRACE/GRACE-FO, GOCE, and static gravity field models: class-based,
% with time-series support, GAX background handling, climatology, and
% the tvANS time-variable anisotropic Wiener filter.
%
% Classes (single point of access)
%   shCoefficients  - One coefficient set: I/O, arithmetic, TN-14, destripe,
%                     Gaussian smoothing, spectra, synthesis, gfct evaluation.
%   shSeries        - Time series of coefficient sets: mean field, climatology,
%                     GAX restore, per-epoch classic chain, tvANS filter,
%                     basin averages with exact deconvolution.
%   shClimatology   - Bias/trend/annual/semi-annual (+extra periods) model:
%                     evaluate at epochs, extract components with 1-sigma
%                     significance, amplitude/phase maps.
%
% New in v2.1
%   shCoefficients.analysis - SH ANALYSIS: gridded or scattered data ->
%                     Stokes coefficients (rings: exact fast per-order
%                     solver; scattered: least squares + Kaula).
%   addDegree1      - TN-13 geocenter insertion (with shLowLevel.readTN13),
%                     completing the GSM + TN-14 + TN-13 chain.
%   filter(..., Blocks="auto") - block-diagonal tvANS (identical results,
%                     Lmax ~ 120 tractable) + posterior sigma stacks
%                     (tsF.sigmaCs/Ss) and basin-average sigmas.
%   climatology(Periods=[161/365.25, 3.66, 7.48]) - S2/K2/K1 tidal alias
%                     terms; per-coefficient significance on components.
%   synthesis(..., Method="fft") - FFT-along-longitude + latitude parity
%                     trick; Legendre recursion stable to degree 2190.
%   dropNaN/select  - GRACE <-> GRACE-FO gap handling.
%   ICGEM 2.0       - piecewise gfct files; shLowLevel.readSINEX covariances/NEQs.
%
% Quick start
%   ts  = shSeries.read("GSM-2_*.gfc");
%   ts  = ts.restore(shSeries.read("GAD-2_*.gfc"));   % add background
%   dts = ts - ts.mean;                               % anomalies
%   [clim, resid] = dts.climatology;                  % seasonal model
%   tsA = dts.destripe(minOrder=6).gaussian(300);     % classic chain
%   [tsF, op] = dts.filter("tvANS");                  % optimal chain
%   avg = tsF.basinAverage(B, Deconvolve=true, Op=op);
%
% Internals (package +shLowLevel, stable but not the primary API)
%   shLowLevel.legendreALF, shLowLevel.legendreCached, shLowLevel.shSynthesis, shLowLevel.shDestripe,
%   shLowLevel.shGaussianWeights, shLowLevel.shDegreeRMS, shLowLevel.readTN14, shLowLevel.writeGFC,
%   shLowLevel.synthesisMatrix, shLowLevel.tvANSFilter, shLowLevel.basinDeconvolve,
%   shLowLevel.resolutionMap, shLowLevel.shIndex, shLowLevel.shAnalysisGrid, shLowLevel.readTN13,
%   shLowLevel.readLoveNumbers (v2.5), shLowLevel.fetchTN (v2.5), setup_shAnalysis (v2.5),
%   shLowLevel.readSINEX, shLowLevel.icgemDate2Year, shLowLevel.kernelFactors, ...
%
% Backward compatibility (folder compat/, add to path if needed)
%   shReadGFC, shSynthesis, shDestripe, shGaussianFilter, shGaussianWeights,
%   shDegreeRMS, shSpectralCrossover, shEvalGFCT, legendreALF,
%   plotSHSpectrum, plotSHCoeffTriangle - v1 signatures, delegating to +shLowLevel.
%
% Tests and documentation
%   runAllTests     - full validation suite (correctness, contract,
%                     robustness, legacy v1, performance).
%   doc shAnalysis  - HTML overview (html/ subfolder on the path).
%
% New in v2.2 (all Python-validated before MATLAB implementation)
%   shLowLevel.basinKernel          - basin kernels from polygons/masks/functions
%                              (buffering + spectral taper)
%   shLowLevel.slepianBasis         - Slepian localization (regional analysis)
%   shLowLevel.readDDK, shLowLevel.applyDDK, g.applyDDK, ts.applyDDK
%                            - DDK anisotropic decorrelation filters
%                              (native binary Wbd support + Nmax
%                              truncation, v2.2.1)
%   ts.removeGIA, clim.removeGIA
%                            - GIA correction (user-supplied rate models)
%   ts.climatology(ARCorrect=true)
%                            - AR(1)-corrected trend/coefficient sigmas
%   shLowLevel.mcPropagate          - Monte-Carlo sigma propagation (independent
%                              or full SINEX covariance)
%   shLowLevel.readMascon           - JPL/CSR mascon netCDF reader
%   shLowLevel.geodetic2geocentric, shLowLevel.geocentric2geodetic, and
%   LatType="geodetic" in synthesis/analysis entry points
%   shLowLevel.shSynthesis(..., 'MaxMemGB', 4) - latitude-banded streaming makes
%                              nmax 2190 real on ordinary RAM
%   shLowLevel.tvANSFilter(..., VCEBands=[0 16 33 61]) - per-order-band VCE
%   shLowLevel.readSINEX(..., Only="estimate") - streaming estimate block of
%                              multi-100-MB SINEX in seconds
%
% New in v2.3 (kernels and deformation; Python-validated)
%   New synthesis/analysis quantities in shLowLevel.kernelFactors:
%     'surface_density'     [kg/m^2]  R rho_ave/3 (2n+1)/(1+kn)
%     'bottom_pressure'     [Pa]      g0 * surface density (OBP)
%     'deformation_up'      [m]       R hn/(1+kn) - elastic vertical
%     'gravity_gradient_rr' [1/s^2]   (GM/R^3)(n+1)(n+2)  (T_rr)
%   'Height' option: upward continuation (R/r)^p to satellite altitude
%   shLowLevel.shSynthesisDeformation / g.deformation - up/north/east elastic
%     load deformation (GNSS comparison); grid and station-list modes;
%     exact dPbar/dphi via shLowLevel.legendreALFDeriv (frozen validated
%     identity, ~5e-10 against numerical gradients)
%
% New in v2.4 (visualization + computation extensions; Python-validated)
%   Visualization: shLowLevel.plotSHMap / g.map (divergent map, coastlines,
%     Hammer projection), shLowLevel.plotBasinSeries (sigma band, gaps, trend),
%     triangle difference mode (RefC/RefS), shLowLevel.plotCovariance (block
%     boundaries), plotSHSpectrum Kaula + crossover overlays,
%     shLowLevel.writeAnimation (MP4, fixed robust scale)
%   shLowLevel.errorMap - analytic sigma maps from full covariances (chol,
%     cross-checked vs mcPropagate)
%   shLowLevel.eofAnalysis - EOF/PCA modes as shCoefficients + unit PCs
%   shLowLevel.basinScaling - forward-modelled gain factors (struct | matrix |
%     handle operators)
%   shSeries.trendBreaks - hinge breakpoints with F-test (betainc-based
%     p-values, no toolbox; null calibrated at 5.1%/5%)
%   shLowLevel.shFanFilter / g.fan / ts.fan - Han degree x order Gaussian
%   shLowLevel.combineCenters - per-(center,month) VCE combination with
%     partial redundancies, robust block-median option, inter-center
%     correlation diagnostic, posterior sigmas
%   shLowLevel.shSynthesisGradientTensor - full NEU tensor at altitude
%     (Laplace trace 7e-16 self-check); legendreALFDeriv 2nd derivative
%   shLowLevel.seaLevelFingerprint - elastic sea-level equation (mass exactly
%     conserved; near-field fall / far-field excess validated)
%   shLowLevel.writeGrid - CF-style netCDF export (roundtrip-tested)
%   Reference systems: shLowLevel.normalFieldCS (even zonals COMPUTED from
%     WGS84/GRS80 defining constants, matches NIMA TR8350.2 to all
%     published digits), g.subtractNormalField (auto-rescaled),
%     shLowLevel.rescaleGMR / g.toReference (exact (GM,R) conversion; field
%     invariance verified); arithmetic guards GM/R mismatches
%
%   demo_shAnalysis - selectable case registry (D01..D16, "list"/"core"/
%     "all"/subset, Visible= for headless runs); rendered gallery with
%     explanations in the workflow guide PDF
%   shLowLevel.fetchITSG / shSeries.fromFolder - on-demand ITSG download
%     (websave, release routing, gap-aware) + folder loading; demos
%     D01-D06 use real data when present, D12-D14 stay synthetic by
%     design (known-truth recovery)
%   shLowLevel.pctile - percentile without the Statistics Toolbox (an audit
%     found prctile had slipped into two plot functions; fixed)
%
% New in v2.5 (uncertainty exactness, Love numbers, real provider files)
%   Constrained tvANS posterior sigma EXACT (diag of (W-I)S(W-I)'+sWNW';
%   the old formula underestimated, not "upper bound"). Basin sigma now
%   carries the deterministic-fit parameter uncertainty (leverage x
%   resVar; MC 1.18 -> 1.00) and the exact constrained noise covariance.
%   ARCorrect applies the Kendall r1 bias correction (MC 1.065 -> 1.028).
%   shLowLevel.readLoveNumbers: layouts, sparse pchip interpolation, degree-1
%   frame conversion (Blewitt 2003; CE->CF/CM, CF<->CM). Real CSR/JPL
%   TN-13 fixtures (256 months each, cross-provider corr 0.995); all 8
%   DDK Wbd files parsed + gain-ordering validated (discovery test).
%   Robust cache-speedup performance test. setup_shAnalysis: one-call
%   path setup (temporary/permanent) + cumulative download levels
%   (core TN files via new shLowLevel.fetchTN with verify-by-parse, DDK
%   filters, ITSG starter months); DryRun plan mode, offline-safe.
%   Real GSFC TN-14 ships as fixture (258 windows, pinned test).
%
% New in v2.4.2 (license, attribution, fixes)
%   MIT LICENSE file (Matthias Weigelt); every help text and doc page
%   closes with the attribution statement. plotSHCoeffTriangle diff mode:
%   sine sectorals S_nn no longer blanked (S wing was one column narrow).
%   shLowLevel.listICGEM(Type="temporal"): full ~70-series catalogue (was 3
%   release-note rows) + Series= option listing single downloadable files.
%
% New in v2.4.1 (plots, spectra, data management)
%   Triangle plots without the center line; spectrum plots on a LINEAR
%     degree axis (log y kept)
%   plotSHSpectrum Quantity= amplitude|rms|variance|cumamplitude|
%     cumrms|cumvariance; order-domain plots via shLowLevel.shOrderRMS (the
%     striping axis; marginals tested against shDegreeRMS); cumulative
%     variance/RMS fields added to shDegreeRMS
%   shLowLevel.version - toolbox metadata (name, version, date, root) from
%   Contents.m; shLowLevel.fetch* now take Update=true for a safe, parse-verified
%   refresh of existing files (TN-13/TN-14 grow monthly upstream).
%   v3.0.0 - BREAKING: the function namespace is renamed +shx ->
%   +shLowLevel (all calls shLowLevel.*; error identifiers
%   'shLowLevel:*'); compat/ and the v1 legacy suite are gone for good.
%   New: shLowLevel.listITSG (server catalogue, ends silent release
%   mixing), fetchITSG Release=/Catalog=/"all", numbered ICGEM catalogue
%   with index/"all" selection, setup_shAnalysis DataFolder= and
%   FetchITSG="all", shLowLevel.fetchLoveNumbers (GROOPS set),
%   shLowLevel.synthesisPoints (pointwise quantities at lat/lon/r).
%   New in v2.7.0: shLowLevel.poleTideConvert (IERS2010 <-> IERS2018 mean-pole
%   conventions for C21/S21, solid + ocean, Python-validated - mixing
%   conventions biases C21/S21 trends at the mm-EWH level); Proxy= in
%   all fetchers (per-call proxy via matlab.net.http); provenance JSON
%   sidecars on writeGFC/writeGrid/writeAnimation (Sidecar=false to
%   disable); writeGrid netCDF is CF-1.8 complete.
%   New in v2.6.0 - the comparison suite: shLowLevel.compareSolutions and
%   shLowLevel.compareSeries (also g.compare / ts.compare) aggregate the
%   standard metric set over the primitives shLowLevel.diffSpectrum (difference
%   degree amplitude, degree correlation, agreement crossover),
%   shLowLevel.spatialStats (area-weighted bias/RMSD/pattern stats, Taylor-
%   ready), shLowLevel.nashSutcliffe, shLowLevel.effectiveCorr (AR(1)-corrected
%   significance), shLowLevel.threeCorneredHat (per-solution noise, N >= 3)
%   and shLowLevel.taylorDiagram. All numerics Python-validated.
%   Also in v2.5.1: applyDDK truncates Lmax-120 filter blocks to the
%   field's nmax (standard n60/n96 use; previously errored), and
%   shLowLevel.pctile accepts vector percentiles.
%   shLowLevel.dataFolder - persistent user data folder (getpref/setpref);
%     all fetchers store beneath it
%   shLowLevel.fetchDDK + shLowLevel.readDDK("DDK<n>") - all eight released DDK
%     filters on demand (mapping in shLowLevel.ddkNames, repository-verified);
%     DDK3 still ships in tests/test_data
%   shLowLevel.fetchITSG Product="daily" - ITSG daily Kalman solutions (one
%     .gfc per day, n40, formal errors; -> dataFolder/itsg_daily; daily
%     filename epochs parse automatically; real fixture file shipped)
%   shLowLevel.listICGEM / shLowLevel.fetchICGEM - ICGEM static-model catalogue
%     (fixture-tested parser) and .gfc download by name; temporal
%     section returns series roots with an honest JS/ZIP note
%
% Claude (Fable 5), 2026-08-07 (v2.1 through v2.4.1 same day).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
