% shAnalysis - Spherical harmonic analysis toolbox
% Version 3.5.1 (R2026a-compatible) 11-Aug-2026
%
% The line above is what ver('shAnalysis') reports as the product name:
% keep it a SHORT name, not a sentence and not a version string (pinned
% by testVersionMetadataIsConsistent).
%
% Spherical harmonic (Stokes coefficient) analysis for GRACE/GRACE-FO,
% GOCE and static gravity field models: class-based, with time-series
% support, GAX background handling, climatology, and the tvANS
% time-variable anisotropic Wiener filter.
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
% Backward compatibility
%   NONE. v3.0.0 removed compat/ and the v1 legacy suite for good. The v1
%   function names live on inside the package: call shLowLevel.shReadGFC,
%   shLowLevel.shSynthesis, shLowLevel.shDestripe etc., or use the classes.
%
% Tests and documentation
%   runAllTests     - full validation suite (correctness, contract,
%                     robustness, performance).
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
%   New in v3.1.0: shLowLevel.standardChain - the canonical pipeline
%   (read -> TN-14 -> degree-1 -> optional GIA -> filter) as a single,
%   correctly ordered entry point with a provenance report; and
%   shLowLevel.designFilter - DDK-class anisotropic filters
%   W = (N + a*inv(S))^-1 N built from YOUR sigmas/covariances, in the
%   readDDK block format so they drop into applyDDK unchanged.
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
%     section returns series roots (superseded in v3.1.1 - the series are
%     downloadable now, see below)
%
% New in v3.5.1 (the open-ocean noise metric)
%   shLowLevel.oceanRMS - area-weighted RMS of a field over the open
%     ocean (> 1000 km from any coast), the standard GRACE noise metric
%     and the natural value for leakageCorrect's NoiseLevel. The ocean
%     mask is USER-SUPPLIED, like Love numbers: base MATLAB's
%     'coastlines' is absent from some installations and a wrong mask
%     silently changes the number. Area weighting is on by default -
%     unweighted over-counts the polar rows, which white noise hides.
%     Python-validated in tools/dev/validate_oceanrms.py.
%
% New in v3.5.0 (stopping is regularisation)
%   leakageCorrect solves an ill-posed problem and SEMICONVERGES: the
%     error against the truth falls, then RISES, while the residual
%     keeps shrinking. Iterating until nothing changes is therefore the
%     wrong rule - on the reference problem the final solution is 361x
%     worse than the best, and a step tolerance of 1e-3 stops 89x past
%     the optimum (1e-4 never triggers at all).
%   NoiseLevel= enables the discrepancy principle (Morozov): stop once
%     the residual reaches the noise level of the data. Lands within 2%
%     of the optimum. Tau= (1.2) is the safety factor. Without it the
%     run is unregularised, info.stoppedBy says so, and a warning is
%     raised. Python-validated in tools/dev/validate_stopping.py.
%
% New in v3.4.1 (the GravIS comparison closed to 1.6 %)
%   leakageCorrect: the Filter= you declare must be the filter the
%     FORWARD MODEL sees, not the one your input file carries. GravIS
%     ice basin averages start from UNFILTERED coefficients and still
%     apply a Wiener filter (~Gaussian, 4 deg latitude half-width)
%     inside the inversion; declaring "none" instead of "gauss445"
%     moved the Greenland trend from -234.9 to -227.4 Gt/yr against a
%     published -231.1. Documented in the guide and on the topic page.
%   designFilter: VDK (Horvath et al. 2018) is DDK rebuilt from each
%     MONTH's error covariance - i.e. a per-epoch usage of designFilter
%     with Noise=. What is missing is the data, not the method: monthly
%     full covariances are not distributed with Level-2 products.
%
% New in v3.4.0 (tidal alias removal; the GravIS validation written up)
%   shSeries.removeAlias - fit a tidal-alias harmonic (S2, 161 d by
%     default) jointly with bias/trend/annual/semi-annual and subtract
%     ONLY that harmonic, as GravIS does for Level-2B. A fixed 100 deg
%     phase offset applies across the GRACE/GRACE-FO boundary (Landerer
%     et al. 2020); ignoring it mis-estimates the amplitude badly.
%   Guide chapter "Validation against GravIS": the whole comparison with
%     the executed code, the results table and an honest account of the
%     remaining 4%.
%
% New in v3.3.1 (leakage masks, validated against GravIS)
%   leakageCorrect: THE MASK MUST COVER EVERY REGION THAT CAN HOLD MASS,
%     not only the one being measured - a target-only mask forces
%     neighbouring signal into the target and biases the result high.
%     Quantified against the published GravIS Greenland series: a
%     Greenland-only mask overshoots by 12%, a union mask including the
%     Canadian Arctic, Iceland and Svalbard by 5%. Documentation only;
%     no API change was needed, Mask= already accepts any union.
%
% New in v3.3.0 (GRAVIS Level-2B reader)
%   shLowLevel.readSHM reads the GRAVIS/GRACE SHM format - YAML header
%     plus GRCOF2 (field) or GRDOTA (rate, 1/yr) records, gzip
%     transparent, GM/R from the header. GravIS Level-2B products are a
%     major public data source the toolbox previously could not ingest
%     at all; a GRDOTA result feeds standardChain(GIA=) directly.
%   Fixtures for both record types plus the GravIS Greenland drainage
%     basins ship in tests/test_data (CC-BY-4.0, see NOTICE).
%
% New in v3.2.2 (fixes found by running on a full 24-year series)
%   standardChain no longer stops when the TN-13/TN-14 tables trail the
%     solutions. They always do - a provider publishes the correction
%     weeks after the monthly field - so the newest months of a fresh
%     series were routinely uncovered and the chain simply errored.
%     Uncovered epochs are dropped and recorded in the report;
%     OnMissing="error" restores the old behaviour, Tolerance= is now
%     settable.
%   shSeries.read/fromFolder no longer mistake writeGFC's
%     "<file>.gfc.provenance.json" sidecars for solutions: a folder
%     written BY the toolbox could not be read back BY the toolbox.
%   shCoefficients.write gained Sidecar=, which only the low-level
%     writeGFC had.
%
% New in v3.2.1 (scientific regression suite - roadmap item 9)
%   tests/testScience.m checks the toolbox against values published
%   OUTSIDE it - the WGS84/GRS80 defining constants, the Gegout97 (PREM)
%   load Love numbers, the closed-form EWH kernel of Wahr et al. (1998),
%   the GRACE-minus-SLR C20 discrepancy, geocenter amplitudes from the
%   providers' technical notes, and the Gaussian half-weight relation.
%   The other suites verify that the code does what the code intends and
%   stay green if a formula is CONSISTENTLY wrong; these do not.
%   A trend regression needs a monthly series, which is too large to
%   ship: set SHX_SERIES_FOLDER to opt in, otherwise that one test is
%   filtered (it deliberately does not touch the persistent data folder).
%
% New in v3.2.0 (leakage correction - roadmap item 8)
%   shLowLevel.leakageCorrect - iterative forward modelling: recover the
%     mass field that, pushed through the SAME chain the data saw,
%     reproduces the observation. Mask= confines the solution to where
%     mass can exist, which is the well-conditioned variant and removes
%     leakage instead of redistributing it. Convergence is judged on the
%     change of the SOLUTION, since a masked problem is inconsistent and
%     its residual floors above any tolerance.
%   shLowLevel.gridScaling - per-pixel scaling factors k from a model
%     series pushed through the same chain. Amplitude-invariant (k is a
%     property of the model's PATTERN) and NaN where the model carries
%     no signal, so a model that does not reach your region is visible
%     rather than silently multiplying data by noise.
%   Both Python-validated first (tools/dev/validate_leakage.py).
%
% New in v3.1.4 (robust safe swaps)
%   shLowLevel.safeMove - the fetchers' final .part -> target swap, with
%   verification and exponential-backoff retries. Windows antivirus and
%   cloud-sync clients hold a freshly written file for a moment, which
%   made a plain movefile report a download failure for a file that had
%   downloaded perfectly. All six fetchers use it.
%
% New in v3.1.3 (documentation sync gate)
%   tools/doc_sync_audit.py gates the three documentation sources against
%   each other in CI: doc snippets against the parsed contracts (option
%   value TYPES included), API-reference freshness, helptoc reachability,
%   version consistency, help placeholders and defaults, narrative
%   coverage. tools/dev/README.md documents the toolchain.
%
% New in v3.1.2 (documentation coverage)
%   Narrative help pages for the v3.x features that previously existed
%   only in the generated API reference: shLowLevel.standardChain,
%   shLowLevel.designFilter and the ICGEM time-series workflow get their
%   own pages; poleTideConvert, synthesisPoints, listITSG and
%   fetchLoveNumbers get sections on the matching topic pages. Workflow
%   guide Edition 5. Every name-value option now carries a real
%   description in the in-file help (58 pointed at the code instead).
%
% New in v3.1.1 (ICGEM time series and real-world file layouts)
%   shLowLevel.fetchICGEM Type="temporal" - a temporal-catalogue row now
%     downloads the WHOLE monthly series to
%     <dataFolder>/icgem/series/<group>_<center>_<series>/.
%     Mode="auto" (default) takes the server's whole-series ZIP in ONE
%     request and falls back to resumable per-file fetching; "archive"
%     and "files" force either. Files= filters, FileList= injects a
%     catalogue (offline mirrors, subsets), Pause=/Retries= carry the
%     rate-limit discipline, info.mode reports what actually happened.
%     The result folder feeds shSeries.fromFolder and standardChain.
%   shLowLevel.shReadGFC - group-wise bulk parsing for ALL files, static
%     and variable (GRGS mean fields: 73.6 MB / 674k terms in ~7 s,
%     previously a multi-minute stall; EIGEN-6C4 177.7 MB / n2190 in
%     ~5 s). FORTRAN D-exponents are read correctly (str2double returns
%     NaN for them - such files were silently corrupted above the
%     degree at which providers switch notation). The ICGEM 2.0
%     acos/asin column order is ... t0 t1 period (period LAST, verified
%     against real CNES/GRGS files); ragged groups (EIGEN-5S/5C) are
%     subgrouped by width with an n/m sanity net.
%   ts.filter("tvANS", ...) forwards Shrinkage/VCEMinDegree/VCEBands to
%     shLowLevel.tvANSFilter - the class method is the full single point
%     of access again.
%
% Claude (Fable 5), 2026-08-07 (v2.1 through v2.4.1 same day);
%   v3.1.1 documentation sync 2026-08-11 by Claude (Opus 5).
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).
