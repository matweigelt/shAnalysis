% shAnalysis - Spherical harmonic analysis toolbox
% Version 3.25.0 (R2026a-compatible) 18-Aug-2026
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
% New in v3.25.0 (relocatable test fixtures)
%   - tests/shxTestDataDir resolves the fixture folder: env
%     SHX_TESTDATA_FOLDER first, preference TestDataFolder second,
%     tests/test_data as the in-repo fallback - fixtures can live on a
%     data drive while CI keeps the repo copy unchanged. All direct
%     test_data paths in the four suites route through it; new
%     setup_shAnalysis option TestDataFolder follows the
%     SeriesFolder/DDKFolder env+setpref convention.
%
% New in v3.24.1 (gate hardening: mlint_lite R5)
%   - tools/dev/mlint_lite.py rule R5: a name-value struct declared in
%     an arguments block (opts.Field) must appear as the LAST input of
%     the function declaration line - exactly the defect that escaped
%     to CI in the v3.24.0 fixed-lag PR (green gates, every call red).
%     Six-case self-test incl. that escaped file verbatim; repo-wide
%     sweep of 143 .m files clean. No toolbox code changed.
%
% New in v3.24.0 (fixed-lag smoother)
%   - rtsSmoother Lag=L: the estimate at epoch t uses observations up
%     to t+L only - the near-real-time production variant (Kurtenbach's
%     practical choice). Definition-true windowed backward recursion,
%     cost O(T*L) against O(T) for the full pass; Lag=0 returns the
%     filter, Lag >= T-1 the full smoother to rounding (unit-tested,
%     also through the matfile store); the error decays with L at the
%     process-memory rate (Python: validate_kalman_qc.py L1/L2).
%     Closes the NRT loop: daily solutions with bounded latency feed
%     hydroExtremeIndex directly (any shSeries does).
%
% New in v3.23.0 (multi-center Kalman: neqCombine inside kalmanChain)
%   - kalmanChain accepts a string ARRAY of SINEX folders (one per
%     center): files are clustered into epoch groups within Tolerance
%     and each group is combined on the normal-equation level by
%     neqCombine BEFORE the filter update - per-epoch VCE weights when
%     every file carries the STATISTICS block, or fixed NeqWeights.
%     rep.centers and rep.sigma2 (F x K per-epoch variance factors)
%     report the combination. Same-background contract across ALL
%     centers (neqCombine). Single-folder behaviour unchanged.
%
% New in v3.22.0 (Joseph-stabilized solution update)
%   - kalmanFilter solution mode uses the Joseph form
%     (I-KH) Pm (I-KH)' + K R K': identical in exact arithmetic
%     (batch-LSA and Wiener-limit tests unchanged), but immune to the
%     (I-K)Pm cancellation when R is small against the prior - the
%     strong-daily-data regime. Measured against a 50-digit reference:
%     standard 1.3e-6 relative covariance error, Joseph 1.9e-13
%     (tools/dev/validate_kalman_qc.py J1/J2). The NEQ mode was
%     already on the numerically clean information form.
%
% New in v3.21.0 (innovation-based quality control, Kvas Sec. 3.3)
%   - kalmanFilter QC=("none")|"flag"|"reject": chi-square innovation
%     test per epoch, solution form d'S^-1 d and NEQ form
%     u'(N P- N + N)^+ u with dof = rank(N) (identical when N = R^-1;
%     Python-validated). "reject" turns a failing epoch into
%     prediction-only so a blunder never enters the recursion - it
%     would otherwise be dragged through all later epochs. Off by
%     default (rank(N) costs an SVD per NEQ epoch). kalmanChain
%     passes QC/QCAlpha through; rep.qcStat and rep.nRejected.
%   - shLowLevel.chi2Quantile - base-MATLAB chi-square quantile
%     (Wilson-Hilferty; same reason pctile exists). Accuracy MEASURED
%     against scipy over the QC range and stated in the help: <= 3.5%
%     at dof <= 10, essentially exact at the operating point dof = P.
%   - Python pre-validation tools/dev/validate_kalman_qc.py (with the
%     importable core kalman_port_base.py): checks Q1-Q4.
%
% New in v3.20.0 (neqCombine: VCE combination on the NEQ level)
%   - shLowLevel.neqCombine - COST-G-style combination of K normal
%     equations BEFORE solving, with Foerstner/Koch variance components
%     (partial redundancies; the invariant sum(r) = sum(nobs) - P is
%     asserted every iteration). Fixed Weights= skip VCE; VCE requires
%     ltpl/nobs on every contribution and errors loudly otherwise.
%     Accepts raw N/b structs or readSINEX NEQ results directly. The
%     rigorous sibling of combineCenters (which works on solutions).
%   - readSINEX parses +SOLUTION/STATISTICS into snx.stats (keys/values
%     verbatim + nobs/nunk/dof/wsos/vfactor) - exactly what
%     neqCombine's VCE feeds on - and now DELIVERS snx.epoch (decimal
%     year from the estimate REF_EPOCH): the help had promised the
%     field since v2.x but it never existed, a latent kalmanChain
%     NEQ-mode crash on first real use, fixed and pinned by test.
%   - Numerics pre-validated in Python (tools/dev/validate_neqcombine.py,
%     5 checks: fixed == stacked GLS, factor recovery, redundancy
%     invariant, combined beats singles, equal-noise symmetry).
%
% New in v3.19.0 (matfile covariance store: long daily Kalman runs)
%   - kalmanFilter StoreCov="matfile": predicted/filtered covariances
%     stream through a v7.3 MAT with partial I/O - RAM stays flat at
%     one (pP)^2 working copy, so n40 daily years become feasible on
%     ordinary machines. Results identical to StoreCov="full" to the
%     last bit (unit-tested through filter AND smoother, gap included).
%     The file path is returned in filt.covFile; the filter's caller
%     deletes it.
%   - rtsSmoother reads the matfile epoch-by-epoch (read-only; file
%     survives smoothing); "diag" still refuses with an identified
%     error.
%   - kalmanChain StoreCov=("auto")|"full"|"matfile": chain-internal
%     temp file, deleted after smoothing (and in the no-smoother
%     branch - no leak); rep.memGB then reports the disk footprint.
%
% New in v3.18.0 (buildCondFun: Kvas EWH-domain covariance conditioning)
%   - shLowLevel.buildCondFun - builds the CondFun handle for estimateVAR
%     implementing Kvas (2019) Sec. 2.4: spectral covariance -> EWH grid
%     (exact Gauss-Legendre quadrature pair, F*G = I), region-block
%     masking (eq. 2.117, e.g. land/ocean) and the distance taper
%     exp(-psi/psi0) (eq. 2.120), back to the spectral domain. Both
%     weights are PSD kernels, so the conditioned covariance stays PSD
%     (Schur product theorem); a singular empirical Sigma(0) from a
%     short series becomes strictly positive definite - the Yule-Walker
%     stabilization Kvas conditions for. Love numbers user-supplied.
%   - Numerics pre-validated in Python (tools/dev/validate_condfun.py,
%     4 checks: identity at psi0=Inf, PSD chain, rcond rescue of a
%     singular covariance, cross-region zeroing).
%
% New in v3.17.0 (Kalman/VAR module: Kurtenbach/Kvas temporal smoothing)
%   - shLowLevel.estimateVAR - empirical VAR(p) process model from an SH
%     state series via Yule-Walker; Order=1 is exactly Kurtenbach (2012)
%     eqs. (3.84)-(3.85); Shrink and CondFun (Kvas-style covariance
%     conditioning) for stability; companion spectral-radius warning.
%   - shLowLevel.kalmanFilter - forward Kalman filter, solution mode
%     (l = x + v with R) and NEQ mode (information-form update straight
%     from readSINEX normal equations); gaps = prediction-only; data
%     contribution per coefficient (Kurtenbach Sec. 3.3.2).
%   - shLowLevel.rtsSmoother - RTS backward pass; with the stationary
%     initialization identical to the joint least-squares adjustment
%     over all epochs (Kvas Sec. 2.3; unit-tested to machine precision).
%   - shLowLevel.kalmanChain - single point of access: model series ->
%     VAR model -> filter -> smoother -> shSeries with formal sigmas;
%     observation climatology restored in solution mode; daily grids
%     with prediction across missing days via Epochs=.
%   - Numerics pre-validated in Python (tools/dev/validate_kalman.py,
%     8 tests including KF+RTS == batch adjustment at 1e-15).
%
% New in v3.16.2 (final documentation audit)
%   - Guide Part IV input tables now carry a DESCRIPTION column, fed
%     from the help Inputs/Options lines the extractor already parsed
%     (881 documented arguments; coverage 99.5% -> 100% after adding
%     the missing demo_shAnalysis Inputs/Options block).
%
% New in v3.16.1 (CoastBufferKm default 300; iter5 acceptance)
%   - oceanChain CoastBufferKm defaults to 300 km (author decision on
%     the live 6e table): the barystatic baseline moves +1.407 ->
%     +1.557 mm/yr, into the published ~1.6-2.2 band; pass 0 for the
%     pre-v3.16.1 numbers. Acceptance criteria rebaselined.
%
% New in v3.16.0 (uniform fetch layout; fetchITSGSINEX)
%   BREAKING: consistent on-disk layout and naming for the fetch
%   family, decided with the author.
%   - fetchSINEX is now fetchITSGSINEX (source-consistent with
%     fetchITSG/fetchITSGBackground); error IDs moved accordingly.
%   - Default targets: series/itsg/{monthly,daily,sinex,background},
%     series/<group_center_series>, series/GAX/<product> for
%     temporal products; static/ for static ICGEM models; DDK/, TN/
%     unchanged. fetchGAX's dest argument is now optional.
%   - Every fetcher warns LOUDLY (once) when it finds data in a
%     pre-v3.16 location instead of silently re-downloading; move
%     the old folders to the new layout to keep skip-if-present.
%   - run_vdk_series: renamed call sites + fixed a latent
%     shLowLevel.shxVersion() typo (first bridge run caught it;
%     correct name is shLowLevel.version()).
%%   - eofSeparate: North MULTIPLET rule (degenerate leading modes kept
%     as a group) guarded by a median-calibrated Marchenko-Pastur bulk
%     edge; oceanChain caps SeparateCirculation at 3 modes (the MP
%     edge is liberal under correlated GRACE residuals - the uncapped
%     chain kept 10). Real trigger: the grown 252-month series
%     degenerated the leading pair (gap 0.44e9 vs dl 0.53e9, nKeep 0).
%   - Machine acceptance (bridge restored): suite 225 tests green,
%     section 5 PASS (3 modes, sigMon 0.0148, trend +1.407), 6c PASS
%     (Amazon droughts 2015/16 > 2010 > 2005: -3.05/-1.43/-1.18),
%     6d PASS (daily path), 6e live-baselined: buffer 300/500 km
%     moves the barystatic trend +1.407 -> +1.557/+1.647 mm/yr INTO
%     the published band; removal(1)+300 km +1.524 (polar-cap
%     leakage correctly removed).
%
% New in v3.15.1 (guide: complete theory edition)
%   - Part I extended for readers without prior knowledge: new entry
%     chapter 0 (what GRACE measures, the product vocabulary GSM/GAX/
%     SINEX/daily/mascons, why stripes and low degrees matter) and new
%     chapters 21-26 covering everything since v3.10: VDK/VADER theory
%     and its exact relation to DDK/tvANS, normal equations and SINEX,
%     ocean-mass restoration and OBP with EOF circulation separation,
%     coastal leakage (buffer numbers and the two-sided support
%     separation, including the one-step warning), standardized
%     hydrological indices, daily Kalman NRT monitoring.
%   - New chapter 27: full reference list (19 primary sources).
%   - Four new computed theory figures (Wiener gain family, 1D leakage
%     + POCS convergence, DSI categories with the detrend effect,
%     causal storage deficit) - all from the session-validated
%     experiments, honestly labelled as computed illustrations.
%
% New in v3.15.0 (coastal-leakage controls for the ocean chain)
%   - oceanChain CoastBufferKm=: great-circle mask erosion (public
%     helper shLowLevel.erodeMask). Quantified before building
%     (Python gauss445 point source): 300 km cuts 54% of the coastal
%     leak, 500 km 82%, at < 1% ocean-area loss.
%   - oceanChain RemoveLandLeakage= (LandLeakIter=5): iterative
%     two-sided support separation (Papoulis-Gerchberg style) of the
%     outside-mask sources, subtracted BEFORE the filter. The
%     pre-validation caught that a naive one-step subtraction WORSENS
%     the leak (interior ringing of the outside reconstruction); the
%     iteration cuts an adjacent-source leak of 47% of the ocean
%     mean by 94% (frozen in CI).
%   - Defaults stay OFF: the published +1.41 mm/yr acceptance stays
%     reproducible until the machine run re-baselines it (accept 6e
%     prints the buffer/removal trend table against ~1.6-2.2).
%
% New in v3.14.0 (flood/drought indices per cell or basin)
%   - shLowLevel.hydroExtremeIndex: GRACE-DSI (Zhao et al. 2017,
%     standardized anomaly per calendar month, 11 USDM-style
%     categories), WSDI (Sinha et al. 2017), and the causal Reager &
%     Famiglietti (2009) storage deficit - with PrecipGrid= the full
%     Flood Potential Index (multi-month flood lead times, Reager et
%     al. 2014). Works on twsChain grids, basin series, or any stack.
%   - Detrend policy QUANTIFIED by the Python pre-validation: an
%     exceptional-drought month under a 0.5 cm/yr trend weakens from
%     -2.41 to -1.44 (out of class) undetrended, 45% of the late
%     decade turns spuriously wet - Detrend="linear" is the DSI/WSDI
%     default; StorageDeficit keeps the physical Reager convention.
%   - Robust sigma option (1.4826*MAD): one corrupt month inflates
%     the classical sigma sixfold, the MAD not at all (measured).
%   - Stage 2, daily: on ITSG daily Kalman solutions the index
%     switches (auto) to a day-of-year climatology - per-DOY means,
%     sigma from residuals in a circular 31-day window (Dec-Jan wrap)
%     with the sqrt(n/(n-1)) correction; raw-value window sigma leaks
%     seasonality (1.99 vs 1.5, measured). Tracks short-lived floods
%     monthly fields miss (Gouweleeuw et al. 2018).
%   - Guide: "Offline by design" - chains never fetch; local paths
%     plus shipped data/gravis freezes suffice, fetch* provisions
%     once (skip-if-present), three-layer path resolution documented.
%
% New in v3.13.0 (VDK/VADER decorrelation from monthly SINEX)
%   - shLowLevel.vdkApply: the Horvath et al. (2018) filter
%     (N + alpha*M)^-1 * N * x with the formal monthly normal-equation
%     matrix from the ITSG SINEX; never forms the filter matrix (one
%     Cholesky per month). Python-prevalidated: form-equivalence to
%     the tvANS Wiener family 1e-15, closed form for M = c*N, Wiener
%     MSE optimality on synthetics.
%   - shLowLevel.signalVarianceKaula: cyclostationary Kaula a*l^b per
%     calendar month from a pre-filtered series, with the EXACT
%     log-chi-square bias correction 0.5*(psi(k/2)-log(k/2)), k=2l+1
%     (without it the intercept biases 6% low - found by the Python
%     pre-validation).
%   - tools/dev/run_vdk_series.m: resumable batch driver for the full
%     series on another machine (~460 MB SINEX/month): skip-present
%     downloads, skip-existing outputs, cached signal model,
%     provenance log, built-in closed-form self-test and the paper's
%     acceptance probes.
%   - tvANSFilter help: the exact algebraic relation and the honest
%     division of labour between tvANS and VDK.
%
% New in v3.12.0 (fetch family: SINEX and ITSG background models)
%   - shLowLevel.fetchSINEX: ITSG monthly normal-equation SINEX from
%     TU Graz (the only public per-month SINEX source; COST-G
%     distributes none). Release routing as in fetchITSG; months are
%     REQUIRED - one n96 month is ~460 MB gzipped (verified live),
%     the full series ~120 GB, so no "all" convenience exists.
%   - shLowLevel.fetchITSGBackground: ITSG monthly background-model
%     means as .gfc (~1 MB); both eras carry dealiasing and the tide
%     models, the GRACE era adds atmosphere/ocean splits, c20,
%     degree-1, GIA, hydrology (palettes verified live). Explicitly
%     NOT the AOD1B GAX split - dealiasing ~ GAC, no GAD substitute.
%   - Naming now consistent: fetchITSG/fetchICGEM for series,
%     fetchGAX/fetchSINEX/fetchITSGBackground for products; all on
%     one shared robust loop (fetchFileSet: caps, pauses, 429 retry,
%     websave fallback) with shared month validation.
%
% New in v3.11.0 (obpChain; residual circulation separation; fetch fallback)
%   - shLowLevel.obpChain: GravIS-style ocean-bottom-pressure FIELDS
%     (GSM + GAD, filter after GAD, anomalies vs the GravIS window
%     2002/04-2020/03, land NaN). Bottom pressure keeps the air
%     column; oceanChain removes it - same core, one step apart.
%   - shLowLevel.eofSeparate + oceanChain SeparateCirculation=: the
%     residual ocean circulation is split from noise (area-weighted
%     EOF, North 1982 rule; Python-prevalidated incl. the
%     Marchenko-Pastur separability limit); sigMon becomes the
%     de-circulated noise RMS.
%   - httpFetch WebsaveFallback (default on): transport-level failures
%     (no HTTP status ever received) fall back to ONE websave attempt;
%     received statuses never fall back - the badStatus contract holds.
%   - The GAD filter-order note (GAD before the filter, GravIS order;
%     ocean mean filter-invariant, fields not) is now stated in both
%     ocean chain helps. addGADFolder is a shared private helper.
%   - tools/dev/machine_accept_v3110.m: sectioned acceptance script
%     for the next bridge session (suite, fetchGAX live, obpChain and
%     EOF acceptance, optional GravIS OBP grid cross-check).
%
% New in v3.10.1 (fetch robustness against ICGEM rate limiting)
%   - fetchGAX: MaxFailures cap (a machine-side failure no longer walks
%     through every remaining file at full retry cost), polite PauseSec
%     between downloads, and a single capped retry on HTTP 429 - the
%     ICGEM server rate-limits bursts (observed live on the acceptance
%     machine; the Retry-After-aware httpFetch waited so faithfully on
%     429 storms that calls exceeded the MCP bridge window).
%   - Full-restoration acceptance (COST-G RL02.1, 252 months, GAD/GAA
%     251 covered): annual amplitude 8.0 -> 9.1 mm (published
%     manometric range), trend +1.41 +- 0.03 mm/yr (moved 0.08 -
%     AOD1B is trend-free by construction), residual RMS 1.21 ->
%     1.49 cm (GAD restores real sub-annual variability).
%
% New in v3.10.0 (fetchGAX; full ocean restoration; API completeness)
%   - shLowLevel.fetchGAX downloads the AOD1B monthly means (GAA/GAB/
%     GAC/GAD) as ICGEM-converted .gfc from the GFZ series pages - 163
%     GRACE-era + 88 GRACE-FO files per product, centre-independent, so
%     one fetch serves every GSM series including COST-G. Listing and
%     download run through the Retry-After-aware fetch layer.
%   - oceanChain full restoration (Chambers & Willis 2010): GADFolder=
%     adds the model ocean signal on the coefficient level, the new
%     GAAFolder= subtracts the per-epoch ocean mean of the atmospheric
%     product; coverage is reported (nGadRestored/nGaaApplied), never
%     silent.
%   - API-table completeness: blank argument descriptions reduced from
%     616 to 0 - a toolbox-convention lexicon fills shared names as a
%     FINAL extraction pass (help text always wins), and 26 functions
%     received real new Inputs/Options help lines.
%   - chains_flow diagram regenerated with the ocean lane and the
%     fetchGAX feed (elbow routing, geometrically verified collision-
%     free including pads).
%
% New in v3.9.0 (ocean chain; Slepian application; gfc dot terms)
%   - shLowLevel.oceanChain (roadmap item 7): the ocean sibling of the
%     validated chains - gravisL2B corrections + GIA (ON: the ocean-
%     floor correction is a measured +0.89 mm/yr lever), gauss445,
%     area-weighted ocean-mean series, trend/annual fit, and the honest
%     residual noise proxy sigMon. Acceptance run (COST-G RL02.1, 252
%     months, |lat|<=66 crude-continent mask): +1.49 +/- 0.02 mm/yr
%     (+312 Gt/yr), annual amp 8.0 mm, sigMon 0.0121 m; the ocean MEAN
%     is filter-invariant to 0.01 mm/yr while unfiltered pixel
%     residuals are stripe-dominated (1.62 m) - the guide-V6 "trend
%     ocean RMS is mostly real signal" finding, now an API. GAD/GAA
%     restoration is a declared limitation (GADFolder= adds GAD where
%     available; AOD1B carries no secular trend by construction).
%   - shLowLevel.slepianProject (roadmap item 8): the application half
%     of the existing slepianBasis - project coefficient series onto
%     the ~Shannon leading tapers and back; regional analysis estimates
%     ~N coefficients instead of P. Cross-validated: the Gauss-Legendre
%     kernel reproduces an independent Python ring-quadrature reference
%     (30-deg cap, lambda_1 = 0.999981).
%   - ICGEM 'dot' secular lines (roadmap item 8) are read as 'trnd'
%     synonyms; files carrying them route through the line parser so
%     the bulk fast-path's group semantics stay untouched.
%
% New in v3.8.10 (name-value convention; quantity "none"; setup prefs)
%   - Name-value convention, codified and regression-tested: canonical
%     spelling is Capitalized for arguments-block options (Filter=,
%     OceanMask=) and lowercase for the legacy inputParser names and all
%     quantity strings ('ewh', 'geoid'); BOTH mechanisms tolerate any
%     casing, and testNVCasingToleranceAndConvention pins that tolerance
%     so refactorings cannot silently break existing scripts.
%   - quantity "none": dimensionless passthrough (kernel 1, no GM, R or
%     kn enter) in kernelFactors and every synthesis front door - the
%     clean way to synthesize raw coefficient fields such as kernels or
%     masks. Retires the documented workaround 'geoid' with GM = R = 1.
%   - setup_shAnalysis persists data locations: SeriesFolder, MasconFile
%     and the new GravisFolder/DDKFolder options are exported as SHX_*
%     env vars for the session AND stored via setpref('shAnalysis', ...)
%     so they survive restarts; the opt-in tests read env first, then
%     the preference.
%
% New in v3.8.9 (Retry-After fetchers; API-table completeness)
%   - shLowLevel.httpFetch / httpRetryDelay: every toolbox download now
%     retries 429/5xx with exponential backoff + jitter and honours the
%     server's Retry-After header (matlab.net.http - websave discards
%     the headers, which is why this was never possible before). All
%     fetchers route through the shared private webFetch; the HTML
%     listing reads retry with backoff. Motivated by the GravIS portal
%     503 bursts hit live during the V7/V8 validation.
%   - Opt-in real-data tests accept persistent MATLAB preferences as an
%     alternative to SHX_* env vars: setpref("shAnalysis",
%     "GravisFolder"|"SeriesFolder"|"DDKFolder", path).
%   - Guide V9 adds the standardChain flow diagram (the generic TN-based
%     core) next to the gravisL2B chain diagram.
%   - API reference tables: sizes are now always populated (help-declared
%     size, else n x 1 / n x m per convention) and descriptions are
%     merged from the help Inputs/Options sections plus toolbox-wide
%     convention texts; ~120 remaining blank descriptions across ~35
%     minor entries are on the roadmap.
%
% New in v3.8.8 (validated chains as one-call methods)
%   The executed GravIS validations (guide V1-V8) are now toolbox
%   methods with exchangeable inputs and full provenance reports:
%   - shLowLevel.gravisL2B        the shared GravIS Level-2B correction
%       core: aux-table C20/C30/C21/S21 + degree 1, NFIL mean, optional
%       GIA rate; uncovered epochs dropped and recorded
%   - shLowLevel.greenlandChain   -225.6 Gt/yr with the tested defaults
%       (published -231.1; guide V1-V6), 32 s end to end
%   - shLowLevel.antarcticaChain  -125.7 Gt/yr + the 25-basin table
%       (read V7 before interpreting the 14% vs the AWI joint product)
%   - shLowLevel.twsChain         river-basin TWS means in cm EWH;
%       amplitude ratio 1.001 and 0.032 cm/yr trend RMS vs GravIS
%   The sigma_trend noise policy requires a USER-SUPPLIED OceanMask
%   (oceanRMS contract) - the mask must be false over land (a bare
%   latitude band verifiably shifts Greenland by 5 Gt/yr).
%   The chains are turnkey: data/gravis ships frozen (2026-08-12) CC-BY
%   copies of the GravIS aux tables, mean, GIA rate and basin polygons
%   (shLowLevel.gravisDataFolder); guide chapter V9 documents the flow
%   with a diagram and runnable example calls, mirrored in the help
%   Example blocks and the HTML API reference. The executable audit
%   evidence behind V1-V8 is collected under tools/audit.
%
% New in v3.8.7 (validation: Antarctica + terrestrial water storage)
%   Guide chapters V7/V8 extend the GravIS validation beyond Greenland,
%   executed end-to-end on public data (2026-08-12):
%   - Antarctica: the unchanged chain gives -125.7 Gt/yr against the
%     AWI joint-basin -146.9 (span 2002-04..2023-02, n=217). The GIA
%     lever (+52.2 Gt/yr) and every stopping pathology reproduce
%     exactly; the official TU Dresden kernel GRID product shows the
%     SAME per-basin damping pattern (corr 0.905, RMS 8.2 vs this
%     chain's 5.9 Gt/yr) - the 14% is the distance between map-based
%     estimators and the Sasgen forward-modelling class.
%   - TWS: eleven river basins, full 252-month span, DDK3 declared for
%     GravIS' VDK blend: amplitude ratio median 1.001 (0.98..1.04);
%     trend RMS 0.225 -> 0.032 cm/yr once the ICE-6G_D GIA rate is
%     subtracted - the GravIS TWS product carries an undocumented GIA
%     correction (their Technical Note omits it).
%   - evalMask now rejects polygon vertices with |lat| > 90 with an
%     identified error: GeoJSON is [lon lat], the toolbox convention is
%     [lat lon], and a swapped input used to build a silent phantom
%     basin with plausible area (caught in the wild by the v3.8.6
%     zeroArea warning during this validation).
%
% New in v3.8.6 (independent audit: 20 findings fixed)
%   An adversarial audit (2026-08-12) reproduced the GravIS chain
%   end-to-end and attacked suite, docs, numerics and error handling.
%   Fixed here, each with a regression test:
%   - readMascon: attribute lookup is case-insensitive (CSR writes
%     'Units'), unrecognized time units ERROR instead of returning raw
%     days as decimal years, grids orient by dimension NAME (a square
%     grid was un-orientable by size), duplicate help Outputs removed.
%   - 22 shipped-fixture assumeTrue(isfile) -> verifyTrue: a renamed
%     fixture now FAILS ~20 tests instead of silently filtering them
%     (the testStandardChain mechanism, suite-wide this time).
%   - testPublishedTrendOnARealSeries now requires Greenland to LOSE
%     mass - it passed a sign-flipped trend before (demonstrated).
%   - shReadGFC warns on truncated files (record count vs the header's
%     max_degree triangle) and on corrupt values parsed to NaN.
%   - kernelFactors validates kn/hn length and 1+kn ~= 0 (CM-frame
%     k1 = -1 gave a silent Inf kernel).
%   - leakageCorrect rejects NaN inside Mask/ResidualRegion with the
%     cause (it used to report 'diverged'); NaN outside stays allowed.
%   - basinKernel warns on zero-area regions; gridScaling gained Quiet=.
%   - fetchLoveNumbers accepts the GROOPS folder's own headerless
%     single-column files (ak135); shCoefficients.read falls back to
%     the filename stem when the header has no modelname, so
%     shSeries.names is populated for COST-G and friends.
%   - sensitivityKernel help: the 2-16% margin and its ordering are
%     properties of the 1-D validation; on 2-D basins the audit
%     measured 1-42% with the ordering reversed. (The full-covariance
%     Noise path was CONFIRMED correct and is now pinned by a test.)
%   - testEigTrickEquivalence exercised zero toolbox code; it now
%     reconstructs tvANSFilter's output from its own op struct.
%   - Three HTML pages carried content after </html>; repaired, and
%     doc_sync_audit now checks document structure.
%   - README/Contents pin the GravIS span dates next to the 2.8%.
%
% New in v3.8.5 (documentation audit)
%   html/shAnalysis.html - the page a user lands on from `doc
%   shAnalysis` - said v3.0.0 for eight releases. Nothing checked it:
%   the version gate only looked at pages that already carried one. It
%   now carries a stamp AND is gated, with a negative test confirming
%   the gate fires.
%   NOTICE listed the GFZ TN-13 fixture as RL06.3.txt; it ships as
%   RL06_3.txt. Same typo that kept testStandardChain filtered.
%
% New in v3.8.4 (a test that was never running)
%   testStandardChain asked for TN-13_GEOC_GFZ_RL06.3.txt while the
%   shipped fixture is RL06_3 (underscore). assumeTrue turned the typo
%   into a filtered test - indistinguishable from a passing one - so the
%   canonical pipeline went unverified from v3.1.0 to v3.8.3. Fixture
%   checks are verifyTrue now, and the newly-running test immediately
%   exposed a RelTol of its own that had never been exercised. Guide
%   gallery figures: checked, all ten are real, the "placeholder"
%   roadmap item was stale.
%
% New in v3.8.3 (opt-in test data, documented and settable)
%   setup_shAnalysis SeriesFolder= and MasconFile= export
%   SHX_SERIES_FOLDER and SHX_MASCON_FILE for the session (startup.m to
%   keep them) and report them in summary.testData. Documented in the
%   guide's new Installation chapter, on the setup page and in the
%   README. Both variables are optional: without them the two tests
%   report as FILTERED rather than passing silently.
%   SHX_SERIES_FOLDER stays deliberately unwired from the persistent
%   data folder - a routine acceptance run must never reach for a
%   network or archive drive.
%
% New in v3.8.2 (readMascon verified against a real product)
%   The synthetic fixture writes lwe_thickness(lon, lat, time); the real
%   JPL RL06.3Mv04 file writes (time, lat, lon). A size check passes
%   either way, so the opt-in test (SHX_MASCON_FILE) checks GEOGRAPHY:
%   Greenland must lose mass over the GRACE era. Verified: -124 cm, with
%   the largest change at 76 N / 291 E, northwest Greenland.
%   Also documented: mascon latitudes are GEODETIC by convention while
%   this toolbox synthesises on GEOCENTRIC ones - 0.19 deg apart near
%   45 deg, about 21 km or 0.4 of a mascon cell.
%
% New in v3.8.1 (README rewritten; the COST-G SINEX question settled)
%   The README still described compat/ and tests/legacy - both removed in
%   v3.0.0 - and its newest section was "New in v2.5", nine releases
%   behind. Rewritten against the current state.
%   COST-G SINEX fetcher: DROPPED from the roadmap, no such source
%   exists. COST-G combines on the normal-equation level internally and
%   publishes combined SOLUTIONS through ICGEM, but does not distribute
%   per-month SINEX or covariances; Dahle et al. (2025) state the
%   variance-covariance information is not available for the combined
%   solutions. Use ITSG (which does publish monthly SINEX) or the GravIS
%   VDK-filtered Level-2B products. Documented on the readSINEX page.
%
% New in v3.8.0 (Live Script tutorials - roadmap item 10)
%   shLowLevel.makeTutorials writes one Live Script per demo case into
%     tutorials/, GENERATED from the demo registry rather than hand
%     written: demo_shAnalysis is the single source of truth, and a
%     parallel set of .mlx (a binary zip) would drift invisibly.
%     Conversion uses an internal MATLAB API; when it is unavailable the
%     .m files are still written and a warning says how to convert by
%     hand.
%   CI: timeout-minutes on the job (20) and on the MATLAB setup step
%     (8). The MathWorks setup step hung four times in one session,
%     taking 5-10 minutes against its usual 1.2, each time needing a
%     manual cancel; now it fails fast and visibly.
%
% New in v3.7.0 (tailored sensitivity kernels)
%   shLowLevel.sensitivityKernel - a basin kernel that trades leakage
%     against propagated noise EXPLICITLY (Swenson & Wahr 2002; the
%     construction behind the ESA CCI and GravIS gridded ice products),
%     subject to unit response over the basin. Closed form via one
%     Lagrange multiplier: no iteration and no stopping criterion, which
%     is its advantage over the forward-modelling route. Sweep Alpha,
%     plot leakage against noise, take the L-curve corner.
%     Python-validated in tools/dev/validate_senskernel.py.
%
% New in v3.6.0 (self-consistent degree 1)
%   shLowLevel.estimateDegree1 - geocentre coefficients from GRACE
%     itself plus an ocean model (Swenson, Chambers & Wahr 2008), the
%     method GravIS and geogravL3 use. Removes the dependency on someone
%     having published a TN-13 series for your exact Level-2 product,
%     and returns the same struct layout - sigmas included - so it drops
%     into addDegree1 unchanged. kn is REQUIRED (degree-1 Love numbers
%     are frame dependent). Python-validated in validate_degree1.py.
%
% New in v3.5.2 (the discrepancy principle, made to actually fire)
%   leakageCorrect ResidualRegion= (default: the Mask). The residual
%     must be measured WHERE THE MODEL IS RESPONSIBLE: with a regional
%     mask the data still contains every other mass signal on the globe,
%     so a global residual never reaches the noise level and the
%     principle silently never fires. info.residualRMSGlobal exposes it.
%   Guide: the GravIS comparison rerun with a real stopping criterion.
%     -224.6 Gt/yr against a published -231.1 (2.8 %; span
%     2002-04..2023-02, n = 217 - quote the dates with the number,
%     the reference moves ~10 Gt/yr over two years of end date),
%     replacing an
%     earlier -227.4 that came from stopping at an arbitrary iteration.
%     Also: the noise level must match the quantity inverted - the
%     open-ocean RMS of a TREND field is mostly barystatic sea level,
%     not error; the trend's noise is sigma_monthly / sqrt(Sxx).
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
