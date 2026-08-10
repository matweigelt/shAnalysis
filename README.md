# shAnalysis

![CI](https://github.com/matweigelt/shAnalysis/actions/workflows/ci.yml/badge.svg)
![MATLAB](https://img.shields.io/badge/MATLAB-base%20only%2C%20R2026a%20tested-orange)
![Platforms](https://img.shields.io/badge/platforms-Windows%20%7C%20Linux-lightgrey)
![Latest tag](https://img.shields.io/github/v/tag/matweigelt/shAnalysis?label=version)
![License](https://img.shields.io/github/license/matweigelt/shAnalysis)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21871298.svg)](https://doi.org/10.5281/zenodo.21871298)

Class-based MATLAB toolbox for spherical harmonic (Stokes coefficient)
analysis of GRACE/GRACE-FO, GOCE, and static gravity field models.
Merges the v1 function toolbox (reading, spectra, destriping, Gaussian
smoothing, synthesis) with the tvANS time-variable anisotropic Wiener
filter into three user-facing classes.

Requires MATLAB R2021a or newer (`arguments` name=value syntax).

## Documentation

The complete **[workflow & theory guide (PDF, 75 pages)](docs/shAnalysis_workflow_guide.pdf)**
covers the processing chain, filter theory, uncertainty propagation, and the
generated API reference (Part IV) with typed inputs/outputs and executable
examples for every public function. HTML help for MATLAB's `doc` browser
ships in `html/`.

## Layout

```
shAnalysis/
  shCoefficients.m   one coefficient set (value class, immutable)
  shSeries.m         time series of coefficient sets
  shClimatology.m    bias/trend/annual/semi-annual model
  runAllTests.m      full validation suite
  +shx/              numerical internals (stable, package-qualified)
  compat/            v1 function signatures (add to path only if needed)
  tests/             correctness / contract / robustness / performance
  tests/legacy/      unmodified v1 suite, runs against compat/
  html/              doc pages (overview, classes, v2.1 function pages)
  tests/test_data/   fixtures incl. two real ITSG monthly solutions
                     (GRACE 2008-04, GRACE-FO 2025-12, TU Graz)
```

## Quick start

```matlab
cd shAnalysis
setup_shAnalysis(Permanent = true, Download = "core")   % path + TN files
% later, monthly refresh of the growing TN files (safe, parse-verified swap):
setup_shAnalysis(Download = "core", Update = true)      % or shx.fetchTN(Update = true)
% levels: "none" (path only) | "core" (TN-13/14) | "filters" (+DDK)
%         | "starter" (+ITSG months); DryRun=true shows the plan
ts  = shSeries.read("GSM-2_*.gfc");     % monthly solutions, epoch-sorted
gad = shSeries.read("GAD-2_*.gfc");
ts  = ts.restore(gad);                  % GSM+GAD, epoch-matched
dts = ts - ts.mean;                     % anomalies about the mean field

[clim, resid] = dts.climatology(Robust=true);   % shClimatology + residuals
A = clim.amplitudeMap("annual", -89:89, 0:359, "ewh", kn);

tsA = dts.destripe(minOrder=6).gaussian(300);   % classic chain
[tsF, op] = dts.filter("tvANS", SignalMode="inhomogeneous");
avg = tsF.basinAverage(B, Deconvolve=true, Op=op);   % unbiased basin series
```

## Running the tests

```matlab
results = runAllTests;        % writes tests/runAllTests_results.txt
```
The complete console output plus a per-test table always land in
`tests/runAllTests_results.txt` (fixed name, overwritten) - upload that
file for review; pasted console text gets truncated.

Note on the legacy suite: one v1 test (kernel ratios) carries a documented
patch - the v2.1 scaled Legendre engine returns EXACT zeros at symmetry
points (e.g. Pbar_30 at the equator), turning zonal ratio samples there
into 0/0; the test now masks zero-field points and separately verifies
they are consistently zero. The v1 engine's ~6e-17 colatitude rounding
fuzz had masked this.



## New in v2.5

**Uncertainty exactness.** (a) With hard constraints the tvANS posterior
sigma is now the EXACT diagonal of (W-I)S(W-I)' + s*WNW' (O(P^2 q) in the
eig basis; Python-validated to 2e-15; constrained-direction identity
Ac'*Cov*Ac = s*Ac'*N*Ac exact). Honest correction: the earlier formula
UNDERESTIMATED along constrained directions (ratios down to ~0.72) - the
old "upper bound" note had the direction wrong. (b) The basin
deconvolution sigma now carries the OLS parameter uncertainty of the
restored deterministic part (design leverage x per-coefficient residual
variance carried on the operator; MC: emp/pred 1.18 -> 1.00, 1.26 -> 1.08
at seasonal leverage peaks), and with constraints the deconvolution noise
covariance uses the exact constrained operator. (c) `ARCorrect` applies
the Kendall first-order r1 bias correction r1' = r1 + (1+3r1)/T (MC:
1.065 -> 1.028 at phi=0.6/T=120; 1.152 -> 1.079 at T=60).

**`shx.readLoveNumbers` (new).** Load-Love-number reader: 2-column and
headered layouts, commented headers, `Columns=` override for the
ambiguous headerless 4-column case (Farrell h l k vs k h l), sparse-table
pchip interpolation in log(1+n) (validated 1.3e-4), and degree-1
reference-frame conversion after Blewitt 2003 (CE->CF/CM, CF<->CM;
reproduces the published PREM CF values and k1(CM) = -1 exactly; ->CE
refused as non-recoverable). Still no hardcoded physics - a reader only.

**Real provider files.** CSR and JPL TN-13 RL06.3 ship as fixtures
(256 paired months each, 2002-04..2026-04; cross-provider C10
correlation vs GFZ 0.995) with a chain test on the real ITSG month. All
eight released DDK Wbd binaries were parsed and validated (241 blocks,
monotone gain ordering DDK1..DDK8); a discovery test re-checks whatever
subset is present locally (test_data + dataFolder/DDK). COST-G SINEX and
real mascon files were not openly reachable (AIUB/CSR 503, PODAAC auth) -
retry list.

**Continuous integration.** `ci/shanalysis-ci.yml` is a ready GitHub
Actions workflow: copy it to `<repo-root>/.github/workflows/`. On a *public*
repository, MathWorks' `setup-matlab` action licenses base MATLAB
automatically on GitHub-hosted runners - no token, no cost (and shAnalysis
needs nothing beyond base MATLAB). On a *private* repository, request a
batch licensing token via the MathWorks Batch Licensing Pilot and store it
as the `MATLAB_BATCH_TOKEN` secret; alternatively register your own machine
as a self-hosted runner and use the local installation. The suite is
CI-safe: all network tests are offline (skip/mirror patterns) and fixtures
ship in `tests/test_data`.

**Comparison suite (v2.6.0).** `shx.compareSolutions(g1, g2)` and
`shx.compareSeries({ts1, ts2, ...})` (or `g.compare(g2)` / `ts.compare(...)`)
return the standard metric set for comparing solutions: difference degree
amplitude with agreement crossover, degree correlation, area-weighted spatial
statistics, chi2/dof error realism against formal sigmas, NSE, AR(1)-corrected
correlation significance, trend/annual component differences with combined
sigmas, and (N >= 3) three-cornered-hat noise levels per solution.
`Plot = true` adds 4-panel figures incl. a Taylor diagram
(`shx.taylorDiagram` also stands alone). Each metric is its own function for
single-test use: `diffSpectrum`, `spatialStats`, `nashSutcliffe`,
`effectiveCorr`, `threeCorneredHat`.

**`setup_shAnalysis` (new).** One-call setup: path (temporary or
`Permanent=true` via savepath, with a startup.m fallback warning on
managed machines) plus cumulative download levels - `"core"` (TN-14 +
TN-13 GFZ/CSR/JPL via the new `shx.fetchTN`, every file verified by
parse before acceptance), `"filters"` (+DDK subset), `"starter"` (+ITSG
months). Idempotent (existing files skipped and re-verified), per-level
failure isolation, `DryRun=true` returns the plan with zero side
effects (contract-tested). The real GSFC TN-14 file ships as a fixture
(258 windows 2002.25-2026.41, pinned test).

**Robust performance test.** `testCacheSpeedup` uses a Legendre-dominated
workload (n120, 361 rings), min-of-N timing and a 0.75 margin - immune to
machine load (the v2.4.2 run failed only because the host was ~2x loaded).

**Python figure port rebuilt** (`shx_port`): scaled Legendre, gfc/Wbd
readers, spectra, filter gains, triangle rendering with the v2.4.2 masks;
self-validates on import (sum rule, parity, filter identities). Guide
figures D01/D04 regenerated from the real shipped data; new D04b figure
shows the fixed full-sine-wing difference triangle on real data.

## New in v2.4.2

**License & attribution.** The toolbox now ships an MIT `LICENSE`
(copyright Matthias Weigelt; the bundled DDK3 Wbd file carries its own MIT
license, ITSG/GFZ test fixtures remain their providers' property). Every
function help text and every HTML doc page closes with
*"Developed by Matthias Weigelt with the help of Claude (Fable 5)."*

**Coefficient triangle (bug fix).** The `RefC`/`RefS` difference branch of
`plotSHCoeffTriangle` masked the sine wing with `tril(...,-1)`, blanking
the sine **sectorals** S_nn and rendering the left wing one column
narrower than the C wing in every difference triangle (demo D04/D07
figures). Fixed: the sine mask is now the full lower triangle n >= m >= 1
including sectorals; a rendered-CData regression test pins every cell.

**ICGEM temporal catalogue (bug fix + extension).** `shx.listICGEM(Type=
"temporal")` previously matched only `getseries` links in the static HTML
- on the current page these are three centers' release notes, so the table
had 3 rows. It now parses the `/sp/` series-page catalogue: **all ~70+
series of all centers** (groups 01_GRACE, 02_COST-G, 03_other, 04_SLR)
with columns `group, center, series, path, url, zip`. New `Series=`
option: pass a catalogue `path` to list the **individual files** of one
series (single monthly `.gfc` files download directly via `websave`) -
the earlier "temporal browsing is JavaScript-only" caveat is obsolete.
Both parsers are fixture-tested offline (`icgem_temporal_fixture.html`,
new `icgem_series_fixture.html`, both real captures of 2026-08-07).

## New in v2.4.1

Plots: triangle center line removed; spectrum plots use a **linear** degree axis (log y).
`plotSHSpectrum` gained `Quantity=` (amplitude | rms | variance | cumamplitude | cumrms |
cumvariance) and order-domain support via **`shx.shOrderRMS`** (the striping axis; degree/order
marginals verified to sum to the same total). `shx.shDegreeRMS` carries `cumVariance`/`cumRMS` now.

Data management: **`shx.dataFolder`** (persistent user folder; all fetchers store beneath it),
**`shx.fetchDDK`** + `shx.readDDK("DDK<n>")` for **all eight** released DDK filters
(DDK3 still ships in the package), **`shx.listICGEM`** (180+ static models as a table;
fixture-tested parser) and **`shx.fetchICGEM("<name>")`** for direct .gfc downloads.
`shx.fetchITSG` now targets `dataFolder/itsg_series` by default and gained
**`Product="daily"`**: the ITSG daily Kalman-smoother solutions (one ICGEM `.gfc` per day,
n40 only, formal errors; ~365 files/~30 MB per year into `dataFolder/itsg_daily`; daily
filename epochs parse automatically — a real daily file ships as test fixture).

## New in v2.4

**Visualization** — `shx.plotSHMap`/`g.map` (divergent scale, coastlines, analytic Hammer
projection, no Mapping Toolbox), `shx.plotBasinSeries` (1σ band, mission-gap patches,
AR(1)-trend annotation), signed **triangle differences** (`'RefC'/'RefS'`),
`shx.plotCovariance` (block-boundary overlay), `plotSHSpectrum` Kaula + crossover
overlays, `shx.writeAnimation` (MP4, fixed robust color scale).

**Computation** —
- `shx.errorMap`: analytic σ-maps from full covariances (Cholesky; matches `mcPropagate`).
- `shx.eofAnalysis`: EOF modes as `shCoefficients`, unit-variance PCs, variance fractions.
- `shx.basinScaling`: forward-modelled gain factors; operators as tvANS struct, matrix, or handle.
- `shSeries.trendBreaks`: continuous piecewise-linear trends + F-test per coefficient
  (p-values via `betainc`, base MATLAB; null rejection calibrated 5.1% @ 5%).
- `shx.shFanFilter` / `g.fan` / `ts.fan`: Han degree×order Gaussian.
- `shx.combineCenters`: **multi-center VCE** with one factor per (center, month) and
  partial redundancies r_c = P − tr(H·w_cN_c⁻¹); robust block-median option;
  `info.interCenterCorr` honesty diagnostic; posterior sigmas. Combine GSM **before**
  TN-14/TN-13. Common-mode errors remain invisible — posterior sigmas are a lower bound.
- `shx.shSynthesisGradientTensor`: all six NEU components at altitude; Laplace trace
  ≤ 7e-16 built in as self-check; exact second Legendre derivatives (stencil²).
- `shx.seaLevelFingerprint`: elastic sea-level equation on the exact quadrature pair;
  mass conserved to machine precision; classic near-field-fall/far-field-excess pattern.
  Elastic only, no rotation, fixed coastlines (documented).
- `shx.writeGrid`: CF-style netCDF (2-D/3-D), roundtrip-tested.
- **Reference systems**: `shx.normalFieldCS` computes the even zonals from the WGS84/GRS80
  **defining constants** (Heiskanen-Moritz closed form; matches NIMA TR8350.2 to all published
  digits — no coefficient tables). `g.subtractNormalField()` rescales the ellipsoid field to the
  model's (GM, R) first — WGS84 GM/a **differ** from the ICGEM values, skipping this costs mm.
  `shx.rescaleGMR`/`g.toReference` convert between conventions exactly (C′=C·(GM₁/GM₂)(R₁/R₂)ⁿ;
  field invariance verified); `+`/`-` now point to `toReference` on GM/R mismatch.
  Permanent-tide convention is never converted silently (documented).
- **Demo registry**: `demo_shAnalysis("list" | "core" | "all" | ["D04","D15"], Visible=false)` —
  16 selectable cases covering every capability (D01 diagnostics … D16 export); real test-data
  files used when present, synthetic fallbacks everywhere. The workflow guide gained a
  **Demo Gallery** chapter: case table plus nine rendered figures with explanations
  (generated with the Python validation port; the DDK3 panel reads the real Wbd binary).
- **Real data on demand**: `shx.fetchITSG(2010:2016)` downloads ITSG monthly solutions
  (websave, automatic release routing, mission-gap aware) into `tests/test_data/itsg_series`;
  `shSeries.fromFolder(dest)` loads them. Demo data policy: D01–D04 use the shipped real
  files (incl. the real 2008→2025 GRACE-FO−GRACE difference with real stripes), D05/D06
  switch to real series once fetched, D12–D14 stay synthetic **by design** (recovery of
  known truth). Gallery figures G1/G2 now show the real ITSG month and the real mass-change map.
- `shx.pctile`: percentile in base MATLAB (audit caught `prctile`/Statistics Toolbox in two
  plot functions; removed).

## New in v2.3

New quantities in `kernelFactors` (single source of truth for synthesis AND analysis):

- **`surface_density`** [kg/m²] — `R·ρ_ave/3·(2n+1)/(1+k'_n)`; exactly `ρ_w ×` the EWH kernel.
- **`bottom_pressure`** [Pa] — `g₀ ×` surface density (`g₀ = GM/R²`); the ocean-bottom-pressure product.
- **`deformation_up`** [m] — `R·h'_n/(1+k'_n)`; elastic vertical load deformation (GRACE↔GNSS uplift).
- **`gravity_gradient_rr`** [1/s²] — `(GM/R³)(n+1)(n+2)`; radial gravity gradient T_rr (1 E = 1e-9 s⁻²).
- **`Height=`** — upward continuation `(R/r)^p` to satellite altitude for potential-type
  quantities (validated against −d/dr and d²/dr² of the continued potential); surface
  quantities reject it with `shSynthesis:heightInvalid`.

Horizontal deformation is not a degree factor — it needs the SH gradient:

- **`shx.shSynthesisDeformation` / `g.deformation`** — up/north/east elastic load
  deformation (Wahr et al. 1998), grid or station-list (`Mode="points"`) evaluation,
  `LatType="geodetic"` for GNSS coordinates. Exact ∂P̄/∂φ via `shx.legendreALFDeriv`
  (frozen identity, calibrated per (n,m) against numerical differentiation; north/east
  validated to ~5e-10). East is NaN at the poles; degree 0 excluded, degree 1 opt-in
  (consistent geocenter/frame handling required).

All Love numbers (k', h', l') remain strictly user-supplied from one consistent loading model.

## New in v2.2

Scientific chain completions:

- **GIA correction** — `ts.removeGIA(gia)` / `clim.removeGIA(gia)` with user-supplied
  rate models (gfc); both routes agree exactly on the trend. Model treated as exact
  (documented; compare several models for a spread estimate).
- **DDK filters** — `shx.readDDK` reads the released binary `Wbd_*` files natively
  (Rietbroek/Kusche BIN format, endian autodetect, verified to 4.4e-16 against the
  documented DDK3 reference values; `Nmax=` truncation for n96 series), plus `.mat`
  containers and an ASCII exchange format; `applyDDK` on coefficients and series.
  The real DDK3 file (`Wbd_2-120.a_1d12p_4`, MIT-licensed, from
  github.com/strawpants/GRACE-filter) ships in `tests/test_data`.
- **Basin kernels** — `shx.basinKernel(idx, polygon|mask|fun, BufferKm=, TaperKm=)`
  replaces manual `Y'*(w.*mask)`; exact great-circle buffering, Jekeli spectral taper.
- **Slepian localization** — `shx.slepianBasis`: concentration eigenproblem, Shannon
  number, orthonormal tapers; the well-posed alternative to Kaula-regularized regional fits.
- **AR(1)-corrected significance** — `ts.climatology(ARCorrect=true)` inflates
  coefficient sigmas by `sqrt((1+r1)/(1-r1))`; MC-validated (white-noise sigmas
  underestimate trend scatter 2.0x at phi=0.6; corrected ratio 1.06).
- **Geodetic latitude** — `shx.geodetic2geocentric` / inverse, plus
  `LatType="geodetic"` in `synthesis`, `analysis`, `shAnalysisGrid`. Mascon and map
  grids are geodetic; forgetting the conversion biases mid-latitudes by ~21 km.

Scale and statistics:

- **Memory-banded synthesis** — `'MaxMemGB'` (default 4): latitude-banded Legendre
  streaming; nmax 2190 on ordinary RAM instead of a 7 GB stack (chunked == monolithic
  exactly).
- **Per-order-band VCE** — `tvANSFilter(..., VCEBands=[0 16 33 61])`: banded monthly
  noise factors (block path), consistently honored in filtering, posterior sigmas,
  and basin deconvolution.
- **Monte-Carlo propagation** — `shx.mcPropagate(fun, g, Cov=, Idx=)`: empirical
  sigmas for any functional, from per-coefficient sigmas or full SINEX covariances.
- **Streaming SINEX** — `readSINEX(..., Only="estimate")` reads the estimate block of
  a 460 MB gz NEQ SINEX in seconds (Java gz stream; gunzip fallback without JVM).
- **Mascon reader** — `shx.readMascon` (base-MATLAB `ncread`), JPL/CSR layout
  auto-detection, "days since" → decimal years. Comparison workflow in the guide.
- Property-based synthesis↔analysis roundtrip tests; performance log
  `tests/perf_log.csv` (date, release, benchmark, seconds — machine-local record).

Full theory, best practices, and step-by-step guides: `shAnalysis_workflow_guide.pdf`.

## New in v2.1

```matlab
% SH ANALYSIS: gridded data -> Stokes coefficients (inverse of synthesis)
gHat = shCoefficients.analysis(grid, lat, lon, 60, quantity="ewh", kn=kn);
%   ring grids: exact fast per-order solver (roundtrip ~1e-15)
%   scattered points: chunked least squares, Kaula=1 for under-determined sampling

% complete processing chain incl. geocenter
g = shCoefficients.read("GSM-2_...gfc").applyTN14("TN-14.txt").addDegree1("TN-13.txt");

% gap-aware series, tidal alias periods, significance
ts   = ts.dropNaN();                                  % GRACE<->GRACE-FO gap
clim = ts.climatology(Periods=[161/365.25 3.66 7.48]);% S2, K2, K1 aliases
tr   = clim.trend();                                  % tr.sigmaC: 1-sigma per coeff

% block-diagonal tvANS: identical results, Lmax ~ 120 tractable; uncertainties
[tsF, op, info] = ts.filter("tvANS", Blocks="auto");  % tsF.sigmaCs/Ss posterior 1-sigma
[avg, out] = tsF.basinAverage(B, Deconvolve=true, Op=op);  % out.sigma (K x T)

% fast synthesis + deep expansions
grid = g.synthesis(-90:90, 0:359, Method="fft");      % 5-20x at 1 deg
P = shx.legendreALF(2190, deg2rad(80));               % scaled recursion, stable

% released covariances / normal equations (ITSG/COST-G SINEX)
snx = shx.readSINEX("solution.snx", Output="covariance", Index=shx.shIndex(60));
[tsF, op] = ts.filter("tvANS", NoiseCov=snx.M);
```

Single fields work the same way:

```matlab
g  = shCoefficients.read("GSM-2_2024032-2024060_GRFO_UTCSR_BA01_0600.gfc");
g  = g.applyTN14("TN-14_C30_C20_SLR_GSFC.txt");
d  = g - g0;
[grid, lat, lon] = d.destripe.gaussian(300).synthesis(-89:89, 0:359, ...
                       quantity="ewh", kn=kn);
```

Every object carries a provenance log (`.history`, shown by `disp`,
written into exported gfc headers).

## Validation

`runAllTests` runs five suites: numerical correctness (golden values,
analytic identities, Python-cross-validated formulas), API contract
(error identifiers, immutability, output dimensions), robustness (NaN
stacks, missing GAX epochs, degenerate basins, corrupted files, poles),
the complete unmodified v1 suite against the compat layer, and
performance benchmarks (including a Legendre-cache speedup assertion).

All tvANS numerics (Wiener MSE identity, eigen-trick equivalence,
constraint identities, unbiased basin deconvolution) were additionally
cross-validated in Python (numpy/scipy) before the MATLAB port.

v2.1 additions, Python-validated before implementation:
posterior variance formula diag((I-W)S) = (V.^2)(s*lam/(lam+s)) to 3.5e-16;
basin-average covariance vs Monte Carlo to 0.7% (4e4 samples);
block-diagonal path identical to the full eigendecomposition (1.8e-15);
FFT synthesis vs direct sum 2.9e-13 (incl. aliased-order bin accumulation
and nonzero start longitude); latitude parity relation exact; scaled
Legendre recursion: sum rule at n=2190 to 3e-11 (lat 89.99 deg), agreement
with the plain recursion 8.4e-15 at nmax=120; ring-grid analysis roundtrip
8.7e-15 incl. quantity kernels/weights/lon0; scattered-LS roundtrip 6.4e-15;
OLS sigma Monte-Carlo ratios 0.99-1.02.

## Known limitations / roadmap

* Both provider-format parsers are verified against REAL data: SINEX
  against an ITSG-Grace2018 n96 monthly SINEX (TU Graz; CODE=degree,
  SOLN=order, PT='--'; truncated 12x12-NEQ fixture in test_data), and
  TN-13 against a real GFZ RL06.3 geocenter file (256 months 2002-2026,
  full fixture in test_data, incl. a real-data degree-1 chain test on
  the ITSG 2008-04 solution and hard expected values).
* Posterior sigmas ignore the hard-constraint projection (slight upper
  bound along constrained directions); basin sigmas treat the unfiltered
  deterministic part as exact.
* ICGEM 1.0 behavior change (documented): 8-column acos/asin lines are
  now read as EIGEN-style "period in the last column"; >=9-column lines
  keep the pre-v2.1 t0/t1 reading bit-for-bit.
* No load Love numbers, GIA models, or other physical tables are
  hardcoded anywhere - EWH synthesis and amplitude maps require an
  explicit `kn`; TN-14 tables are parsed from user-supplied files.
* `shx.readTN14` assumes the standard 10-column TN-14 layout; TN-13
  (degree-1) parsing is on the roadmap.
* Empirical noise covariances (tvANS default) are signal-contaminated at
  low degrees; supply a full covariance via `NoiseCov` when available.
  A SINEX/ITSG covariance reader is on the roadmap.
* FFT-along-longitude synthesis (regular global grids) planned for v2.1.

Claude (Fable 5), 2026-08-07.
