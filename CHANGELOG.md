# Changelog

All notable changes to shAnalysis. The version line in `Contents.m`
(read by `shLowLevel.version` and MATLAB's `ver`) is the single source of truth.

## [3.6.0] - 2026-08-11

### Added
- `shLowLevel.estimateDegree1`: geocentre coefficients (C10, C11, S11)
  estimated from GRACE itself plus an ocean model, following Swenson,
  Chambers & Wahr (2008) - the method both GravIS and geogravL3 use.
  This removes the dependency on somebody having published a TN-13
  series for your exact Level-2 product and release, and guarantees the
  degree-1 terms are consistent with the solutions actually processed.
  The output struct matches `readTN13` exactly, sigmas included, so it
  drops into `addDegree1` unchanged. `kn` is required: degree-1 load
  Love numbers are frame dependent and are never assumed.
- Python-validated (`tools/dev/validate_degree1.py`): exact recovery of
  a known geocentre from a synthetic land/ocean world, 1.5-4% error at
  2 mm noise, a demonstrable bias when the ocean model is omitted, and
  a rising condition number on a degenerate ocean domain.

### Notes on three things the acceptance machine caught
- A CONSTANT nuisance column is carried beside the three patterns. The
  ocean residual generally has a degree-0 component (the observed field
  is degrees 2+, and any mass imbalance between data and ocean model is
  an offset), and without a constant to absorb it that offset lands on
  C10, which has a large ocean mean: a **153% error** before the fix.
- The design matrix is column-equilibrated before `cond` is reported.
  The constant column has a completely different scale from the basis
  patterns, so an unscaled condition number reports the units rather
  than the geometry - 3.8e7 on a problem whose geometry is fine. After
  equilibration: 1.3 for a global ocean, 58 for a polar-only one.
- `addDegree1` REQUIRES sigma fields, so the first version was not the
  drop-in its own help text claimed. Formal sigmas now come from the
  ocean-fit residual propagated through the least-squares solution.

## [3.5.2] - 2026-08-11

### Fixed
- `leakageCorrect` gained `ResidualRegion=` (default: the `Mask`). The
  discrepancy principle added in 3.5.0 measured the residual GLOBALLY,
  and with a regional mask that can never reach the noise level: the
  model describes one region while the data contains Antarctica,
  Alaska and global hydrology. On the GravIS Greenland case the
  residual was 0.0085 globally against 0.0008 near the mask, with a
  noise level of 0.0024 - so the principle silently never fired and
  every run stopped at `MaxIter`. `info.residualRMSGlobal` now exposes
  the gap; a large one means unmodelled signal elsewhere.

### Measured
- The GravIS Greenland comparison, rerun with a real stopping criterion:
  **-224.6 Gt/yr** against a published -231.1, i.e. **2.8%**. This
  replaces -227.4 (1.6%), which came from stopping at an arbitrary
  iteration count. The closer number was a coincidence of where the
  iteration happened to be; a result that depends on an arbitrary
  choice is not a result.
- The noise level must match the quantity being inverted. The
  open-ocean RMS of a TREND field is 0.0024 m/yr, but most of that is
  real barystatic sea-level rise, not error - using it stops after 4
  iterations and gives -201.0 Gt/yr, 13% low. The trend's noise is the
  monthly noise propagated through the fit, sigma_monthly / sqrt(Sxx) =
  0.0115 / sqrt(8353) = 0.000125 m/yr, nineteen times smaller. That
  stops at 82 iterations and gives the -224.6 above.

## [3.5.1] - 2026-08-11

### Added
- `shLowLevel.oceanRMS`: area-weighted RMS of a field over the OPEN
  ocean - ocean points more than `MinDistanceKm` (default 1000) from the
  nearest non-ocean point. Far from land a GRACE field should contain
  almost no real signal, so what remains is error; this is the standard
  noise metric processing centres quote, and the value `leakageCorrect`
  needs for `NoiseLevel`. Without it, v3.5.0's discrepancy principle had
  no natural input.
- The ocean mask is USER-SUPPLIED, deliberately. Base MATLAB's
  `coastlines` data set is not present in every installation (it is
  absent on the acceptance machine, which is why `plotSHMap` carries a
  fallback), and a wrong mask silently changes the number.
- Area weighting is on by default: an unweighted RMS over a lat/lon grid
  over-counts the polar rows. On white noise the weighted and unweighted
  values agree, which is exactly why the difference goes unnoticed until
  the field has structure - the test pins both cases.
- Python-validated (`tools/dev/validate_oceanrms.py`): the erosion moves
  the boundary by exactly d/R against an analytic spherical cap, and the
  weighted mean of cos^2 over the sphere comes out 0.6666 against 2/3.

## [3.5.0] - 2026-08-11

### Fixed (correctness of a shipped result)
- `leakageCorrect` solves an ill-posed inverse problem and therefore
  SEMICONVERGES: the error against the truth falls, reaches a minimum,
  then rises as the iteration begins to fit noise - while the residual
  keeps shrinking. The previous stopping rule chased a small solution
  step, which is exactly wrong. Measured on the reference problem
  (`tools/dev/validate_stopping.py`): the final solution is **361x**
  worse than the best one, `Tol = 1e-3` stops **89x** past the optimum,
  and the shipped default `Tol = 1e-4` never triggers at all.
- `NoiseLevel=` enables the discrepancy principle (Morozov): stop once
  the residual reaches the noise level of the data, since fitting more
  closely than the noise is fitting noise. It lands within 2% of the
  optimum. `Tau=` (1.2) is the safety factor; 1.2-1.5 is best, 1.0
  slightly overfits. The natural estimate for `NoiseLevel` is the
  open-ocean RMS of the field, which is what processing centres use.
- On a noisy version of the reference disc, verified on the acceptance
  machine: the discrepancy principle stopped after 2 iterations and was
  **4.3x closer to the truth** than the unregularised run, which took
  239.
- `info.stoppedBy` ("discrepancy" | "step" | "maxIter" | "zeroField")
  plus `residualRMS` and `noiseLevel` make the regularisation state
  inspectable. Without `NoiseLevel` a warning states that the result
  depends on `MaxIter`.

## [3.4.1] - 2026-08-11

### Documented
- The GravIS Greenland comparison closes to **1.6%** once the filter is
  declared correctly. Dahle et al. (2025, Sect. 2.2.1) specify the
  basin-average method: spectral masking, then a Wiener optimal filter
  approximately equivalent to a Gaussian of 4 deg latitude half-width,
  then conversion to surface mass, then least-squares adjustment. The
  forward model is filtered IDENTICALLY to the observations. The inputs
  are genuinely unfiltered Level-2B coefficients, so `Filter="none"`
  looked right and was not: declaring `Filter="gauss445"` moved the
  trend from -234.9 to -227.4 Gt/yr against a published -231.1, turning
  a 1.6% overshoot into a 1.6% undershoot. It also drifts less with
  continued iteration (2.0 Gt/yr over 300->900 steps against 2.5),
  because the filter regularises the inverse problem.
- The same paper specifies the mask: 1 until 200 km outside the
  grounding line, tapering to 0 by 600 km - `basinKernel(BufferKm=200,
  TaperKm=400)`, not a hard outline.
- **VDK**: it is DDK with the filter matrix rebuilt from each MONTH's
  formal error covariance, i.e. a per-epoch usage of `designFilter` with
  `Noise=`. What the toolbox lacks is the data, not the method - monthly
  full covariances are not distributed with Level-2 products, and for
  COST-G they are not available at all. Documented on the designFilter
  page with the two practical routes (GravIS VDK-filtered Level-2B
  products, or a SINEX normal-equation covariance). Note the GravIS ice
  basin averages do not use VDK at all.

## [3.4.0] - 2026-08-11

### Added
- `shSeries.removeAlias`: fits a tidal-alias harmonic (S2, 161 days by
  default) jointly with bias, trend, annual and semi-annual and
  subtracts only that harmonic - the correction GravIS applies to its
  Level-2B products. Across the GRACE/GRACE-FO boundary the nodal
  planes are not aligned, so a FIXED 100 degree phase offset applies to
  the later mission (Landerer et al. 2020). The offset adds no free
  parameters and ignoring it mis-estimates the amplitude badly
  (validated in `tools/dev/validate_alias.py`).
- Guide chapter "Validation against GravIS": the complete worked
  comparison against the published Greenland basin averages - the
  reference and why it must be recomputed over the span actually used,
  all four Level-2B corrections, the basin average, the results table,
  and what each step is worth.

### Measured
- Removing the S2 alias moves a twenty-year Greenland trend by
  **0.2 Gt/yr**, i.e. nothing: a 161-day harmonic is very nearly
  orthogonal to a linear trend over two decades. It remains a real
  correction for the monthly series. Recorded because it was expected to
  explain part of the residual gap and did not.

## [3.3.1] - 2026-08-11

### Fixed (documentation)
- `leakageCorrect`: the mask must cover EVERY region that can hold mass,
  not only the one being measured. A target-only mask tells the
  iteration that mass exists nowhere else, so neighbouring signal is
  forced into the target and the result is biased high. Validated
  against the published GravIS Greenland basin-average series (COST-G
  RL02, -231.1 Gt/yr over the matching span, all four GravIS
  corrections applied to the same Level-2 input):

  | mask | Gt/yr | vs published |
  |---|---|---|
  | none (naive grid integral) | -170.3 | -26% |
  | Greenland only | -259.8 | +12% |
  | union with Canadian Arctic, Iceland, Svalbard | -242.5 | +5% |

  The neighbours absorb -112 Gt/yr in the union solution, and it
  converges faster. No API change: `Mask=` already accepts any union -
  the error was in how it was used, which is exactly why it needed
  saying in the help rather than fixing in code.

## [3.3.0] - 2026-08-11

### Added
- `shLowLevel.readSHM`: the GRAVIS/GRACE SHM format - a YAML header
  followed by `GRCOF2` (a field) or `GRDOTA` (a rate in 1/yr) records,
  gzip transparent, GM and the reference radius taken from the header.
  GravIS Level-2B products - monthly solutions, mean fields, the
  ICE-6G_D (VM5a) GIA model - are a major public data source the
  toolbox could not read at all. A `GRDOTA` result feeds
  `standardChain(GIA=)` directly, which is what it always wanted.
  Python-validated first (`tools/dev/validate_shm.py`).
- Fixtures for both record types, and the GravIS Greenland drainage
  basin geometries, under CC-BY-4.0 with attribution in NOTICE.

## [3.2.2] - 2026-08-11

### Fixed
- `standardChain` stopped with "No TN-13 entry within 0.050 yr" whenever
  the correction tables trailed the solutions - which they always do,
  since TN-13/TN-14 are published weeks after the monthly fields. On a
  current 257-month ITSG series the two newest epochs were uncovered and
  the whole chain failed. Uncovered epochs are now dropped and the fact
  recorded in `rep.steps`; `OnMissing="error"` restores the old
  behaviour and `Tolerance=` is settable. Dropping beats the
  alternative of keeping them, which would splice UNCORRECTED months
  onto corrected ones - a step in the series exactly where people look
  for the newest signal.
- `shSeries.read` / `fromFolder` treated `writeGFC`'s
  `<file>.gfc.provenance.json` sidecars as solutions: the read pattern
  has to be loose enough for `.gfc.gz`, so `*.gfc*` matched them and a
  folder written BY the toolbox could not be read back BY the toolbox.
- `shCoefficients.write` gained `Sidecar=`, which only the low-level
  `writeGFC` had, so class-API users could not turn the sidecar off.

### Notes
- All three were found by running the toolbox against a full 24-year
  ITSG series (257 monthly n96 solutions, 2002-2026) rather than against
  fixtures. The suite is fixture-bound by necessity; the failures it
  cannot see are the ones that need real data of real length.

## [3.2.1] - 2026-08-11

### Added
- `tests/testScience.m`: regression against values published OUTSIDE the
  toolbox - WGS84 and GRS80 J2 from the defining constants, the Gegout97
  (PREM) load Love numbers including the CM-frame `k'_1 = 0`, the
  closed-form EWH kernel of Wahr et al. (1998) evaluated with the
  shipped Love numbers, the GRACE-minus-SLR C20 discrepancy, TN-13
  geocenter amplitudes for all three providers, and the Gaussian
  half-weight relation. The existing suites check that the code does
  what the code intends, and stay green if a formula is consistently
  wrong; these have somewhere to fail. Verified that the kernel test
  catches both a missing `rho_water` (factor 1000) and a `k_n` sign
  error.
- The trend regression is opt-in via `SHX_SERIES_FOLDER`, because a
  trend needs a monthly series that is far too large to ship. It
  deliberately does NOT read the persistent data folder, so a routine
  acceptance run never touches a network or archive drive.

### Fixed
- Test cleanups no longer throw from `onCleanup` destructors. `d =
  tempname` invents a NAME, so a test whose code path did not reach the
  `mkdir` left `rmdir` with nothing to remove: every CI log carried a
  destructor warning that MATLAB cannot attach to any test, and which
  can mask the real failure that ended the test early. All 30 sites now
  use `tests/rmIfFolder.m`, and one genuinely dead variable in
  `testRobustness` (a `tempname` that was never used for anything) is
  gone.

## [3.2.0] - 2026-08-11

### Added
- `shLowLevel.leakageCorrect`: iterative forward modelling of filter
  leakage and attenuation. Finds the mass field whose filtered image
  matches the observation, using only the filter - no external signal
  model. `Mask=` confines the solution to where mass can exist, which is
  both better conditioned and the variant that removes leakage rather
  than redistributing it: on a synthetic disc with an exact mask it
  recovers 0.9998 of the truth and leaves exactly zero outside.
- `shLowLevel.gridScaling`: per-pixel scaling factors from a model
  series pushed through the same chain. Amplitude-invariant to 1e-15 -
  the factors describe the model's PATTERN, not its scale.
- Both implemented and validated in Python
  (`tools/dev/validate_leakage.py`, a compact numpy SH port whose
  analysis/synthesis round trip is exact to 4e-16 and whose Gaussian
  weights match `shGaussianWeights` to 2.7e-11) BEFORE the MATLAB port.
  The validation caught two design errors before any MATLAB was written:
  a scaling test that was exact by construction and therefore validated
  nothing, and a test cap so large the filter barely touched it.

### Two contracts that are easy to get wrong
- Convergence is judged on the CHANGE OF THE SOLUTION, not the residual.
  A masked problem is inconsistent - no field confined to the region
  reproduces an observation with energy outside it - so the residual
  plateaus at a nonzero floor while the solution is converged. The first
  implementation stopped on the residual and reported failure on exactly
  the case the mask exists for.
- `gridScaling` returns NaN where the model carries no signal instead of
  a ratio of two numerical zeros, and `info.kMedian` summarises the
  pixels the model gives mass to. Over all finite pixels the median is
  dominated by pure-leakage pixels, whose `k = 0` is correct but is not
  the basin scaling factor anyone is looking for.

## [3.1.4] - 2026-08-11

### Added
- `shLowLevel.safeMove`: the `.part` -> target swap every fetcher ends
  with, now with verification and exponential-backoff retries
  (0.25/0.5/1/2/4 s). Windows antivirus and cloud-sync clients open a
  just-written file to scan it and hold it briefly; a plain `movefile`
  hits that window and reports a download failure for a file that
  downloaded perfectly. Each attempt checks that the destination exists
  and the source is gone rather than trusting the status flag, because
  synced and network volumes have been observed reporting success early.
  All six call sites in fetchICGEM (2), fetchITSG, fetchTN,
  fetchLoveNumbers and fetchDDK go through it.
- `shLowLevel:safeMove:locked` names the usual cause and the two fixes
  (exclude the data folder from on-access scanning, or move it off the
  synced drive) instead of surfacing a bare OS message.

## [3.1.3] - 2026-08-11

### Added
- `tools/doc_sync_audit.py`, a sixth quality gate, wired into CI before
  the MATLAB setup step. The five existing gates were all green while
  the API reference listed 5 of 12 `fetchICGEM` options, eleven help
  pages were unreachable from the Help browser, the v3.1.1 tag reported
  3.1.0 and the guide advertised a call that threw - none of that was
  detectable automatically. The gate checks doc snippets against the
  parsed contracts (including option VALUE types, which is the class a
  name-only check misses), API-reference freshness, `helptoc.xml`
  reachability, version consistency, help placeholders and defaults, and
  a narrative-coverage floor. Every check is verified by deliberately
  breaking the repository and confirming it fires.
- `tools/dev/README.md`: gate order, regeneration order, the SIGPIPE
  trap, and the guide-asset recovery note.

### Fixed
- `fetchITSG` help claimed `Nmax (96)` unconditionally. The default is
  NaN, which resolves to 96 for monthly and 40 for daily solutions, and
  120 is also accepted for monthly - the help listed neither.
- `seaLevelFingerprint` documented `Tol` as `1e-8 * |eustatic|`, mixing
  the option's value with its effect; the default is 1e-8 and the
  tolerance is relative.
- `vceRescale` help wrote `round(2/3*Lmax)` where the arguments block
  says `round(2/3 * idx.Lmax)`.
- `shCoefficients.mtimes`/`uminus` and `shClimatology.withNote`/
  `fromCoef` were documented only in the generated API reference.

## [3.1.2] - 2026-08-11

### Added
- Narrative documentation for the v3.x features that existed only in the
  generated API reference: new help pages `standardChain.html`,
  `designFilter.html` and `icgemSeries.html`, plus sections for
  `poleTideConvert` (referenceSystems), `synthesisPoints` (shSynthesis),
  `listITSG` (fetchITSG) and `fetchLoveNumbers` (readLoveNumbers). All
  reachable from `helptoc.xml` (46/46 pages).
- Workflow guide Edition 5: the ICGEM time-series workflow and why
  `Mode="auto"` fetches one archive rather than 300 files, the standard
  chain and why its order is not negotiable, custom filter design, and a
  section on what real provider files actually contain (FORTRAN
  D-exponents, the ICGEM 2.0 column order, ragged record groups).

### Fixed
- The guide's figures were read from a container path that does not
  survive a session, so the PDF was not reproducible and `make_figs.py`
  could only regenerate 3 of the 10. The assets were recovered from the
  shipped PDF and now live in `tools/dev/guide_assets/` in the repo.
- 58 name-value options across 23 entities were "documented" as `see
  arguments block`: `help` showed nothing while the generated API
  reference showed size, type and default. All replaced with real
  descriptions; pinned to zero by `testHelpHasNoPlaceholders`.
- Four help defaults disagreed textually with their arguments blocks
  (`synthesisMatrix.NLat/NLon`, `tvANSFilter.VCEMinDegree`,
  `vceRescale.MinDegree`).
- `fetchICGEM` Outputs help never mentioned `info.mode`, added in the
  series work, nor that bulk selections concatenate their file lists.

## [3.1.1] - 2026-08-11

### Added
- `fetchICGEM` fetches TIME SERIES: temporal-catalogue rows (from
  `listICGEM(Type="temporal")`) download a whole monthly series into
  `<dataFolder>/icgem/series/<group>_<center>_<series>/` - per file
  (resumable, throttled, each file verified before swap; `Files=`
  filter, `FileList=` injection for subsets/mirrors). DEFAULT is the
  server's whole-series ZIP in ONE request (`Mode="auto"`; the rate
  limiter punishes hundreds of per-file requests), with automatic
  per-file fallback; `Mode="archive"|"files"` force either. `Type=`
  mirrors listICGEM so numeric/name selection works on either
  catalogue; the result folder feeds `shSeries.fromFolder` and
  `standardChain` directly.

### Fixed
- `fetchICGEM` bulk selection of temporal series: the multi-selection
  loop assumed one file per row and crashed assigning a 283-file
  series into `strings(1,K)` - AFTER the download had succeeded. File
  outputs are now collected per selection and concatenated; `info`
  carries a uniform `mode` field ("model"|"archive"|"files") across
  all paths; the summary reports selections and files separately.
- `shSeries.filter("tvANS", ...)` accepts and forwards `Shrinkage`,
  `VCEMinDegree` and `VCEBands`. The workflow guide has advertised
  `ts.filter("tvANS", Blocks="auto", VCEBands=[0 16 33 61])` since
  Edition 2; the class method rejected it (`MATLAB:TooManyInputs`), so
  banded VCE was reachable only through the low-level entry point. The
  defaults stay in `shLowLevel.tvANSFilter` - `[]` means "do not
  override" - so there is exactly one home for each default.
- Documentation metadata brought back in sync with the code: `Contents.m`
  is at 3.1.1 (the v3.1.1 tag shipped a 3.1.0 version line, so
  `ver('shAnalysis')` under-reported), `CITATION.cff` follows, and the
  `[3.1.0]`/`[3.0.1]`/`[2.7.0]` sections are dated instead of
  "Unreleased". `Contents.m` no longer advertises the `compat/` folder
  and the legacy v1 suite, both removed in v3.0.0, and its first line is
  now a short product name (MATLAB's `ver` printed a truncated sentence).
- `apiReference.html` regenerated: it predated the series work and listed
  5 of the 12 `fetchICGEM` options. `make_apiref.py` now parses the
  version from `Contents.m` instead of carrying a hardcoded string.
- `helptoc.xml` reaches all 43 help pages (11 were unreachable from the
  Help browser) and its labels no longer name superseded versions;
  `shAnalysis.html` is a v3 page again instead of describing `compat/`.
- `shCoefficients.disp`, `shSeries.disp` and `shClimatology.disp` have
  help text; `help_audit.py` no longer exempts `disp`.
- `shReadGFC` bulk parser survives ragged record groups: EIGEN-5S/5C
  carry a single gfc line with a trailing epoch among thousands of
  uniform ones, which silently misaligned the rectangular sscanf and
  requested a (date x date) coefficient matrix. Groups are now checked
  for rectangularity (rows out == lines in), split by numeric width
  when ragged, and guarded by an n/m integer/bound sanity net; dirty
  groups still fall back to the line parser.

## [3.1.0] - 2026-08-10

### Fixed
- `shReadGFC` ICGEM 2.0 acos/asin column order corrected to the
  real-world layout `... sigC sigS t0 t1 period` (verified against
  CNES_GRGS.RL05MF). The previous parser assumed period-first and
  MIS-PARSED such files (period = a date code, t1 = 1.0) - a
  correctness bug independent of, and older than, the speed problem.
- `shReadGFC` parses variable-term files (GRGS mean fields, ICGEM 2.0)
  in bulk too: group-wise sscanf per record key, vectorized epoch
  conversion, vectorized variableTerms assembly (the per-element
  struct growth was O(N^2)). `icgemDate2Year` is pure arithmetic now -
  the datetime/calendarDuration construction cost ~1 ms per call and
  was invoked per record. A 73 MB CNES_GRGS mean field verification
  went from tens of minutes (the second reported "stall") to seconds;
  dirty records still fall back to the authoritative line parser,
  which itself now assembles variableTerms vectorized.
- `shReadGFC` fast path accepts FORTRAN D-exponents (EIGEN-6C4 switches
  from e to D serialization at degree 371): normalized to E after the
  key strip, matching str2double semantics. Previously such files fell
  back to the legacy line parser - minutes for a 178 MB model, which
  presented as fetchICGEM "stalling" during download verification.
  CI exposed the deeper bug: str2double returns NaN for D-exponents,
  so the legacy parser had silently corrupted such files to NaN above
  the switch degree - both paths now normalize and agree.
- `fetchICGEM` bulk mode ("all", index vectors) no longer stalls after
  the first file: the ICGEM server rate-limits rapid requests (HTTP
  429 / tarpit). Bulk fetching now throttles between models (`Pause=`,
  3 s), retries 429/timeouts with 30/60/120 s backoff (`Retries=`),
  prints a [k/K] counter with per-file size and timing, continues past
  failed models with a summary instead of aborting, and honours
  `Quiet=`. Interrupted runs resume (present files are skipped).

### Changed
- `shReadGFC` is ~100x faster on large static files: bulk parsing
  replaces the per-line fgetl/str2double loop (28.6 s -> 0.21 s at
  n720 measured on PCWIN64; EGM2008-class n2190 in 1.7 s instead of
  minutes). Variable-term files (gfct/trnd/dot/acos/asin) keep the
  proven line parser; both paths verified struct-identical on all
  fixtures and guarded by an equivalence test.

### Added
- `shLowLevel.standardChain`: the canonical GRACE post-processing chain
  (read folder -> TN-14 C20/C30 -> degree-1 -> optional GIA trend ->
  Gaussian/DDK/custom-W filter) as one correctly ordered call with a
  step-by-step provenance report; explicit-file overrides for
  reproducible pipelines.
- `shLowLevel.designFilter`: DDK-class anisotropic filter
  W = (N + Alpha*inv(S))^-1 N per order block, from formal sigmas +
  Kaula rule or full covariances (buildNoiseCov/buildSignalCov);
  output in the readDDK block format (drops into applyDDK and the
  class methods unchanged). Python-validated with pinned references.
- Workflow guide version stamp is generated from Contents.m (the PDF
  title claimed v2.5 since the beginning).

## [3.0.1] - 2026-08-10

### Added
- `tools/help_audit.py`: documentation completeness gate, run in CI
  before the test suite. Every public function and method must document
  every positional input, every option WITH its default, and every
  output with size/type annotations, plus an example. 328 findings
  fixed across the tree to reach zero; the gate keeps it there.

## [3.0.0] - 2026-08-10 - DOI: 10.5281/zenodo.21876348

### BREAKING
- Namespace renamed `+shx` -> `+shLowLevel`: every call site changes
  from `shx.fun` to `shLowLevel.fun`; error identifiers from `shx:*`
  to `shLowLevel:*`. Mechanical migration: replace `shx.` and `shx:`.
- `compat/` (v1 wrappers) and the legacy v1 suite are removed for good.

### Added
- `shLowLevel.listITSG`: numbered catalogue of the TU Graz server
  (releases x monthly n60/n96/n120 x daily), local-mirror mode.
- `fetchITSG`: `Release=` override (no more silent 2018/operational
  mixing), `months = "all"`, `Catalog=` selection by catalogue row,
  `BaseURL=` mirror, n120.
- `listICGEM` numbered; `fetchICGEM` accepts row indices, name vectors,
  and `"all"`.
- `setup_shAnalysis`: `DataFolder=` (applied before any fetcher) and
  `FetchITSG = "all"`.
- `shLowLevel.fetchLoveNumbers`: load/deformation Love numbers from the
  GROOPS collection (TU Graz), safe-swap + parser.
- `shLowLevel.synthesisPoints`: pointwise potential/disturbance/anomaly
  at arbitrary (lat, lon, r) incl. upward continuation; Python-pinned.
- Generated API reference page `html/apiReference.html` (inputs,
  options with defaults, outputs for every public function/method).

## [2.7.0] - 2026-08-10

### Added
- `shLowLevel.poleTideConvert`: C21/S21 conversion between the IERS2010 and
  IERS2018 (secular pole) mean-pole conventions, solid + ocean pole
  tide, overridable coefficients; raw (C, S, epoch) and object forms.
  Python-validated (identities exact, dS21 trend to < 1% of the
  dominant term).
- `Proxy=` option in all four fetchers and `setup_shAnalysis`:
  per-call proxy via matlab.net.http (base MATLAB) for institutional
  networks; websave path continues to honour MATLAB Web Preferences.
- Provenance JSON sidecars (`<file>.provenance.json`: tool version,
  MATLAB release, metadata) on `writeGFC`, `writeGrid`,
  `writeAnimation`; `Sidecar = false` disables.
- `writeGrid` netCDF output is CF-1.8 complete (standard_name, axis,
  coordinates, long_name, Conventions/title/source/history globals);
  the source stamp is dynamic (was hardcoded "v2.4").

## [2.6.0] - 2026-08-10 - DOI: 10.5281/zenodo.21871299

### Added
- Comparison suite: `shLowLevel.compareSolutions` / `shLowLevel.compareSeries`
  aggregators (also `g.compare` / `ts.compare`) over new standalone
  metrics `diffSpectrum`, `spatialStats`, `nashSutcliffe`,
  `effectiveCorr` (AR(1)-corrected significance), `threeCorneredHat`
  (N >= 3), and `shLowLevel.taylorDiagram`. All numerics Python-validated.
- GitHub Actions CI (`ci/shanalysis-ci.yml`); full suite green on
  Linux and Windows.

### Fixed
- `writeAnimation`: platform-aware VideoWriter profile (.avi portable;
  MPEG-4 unavailable on Linux now errors clearly).
- `shDegreeRMS` help documented fields that did not exist
  (`n`/`amp` -> `degree`/`degAmplitude`/...).
- Demo provenance line now generated from `shLowLevel.version`.

## [2.5.1] - 2026-08-08

### Added
- `Update=` option in all four fetchers (safe, parse-verified swap;
  TN-13/TN-14 grow monthly upstream); `BaseURL=` mirror-folder mode in
  `fetchTN`; forwarded through `setup_shAnalysis(Update=true)`.
- `shLowLevel.version`: toolbox metadata parsed from `Contents.m`.

### Fixed
- `applyDDK` truncates Lmax-120 filter blocks to the field nmax
  (standard n60/n96 use; previously errored).
- `shLowLevel.pctile` accepts vector percentiles.
- `setCoefficient` initializes both sigma stacks (pairing invariant);
  a lone sigmaC used to break shSeries stacking.

## [2.5.0] - 2026-08-07

### Added
- Exact constrained tvANS posterior sigmas (the pre-2.5 formula
  UNDERESTIMATED; documented); basin sigmas include the deterministic
  model contribution.
- Kendall AR(1) correction in `fitDeterministicModel`.
- `shLowLevel.readLoveNumbers` (layouts, frames CE/CF/CM, pchip in log(1+n)).
- `setup_shAnalysis` one-call setup; `shLowLevel.fetchTN`.
- Generated 69-page API reference (Part IV) with executable examples
  (`docs/apiExamples.json`, `testAPIExamplesRun`).
- Real-file validation: CSR/JPL/GFZ TN-13, GSFC TN-14, DDK Wbd
  binaries, ITSG monthly/daily, 460 MB operational SINEX streaming.

## [2.4.2] - earlier

- Full-sine-wing triangle difference plots; robust cache-speedup test;
  DDK reader; ICGEM catalogue; persistent data folder.
