# Changelog

All notable changes to shAnalysis. The version line in `Contents.m`
(read by `shLowLevel.version` and MATLAB's `ver`) is the single source of truth.

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
