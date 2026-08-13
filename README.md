# shAnalysis

![CI](https://github.com/matweigelt/shAnalysis/actions/workflows/ci.yml/badge.svg)
![MATLAB](https://img.shields.io/badge/MATLAB-base%20only%2C%20R2026a%20tested-orange)
![Platforms](https://img.shields.io/badge/platforms-Windows%20%7C%20Linux-lightgrey)
![Latest tag](https://img.shields.io/github/v/tag/matweigelt/shAnalysis?label=version)
![License](https://img.shields.io/github/license/matweigelt/shAnalysis)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21871298.svg)](https://doi.org/10.5281/zenodo.21871298)

Class-based MATLAB toolbox for spherical harmonic (Stokes coefficient)
analysis of GRACE/GRACE-FO, GOCE, and static gravity field models.
Three user-facing classes cover the whole chain - reading, spectra,
filtering, synthesis, climatology, basin averages - and every numerical
building block is callable directly from the `+shLowLevel` package.

Requires MATLAB R2021a or newer (`arguments` name=value syntax) and
**nothing beyond base MATLAB**: no toolboxes.

Validated against the published GravIS Greenland ice-mass series: 2.8%
on the mass trend over 2002-04..2023-02 (217 monthly solutions - the
span of the auxiliary correction tables at the time of the
reproduction), matching their documented processing step for step
(see the guide's *Validation against GravIS* chapter). The span dates
matter: the reference trend moves by ~10 Gt/yr over two years of end
date, so the percentage is only meaningful together with them.
An independent 2026-08 audit reproduced the chain end-to-end from
public data (headline -225.7 vs -224.6 Gt/yr; the wrong-noise control
-201.0 exactly). The same chain reproduces the GravIS terrestrial
water storage product over eleven major river basins to 0.032 cm/yr
trend RMS and a 1.001 median amplitude ratio (DDK3 + ICE-6G_D GIA,
full 252-month span), and locates the Antarctic one-map-inversion
limit at 14% against the AWI joint-basin product - with a 0.905
difference-pattern correlation to the official TU Dresden kernel grid,
i.e. a method-class distance, not a chain defect (guide chapters
V7/V8). All three validations are available as one-call methods with
exchangeable inputs: `shLowLevel.greenlandChain`,
`shLowLevel.antarcticaChain` and `shLowLevel.twsChain`, built on the
shared `shLowLevel.gravisL2B` correction core. All downloads retry
politely on 429/5xx honouring the server's Retry-After header
(`shLowLevel.httpFetch`). Data locations persist across sessions via
`setup_shAnalysis(SeriesFolder=..., GravisFolder=..., DDKFolder=...)`
(setpref-backed); `quantity="none"` synthesizes raw coefficient
fields. `shLowLevel.oceanChain` closes the chain family with the
barystatic ocean series (measured GIA lever +0.89 mm/yr), and
`slepianProject` turns the Slepian basis into well-posed regional
analysis. `shLowLevel.fetchGAX` completes the ocean story:
GAD/GAA from 2002 on, one call, full Chambers-&-Willis restoration in
`oceanChain`. The ocean lane now forks: `obpChain` keeps the air
column and delivers GravIS-style bottom-pressure fields, while
`eofSeparate` splits residual circulation from noise (North rule) so
`oceanChain`'s noise proxy stops blaming the ocean for being dynamic. The fetch family is complete and consistently named:
series via `fetchITSG`/`fetchICGEM`, products via `fetchGAX`,
`fetchSINEX` (ITSG monthly normals - the only public per-month SINEX,
~460 MB each) and `fetchITSGBackground`. Those normals feed `vdkApply` - the
VDK/VADER decorrelation of Horvath et al. (2018) with true monthly
covariance structure; `run_vdk_series.m` drives the full-series batch. And the hydrology story
closes the loop from field to warning: `hydroExtremeIndex` delivers
GRACE-DSI drought categories and the Reager flood-predisposition
deficit for every cell or basin.

## Documentation

The complete **[workflow & theory guide (PDF)](docs/shAnalysis_workflow_guide.pdf)**
covers the processing chain, filter theory, uncertainty propagation, worked
recipes, and a full chapter reproducing the published GravIS Greenland
series step by step. The generated API reference lists typed inputs and
outputs with an executable example for every public function. HTML help
for MATLAB's `doc` browser ships in `html/`, and
`shLowLevel.makeTutorials` writes Live Scripts you can open and run.

## Layout

```
shAnalysis/
  shCoefficients.m   one coefficient set (value class, immutable)
  shSeries.m         time series of coefficient sets
  shClimatology.m    bias/trend/annual/semi-annual model
  runAllTests.m      full validation suite
  +shLowLevel/       numerical internals (stable, package-qualified)
  demo_shAnalysis.m  16 selectable demonstration cases
  tests/             correctness / contract / robustness / science /
                     performance
  tests/test_data/   fixtures incl. two real ITSG monthly solutions and
                     GravIS Level-2B samples (CC-BY-4.0, see NOTICE)
  html/              doc pages for MATLAB's `doc` browser
  tools/             CI gates; tools/dev/ the doc generators (see its README)
```

## Quick start

```matlab
cd shAnalysis
setup_shAnalysis(Permanent = true, Download = "core")   % path + TN files
shLowLevel.makeTutorials(Cases = "core");   % Live Scripts to open and run
% later, monthly refresh of the growing TN files (safe, parse-verified swap):
setup_shAnalysis(Download = "core", Update = true)      % or shLowLevel.fetchTN(Update = true)
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

## Opt-in test data

Two checks need data too large to ship. Point `setup_shAnalysis` at yours:

```matlab
setup_shAnalysis(SeriesFolder = "D:/grace/itsg", ...
                 MasconFile = "D:/grace/GRCTellus_JPL.nc");
```

| Variable | Enables | Data |
|---|---|---|
| `SHX_SERIES_FOLDER` | the trend regression in `testScience` | monthly `.gfc` files, > 24 epochs, GSM level |
| `SHX_MASCON_FILE` | `readMascon` against a real product | CSR (no login), GSFC, or JPL PO.DAAC |

Neither is required. Without them both tests report as filtered rather
than passing silently. They exist because the fixture suite cannot see
the class of failure that only appears on real data of real length -
three bugs surfaced the first time the toolbox ran against a full
24-year series.

## Running the tests

```matlab
results = runAllTests;        % writes tests/runAllTests_results.txt
```
The complete console output plus a per-test table always land in
`tests/runAllTests_results.txt` (fixed name, overwritten) - upload that
file for review; pasted console text gets truncated.

Five suites run: numerical correctness, API contract, robustness,
**science** (regression against values published outside this toolbox -
defining constants, load Love numbers, the closed-form EWH kernel,
geocentre amplitudes) and performance. Set `SHX_SERIES_FOLDER` to a
folder of monthly `.gfc` files to enable the opt-in trend regression;
the fixture suite cannot see the class of failure that only appears on
a real series of real length.

## What is new

The **[CHANGELOG](CHANGELOG.md)** carries the full history with the
reasoning behind each change. The headlines since v3.0.0:

**Leakage correction (v3.2.0+).** `shLowLevel.leakageCorrect` recovers the
mass field whose filtered image matches the observation;
`shLowLevel.gridScaling` gives per-pixel scaling factors from a model
series; `shLowLevel.sensitivityKernel` (v3.7.0) trades leakage against
propagated noise *explicitly* - closed form, no iteration, no stopping
criterion.

Two rules that cost real accuracy when broken, both measured:

* **The mask must cover every region that can hold mass**, not only the
  one you are measuring. A Greenland-only mask forces the Canadian
  Arctic and Svalbard into Greenland: +12% instead of +4%.
* **Declare the filter the forward model sees**, not the one your input
  file carries. GravIS ice basin averages start from *unfiltered*
  coefficients and still apply a Wiener filter inside the inversion;
  declaring `"none"` instead of `"gauss445"` moved the trend by 3%.

**Stopping is regularisation (v3.5.0+).** `leakageCorrect` solves an
ill-posed problem and semiconverges: the error against the truth falls,
then *rises*, while the residual keeps shrinking. Pass `NoiseLevel` -
from `shLowLevel.oceanRMS`, the open-ocean noise metric - and the run
stops on the discrepancy principle. `info.stoppedBy` says whether the
answer was regularised at all.

**Self-consistent geocentre (v3.6.0).** `shLowLevel.estimateDegree1`
derives C10/C11/S11 from the data plus an ocean model (Swenson et al.
2008), so you no longer depend on somebody having published a TN-13
series for your exact Level-2 product. Output layout identical to
`readTN13`.

**GRAVIS Level-2B (v3.3.0).** `shLowLevel.readSHM` reads the SHM format
(`GRCOF2` fields, `GRDOTA` rate models) used by GravIS mean fields, GIA
models and monthly solutions - gzip transparent, GM/R from the header.

**The standard chain (v3.1.0).** `shLowLevel.standardChain` runs
read -> TN-14 -> degree 1 -> GIA -> filter in the one correct order and
returns a provenance report. Correction tables always trail the
solutions, so uncovered epochs are dropped and recorded rather than
raising an error.

**Documentation is gated (v3.1.2+).** Six CI gates run before the test
suite; `tools/doc_sync_audit.py` checks the in-file help, `html/` and
the workflow guide against the source and against each other - every
documented call is contract-checked, including option value *types*.

## Data sources

| What | Where | Toolbox route |
|---|---|---|
| ITSG monthly / daily solutions | ftp.tugraz.at | `shLowLevel.listITSG`, `fetchITSG` |
| ICGEM static models and time series | icgem.gfz.de | `shLowLevel.listICGEM`, `fetchICGEM` |
| TN-13 (degree 1), TN-14 (C20/C30) | GFZ ISDC | `shLowLevel.fetchTN` |
| DDK filter matrices | upstream repository | `shLowLevel.fetchDDK` |
| Load Love numbers (GROOPS set) | ftp.tugraz.at | `shLowLevel.fetchLoveNumbers` |
| GravIS Level-2B and Level-3 | gravis.gfz.de, ISDC | `shLowLevel.readSHM` (manual download) |

**On COST-G covariances.** COST-G publishes its combined *solutions*
through ICGEM, and combines on the normal-equation level internally, but
does not distribute per-month SINEX normal equations or covariance
matrices as a public product - Dahle et al. (2025) note the
variance-covariance information is not available for the combined
solutions. There is therefore no COST-G SINEX endpoint to fetch. For a
full covariance use **ITSG**, which does publish monthly SINEX normal
equations (`shLowLevel.readSINEX`), or the GravIS **VDK-filtered**
Level-2B products if what you want is the filtered result rather than
the covariance itself.

## History

Releases before v3.0.0 are described in the CHANGELOG. v3.0.0 was a
breaking change: the package was renamed `+shx` -> `+shLowLevel`, and the
`compat/` folder with the v1 function API was removed together with the
v1 legacy test suite. The v1 function names live on inside the package -
call `shLowLevel.shReadGFC(...)` instead of `shReadGFC(...)`, or use the
classes.

## Citing

See [CITATION.cff](CITATION.cff). Bundled third-party fixtures and their
licences are listed in [NOTICE](NOTICE).

Developed by Matthias Weigelt with the help of Claude.
