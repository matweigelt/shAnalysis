"""outdocs_inject.py - author + inject "%   Outputs" blocks into help text.

For every public function/method whose help lacks an Outputs section,
insert a house-style block (name, size, type, description) right before
the attribution line of that help block. Idempotent: entries whose help
already contains an Outputs heading are skipped.

Developed by Matthias Weigelt with the help of Claude (Fable 5),
2026-08-07 (v2.5).
"""
import re

ROOT = "/home/claude/shx_build/shAnalysis"
n1 = "(nmax+1 x nmax+1)"

OUTSPEC = {
# ------------------------------------------------------------------ shx
"shx.applyDDK": [
    ("Cf", n1 + " double   DDK-filtered cosine coefficients"),
    ("Sf", n1 + " double   DDK-filtered sine coefficients")],
"shx.basinDeconvolve": [
    ("avg", "(K x T) double   deconvolved basin averages"),
    ("out", "struct: c (K x T) deconvolved coefficients, sigma (K x T) "
            "1-sigma incl. deterministic term (v2.5), avgNaive, attn, "
            "condA")],
"shx.buildNoiseCov": [
    ("N", "(P x P) double or block struct   noise covariance (order/"
          "parity block-diagonal for Assemble='blocks')"),
    ("info", "struct: mode, shrinkage, blocks metadata")],
"shx.buildSignalCov": [
    ("S", "(P x P) double   signal covariance in idx ordering"),
    ("info", "struct: mode, iterations, degree-variance model")],
"shx.csFromVec": [
    ("C", n1 + " double   cosine coefficients scattered from x"),
    ("S", n1 + " double   sine coefficients scattered from x")],
"shx.dataFolder": [
    ("folder", "(1 x 1) string   current data folder (created on "
               "demand)")],
"shx.ddkNames": [
    ("names", "(1 x 8) string   Wbd file names for DDK1..DDK8")],
"shx.icgemDate2Year": [
    ("t", "(same size as input) double   decimal years")],
"shx.mcPropagate": [
    ("out", "struct: mean, sigma, samples (N x Q), q16/q50/q84 "
            "percentiles of the propagated functional")],
"shx.opApply": [
    ("xf", "(P x cols(x)) double   W_t * x (or W_t' * x in 'transp' "
           "mode) without forming the P x P matrix")],
"shx.parseGraceFilename": [
    ("meta", "struct: center, product, release, epochStart/Stop "
             "(decimal years), raw tokens")],
"shx.plotCovariance": [
    ("h", "(1 x 1) graphics handle   image of the covariance/"
          "correlation structure")],
"shx.plotSHCoeffTriangle": [
    ("h", "(1 x 1) graphics handle   butterfly triangle image (abs "
          "log10 or signed diff with RefC/RefS)")],
"shx.plotSHSpectrum": [
    ("h", "(1 x 1) graphics handle   log-log spectrum axes")],
"shx.readDDK": [
    ("W", "struct: nmax (1 x 1), name, blocks (1 x 241) struct array "
          "with m, cs, n (degrees), M (square block matrix)")],
"shx.readMascon": [
    ("mas", "struct: lat (I x 1), lon (J x 1), t (T x 1 decimal "
            "years), ewh (I x J x T) [m], units/meta from the netCDF")],
"shx.readSINEX": [
    ("snx", "struct: idx (shIndex-compatible), est (P x 1), M (P x P "
            "covariance or normal-equation matrix per Output=), "
            "sigma (P x 1), header meta")],
"shx.readTN14": [
    ("tn", "struct: epoch/epochStart/epochStop (M x 1 decimal years), "
           "C20, sigmaC20, C30, sigmaC30 (M x 1; C30 NaN before "
           "2012.16)")],
"shx.resolutionMap": [
    ("res", "(nlat x nlon) double   local effective resolution [km] of "
            "the stored filter operator")],
"shx.shDegreeRMS": [
    ("spec", "struct: n (nmax+1 x 1), amp/rms/var (nmax+1 x 1) per "
             "Quantity, err (nmax+1 x 1) when sigmas present")],
"shx.shDestripe": [
    ("Cf", n1 + " double   destriped cosine coefficients"),
    ("Sf", n1 + " double   destriped sine coefficients")],
"shx.shEvalGFCT": [
    ("C", n1 + " double   coefficients evaluated at the epoch"),
    ("S", n1 + " double   sine coefficients at the epoch")],
"shx.shFanFilter": [
    ("Cf", n1 + " double   fan-filtered cosine coefficients"),
    ("Sf", n1 + " double   fan-filtered sine coefficients"),
    ("Wnm", n1 + " double   separable degree x order gain")],
"shx.shGaussianFilter": [
    ("Cf", n1 + " double   smoothed cosine coefficients"),
    ("Sf", n1 + " double   smoothed sine coefficients"),
    ("Wn", "(nmax+1 x 1) double   Jekeli degree weights applied")],
"shx.shGaussianWeights": [
    ("Wn", "(nmax+1 x 1) double   Jekeli (1981) weights, W(1) = 1, "
           "monotone non-increasing")],
"shx.shIndex": [
    ("idx", "struct: n/m/cs (P x 1), P (1 x 1), nmax, MinDegree - the "
            "canonical vector ordering of the toolbox")],
"shx.shReadGFC": [
    ("C", n1 + " double   cosine coefficients (static part)"),
    ("S", n1 + " double   sine coefficients"),
    ("header", "struct: GM, R, nmax, tide system, errors flag, "
               "variableTerms for gfct files")],
"shx.shSpectralCrossover": [
    ("nc", "(1 x 1) double   first degree where the error spectrum "
           "exceeds the signal (NaN if none)")],
"shx.shSynthesis": [
    ("grid", "(nlat x nlon) double   synthesized field in the "
             "requested quantity units"),
    ("lat", "(1 x nlat) double   geocentric latitudes used"),
    ("lon", "(1 x nlon) double   longitudes used"),
    ("P", "struct   reusable Legendre/cache handle for repeat calls")],
"shx.shSynthesisDeformation": [
    ("up", "(nlat x nlon | npts) double   vertical elastic deformation "
           "[m]"),
    ("north", "same size   horizontal north component [m]"),
    ("east", "same size   horizontal east component [m]")],
"shx.vceRescale": [
    ("s", "(T x 1) or (nBands x T) double   monthly (band-wise) VCE "
          "noise variance factors")],
"shx.vecFromCS": [
    ("x", "(P x 1) double   coefficients gathered in idx ordering")],
"demo_shAnalysis": [
    ("reg", "(D x 3) table   demo registry: id, title, exercised API")],
# --------------------------------------------------- shCoefficients
"shCoefficients.shCoefficients": [
    ("obj", "(1 x 1) shCoefficients   immutable coefficient set with "
            "GM, R, epoch, sigmas, history")],
"shCoefficients.plus": [
    ("out", "(1 x 1) shCoefficients   sum; sigmas RSS; epochs must "
            "match")],
"shCoefficients.minus": [
    ("out", "(1 x 1) shCoefficients   difference; sigmas RSS")],
"shCoefficients.uminus": [
    ("out", "(1 x 1) shCoefficients   negated coefficients")],
"shCoefficients.mtimes": [
    ("out", "(1 x 1) shCoefficients   scalar-scaled; sigmas scale by "
            "|a|")],
"shCoefficients.times": [
    ("out", "(1 x 1) shCoefficients   elementwise-scaled coefficients")],
"shCoefficients.truncate": [
    ("out", "(1 x 1) shCoefficients   truncated to the new nmax")],
"shCoefficients.setCoefficient": [
    ("out", "(1 x 1) shCoefficients   copy with C/S(n+1, m+1) "
            "replaced")],
"shCoefficients.applyTN14": [
    ("out", "(1 x 1) shCoefficients   C20 (and C30 when available) "
            "replaced by the SLR values, sigmas updated")],
"shCoefficients.addDegree1": [
    ("out", "(1 x 1) shCoefficients   degree-1 row set from the TN-13 "
            "record nearest to obj.epoch")],
"shCoefficients.destripe": [
    ("out", "(1 x 1) shCoefficients   Swenson-Wahr decorrelated")],
"shCoefficients.gaussian": [
    ("out", "(1 x 1) shCoefficients   Gaussian-smoothed; sigmas scaled "
            "by Wn")],
"shCoefficients.degreeRMS": [
    ("spec", "struct: n, amp/rms/var, err   degree spectrum incl. "
             "formal-error curve")],
"shCoefficients.crossover": [
    ("nc", "(1 x 1) double   signal/error crossover degree (NaN if "
           "none within nmax)")],
"shCoefficients.spectrum": [
    ("h", "(1 x 1) graphics handle   spectrum plot")],
"shCoefficients.triangle": [
    ("h", "(1 x 1) graphics handle   coefficient triangle plot")],
"shCoefficients.synthesis": [
    ("grid", "(nlat x nlon) double   synthesized field"),
    ("lat", "(1 x nlat) double   geocentric latitudes"),
    ("lon", "(1 x nlon) double   longitudes")],
"shCoefficients.toReference": [
    ("out", "(1 x 1) shCoefficients   rescaled to the given GM/R")],
"shCoefficients.subtractNormalField": [
    ("out", "(1 x 1) shCoefficients   disturbing field (even zonals of "
            "the normal field removed)")],
"shCoefficients.map": [
    ("h", "(1 x 1) graphics handle   synthesized map plot")],
"shCoefficients.fan": [
    ("out", "(1 x 1) shCoefficients   fan-filtered")],
"shCoefficients.deformation": [
    ("up", "(nlat x nlon | npts) double   vertical deformation [m]"),
    ("north", "same size   north component [m]"),
    ("east", "same size   east component [m]")],
"shCoefficients.applyDDK": [
    ("out", "(1 x 1) shCoefficients   DDK-decorrelated")],
"shCoefficients.evalAt": [
    ("out", "(1 x 1) shCoefficients   gfct variable terms evaluated at "
            "the epoch; static models pass through unchanged")],
"shCoefficients.vec": [
    ("x", "(P x 1) double   coefficients in idx ordering")],
"shCoefficients.read": [
    ("obj", "(1 x 1) shCoefficients   parsed gfc/gfct(.gz) with GM, R, "
            "sigmas, variable terms, Epoch= as given")],
"shCoefficients.fromVec": [
    ("obj", "(1 x 1) shCoefficients   rebuilt from x with sizes/GM/R/"
            "epoch of the template")],
"shCoefficients.analysis": [
    ("obj", "(1 x 1) shCoefficients   Stokes coefficients estimated "
            "from the grid (exact on ring grids; Kaula for scattered "
            "points)")],
# --------------------------------------------------------- shSeries
"shSeries.shSeries": [
    ("obj", "(1 x 1) shSeries   epoch-sorted monthly stack with "
            "sigmas and history")],
"shSeries.at": [
    ("g", "(1 x 1) shCoefficients   month k with epoch and sigmas")],
"shSeries.mean": [
    ("g", "(1 x 1) shCoefficients   temporal mean field (omits NaN "
          "months)")],
"shSeries.climatology": [
    ("clim", "(1 x 1) shClimatology   fitted bias/trend/annual/"
             "semiannual (+Periods=) with coefficient sigmas"),
    ("resid", "(1 x 1) shSeries   residual series about the fit")],
"shSeries.trendBreaks": [
    ("out", "struct: trends per segment (shCoefficients), F/p per "
            "break, segment epochs")],
"shSeries.fan": [
    ("out", "(1 x 1) shSeries   fan-filtered per month")],
"shSeries.applyTN14": [
    ("out", "(1 x 1) shSeries   C20/C30 replaced epoch-matched")],
"shSeries.addDegree1": [
    ("out", "(1 x 1) shSeries   degree 1 completed epoch-matched")],
"shSeries.removeGIA": [
    ("out", "(1 x 1) shSeries   GIA trend removed about the series "
            "mean epoch")],
"shSeries.applyDDK": [
    ("out", "(1 x 1) shSeries   DDK-decorrelated per month")],
"shSeries.restore": [
    ("out", "(1 x 1) shSeries   background model added back "
            "epoch-matched (e.g. GSM + GAD)")],
"shSeries.minus": [
    ("out", "(1 x 1) shSeries   per-month difference (series or "
            "single field)")],
"shSeries.destripe": [
    ("out", "(1 x 1) shSeries   destriped per month")],
"shSeries.gaussian": [
    ("out", "(1 x 1) shSeries   Gaussian-smoothed per month")],
"shSeries.truncate": [
    ("out", "(1 x 1) shSeries   truncated per month")],
"shSeries.dropNaN": [
    ("out", "(1 x 1) shSeries   gap months removed")],
"shSeries.select": [
    ("out", "(1 x 1) shSeries   months selected by logical/index "
            "mask")],
"shSeries.filter": [
    ("out", "(1 x 1) shSeries   tvANS-filtered series with sigmaCs/Ss "
            "posterior stacks (exact incl. constraints, v2.5)"),
    ("op", "struct   stored linear operator (V/Ut/lam/s or blocks, "
           "model, detLeverage/detResVar) for deconvolution and "
           "resolution maps"),
    ("info", "struct: sigmaXfres (P x T), sigmaNote, VCE diagnostics")],
"shSeries.basinAverage": [
    ("avg", "(K x T) double   basin averages (deconvolved when "
            "Deconvolve=true)"),
    ("out", "struct: sigma (K x T), c, attn, condA (deconvolution "
            "path)")],
"shSeries.read": [
    ("obj", "(1 x 1) shSeries   wildcard-read, epoch-sorted stack")],
"shSeries.fromFolder": [
    ("obj", "(1 x 1) shSeries   all gfc(.gz) files of a folder, "
            "epoch-sorted")],
# ----------------------------------------------------- shClimatology
"shClimatology.withNote": [
    ("obj", "(1 x 1) shClimatology   copy with the note appended to "
            "its history")],
"shClimatology.removeGIA": [
    ("obj", "(1 x 1) shClimatology   GIA rate removed from the trend "
            "component")],
"shClimatology.eval": [
    ("g", "(1 x 1) shCoefficients   model field evaluated at the "
          "epoch")],
"shClimatology.bias": [
    ("g", "(1 x 1) shCoefficients   bias component with sigmas")],
"shClimatology.trend": [
    ("g", "(1 x 1) shCoefficients   trend component [units/yr] with "
          "sigmas")],
"shClimatology.cosAnnual": [
    ("g", "(1 x 1) shCoefficients   cosine annual component")],
"shClimatology.sinAnnual": [
    ("g", "(1 x 1) shCoefficients   sine annual component")],
"shClimatology.cosSemiannual": [
    ("g", "(1 x 1) shCoefficients   cosine semiannual component")],
"shClimatology.sinSemiannual": [
    ("g", "(1 x 1) shCoefficients   sine semiannual component")],
"shClimatology.periodic": [
    ("gc", "(1 x 1) shCoefficients   cosine component of the k-th "
           "extra period"),
    ("gs", "(1 x 1) shCoefficients   sine component")],
"shClimatology.amplitudeMap": [
    ("A", "(nlat x nlon) double   amplitude of the harmonic in the "
          "requested quantity"),
    ("lat", "(1 x nlat) double   latitudes"),
    ("lon", "(1 x nlon) double   longitudes"),
    ("phase", "(nlat x nlon) double   phase [rad] of the harmonic")],
"shClimatology.fromCoef": [
    ("obj", "(1 x 1) shClimatology   rebuilt from a coefficient table "
            "with the series' layout")],
}


def format_block(specs, indent):
    L = [indent + "%   Outputs"]
    for name, desc in specs:
        L.append(indent + "%%     %-10s %s" % (name, desc))
    return [x.replace("%%", "%") for x in L]


def inject_function_file(path, key):
    if key not in OUTSPEC:
        return "no-spec"
    src = open(path).read().split("\n")
    # first function line
    fi = next(i for i, L in enumerate(src) if re.match(r"\s*function\s", L))
    j = fi + 1
    while j < len(src) and src[j].lstrip().startswith("%"):
        j += 1
    help_lines = src[fi + 1:j]
    if any(re.match(r"\s*%\s*Outputs?\s*$", L) for L in help_lines):
        return "present"
    # insert before the attribution line (always present), else at end
    ins = j
    for k in range(fi + 1, j):
        if "Developed by Matthias Weigelt" in src[k]:
            ins = k
            break
    block = format_block(OUTSPEC[key], "") + ["%"]
    src[ins:ins] = block
    open(path, "w").write("\n".join(src))
    return "injected"


def inject_classdef(path, cname):
    src = open(path).read().split("\n")
    changed = 0
    i = 0
    while i < len(src):
        fm = re.match(r"^\s{4}function\s+(?:(?:\[[^\]]*\]|\w+)\s*=\s*)?"
                      r"(\w+)\s*\(", src[i])
        if fm:
            key = cname + "." + fm.group(1)
            j = i + 1
            while j < len(src) and src[j].lstrip().startswith("%"):
                j += 1
            help_lines = src[i + 1:j]
            has = any(re.match(r"\s*%\s*Outputs?\s*$", L)
                      for L in help_lines)
            if key in OUTSPEC and not has and help_lines:
                block = format_block(OUTSPEC[key], "        ")
                src[j:j] = block
                changed += 1
                i = j + len(block)
                continue
        i += 1
    if changed:
        open(path, "w").write("\n".join(src))
    return changed


def main():
    import glob
    import os
    stats = {"injected": 0, "present": 0, "no-spec": 0}
    for f in glob.glob(os.path.join(ROOT, "+shx", "*.m")):
        key = "shx." + os.path.basename(f)[:-2]
        stats[inject_function_file(f, key)] += 1
    for f in ["demo_shAnalysis.m"]:
        stats[inject_function_file(os.path.join(ROOT, f), f[:-2])] += 1
    for f, cn in [("shCoefficients.m", "shCoefficients"),
                  ("shSeries.m", "shSeries"),
                  ("shClimatology.m", "shClimatology")]:
        n = inject_classdef(os.path.join(ROOT, f), cn)
        print(f, "methods injected:", n)
    print(stats)


if __name__ == "__main__":
    main()
