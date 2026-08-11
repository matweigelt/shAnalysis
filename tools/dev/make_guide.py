"""make_guide.py - build shAnalysis_workflow_guide.pdf (Edition 5, v3.1.1).

Rebuilt after the original builder was lost; content reconstructed from the
shipped Edition-2 PDF, updated for v2.4.2, with all display equations
rendered as matplotlib mathtext (Computer Modern) images instead of the
former flattened inline text. Figures live in tools/dev/guide_assets/ IN
THE REPO - they were lost once when the container path they used to be
read from vanished, and were recovered from the shipped PDF.

Snippets are declared via code(\"\"\"...\"\"\") so check_api.py / mlint_lite can
lint every guide code block (gate 2).

Developed by Matthias Weigelt with the help of Claude (Fable 5), 2026-08-07.
"""
import hashlib
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from reportlab.lib import colors  # noqa: E402
from reportlab.lib.pagesizes import A4  # noqa: E402
from reportlab.lib.styles import ParagraphStyle  # noqa: E402
from reportlab.lib.units import cm  # noqa: E402
from reportlab.platypus import (BaseDocTemplate, Frame, Image,  # noqa: E402
                                PageBreak, PageTemplate, Paragraph, Spacer,
                                Table, TableStyle, XPreformatted)

matplotlib.rcParams["mathtext.fontset"] = "cm"

OUT = os.environ.get("SHX_GUIDE_OUT",
                     "/tmp/shx_git/docs/shAnalysis_workflow_guide.pdf")
# figures live IN THE REPO (they were lost once when a container path
# vanished; recovered from the shipped PDF and committed since v3.1.2)
ASSETS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "guide_assets")
EQDIR = "/tmp/guide_eq"
os.makedirs(EQDIR, exist_ok=True)

import re as _re
_cm = open("/tmp/shx_git/Contents.m").read()
VERSION = "v" + _re.search(r"% Version (\S+)", _cm).group(1)
STAMP = ("shAnalysis %s - Workflow & Theory Guide - Developed by "
         "Matthias Weigelt with the help of Claude, 2026-08-11"
         % VERSION)

PAGE_W, PAGE_H = A4
MARGIN = 1.9 * cm
TEXT_W = PAGE_W - 2 * MARGIN

# ----------------------------------------------------------------- styles
S = {}
S["body"] = ParagraphStyle("body", fontName="Helvetica", fontSize=9.3,
                           leading=12.3, spaceAfter=4)
S["bullet"] = ParagraphStyle("bullet", parent=S["body"], leftIndent=14,
                             bulletIndent=4, spaceAfter=2)
S["h1"] = ParagraphStyle("h1", fontName="Helvetica-Bold", fontSize=15,
                         leading=18, spaceBefore=10, spaceAfter=6,
                         textColor=colors.HexColor("#1a3a5c"))
S["h2"] = ParagraphStyle("h2", fontName="Helvetica-Bold", fontSize=11.5,
                         leading=14, spaceBefore=9, spaceAfter=3,
                         textColor=colors.HexColor("#1a3a5c"))
S["h3"] = ParagraphStyle("h3", fontName="Helvetica-Bold", fontSize=9.8,
                         leading=12, spaceBefore=7, spaceAfter=2)
S["title"] = ParagraphStyle("title", fontName="Helvetica-Bold", fontSize=21,
                            leading=26, spaceAfter=4,
                            textColor=colors.HexColor("#1a3a5c"))
S["subtitle"] = ParagraphStyle("subtitle", fontName="Helvetica", fontSize=11,
                               leading=14, spaceAfter=14,
                               textColor=colors.HexColor("#444444"))
S["code"] = ParagraphStyle("code", fontName="Courier", fontSize=7.8,
                           leading=9.6)
S["caption"] = ParagraphStyle("caption", parent=S["body"], fontSize=8.6,
                              leading=11, textColor=colors.HexColor("#333333"),
                              spaceBefore=2, spaceAfter=10)
S["cell"] = ParagraphStyle("cell", parent=S["body"], fontSize=8.6,
                           leading=10.8, spaceAfter=0)
S["cellb"] = ParagraphStyle("cellb", parent=S["cell"],
                            fontName="Helvetica-Bold")

TBL = TableStyle([
    ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#b7c4d1")),
    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e8eef4")),
    ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
    ("TOPPADDING", (0, 0), (-1, -1), 3),
    ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
])
TBL_PLAIN = TableStyle([
    ("VALIGN", (0, 0), (-1, -1), "TOP"),
    ("TOPPADDING", (0, 0), (-1, -1), 2),
    ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
])
CODEBOX = TableStyle([
    ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#f4f6f8")),
    ("BOX", (0, 0), (-1, -1), 0.4, colors.HexColor("#c9d3dc")),
    ("LEFTPADDING", (0, 0), (-1, -1), 7),
    ("RIGHTPADDING", (0, 0), (-1, -1), 7),
    ("TOPPADDING", (0, 0), (-1, -1), 5),
    ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
])


# ------------------------------------------------------------- equations
def eq(tex, fontsize=11.5, pad=0.35):
    """Render display math via matplotlib mathtext -> centered Image."""
    key = hashlib.md5((tex + str(fontsize)).encode()).hexdigest()[:12]
    png = os.path.join(EQDIR, key + ".png")
    if not os.path.exists(png):
        fig = plt.figure(figsize=(0.01, 0.01))
        fig.text(0, 0, "$%s$" % tex, fontsize=fontsize)
        fig.savefig(png, dpi=300, bbox_inches="tight",
                    pad_inches=0.02, transparent=True)
        plt.close(fig)
    from PIL import Image as PILImage
    w, h = PILImage.open(png).size
    scale = 72.0 / 300.0
    wpt, hpt = w * scale, h * scale
    if wpt > TEXT_W - 20:
        f = (TEXT_W - 20) / wpt
        wpt, hpt = wpt * f, hpt * f
    img = Image(png, width=wpt, height=hpt)
    img.hAlign = "CENTER"
    return [Spacer(1, pad * 10), img, Spacer(1, pad * 10)]


def ieq(tex, fontsize=8.8):
    """Small equation image for table cells (left-aligned)."""
    key = hashlib.md5(("i" + tex + str(fontsize)).encode()).hexdigest()[:12]
    png = os.path.join(EQDIR, key + ".png")
    if not os.path.exists(png):
        fig = plt.figure(figsize=(0.01, 0.01))
        fig.text(0, 0, "$%s$" % tex, fontsize=fontsize)
        fig.savefig(png, dpi=300, bbox_inches="tight",
                    pad_inches=0.02, transparent=True)
        plt.close(fig)
    from PIL import Image as PILImage
    w, h = PILImage.open(png).size
    scale = 72.0 / 300.0
    img = Image(png, width=w * scale, height=h * scale)
    img.hAlign = "LEFT"
    return img


SNIPPETS = []


def code(src):
    """Register + typeset a MATLAB snippet (linted by check_api/mlint)."""
    SNIPPETS.append(src)
    pre = XPreformatted(src.strip("\n"), S["code"])
    t = Table([[pre]], colWidths=[TEXT_W])
    t.setStyle(CODEBOX)
    return [Spacer(1, 3), t, Spacer(1, 5)]


def para(txt, style="body"):
    return Paragraph(txt, S[style])


def bullets(items):
    return [Paragraph("&bull; " + it, S["bullet"]) for it in items]


def fig(png, caption, width=None):
    from PIL import Image as PILImage
    p = os.path.join(ASSETS, png)
    w, h = PILImage.open(p).size
    wpt = width or min(TEXT_W, w * 72.0 / 271.0)
    img = Image(p, width=wpt, height=wpt * h / w)
    img.hAlign = "CENTER"
    return [Spacer(1, 6), img, Paragraph(caption, S["caption"])]


def tbl(header, rows, widths, style=None):
    data = [[Paragraph(c, S["cellb"]) for c in header]]
    for r in rows:
        data.append([c if not isinstance(c, str)
                     else Paragraph(c, S["cell"]) for c in r])
    t = Table(data, colWidths=widths, repeatRows=1)
    t.setStyle(style or TBL)
    return [Spacer(1, 3), t, Spacer(1, 6)]


# ------------------------------------------------------------ page frame
def on_page(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(colors.HexColor("#666666"))
    canvas.drawString(MARGIN, 1.05 * cm, STAMP)
    canvas.drawRightString(PAGE_W - MARGIN, 1.05 * cm, "p. %d" % doc.page)
    canvas.restoreState()


doc = BaseDocTemplate(OUT, pagesize=A4,
                      leftMargin=MARGIN, rightMargin=MARGIN,
                      topMargin=1.6 * cm, bottomMargin=1.7 * cm,
                      title="shAnalysis Workflow & Theory Guide",
                      author="Matthias Weigelt with Claude (Fable 5)")
frame = Frame(MARGIN, 1.7 * cm, TEXT_W, PAGE_H - 3.3 * cm, id="main")
doc.addPageTemplates([PageTemplate(id="page", frames=[frame],
                                   onPage=on_page)])

story = []

# ================================================================= title
story += [Spacer(1, 30),
          para("shAnalysis Workflow &amp; Theory Guide", "title"),
          para("GRACE / GRACE-FO spherical-harmonic processing in base "
               "MATLAB &mdash; theory, best practices, and step-by-step "
               "recipes", "subtitle")]
story += tbl(["", ""], [
    ["Toolbox version", "shAnalysis " + VERSION],
    ["Edition", "5 (adds: ICGEM time-series download, the standard chain "
                "and custom filter design as single entry points, and what "
                "real provider files actually look like &mdash; FORTRAN "
                "D-exponents, the ICGEM 2.0 column order, ragged record "
                "groups)"],
    ["Scope", "ICGEM I/O, spectral diagnostics, synthesis &amp; analysis, "
              "destriping, Gaussian / DDK / tvANS filtering, VCE, basin "
              "averages &amp; deconvolution, climatology &amp; AR(1) "
              "significance, GIA, degree-1/C20 chain, Slepian localization, "
              "SINEX covariances, mascon comparison, load deformation, "
              "gradient tensor, multi-center combination, sea-level "
              "fingerprints, Love-number I/O incl. frame conversion, "
              "ICGEM series download, pole-tide conventions, pointwise "
              "synthesis at satellite altitude"],
    ["Dependencies", "base MATLAB only (no toolboxes); Java optionally for "
                     "gz streaming"],
    ["Validation", "every numerical method independently validated in "
                   "Python (numpy/scipy) before MATLAB implementation; "
                   "169 unit tests"],
    ["License", "MIT (see LICENSE; bundled DDK3 Wbd file MIT from its "
                "upstream repository; ITSG/GFZ fixtures remain their "
                "providers' property)"],
    ["Provenance", "Developed by Matthias Weigelt with the help of "
                   "Claude (Fable 5), 2026-08-07"],
], [3.4 * cm, TEXT_W - 3.4 * cm],
    TableStyle([("GRID", (0, 0), (-1, -1), 0.4,
                 colors.HexColor("#b7c4d1")),
                ("BACKGROUND", (0, 0), (0, -1),
                 colors.HexColor("#e8eef4")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("TOPPADDING", (0, 0), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3)]))
story += [Spacer(1, 8),
          para("<b>How to read this guide.</b> Part I develops the theory "
               "exactly as implemented (equation forms match the code, "
               "including storage conventions and numerically validated "
               "tolerances). Part II condenses the best practices the "
               "implementation enforces or enables. Part IV is the "
               "complete API reference &mdash; every public function and "
               "method with typed, dimensioned inputs/outputs, defaults "
               "and a real-data example, generated from the source "
               "arguments blocks at build time. Part III gives "
               "copy-paste step-by-step recipes for the standard products. "
               "The HTML reference (<font face='Courier'>doc shAnalysis"
               "</font>) documents every function signature; this guide "
               "explains why and in which order."),
          PageBreak()]

# ========================================================= Part I theory
story += [para("Part I &mdash; Theory as implemented", "h1")]

story += [para("1. Spherical harmonics, normalization, storage", "h2"),
          para("The disturbing potential is expanded in fully "
               "(4&pi;-)normalized real spherical harmonics:")]
story += eq(r"V(r,\vartheta,\lambda)=\frac{GM}{R}\sum_{n=0}^{N}"
            r"\left(\frac{R}{r}\right)^{n+1}\sum_{m=0}^{n}"
            r"\left[\bar{C}_{nm}\cos m\lambda+\bar{S}_{nm}\sin m\lambda"
            r"\right]\bar{P}_{nm}(\cos\vartheta)")
story += [para("with the normalization &int;|Y<sub>nm</sub>|&sup2; "
               "d&Omega; = 4&pi;, the ICGEM convention. Throughout the "
               "toolbox coefficients live in lower-triangular matrices "
               "with 1-based indexing <font face='Courier'>C(n+1, m+1)"
               "</font>; vectorized orderings are always mediated by "
               "<font face='Courier'>shLowLevel.shIndex</font> (fields n, m, cs, "
               "P; MinDegree defaults to 2 for GRACE work, 0 for "
               "analysis). GM and R ride with the data (read from the gfc "
               "header); the defaults (3.986004415e14, 6378136.3) are "
               "overridable, never baked-in physics. Load Love numbers "
               "k<sub>n</sub> for water-storage kernels must be supplied "
               "explicitly &mdash; the toolbox refuses to guess them.")]

story += [para("2. Legendre functions: scaled Holmes&ndash;Featherstone "
               "recursion", "h2"),
          para("Fully normalized associated Legendre functions are "
               "computed by the standard order-diagonal seed plus degree "
               "recursion. The naive diagonal seed underflows near the "
               "poles at high degree; <font face='Courier'>"
               "shLowLevel.legendreALF</font> therefore implements the globally "
               "scaled Holmes&ndash;Featherstone variant (scale "
               "10<super>-280</super>, rescaling on the fly), verified "
               "stable to degree 2190: the sum rule and parity")]
story += eq(r"\sum_{m=0}^{n}\bar{P}_{nm}^{\,2}(\cos\vartheta)=2n+1"
            r"\ \ (\sim 10^{-11}\ \mathrm{at}\ n=2190),\qquad"
            r"\bar{P}_{nm}(-\varphi)=(-1)^{n+m}\,\bar{P}_{nm}(\varphi)")
story += [para("hold across all latitudes; the parity relation halves the "
               "recursion cost on hemispherically symmetric grids (the "
               "parity trick).")]

story += [para("3. Synthesis", "h2"),
          para("Spatial synthesis evaluates, per latitude, the order sums")]
story += eq(r"A_m(\varphi)=\sum_{n=m}^{N}f_n\,\bar{C}_{nm}\,"
            r"\bar{P}_{nm}(\varphi),\qquad B_m(\varphi)\ \mathrm{"
            r"analogously\ with}\ \bar{S}_{nm},")
story += [para("then the longitude sum")]
story += eq(r"g(\varphi,\lambda)=\sum_{m=0}^{N}\left[A_m(\varphi)"
            r"\cos m\lambda+B_m(\varphi)\sin m\lambda\right].")
story += [para("For uniform full-circle longitude grids this last sum is "
               "an inverse FFT after folding aliased orders into "
               "mod(m, n<sub>lon</sub>) bins &mdash; exact, validated "
               "against the direct product to 3e-13 (<font face='Courier'>"
               "Method=\"auto\"|\"fft\"|\"direct\"</font>)."),
          para("<b>Memory model (v2.2).</b> The Legendre stack "
               "(n<sub>max</sub>+1)&sup2; &times; n<sub>lat</sub> &times; "
               "8 bytes is the memory driver: 7 GB at n<sub>max</sub> 2190 "
               "over 181 latitudes. Above the 'MaxMemGB' budget (default "
               "4), synthesis streams latitude bands: compute the Legendre "
               "block, fold into A<sub>m</sub>/B<sub>m</sub>, discard. "
               "Chunked equals monolithic exactly (Python-validated, max "
               "deviation 0.0); the class-level Legendre cache is bypassed "
               "automatically when the stack would not fit.")]

story += [para("4. Quantity kernels", "h2"),
          para("A single degree factor f<sub>n</sub> maps Stokes "
               "coefficients to the target quantity (<font face='Courier'>"
               "shLowLevel.kernelFactors</font>):")]
story += tbl(["quantity", "f<sub>n</sub>", "unit"], [
    ["geoid", ieq(r"R"), "m"],
    ["potential", ieq(r"GM/R"), "m&sup2;/s&sup2;"],
    ["gravity disturbance", ieq(r"\frac{GM}{R^2}\,(n+1)"), "m/s&sup2;"],
    ["gravity anomaly", ieq(r"\frac{GM}{R^2}\,(n-1)"), "m/s&sup2;"],
    ["EWH (water)",
     ieq(r"R\;\frac{\rho_{\mathrm{ave}}}{3\rho_w}\;\frac{2n+1}{1+k_n}"),
     "m"],
    ["surface density *",
     ieq(r"R\;\frac{\rho_{\mathrm{ave}}}{3}\;\frac{2n+1}{1+k'_n}"),
     "kg/m&sup2;"],
    ["bottom pressure *",
     ieq(r"g_0R\,\frac{\rho_{\mathrm{ave}}}{3}\,\frac{2n+1}{1+k'_n},"
         r"\quad g_0=\frac{GM}{R^2}"), "Pa"],
    ["vertical deformation *", ieq(r"R\;\frac{h'_n}{1+k'_n}"), "m"],
    ["gravity gradient T<sub>rr</sub> *",
     ieq(r"\frac{GM}{R^3}\,(n+1)(n+2)"), "1/s&sup2;"],
], [4.6 * cm, TEXT_W - 4.6 * cm - 2.3 * cm, 2.3 * cm])
story += [para("(* new in v2.3.) The 'Height' option continues "
               "potential-type kernels upward by (R/r)<super>p</super> "
               "(p = n+1 potential, n+2 anomaly/disturbance, n+3 "
               "T<sub>rr</sub>; validated against -d/dr and "
               "d&sup2;/dr&sup2; of the continued potential) &mdash; "
               "evaluate directly at satellite altitude. Surface "
               "quantities reject Height by contract."),
          para("The EWH kernel divides by (1+k<sub>n</sub>): the observed "
               "potential contains the solid-Earth loading response, which "
               "must be removed to recover surface mass. k<sub>n</sub> is "
               "always user-supplied (e.g. Wahr et al. 1998 or your "
               "loading-theory set of choice).")]

story += [para("5. Elastic load deformation (v2.3)", "h2"),
          para("A surface load deforms the elastic Earth; the displacement "
               "at the surface, expressed through the load-induced "
               "(residual) Stokes coefficients, is (Wahr et al. 1998):")]
story += eq(r"u_{\mathrm{up}}=R\sum_{n,m}\frac{h'_n}{1+k'_n}\,"
            r"\Delta\bar{C}_{nm}\,\bar{Y}_{nm},\qquad"
            r"u_{\mathrm{north}}=R\sum_{n,m}\frac{l'_n}{1+k'_n}\,"
            r"\Delta\bar{C}_{nm}\,\frac{\partial\bar{Y}_{nm}}"
            r"{\partial\varphi},\qquad"
            r"u_{\mathrm{east}}=R\sum_{n,m}\frac{l'_n}{1+k'_n}\,"
            r"\Delta\bar{C}_{nm}\,\frac{1}{\cos\varphi}\,"
            r"\frac{\partial\bar{Y}_{nm}}{\partial\lambda}", 10.5)
story += [para("The vertical is a pure degree factor ('deformation_up'); "
               "the horizontals need the exact spherical-harmonic "
               "gradient. <font face='Courier'>shLowLevel.legendreALFDeriv</font> "
               "provides d<font face='Courier'>Pbar</font>/d&phi; through "
               "a frozen identity &mdash; the (n,m) coefficient pattern "
               "(with &radic;2 corrections at m=0 and m=1 from the "
               "normalization of Pbar<sub>n0</sub>/Pbar<sub>n1</sub>) was "
               "calibrated by least squares against numerical "
               "differentiation and then frozen; residuals sit at the "
               "finite-difference floor (~1e-9). <font face='Courier'>"
               "shLowLevel.shSynthesisDeformation</font> / <font face='Courier'>"
               "g.deformation</font> evaluate all three components on "
               "grids or GNSS station lists (Mode=\"points\", "
               "LatType=\"geodetic\"); north/east validated point-wise "
               "against brute-force gradients to ~5e-10. Degree 0 is "
               "excluded; degree 1 is opt-in and only meaningful with "
               "consistent geocenter (CF/CE frame) handling; all three "
               "Love-number sets must come from ONE loading model. East "
               "is NaN at the poles (undefined direction).")]

story += [para("6. Analysis (the inverse problem)", "h2"),
          para("<font face='Courier'>shLowLevel.shAnalysisGrid</font> / "
               "<font face='Courier'>shCoefficients.analysis</font> "
               "estimate coefficients from data. Two regimes:")]
story += bullets([
    "<b>Ring grids</b> (uniform full-circle longitudes, arbitrary "
    "latitude rings): an FFT along longitude decouples the normal "
    "equations per order; each order solves a small least-squares over "
    "latitude rings. Exact recovery of band-limited input to ~8.7e-15 "
    "when the sampling theorem is met (n<sub>lat</sub> &ge; "
    "n<sub>max</sub>+1, n<sub>lon</sub> &ge; 2n<sub>max</sub>+1).",
    "<b>Scattered points</b>: chunked accumulation of full normal "
    "equations. Under-determined or ill-conditioned samplings (polar "
    "gaps, regional patches) need Kaula regularization: a diagonal prior"])
story += eq(r"\sigma^2_{\mathrm{prior}}(n)=\left(\frac{K}{n^{2}}"
            r"\right)^{2}")
story += [para("stabilizes the inversion at the cost of damping "
               "unconstrained coefficients toward zero. For regional data "
               "the mathematically cleaner route is Slepian localization "
               "(Section 12) &mdash; estimate only the ~N concentrated "
               "linear combinations the data can actually see.")]

story += [para("7. Destriping and Gaussian smoothing", "h2"),
          para("GRACE monthly solutions carry correlated north&ndash;south "
               "stripes: within one order and parity, coefficients of "
               "successive degrees correlate. The classic "
               "Swenson&ndash;Wahr destriper (<font face='Courier'>"
               "shLowLevel.shDestripe</font>) fits and removes a low-order "
               "polynomial across same-parity degree series per order "
               "(defaults: orders &ge; 8, polynomial order 2, window 5, "
               "all configurable). Gaussian smoothing "
               "(<font face='Courier'>shLowLevel.shGaussianFilter</font>) "
               "multiplies degree-wise by Jekeli weights from the "
               "recursion")]
story += eq(r"b=\frac{\ln 2}{1-\cos(r/R)},\qquad W_0=1,\qquad"
            r"W_1=\frac{1+e^{-2b}}{1-e^{-2b}}-\frac{1}{b},\qquad"
            r"W_{n+1}=-\frac{2n+1}{b}\,W_n+W_{n-1}")
story += [para("with the half-response radius r as the single tuning "
               "knob. Both are provided primarily as the baseline against "
               "which tvANS and DDK are judged.")]

story += [para("8. DDK anisotropic decorrelation (v2.2)", "h2"),
          para("The DDK filters (Kusche 2007; Kusche et al. 2009: "
               "DDK1&ndash;DDK8, ordered strong&rarr;weak) replace the "
               "scalar Gaussian weight by a full matrix per (order, C/S) "
               "block, derived from a GRACE error covariance and a signal "
               "prior &mdash; anisotropic, order-dependent smoothing that "
               "suppresses stripes with less signal attenuation than an "
               "equivalent Gaussian. <font face='Courier'>shLowLevel.readDDK"
               "</font> parses the released binary Wbd files natively "
               "(Rietbroek/Kusche BIN format; endian autodetect; verified "
               "to 4.4e-16 against the repository-documented DDK3 "
               "reference values; DDK1..8 = a_1d14p_4, 1d13, 1d12, 5d11, "
               "1d11, 5d10, 1d10, 5d9; Nmax= truncates an Lmax-120 filter "
               "to n96 series), plus .mat containers and an ASCII "
               "exchange layout; <font face='Courier'>applyDDK</font> "
               "applies block-wise c &rarr; M&middot;c, passing uncovered "
               "degrees through unchanged. Formal errors are invalidated "
               "by the filter (coefficients become correlated); propagate "
               "uncertainties with <font face='Courier'>shLowLevel.mcPropagate"
               "</font>."),
          para("<b>Note.</b> The real DDK3 file (Wbd_2-120.a_1d12p_4, "
               "MIT-licensed, from github.com/strawpants/GRACE-filter) "
               "ships in tests/test_data; the unit tests pin the parser "
               "to the repository-documented pack values. Degrees below "
               "the filter's Lmin (2) pass through unfiltered, as in the "
               "reference implementation.")]

story += [para("9. The tvANS filter (time-variable anisotropic "
               "noise/signal Wiener filter)", "h2"),
          para("The toolbox's flagship decorrelation filter. Model per "
               "month t: x = signal + noise with signal covariance S "
               "(estimated from the data themselves) and noise covariance "
               "s<sub>t</sub>N (a shape N, e.g. from a SINEX covariance "
               "or the built-in order/parity model, scaled monthly by VCE "
               "factors s<sub>t</sub>). The Wiener estimate of the "
               "residual field is")]
story += eq(r"\hat{x}_t=W_t\,x_t,\qquad W_t=S\,(S+s_tN)^{-1}")
story += [para("Implementation: one generalized symmetric "
               "eigendecomposition gives, for every month simultaneously,")]
story += eq(r"SU=NU\Lambda\ \ \Rightarrow\ \ W_t=V\,\mathrm{diag}\!"
            r"\left(\frac{\lambda_i}{\lambda_i+s_t}\right)U^{\mathsf{T}},"
            r"\qquad V=U^{-\mathsf{T}}")
story += [para("&mdash; O(P&sup2;) per month instead of a fresh P&sup3; "
               "solve. The deterministic model (bias/trend/annual/"
               "semi-annual) is removed first and restored after "
               "filtering; VCE factors are estimated from high-degree "
               "residuals via a robust median &chi;&sup2; estimator "
               "(<font face='Courier'>shLowLevel.vceRescale</font>). With a "
               "block-diagonal N (order/parity blocks), the whole "
               "machinery runs per block (Blocks=\"auto\"), tractable to "
               "Lmax &asymp; 120. External noise models join the block "
               "path when they are block-diagonal in the (order, C/S, "
               "parity) partition (verified, not assumed &mdash; v2.2.2); "
               "dense covariances fall back to the full path (auto) or "
               "error loudly (on)."),
          para("<b>Per-order-band VCE (v2.2).</b> One s<sub>t</sub> per "
               "month is coarse: striping strength varies with order. "
               "VCEBands=[0 16 33 61] estimates independent monthly "
               "factors per order band (block path only &mdash; band-wise "
               "scaling breaks the single global eigendecomposition, "
               "which the code asserts). The per-block factors are "
               "honored consistently in the filter application, the "
               "posterior sigmas, and the basin deconvolution covariance "
               "(identity of the banded block filter with the directly "
               "constructed W validated to 5.6e-16)."),
          para("<b>Posterior uncertainty.</b> The filtered-residual "
               "variance")]
story += eq(r"\mathrm{diag}\left((I-W_t)\,S\right)=(V\circ V)\;"
            r"\frac{s_t\lambda}{\lambda+s_t}")
story += [para("(validated 3.5e-16) is exposed as per-coefficient/month "
               "sigma stacks (<font face='Courier'>tsF.sigmaCs/Ss"
               "</font>). With hard constraints (v2.5) the EXACT diagonal "
               "of (W&minus;I)S(W&minus;I)' + s&middot;WNW' is used "
               "(O(P&sup2;q) in the eig basis, validated 2e-15 against "
               "the brute-force covariance; the constrained-direction "
               "identity Ac'&middot;Cov&middot;Ac = s&middot;Ac'NAc holds "
               "exactly). Honest correction: the earlier formula "
               "UNDERESTIMATED along constrained directions (ratios down "
               "to ~0.72) &mdash; the pre-v2.5 'upper bound' note had the "
               "direction wrong.")]

story += [para("10. Basin averages and deconvolution", "h2"),
          para("A basin kernel b (Section 11) gives the filtered basin "
               "series a<sub>t</sub> = b<super>T</super>x&#770;"
               "<sub>t</sub> / (b<super>T</super>b). Because filtering "
               "attenuates and leaks, <font face='Courier'>"
               "shLowLevel.basinDeconvolve</font> solves the small linear system "
               "that undoes the filter's action on a set of basin kernels "
               "B:")]
story += eq(r"c_t=A^{-1}B^{\mathsf{T}}\hat{x}_t,\quad "
            r"A=B^{\mathsf{T}}W_tB\ \ (K\times K);\qquad"
            r"\mathrm{cov}(c_t)=A^{-1}\left[G_v^{\mathsf{T}}\,"
            r"\mathrm{diag}(s_{\mathrm{vec}}\,g^2)\,G_v\right]"
            r"A^{-\mathsf{T}}")
story += [para("where g are the spectral gains and s<sub>vec</sub> the "
               "(possibly band-wise) monthly noise factors &mdash; the "
               "v2.2 form is exactly the v2.1 formula when bands are off. "
               "Ill-conditioned kernel sets accept a relative ridge "
               "(&lambda; = Ridge&middot;||A||), with the covariance "
               "consistently using the regularized inverse. Since v2.5 "
               "the basin sigma also carries the OLS parameter "
               "uncertainty of the restored deterministic part "
               "(design leverage h<sub>t</sub> &times; per-coefficient "
               "residual variance, carried on the operator; Monte-Carlo: "
               "emp/pred 1.18 without the term, 1.00 with it, 1.26 "
               "&rarr; 1.08 at seasonal leverage peaks), and with "
               "constraints the deconvolution noise covariance uses the "
               "exact constrained operator. The deterministic-part "
               "errors remain correlated across epochs through the "
               "shared fit (pointwise sigma; documented).")]

story += [para("11. Basin kernels: buffering and tapering (v2.2)", "h2"),
          para("<font face='Courier'>shLowLevel.basinKernel</font> builds "
               "b = Y<super>T</super>(w&middot;mask) by exact "
               "Gauss&ndash;Legendre quadrature from a polygon, mask, or "
               "f(lat,lon) handle. Two leakage controls:")]
story += bullets([
    "<b>BufferKm</b>: grow/shrink the region by exact great-circle "
    "distance before quadrature (trade leakage-in against leakage-out; "
    "computed in latitude chunks, practical to Lmax ~ 96).",
    "<b>TaperKm</b>: multiply the kernel spectrally by Jekeli Gaussian "
    "weights &mdash; soft edges suppress the ringing of the truncated "
    "indicator. Rule of thumb: taper radius &asymp; the coefficient "
    "filter's half-response radius."])
story += [para("The quadrature area fraction is returned; b(1) equals it "
               "exactly (Y<sub>00</sub> = 1), which the tests pin down.")]

story += [para("12. Slepian localization (v2.2)", "h2"),
          para("For a region R, the concentration problem maximizes the "
               "energy ratio over band-limited g = Yx:")]
story += eq(r"\lambda=\frac{\int_{R}g^{2}\,d\Omega}"
            r"{\int_{\Omega}g^{2}\,d\Omega}\ \rightarrow\ \max,"
            r"\qquad K=Y^{\mathsf{T}}\mathrm{diag}(w\cdot\mathrm{mask})"
            r"\,Y,\qquad N_{\mathrm{Shannon}}=\sum_i\lambda_i="
            r"\frac{A_R}{4\pi}\,P")
story += [para("Eigenvectors (Slepian tapers) are orthonormal; "
               "eigenvalues are concentrations in [0,1]; the Shannon "
               "number (exact by the addition theorem &mdash; the trace "
               "identity the tests verify) counts the tapers the region "
               "supports. Regional analysis then estimates ~N Slepian "
               "coefficients a = G<super>T</super>x instead of P Stokes "
               "coefficients &mdash; well-posed by construction, no Kaula "
               "prior needed. Validated on the 30&deg; polar cap: "
               "&lambda; bounds to 4e-16, orthonormality 3e-15.")]

story += [para("13. Climatology, alias periods, AR(1) significance "
               "(v2.2)", "h2"),
          para("<font face='Courier'>ts.climatology()</font> fits bias, "
               "trend, annual and semi-annual terms (plus optional extra "
               "periods) per coefficient. GRACE tidal aliasing puts S2 at "
               "161 d, K2 at 3.66 yr, K1 at 7.48 yr: "
               "Periods=[161/365.25, 3.66, 7.48] keeps these out of trend "
               "and semi-annual estimates. Per-coefficient 1-sigma "
               "significance comes from the OLS covariance (weighted and "
               "robust variants available) and rides on every component "
               "accessor, enabling significance-masked maps."),
          para("<b>AR(1) correction.</b> GRACE residuals are temporally "
               "correlated; white-noise sigmas are optimistic. "
               "ARCorrect=true estimates the lag-1 residual "
               "autocorrelation r<sub>1</sub> per coefficient and "
               "inflates sigmas by")]
story += eq(r"\sigma_{\mathrm{AR(1)}}=\sigma\;\sqrt{\frac{1+r_1}{1-r_1}}")
story += [para("&mdash; the standard first-order effective-sample-size "
               "correction, since v2.5 applied to the Kendall bias-"
               "corrected r<sub>1</sub>' = r<sub>1</sub> + "
               "(1+3r<sub>1</sub>)/T. Monte-Carlo validation: uncorrected "
               "sigmas underestimate the empirical trend scatter by "
               "2.00&times; (&phi; = 0.6, T = 120); with the raw-"
               "r<sub>1</sub> inflation the ratio is 1.065, with the "
               "bias-corrected r<sub>1</sub> 1.028 (1.152 &rarr; 1.079 "
               "at T = 60). Still first-order and slightly conservative "
               "&mdash; the residual few percent stems from the "
               "regression absorbing low-frequency noise.")]

story += [para("14. GIA correction (v2.2)", "h2"),
          para("Secular Stokes trends mix present-day mass change with "
               "glacial isostatic adjustment. <font face='Courier'>"
               "ts.removeGIA(gia)</font> subtracts rate&middot;"
               "(t - T<sub>0</sub>) per epoch; <font face='Courier'>"
               "clim.removeGIA(gia)</font> corrects a fitted trend "
               "directly &mdash; both routes agree identically (tested to "
               "1e-13). GIA models (ICE-6G_D, Caron et al., A/W13, ...) "
               "are user-supplied gfc rate files, never embedded. The "
               "model is treated as exact: its uncertainty, often the "
               "dominant regional trend error (Laurentia, Fennoscandia, "
               "Antarctica), is not propagated &mdash; difference several "
               "models as a spread estimate.")]

story += [para("15. Geocentric vs. geodetic latitude (v2.2)", "h2"),
          para("All SH mathematics here runs in geocentric latitude; "
               "maps, GIS products and mascon grids are geodetic:")]
story += eq(r"\tan\varphi'=(1-f)^{2}\,\tan\varphi")
story += [para("with the WGS84 flattening as overridable default. The "
               "difference peaks at 0.192&deg; (~21 km) near 45&deg; "
               "&mdash; large enough to visibly bias basin averages and "
               "point comparisons. Synthesis, analysis and "
               "shAnalysisGrid accept LatType=\"geodetic\" and convert "
               "internally.")]

story += [para("16. The degree-1 / C20 completion chain", "h2"),
          para("GRACE does not observe geocenter motion (degree 1), and "
               "its C20 is weak: the standard product chain is GSM + "
               "TN-14 (C20/C30 from SLR) + TN-13 (degree-1). "
               "<font face='Courier'>shLowLevel.readTN14</font> and "
               "<font face='Courier'>shLowLevel.readTN13</font> parse the "
               "official technical notes (TN-13 verified against the real "
               "GFZ RL06 file: 256 paired months 2002-04 to 2026-04, "
               "including the 2017.4&ndash;2018.4 gap); "
               "<font face='Courier'>ts.addDegree1</font> and the "
               "C20/C30 replacement close the chain with epoch matching "
               "and provenance notes in the history.")]

story += [para("17. Normal field and (GM, R) conventions (v2.4)", "h2"),
          para("Geoid maps need the DISTURBING field: subtract the "
               "ellipsoidal normal potential. <font face='Courier'>"
               "shLowLevel.normalFieldCS</font> computes the even zonals from "
               "the WGS84/GRS80 defining constants (GM, a, f, &omega;) "
               "via the Heiskanen&ndash;Moritz closed form "
               "(q<sub>0</sub> series &rarr; J<sub>2</sub> &rarr; "
               "J<sub>2n</sub>) &mdash; no coefficient tables; validated "
               "against NIMA TR8350.2 to all published digits. Crucially, "
               "WGS84 GM (3.986004418e14) and a (6378137.0) DIFFER from "
               "the ICGEM conventions (3.986004415e14, 6378136.3): "
               "<font face='Courier'>g.subtractNormalField</font> "
               "rescales the ellipsoid field to the model's (GM, R) "
               "first via")]
story += eq(r"\bar{C}'_{nm}=\bar{C}_{nm}\;\frac{GM_1}{GM_2}"
            r"\left(\frac{R_1}{R_2}\right)^{n}")
story += [para("(<font face='Courier'>shLowLevel.rescaleGMR</font> / "
               "<font face='Courier'>g.toReference</font>; the physical "
               "field is invariant &mdash; verified by synthesizing the "
               "potential at a common radius from both representations). "
               "Skipping the rescaling costs millimetres. The arithmetic "
               "operators refuse GM/R-mismatched operands "
               "(shCoefficients:constantsMismatch) instead of silently "
               "mixing conventions. Permanent tide: the ellipsoid is "
               "tide-free by construction; the model's C20 tide system "
               "(ICGEM header) is NOT converted silently &mdash; "
               "zero-tide vs tide-free differ by ~4.2e-9 in "
               "Cbar<sub>20</sub>.")]

story += [para("18. Multi-center combination (v2.4)", "h2"),
          para("Several centers observe the same monthly field with "
               "different noise. One variance factor per (center, month) "
               "&mdash; the COST-G granularity &mdash; estimated by "
               "Foerstner iteration with PARTIAL redundancies:")]
story += eq(r"\hat{x}=\left(\sum_c w_cN_c^{-1}\right)^{-1}"
            r"\sum_c w_cN_c^{-1}x_c,\qquad "
            r"s_c^{2}=\frac{v_c^{\mathsf{T}}N_c^{-1}v_c}{r_c},\qquad "
            r"r_c=P-\mathrm{tr}\!\left(H\,w_cN_c^{-1}\right)")
story += [para("Without r<sub>c</sub> the factors bias low "
               "(v<sub>c</sub> is correlated with x&#770;). Noise shapes "
               "N<sub>c</sub> come from each center's own climatology "
               "residuals as (order, C/S, parity) blocks; the whole "
               "iteration runs block-wise. Robust=true replaces the "
               "global ratio by the median of per-block unbiased "
               "estimates q<sub>b</sub>/r<sub>b</sub>. The empirical "
               "shapes are trace-normalized, so the factors carry each "
               "center's absolute noise power and are directly comparable "
               "(weights invariant; v2.4.1). Two limits are structural: "
               "common-mode errors (shared K-band data, AOD1B, tides) are "
               "invisible to inter-center VCE &mdash; posterior sigmas "
               "are a LOWER bound &mdash; and inter-center correlations "
               "violate independence; info.interCenterCorr quantifies the "
               "optimism. Combine GSM before TN-14/TN-13, then apply the "
               "replacements once.")]

story += [para("19. Elastic sea-level fingerprints (v2.4)", "h2"),
          para("A melting ice mass does not raise the ocean uniformly: "
               "the ocean surface follows the perturbed geoid minus the "
               "deformed crust, under global mass conservation (Farrell "
               "&amp; Clark 1976, elastic limit):")]
story += eq(r"S=\mathcal{O}\,(N-U+\kappa),\qquad "
            r"\int_{\mathrm{ocean}}\rho_w\,S\;dA="
            r"-\int_{\mathrm{land}}\sigma\;dA")
story += [para("<font face='Courier'>shLowLevel.seaLevelFingerprint</font> "
               "iterates this on the toolbox's exact Gauss&ndash;Legendre "
               "quadrature pair. Validation: mass conserved to machine "
               "precision, ~11 iterations, and the classic pattern "
               "&mdash; sea level FALLS near the melting mass (lost "
               "self-attraction plus rebound; validation run: near field "
               "-0.3 mm vs far field +0.8 mm about a +0.6 mm eustatic "
               "mean per unit load). Elastic only, no rotational "
               "feedback, fixed coastlines.")]

story += [para("20. Gravity gradient tensor (v2.4)", "h2"),
          para("<font face='Courier'>shLowLevel.shSynthesisGradientTensor"
               "</font> synthesizes all six NEU components at altitude; "
               "second Legendre derivatives come from the frozen "
               "first-derivative stencil applied twice (exact identity). "
               "Every call self-checks the Laplace invariant "
               "trace(G) = 0, validated at 7e-16 of max|G|; "
               "G<sub>uu</sub> is bit-consistent with the "
               "'gravity_gradient_rr' kernel route. 1 Eotvos = 1e-9 "
               "s<super>-2</super>."),
          PageBreak()]

# ================================================== Part II best practices
story += [para("Part II &mdash; Best practices", "h1"),
          para("Processing chain order", "h2"),
          para("The order below is not arbitrary; each step assumes the "
               "previous ones:")]
story += bullets([
    "<b>1.</b> Read GSM monthlies (shSeries.read, gz transparent, ICGEM "
    "1.0/2.0 incl. gfct).",
    "<b>2.</b> Replace C20 (and C30 in the accelerometer-degraded years) "
    "from TN-14; add degree-1 from TN-13 (addDegree1). Do this BEFORE "
    "filtering &mdash; filters see the stripes, not the low-degree "
    "completions.",
    "<b>3.</b> Bridge the GRACE&harr;GRACE-FO gap explicitly: dropNaN / "
    "select; never fit a climatology across NaNs.",
    "<b>4.</b> Filter: tvANS (or DDK for benchmarking; Gaussian as "
    "baseline).",
    "<b>5.</b> Remove GIA (trend work only).",
    "<b>6.</b> Fit climatology with tidal alias periods and "
    "ARCorrect=true; or form basin series with deconvolution.",
    "<b>7.</b> Synthesize maps / averages; propagate sigmas (analytic "
    "where derived, mcPropagate everywhere else)."])

story += [para("Choosing a filter", "h2")]
story += tbl(["situation", "recommendation"], [
    ["global EWH maps, standard use",
     "tvANS (Blocks=\"auto\"); Gaussian 300 km only as sanity baseline"],
    ["cross-paper comparability",
     "add DDK3-DDK6 variants via readDDK (binary Wbd, Nmax=96) + applyDDK"],
    ["strong order-dependent striping",
     "tvANS with VCEBands=[0 16 33 61] (v2.2)"],
    ["small basins (&lt; ~200,000 km&sup2;)",
     "expect leakage; taper kernels, consider basinDeconvolve with "
     "several neighboring kernels"],
    ["regional studies",
     "Slepian basis instead of filtering a global field (Part III, G4)"],
], [6.3 * cm, TEXT_W - 6.3 * cm])

story += [para("Uncertainty accounting &mdash; what is honest and what "
               "is approximate", "h2")]
story += bullets([
    "Formal per-coefficient sigmas (readers) &rarr; exact as published.",
    "tvANS posterior sigmas &rarr; exact for the unconstrained AND (v2.5) "
    "the constrained filter.",
    "Basin sigmas &rarr; noise part consistent with (banded) VCE and the "
    "exact (constrained) operator; deterministic-fit parameter "
    "uncertainty included since v2.5 (pointwise; epoch-correlated "
    "through the shared fit).",
    "Climatology sigmas &rarr; OLS/weighted/robust; AR(1) correction "
    "first-order with Kendall-corrected r<sub>1</sub> (v2.5; residual "
    "~3% optimism from regression absorption).",
    "GIA &rarr; model treated exact; use model spread.",
    "Multi-center combination &rarr; posterior sigmas are a LOWER bound "
    "(common mode invisible; see info.interCenterCorr).",
    "DDK-filtered fields &rarr; no analytic sigmas; use mcPropagate.",
    "Anything else &rarr; mcPropagate: independent-sigma or full SINEX "
    "covariance sampling validates every analytic formula empirically "
    "(the toolbox's own tests do exactly this)."])

story += [para("Validation habits the toolbox is built around", "h2")]
story += bullets([
    "Every numerical method was validated in Python (numpy/scipy) before "
    "the MATLAB implementation; the validation numbers are quoted in the "
    "help texts. Keep the habit for your own extensions.",
    "Ring-grid analysis&harr;synthesis roundtrips are exact to ~1e-14: "
    "use them as a self-test on your grids (Part III, G1 step 5).",
    "Real provider files beat synthetic fixtures: the test suite "
    "auto-discovers any *.gfc / *.snx(.gz) you drop into tests/test_data "
    "and sanity-checks them (large SINEX are streamed, v2.2).",
    "runAllTests is the gate: run it after every toolbox update; "
    "perf_log.csv accumulates machine-local timing history."])

story += [para("Common silent errors this toolbox guards against", "h2")]
story += tbl(["error", "guard"], [
    ["geodetic latitudes fed to SH math (~21 km bias)",
     "LatType=\"geodetic\" + conversion helpers (v2.2)"],
    ["hardcoded Love numbers", "kn is a required explicit input for EWH"],
    ["climatology across the 2017-2018 gap NaNs",
     "assertClean errors loudly; dropNaN/select"],
    ["trend maps without GIA", "removeGIA one-liner (v2.2)"],
    ["optimistic trend sigmas", "ARCorrect=true (v2.2)"],
    ["tidal alias leakage into trends", "Periods=[161/365.25, 3.66, 7.48]"],
    ["zero-variance rows collapsing the noise floor",
     "positive-diagonal floor + shLowLevel:buildNoiseCov:allZeroVariance"],
    ["8-column ICGEM 1.0 acos/asin ambiguity",
     "EIGEN-style period read; documented in README"],
    ["sine sectorals S<sub>nn</sub> blanked in difference triangles",
     "fixed v2.4.2; rendered-CData regression test"],
], [7.6 * cm, TEXT_W - 7.6 * cm])
story += [PageBreak()]

# ===================================================== Part III guides
story += [para("Part III &mdash; Step-by-step guides", "h1")]

story += [para("G1. Monthly EWH anomaly maps from ITSG files", "h2")]
story += code("""
% 1. read the monthly GSM series (gz ok, ICGEM 1.0/2.0)
ts = shSeries.read("data/ITSG-Grace2018_n96_*.gfc");

% 2. low-degree completion BEFORE filtering
tn14 = shLowLevel.readTN14("TN-14_C30_C20_SLR_GSFC.txt");
ts   = ts.applyTN14(tn14);                  % C20 (+C30 where flagged)
tn13 = shLowLevel.readTN13("TN-13_GEOC_GFZ_RL06.txt");
ts   = ts.addDegree1(tn13);

% 3. bridge the mission gap explicitly
ts = ts.dropNaN();

% 4. anomalies w.r.t. a mean field, then filter
ts  = ts - ts.mean();                       % anomalies
tsF = ts.filter("tvANS", Blocks="auto");    % or .gaussian(300)

% 5. EWH synthesis (Love numbers user-supplied!)
kn  = readmatrix("loadLoveNumbers.txt");    % your set, degrees 0..nmax
lat = -89.5:1:89.5; lon = 0.5:1:359.5;
ewh = tsF.at(24).synthesis(lat, lon, quantity="ewh", kn=kn);  % [m]

% self-test: analysis roundtrip on a ring grid should be ~1e-14
gChk = shCoefficients.analysis(ewh, lat, lon, tsF.nmax, ...
          quantity='ewh', kn=kn);
""")
story += [para("Interpretation: values are equivalent water height "
               "anomalies relative to the removed mean, in meters. The "
               "FFT synthesis path engages automatically on this uniform "
               "grid; check <font face='Courier'>tsF.at(24).history"
               "</font> for the full provenance chain.")]

story += [para("G2. Basin time series with uncertainties", "h2")]
story += code("""
idx = shLowLevel.shIndex(ts.nmax, MinDegree=2);
% kernel from a polygon [lat lon], buffered + tapered (v2.2)
amazon = readmatrix("amazon_polygon.txt");
[b, info] = shLowLevel.basinKernel(idx, amazon, BufferKm=150, TaperKm=350);

% filtered series + filter operator
[tsF, op] = ts.filter("tvANS", Blocks="auto", VCEBands=[0 16 33 61]);

% deconvolved basin average with 1-sigma (banded-VCE-consistent, v2.2)
[avg, out] = shLowLevel.basinDeconvolve(b, op);
shLowLevel.plotBasinSeries(op.tYears, avg(1,:)', out.sigma(1,:)', Units="cm");
""")
story += [para("The deconvolution undoes the filter's attenuation of the "
               "kernel; the sigma comes from the posterior covariance "
               "with per-band noise factors. For several basins pass a "
               "kernel matrix B &mdash; the joint solve also "
               "redistributes mutual leakage. (Edition-2 erratum fixed: "
               "the signature is <font face='Courier'>basinDeconvolve(B, "
               "op)</font>, epochs live in <font face='Courier'>"
               "op.tYears</font>, averages in the first output.)")]

story += [para("G3. Significance-masked trend maps (GIA- and "
               "AR(1)-honest)", "h2")]
story += code("""
gia = shCoefficients.read("ICE-6G_D_rates.gfc");       % rate [1/yr]
tsG = tsF.removeGIA(gia);
clim = tsG.climatology(Periods=[161/365.25, 3.66, 7.48], ...
                       ARCorrect=true);
tr   = clim.trend();                        % shCoefficients WITH sigmas
kn   = readmatrix("loadLoveNumbers.txt");
rate = tr.synthesis(lat, lon, quantity="ewh", kn=kn);         % [m/yr]
mc   = shLowLevel.mcPropagate(@(g) g.synthesis(lat, lon, ...
          quantity="ewh", kn=kn), tr, N=300, Seed=1);
sigma = mc.sigma;
mask = abs(rate) > 2 * sigma;               % 2-sigma significance
rate(~mask) = NaN; imagesc(lon, lat, rate);
""")
story += [para("Compare against a second GIA model (Caron) and show the "
               "difference as the model-spread contribution to the error "
               "budget &mdash; usually dominant over formal errors in "
               "formerly glaciated regions.")]

story += [para("G4. Regional analysis with Slepian functions", "h2")]
story += code("""
idx = shLowLevel.shIndex(60, MinDegree=0);
region = @(la,lo) double(la>-90 & la<-60);      % Antarctica-ish belt
[G, lam, infoS] = shLowLevel.slepianBasis(idx, region);
fprintf("Shannon %.1f -> keeping %d tapers\\n", infoS.shannon, numel(lam));

x  = shLowLevel.vecFromCS(g.C, g.S, idx);              % one month, idx order
a  = G' * x;                                    % Slepian coefficients
xR = G * a;                                     % regional part of x
[Cr, Sr] = shLowLevel.csFromVec(xR, idx);
gR = shCoefficients(Cr, Sr, GM=g.GM, R=g.R);
""")
story += [para("Estimating a (dimension &asymp; Shannon number) instead "
               "of all P Stokes coefficients is well-posed on regional "
               "data by construction; concentrations &lambda; near 1 "
               "tell you which tapers the region genuinely constrains.")]

story += [para("G5. Working with released covariances (SINEX)", "h2")]
story += code("""
idx = shLowLevel.shIndex(96, MinDegree=2);
% full read (monthly solutions to ~30 MB): covariance in idx order
snx = shLowLevel.readSINEX("ITSG-..._n96_2008-04.snx.gz", ...
                    Output="covariance", Index=idx);
N   = snx.M;                                % noise shape for tvANS
[tsF, op] = ts.filter("tvANS", NoiseCov=N, Blocks="off");

% multi-100-MB NEQ SINEX: stream just the estimates in seconds (v2.2)
est = shLowLevel.readSINEX("ITSG-..._operational_n96_2020-01.snx.gz", ...
                    Only="estimate");

% propagate the full covariance through ANY functional (v2.2)
out = shLowLevel.mcPropagate(@(gs) b' * shLowLevel.vecFromCS(gs.C,gs.S,idx), ...
                      g0, Cov=snx.M, Idx=idx, N=2000);
""")

story += [para("G6. Comparing against mascons", "h2")]
story += code("""
mas = shLowLevel.readMascon("GRCTellus.JPL...MSCNv03CRI.nc");
% mascon lats are GEODETIC -> one switch avoids the 21-km bias (v2.2)
kn = readmatrix("loadLoveNumbers.txt");
k  = find(abs(mas.epoch - 2010.29) < 0.05, 1);     % match an epoch
ewhSH = tsF.at(k).synthesis(mas.lat', mas.lon', ...
          quantity="ewh", kn=kn, LatType="geodetic");  % [m]
d = 100*ewhSH - mas.ewh(:,:,k);                    % mascon units: cm
fprintf("RMS diff %.2f cm\\n", sqrt(mean(d(:).^2, 'omitnan')));
""")
story += [para("Expect coastal and small-basin differences: mascons "
               "carry their own regularization. Compare basin-integrated "
               "series rather than pointwise fields for a fair "
               "statement; use the same GIA convention on both sides "
               "(JPL mascons ship GIA-corrected &mdash; remove GIA from "
               "the SH side too).")]

story += [para("G7. GNSS loading comparison (v2.3)", "h2")]
story += code("""
% residual monthly series (mean removed), filtered
ts  = ts - ts.mean();                       % anomalies
tsF = ts.filter("tvANS", Blocks="auto");

% load Love numbers k', h', l' - ONE consistent model, user-supplied;
% v2.5 reader: layouts, sparse pchip interpolation, CF-frame conversion
LN = shLowLevel.readLoveNumbers("loadLoveNumbers_PREM.txt", ...
    Columns="n k h l", MaxDegree=tsF.nmax, Interp="pchip", ...
    InFrame="CE", OutFrame="CF");               % match TN-13 (CF)
kn = LN.kn; hn = LN.hn; ln = LN.ln;

% GNSS stations (geodetic coordinates!)
sta = readmatrix("stations.txt");               % [lat lon]
T = tsF.nEpochs;
up = zeros(size(sta,1), T); no = up; ea = up;
for k = 1:T
    [up(:,k), no(:,k), ea(:,k)] = tsF.at(k).deformation( ...
        sta(:,1)', sta(:,2)', kn=kn, hn=hn, ln=ln, ...
        Mode="points", LatType="geodetic");
end
% compare up(i,:) against the GNSS vertical of station i (same epochs);
% typical GRACE-elastic amplitudes: mm to ~1 cm vertical, ~1/3 horizontal
""")
story += [para("Degree 1 is excluded by default; include it only when "
               "both the GRACE side (TN-13 in the CF frame) and the Love "
               "numbers use the same frame convention. The elastic "
               "prediction excludes GIA and poroelastic effects &mdash; "
               "detrend both series or remove GIA explicitly before "
               "comparing trends.")]

story += [para("G7b. Data management (v2.4.1, temporal catalogue "
               "v2.4.2)", "h2")]
story += code("""
shLowLevel.dataFolder("D:/geodata/shLowLevel");      % once; persists (getpref)
shLowLevel.fetchDDK(1:8);                     % all released DDK filters
W5 = shLowLevel.readDDK("DDK5");              % by name, from the data folder
T = shLowLevel.listICGEM();                   % 180+ static models as a table
f = shLowLevel.fetchICGEM("EGM2008");         % -> dataFolder/icgem/EGM2008.gfc
g = shCoefficients.read(f);
shLowLevel.fetchITSG(2010:2016);              % -> dataFolder/itsg_series
shLowLevel.fetchITSG("2008-04", Product="daily"); % Kalman daily n40 (v2.4.1)

% FULL temporal catalogue: all ~70+ series, all centers
Tt = shLowLevel.listICGEM(Type="temporal");
% fetch a WHOLE monthly series (v3.1.1) - one row of that catalogue
fs = shLowLevel.fetchICGEM(17, Type="temporal");
ts = shSeries.fromFolder(fileparts(fs(1)));   % straight into a series
""")
story += [para("The ICGEM static catalogue is parsed from icgem.gfz.de "
               "(fixture-tested offline). The temporal catalogue lists "
               "every series of every center &mdash; groups 01_GRACE, "
               "02_COST-G, 03_other, 04_SLR &mdash; with columns group, "
               "center, series, path, url and zip. Since v3.1.1 a "
               "temporal row is not just browsable but fetchable: "
               "<font face='Courier'>fetchICGEM(idx, Type=\"temporal\")</font> "
               "downloads the whole series into "
               "<font face='Courier'>&lt;dataFolder&gt;/icgem/series/"
               "&lt;group&gt;_&lt;center&gt;_&lt;series&gt;/</font>, which "
               "<font face='Courier'>shSeries.fromFolder</font> and "
               "<font face='Courier'>shLowLevel.standardChain</font> consume "
               "directly. For ITSG monthlies "
               "<font face='Courier'>shLowLevel.fetchITSG</font> remains the "
               "convenient route.")]

story += [para("Mode: one request, not three hundred", "h3")]
story += [para("The default <font face='Courier'>Mode=\"auto\"</font> takes "
               "the server's whole-series ZIP in a SINGLE request and "
               "unpacks it. This is not an optimization but a courtesy "
               "requirement: icgem.gfz.de rate-limits, and fetching a "
               "283-file series one file at a time earns HTTP 429 and a "
               "tarpit (field-observed as &quot;too many connections&quot;). "
               "If the archive is unavailable the fetcher falls back to "
               "per-file mode automatically &mdash; resumable, throttled "
               "with <font face='Courier'>Pause=</font> and "
               "<font face='Courier'>Retries=</font> exponential backoff, "
               "each file verified by parse before it replaces anything. "
               "<font face='Courier'>Mode=\"archive\"|\"files\"</font> force "
               "either path; <font face='Courier'>info.mode</font> reports "
               "which one actually ran, so a script can tell a fresh "
               "download from a no-op. <font face='Courier'>Files=</font> "
               "filters the series, and "
               "<font face='Courier'>FileList=</font> injects a catalogue "
               "table for offline mirrors and subsets.")]
story += code("""
% forced per-file, filtered to a window, gentle on the server
fs = shLowLevel.fetchICGEM(17, Type = "temporal", Mode = "files", ...
        Files = "*2008*.gfc", Pause = 3, Retries = 3);
[ts, rep] = shLowLevel.standardChain(fileparts(fs(1)), Filter = "DDK3");
""")

story += [para("G8. Multi-center combination (v2.4)", "h2")]
story += code("""
shLowLevel.fetchITSG(2010:2016);              % download once (websave)
tsI = shSeries.fromFolder("ITSG/");    % GSM level, same nmax
tsC = shSeries.fromFolder("CSR/");
tsG = shSeries.fromFolder("GFZ/");
[tsComb, info] = shLowLevel.combineCenters({tsI, tsC, tsG}, Robust=true);

figure; plot(info.epochs, info.weights');     % weight time series
legend(info.centers); ylabel("relative weight");
figure; imagesc(info.interCenterCorr);        % honesty diagnostic
% THEN degree-1 / C20-C30 replacements, once:
tsComb = tsComb.applyTN14("TN-14.txt");
tsComb = tsComb.addDegree1("TN-13.txt");
""")
story += [para("Watch the weight series for regime changes "
               "(single-accelerometer era, GRACE-FO transition) and "
               "treat posterior sigmas as lower bounds &mdash; the "
               "common mode is invisible.")]

story += [para("G9. Sea-level fingerprint (v2.4)", "h2")]
story += code("""
idx = shLowLevel.shIndex(96, MinDegree=0);
LN = readmatrix("loadLoveNumbers.txt"); kn = LN(:,2); hn = LN(:,3);
greenland = @(la,lo) -280 * inGreenland(la, lo);   % kg/m^2/yr, say
oceanF = @(la,lo) oceanMaskFun(la, lo);            % 1 over ocean
[S, grid, info] = shLowLevel.seaLevelFingerprint(greenland, oceanF, idx, ...
    kn=kn, hn=hn);
shLowLevel.plotSHMap(info.S2D / info.eustatic, grid.latDeg, grid.lonDeg, ...
    Units="S / eustatic", Title="Greenland fingerprint");
""")
story += [para("Plot S normalized by the eustatic mean: near-field "
               "values drop below zero, far-field plateaus around "
               "1.1&ndash;1.3 &mdash; if not, check the ocean mask and "
               "that idx uses MinDegree = 0.")]

story += [para("G10. The standard chain in one call", "h2")]
story += code("""
% the whole post-processing chain, in the ONE correct order
gia = shCoefficients.read("ICE-6G_D_trend.gfc");   % a TREND field [1/yr]
[ts, rep] = shLowLevel.standardChain("D:/grace/monthly", ...
    TN14File = "TN-14_C30_C20_SLR_GSFC.txt", ...   % TN14 itself is a FLAG
    Degree1 = "CSR", ...                           % provider, not a path
    GIA = gia, Filter = "DDK3");
disp(rep.steps')          % what was applied, in order, with provenance
""")
story += [para("Order is not a matter of taste here. The SLR C20/C30 "
               "replacement and the degree-1 restoration must happen on "
               "the UNFILTERED coefficients &mdash; a filter mixes degrees, "
               "so replacing a coefficient afterwards inserts an "
               "unfiltered value into a filtered field and the low "
               "degrees no longer mean anything consistent. GIA is a "
               "trend in the same (corrected) coefficients, so it comes "
               "next; filtering is last. "
               "<font face='Courier'>standardChain</font> encodes exactly "
               "that order and returns REP, a provenance record naming "
               "every step, the files it used and the toolbox version "
               "that ran &mdash; the thing you paste into a paper's data "
               "section. Every stage is optional: leave "
               "<font face='Courier'>GIA=</font> out and the report says "
               "so, rather than silently skipping a step you assumed had "
               "happened.")]
story += [para("Two contracts worth reading before the first call: "
               "<font face='Courier'>TN14</font> is a logical FLAG and the "
               "path goes in <font face='Courier'>TN14File</font>, while "
               "<font face='Courier'>Degree1</font> names the TN-13 "
               "provider (\"GFZ\" | \"CSR\" | \"JPL\" | \"none\") with the "
               "path in <font face='Courier'>Degree1File</font>; and "
               "<font face='Courier'>GIA</font> takes an "
               "<font face='Courier'>shCoefficients</font> TREND field in "
               "[1/yr], not a filename &mdash; it is applied as "
               "(t &minus; GIAEpoch) &times; GIA. Left at their defaults "
               "the files come from the persistent data folder that "
               "<font face='Courier'>setup_shAnalysis</font> populates.")]

story += [para("G11. Designing your own anisotropic filter", "h2")]
story += code("""
% a DDK-class filter built from YOUR error covariance
g  = shCoefficients.read("ITSG-Grace2018_n60_2008-04.gfc");
W  = shLowLevel.designFilter(g.sigmaC, g.sigmaS, Kaula = 1e-5);
[Cf, Sf] = shLowLevel.applyDDK(g.C, g.S, W);   % same block format as readDDK
""")
story += [para("The released DDK filters are built from a specific "
               "center's normal equations. If you have your own sigmas "
               "or a released covariance, "
               "<font face='Courier'>designFilter</font> builds the same "
               "kind of operator, W = (N + a S<sup>-1</sup>)<sup>-1</sup> N, "
               "in the block format "
               "<font face='Courier'>readDDK</font> produces &mdash; so it "
               "drops into <font face='Courier'>applyDDK</font>, "
               "<font face='Courier'>g.applyDDK</font> and "
               "<font face='Courier'>standardChain(Filter=W)</font> "
               "unchanged. Note that <font face='Courier'>Noise=</font> "
               "expects a COVARIANCE, not sigmas: passing standard "
               "deviations where a variance is expected is the classic "
               "way to get a filter that looks plausible and smooths by "
               "the square root of what you intended.")]

story += [para("G12. Correcting leakage", "h2")]
story += code("""
% forward modelling: only the FILTER has to be known
ewh = g.synthesis(-89:89, 0:359, quantity = "ewh", kn = kn);
[m, info] = shLowLevel.leakageCorrect(ewh, -89:89, 0:359, ...
    Filter = "gauss300", Mask = basinMask, Gain = 2);

% or scaling factors, if you trust a model's spatial PATTERN
k = shLowLevel.gridScaling(modelEWH, -89:89, 0:359, Filter = "gauss300");
corrected = k .* ewh;                 % NaN where the model is blind
""")
story += [para("A filter removes real signal and spreads the rest across "
               "basin boundaries; a 6&deg; disc under a 500 km Gaussian "
               "keeps about half its peak. The two corrections differ in "
               "what you have to BELIEVE, which is how to choose between "
               "them. Forward modelling "
               "(<font face='Courier'>leakageCorrect</font>) needs only "
               "the filter and, usefully, a mask saying where mass can "
               "exist &mdash; it then has to explain the observation with "
               "mass inside that region, so leakage outside is removed "
               "rather than redistributed. Scaling factors "
               "(<font face='Courier'>gridScaling</font>) need a model "
               "series, and inherit the correctness of its spatial "
               "pattern: the factors are invariant to the model&rsquo;s "
               "amplitude but not to its shape.")]
story += [para("Two behaviours worth knowing before you read the output. "
               "Convergence is judged on the change of the SOLUTION, not "
               "on the residual: with a mask the problem is inconsistent, "
               "so the residual settles at a floor while the solution is "
               "converged. And <font face='Courier'>gridScaling</font> "
               "returns NaN where the model carries no signal instead of "
               "a ratio of two numerical zeros &mdash; a model that does "
               "not reach your basin shows up as missing coverage rather "
               "than as noise multiplying your data.")]

story += [para("G13. What real provider files look like", "h2")]
story += [para("Format specifications describe what files ought to "
               "contain. The reader is written against what they "
               "actually contain, and three separate field failures are "
               "worth knowing about, because they are silent &mdash; you "
               "get numbers, just not the right ones.")]
story += bullets([
    "<b>FORTRAN D-exponents.</b> Providers switch from "
    "<font face='Courier'>1.23E-10</font> to "
    "<font face='Courier'>1.23D-10</font> partway through large files "
    "(EIGEN-6C4 does it above degree ~370). MATLAB's "
    "<font face='Courier'>str2double</font> returns NaN for those, so a "
    "line-by-line parser quietly corrupts the high degrees of exactly "
    "the models where the high degrees are the point. The reader "
    "normalizes D/d to E/e before conversion.",
    "<b>ICGEM 2.0 column order.</b> The real-world layout of an "
    "<font face='Courier'>acos</font>/<font face='Courier'>asin</font> "
    "line is <font face='Courier'>... sigC sigS t0 t1 period</font> "
    "&mdash; period LAST, verified against CNES/GRGS files. A "
    "period-first assumption mis-parses every such file and produces a "
    "time-variable model whose seasonal terms are wrong without ever "
    "raising an error.",
    "<b>Ragged record groups.</b> EIGEN-5S and 5C carry a single "
    "<font face='Courier'>gfc</font> line with a trailing epoch among "
    "thousands of uniform ones. Bulk "
    "<font face='Courier'>sscanf(txt, '%f', [nc Inf])</font> does not "
    "complain about that &mdash; it silently re-flows the whole group by "
    "one column and asks for a coefficient matrix indexed by a date. "
    "The reader checks rows-out == lines-in per group, subgroups by "
    "numeric width when ragged, and keeps an n/m sanity net.",
])
story += [para("The practical rule: when a new provider or a new release "
               "enters your pipeline, read one file and check "
               "<font face='Courier'>g.C(3,1)</font> and the highest "
               "degrees against the provider's own published values "
               "before you trust a thousand of them. Bulk parsing is "
               "fast &mdash; EIGEN-6C4 at 177.7 MB and degree 2190 in "
               "about 5 s, a 73.6 MB GRGS mean field with 674k variable "
               "terms in about 7 s &mdash; but fast and correct are "
               "different properties.")]
story += [PageBreak()]

# ======================================================== demo gallery
story += [para("Demo gallery", "h1"),
          para("<font face='Courier'>demo_shAnalysis</font> is a case "
               "registry: <font face='Courier'>demo_shAnalysis(\"list\")"
               "</font> prints the table below, <font face='Courier'>"
               "demo_shAnalysis(\"all\")</font> runs everything, "
               "<font face='Courier'>demo_shAnalysis([\"D04\",\"D15\"])"
               "</font> runs a selection, and Visible=false keeps figures "
               "hidden for smoke tests. Real files in tests/test_data "
               "(ITSG month, DDK3 Wbd) are used when present; every case "
               "falls back to synthetic data. Love numbers in the demos "
               "are SYNTHETIC &mdash; real work requires a real loading "
               "model.")]
story += tbl(["ID", "Case", "Exercises"], [
    ["D01", "Read &amp; spectral diagnostics",
     "shReadGFC, shDegreeRMS, triangle, spectrum"],
    ["D02", "Synthesis quantities &amp; maps",
     "shSynthesis, kernelFactors, plotSHMap, Height"],
    ["D03", "Normal field &rarr; geoid",
     "normalFieldCS, subtractNormalField, toReference"],
    ["D04", "Filter comparison",
     "gaussian, fan, destripe, DDK, triangle diff"],
    ["D05", "Climatology &amp; basin series",
     "climatology, basinKernel, plotBasinSeries"],
    ["D06", "tvANS pipeline", "buildNoiseCov, tvANSFilter, basinDeconvolve"],
    ["D07", "Analysis (inverse)", "shAnalysisGrid, Kaula"],
    ["D08", "Basin tools", "basinKernel taper, deconvolve vs scaling"],
    ["D09", "Uncertainty", "errorMap, mcPropagate, plotCovariance"],
    ["D10", "Load deformation", "shSynthesisDeformation, Mode=points"],
    ["D11", "Gradient tensor", "shSynthesisGradientTensor"],
    ["D12", "EOF analysis", "eofAnalysis"],
    ["D13", "Trend breakpoints", "trendBreaks, F-test"],
    ["D14", "Multi-center combination", "combineCenters, VCE weights"],
    ["D15", "Sea-level fingerprint", "seaLevelFingerprint, evalMask"],
    ["D16", "Export", "writeGrid netCDF, writeAnimation MP4"],
], [1.1 * cm, 5.6 * cm, TEXT_W - 6.7 * cm])
story += [para("<b>Note.</b> Provenance of the figures below: rendered "
               "with the toolbox's Python validation port &mdash; the "
               "same algorithms that gate every MATLAB implementation, "
               "numerically equivalent by construction. D01/D02 use the "
               "REAL ITSG files shipped in tests/test_data (and D04 the "
               "real DDK3 Wbd binary); the remaining figures use "
               "synthetic data with known ground truth, because recovery "
               "demonstrations (VCE factors, EOF modes, F-test "
               "calibration) require a truth to recover. "
               "<font face='Courier'>shLowLevel.fetchITSG(years)</font> extends "
               "the real-data demos to full monthly series. The MATLAB "
               "demo cases produce the interactive originals.")]

story += fig("d01.png",
             "<b>D01 &mdash; REAL DATA:</b> the shipped ITSG-Grace2018 "
             "2008-04 month (rendered by the rebuilt Python port, v2.5). "
             "Coefficient triangle (apex up; C00/C20 removed for "
             "display) and the degree-amplitude spectrum with the REAL "
             "formal errors from the gfc file. At this n60 truncation "
             "the formal errors stay ~3 orders below the signal - the "
             "signal/error crossover (the resolution limit) lies beyond "
             "the truncation degree; Kaula decay as the yardstick.")
story += fig("f-002.png",
             "<b>D02 &mdash; REAL DATA.</b> Left: the disturbing geoid of "
             "2008-04 after subtracting the WGS84 normal field computed "
             "from defining constants and rescaled to the ITSG (GM, R) "
             "&mdash; the classic &plusmn;100 m map (Indian Ocean low, "
             "New Guinea high). Right: the real 2008-2025 mass change "
             "(GRACE-FO minus GRACE, Gaussian 350 km, EWH): Greenland "
             "and West Antarctic mass loss stand out; residual "
             "meridional stripes at low latitudes are the real GRACE "
             "error structure. Honest caveats: synthetic Love numbers, "
             "no GIA correction, degree 1 absent &mdash; patterns are "
             "real, amplitudes approximate.")
story += fig("d04_gains.png",
             "<b>D04 &mdash; Filter gain over the (degree, order) "
             "triangle.</b> Gaussian: isotropic in degree. Fan (Han): "
             "the order Gaussian additionally damps high orders &mdash; "
             "the striping direction. DDK3 (diagonal of the real Wbd "
             "blocks): anisotropic and order-dependent &mdash; "
             "near-zonal signal survives to n~35 while high orders are "
             "suppressed hard; exactly the structure a GRACE error "
             "covariance demands.")
story += fig("d04_diff.png",
             "<b>D04b &mdash; REAL DATA, v2.4.2 fix:</b> the signed "
             "difference triangle (Gaussian 350 &minus; raw) on the real "
             "GRACE-FO &minus; GRACE field at n40 &mdash; the removed "
             "signal concentrates at the low orders where the real "
             "stripes and secular signal live. The sine wing now spans "
             "the full n &ge; m &ge; 1 support including the sectorals "
             "S<sub>nn</sub>: before v2.4.2 the difference mode blanked "
             "the sectoral edge and every diff triangle rendered the "
             "left wing one column narrower than the right (regression-"
             "tested on the rendered CData since).", width=330)
story += fig("f-006.png",
             "<b>D05/D06 &mdash; The standard basin figure</b> from "
             "plotBasinSeries: 1-sigma band, grey mission-gap patch, "
             "AR(1)-corrected trend annotation. With tvANS + "
             "basinDeconvolve the band comes from propagated posterior "
             "covariances instead of an empirical scatter.")
story += fig("f-008.png",
             "<b>D15 &mdash; Sea-level fingerprint</b> of a polar-cap "
             "melt, normalized by the eustatic mean: the near field "
             "FALLS (blue; lost self-attraction plus crustal rebound), "
             "the far field exceeds eustatic by ~20-30%. Mass conserved "
             "to machine precision &mdash; the printed residual is the "
             "built-in check.")
story += fig("f-010.png",
             "<b>D14 &mdash; Multi-center VCE:</b> per-(center, month) "
             "weights track the time-varying noise levels; the "
             "inter-center residual correlation matrix is the honesty "
             "diagnostic &mdash; VCE assumes independence, GRACE centers "
             "share data and background models.")
story += fig("f-012.png",
             "<b>D10 &mdash; Elastic load deformation:</b> up, north, "
             "east from the same residual field. Horizontals are ~1/3 "
             "of the vertical (validation run: 0.38) and require the "
             "exact SH gradient (legendreALFDeriv); east is undefined "
             "at the poles.")
story += fig("f-014.png",
             "<b>D11 &mdash; All six gradient tensor components</b> at "
             "250 km in Eotvos. The Laplace invariant G<sub>nn</sub> + "
             "G<sub>ee</sub> + G<sub>uu</sub> = 0 holds at roundoff "
             "(4e-15 of max|G| here) &mdash; a full-chain check of first "
             "and second Legendre derivatives, trig derivatives and "
             "frame terms at once.")
story += fig("f-016.png",
             "<b>D13 &mdash; Trend breakpoint:</b> continuous "
             "piecewise-linear hinge model vs the no-break fit; the "
             "F-test (p via betainc, base MATLAB) quantifies the break. "
             "Null rejection is calibrated at the nominal level (5.1% "
             "at 5%).")
story += [PageBreak()]

# ============================================================= appendix

# ===================================================== Part IV API reference
import sys as _sys
if "/home/claude" not in _sys.path:
    _sys.path.insert(0, "/home/claude")
from api_data import API as _API


def _esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;")
             .replace(">", "&gt;"))


_INSTVAR = {"shCoefficients": "g", "shSeries": "ts", "shClimatology": "clim"}


def api_entry(e, qual):
    out = [para("<b>%s</b>" % _esc(qual), "h3")]
    # call-form signature
    if e["kind"] == "method" and not e.get("static"):
        callee = _INSTVAR[e["cls"]] + "." + e["name"]
        ins = [a for a in e["inputs"] if a not in ("obj",)]
    elif e["kind"] == "method":
        callee = e["cls"] + "." + e["name"]
        ins = e["inputs"]
    else:
        callee = qual
        ins = e["inputs"]
    nvargs = [a for a in e.get("args", []) if a["nv"]]
    ins = [a for a in ins if a != "opts"] +           (["Name=Value"] if nvargs else [])
    outs = e["outputs"]
    sig = (("[" + ", ".join(outs) + "] = ") if len(outs) > 1 else
           (outs[0] + " = " if outs else "")) + callee +           "(" + ", ".join(ins) + ")"
    out.append(para("<font face='Courier' size='8.2'>%s</font>" % _esc(sig)))
    d = (e.get("h1", "") + " " + e.get("desc", "")).strip()
    if len(d) > 520:
        d = d[:517] + "..."
    if d:
        out.append(para(_esc(d)))
    # inputs table from the arguments block (ground truth)
    rows = []
    for a in e.get("args", []):
        typ = a["cls"] or "&mdash;"
        if a.get("allowed"):
            typ += " &isin; " + _esc(a["allowed"])
        rows.append([a["name"], a["size"] or "&mdash;", typ,
                     _esc(a["default"]) if a["default"] else "required",
                     "name-value" if a["nv"] else "positional"])
    if not rows:
        rows = [[i, "&mdash;", "&mdash;", "required", "positional"]
                for i in ins if i != "Name=Value"]
    if rows:
        out += tbl(["input", "size", "type", "default", "kind"], rows,
                   [0.16 * TEXT_W, 0.11 * TEXT_W, 0.33 * TEXT_W,
                    0.24 * TEXT_W, 0.16 * TEXT_W])
    # outputs
    orow = e.get("outdocs") or [[o, ""] for o in outs]
    if orow:
        out += tbl(["output", "description"],
                   [[_esc(a), _esc(b)] for a, b in orow],
                   [0.18 * TEXT_W, 0.82 * TEXT_W])
    if e.get("example"):
        out.append(para("<i>Example</i>"))
        out += code(e["example"])
    return out


story += [PageBreak(), para("Part IV &mdash; API reference", "h1"),
          para("Generated from the source: the <b>arguments blocks are the "
               "ground truth</b> for input names, dimensions, classes and "
               "defaults (default 'required' = mandatory positional input); "
               "name-value options are passed as Name=Value. Conventions "
               "throughout: coefficient triangles are addressed C(n+1, m+1); "
               "vectors follow shLowLevel.shIndex ordering (MinDegree = 2 unless "
               "stated); latitudes are geocentric degrees; radii in km; "
               "epochs in decimal years. Every entry closes with a "
               "real-data example from the canonical GRACE chain; all "
               "examples are lint- and API-gated at build time.")]

for _cls in _API["classes"]:
    story += [para(_cls["name"] + " class", "h2"),
              para(_esc((_cls["h1"] + " " + _cls["desc"]).strip()[:600]))]
    if _cls["props"]:
        _rows = []
        for _p in _cls["props"]:
            _nm, _sz, _ty, _df, _cm = (_p + [""] * 5)[:5]
            _rows.append([_esc(_nm), _esc(_sz) or "&mdash;",
                          _esc(_ty) or "&mdash;",
                          _esc(_df) if _df else "&mdash;", _esc(_cm)])
        story += tbl(["property", "size", "type", "default", "notes"],
                     _rows, [0.18 * TEXT_W, 0.12 * TEXT_W, 0.2 * TEXT_W,
                             0.3 * TEXT_W, 0.2 * TEXT_W])
    for _m in _cls["methods"]:
        story += api_entry(_m, _cls["name"] + "." + _m["name"] +
                           (" (static)" if _m.get("static") else ""))

_CATS = [
 ("Readers, writers &amp; data management",
  ["shReadGFC", "shEvalGFCT", "readTN13", "readTN14", "readSINEX",
   "readDDK", "readMascon", "readLoveNumbers", "listICGEM", "fetchITSG",
   "fetchICGEM", "fetchDDK", "fetchTN", "dataFolder", "ddkNames",
   "parseGraceFilename", "icgemDate2Year", "writeGFC", "writeGrid",
   "writeAnimation"]),
 ("Legendre, spectral tools &amp; indexing",
  ["legendreALF", "legendreALFDeriv", "legendreCached", "ylm", "shIndex",
   "vecFromCS", "csFromVec", "shDegreeRMS", "shOrderRMS",
   "shSpectralCrossover", "gaussLegendre"]),
 ("Synthesis &amp; analysis",
  ["shSynthesis", "synthesisMatrix", "shSynthesisDeformation",
   "shSynthesisGradientTensor", "shAnalysisGrid", "kernelFactors",
   "normalFieldCS", "rescaleGMR", "geodetic2geocentric",
   "geocentric2geodetic", "evalMask"]),
 ("Filtering &amp; operators",
  ["shGaussianWeights", "shGaussianFilter", "shFanFilter", "shDestripe",
   "applyDDK", "tvANSFilter", "buildNoiseCov", "buildSignalCov",
   "vceRescale", "opApply", "resolutionMap"]),
 ("Basins, uncertainty &amp; analysis tools",
  ["basinKernel", "basinDeconvolve", "basinScaling",
   "fitDeterministicModel", "mcPropagate", "errorMap", "combineCenters",
   "eofAnalysis", "slepianBasis", "seaLevelFingerprint", "removeGIA",
   "pctile"]),
 ("Comparison &amp; validation metrics (v2.6.0)",
  ["compareSolutions", "compareSeries", "diffSpectrum", "spatialStats",
   "nashSutcliffe", "effectiveCorr", "threeCorneredHat",
   "taylorDiagram"]),
 ("Plotting",
  ["plotSHMap", "plotSHSpectrum", "plotSHCoeffTriangle",
   "plotBasinSeries", "plotCovariance"]),
]
_bynm = {e["name"]: e for e in _API["shLowLevel"]}
_used = set()
story += [para("Package +shLowLevel", "h2")]
for _t, _names in _CATS:
    story += [para(_t, "h3"), Spacer(1, 2)]
    for _n in _names:
        if _n in _bynm:
            story += api_entry(_bynm[_n], "shLowLevel." + _n)
            _used.add(_n)
_rest = [n for n in sorted(_bynm) if n not in _used]
if _rest:
    story += [para("Utilities", "h3")]
    for _n in _rest:
        story += api_entry(_bynm[_n], "shLowLevel." + _n)

story += [para("Toolbox root", "h2")]
for _e in _API["root"]:
    story += api_entry(_e, _e["name"])


story += [para("Appendix", "h1"),
          para("A. Error-identifier conventions", "h2"),
          para("Namespaces: <font face='Courier'>shLowLevel:&lt;function&gt;:"
               "&lt;what&gt;</font> for package functions, "
               "<font face='Courier'>shCoefficients:/shSeries:/"
               "shClimatology:</font> for class methods. Every guard has "
               "a unit test; every message names the offending quantity. "
               "Instrumentation asserts (nonFiniteVCE, nonFiniteEig, "
               "nonFiniteOperator, allZeroVariance) remain in the code "
               "as a permanent safety net after the v2.1 debugging "
               "campaign."),
          para("B. Known limitations (documented, tested where "
               "applicable)", "h2")]
story += bullets([
    "TN-13 verified against the real GFZ, CSR and JPL RL06.3 files "
    "(v2.5; 256 paired months each, cross-provider C10 correlation "
    "0.995); SINEX against ITSG monthly files and synthetic fixtures "
    "&mdash; COST-G/other SINEX dialects may need parser tweaks (drop "
    "files into tests/test_data, the discovery tests pick them up).",
    "DDK: all eight released Wbd binaries parsed and gain-ordering "
    "validated (v2.5); a discovery test re-checks whatever subset is "
    "present locally. Other BIN dialects untested.",
    "readMascon validated against a synthetic JPL-style fixture; real "
    "product files were not openly reachable at build time (auth/503) "
    "&mdash; verify against yours (units, GIA convention, CRI flag).",
    "GIA models treated as exact; AR(1) correction is first-order "
    "(Kendall-corrected r<sub>1</sub> since v2.5, ~3% residual optimism).",
    "Multi-center combination: common-mode errors invisible &rarr; "
    "posterior sigmas a lower bound.",
    "Basin sigma: deterministic-part errors correlated across epochs "
    "(pointwise 1-sigma reported).",
    "Streaming gz SINEX needs the JVM (fallback: gunzip to temp).",
    "ICGEM catalogue parsers are fixture-tested against real page "
    "captures (2026-08-07); a future site redesign would need refreshed "
    "fixtures (shLowLevel:listICGEM:parseFailed errors loudly)."])
story += [para("C. License &amp; provenance", "h2"),
          para("The toolbox is released under the MIT license (LICENSE "
               "in the package root; copyright Matthias Weigelt). The "
               "bundled DDK3 Wbd file carries its own MIT license from "
               "its upstream repository; ITSG/GFZ files in "
               "tests/test_data remain the property of their providers "
               "and are included as test fixtures only."),
          para("All code, tests, documentation and this guide: developed "
               "by Matthias Weigelt with the help of Claude (Fable 5), "
               "2026-08-07 (v1 through v2.4.2 in working sessions the "
               "same week). Every numerical method was independently "
               "validated in Python before MATLAB implementation; MATLAB "
               "execution itself happens on the user's side with "
               "runAllTests as the acceptance gate (118/118 green at "
               "v2.1; 134 at v2.2.2; 167 at v2.4.1; 168 at v2.4.2; 179 "
               "at v2.5 &mdash; v2.5 adds the exact constrained "
               "posterior, the basin deterministic-sigma term, the "
               "Kendall AR(1) correction, shLowLevel.readLoveNumbers, "
               "setup_shAnalysis/shLowLevel.fetchTN with verified staged "
               "downloads, real CSR/JPL TN-13 + GSFC TN-14 fixtures "
               "and the DDK-pack validation).")]

# ---------------------------------------------------------------- build
doc.build(story)
print("built %s (%d snippets registered)" % (OUT, len(SNIPPETS)))
