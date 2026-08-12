"""api_extract.py - build the shAnalysis API reference from the source.

Parses every public function (+shLowLevel/*.m, toolbox root) and every public
method/property of the three classes. The arguments blocks are the
ground truth for input names, dimensions, types and defaults; the help
text supplies descriptions, output tables and examples. Writes
/home/claude/api_data.py for make_guide.py to render.

Developed by Matthias Weigelt with the help of Claude (Fable 5),
2026-08-07 (v2.5).
"""
import glob
import os
import re

ROOT = "/home/claude/shx_build/shAnalysis"


def _join_continuations(lines):
    out, buf = [], ""
    for L in lines:
        s = L.rstrip()
        core = s.split("%")[0] if not s.lstrip().startswith("%") else s
        if buf:
            s = buf + " " + s.lstrip()
            core = s.split("%")[0]
            buf = ""
        if core.rstrip().endswith("..."):
            buf = core.rstrip()[:-3].rstrip()
            continue
        out.append(s)
    if buf:
        out.append(buf)
    return out


def parse_arguments_block(lines, i):
    """Parse one arguments...end block starting at index i ('arguments')."""
    args = []
    i += 1
    while i < len(lines):
        s = lines[i].strip()
        if s == "end":
            break
        if s and not s.startswith("%"):
            s = s.split("%")[0].strip()        # drop %#ok<...> etc.
            m = re.match(
                r"([\w.]+)\s*"                 # name (may be opts.X)
                r"(\([^)]*\))?\s*"             # (size)
                r"([a-zA-Z]\w*)?\s*"           # class
                r"(\{.*\})?\s*"                # {validators} (may nest cells)
                r"(?:=\s*(.+))?$", s)
            if m:
                name, size, cls, vals, dflt = m.groups()
                nv = name.startswith("opts.")
                entry = dict(
                    name=name.split(".", 1)[1] if nv else name,
                    nv=nv,
                    size=(size or "").strip("()").replace(",", " x "),
                    cls=cls or "",
                    validators=(vals or "").strip("{}"),
                    default=(dflt or "").strip())
                mm = re.search(r"mustBeMember\([^,]+,\s*(\[[^\]]*\]|\{[^}]*\})",
                               entry["validators"])
                if mm:
                    entry["allowed"] = mm.group(1)
                args.append(entry)
        i += 1
    return args, i


def parse_help(block):
    """From a contiguous %-comment block: H1, description, outputs, example."""
    lines = [re.sub(r"^\s*%\s?", "", L.rstrip()) for L in block]
    h1 = lines[0].strip() if lines else ""
    h1 = re.sub(r"^[A-Z_0-9]+\s+", "", h1)     # strip NAME prefix
    desc, i = [], 1
    while i < len(lines):
        s = lines[i].strip()
        if s == "" and desc:
            break
        if re.match(r"(Inputs?|Options?|Outputs?|Example|Recognized|Notes?"
                    r"|Syntax|Usage)\b", s):
            break
        if s:
            desc.append(s)
        i += 1
    # outputs section: "  Outputs" then indented "name   description"
    outputs = []
    m = re.search(r"^\s*Outputs?\s*$", "\n".join(lines), re.M)
    if m:
        idx = [k for k, L in enumerate(lines)
               if re.match(r"^\s*Outputs?\s*$", L)][0] + 1
        while idx < len(lines):
            s = lines[idx]
            if s.strip() == "" or re.match(
                    r"^\s*(Inputs?|Options?|Example|Notes?|Developed)", s):
                break
            mo = re.match(r"^\s{1,8}([\w.]+(?:,\s*[\w.]+)*)\s{2,}(.+)$", s)
            if mo:
                outputs.append([mo.group(1), mo.group(2).strip()])
            elif outputs:
                outputs[-1][1] += " " + s.strip()
            idx += 1
    # example block
    example = ""
    idx = None
    for k, L in enumerate(lines):
        if re.match(r"^\s*Examples?\b", L.strip()):
            idx = k + 1
            break
    if idx is not None:
        ex = []
        while idx < len(lines):
            s = lines[idx]
            if re.match(r"^\s*(Developed|Outputs?|Inputs?|Options?"
                        r"|Error identifiers)", s) or (
                    s.strip() == "" and ex and lines[idx - 1].strip() == ""):
                break
            ex.append(s.rstrip())
            idx += 1
        example = "\n".join(x for x in ex).strip("\n")
    # Inputs/Options sections: "name (size/type) description", continuation
    # lines indented deeper. Collected for merging into arguments-block args.
    iodocs = {}
    for hdr in ("Inputs?", "Options?"):
        for k, L in enumerate(lines):
            if re.match(r"^\s*" + hdr + r"\s*$", L):
                idx = k + 1
                while idx < len(lines):
                    t = lines[idx]
                    if t.strip() == "" or re.match(
                            r"^\s*(Inputs?|Options?|Outputs?|Example|Notes?"
                            r"|Developed|Error identifiers)", t):
                        break
                    mo = re.match(r"^\s{1,8}([\w.]+)\s+(.*)$", t)
                    if mo and not t.startswith("        "):
                        nm = mo.group(1)
                        rest = mo.group(2).strip()
                        sz = ""
                        ms = re.match(r"\(([^)]{1,24})\)\s*(.*)$", rest)
                        if ms:
                            sz = ms.group(1); rest = ms.group(2).strip()
                        iodocs[nm] = [sz, rest]
                    elif iodocs:
                        last = list(iodocs)[-1]
                        iodocs[last][1] += " " + t.strip()
                    idx += 1
    return h1, " ".join(desc), outputs, example, iodocs


def leading_help(lines, i):
    """Contiguous comment block starting at line i."""
    block = []
    while i < len(lines) and lines[i].lstrip().startswith("%"):
        block.append(lines[i])
        i += 1
    return block, i



def merge_iodocs(args, iodocs):
    """Merge help Inputs/Options docs into arguments-block args and fill
    missing sizes: help-declared size first, then n x 1 for vector-ish
    entries, n x m otherwise (per the API-table convention)."""
    for a in args:
        low = {k.lower(): v for k, v in iodocs.items()}
        doc = (iodocs.get(a["name"]) or iodocs.get("opts." + a["name"])
               or low.get(a["name"].lower())
               or low.get("opts." + a["name"].lower()))
        if doc:
            sz, de = doc
            if de and not a.get("desc"):
                a["desc"] = de
            if not a["size"] and sz:
                szc = sz.replace(",", " x ").replace(":", "n")
                if re.fullmatch(r"[\w. +*-]+( x [\w. +*-]+)*", szc):
                    a["size"] = szc
        if not a["size"]:
            de = (a.get("desc") or "") + " " + a.get("validators", "")
            if a.get("cls") in ("string", "char", "logical") or \
               re.search(r"mustBeScalar|\(1 x 1\)", de):
                a["size"] = "1 x 1"
            elif re.search(r"vector|column|list of|1 x N|N x 1", de, re.I):
                a["size"] = "n x 1"
            else:
                a["size"] = "n x m"
        if not a.get("desc"):
            key = a["name"].lower()
            stem = re.sub(r"(vec|deg|grid|0|1)$", "", key)
            a["desc"] = (CONVENTION_DESCS.get(key)
                         or CONVENTION_DESCS.get(stem, ""))
    return args


# Toolbox-wide fixed meanings (C(n+1,m+1) indexing, geocentric latitudes,
# user-supplied Love numbers). Used only when a help section adds nothing.
CONVENTION_DESCS = {
    # --- v3.10.0 completeness pass: shared toolbox conventions.
    # Defaults in the tables come from the arguments blocks, never from
    # here - these entries carry semantics only.
    "ln": "load Love numbers l_n (horizontal), same layout as kn",
    "quantity": "output functional: 'ewh', 'geoid', 'potential', "
        "'gravity_anomaly', 'gravity_disturbance', "
        "'gravity_gradient_rr' or 'none' (dimensionless passthrough)",
    "latdeg": "geocentric latitudes [deg]",
    "londeg": "longitudes [deg], [0, 360)",
    "latvec": "geocentric latitude vector of the grid [deg]",
    "lonvec": "longitude vector of the grid [deg], [0, 360)",
    "idx": "coefficient index from shLowLevel.shIndex (mind its "
        "MinDegree = 2 default)",
    "grid": "spatial grid, (nLat x nLon) or (nLat x nLon x T)",
    "nmax": "maximum spherical-harmonic degree",
    "nmin": "minimum degree included",
    "mindegree": "lowest degree carried by the index/operation",
    "epoch": "epoch as decimal year",
    "epochs": "epochs as decimal years, (T x 1)",
    "tyears": "time stamps as decimal years, (T x 1)",
    "rho_ave": "mean Earth density [kg/m^3] (default 5517, overridable)",
    "rho_water": "water density [kg/m^3] (default 1000, overridable)",
    "sigmac": "1-sigma uncertainties of C, same indexing",
    "sigmas": "1-sigma uncertainties of S, same indexing",
    "height": "evaluation height above the reference radius [m]",
    "weights": "per-observation or per-cell weights",
    "kaula": "Kaula-rule regularization scale (0 disables)",
    "mask": "logical region mask on the working grid",
    "names": "display names, one per series/solution",
    "plot": "produce the diagnostic figure",
    "ax": "target axes handle ([] creates a new figure)",
    "clim": "color limits [lo hi] ([] = automatic)",
    "coast": "draw coastlines",
    "units": "unit label used for annotation",
    "projection": "map projection name",
    "title": "figure title",
    "colormap": "colormap name or array",
    "dest": "destination folder (created if absent)",
    "baseurl": "server base URL - override for mirrors or testing",
    "timeout": "per-request timeout [s]",
    "proxy": "proxy server URL ('' = direct)",
    "update": "re-download/overwrite files that already exist",
    "quiet": "suppress progress output",
    "release": "product release identifier",
    "catalog": "pre-fetched catalogue table (skips the listing call)",
    "filename": "output/input file path",
    "sidecar": "write the metadata sidecar next to the file",
    "robust": "iteratively reweighted (Huber) estimation",
    "huberk": "Huber tuning constant",
    "maxiter": "iteration cap",
    "tol": "convergence tolerance",
    "oversample": "quadrature-grid refinement factor (boundary "
        "resolution grows linearly, cost quadratically)",
    "bufferkm": "outward buffer of the region boundary [km]",
    "taperkm": "cosine taper width at the region edge [km]",
    "chunksize": "rows processed per block (memory/speed trade-off)",
    "usecache": "reuse the Legendre cache between calls",
    "maxmemgb": "memory ceiling for the operation [GB]",
    "lattype": "'geocentric' (native) or 'geodetic' (converted via "
        "Flattening)",
    "flattening": "flattening used for geodetic-latitude conversion "
        "(default 1/298.257223563)",
    "t0": "reference epoch of the fit [decimal years]",
    "periods": "periodic components to fit [yr]",
    "arcorrect": "inflate sigmas for lag-1 autocorrelation",
    "breaks": "break epochs for piecewise terms [decimal years]",
    "noisecov": "noise covariance from shLowLevel.buildNoiseCov",
    "constraints": "linear constraint spec applied to the estimate",
    "blocks": "order-block structure exploited by the solver",
    "seed": "random seed for reproducibility",
    "keepsamples": "return the raw Monte-Carlo samples",
    "fun": "function handle mapping a coefficient set to the target "
        "quantity",
    "cov": "coefficient covariance (full or per-coefficient)",
    "n": "number of samples/realizations",
    "basin": "basin polygon [lat lon] in degrees or mask",
    "matchtolerance": "epoch matching tolerance [yr]",
    "tslist": "cell array of shSeries to compare",
    "tscell": "cell array of shSeries, one per centre",
    "tsmodel": "model series the scaling is derived from",
    "allowmissing": "tolerate epochs absent from some centres",
    "tolerance": "epoch matching tolerance [yr]",
    "region": "region as polygon [lat lon] deg, mask, or "
        "@(lat, lon) handle (see evalMask)",
    "which": "selection of items to act on",
    "system": "normal-field system, e.g. 'GRS80' or 'WGS84'",
    "a": "semi-major axis of the normal ellipsoid [m]",
    "f": "flattening of the normal ellipsoid",
    "omega": "angular velocity of the normal ellipsoid [rad/s]",
    "gm1": "GM the coefficients are currently scaled to [m^3/s^2]",
    "r1": "R the coefficients are currently scaled to [m]",
    "gm2": "target GM [m^3/s^2]",
    "r2": "target R [m]",
    "columns": "column selection/order of the input table",
    "maxdegree": "truncate the table at this degree",
    "interp": "fill gaps by interpolation across degree",
    "inframe": "reference frame of the input (CM/CE/CF)",
    "outframe": "reference frame of the output (CM/CE/CF)",
    "gapthreshold": "gap length that breaks the plotted line [yr]",
    "trend": "overlay the fitted trend line",
    "label": "series label used in the legend",
    "months": "month selection, 'YYYY-MM' strings",
    "product": "product identifier to download",
    "pattern": "filename glob the folder is scanned with",
    "truncate": "truncate solutions at this degree",
    "xres": "residual coefficient stack after the deterministic fit",
    "fullcov": "assemble the full covariance (memory-heavy)",
    "shrinkage": "shrink off-diagonal covariance toward diagonal",
    "assemble": "return the assembled matrix instead of factors",
    "mode": "operating mode of the routine (see the function help)",
    "signalmode": "signal-covariance construction mode",
    "nitersignal": "signal/noise re-estimation iterations",
    "vcemindegree": "lowest degree entering variance-component "
        "estimation",
    "vcebands": "degree bands for variance-component estimation",
    "floorrel": "relative floor applied to the signal spectrum",
    "mapsmooth": "spatial smoothing radius of the variance map [km]",
    "degvar": "degree-variance model of the signal",
    "noise": "noise degree-variance model or level",
    "signal": "signal degree-variance model",
    "alpha": "filter strength/regularization parameter",
    "op": "linear operator/kernel the map is derived from",
    "t": "target point or epoch of the evaluation",
    "naz": "number of azimuth samples",
    "psimax": "maximum spherical distance evaluated [deg]",
    "npsi": "number of spherical-distance samples",
    "loadregion": "load region as polygon/mask/handle (see evalMask)",
    "ocean": "ocean function as polygon/mask/handle (see evalMask)",
    "loadvalue": "load amplitude assigned to the region [m EWH]",
    "b": "basin-kernel coefficient vector from basinKernel",
    "permonth": "return one factor per month instead of one overall",
    "g1": "first coefficient set",
    "g2": "second coefficient set",
    "in": "input the object is constructed from (see the class help)",
    "producttype": "product type string, e.g. 'gravity_field'",
    "tidesystem": "tide system, e.g. 'zero_tide'",
    "header": "raw header key/value struct",
    "variableterms": "gfct time-variable terms (trnd/acos/asin)",
    "history": "provenance strings carried along",
    "modelname": "model name written to the header",
    "comment": "free-text comment written to the header",
    "cval": "cosine coefficient value",
    "sval": "sine coefficient value",
    "m": "spherical-harmonic order",
    "x": "coefficient stack, (P x T) in index ordering",
    "ts": "shSeries the operation runs on",
    "framerate": "animation frame rate [1/s]",
    "description": "free-text description written to the sidecar",
    "obj": "the object the method is called on",
    "gm": "gravitational constant times Earth mass [m^3/s^2] "
          "(default 3.986004415e14, overridable)",
    "r": "reference radius [m] (default 6378136.3, overridable)",
    "kn": "load Love numbers, (degree, kn) table or column vector - "
          "always user-supplied, no reference frame is assumed",
    "hn": "load Love numbers h_n, same layout as kn",
    "c": "cosine Stokes coefficients, C(n+1, m+1) indexing",
    "s": "sine Stokes coefficients, S(n+1, m+1) indexing",
    "cs": "stacked cosine coefficients, C(n+1, m+1, T)",
    "ss": "stacked sine coefficients, S(n+1, m+1, T)",
    "lat": "geocentric latitudes [deg]",
    "latdeg": "geocentric latitudes [deg]",
    "lon": "longitudes [deg], lon in [0, 360)",
    "londeg": "longitudes [deg], lon in [0, 360)",
    "nmax": "maximum spherical harmonic degree",
    "nmin": "minimum spherical harmonic degree",
    "n": "spherical harmonic degree",
    "m": "spherical harmonic order",
    "idx": "index struct from shLowLevel.shIndex (fields n, m, cs, P, "
           "Lmax, minDegree, pos)",
    "epoch": "epoch [decimal years]",
    "epochs": "epochs [decimal years]",
    "quiet": "suppress progress output",
    "quantity": "functional of the field, see shLowLevel.kernelFactors",
    "folder": "folder path",
    "grid": "field values on the (lat, lon) grid",
    "w": "filter container from shLowLevel.readDDK / designFilter, or "
         "degree weights",
    "ts": "shSeries object",
    "tn": "correction-table struct as returned by the matching reader",
}


def parse_function_file(path):
    raw = open(path).read().split("\n")
    lines = raw
    m = None
    for i, L in enumerate(lines):
        m = re.match(r"\s*function\s+(?:(\[[^\]]*\]|\w+)\s*=\s*)?(\w+)\s*"
                     r"\(([^)]*)\)", L)
        if m:
            break
    if not m:
        return None
    outs = (m.group(1) or "").strip("[]")
    outs = [o.strip() for o in outs.split(",") if o.strip()]
    name = m.group(2)
    ins = [a.strip() for a in m.group(3).split(",") if a.strip()]
    hb, j = leading_help(lines, i + 1)
    h1, desc, outdocs, example, iodocs = parse_help(hb)
    args = []
    jl = _join_continuations(lines[j:j + 200])
    for k, L in enumerate(jl):
        if L.strip() == "arguments":
            args, _ = parse_arguments_block(jl, k)
            break
        if re.match(r"^\s*function\s", L):
            break                    # never cross into local functions
    args = merge_iodocs(args, iodocs)
    return dict(kind="function", name=name, inputs=ins, outputs=outs,
                h1=h1, desc=desc, outdocs=outdocs, example=example,
                args=args)


def parse_classdef(path):
    src = open(path).read()
    lines = src.split("\n")
    cname = re.search(r"^classdef\s+(\w+)", src, re.M).group(1)
    ci = [i for i, L in enumerate(lines) if L.startswith("classdef")][0]
    hb, _ = leading_help(lines, ci + 1)
    ch1, cdesc, _, _, _ = parse_help(hb)
    props = []
    for pm in re.finditer(r"^properties\s*(\([^)]*\))?\s*$(.*?)^end\s*$",
                          src, re.M | re.S):
        attrs = (pm.group(1) or "").lower()
        # SetAccess=private is fine (immutable style); exclude only
        # blocks the USER cannot read: GetAccess private/hidden
        if re.search(r"getaccess\s*=\s*private", attrs) or "hidden" in attrs:
            continue
        for L in pm.group(2).split("\n"):
            code_part, _, cmt = L.partition("%")
            mo = re.match(
                r"^\s{4}(\w+)\s*(\([^)]*\))?\s*([a-zA-Z]\w*)?"
                r"\s*(?:\{[^}]*\})?\s*(?:=\s*(.+))?\s*$", code_part)
            if mo and mo.group(1) not in ("end", "properties"):
                props.append([mo.group(1),
                              (mo.group(2) or "").strip("()").replace(",", " x "),
                              (mo.group(3) or ""),
                              (mo.group(4) or "").strip(),
                              cmt.strip()])
    # methods blocks are at column 0 in this codebase; a block runs to
    # the next column-0 methods/properties/events/end keyword. Functions
    # are declared at 4-space indent; deeper indents are nested locals.
    starts = []
    for i, L in enumerate(lines):
        m = re.match(r"^methods\s*(\(([^)]*)\))?\s*$", L)
        if m:
            starts.append((i, (m.group(2) or "").lower()))
    delim = [i for i, L in enumerate(lines)
             if re.match(r"^(methods|properties|events|end)\b", L)]
    methods = []
    for s, attrs in starts:
        if "private" in attrs or "hidden" in attrs:
            continue
        static = "static" in attrs
        epos = min([d for d in delim if d > s], default=len(lines))
        i = s + 1
        while i < epos:
            fm = re.match(r"^\s{4}function\s+(?:(\[[^\]]*\]|\w+)\s*=\s*)?"
                          r"(\w+)\s*\(([^)]*)\)", lines[i])
            if fm:
                outs = (fm.group(1) or "").strip("[]")
                outs = [o.strip() for o in outs.split(",") if o.strip()]
                mname = fm.group(2)
                ins = [a.strip() for a in fm.group(3).split(",")
                       if a.strip()]
                hb, j = leading_help(lines, i + 1)
                h1, desc, outdocs, example, iodocs = parse_help(hb)
                jl = _join_continuations(lines[j:j + 120])
                margs = []
                for k, L2 in enumerate(jl):
                    if L2.strip() == "arguments":
                        margs, _ = parse_arguments_block(jl, k)
                        break
                    if re.match(r"^\s*function\s", L2):
                        break
                margs = merge_iodocs(margs, iodocs)
                methods.append(dict(
                    kind="method", cls=cname, name=mname, static=static,
                    inputs=ins, outputs=outs, h1=h1, desc=desc,
                    outdocs=outdocs, example=example, args=margs))
            i += 1
    return dict(name=cname, h1=ch1, desc=cdesc, props=props,
                methods=methods)


def main():
    entries = {"classes": [], "shLowLevel": [], "root": []}
    for f in ["shCoefficients.m", "shSeries.m", "shClimatology.m"]:
        entries["classes"].append(parse_classdef(os.path.join(ROOT, f)))
    for f in sorted(glob.glob(os.path.join(ROOT, "+shLowLevel", "*.m"))):
        e = parse_function_file(f)
        if e:
            entries["shLowLevel"].append(e)
    for f in ["setup_shAnalysis.m", "runAllTests.m", "demo_shAnalysis.m"]:
        e = parse_function_file(os.path.join(ROOT, f))
        if e:
            entries["root"].append(e)
    # coverage report
    miss_ex, miss_desc = [], []
    for e in entries["shLowLevel"] + entries["root"]:
        if not e["example"]:
            miss_ex.append(e["name"])
        if not (e["h1"] or e["desc"]):
            miss_desc.append(e["name"])
    nm = 0
    for c in entries["classes"]:
        for m in c["methods"]:
            nm += 1
            if not m["example"]:
                miss_ex.append(c["name"] + "." + m["name"])
            if not (m["h1"] or m["desc"]):
                miss_desc.append(c["name"] + "." + m["name"])
    print(f"classes: {len(entries['classes'])} with {nm} public methods; "
          f"shLowLevel: {len(entries['shLowLevel'])}; root: {len(entries['root'])}")
    print(f"missing example: {len(miss_ex)}")
    print("  " + ", ".join(miss_ex))
    print(f"missing description: {len(miss_desc)}")
    if miss_desc:
        print("  " + ", ".join(miss_desc))
    # v3.10.0: convention fallback as a FINAL pass over every entity,
    # whatever extraction path built it - help text always wins, the
    # lexicon only fills what stayed blank.
    def _fallback(a):
        if isinstance(a, dict) and not (a.get("desc") or "").strip():
            key = str(a.get("name", "")).lower()
            stem = re.sub(r"(vec|deg|grid|0|1)$", "", key)
            a["desc"] = (CONVENTION_DESCS.get(key)
                         or CONVENTION_DESCS.get(stem, ""))
    def _walk_entity(e):
        for k in ("args", "options", "inputs", "outputs", "props",
                  "properties"):
            for a in e.get(k, []) or []:
                _fallback(a)
        for meth in e.get("methods", []) or []:
            _walk_entity(meth)
    for _grp in entries.values() if isinstance(entries, dict) else [entries]:
        for _e in _grp:
            if isinstance(_e, dict):
                _walk_entity(_e)
    with open("/home/claude/api_data.py", "w") as f:
        f.write("# generated by api_extract.py - do not edit\n")
        f.write("API = " + repr(entries) + "\n")
    print("wrote /home/claude/api_data.py")


if __name__ == "__main__":
    main()
