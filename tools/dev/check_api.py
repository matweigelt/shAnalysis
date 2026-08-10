"""Cross-check API usage in demo + doc snippets against the code, plus
repo-wide inputParser consistency. Recreated for the v2.4.1 session per the
handoff spec; guards against the paid-for bug classes (invented methods
tvANS/meanField, wrong positional arity, orphaned p.Results.X).

  C1  obj.method( in demo + snippets must exist on the mapped class:
        ts*   -> shSeries      g*    -> shCoefficients
        clim* -> shClimatology
      (methods AND properties count; get./set. stripped)
  C2  shLowLevel.fun( positional-arity check against +shLowLevel signatures
      (continuations folded first; args before the first Name=value pair;
       varargin signatures skipped; 'opts' never positional)
  C3  every p.Results.X in the repo needs a matching addParameter/
      addOptional/addRequired for X in the same file

Snippet sources: demo_shAnalysis.m and <pre> blocks in html/*.html
(the guide-builder snippets re-attach when make_guide.py is rebuilt).
Claude (Fable 5), 2026-08-07.
"""
import glob
import html as htmllib
import os
import re
import sys

ROOT = "/home/claude/shx_build/shAnalysis"

CLASS_FILES = {"shSeries": "shSeries.m",
               "shCoefficients": "shCoefficients.m",
               "shClimatology": "shClimatology.m"}

PREFIX_MAP = [("clim", "shClimatology"),   # longest prefix first
              ("ts", "shSeries"),
              ("g", "shCoefficients")]

# names that are base-MATLAB struct/field or builtin patterns, not methods
IGNORE_METHODS = {"Results", "sigma", "C", "S"}


def fold(txt):
    return re.sub(r"\.\.\.[^\n]*\n", " ", txt)


def strip_comments(txt):
    return re.sub(r"%[^\n]*", "", txt)


def class_api(path):
    src = open(path, encoding="utf-8").read()
    api = set()
    for m in re.finditer(r"^\s*function\s+(?:\[[^\]]*\]\s*=\s*|\w+\s*=\s*)?"
                         r"([\w\.]+)\s*[\(\n,]", fold(src), re.M):
        name = m.group(1)
        name = re.sub(r"^(get|set)\.", "", name)
        api.add(name)
    in_props = False
    for line in src.split("\n"):
        s = line.strip()
        if re.match(r"properties\b", s):
            in_props = True
            continue
        if in_props:
            if s == "end":
                in_props = False
                continue
            mm = re.match(r"([A-Za-z_]\w*)\b", s)
            if mm:
                api.add(mm.group(1))
    return api


def pkg_signatures():
    sig = {}
    for p in sorted(glob.glob(os.path.join(ROOT, "+shLowLevel", "*.m"))):
        head = fold(open(p, encoding="utf-8").read()[:2000])
        mm = re.search(r"function\s+(?:\[[^\]]*\]\s*=\s*|\w+\s*=\s*)?"
                       r"(\w+)\s*\(([^)]*)\)", head)
        if not mm:
            continue
        name = mm.group(1)
        params = [q.strip() for q in mm.group(2).split(",") if q.strip()]
        if "varargin" in params:
            continue
        sig[name] = sum(1 for q in params if q != "opts")
    return sig


def balanced_call(text, start):
    depth = 0
    for j in range(start, len(text)):
        if text[j] == "(":
            depth += 1
        elif text[j] == ")":
            depth -= 1
            if depth == 0:
                return text[start + 1:j]
    return None


def toplevel_args(s):
    args, depth, cur = [], 0, ""
    for ch in s:
        if ch in "([{":
            depth += 1
        if ch in ")]}":
            depth -= 1
        if ch == "," and depth == 0:
            args.append(cur)
            cur = ""
        else:
            cur += ch
    if cur.strip():
        args.append(cur)
    return args


def scan_methods(name, text, apis):
    bad = []
    txt = strip_comments(fold(text))
    for mm in re.finditer(r"\b([a-zA-Z]\w*)\.([a-zA-Z]\w*)\s*\(", txt):
        var, meth = mm.group(1), mm.group(2)
        if var == "shLowLevel" or meth in IGNORE_METHODS:
            continue
        cls = None
        for pre, c in PREFIX_MAP:
            if var.lower().startswith(pre):
                cls = c
                break
        if cls is None:
            continue
        if meth not in apis[cls]:
            ln = text[:mm.start()].count("\n") + 1
            bad.append((name, ln, f"{var}.{meth}( -- no method/property "
                                  f"'{meth}' on {cls}"))
    return bad


def scan_pkg(name, text, sig):
    bad = []
    txt = strip_comments(fold(text))
    for mm in re.finditer(r"\bshx\.(\w+)\s*\(", txt):
        fn = mm.group(1)
        if fn not in sig:
            if os.path.exists(os.path.join(ROOT, "+shLowLevel", fn + ".m")):
                continue                       # varargin or unparsed
            ln = text[:mm.start()].count("\n") + 1
            bad.append((name, ln, f"shLowLevel.{fn} -- no such package function"))
            continue
        inner = balanced_call(txt, mm.end() - 1)
        if inner is None:
            continue
        pos = 0
        for a in toplevel_args(inner):
            if re.match(r"\s*\w+\s*=[^=]", a):
                break
            pos += 1
        if pos > sig[fn]:
            ln = text[:mm.start()].count("\n") + 1
            bad.append((name, ln, f"shLowLevel.{fn}: {pos} positional args, "
                                  f"signature takes {sig[fn]}"))
    return bad


def scan_inputparser(path):
    bad = []
    src = fold(open(path, encoding="utf-8").read())
    declared = set(re.findall(
        r"add(?:Parameter|Optional|Required)\s*\(\s*\w+\s*,\s*['\"](\w+)['\"]",
        src))
    if not declared and "p.Results." not in src:
        return bad
    for mm in re.finditer(r"\bp\.Results\.(\w+)", src):
        if mm.group(1) not in declared:
            ln = src[:mm.start()].count("\n") + 1
            bad.append((path, ln, f"p.Results.{mm.group(1)} has no "
                                  f"addParameter/Optional/Required"))
    return bad


def html_snippets():
    out = []
    for p in sorted(glob.glob(os.path.join(ROOT, "html", "*.html"))):
        src = open(p, encoding="utf-8", errors="replace").read()
        for i, b in enumerate(re.findall(r"<pre[^>]*>(.*?)</pre>",
                                         src, re.S)):
            out.append((f"{os.path.basename(p)}#pre{i}",
                        htmllib.unescape(re.sub(r"<[^>]+>", "", b))))
    return out


def main():
    apis = {c: class_api(os.path.join(ROOT, f))
            for c, f in CLASS_FILES.items()}
    sig = pkg_signatures()
    bad = []

    demo = open(os.path.join(ROOT, "demo_shAnalysis.m"),
                encoding="utf-8").read()
    bad += scan_methods("demo_shAnalysis.m", demo, apis)
    bad += scan_pkg("demo_shAnalysis.m", demo, sig)

    for name, snip in html_snippets():
        bad += scan_methods(name, snip, apis)
        bad += scan_pkg(name, snip, sig)

    for p in sorted(glob.glob(os.path.join(ROOT, "**", "*.m"),
                              recursive=True)):
        bad += scan_inputparser(p)

    for name, ln, msg in bad:
        print(f"{name}:{ln}: {msg}")
    n_api = sum(len(a) for a in apis.values())
    print(f"check_api: {n_api} class API names, {len(sig)} shLowLevel signatures, "
          f"{len(bad)} finding(s)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
