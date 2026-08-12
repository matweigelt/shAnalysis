#!/usr/bin/env python3
"""doc_sync_audit: keep the THREE documentation sources in sync.

shAnalysis documents itself in three places, and before v3.1.2 they had
drifted apart without a single gate noticing - all five existing gates
were green while the API reference listed 5 of 12 fetchICGEM options,
eleven help pages were unreachable from the Help browser, the tagged
release reported the previous version, and the workflow guide advertised
a call that threw. Each of those is checked here.

  1. the in-file help text        (help shLowLevel.foo)
  2. html/                        (doc shAnalysis, incl. the generated
                                   apiReference.html)
  3. docs/shAnalysis_workflow_guide.pdf, authored in tools/dev/make_guide.py

Checks, in order of how badly they bite:

  SNIPPET   every call in an html <pre> block or a guide code block must
            match the parsed contract: the function exists, positional
            count in range, option names declared, and option VALUES of
            the declared type (a string literal into a logical option, or
            a value outside a mustBeMember set, is a snippet that throws).
  APIREF    every public entity and every one of its name-value options
            appears in html/apiReference.html - i.e. it was regenerated.
  TOC       html/helptoc.xml reaches every page and points at no ghosts.
  VERSION   Contents.m (the single source of truth) agrees with
            CITATION.cff, the top CHANGELOG section and the API
            reference; no "Unreleased" section survives a release.
  HELP      no option is "documented" as a pointer to the code, and every
            documented default matches the arguments block.
  COVERAGE  every public entity is mentioned in at least one narrative
            source, not only in the generated reference.

Stdlib only; run from the repository root. Nonzero exit on any finding.

Developed by Matthias Weigelt with the help of Claude (Opus 5),
2026-08-11 (v3.1.2).
"""
import html as _html
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = os.path.join(ROOT, "+shLowLevel")
HTML = os.path.join(ROOT, "html")
GUIDE_SRC = os.path.join(ROOT, "tools", "dev", "make_guide.py")
ROOT_FUNCS = ["setup_shAnalysis.m", "runAllTests.m", "demo_shAnalysis.m"]
CLASSES = ["shCoefficients.m", "shSeries.m", "shClimatology.m"]

# "Name=Value" is MATLAB's own signature placeholder, not a real option
PLACEHOLDER_OPTS = {"Name=Value", "NameValue"}
# prose in narrative pages writes calls with a literal ellipsis
ELLIPSIS = "..."


# --------------------------------------------------------------- parsing
def join_cont(lines):
    out, buf = [], ""
    for L in lines:
        s = L.rstrip("\n")
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


def parse_arguments(lines, start):
    """Parse one arguments...end block; returns a list of option dicts."""
    args, i = [], start + 1
    while i < len(lines):
        s = lines[i].strip()
        if s == "end":
            break
        if s and not s.startswith("%"):
            s = s.split("%")[0].strip()
            m = re.match(r"([\w.]+)\s*"          # name (may be opts.X)
                         r"(\([^)]*\))?\s*"      # (size)
                         r"([a-zA-Z]\w*)?\s*"    # class
                         r"(\{.*\})?\s*"         # {validators}
                         r"(?:=\s*(.+))?$", s)
            if m:
                name, size, cls, vals, dflt = m.groups()
                nv = name.startswith("opts.")
                members = None
                if vals:
                    mm = re.search(r"mustBeMember\s*\([^,]+,\s*(.+?)\)\s*\}?\s*$",
                                   vals)
                    if mm:
                        # split the list - a '""' member defeats a pair regex
                        inner = mm.group(1).strip().strip("[]{}")
                        members = {t.strip().strip("'\"")
                                   for t in split_args(inner)}
                        members = members or None
                args.append(dict(name=name.split(".", 1)[1] if nv else name,
                                 nv=nv, cls=cls or "",
                                 members=members,
                                 default=(dflt or "").strip()))
        i += 1
    return args


def parse_one(lines, i0):
    """Parse the function starting at/after i0."""
    fre = re.compile(r"^\s*function\s+(?:(\[[^\]]*\]|\w+)\s*=\s*)?"
                     r"(\w+)\s*\(([^)]*)\)")
    i = i0
    while i < len(lines) and not fre.match(lines[i]):
        i += 1
    if i >= len(lines):
        return None
    m = fre.match(lines[i])
    outs = re.findall(r"\w+", m.group(1)) if m.group(1) else []
    ins = [a.strip() for a in m.group(3).split(",") if a.strip()]
    j = i + 1
    help_lines = []
    while j < len(lines) and lines[j].lstrip().startswith("%"):
        help_lines.append(lines[j].lstrip().lstrip("%").rstrip("\n"))
        j += 1
    args = []
    jl = join_cont(lines[j:j + 250])
    for k, L in enumerate(jl):
        if re.match(r"^\s*function\s", L):
            break
        if L.strip() == "arguments":
            args = parse_arguments(jl, k)
            break
    return dict(name=m.group(2), outs=outs, ins=ins,
                help="\n".join(help_lines), args=args)


def entities():
    """All public functions and methods, keyed by qualified name."""
    out = {}
    for f in sorted(os.listdir(PKG)):
        if f.endswith(".m"):
            e = parse_one(open(os.path.join(PKG, f), encoding="utf-8",
                               errors="replace").readlines(), 0)
            if e:
                out["shLowLevel." + e["name"]] = e
    for f in ROOT_FUNCS:
        p = os.path.join(ROOT, f)
        if os.path.isfile(p):
            e = parse_one(open(p, encoding="utf-8",
                               errors="replace").readlines(), 0)
            if e:
                out[e["name"]] = e
    for cls in CLASSES:
        p = os.path.join(ROOT, cls)
        if not os.path.isfile(p):
            continue
        lines = open(p, encoding="utf-8", errors="replace").readlines()
        i, seen = 0, set()
        while i < len(lines):
            mb = re.match(r"^(\s*)methods\b(.*)$", lines[i])
            if not mb:
                i += 1
                continue
            ind, attrs = mb.group(1), mb.group(2)
            public = not (re.search(r"Access\s*=\s*'?(private|protected)",
                                    attrs) or "Hidden" in attrs)
            j = i + 1
            endre = re.compile("^" + ind + r"end\b")
            while j < len(lines) and not endre.match(lines[j]):
                if public and re.match(r"^\s*function\s", lines[j]):
                    e = parse_one(lines, j)
                    if e and e["name"] not in seen:
                        seen.add(e["name"])
                        out[cls[:-2] + "." + e["name"]] = e
                j += 1
            i = j + 1
    return out


# -------------------------------------------------------- snippet checks
def split_args(s):
    out, d, cur = [], 0, ""
    for ch in s:
        if ch in "([{":
            d += 1
        if ch in ")]}":
            d -= 1
        if ch == "," and d == 0:
            out.append(cur)
            cur = ""
        else:
            cur += ch
    if cur.strip():
        out.append(cur)
    return [a.strip() for a in out if a.strip()]


def sig_of(e):
    ins = e["ins"]
    varargin = "varargin" in ins
    args = e.get("args", [])
    pos = [a for a in args if not a["nv"] and a["name"] not in ("obj",)]
    nv = {a["name"] for a in args if a["nv"]}
    byname = {a["name"]: a for a in args if a["nv"]}
    if args:
        minp = len([a for a in pos if not a["default"]])
        maxp = 999 if varargin else len(pos)
        return minp, maxp, (None if varargin else nv), byname
    ins2 = [i for i in ins if i not in ("obj", "opts", "varargin")]
    return (1 if ins2 else 0), (999 if varargin else len(ins2)), None, {}


def check_value(fn, opt, spec, where, probs):
    """The class that a name-only check misses: a literal of the wrong TYPE."""
    val = spec.strip()
    a = opt
    # narrative pages list the alternatives: Method="auto"|"rings"|"ls"
    if "|" in val and a["members"]:
        for alt in re.findall(r"['\"]([^'\"]*)['\"]", val):
            if alt not in a["members"]:
                probs.append("%s: %s(..., %s=%s) - not one of {%s}"
                             % (where, fn, a["name"], alt,
                                ", ".join(sorted(a["members"]))))
        return
    lit = re.match(r"^(['\"])(.*)\1$", val)
    if lit:
        text = lit.group(2)
        if a["cls"] == "logical":
            probs.append("%s: %s(..., %s=%s) - %s is a LOGICAL flag, not a "
                         "path/name" % (where, fn, a["name"], val, a["name"]))
        elif a["cls"] == "double":
            probs.append("%s: %s(..., %s=%s) - %s is numeric"
                         % (where, fn, a["name"], val, a["name"]))
        elif a["members"] and text not in a["members"]:
            probs.append("%s: %s(..., %s=%s) - not one of {%s}"
                         % (where, fn, a["name"], val,
                            ", ".join(sorted(a["members"]))))
    elif a["cls"] == "logical" and re.match(r"^[\d.]+e?[-+]?\d*$", val):
        if val not in ("0", "1"):
            probs.append("%s: %s(..., %s=%s) - %s is a logical flag"
                         % (where, fn, a["name"], val, a["name"]))


def audit_snippet(where, src, sig, meth, probs):
    src = clean_matlab(src)
    for mo in re.finditer(r"([\w.]+)\s*\(", src):
        fn = mo.group(1)
        i, d, buf = mo.end(), 1, ""
        while i < len(src) and d:
            ch = src[i]
            if ch == "(":
                d += 1
            if ch == ")":
                d -= 1
            if d:
                buf += ch
            i += 1
        cands = []
        if fn in sig:
            cands = [(fn, sig[fn])]
        else:
            mm = re.match(r"^(\w+)\.(\w+)$", fn)
            if mm and mm.group(1) != "shLowLevel" and mm.group(2) in meth:
                cands = [(fn, s) for s in meth[mm.group(2)]]
        if not cands:
            continue
        trials = []
        for name, (minp, maxp, nvs, byname) in cands:
            local = []
            raw = split_args(buf)
            elided = ELLIPSIS in raw
            args = [a for a in raw if a != ELLIPSIS]
            npos, used = 0, []
            for a in args:
                m = re.match(r"^([A-Za-z]\w*)\s*=\s*([^=].*)$", a, re.S)
                if m and not a.startswith("@"):
                    used.append((m.group(1), m.group(2)))
                elif not used:
                    npos += 1
            if not elided and not (minp <= npos <= maxp):
                local.append("%s: %s(...) uses %d positional args, contract "
                             "allows [%d,%d]" % (where, name, npos, minp, maxp))
            for nm, val in used:
                if nm + "=" + val.strip() in PLACEHOLDER_OPTS:
                    continue
                if nvs is not None and nm not in nvs:
                    local.append("%s: %s(..., %s=...) - option not declared "
                                 "in the arguments block" % (where, name, nm))
                elif nm in byname:
                    check_value(name, byname[nm], val, where, local)
            trials.append(local)
        if all(t for t in trials):        # every candidate rejects it
            probs.extend(trials[0])


def clean_matlab(src):
    """Strip comments and join continuations, honouring string literals.

    Without this a trailing comment on a continued line ("..., ...  % note")
    leaves the NEXT argument glued behind an ellipsis, where the name=value
    regex no longer matches it - so the argument is silently dropped and
    never checked. That defeated the type check in testing.
    """
    out = []
    for line in src.split("\n"):
        q = None          # None | "'" | '"'
        cut = len(line)
        i = 0
        while i < len(line):
            ch = line[i]
            if q:
                if ch == q:
                    if i + 1 < len(line) and line[i + 1] == q:
                        i += 1          # doubled quote inside the literal
                    else:
                        q = None
            elif ch in "'\"":
                # a quote after an operand is a transpose, not a literal
                if ch == "'" and i and (line[i - 1].isalnum()
                                        or line[i - 1] in ")]}._"):
                    pass
                else:
                    q = ch
            elif ch == "%":
                cut = i
                break
            i += 1
        out.append(line[:cut].rstrip())
    joined, buf = [], ""
    for line in out:
        if line.endswith("..."):
            buf += line[:-3].rstrip() + " "
            continue
        joined.append(buf + line.strip() if buf else line)
        buf = ""
    if buf:
        joined.append(buf)
    return "\n".join(joined)


# defaults that are SENTINELS: the help may document the EFFECT instead
# (T0 (mean(tYears)) is better documentation than T0 (NaN))
SENTINELS = {"NaN", "[]", "", "table()", "struct()", "struct([])",
             "string.empty(0,1)", "double.empty(1,0)", "strings(1,0)",
             "true(0,1)"}


def styled_default(help_text, name):
    """The house style is "Name (default)  description".

    Returns the parenthesised default, or None when the help does not use
    the style for this option. Parentheses are matched by depth - a plain
    regex captures "6), Tol (1e-3" out of a line listing two options.
    """
    for m in re.finditer(r"(?:^|\s)" + re.escape(name) + r"\s*\(", help_text):
        i, d = m.end(), 1
        while i < len(help_text) and d:
            if help_text[i] == "(":
                d += 1
            elif help_text[i] == ")":
                d -= 1
            i += 1
        if d:
            continue
        inner = help_text[m.end():i - 1]
        rest = help_text[i:]
        if rest.startswith("  ") or rest.startswith("\n") or rest == "":
            return inner
    return None


def html_text(path):
    s = open(path, encoding="utf-8", errors="replace").read()
    s = re.sub(r"<style.*?</style>", " ", s, flags=re.S)
    return _html.unescape(re.sub(r"<[^>]+>", " ", s))


def html_snippets(path):
    s = open(path, encoding="utf-8", errors="replace").read()
    for blk in re.findall(r"<pre>(.*?)</pre>", s, re.S):
        yield _html.unescape(re.sub(r"<[^>]+>", "", blk))


def guide_snippets():
    if not os.path.isfile(GUIDE_SRC):
        return
    g = open(GUIDE_SRC, encoding="utf-8").read()
    for blk in re.findall(r'code\(\s*"""(.*?)"""', g, re.S):
        yield blk


# ------------------------------------------------------------- the audit
def main():
    ent = entities()
    probs = []
    sig = {k: sig_of(e) for k, e in ent.items()}
    meth = {}
    for k, e in ent.items():
        if "." in k and not k.startswith("shLowLevel."):
            meth.setdefault(k.split(".", 1)[1], []).append(sig_of(e))

    # --- SNIPPET
    pages = sorted(f for f in os.listdir(HTML) if f.endswith(".html"))
    for f in pages:
        if f == "apiReference.html":
            continue
        for blk in html_snippets(os.path.join(HTML, f)):
            audit_snippet("html/" + f, blk, sig, meth, probs)
    for blk in guide_snippets():
        audit_snippet("guide", blk, sig, meth, probs)

    # --- APIREF: the generated reference must know every option
    apath = os.path.join(HTML, "apiReference.html")
    api = open(apath, encoding="utf-8", errors="replace").read()
    for k, e in sorted(ent.items()):
        if 'id="%s"' % k not in api:
            probs.append("apiReference.html: no entry for %s - regenerate "
                         "with tools/dev/make_apiref.py" % k)
            continue
        blk = api.split('id="%s"' % k, 1)[1].split("</div>", 1)[0]
        for a in e["args"]:
            if a["nv"] and "<code>%s</code>" % a["name"] not in blk:
                probs.append("apiReference.html: %s is missing option '%s' - "
                             "the page is stale" % (k, a["name"]))

    # --- TOC
    toc = open(os.path.join(HTML, "helptoc.xml"), encoding="utf-8").read()
    targets = set(re.findall(r'target="([^"]+)"', toc))
    for f in set(pages) - targets:
        probs.append("helptoc.xml: %s is unreachable from the Help browser" % f)
    for f in targets - set(pages):
        probs.append("helptoc.xml: target %s does not exist" % f)

    # --- VERSION
    cm = open(os.path.join(ROOT, "Contents.m"), encoding="utf-8").read()
    mv = re.search(r"^%\s*Version\s+(\S+)", cm, re.M)
    if not mv:
        probs.append("Contents.m: no Version line")
    else:
        ver = mv.group(1)
        cff = open(os.path.join(ROOT, "CITATION.cff"), encoding="utf-8").read()
        mc = re.search(r'version:\s*"([^"]+)"', cff)
        if not mc or mc.group(1) != ver:
            probs.append("CITATION.cff: version %s != Contents.m %s"
                         % (mc.group(1) if mc else "?", ver))
        chg = open(os.path.join(ROOT, "CHANGELOG.md"), encoding="utf-8").read()
        mh = re.search(r"^##\s*\[([^\]]+)\]", chg, re.M)
        if not mh or mh.group(1) != ver:
            probs.append("CHANGELOG.md: top section [%s] != Contents.m %s"
                         % (mh.group(1) if mh else "?", ver))
        for openhdr in re.findall(r"^##\s*\[[^\]]+\]\s*-\s*Unreleased.*$",
                                  chg, re.M):
            probs.append("CHANGELOG.md: %s - a tagged release leaves no "
                         "Unreleased section" % openhdr.strip())
        if "API reference (v%s)" % ver not in api:
            probs.append("apiReference.html: title is not v%s - regenerate it"
                         % ver)
        # The overview page is where a user lands from `doc shAnalysis`.
        # It carried "v3.0.0" for eight releases because nothing checked
        # it: the gate only looked at pages that already had a version.
        ov = os.path.join(ROOT, "html", "shAnalysis.html")
        if os.path.isfile(ov):
            ot = open(ov, encoding="utf-8").read()
            mo = re.search(r'class="ver">\s*Version\s+(\S+?)\s', ot)
            if not mo:
                probs.append("shAnalysis.html: no version stamp - the "
                             "overview page must state the version so it "
                             "cannot silently go stale")
            elif mo.group(1) != ver:
                probs.append("shAnalysis.html: version %s != Contents.m %s"
                             % (mo.group(1), ver))
        # the first Contents.m line is what ver() shows as the product NAME
        first = cm.split("\n", 1)[0].lstrip("% ").strip()
        if len(first) >= 60 or re.search(r"\d+\.\d+", first):
            probs.append("Contents.m: line 1 is shown by ver() as the product "
                         "name - keep it short and version-free (%r)" % first)

    # --- HELP
    for k, e in sorted(ent.items()):
        for a in e["args"]:
            if not a["nv"]:
                continue
            segs = re.findall(r"^.*\b" + re.escape(a["name"]) + r"\b.*$",
                              e["help"], re.M)
            if not segs:
                continue
            if any("see arguments block" in s.lower() for s in segs):
                probs.append("%s: option '%s' is documented as a pointer to "
                             "the code" % (k, a["name"]))
                continue
            d = re.sub(r"[\s'\"]+", "", a["default"])
            if not d:
                continue
            # house style is "Name (default)  description" - when the help
            # uses it, compare the parenthesised value EXACTLY. A lenient
            # substring test cannot catch a drifted one-character default
            # like (2) -> (7), because "2" occurs in almost any line.
            styled = styled_default(e["help"], a["name"])
            if styled is not None:
                if (re.sub(r"[\s'\"]+", "", styled) != d
                        and d not in SENTINELS):
                    probs.append("%s: help says %s (%s) but the arguments "
                                 "block says %s"
                                 % (k, a["name"], styled, a["default"]))
            elif not any(d[:20] in re.sub(r"[\s'\"]+", "", s)
                         for s in segs):
                probs.append("%s: option '%s' documented without its "
                             "arguments-block default (%s)"
                             % (k, a["name"], a["default"]))

    # --- COVERAGE: narrative sources, not just the generated reference
    blob = "\n".join(html_text(os.path.join(HTML, f))
                     for f in pages if f != "apiReference.html")
    if os.path.isfile(GUIDE_SRC):
        blob += "\n" + open(GUIDE_SRC, encoding="utf-8").read()
    for k in sorted(ent):
        short = k.split(".")[-1]
        if not re.search(r"\b" + re.escape(short) + r"\b", blob):
            probs.append("%s is documented ONLY in the generated API "
                         "reference - no narrative page or guide section" % k)

    # structural (audit F-6): content after </html> renders by browser
    # error recovery and passes every text-presence check above
    for fn in sorted(os.listdir(HTML)):
        if not fn.endswith(".html"):
            continue
        d = open(os.path.join(HTML, fn), encoding="utf-8").read()
        i = d.rfind("</html>")
        if i < 0:
            probs.append("%s: no closing </html>" % fn)
        elif d[i + 7:].strip():
            probs.append("%s: %d bytes after </html>" % (fn, len(d[i + 7:].strip())))
        if d.count("</html>") > 1:
            probs.append("%s: multiple </html>" % fn)

    if probs:
        print("doc_sync_audit: %d finding(s)" % len(probs))
        for p in probs:
            print("  -", p)
        sys.exit(1)
    print("doc_sync_audit: 0 findings - help, html and guide agree "
          "(%d public entities, %d help pages)" % (len(ent), len(pages)))


if __name__ == "__main__":
    main()
