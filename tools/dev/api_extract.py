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
    return h1, " ".join(desc), outputs, example


def leading_help(lines, i):
    """Contiguous comment block starting at line i."""
    block = []
    while i < len(lines) and lines[i].lstrip().startswith("%"):
        block.append(lines[i])
        i += 1
    return block, i


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
    h1, desc, outdocs, example = parse_help(hb)
    args = []
    jl = _join_continuations(lines[j:j + 200])
    for k, L in enumerate(jl):
        if L.strip() == "arguments":
            args, _ = parse_arguments_block(jl, k)
            break
        if re.match(r"^\s*function\s", L):
            break                    # never cross into local functions
    return dict(kind="function", name=name, inputs=ins, outputs=outs,
                h1=h1, desc=desc, outdocs=outdocs, example=example,
                args=args)


def parse_classdef(path):
    src = open(path).read()
    lines = src.split("\n")
    cname = re.search(r"^classdef\s+(\w+)", src, re.M).group(1)
    ci = [i for i, L in enumerate(lines) if L.startswith("classdef")][0]
    hb, _ = leading_help(lines, ci + 1)
    ch1, cdesc, _, _ = parse_help(hb)
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
                h1, desc, outdocs, example = parse_help(hb)
                jl = _join_continuations(lines[j:j + 120])
                margs = []
                for k, L2 in enumerate(jl):
                    if L2.strip() == "arguments":
                        margs, _ = parse_arguments_block(jl, k)
                        break
                    if re.match(r"^\s*function\s", L2):
                        break
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
    with open("/home/claude/api_data.py", "w") as f:
        f.write("# generated by api_extract.py - do not edit\n")
        f.write("API = " + repr(entries) + "\n")
    print("wrote /home/claude/api_data.py")


if __name__ == "__main__":
    main()
