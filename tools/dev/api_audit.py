"""api_audit.py - gate: every API example must match the parsed contract.

Checks per call in each example: positional-argument count within the
[required, total] range of the arguments block, name-value options among
the declared ones, and the number of requested outputs not exceeding the
declared outputs. varargin functions are unbounded with free options;
functions without arguments blocks allow 1..len(inputs) positionals.

Exit code 1 on any finding. Developed by Matthias Weigelt with the help
of Claude (Fable 5), 2026-08-07 (v2.5).
"""
import re
import sys

sys.path.insert(0, "/home/claude")
exec(open("/home/claude/api_data.py").read())     # -> API

INF = 999


def sig_from(e):
    ins = e["inputs"]
    varargin = "varargin" in ins
    args = e.get("args", [])
    pos = [a for a in args if not a["nv"] and a["name"] not in ("obj",)]
    nv = {a["name"] for a in args if a["nv"]}
    if args:
        minp = len([a for a in pos if not a["default"]])
        maxp = INF if varargin else len(pos)
        nvs = None if varargin else nv
    else:
        ins2 = [i for i in ins if i not in ("obj", "opts", "varargin")]
        minp, maxp = (1 if ins2 else 0), (INF if varargin else len(ins2))
        nvs = None
    return minp, maxp, nvs, len(e["outputs"])


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


def check(fn, argstr, nout, minp, maxp, nv, maxout, where, probs):
    args = split_args(argstr)
    npos, nvused = 0, []
    for a in args:
        m = re.match(r"^([A-Za-z]\w*)\s*=\s*[^=]", a)
        if m and not a.startswith("@"):
            nvused.append(m.group(1))
        elif not nvused:
            npos += 1
    if not (minp <= npos <= maxp):
        probs.append(f"{where}: {fn} positional {npos} not in "
                     f"[{minp},{maxp}]")
    if nv is not None:
        for n in nvused:
            if n not in nv:
                probs.append(f"{where}: {fn} unknown option {n}")
    if nout > maxout > 0:
        probs.append(f"{where}: {fn} requests {nout} outputs, "
                     f"declares {maxout}")


def main():
    sig = {}
    for e in API["shLowLevel"]:
        sig["shLowLevel." + e["name"]] = sig_from(e)
    for e in API["root"]:
        sig[e["name"]] = sig_from(e)
    meth = {}
    for c in API["classes"]:
        for m in c["methods"]:
            meth.setdefault((m["name"], m.get("static", False)),
                            []).append((c["name"],) + sig_from(m))
            if m.get("static"):
                sig[c["name"] + "." + m["name"]] = sig_from(m)
    probs = []

    def audit(name, ex):
        ex = re.sub(r"\.\.\.\s*\n\s*", " ", ex)     # join continuations
        for mo in re.finditer(r"([\w.]+)\s*\(", ex):
            fn = mo.group(1)
            i, d, buf = mo.end(), 1, ""
            while i < len(ex) and d:
                ch = ex[i]
                if ch == "(":
                    d += 1
                if ch == ")":
                    d -= 1
                if d:
                    buf += ch
                i += 1
            # requested outputs: look back for [a, b] = or a =
            pre = ex[:mo.start()].rstrip()
            nout = 1
            mo2 = re.search(r"\[([^\]]*)\]\s*=\s*$", pre)
            if mo2:
                nout = len([x for x in mo2.group(1).split(",")
                            if x.strip()])
            elif not re.search(r"=\s*$", pre):
                nout = 0
            if fn in sig:
                check(fn, buf, nout, *sig[fn], name, probs)
                continue
            mm = re.match(r"^(\w+)\.(\w+)$", fn)
            if mm and mm.group(1) not in ("shLowLevel",):
                key = (mm.group(2), False)
                if key in meth:
                    trial = []
                    for cand in meth[key]:
                        p2 = []
                        check(fn, buf, nout, *cand[1:], name, p2)
                        trial.append(p2)
                    if all(t for t in trial):
                        probs.extend(trial[0])

    for e in API["shLowLevel"]:
        audit("shLowLevel." + e["name"], e["example"])
    for e in API["root"]:
        audit(e["name"], e["example"])
    for c in API["classes"]:
        for m in c["methods"]:
            audit(c["name"] + "." + m["name"], m["example"])
    if probs:
        print("api_audit: %d finding(s)" % len(probs))
        for p in probs:
            print(" -", p)
        sys.exit(1)
    print("api_audit: 0 finding(s)")


if __name__ == "__main__":
    main()
