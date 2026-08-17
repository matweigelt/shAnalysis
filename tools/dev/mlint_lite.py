"""Pragmatic MATLAB syntax checker (recreated for the v2.4.1 session).
Rules, per the shAnalysis handoff spec:
  R1  bracket balance ((), [], {}) per file
  R2  block/end balance: openers (function/if/for/while/switch/try/parfor/
      arguments/methods/properties/classdef) at bracket depth 0 versus 'end'
      at depth 0; 'end' inside brackets is indexing
  R3  "(expression).method"  -> MATLAB parse error
      classification: walk back from the matching '(' -- if the char before
      it is an identifier char or a closing bracket, it is a call/index
      (legal); otherwise the paren wraps an expression (illegal to dot)
  R4  "shx.fun(...).field"   -> parses, fails at runtime (dot on a package
      function result; package functions are never variables)
  R5  name-value struct declared in an arguments block (opts.Field) but
      missing as the LAST input of the function declaration line ->
      MATLAB:functionValidation:MismatchBetweenBlockAndLine at call time.
      Added 2026-08-18 after exactly this escaped to CI (rtsSmoother
      Lag option, PR #66): mlint_lite was green, every call red.
Tokenizes away comments (%, %{..%}), char/string literals (transpose-aware),
and folds ... line continuations. A tripwire, not a parser.
Usage: python3 mlint_lite.py file1.m [file2.m ...]   (also reads .txt snippets)
Claude (Fable 5), 2026-08-07; R5 added 2026-08-18.
"""
import re
import sys

BLOCK_OPENERS = {"function", "if", "for", "while", "switch", "try",
                 "parfor", "arguments", "methods", "properties", "classdef"}
# these are openers ONLY at statement start (they are also legal identifiers)
STMT_ONLY = {"arguments", "methods", "properties", "function", "classdef"}

OPEN, CLOSE = "([{", ")]}"


def strip_tokens(src):
    """Remove comments and string literals, preserving line structure.
    Returns text of identical length/line layout with stripped chars blanked."""
    out = list(src)
    i, n = 0, len(src)
    in_block_comment = 0
    while i < n:
        c = src[i]
        two = src[i:i + 2]
        if in_block_comment:
            if two == "%}":
                out[i] = out[i + 1] = " "
                i += 2
                in_block_comment -= 1
            else:
                if c != "\n":
                    out[i] = " "
                i += 1
            continue
        if two == "%{" and src[:i].split("\n")[-1].strip() == "":
            in_block_comment += 1
            out[i] = out[i + 1] = " "
            i += 2
            continue
        if c == "%":
            while i < n and src[i] != "\n":
                out[i] = " "
                i += 1
            continue
        if c == '"':                                   # string scalar
            out[i] = " "
            i += 1
            while i < n and src[i] != "\n":
                if src[i] == '"':
                    if src[i:i + 2] == '""':
                        out[i] = out[i + 1] = " "
                        i += 2
                        continue
                    out[i] = " "
                    i += 1
                    break
                out[i] = " "
                i += 1
            continue
        if c == "'":
            prev = src[i - 1] if i else ""
            if re.match(r"[\w\)\]\}\.']", prev):       # transpose
                i += 1
                continue
            out[i] = " "
            i += 1
            while i < n and src[i] != "\n":            # char literal
                if src[i] == "'":
                    if src[i:i + 2] == "''":
                        out[i] = out[i + 1] = " "
                        i += 2
                        continue
                    out[i] = " "
                    i += 1
                    break
                out[i] = " "
                i += 1
            continue
        i += 1
    return "".join(out)


def fold_continuations(txt):
    return re.sub(r"\.\.\.[^\n]*\n", " ", txt)


def check_balance(name, txt):
    """R1 + R2 on tokenized text."""
    bad = []
    depth = 0
    opens = ends = 0
    stack = []
    for ln, line in enumerate(txt.split("\n"), 1):
        stmt_start = True
        j = 0
        while j < len(line):
            c = line[j]
            if c in OPEN:
                depth += 1
                stack.append((c, ln))
                stmt_start = False
                j += 1
                continue
            if c in CLOSE:
                depth -= 1
                if depth < 0:
                    bad.append((name, ln, f"unbalanced '{c}'"))
                    depth = 0
                elif stack:
                    o, oln = stack.pop()
                    if OPEN.index(o) != CLOSE.index(c):
                        bad.append((name, ln,
                                    f"'{c}' closes '{o}' from line {oln}"))
                stmt_start = False
                j += 1
                continue
            if c in ";,":
                if depth == 0:
                    stmt_start = True
                j += 1
                continue
            if c.isspace():
                j += 1
                continue
            m = re.match(r"[A-Za-z_]\w*", line[j:])
            if m:
                w = m.group(0)
                if depth == 0:
                    if w == "end":
                        ends += 1
                    elif w in BLOCK_OPENERS:
                        if w in STMT_ONLY:
                            if stmt_start:
                                opens += 1
                        else:
                            opens += 1
                j += len(w)
                stmt_start = False
                continue
            stmt_start = False
            j += 1
    if depth != 0:
        where = stack[-1][1] if stack else "?"
        bad.append((name, where, f"{depth} unclosed bracket(s)"))
    if opens != ends:
        bad.append((name, 0, f"block imbalance: {opens} openers vs "
                             f"{ends} end(s)"))
    return bad


def check_dot_rules(name, txt):
    """R3 + R4 on tokenized, continuation-folded text."""
    bad = []
    # R4 first: shx.fun( ... ).field
    for mm in re.finditer(r"\bshx\.\w+\s*\(", txt):
        j = mm.end() - 1
        depth = 0
        while j < len(txt):
            if txt[j] == "(":
                depth += 1
            elif txt[j] == ")":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        if j < len(txt) and re.match(r"\)\s*\.\s*[A-Za-z_]", txt[j:j + 8]
                                     if len(txt) - j >= 8 else txt[j:]):
            ln = txt[:mm.start()].count("\n") + 1
            bad.append((name, ln, "shx.fun(...).field -- dot on package "
                                  "function result (runtime error)"))
    # R3: ").method" where the matching "(" wraps an expression
    for mm in re.finditer(r"\)\s*\.\s*[A-Za-z_]\w*", txt):
        j = mm.start()
        depth = 0
        while j >= 0:
            if txt[j] == ")":
                depth += 1
            elif txt[j] == "(":
                depth -= 1
                if depth == 0:
                    break
            j -= 1
        if j < 0:
            continue
        k = j - 1
        while k >= 0 and txt[k] in " \t":
            k -= 1
        prev = txt[k] if k >= 0 else ""
        if re.match(r"[\w\)\]\}]", prev):
            # call or index -> legal; but shx.fun( is caught by R4 above
            continue
        if prev == ".":
            # dynamic field name s.(expr) -> legal MATLAB
            continue
        inner = txt[j + 1:mm.start()].strip()
        if re.fullmatch(r"[A-Za-z_][\w\.]*", inner):
            continue                     # (x).f is legal-ish; skip noise
        ln = txt[:mm.start()].count("\n") + 1
        bad.append((name, ln, f"(expression).method parse error: "
                              f"({inner[:40]}...)" if len(inner) > 40 else
                              f"(expression).method parse error: ({inner})"))
    return bad


def check_arguments_blocks(name, txt):
    """R5 on tokenized, continuation-folded text: every struct used for
    name-value entries in an (input) arguments block must be the last
    input in the function declaration it belongs to."""
    bad = []
    lines = txt.split("\n")
    # function declarations: line no -> ordered input list
    funcs = []                                    # (lineno, [inputs])
    for i, ln in enumerate(lines):
        m = re.match(r"\s*function\b(.*)$", ln)
        if not m:
            continue
        rest = m.group(1)
        rest = rest.split("=", 1)[1] if "=" in rest else rest
        pm = re.search(r"\(([^)]*)\)", rest)
        args = ([a.strip() for a in pm.group(1).split(",") if a.strip()]
                if pm else [])
        funcs.append((i, args))
    if not funcs:
        return bad
    # arguments blocks: attach each to the nearest function above
    i = 0
    while i < len(lines):
        m = re.match(r"\s*arguments\b\s*(\(([^)]*)\))?", lines[i])
        if not m:
            i += 1
            continue
        attrs = (m.group(2) or "").lower()
        start = i
        depth = 0                                 # nested block openers
        structs = []
        i += 1
        while i < len(lines):
            ln = lines[i].strip()
            if re.match(r"\bend\b", ln) and depth == 0:
                break
            sm = re.match(r"([A-Za-z_]\w*)\.[A-Za-z_]\w*", ln)
            if sm:
                structs.append(sm.group(1))
            i += 1
        if "output" not in attrs:
            owner = max((f for f in funcs if f[0] < start),
                        key=lambda f: f[0], default=None)
            if owner:
                for sname in dict.fromkeys(structs):
                    if sname not in owner[1]:
                        bad.append((name, start + 1,
                                    f"arguments block declares name-value "
                                    f"struct '{sname}' but the function "
                                    f"line lacks it (Mismatch"
                                    f"BetweenBlockAndLine at call time)"))
                    elif owner[1] and owner[1][-1] != sname:
                        bad.append((name, start + 1,
                                    f"name-value struct '{sname}' must be "
                                    f"the LAST input of the function line"))
        i += 1
    return bad


def lint_text(name, src):
    txt = fold_continuations(strip_tokens(src))
    return (check_balance(name, txt) + check_dot_rules(name, txt)
            + check_arguments_blocks(name, txt))


def main(paths):
    bad = []
    for p in paths:
        with open(p, encoding="utf-8", errors="replace") as f:
            bad += lint_text(p, f.read())
    for name, ln, msg in bad:
        print(f"{name}:{ln}: {msg}")
    print(f"mlint_lite: {len(paths)} file(s), {len(bad)} finding(s)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
