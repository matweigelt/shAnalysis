"""Attribution sweep (v2.4.2): append
    Developed by Matthias Weigelt with the help of Claude (Fable 5).
to (a) the MAIN help block of every .m file in the toolbox and (b) the
footer of every html/ doc page. Existing dated provenance stamps stay
(version history); files already containing the statement are skipped.
Idempotent. Claude (Fable 5), 2026-08-07.
"""
import glob
import os
import re
import sys

ROOT = "/home/claude/shx_build/shAnalysis"
STMT = "Developed by Matthias Weigelt with the help of Claude (Fable 5)."


def find_help_block(lines):
    """Return (insert_at, indent) for the end of the main help block, or None.
    Main help block: contiguous %-lines directly after the first
    function/classdef statement (incl. ... continuations), or the leading
    %-block for pure comment/script files (Contents.m)."""
    i = 0
    n = len(lines)
    first = lines[0].lstrip() if n else ""
    if re.match(r"(function|classdef)\b", first):
        # skip the (possibly continued) declaration
        while i < n and lines[i].rstrip().endswith("..."):
            i += 1
        i += 1
    elif first.startswith("%"):
        i = 0
    else:
        return None
    if i >= n or not lines[i].lstrip().startswith("%"):
        return None
    start = i
    while i < n and lines[i].lstrip().startswith("%"):
        i += 1
    indent = re.match(r"\s*", lines[start]).group(0)
    return i, indent


def process_m(path):
    with open(path, encoding="utf-8") as f:
        src = f.read()
    if "Developed by Matthias Weigelt" in src:   # incl. wrapped variants
        return "skip (present)"
    lines = src.split("\n")
    hb = find_help_block(lines)
    if hb is None:
        return "SKIP (no help block)"
    end, indent = hb
    lines.insert(end, f"{indent}%   {STMT}")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    return "ok"


def process_html(path):
    with open(path, encoding="utf-8") as f:
        src = f.read()
    if "Developed by Matthias Weigelt" in src:
        return "skip (present)"
    foot = f'<p class="provenance">{STMT}</p>\n'
    if "</body>" in src:
        src = src.replace("</body>", foot + "</body>", 1)
    else:
        src = src.rstrip() + "\n" + foot
    with open(path, "w", encoding="utf-8") as f:
        f.write(src)
    return "ok"


def main():
    counts = {}
    skipped = []
    for p in sorted(glob.glob(os.path.join(ROOT, "**", "*.m"),
                              recursive=True)):
        r = process_m(p)
        counts[r] = counts.get(r, 0) + 1
        if r.startswith("SKIP"):
            skipped.append(os.path.relpath(p, ROOT))
    for p in sorted(glob.glob(os.path.join(ROOT, "html", "*.html"))):
        r = "html " + process_html(p)
        counts[r] = counts.get(r, 0) + 1
    print(counts)
    for s in skipped:
        print("  no help block:", s)


if __name__ == "__main__":
    sys.exit(main())
