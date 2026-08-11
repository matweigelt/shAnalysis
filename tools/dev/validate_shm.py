"""validate_shm.py - the GRAVIS/GRACE SHM format, parsed in Python first.

GravIS Level-2B products (mean fields, GIA rate models, monthly
solutions) do NOT use the ICGEM gfc layout. They are:

    <YAML header>
    # End of YAML header
    GRCOF2  n  m  C  S  sigC  sigS  t0  t1  flags     (a FIELD)
    GRDOTA  n  m  Cdot  Sdot                          (a RATE, 1/yr)

GM and the reference radius live in the YAML as
earth_gravity_param/mean_equator_radius. GRDOTA records carry no
sigmas, and the two record types can be told apart only by the keyword -
which is what makes a gfc parser silently useless here rather than
loudly wrong.

Checks the parse against the real files and prints the values pinned in
tests/testCorrectness.m.

Developed by Matthias Weigelt with the help of Claude (Opus 5),
2026-08-11 (v3.3.0).
"""
import gzip
import re
import sys


def read_shm(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", errors="replace") as fh:
        lines = fh.read().split("\n")
    GM = R = None
    kind = None
    C, S, sC, sS = {}, {}, {}, {}
    nmax = 0
    inhead = True
    for i, L in enumerate(lines):
        if inhead:
            m = re.search(r"value\s*:\s*([-\d.Ee+]+)", L)
            if m:
                # the value line belongs to the preceding long_name block
                ctx = " ".join(lines[max(0, i - 4):i])
                if "gravitational constant" in ctx:
                    GM = float(m.group(1))
                elif "equator radius" in ctx:
                    R = float(m.group(1))
            if L.startswith("# End of YAML header"):
                inhead = False
            continue
        t = L.split()
        if not t or t[0] not in ("GRCOF2", "GRDOTA"):
            continue
        kind = t[0] if kind is None else kind
        n, m_ = int(t[1]), int(t[2])
        nmax = max(nmax, n)
        C[(n, m_)] = float(t[3].replace("D", "E"))
        S[(n, m_)] = float(t[4].replace("D", "E"))
        if t[0] == "GRCOF2" and len(t) >= 7:
            sC[(n, m_)] = float(t[5].replace("D", "E"))
            sS[(n, m_)] = float(t[6].replace("D", "E"))
    return dict(kind=kind, GM=GM, R=R, nmax=nmax, C=C, S=S, sC=sC, sS=sS)


def main(mean_path, gia_path):
    M = read_shm(mean_path)
    print("MEAN  kind=%s nmax=%d GM=%.10E R=%.10E records=%d"
          % (M["kind"], M["nmax"], M["GM"] or 0, M["R"] or 0, len(M["C"])))
    print("      C00=%.12E  C20=%.12E  sigC20=%.4E"
          % (M["C"][(0, 0)], M["C"][(2, 0)], M["sC"].get((2, 0), float("nan"))))
    assert M["kind"] == "GRCOF2"
    assert abs(M["C"][(0, 0)] - 1.0) < 1e-15, "C00 of a field must be 1"
    assert abs(M["C"][(2, 0)] + 4.84165e-4) < 1e-8, "C20 out of range"
    assert (M["nmax"] + 1) * (M["nmax"] + 2) // 2 >= len(M["C"])

    G = read_shm(gia_path)
    print("GIA   kind=%s nmax=%d GM=%.10E R=%.10E records=%d"
          % (G["kind"], G["nmax"], G["GM"] or 0, G["R"] or 0, len(G["C"])))
    print("      C20dot=%.12E  C21dot=%.12E  S21dot=%.12E"
          % (G["C"][(2, 0)], G["C"][(2, 1)], G["S"][(2, 1)]))
    assert G["kind"] == "GRDOTA"
    assert not G["sC"], "a rate model carries no sigmas"
    # a GIA rate is tiny per year but must not be zero
    assert 1e-12 < abs(G["C"][(2, 0)]) < 1e-9, G["C"][(2, 0)]
    # degrees 0 and 1 are zero by construction in this model
    for k in [(0, 0), (1, 0), (1, 1)]:
        assert G["C"][k] == 0.0 and G["S"][k] == 0.0, k
    print("      degree 0/1 are zero as expected")
    print("\nvalidate_shm: format understood, both record types parse")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
