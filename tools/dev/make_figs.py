"""make_figs.py - regenerate guide figures from the Python port + real data.

Outputs to /home/claude/guide_assets: d01.png (triangle + spectrum, real
ITSG 2008-04), d04_gains.png (Gaussian / fan / DDK3 gain triangles),
d04_diff.png (NEW: signed-difference triangle on the real GRACE-FO -
GRACE field, full sine wing incl. sectorals - the v2.4.2 fix).

Developed by Matthias Weigelt with the help of Claude (Fable 5), 2026-08-07.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, "/home/claude/shx_port")
import shx_port as sp  # noqa: E402

DATA = "/home/claude/shx_build/shAnalysis/tests/test_data"
OUT = "/home/claude/guide_assets"

g1 = sp.read_gfc(f"{DATA}/ITSG-Grace2018_n60_2008-04.gfc")
g2 = sp.read_gfc(f"{DATA}/ITSG-Grace_operational_n60_2025-12.gfc")

# ------------------------------------------------- D01 triangle + spectrum
C = g1["C"].copy(); S = g1["S"].copy()
Cd = C.copy(); Cd[0, 0] = 0; Cd[2, 0] = 0            # display: drop C00, C20
img = sp.triangle_image(Cd, S)
n, amp, err = sp.degree_rms(C, S, g1["sigmaC"], g1["sigmaS"])
kaula = 1e-5 / np.maximum(n, 1) ** 2
cross = np.argmax((amp[2:] < err[2:])) + 2 if np.any(amp[2:] < err[2:]) else None

fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.5, 2.15), dpi=271)
nmax = C.shape[0] - 1
im = a1.imshow(img, extent=[-nmax - .5, nmax + .5, nmax + .5, -.5],
               aspect="equal", cmap="viridis")
a1.set_xlabel(r"$\leftarrow$ order m ($S_{nm}$)   order m ($C_{nm}$) "
              r"$\rightarrow$", fontsize=6)
a1.set_ylabel("degree n", fontsize=6)
a1.set_title("coefficient triangle: ITSG-Grace2018 2008-04 (n60)",
             fontsize=6.5)
a1.tick_params(labelsize=5.5)
cb = fig.colorbar(im, ax=a1, fraction=.046, pad=.03)
cb.set_label(r"log$_{10}$|coefficient|", fontsize=5.5)
cb.ax.tick_params(labelsize=5)

a2.semilogy(n[2:], amp[2:], "-", color="#1a5ca8", lw=1.2, label="signal")
a2.semilogy(n[2:], err[2:], "-", color="#c8571b", lw=1.2,
            label="formal error (real)")
a2.semilogy(n[2:], kaula[2:], "--", color="#666666", lw=.9,
            label=r"Kaula $10^{-5}/n^2$")
if cross:
    a2.axvline(cross, color="#999999", lw=.7, ls=":")
    a2.text(cross + 1, amp[2] * .5, f"crossover n={cross}", fontsize=5.5,
            color="#555555")
a2.set_xlabel("degree n", fontsize=6)
a2.set_ylabel("degree amplitude", fontsize=6)
a2.set_title("spectrum, Kaula + crossover overlays", fontsize=6.5)
a2.tick_params(labelsize=5.5)
a2.legend(fontsize=5.5, frameon=False)
a2.grid(alpha=.3, lw=.4)
fig.tight_layout()
fig.savefig(f"{OUT}/d01.png", dpi=271, bbox_inches="tight")
plt.close(fig)

# --------------------------------------------------- D04 gain triangles
nmaxG = 60
tri = np.tril(np.ones((nmaxG + 1, nmaxG + 1)))
Wg = sp.gaussian_weights(nmaxG, 350)
Ggauss = np.where(tri > 0, np.tile(Wg[:, None], (1, nmaxG + 1)), np.nan)
Gfan = np.where(tri > 0, sp.fan_gain(nmaxG, 350, 200), np.nan)
Wddk = sp.read_wbd(f"{DATA}/Wbd_2-120.a_1d12p_4")
Gddk = sp.ddk_diag_gain(Wddk)[: nmaxG + 1, : nmaxG + 1]
Gddk = np.where(tri > 0, Gddk, np.nan)

fig, axes = plt.subplots(1, 3, figsize=(6.5, 1.95), dpi=271)
for ax, G, tt in zip(axes, [Ggauss, Gfan, Gddk],
                     ["Gaussian 350 km", "fan 350/200 km",
                      "DDK3 (real Wbd diagonal)"]):
    im = ax.imshow(G, extent=[-.5, nmaxG + .5, nmaxG + .5, -.5],
                   aspect="equal", cmap="magma", vmin=0, vmax=1)
    ax.set_title(tt, fontsize=6.5)
    ax.set_xlabel("order m", fontsize=6)
    ax.tick_params(labelsize=5.5)
axes[0].set_ylabel("degree n", fontsize=6)
cb = fig.colorbar(im, ax=axes, fraction=.02, pad=.02)
cb.set_label("gain", fontsize=5.5)
cb.ax.tick_params(labelsize=5)
fig.savefig(f"{OUT}/d04_gains.png", dpi=271, bbox_inches="tight")
plt.close(fig)

# -------------------------- NEW: diff-mode triangle, real data, v2.4.2
nmaxD = 40
n1 = nmaxD + 1
Cd_ = (g2["C"] - g1["C"])[:n1, :n1]
Sd_ = (g2["S"] - g1["S"])[:n1, :n1]
Wg40 = sp.gaussian_weights(nmaxD, 350)
CG = Wg40[:, None] * Cd_
SG = Wg40[:, None] * Sd_
img = sp.triangle_image(CG, SG, ref=(Cd_, Sd_))
a = np.nanmax(np.abs(img))
fig, ax = plt.subplots(figsize=(6.5, 2.35), dpi=271)
im = ax.imshow(img, extent=[-nmaxD - .5, nmaxD + .5, nmaxD + .5, -.5],
               aspect="equal", cmap="RdBu_r", vmin=-a, vmax=a)
ax.set_xlabel(r"$\leftarrow$ order m ($S_{nm}$)     order m ($C_{nm}$) "
              r"$\rightarrow$", fontsize=6)
ax.set_ylabel("degree n", fontsize=6)
ax.set_title("signed difference: Gaussian 350 $-$ raw, real GRACE-FO $-$ "
             "GRACE field (n40)\nfull sine wing incl. sectorals "
             "$S_{nn}$ (v2.4.2 fix)", fontsize=6.5)
ax.tick_params(labelsize=5.5)
cb = fig.colorbar(im, ax=ax, fraction=.03, pad=.03)
cb.set_label("coefficient difference", fontsize=5.5)
cb.ax.tick_params(labelsize=5)
fig.tight_layout()
fig.savefig(f"{OUT}/d04_diff.png", dpi=271, bbox_inches="tight")
plt.close(fig)

print("figures written: d01.png d04_gains.png d04_diff.png")
