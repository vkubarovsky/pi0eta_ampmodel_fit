"""The De Masi sin(phi) moment: published, re-extracted from the phi
distributions we actually fit, and the model."""
import math, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import datasets as D, amplitudes as A

P = np.load("runs/Rosen_B/fitpar.npy"); LAM = 1.0630
S = {s["key"]: s for s in D.SETS}

def extract(rows):
    n = d = 0.0
    for Q2, xB, mt, _, v, e, eps, _, phi in rows:
        f = math.sin(phi); w = 1/e**2
        n += w*f*v; d += w*f*f
    return n/d, 1/math.sqrt(d)

B = {}
for r in S["bsa_demasi_phi"]["rows"]: B.setdefault(r[7], []).append(r)
rows = []
for b, rr in B.items():
    Q2 = np.mean([r[0] for r in rr]); xB = np.mean([r[1] for r in rr])
    mt = np.mean([r[2] for r in rr]); a, e = extract(rr)
    rows.append((Q2, xB, mt, a, e))
rows.sort()

PUB = sorted((r[0], r[1], r[2], r[4], r[5]) for r in S["bsa_demasi"]["rows"])
EDGES = [(1.0, 1.6), (1.6, 2.1), (2.1, 2.6), (2.6, 3.2), (3.2, 4.2)]
fig, ax = plt.subplots(1, 5, figsize=(19, 4.2), sharey=True)
for a_, (lo, hi) in zip(ax, EDGES):
    sel = [r for r in rows if lo <= r[0] < hi]
    psel = [r for r in PUB if lo <= r[0] < hi]
    if psel:
        a_.errorbar([r[2] for r in psel], [r[3] for r in psel], [r[4] for r in psel],
                    fmt="s", ms=9, mfc="none", color="#c01c28", capsize=3, lw=1.3,
                    label="published moment", zorder=4)
    a_.errorbar([r[2] for r in sel], [r[3] for r in sel], [r[4] for r in sel],
                fmt="o", ms=5, color="k", capsize=2, lw=1.2,
                label=r"re-extracted from $\phi$", zorder=5)
    mt = np.linspace(0.1, 1.55, 120)
    Q2 = np.mean([r[0] for r in sel]); xB = np.mean([r[1] for r in sel])
    y = [A.bsa_sinphi(P, "pi0p", -x, xB, Q2, 5.776) for x in mt]
    y = [np.nan if v is None else v for v in y]
    a_.plot(mt, [LAM*v for v in y], "-", color="#1a5fb4", lw=2.2,
            label=r"model $\times\,\lambda$")
    a_.plot(mt, y, "--", color="#1a5fb4", lw=1.4, alpha=.8, label="model")
    a_.set_title(rf"$Q^2 = {Q2:.2f}$ GeV$^2$,  $x_B = {xB:.2f}$", fontsize=11)
    a_.set_xlabel(r"$|t|$   [GeV$^2$]"); a_.grid(alpha=.25, lw=.5)
    a_.axhline(0, color="k", lw=.8); a_.set_xlim(0, 1.6)
ax[0].set_ylabel(r"$A_{LU}^{\sin\phi}$"); ax[0].set_ylim(-0.04, 0.20)
ax[0].legend(fontsize=8.5, frameon=False, loc="upper right")
fig.suptitle(r"CLAS6 $\pi^0$ beam-spin asymmetry, De Masi PRC 77 042201(R) — "
             r"the paper shows plots only, so the fit uses the $\phi$ distributions; "
             r"this is the moment recovered from them", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig("figures/demasi_moments.png", dpi=150)
print("-> figures/demasi_moments.png")
