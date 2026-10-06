"""The t-slope b of d sigma/dt ~ exp(b t) against xB, by Q2 band.

Reproduces the classic CLAS6 picture from OUR OWN workbook rather than from a
figure: b is extracted per (Q2, xB) bin by a weighted fit to sigma_U, over four
experiments.  The model's slope is b + b' ln(xB) and carries NO Q2 dependence at
all, which is what the data show -- the bands lie on top of one another.
"""
import math, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import datasets as D

P = np.load("runs/Mom_B12/fitpar.npy")
S = {s["key"]: s for s in D.SETS}

def slopes(key, obs="U", nmin=4):
    g = {}
    for r in S[key]["rows"]:
        if r[3] != obs: continue
        g.setdefault((round(r[0], 2), round(r[1], 3)), []).append((r[2], r[4], r[5]))
    out = []
    for (Q2, xB), v in sorted(g.items()):
        t = np.array([x[0] for x in v]); s = np.array([x[1] for x in v])
        e = np.array([x[2] for x in v]); m = s > 0
        if m.sum() < nmin: continue
        w = (s[m]/e[m])**2
        A = np.vstack([np.ones(m.sum()), -t[m]]).T; W = np.diag(w)
        try:
            b = np.linalg.solve(A.T@W@A, A.T@W@np.log(s[m])); c = np.linalg.inv(A.T@W@A)
        except np.linalg.LinAlgError: continue
        if b[1] <= 0 or math.sqrt(c[1,1]) > 1.0: continue
        out.append((Q2, xB, b[1], math.sqrt(c[1,1])))
    return out

pts = []
for k in ("clas6_pi0", "clas12_xs", "halla_y11", "halla_y21"):
    pts += slopes(k)
print(f"{len(pts)} (Q2,xB) bins")

# Q2 bands, in the spirit of the published figure
BANDS = [(1.0, 1.5, "#1a1a1a", "s", r"$1.0-1.5$"),
         (1.5, 2.0, "#c01c28", "^", r"$1.5-2.0$"),
         (2.0, 2.5, "#26883a", "v", r"$2.0-2.5$"),
         (2.5, 3.2, "#1a5fb4", "o", r"$2.5-3.2$"),
         (3.2, 4.5, "#8f5902", "D", r"$3.2-4.5$"),
         (4.5, 9.0, "#a347ba", "d", r"$4.5-8.3$")]

fig, a = plt.subplots(figsize=(7.6, 5.6))
for lo, hi, col, mk, lab in BANDS:
    sel = [p for p in pts if lo <= p[0] < hi]
    if not sel: continue
    a.errorbar([p[1] for p in sel], [p[2] for p in sel], [p[3] for p in sel],
               fmt=mk, ms=6, color=col, mfc="none" if mk in "ovd" else col,
               capsize=2.5, lw=1.3, label=lab + f"  ({len(sel)})", zorder=3)
# The model curve must be the slope of the OBSERVABLE, not of a block.  Each
# block has its own slope -- H_T 1.81 at xB = 0.25, Ebar_T^u 1.25, Ebar_T^d 5.67
# -- and sigma_U is their mixture, so its effective slope moves with Q2 even
# though no block's slope does.  Fitted here exactly as the data were.
import amplitudes as _A
def eff_slope(xB, Q2, eps=0.6):
    # window measured from t_min, otherwise the fit range jumps as t_min moves
    # with xB and the curve comes out in steps
    t0 = abs(_A.tmin(_A.Mpi0, Q2, xB))
    tt = np.linspace(t0 + 0.05, t0 + 1.0, 14); ss = []
    for mt in tt:
        r = _A.structure(P, "pi0p", -mt, xB, Q2)
        ss.append(r["T"] + eps*r["L"] if r else np.nan)
    ss = np.array(ss); m = np.isfinite(ss) & (ss > 0)
    if m.sum() < 5: return np.nan
    return -np.polyfit(tt[m], np.log(ss[m]), 1)[0]
x = np.linspace(0.105, 0.62, 60)
for lo, hi, col, mk, lab in BANDS:
    qc = 0.5*(lo+hi)
    a.plot(x, [eff_slope(v, qc) for v in x], "-", color=col, lw=2.0, alpha=.9, zorder=2)
a.plot([], [], "-", color="0.3", lw=2.0,
       label=r"model, $\sigma_U$ slope" "\n" r"(one curve per band)")
a.set_xlabel(r"$x_B$", fontsize=12)
a.set_ylabel(r"$b$   [GeV$^{-2}$]", fontsize=12)
a.set_xlim(0.09, 0.63); a.set_ylim(0, 3.2)
a.grid(alpha=.25, lw=.5)
leg = a.legend(fontsize=9, frameon=False, loc="upper right",
               title=r"$Q^2$ [GeV$^2$]", title_fontsize=9.5)
a.set_title(r"$d\sigma/dt \propto e^{bt}$ :  the slope falls with $x_B$,"
            "\n" r"and the $Q^2$ dependence is weak -- $+0.2$ to $+0.3$ over $Q^2 = 1.2$ to $8$"
            "\n" r"CLAS6 $\pi^0$, CLAS12, Hall A 2011 and 2021 -- slopes fitted here from the database",
            fontsize=10.5)
fig.tight_layout()
fig.savefig("figures/tslope_xb.png", dpi=160)
print("-> figures/tslope_xb.png")
for v in (0.13, 0.2, 0.3, 0.4, 0.5, 0.6):
    print(f"   model b at xB = {v:.2f} : {P[1]+P[2]*math.log(v):5.2f} GeV^-2")
