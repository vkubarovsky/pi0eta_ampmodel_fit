"""sigma_T and sigma_TT against |t| at one (Q^2, x_B), VPK model + data.

sigma_TT is measured directly, so the data points are on the same footing as
the curve.  sigma_T is not: the experiments give sigma_U = sigma_T + eps sigma_L,
so the left panel also carries the model's sigma_U (dashed) next to the
sigma_U points, and sigma_T alone as the solid band.

    ~/.venv/bin/python3 tdep_plot.py --q2 2.24 --xB 0.332
"""
import argparse
import csv
import pathlib

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import amplitudes as amp
from cov_c1 import draws

HERE = pathlib.Path(__file__).resolve().parent

ap = argparse.ArgumentParser()
ap.add_argument("--q2", type=float, default=2.24)
ap.add_argument("--xB", type=float, default=0.332)
ap.add_argument("--tmax", type=float, default=1.8)
ap.add_argument("--ndraws", type=int, default=1200)
ap.add_argument("--dq2", type=float, default=0.12, help="window for the data")
ap.add_argument("--dxb", type=float, default=0.02)
ap.add_argument("--label", default="VPK model")
ap.add_argument("--out", type=pathlib.Path, default=None)
args = ap.parse_args()
Q2, XB = args.q2, args.xB

p0, pars = draws(args.ndraws)
eps0 = amp.epsilon(XB, Q2, 5.75)
tmin = abs(amp.tmin(amp.Mpi0, Q2, XB))
tt = np.linspace(tmin + 1e-3, args.tmax, 60)

band = {k: np.full((3, len(tt)), np.nan) for k in ("T", "TT", "U")}
med = {k: np.full(len(tt), np.nan) for k in ("T", "TT", "U")}
for j, t in enumerate(tt):
    vals = {"T": [], "TT": [], "U": []}
    for q in pars:
        s = amp.structure(q, "pi0p", -t, XB, Q2)
        if s is None:
            continue
        vals["T"].append(s["T"]); vals["TT"].append(s["TT"])
        vals["U"].append(s["T"] + eps0 * s["L"])
    for k, a in vals.items():
        if a:
            a = np.array(a)
            band[k][:, j] = (np.median(a), np.percentile(a, 16), np.percentile(a, 84))

# ------------------------------------------------------------------- data
data = {}
with open(HERE / "db_csv/cross_sections.csv") as f:
    for r in csv.DictReader(f):
        try:
            if r["meson"] != "pi0" or r["target"] != "p":
                continue
            q2, xb, t = float(r["Q2"]), float(r["xB"]), abs(float(r["t"]))
        except Exception:
            continue
        if abs(q2 - Q2) > args.dq2 or abs(xb - XB) > args.dxb or t > args.tmax:
            continue
        d = data.setdefault(r["exp"], {"t": [], "U": [], "eU": [], "TT": [], "eTT": []})
        d["t"].append(t)
        d["U"].append(float(r["s_u"]))
        d["eU"].append(np.hypot(float(r["stat_U"]), float(r["sys_U"] or 0)))
        d["TT"].append(float(r["s_TT"]))
        d["eTT"].append(np.hypot(float(r["stat_TT"]), float(r["sys_TT"] or 0)))

STYLE = {"CLAS6_y12": ("o", "tab:red", "CLAS6 (in the fit)"),
         "HallA_y11": ("s", "tab:green", "Hall A 2011"),
         "CLAS12_y25": ("^", "tab:orange", "CLAS12 (preliminary, not fitted)")}

# ------------------------------------------------------------------- plot
C = "#1a5fb4"
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.8))
for a, key, lab in ((ax[0], "T", r"$\mathrm{d}\sigma_T/\mathrm{d}|t|$"),
                    (ax[1], "TT", r"$\mathrm{d}\sigma_{TT}/\mathrm{d}|t|$")):
    a.fill_between(tt, band[key][1], band[key][2], color=C, alpha=0.22, lw=0,
                   label=f"{args.label}, 16-84%")
    a.plot(tt, band[key][0], "-", color=C, lw=2, label=f"{args.label}: $\\sigma_{{{key}}}$")
    a.axvline(tmin, color="0.6", ls=":", lw=1)
    a.text(tmin, 0.02, r" $|t|_{\min}$", transform=a.get_xaxis_transform(),
           fontsize=8, color="0.45")
    a.set_xlabel(r"$|t|$   [(GeV/c)$^2$]")
    a.set_ylabel(lab + r"   [nb/(GeV/c)$^2$]")
    a.grid(alpha=.25, lw=.5)
    a.set_xlim(0, args.tmax)

ax[0].plot(tt, band["U"][0], "--", color=C, lw=1.4,
           label=rf"{args.label}: $\sigma_U=\sigma_T+\epsilon\sigma_L$ ($\epsilon$={eps0:.2f})")
ax[0].set_yscale("log")
ax[0].set_ylim(0.6*np.nanmin(band["T"][1]), 2.2*np.nanmax(band["U"][2]))
for exp, d in data.items():
    m, c, lab = STYLE.get(exp, ("v", "0.4", exp))
    ax[0].errorbar(d["t"], d["U"], yerr=d["eU"], fmt=m, color=c, ms=5, lw=1,
                   label=lab + r"  ($\sigma_U$)")
    ax[1].errorbar(d["t"], d["TT"], yerr=d["eTT"], fmt=m, color=c, ms=5, lw=1,
                   label=lab)
ax[1].axhline(0, color="0.6", lw=0.8, ls=":")
ax[0].legend(frameon=False, fontsize=8)
ax[1].legend(frameon=False, fontsize=8)
fig.suptitle(rf"$\gamma^* p \to \pi^0 p'$ at $Q^2$ = {Q2:g} (GeV/c)$^2$, "
             rf"$x_B$ = {XB:g}   —   {args.label}, band = fit covariance "
             f"({args.ndraws} draws, 16–84%)", fontsize=10)
fig.text(0.995, 0.012, "fit C1_with_compass / amp2609", ha="right", va="bottom",
         fontsize=6.5, color="0.45")
fig.tight_layout(rect=[0, 0.02, 1, 0.93])
out = args.out or HERE / "figures" / f"tdep_T_TT_Q2{Q2:g}_xB{XB:g}.png"
fig.savefig(out, dpi=150)
print("->", out)

print(f"\n  |t|min = {tmin:.3f},  eps(CLAS6, 5.75 GeV) = {eps0:.3f}")
print("   |t|  |    sigma_T  [16,84]        |   sigma_TT [16,84]          |  |sTT|/sT")
for t in (0.18, 0.25, 0.35, 0.50, 0.80, 1.20, 1.70):
    if t < tmin:
        continue
    j = int(np.argmin(abs(tt - t)))
    T, TT = band["T"][:, j], band["TT"][:, j]
    print(f"  {t:4.2f} | {T[0]:8.1f} [{T[1]:7.1f},{T[2]:8.1f}] | "
          f"{TT[0]:8.1f} [{TT[1]:7.1f},{TT[2]:8.1f}] | {abs(TT[0])/T[0]:6.3f}")
