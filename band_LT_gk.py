"""sigma_L, sigma_T and R = L/T vs Q^2 at one (x_B, |t|): amplitude-model band
plus the GK 2024 curve, and a mark showing where data actually constrain it.

Same band machinery as band_LT_plot.py (draws from the parameter covariance,
16-84 percentiles), with two additions:

* the GK_2024 curve, precomputed by ~/gk_work/compass/gk_curve.py (that script
  owns libGKPi0; here we only read its JSON).  If the file is missing it is
  produced on the fly.
* the Q^2 interval where the fit is actually constrained at this (x_B, |t|):
  the cross-section points in the fit within |dx_B| < 0.04 and |d|t|| < 0.08.
  Outside it the band is the functional form talking, not the data.

Run with the venv python (pandas + iminuit):

    ~/.venv/bin/python3 band_LT_gk.py --xB 0.2 --t 0.3
"""
import argparse
import csv
import json
import pathlib
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import fitrun as R, datasets as D, amplitudes as amp

HERE = pathlib.Path(__file__).resolve().parent
COMPASS_DIR = pathlib.Path.home() / "gk_work" / "compass"
TAG = "C1_with_compass"

ap = argparse.ArgumentParser()
ap.add_argument("--xB", type=float, default=0.2)
ap.add_argument("--t", type=float, default=0.3, help="|t| in (GeV/c)^2")
ap.add_argument("--q2-min", type=float, default=1.0)
ap.add_argument("--q2-max", type=float, default=6.0)
ap.add_argument("--ndraws", type=int, default=1500)
ap.add_argument("--nq2-gk", type=int, default=11)
ap.add_argument("--no-gk", action="store_true")
ap.add_argument("--label", default="VPK model",
                help="name of the amplitude model in the legend and title")
ap.add_argument("--out", type=pathlib.Path, default=None)
args = ap.parse_args()
XB, MT = args.xB, args.t

# ---------------------------------------------------------------- the band
p0 = np.load(HERE / f"runs/{TAG}/fitpar.npy")
rec = json.load(open(HERE / f"runs/{TAG}/summary.json"))
keys, lamd = rec["fitted"], rec["norms"]
nn = [k for k in keys if k in R.NORMS]
nf = len(R.FREE)


def resid(x):
    q = R.expand(x[:nf]); lam = dict(zip(nn, x[nf:])); r = []
    for k in keys:
        s = D.BY[k]; L = lam.get(k, 1.0)
        for row in s["rows"]:
            v = s["predict"](q, row)
            r.append((row[4] - L * v) / row[5] if v is not None else 5.0)
    for k in nn:
        r.append((lam[k] - 1.0) / R.NORMS[k])
    for xx in (0.05, 0.35, 1.00):
        sl = R.slopes(q, xx)
        for n in R.BL:
            r.append(R.W * min(0.0, sl[n]))
    return np.array(r)


x0 = np.concatenate([p0[R.FREE], [lamd[k] for k in nn]])
f0 = resid(x0)
J = np.zeros((len(f0), len(x0)))
for i in range(len(x0)):
    h = 1e-6 * max(abs(x0[i]), 1e-3); xp = x0.copy(); xp[i] += h
    J[:, i] = (resid(xp) - f0) / h
C = np.linalg.pinv(J.T @ J)[:nf, :nf]
LO, HI = R.bounds()
rng = np.random.default_rng(20260929)
pars = [R.expand(d) for d in
        np.clip(rng.multivariate_normal(x0[:nf], C, size=args.ndraws),
                LO[R.FREE], HI[R.FREE])]

Q2 = np.linspace(args.q2_min, args.q2_max, 26)
band = {k: np.full((3, len(Q2)), np.nan) for k in ("L", "T", "R")}
for j, q2 in enumerate(Q2):
    sL, sT, rr = [], [], []
    for q in pars:
        s = amp.structure(q, "pi0p", -MT, XB, q2)
        if s is None:
            continue
        sL.append(s["L"]); sT.append(s["T"]); rr.append(s["L"] / s["T"])
    if not sL:
        continue
    for key, a in (("L", sL), ("T", sT), ("R", rr)):
        a = np.array(a)
        band[key][:, j] = (np.median(a), np.percentile(a, 16), np.percentile(a, 84))

# ------------------------------------------------------------- GK 2024
gk = None
if not args.no_gk:
    gkf = COMPASS_DIR / "results" / f"gk_curve_xB{XB:g}_t{MT:g}_GK_2024.json"
    if not gkf.exists():
        subprocess.run([sys.executable, str(COMPASS_DIR / "gk_curve.py"),
                        "--xB", str(XB), "--t", str(MT),
                        "--q2-min", str(args.q2_min), "--q2-max", str(args.q2_max),
                        "--nq2", str(args.nq2_gk), "--workers", str(args.nq2_gk)],
                       cwd=COMPASS_DIR, check=True)
    gk = json.loads(gkf.read_text())
    gk_rows = [r for r in gk["rows"] if r["sigmaL"] is not None]

# --------------------------------------------- where the data constrain it
cov = []
with open(HERE / "db_csv/cross_sections.csv") as f:
    for r in csv.DictReader(f):
        try:
            if r["in_fit"] != "yes":
                continue
            xb, t, q2 = float(r["xB"]), abs(float(r["t"])), float(r["Q2"])
        except Exception:
            continue
        if abs(xb - XB) < 0.04 and abs(t - MT) < 0.08:
            cov.append((q2, r["exp"]))
q2cov = [c[0] for c in cov]
exps = sorted({c[1] for c in cov})

# ------------------------------------------------------------------ plot
C_MOD = "#1a5fb4"
C_GK = "#a50f15"
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.8))
for a, key, lab, log in ((ax[0], "L", r"$d\sigma_L/d|t|$   [nb/(GeV/c)$^2$]", True),
                         (ax[1], "T", r"$d\sigma_T/d|t|$   [nb/(GeV/c)$^2$]", True),
                         (ax[2], "R", r"$R = \sigma_L/\sigma_T$", False)):
    m, lo, hi = band[key]
    if q2cov:
        a.axvspan(min(q2cov), max(q2cov), color="0.5", alpha=0.10, lw=0, zorder=0)
    a.fill_between(Q2, lo, hi, color=C_MOD, alpha=0.22, lw=0,
                   label=f"{args.label}, 16-84%")
    a.plot(Q2, m, "-", color=C_MOD, lw=2, label=f"{args.label} median")
    if gk:
        gq = [r["Q2"] for r in gk_rows]
        gy = {"L": [r["sigmaL"] for r in gk_rows],
              "T": [r["sigmaT"] for r in gk_rows]}
        y = gy[key] if key in gy else [r["sigmaL"] / r["sigmaT"] for r in gk_rows]
        a.plot(gq, y, "--o", color=C_GK, lw=1.8, ms=3.5, label="GK 2024")
    if log:
        a.set_yscale("log")
    else:
        a.set_ylim(0, None)
    a.set_xlabel(r"$Q^2$   [(GeV/c)$^2$]")
    a.set_ylabel(lab)
    a.grid(alpha=.25, lw=.5)
    a.set_xlim(args.q2_min, args.q2_max)
ax[0].legend(frameon=False, fontsize=8.5, loc="lower left")
cov_txt = (f"grey: Q$^2$ covered by fitted data at this $(x_B, |t|)$ — "
           f"{len(cov)} points, {', '.join(exps)}" if cov else
           "no fitted cross-section point at this $(x_B, |t|)$: "
           "the whole curve is functional form")
fig.suptitle(r"$\gamma^* p \to \pi^0 p'$ at "
             rf"$x_B = {XB:g}$, $|t| = {MT:g}$ (GeV/c)$^2$" "\n"
             f"band = {args.label} fit uncertainty at fixed functional form "
             f"({args.ndraws} draws from the parameter covariance, 16–84%)\n"
             + cov_txt, fontsize=9.5, y=0.995, va="top")
fig.text(0.995, 0.012, f"fit {TAG} / amp2609   vs   libGKPi0 GK_2024",
         ha="right", va="bottom", fontsize=6.5, color="0.45")
fig.tight_layout(rect=[0, 0.02, 1, 0.88])
fig.subplots_adjust(top=0.80)
out = args.out or HERE / "figures" / f"band_LT_gk_xB{XB:g}_t{MT:g}.png"
fig.savefig(out, dpi=150)
print("->", out)

print(f"\n  Q^2 |      amp sigma_L        |       amp sigma_T       |    amp R"
      f"          |  GK sigma_L  sigma_T     R")
for j, q2 in enumerate(Q2):
    if round(q2, 6) not in [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]:
        continue
    g = ""
    if gk:
        gq = np.array([r["Q2"] for r in gk_rows])
        gl = np.interp(q2, gq, [r["sigmaL"] for r in gk_rows])
        gt = np.interp(q2, gq, [r["sigmaT"] for r in gk_rows])
        g = f" | {gl:9.3f} {gt:9.2f} {gl/gt:7.3f}"
    print(f"  {q2:3.0f} | {band['L'][0,j]:8.3g} [{band['L'][1,j]:7.3g},{band['L'][2,j]:8.3g}] |"
          f" {band['T'][0,j]:7.4g} [{band['T'][1,j]:6.4g},{band['T'][2,j]:7.4g}] |"
          f" {band['R'][0,j]:.3f} [{band['R'][1,j]:.3f},{band['R'][2,j]:.3f}]{g}")
if cov:
    print(f"\n  data coverage at x_B = {XB} +- 0.04, |t| = {MT} +- 0.08 : "
          f"Q^2 {min(q2cov):.2f} - {max(q2cov):.2f}, {len(cov)} points from {', '.join(exps)}")
