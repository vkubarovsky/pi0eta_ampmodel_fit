"""sigma_L, sigma_T and R = L/T vs Q^2 at one (x_B, |t|): amplitude-model bands
plus the GK 2024 curve, and a mark showing where data actually constrain it.

--tags takes one or more run tags and draws a band for each, so the model before
and after the Htilde block can be compared on the same axes.  The covariance is
recomputed per tag from the Jacobian at that tag's own minimum.

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
ap.add_argument("--tags", default="C1_with_compass",
                help="comma-separated run tags; one band each")
ap.add_argument("--labels", default=None,
                help="comma-separated legend labels, one per tag")
ap.add_argument("--out", type=pathlib.Path, default=None)
args = ap.parse_args()
XB, MT = args.xB, args.t

# ---------------------------------------------------------------- the bands
TAGS = [t.strip() for t in args.tags.split(",") if t.strip()]
LABELS = ([l.strip() for l in args.labels.split(",")] if args.labels
          else ["VPK model" if t == "C1_with_compass" else t for t in TAGS])
Q2 = np.linspace(args.q2_min, args.q2_max, 26)


def band_for(tag):
    """Covariance from the Jacobian at THIS tag's minimum, then draws.

    fitrun reads NPAR at import time, so a run with the Htilde slots and one
    without cannot share a process.  R.FREE and R.bounds() are rebuilt here from
    the length of the tag's own parameter vector.
    """
    p0 = np.load(HERE / f"runs/{tag}/fitpar.npy")
    rec = json.load(open(HERE / f"runs/{tag}/summary.json"))
    keys, lamd = rec["fitted"], rec["norms"]
    nn = [k for k in keys if k in R.NORMS]

    npar = len(p0)
    frozen = {23, 28, 29, 30, 31, 32, 33} | ({25, 42} if npar > 37 else set())
    FREE = [i for i in range(npar) if i not in set(R.TIES) | frozen]
    nf = len(FREE)
    BLs = {k: v for k, v in R.BL.items() if max(v) < npar}

    def expand(x):
        q = np.zeros(npar); q[FREE] = x; q[23] = 0.0
        if npar > 42: q[42] = p0[42]
        for t_, s_ in R.TIES.items(): q[t_] = q[s_]
        return q

    def resid(x):
        q = expand(x[:nf]); lam = dict(zip(nn, x[nf:])); r = []
        for k in keys:
            s = D.BY[k]; L = lam.get(k, 1.0)
            for row in s["rows"]:
                v = s["predict"](q, row)
                r.append((row[4] - L * v) / row[5] if v is not None else 5.0)
        for k in nn:
            r.append((lam[k] - 1.0) / R.NORMS[k])
        for xx in (0.05, 0.35, 1.00):
            L_ = np.log(xx)
            for b, bp in BLs.values():
                r.append(R.W * min(0.0, q[b] + q[bp] * L_))
        return np.array(r)

    x0 = np.concatenate([p0[FREE], [lamd[k] for k in nn]])
    f0 = resid(x0)
    J = np.zeros((len(f0), len(x0)))
    for i in range(len(x0)):
        h = 1e-6 * max(abs(x0[i]), 1e-3); xp = x0.copy(); xp[i] += h
        J[:, i] = (resid(xp) - f0) / h
    C = np.linalg.pinv(J.T @ J)
    C = 0.5*(C + C.T)                      # pinv of a numerical J^T J is not
    C = C[:nf, :nf]                        # exactly symmetric; numpy complains
    w = np.linalg.eigvalsh(C)
    if w.min() < 0:
        # a few tiny negative eigenvalues are round-off on near-degenerate
        # directions; clip them rather than let multivariate_normal warn and
        # silently do the same thing
        print(f"   [{tag}] covariance had {int((w < 0).sum())} negative "
              f"eigenvalues, smallest {w.min():.2e} vs largest {w.max():.2e}; "
              f"clipped to zero")
        V = np.linalg.eigh(C)[1]
        C = V @ np.diag(np.clip(w, 0.0, None)) @ V.T
    LOb, HIb = R.bounds()
    LOb = np.resize(LOb, npar) if len(LOb) < npar else LOb[:npar]
    HIb = np.resize(HIb, npar) if len(HIb) < npar else HIb[:npar]
    rng = np.random.default_rng(20260929)
    pars = [expand(d) for d in
            np.clip(rng.multivariate_normal(x0[:nf], C, size=args.ndraws),
                    LOb[FREE], HIb[FREE])]

    out = {k: np.full((3, len(Q2)), np.nan) for k in ("L", "T", "R")}
    for j, q2 in enumerate(Q2):
        sL, sT, rr = [], [], []
        for q in pars:
            sv = amp.structure(q, "pi0p", -MT, XB, q2)
            if sv is None: continue
            sL.append(sv["L"]); sT.append(sv["T"]); rr.append(sv["L"] / sv["T"])
        if not sL: continue
        for key, a in (("L", sL), ("T", sT), ("R", rr)):
            a = np.array(a)
            out[key][:, j] = (np.median(a), np.percentile(a, 16),
                              np.percentile(a, 84))
    return out


BANDS = [band_for(t) for t in TAGS]
band = BANDS[0]

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
C_SET = ["#7a7a7a", "#1a5fb4", "#1b7837"]
LS_SET = ["--", "-", "-."]
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.8))
for a, key, lab, log in ((ax[0], "L", r"$d\sigma_L/d|t|$   [nb/(GeV/c)$^2$]", True),
                         (ax[1], "T", r"$d\sigma_T/d|t|$   [nb/(GeV/c)$^2$]", True),
                         (ax[2], "R", r"$R = \sigma_L/\sigma_T$", False)):
    if q2cov:
        a.axvspan(min(q2cov), max(q2cov), color="0.5", alpha=0.10, lw=0, zorder=0)
    for bi, (bd, lb) in enumerate(zip(BANDS, LABELS)):
        m, lo, hi = bd[key]
        col, ls = C_SET[bi % len(C_SET)], LS_SET[bi % len(LS_SET)]
        a.fill_between(Q2, lo, hi, color=col, alpha=0.20, lw=0,
                       label=f"{lb}, 16-84%")
        a.plot(Q2, m, ls, color=col, lw=2, label=f"{lb} median")
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
             f"bands = fit uncertainty at fixed functional form "
             f"({args.ndraws} draws from the parameter covariance, 16–84%)\n"
             + cov_txt, fontsize=9.5, y=0.995, va="top")
fig.text(0.995, 0.012, f"fits {' , '.join(TAGS)}   vs   libGKPi0 GK_2024",
         ha="right", va="bottom", fontsize=6.5, color="0.45")
fig.tight_layout(rect=[0, 0.02, 1, 0.88])
fig.subplots_adjust(top=0.80)
out = args.out or HERE / "figures" / f"band_LT_gk_xB{XB:g}_t{MT:g}.png"
fig.savefig(out, dpi=150)
print("->", out)

for bd, lb in zip(BANDS, LABELS):
    print(f"\n=== {lb} ===")
    print(f"  Q^2 |        sigma_L          |        sigma_T          |"
          f"        R            |  GK sigma_L  sigma_T     R")
    for j, q2 in enumerate(Q2):
        if round(q2, 6) not in [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]:
            continue
        g = ""
        if gk:
            gq = np.array([r["Q2"] for r in gk_rows])
            gl = np.interp(q2, gq, [r["sigmaL"] for r in gk_rows])
            gt = np.interp(q2, gq, [r["sigmaT"] for r in gk_rows])
            g = f" | {gl:9.3f} {gt:9.2f} {gl/gt:7.3f}"
        print(f"  {q2:3.0f} | {bd['L'][0,j]:8.3g} [{bd['L'][1,j]:7.3g},{bd['L'][2,j]:8.3g}] |"
              f" {bd['T'][0,j]:7.4g} [{bd['T'][1,j]:6.4g},{bd['T'][2,j]:7.4g}] |"
              f" {bd['R'][0,j]:.3f} [{bd['R'][1,j]:.3f},{bd['R'][2,j]:.3f}]{g}")
if cov:
    print(f"\n  data coverage at x_B = {XB} +- 0.04, |t| = {MT} +- 0.08 : "
          f"Q^2 {min(q2cov):.2f} - {max(q2cov):.2f}, {len(cov)} points from {', '.join(exps)}")
