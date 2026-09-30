"""t-slopes of the VPK model: the fitted parameters and the effective b of
sigma_L, sigma_T and R = sigma_L/sigma_T.

Two different things are printed, and they are not the same number:

1. the slope parameters of the flavour GFFs, which enter the *amplitude* as
   exp(b_eff * t) with b_eff = b + b' ln x_B.  A cross section is the square of
   an amplitude, so the GFF contributes 2*b_eff to the slope of sigma.

2. the effective slope actually seen in a cross section,
   b = -d ln sigma / d|t|, fitted over a |t| window.  It differs from 2*b_eff
   because of the explicit t' = t - t_min factors in the model:
       sigma_L  ~ t' * exp(2 b_L t)            (vanishes at t_min)
       sigma_T  ~ |H_T|^2 exp(2 b_HT t)  +  t' |Ebar_T|^2 exp(2 b_ET t)
   so b runs with |t| and with the H_T / Ebar_T mixture.

Uncertainties come from the same covariance draws as band_LT_gk.py.

    ~/.venv/bin/python3 tslope.py --xB 0.15 0.2 0.3 0.4 --q2 2.0
"""
import argparse
import json
import pathlib

import numpy as np

import fitrun as R, datasets as D, amplitudes as amp

HERE = pathlib.Path(__file__).resolve().parent
TAG = "C1_with_compass"
GK_JSON = (pathlib.Path.home() / "gk_work/compass/results"
           / "LT_GK_2024_t015-02-03-04-05.json")

ap = argparse.ArgumentParser()
ap.add_argument("--xB", type=float, nargs="+", default=[0.15, 0.2, 0.3, 0.4])
ap.add_argument("--q2", type=float, default=2.0)
ap.add_argument("--t-lo", type=float, default=0.15)
ap.add_argument("--t-hi", type=float, default=0.50)
ap.add_argument("--ndraws", type=int, default=600)
ap.add_argument("--tprime", action="store_true",
                help="the window is t' = |t| - |t|min, so every x_B is fitted "
                     "at the same distance from threshold")
args = ap.parse_args()

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
Cfull = np.linalg.pinv(J.T @ J)
C = Cfull[:nf, :nf]
err = np.zeros(len(p0))
for k, idx in enumerate(R.FREE):
    err[idx] = np.sqrt(max(C[k, k], 0.0))

# ------------------------------------------------- 1. the fitted parameters
GFF = [("H_T^u", 1, 2, 3), ("H_T^d", 5, 6, 7),
       ("Ebar_T^u", 9, 10, 11), ("Ebar_T^d", 14, 12, 27),
       ("L (twist-2)", 16, 21, 17)]
print(f"VPK model {TAG}: amplitude slope parameters,  b_eff(x_B) = b + b' ln x_B\n")
bp_lab = "b'"
print(f"  {'GFF':11s} {'b':>16s} {bp_lab:>16s}   b_eff at x_B = "
      + "  ".join(f"{x:g}" for x in args.xB))
for name, ib, ibp, _ in GFF:
    b, bp = p0[ib], p0[ibp]
    eb = f"{b:7.3f} +- {err[ib]:5.3f}" if err[ib] else f"{b:7.3f}  (fixed)"
    ebp = f"{bp:7.3f} +- {err[ibp]:5.3f}" if err[ibp] else f"{bp:7.3f}  (fixed)"
    vals = "  ".join(f"{b + bp*np.log(x):6.2f}" for x in args.xB)
    print(f"  {name:11s} {eb:>16s} {ebp:>16s}    {vals}")
print("\n  (these are amplitude slopes in GeV^-2; a cross section gets 2*b_eff "
      "from the GFF,\n   plus the t' factors, which is why the numbers below differ)")

# --------------------------------------- 2. effective slopes of the sigmas
LO, HI = R.bounds()
rng = np.random.default_rng(20260930)
pars = [R.expand(d) for d in
        np.clip(rng.multivariate_normal(x0[:nf], C, size=args.ndraws),
                LO[R.FREE], HI[R.FREE])]
tt = np.linspace(args.t_lo, args.t_hi, 8)


def slope(p, xB, Q2, key):
    y, x = [], []
    off = abs(amp.tmin(amp.Mpi0, Q2, xB)) if args.tprime else 0.0
    for t in tt + off:
        s = amp.structure(p, "pi0p", -t, xB, Q2)
        if s is None:
            continue
        v = s["L"] if key == "L" else s["T"] if key == "T" else s["L"] / s["T"]
        if v <= 0:
            continue
        y.append(np.log(v)); x.append(t)
    if len(x) < 4:
        return np.nan
    return -np.polyfit(x, y, 1)[0]


var = "t'" if args.tprime else "|t|"
print(f"\neffective slope b = -d ln sigma / d|t| [GeV^-2], "
      f"fitted over {var} = {args.t_lo}-{args.t_hi}, Q^2 = {args.q2:g}\n")
print(f"  {'x_B':>5} | {'sigma_L':>22} | {'sigma_T':>22} | {'R = L/T':>20} | GK 2024 L / T / R")
gk = json.loads(GK_JSON.read_text()) if GK_JSON.exists() else None
for xB in args.xB:
    out = {}
    for key in ("L", "T", "R"):
        v = np.array([slope(p, xB, args.q2, key) for p in pars])
        v = v[np.isfinite(v)]
        out[key] = (np.median(v), np.percentile(v, 16), np.percentile(v, 84)) \
            if len(v) else (np.nan,) * 3
    g = "        -"
    if gk and any(abs(r["xB"] - xB) < 1e-9 for r in gk["rows"]):
        rows = [r for r in gk["rows"]
                if abs(r["xB"] - xB) < 1e-9 and abs(r["Q2"] - args.q2) < 0.35]
        if rows:
            r = min(rows, key=lambda r: abs(r["Q2"] - args.q2))
            sel = [(t, r["sigmaL"][i], r["sigmaT"][i])
                   for i, t in enumerate(gk["t"]) if t >= r["tmin"]]
            if len(sel) >= 4:
                x = [s[0] for s in sel]
                bl = -np.polyfit(x, np.log([s[1] for s in sel]), 1)[0]
                bt = -np.polyfit(x, np.log([s[2] for s in sel]), 1)[0]
                g = f"{bl:6.2f} /{bt:6.2f} /{bl-bt:6.2f}"
    def f(k):
        m, lo, hi = out[k]
        return f"{m:6.2f} [{lo:5.2f},{hi:5.2f}]"
    print(f"  {xB:5g} | {f('L'):>22} | {f('T'):>22} | {f('R'):>20} | {g}")
print("\n  R's slope is b_L - b_T by construction; a positive value means R "
      "falls with |t|.")
