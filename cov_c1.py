"""Parameter covariance of a fit run, and draws from it.

Factored out of band_LT_plot.py / band_LT_gk.py so the band machinery is
written once.  Import needs the venv python (pandas + iminuit).

    from cov_c1 import draws
    pars = draws(n=1500)                 # list of expanded 37-slot vectors
"""
import json
import pathlib

import numpy as np

import fitrun as R, datasets as D

HERE = pathlib.Path(__file__).resolve().parent


def covariance(tag="C1_with_compass"):
    """(p0, C, err) - best fit, covariance of the free model parameters, and
    the per-slot 1-sigma errors (zero for slots that were not free)."""
    p0 = np.load(HERE / f"runs/{tag}/fitpar.npy")
    rec = json.load(open(HERE / f"runs/{tag}/summary.json"))
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
    err = np.zeros(len(p0))
    for k, idx in enumerate(R.FREE):
        err[idx] = np.sqrt(max(C[k, k], 0.0))
    return p0, C, err


def draws(n=1500, tag="C1_with_compass", seed=20260930):
    """n parameter vectors sampled from N(p0, C), clipped to the fit bounds."""
    p0, C, _ = covariance(tag)
    LO, HI = R.bounds()
    rng = np.random.default_rng(seed)
    x = np.clip(rng.multivariate_normal(p0[R.FREE], C, size=n),
                LO[R.FREE], HI[R.FREE])
    return p0, [R.expand(d) for d in x]
