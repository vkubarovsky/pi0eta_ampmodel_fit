"""Error bands on sigma_L, sigma_T and R = L/T from the parameter covariance.

Linear propagation is wrong for a ratio, so the parameters are sampled from
N(p, C) instead and the band read off the percentiles.  C comes from the
Jacobian at the minimum, the same one the report's error column uses.

Two parameters sit on a bound (b of H_T^u and of the longitudinal block, both
at zero, held there by the requirement that form factors fall with |t|).  A
Gaussian sample walks through that wall, so every draw is clipped back into the
allowed box -- which makes the band correct on the inside and conservative on
the side where the bound bites.
"""
import json, sys
import numpy as np
import fitrun as R, datasets as D, amplitudes as amp

TAG = "C1_with_compass"
p0 = np.load(f"runs/{TAG}/fitpar.npy")
rec = json.load(open(f"runs/{TAG}/summary.json"))
keys, lamd = rec["fitted"], rec["norms"]
nn = [k for k in keys if k in R.NORMS]
nf = len(R.FREE)

def resid(x):
    q = R.expand(x[:nf]); lam = dict(zip(nn, x[nf:])); r = []
    for k in keys:
        s = D.BY[k]; L = lam.get(k, 1.0)
        for row in s["rows"]:
            v = s["predict"](q, row)
            r.append((row[4] - L*v)/row[5] if v is not None else 5.0)
    for k in nn: r.append((lam[k] - 1.0)/R.NORMS[k])
    for xx in (0.05, 0.35, 1.00):
        sl = R.slopes(q, xx)
        for n in R.BL: r.append(R.W*min(0.0, sl[n]))
    return np.array(r)

x0 = np.concatenate([p0[R.FREE], [lamd[k] for k in nn]])
f0 = resid(x0)
J = np.zeros((len(f0), len(x0)))
for i in range(len(x0)):
    h = 1e-6*max(abs(x0[i]), 1e-3); xp = x0.copy(); xp[i] += h
    J[:, i] = (resid(xp) - f0)/h
C = np.linalg.pinv(J.T @ J)[:nf, :nf]          # free model parameters only
LO, HI = R.bounds()

rng = np.random.default_rng(20260929)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
draws = rng.multivariate_normal(x0[:nf], C, size=N)
draws = np.clip(draws, LO[R.FREE], HI[R.FREE])
pars = [R.expand(d) for d in draws]

xB, mt, E = 0.2, 0.3, 10.6
print(f"pi0 p,  xB = {xB},  -t = {mt} GeV^2,  {N} draws from the covariance\n")
print(f"{'Q2':>5} | {'sigma_L [nb/GeV^2]':>26} | {'sigma_T':>26} | {'R = L/T':>22}")
print(f"{'':>5} | {'median   [16%, 84%]':>26} | {'median   [16%, 84%]':>26} | {'median   [16%, 84%]':>22}")
print("-"*88)
for Q2 in (1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0):
    sL, sT, Rr = [], [], []
    s = amp.structure(p0, "pi0p", -mt, xB, Q2)
    for q in pars:
        t = amp.structure(q, "pi0p", -mt, xB, Q2)
        if t is None: continue
        sL.append(t["L"]); sT.append(t["T"]); Rr.append(t["L"]/t["T"])
    if not sL: print(f"{Q2:5.1f} |  вне области"); continue
    f = lambda a: (np.median(a), np.percentile(a, 16), np.percentile(a, 84))
    (mL,loL,hiL), (mT,loT,hiT), (mR,loR,hiR) = f(np.array(sL)), f(np.array(sT)), f(np.array(Rr))
    print(f"{Q2:5.1f} | {mL:8.3g} [{loL:7.3g},{hiL:8.3g}] | {mT:8.3g} [{loT:7.3g},{hiT:8.3g}] |"
          f" {mR:6.3f} [{loR:5.3f},{hiR:6.3f}]")
