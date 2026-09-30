"""sigma_L, sigma_T and R = L/T against Q^2, with the fit's error band."""
import json, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import fitrun as R, datasets as D, amplitudes as amp

TAG = "C1_with_compass"
p0 = np.load(f"runs/{TAG}/fitpar.npy")
rec = json.load(open(f"runs/{TAG}/summary.json"))
keys, lamd = rec["fitted"], rec["norms"]
nn = [k for k in keys if k in R.NORMS]; nf = len(R.FREE)

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

x0 = np.concatenate([p0[R.FREE], [lamd[k] for k in nn]]); f0 = resid(x0)
J = np.zeros((len(f0), len(x0)))
for i in range(len(x0)):
    h = 1e-6*max(abs(x0[i]), 1e-3); xp = x0.copy(); xp[i] += h
    J[:, i] = (resid(xp) - f0)/h
C = np.linalg.pinv(J.T @ J)[:nf, :nf]
LO, HI = R.bounds()
rng = np.random.default_rng(20260929)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
pars = [R.expand(d) for d in
        np.clip(rng.multivariate_normal(x0[:nf], C, size=N), LO[R.FREE], HI[R.FREE])]

XB, MT = 0.2, 0.3
Q2 = np.linspace(1.0, 6.0, 26)
band = {k: np.zeros((3, len(Q2))) for k in ("L", "T", "R")}
for j, q2 in enumerate(Q2):
    sL, sT, rr = [], [], []
    for q in pars:
        s = amp.structure(q, "pi0p", -MT, XB, q2)
        if s is None: continue
        sL.append(s["L"]); sT.append(s["T"]); rr.append(s["L"]/s["T"])
    for key, a in (("L", sL), ("T", sT), ("R", rr)):
        a = np.array(a)
        band[key][:, j] = (np.median(a), np.percentile(a, 16), np.percentile(a, 84))

C_MOD, C_BAND = "#1a5fb4", "#1a5fb4"
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))
for a, key, lab, log in ((ax[0], "L", r"$d\sigma_L/d|t|$   [nb/(GeV/c)$^2$]", True),
                         (ax[1], "T", r"$d\sigma_T/d|t|$   [nb/(GeV/c)$^2$]", True),
                         (ax[2], "R", r"$R = \sigma_L/\sigma_T$", False)):
    m, lo, hi = band[key]
    a.fill_between(Q2, lo, hi, color=C_BAND, alpha=0.22, lw=0)
    a.plot(Q2, m, "-", color=C_MOD, lw=2)
    if log: a.set_yscale("log")
    else: a.set_ylim(0, None)
    a.set_xlabel(r"$Q^2$   [(GeV/c)$^2$]")
    a.set_ylabel(lab)
    a.grid(alpha=.25, lw=.5)
    a.set_xlim(1, 6)
fig.suptitle(r"VPK amplitude model — $\gamma^* p \to \pi^0 p'$ at "
             rf"$x_B = {XB}$, $|t| = {MT}$ (GeV/c)$^2$" "\n"
             "band = fit uncertainty at fixed functional form "
             f"({N} draws from the parameter covariance, 16–84%)", fontsize=10)
# The run tag is what ties this figure to one parameter vector.  It belongs on
# the figure, small, not in the title where it says nothing to a reader.
fig.text(0.995, 0.012, f"fit {TAG} / amp2609", ha="right", va="bottom",
         fontsize=6.5, color="0.45")
fig.tight_layout(rect=[0, 0.02, 1, 0.90])
fig.savefig("figures/band_LT_xB0.2_t0.3.png", dpi=150)
print("-> figures/band_LT_xB0.2_t0.3.png")
for j, q2 in enumerate(Q2):
    if q2 in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0):
        print(f"  Q2={q2:.0f}  R = {band['R'][0,j]:.3f} [{band['R'][1,j]:.3f}, {band['R'][2,j]:.3f}]")
