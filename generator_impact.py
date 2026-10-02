"""What changes for the GENERATOR if the Htilde fit replaces the old one.

The generator samples sigma_U = sigma_T + eps sigma_L for the rate and the phi
shape from sigma_LT, sigma_TT.  chi2 is not the quantity that matters here: the
cross section over the campaign box is.  Prints the ratio new/old of the rate and
of the two phi modulations, and the change in the t' spectrum, which is where
the forward behaviour of sigma_L actually shows up.
"""
import math, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import amplitudes as A

OLD = np.load("runs/old_ctl1/fitpar.npy")
NEW = np.load("runs/Htil_free2/fitpar.npy")
EB  = 10.6
Q2G = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]
XBG = [0.15, 0.20, 0.25, 0.30, 0.35, 0.45, 0.55]
TPG = [0.005, 0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0]

def obs(p, ch, Q2, xB, d):
    mM = A.Meta if ch in ("etap", "etan") else A.Mpi0
    W2 = Q2*(1/xB - 1) + A.Mp**2
    if W2 < 1.8**2: return None
    eps = A.epsilon(xB, Q2, EB)
    if not (0.0 < eps < 1.0): return None
    s = A.structure(p, ch, A.tmin(mM, Q2, xB) - d, xB, Q2)
    if s is None: return None
    s0 = s["T"] + eps*s["L"]
    if s0 <= 0: return None
    return dict(s0=s0, cos=math.sqrt(2*eps*(1+eps))*s["LT"]/s0,
                cos2=eps*s["TT"]/s0, L=s["L"], T=s["T"])

for ch in ("pi0p", "etap", "pi0n"):
    r = []; dc = []; dc2 = []
    for Q2 in Q2G:
        for xB in XBG:
            for d in TPG:
                a, b = obs(OLD, ch, Q2, xB, d), obs(NEW, ch, Q2, xB, d)
                if a is None or b is None: continue
                r.append(b["s0"]/a["s0"]); dc.append(b["cos"]-a["cos"])
                dc2.append(b["cos2"]-a["cos2"])
    r = np.array(r); dc = np.array(dc); dc2 = np.array(dc2)
    q = lambda v, f: float(np.quantile(v, f))
    print(f"{ch}:  n={len(r)}")
    print(f"   rate sigma_U  new/old : median {np.median(r):.3f}   "
          f"5-95% [{q(r,.05):.3f}, {q(r,.95):.3f}]   min {r.min():.3f} max {r.max():.3f}")
    print(f"   cos phi  coefficient  : shift median {np.median(dc):+.4f}   "
          f"5-95% [{q(dc,.05):+.4f}, {q(dc,.95):+.4f}]")
    print(f"   cos2phi  coefficient  : shift median {np.median(dc2):+.4f}   "
          f"5-95% [{q(dc2,.05):+.4f}, {q(dc2,.95):+.4f}]")

# the t' spectrum: the shape the generator actually samples
fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.0))
tp = np.logspace(math.log10(0.003), math.log10(2.0), 120)
for j, (Q2, xB) in enumerate([(1.5, 0.15), (2.5, 0.25), (4.0, 0.40)]):
    a = ax[j]
    for p, c, ls, lab in ((OLD, "0.45", "--", r"no $\tilde H$"),
                          (NEW, "#1a5fb4", "-", r"with $\tilde H$")):
        y = [obs(p, "pi0p", Q2, xB, d) for d in tp]
        a.plot(tp, [v["s0"] if v else np.nan for v in y], ls, color=c, lw=1.8, label=lab)
    a.set_xscale("log"); a.set_yscale("log"); a.grid(alpha=.3)
    a.set_xlabel(r"$-t'$  [GeV$^2$]")
    if j == 0: a.set_ylabel(r"$\sigma_U=\sigma_T+\epsilon\sigma_L$  [nb/GeV$^2$]")
    a.set_title(rf"$Q^2$={Q2}, $x_B$={xB}", fontsize=10)
    a.legend(fontsize=8)
    r = [(obs(NEW,"pi0p",Q2,xB,d) or {}).get("s0", np.nan) /
         ((obs(OLD,"pi0p",Q2,xB,d) or {}).get("s0", np.nan)) for d in tp]
    tw = a.twinx(); tw.plot(tp, r, ":", color="#c2621b", lw=1.4)
    tw.set_ylim(0.6, 1.8); tw.tick_params(labelsize=8, colors="#c2621b")
    if j == 2: tw.set_ylabel("new / old", color="#c2621b", fontsize=9)
fig.suptitle(r"What the generator samples: $\sigma_U$ vs $-t'$ for $\pi^0p$"
             "   (dotted orange, right axis: the ratio)", fontsize=11)
fig.tight_layout(); fig.savefig("figures/generator_impact.png", dpi=150)
print("\n-> figures/generator_impact.png")
