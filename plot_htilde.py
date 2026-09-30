"""sigma_L, sigma_T and R = L/T versus -t': the old longitudinal sector against
the corrected one.  The old model had no <Htilde>: T00p = rho_nf*|T00m| carried
the sqrt(-t') of the FLIP amplitude, so sigma_L died linearly at the forward
peak.  Angular momentum says the non-flip amplitude survives t' -> 0."""
import sys, json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import amplitudes as A

XB, Q2, CH = 0.2, 2.5, "pi0p"
RUNS = [("runs/old_ctl1/fitpar.npy", "no $\\tilde H$ (old)", "0.45", "--"),
        ("runs/Htil_xi0_rs/fitpar.npy",        "with $\\tilde H$, $n_x=0$",     "C3", "-"),
        ("runs/Htil_free2/fitpar.npy",        "with $\\tilde H$, explicit $\\xi$ (converged)", "C0", "-.")]

tmin = A.tmin(A.Mpi0, Q2, XB)
tp = np.concatenate([np.linspace(0.0005, 0.05, 60), np.linspace(0.05, 1.2, 140)])
fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.0))
for path, lab, col, ls in RUNS:
    try: p = np.load(path)
    except FileNotFoundError: print("missing", path); continue
    L, T = [], []
    for d in tp:
        s = A.structure(p, CH, tmin - d, XB, Q2)
        L.append(s["L"]); T.append(s["T"])
    L, T = np.array(L), np.array(T)
    ax[0].plot(tp, L, ls, color=col, label=lab, lw=1.8)
    ax[1].plot(tp, T, ls, color=col, label=lab, lw=1.8)
    ax[2].plot(tp, L/T, ls, color=col, label=lab, lw=1.8)
for a, yl in zip(ax, [r"$\sigma_L$  (nb/GeV$^2$)", r"$\sigma_T$  (nb/GeV$^2$)",
                      r"$R=\sigma_L/\sigma_T$"]):
    a.set_xlabel(r"$-t'$  (GeV$^2$)"); a.set_ylabel(yl)
    a.set_xscale("log"); a.grid(alpha=.3); a.legend(fontsize=8)
ax[0].set_yscale("log")
fig.suptitle(rf"$\pi^0 p$,  $x_B={XB}$,  $Q^2={Q2}$ GeV$^2$ — the forward limit of $\sigma_L$",
             fontsize=11)
fig.tight_layout()
fig.savefig("figures/htilde_forward.png", dpi=150)
print("-> figures/htilde_forward.png")
for path, lab, _, _ in RUNS:
    try: p = np.load(path)
    except FileNotFoundError: continue
    print(f"\n{lab}")
    for d in (0.002, 0.01, 0.05, 0.3, 1.0):
        s = A.structure(p, CH, tmin - d, XB, Q2)
        print(f"   -t'={d:5.3f}  sigL={s['L']:9.4g}  R={s['L']/s['T']:8.5f}")
