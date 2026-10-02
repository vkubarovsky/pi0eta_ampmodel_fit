"""R = sigma_L/sigma_T as Q^2 -> 0.  An EXACT constraint, not an asymptotic one.

A real photon has no longitudinal polarisation state, so sigma_L must vanish at
Q^2 = 0.  Current conservation fixes the rate.  With q = (q0, 0, 0, |q|),

    eps_L = (|q|, 0, 0, q0)/sqrt(Q^2),      M_L = eps_L . J
    q.J = 0  =>  J0 = (|q|/q0) Jz
    M_L = (|q| J0 - q0 Jz)/sqrt(Q^2) = Jz (|q|^2 - q0^2)/(q0 sqrt(Q^2))
        = Jz sqrt(Q^2)/q0

so M_L ~ sqrt(Q^2), sigma_L ~ Q^2, and at fixed W and t

    R = sigma_L/sigma_T  ~  Q^2  ->  0 .

In our parametrisation T00 ~ sqrt(A) (Q^2)^{nQ/2} with A ~ Q^-4 at fixed W, so
the limit requires nQ_L -> +3.  The fit gives +0.688 (Htilde) and the old
longitudinal block gave +0.396: both blow up instead of vanishing.  Neither
model was ever asked about this region -- the data start at Q^2 = 1 -- but the
radiative-correction integral evaluates the structure functions at the SHIFTED
Q2t, clamped only at 1e-6 (exclurad_py/core/radiative.py:228).
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import amplitudes as A

W, T = 2.2, -0.3
Mp = A.Mp
xB = lambda Q2: Q2/(Q2 + W*W - Mp*Mp)
OLD = np.load("runs/old_ctl1/fitpar.npy")
NEW = np.load("runs/Htil_free2/fitpar.npy")
Q = np.logspace(-2.2, 0.85, 200)

fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.3))
for p, c, ls, lab in ((OLD, "0.45", "--", r"VPK old (no $\tilde H$)"),
                      (NEW, "#1a5fb4", "-", r"VPK new (with $\tilde H$)")):
    def rr(q):
        s_ = A.structure(p, "pi0p", T, xB(q), q)
        return np.nan if s_ is None else s_["L"]/s_["T"]
    R = np.array([rr(q) for q in Q])
    ax[0].plot(Q, R, ls, color=c, lw=2, label=lab)
    ax[1].plot(Q, np.gradient(np.log(R), np.log(Q)), ls, color=c, lw=2, label=lab)
ref = Q/Q[120]*0.06
ax[0].plot(Q, ref, ":", color="#1b7837", lw=2.2, label=r"required: $R\propto Q^2$")
ax[0].set_yscale("log"); ax[0].set_ylim(1e-3, 30)
ax[0].set_ylabel(r"$R=\sigma_L/\sigma_T$")
ax[1].axhline(1.0, color="#1b7837", ls=":", lw=2.2,
              label=r"required: $d\ln R/d\ln Q^2=+1$")
ax[1].axhline(0.0, color="k", lw=0.8)
ax[1].set_ylabel(r"$d\ln R\,/\,d\ln Q^2$"); ax[1].set_ylim(-2.0, 1.6)
for a in ax:
    a.set_xscale("log"); a.set_xlabel(r"$Q^2$   [GeV$^2$]")
    a.axvspan(1.0, 7.0, color="0.5", alpha=0.12, lw=0)
    a.grid(alpha=.25); a.legend(fontsize=8, loc="best")
    a.set_xlim(Q[0], Q[-1])
ax[0].text(2.0, 2e-3, "fitted\nrange", fontsize=8, color="0.35", ha="center")
fig.suptitle(rf"The photoproduction limit at fixed $W={W}$ GeV, $t={T}$ GeV$^2$ "
             r"— a real photon has no longitudinal state, so $\sigma_L\to0$",
             fontsize=10.5)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("figures/q2_zero_limit.png", dpi=150)
print("-> figures/q2_zero_limit.png")
print(f"{'Q2':>8} {'R old':>10} {'R new':>10}   (required: R -> 0 as Q^2)")
for q in (1.0, 0.3, 0.1, 0.03, 0.01):
    r = []
    for p in (OLD, NEW):
        s_ = A.structure(p, "pi0p", T, xB(q), q)
        r.append(np.nan if s_ is None else s_["L"]/s_["T"])
    print(f"{q:8.3f} {r[0]:10.4f} {r[1]:10.4f}")
