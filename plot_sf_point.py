"""sigma_T, sigma_L, sigma_LT and R = sigma_L/sigma_T versus |t| at one (Q2, xB),
proton and neutron, for the model before and after the Htilde block.

Drawn after the Rosenbluth data lost from the workbook were recovered
(halla_rosenbluth.py), so the two fits here are the ones that saw them:
amp2609 refitted = Rosen_old_A, with Htilde = Rosen_B.
"""
import sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import amplitudes as A

Q2, XB = 2.5, 0.3
OLD = np.load("runs/Rosen_old_A/fitpar.npy")
NEW = np.load("runs/Rosen_B/fitpar.npy")
MT = np.linspace(0.08, 1.6, 300)

COLS = (("T",  r"$\sigma_T$"),
        ("L",  r"$\sigma_L$"),
        ("LT", r"$\sigma_{LT}$"),
        ("R",  r"$R = \sigma_L/\sigma_T$"))

def val(p, ch, mt, key):
    s = A.structure(p, ch, -mt, XB, Q2)
    if s is None: return np.nan
    return s["L"] / s["T"] if key == "R" else s[key]

fig, ax = plt.subplots(2, 4, figsize=(17.0, 7.0), sharex=True)
for row, ch, chlab in ((0, "pi0p", r"$\pi^0 p$"), (1, "pi0n", r"$\pi^0 n$")):
    for col, (key, lab) in enumerate(COLS):
        a = ax[row][col]
        for p, c, ls, name in ((OLD, "0.45", "--", "amp2609"),
                               (NEW, "#1a5fb4", "-", r"with $\tilde H$")):
            a.plot(MT, [val(p, ch, mt, key) for mt in MT], ls,
                   color=c, lw=2, label=name)
        a.axhline(0, color="k", lw=0.8)
        a.grid(alpha=.25, lw=.5)
        if key in ("T", "L"):
            a.set_yscale("symlog", linthresh=1.0)
        if key == "R":
            a.set_ylim(0, None)
            a.set_ylabel(lab)
        else:
            a.set_ylabel(rf"{lab}   [nb/GeV$^2$]" if col == 0 else lab)
        if row == 1: a.set_xlabel(r"$|t|$   [GeV$^2$]")
        if row == 0 and col == 0: a.legend(fontsize=9, frameon=False)
        a.text(0.97, 0.93, chlab, transform=a.transAxes, ha="right", va="top",
               fontsize=13)
fig.suptitle(rf"$Q^2 = {Q2}$ GeV$^2$, $x_B = {XB}$ — the two fits that saw the "
             r"recovered Rosenbluth data", fontsize=11.5)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig("figures/sf_Q2.5_xB0.3.png", dpi=150)
print("-> figures/sf_Q2.5_xB0.3.png")

hdr = (f"\n{'ch':6} {'|t|':>6} | {'sigT old':>9} {'sigT new':>9} | {'sigL old':>9} "
       f"{'sigL new':>9} | {'sigLT old':>9} {'sigLT new':>9} | {'R old':>7} {'R new':>7}")
print(hdr)
for ch in ("pi0p", "pi0n"):
    for mt in (0.1, 0.2, 0.4, 0.8, 1.5):
        a_, b_ = A.structure(OLD, ch, -mt, XB, Q2), A.structure(NEW, ch, -mt, XB, Q2)
        if a_ is None or b_ is None: continue
        print(f"{ch:6} {mt:6.2f} | {a_['T']:9.2f} {b_['T']:9.2f} | {a_['L']:9.2f} "
              f"{b_['L']:9.2f} | {a_['LT']:9.2f} {b_['LT']:9.2f} | "
              f"{a_['L']/a_['T']:7.3f} {b_['L']/b_['T']:7.3f}")
