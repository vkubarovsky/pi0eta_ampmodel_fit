"""exp[(b + b' ln xB) t] for every GFF block, xB = 0.1 - 0.5.

The t-dependence of each block, stripped of its normalisation, so the only thing
on show is how fast it falls.  <Etilde> in amp2610 is the outlier: b' = -0.30
against -1.0 to -1.3 everywhere else, so it flattens instead of steepening
towards small xB, and past the fitted |t| = 1.8 it stops falling at all.
"""
import math, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NEW = np.load("runs/Htil_free2/fitpar.npy")
OLD = np.load("runs/old_ctl1/fitpar.npy")
BLOCKS = [
    (r"$\langle\tilde E\rangle$  amp2610", NEW[16], NEW[21], "#c0392b", "-",  2.6),
    (r"L block  amp2609 (was $\tilde E$'s role)", OLD[16], OLD[21], "0.35", "--", 1.8),
    (r"$\langle\tilde H\rangle$  amp2610", NEW[38], NEW[39], "#1a5fb4", "-", 1.8),
    (r"$\langle H_T\rangle$  amp2610",     NEW[1],  NEW[2],  "#1b7837", "-.", 1.5),
    (r"$\langle\bar E_T\rangle$  amp2610", NEW[9],  NEW[10], "#8e44ad", ":",  1.8),
]
XB = [0.1, 0.2, 0.3, 0.4, 0.5]
MT = np.linspace(0.0, 3.0, 400)

fig, ax = plt.subplots(1, 5, figsize=(16.0, 3.9), sharey=True)
for a, x in zip(ax, XB):
    L = math.log(x)
    a.axvspan(1.8, 3.0, color="0.85", zorder=0)
    a.axvline(2.5, color="0.55", lw=1, ls=":")
    for lab, b, bp, col, ls, lw in BLOCKS:
        s = b + bp*L
        a.plot(MT, np.exp(-s*MT), ls, color=col, lw=lw,
               label=f"{lab}   $b+b'\\ln x_B$ = {s:.2f}" if x == XB[0] else None)
    a.set_yscale("log"); a.set_ylim(1e-4, 1.5)
    a.set_xlim(0, 3.0); a.grid(alpha=.25, lw=.5)
    a.set_xlabel(r"$|t|$   [GeV$^2$]")
    a.set_title(rf"$x_B = {x}$", fontsize=11)
ax[0].set_ylabel(r"$\exp[(b + b'\ln x_B)\,t]$")
ax[2].text(2.4, 2.5e-4, "outside\nfit validity", fontsize=7.5, color="0.35",
           ha="center")
h, l = ax[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=5, fontsize=8.0, frameon=False,
           bbox_to_anchor=(0.5, -0.03))
fig.suptitle(r"t-dependence of each GFF block, normalisation removed — "
             r"grey band is beyond the fitted $|t| = 1.8$, dotted line is the "
             r"RG-A card edge $|t| = 2.5$", fontsize=10.5)
fig.tight_layout(rect=[0, 0.09, 1, 0.93])
fig.savefig("figures/tslopes_xB.png", dpi=150, bbox_inches="tight")
print("-> figures/tslopes_xB.png")
print()
print(f"{'xB':>5} " + "".join(f"{n.split()[0][-12:]:>14}" for n, *_ in BLOCKS))
print(f"{'':5} " + "".join(f"{'slope':>14}" for _ in BLOCKS))
for x in XB:
    L = math.log(x)
    print(f"{x:5.2f} " + "".join(f"{b+bp*L:14.3f}" for _, b, bp, *_ in BLOCKS))
print()
print("suppression at |t| = 2.5 relative to |t| = 0:")
print(f"{'xB':>5} " + "".join(f"{n.split()[0][-12:]:>14}" for n, *_ in BLOCKS))
for x in XB:
    L = math.log(x)
    print(f"{x:5.2f} " + "".join(f"{math.exp(-(b+bp*L)*2.5):14.2e}"
                                 for _, b, bp, *_ in BLOCKS))
