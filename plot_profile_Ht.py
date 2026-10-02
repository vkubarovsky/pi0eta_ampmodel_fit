"""The chi2 profile in N_Htilde and what it allows for R = sigma_L/sigma_T.

There is no sigma_L measurement in the database: the longitudinal sector enters
only through eps*sigma_L inside sigma_U and through sigma_LT, sigma_LT'.  On this
29-parameter surface the refit at a fixed N does not always reach the profile
minimum, so the resolution of the scan is set by OPTIMISER NOISE, not by the
curvature of chi2.  That noise is measured here from the spread between
independent continuation chains that visited the same node, and drawn as a band;
quoting a Delta chi2 = 1 interval below it would be fiction.
"""
import glob, json, re, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = json.load(open("runs/scan_Ht/profile.json"))
N  = np.array([r["N"] for r in rows]); c = np.array([r["chi2"] for r in rows])
R0 = np.array([r["R0"] for r in rows]); R3 = np.array([r["R3"] for r in rows])
d  = c - c.min()

# noise floor: how far apart do independent chains land at the same N?
seen = {}
for f in glob.glob("runs/scan_Ht_*/summary_N*.json"):
    n = float(re.search(r"summary_N(-?[\d.]+)\.json", f).group(1))
    seen.setdefault(n, []).append(json.load(open(f))["chi2_fitted"])
spread = sorted(max(v) - min(v) for v in seen.values() if len(v) > 1)
FLOOR = float(np.median(spread))
print("chain spreads at shared nodes:", [round(s, 1) for s in spread])
print(f"noise floor (median spread) = {FLOOR:.1f} units")

fig, ax = plt.subplots(1, 2, figsize=(10.6, 4.2))
for a in ax:
    a.axhspan(0, FLOOR, color="0.85", zorder=0)
    a.axhline(FLOOR, color="0.5", lw=1, ls="--")
ax[0].plot(N, d, "o-", color="C3", lw=1.6)
ax[0].text(N.max(), FLOOR, " optimiser noise", va="bottom", ha="right",
           fontsize=8, color="0.35")
ax[0].set_xlabel(r"$N_{\tilde H}$"); ax[0].set_ylabel(r"$\Delta\chi^2$")
ax[0].set_ylim(-3, 160); ax[0].grid(alpha=.3)
ax[0].set_title(r"profile in $N_{\tilde H}$, all else refitted", fontsize=10)

o = np.argsort(R0)
ax[1].plot(R0[o], d[o], "o-", color="C0", lw=1.6, label=r"$-t'=0.002$ (forward)")
o3 = np.argsort(R3)
ax[1].plot(R3[o3], d[o3], "s-", color="C2", lw=1.6, label=r"$-t'=0.30$")
ax[1].axvline(0.0009, color="0.3", lw=1.2, ls=":")
ax[1].text(0.0009, 150, " old model", rotation=90, va="top", fontsize=8, color="0.3")
ax[1].set_xlabel(r"$R=\sigma_L/\sigma_T$"); ax[1].set_ylabel(r"$\Delta\chi^2$")
ax[1].set_ylim(-3, 160); ax[1].grid(alpha=.3); ax[1].legend(fontsize=8)
ax[1].set_title(r"$\pi^0p$, $x_B=0.2$, $Q^2=2.5$ GeV$^2$", fontsize=10)
fig.tight_layout(); fig.savefig("figures/htilde_profile.png", dpi=150)
print("-> figures/htilde_profile.png")

ok = d <= FLOOR
print(f"\nwithin the noise floor:  N_Htil {N[ok].min():g}-{N[ok].max():g},"
      f"  R_fwd {R0[ok].min():.3f}-{R0[ok].max():.3f},"
      f"  R(0.3) {R3[ok].min():.3f}-{R3[ok].max():.3f}")
out = d > FLOOR
print("clearly outside:", ", ".join(f"N={n:g}(+{dd:.0f})" for n, dd in zip(N[out], d[out])))
