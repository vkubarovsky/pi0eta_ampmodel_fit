"""Hall-A pi0 off proton and neutron at Q2 = 1.75, xB = 0.36 -- the kinematic
point of the CERN-talk slide -- with the data and the two model fits.

Proton:  PRL 117 262001 (2016), Rosenbluth separated, sigma_T / sigma_LT / sigma_TT.
Neutron: PRL 118 222002 (2017), sigma_LT / sigma_TT, and sigma_U at the two beam
         energies (the measured quantity; its derived sigma_T is not fitted).

Curves: amp2609 refitted on the recovered Rosenbluth data (Rosen_old_A) and the
Htilde model (Rosen_B).  Both saw exactly these points.
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import amplitudes as A, datasets as D

Q2, XB = 1.75, 0.36
OLD = np.load("runs/Rosen_old_A/fitpar.npy")
NEW = np.load("runs/Rosen_B/fitpar.npy")
MT  = np.linspace(0.12, 0.42, 200)
CP, CN = "#1a5fb4", "#c01c28"          # proton blue, neutron red (as on the slide)

SET = {s["key"]: s for s in D.SETS}

def pick(key, obs, q2=Q2):
    """(|t|, value, error) for one observable at this Q2."""
    r = [(x[2], x[4], x[5]) for x in SET[key]["rows"]
         if x[3] == obs and abs(x[0] - q2) < 1e-6]
    return np.array(sorted(r)).T if r else np.zeros((3, 0))

def curve(p, ch, key, eps=None):
    out = []
    for mt in MT:
        s = A.structure(p, ch, -mt, XB, Q2)
        if s is None: out.append(np.nan)
        elif key == "U": out.append(s["T"] + eps*s["L"])
        else: out.append(s[key])
    return np.array(out)

fig, ax = plt.subplots(1, 4, figsize=(18.0, 4.8))

# ---- 1. sigma_TT, the slide's panel ----------------------------------------
a = ax[0]
for ch, c, key, lab in (("pi0p", CP, "halla_y16", r"$\pi^0/$proton"),
                        ("pi0n", CN, "halla_n",   r"$\pi^0/$neutron")):
    t, v, e = pick(key, "TT")
    a.errorbar(t, v, e, fmt="o" if ch == "pi0p" else "^", color=c, ms=7,
               capsize=3, lw=1.4, label=lab, zorder=5)
    a.plot(MT, curve(OLD, ch, "TT"), "--", color=c, lw=1.6, alpha=.75)
    a.plot(MT, curve(NEW, ch, "TT"), "-",  color=c, lw=2.2)
a.set_ylabel(r"$\sigma_{TT}$   [nb/GeV$^2$]")
a.legend(fontsize=10, frameon=False, loc="lower left")

# ---- 2. sigma_LT ------------------------------------------------------------
a = ax[1]
for ch, c, key in (("pi0p", CP, "halla_y16"), ("pi0n", CN, "halla_n")):
    t, v, e = pick(key, "LT")
    a.errorbar(t, v, e, fmt="o" if ch == "pi0p" else "^", color=c, ms=7,
               capsize=3, lw=1.4, zorder=5)
    a.plot(MT, curve(OLD, ch, "LT"), "--", color=c, lw=1.6, alpha=.75)
    a.plot(MT, curve(NEW, ch, "LT"), "-",  color=c, lw=2.2)
a.set_ylabel(r"$\sigma_{LT}$   [nb/GeV$^2$]")

# ---- 3. sigma_T (p, separated) and sigma_U (n, measured at eps = 0.79) ------
a = ax[2]
t, v, e = pick("halla_y16", "T")
a.errorbar(t, v, e, fmt="o", color=CP, ms=7, capsize=3, lw=1.4,
           label=r"$\sigma_T$  proton", zorder=5)
a.plot(MT, curve(OLD, "pi0p", "T"), "--", color=CP, lw=1.6, alpha=.75)
a.plot(MT, curve(NEW, "pi0p", "T"), "-",  color=CP, lw=2.2)
EPS_N = 0.791
rows = [(x[2], x[4], x[5]) for x in SET["halla_n_U"]["rows"]
        if abs(x[6] - EPS_N) < 1e-3]
t, v, e = np.array(sorted(rows)).T
a.errorbar(t, v, e, fmt="^", color=CN, ms=7, capsize=3, lw=1.4,
           label=r"$\sigma_U$  neutron, $\epsilon = 0.79$", zorder=5)
a.plot(MT, curve(OLD, "pi0n", "U", EPS_N), "--", color=CN, lw=1.6, alpha=.75)
a.plot(MT, curve(NEW, "pi0n", "U", EPS_N), "-",  color=CN, lw=2.2)
a.set_yscale("log")
a.set_ylabel(r"$\sigma_T$ (p),  $\sigma_U$ (n)   [nb/GeV$^2$]")
a.legend(fontsize=9, frameon=False, loc="lower right")

# ---- 4. neutron / proton in sigma_TT ---------------------------------------
a = ax[3]
tp, vp, ep = pick("halla_y16", "TT")
tn, vn, en = pick("halla_n",   "TT")
rt, rv, re = [], [], []
for i in range(len(tp)):                       # the three shared |t| bins
    j = int(np.argmin(abs(tn - tp[i])))
    if abs(tn[j] - tp[i]) > 0.01: continue
    r = vn[j] / vp[i]
    rt.append(tp[i]); rv.append(r)
    re.append(abs(r)*np.hypot(en[j]/vn[j], ep[i]/vp[i]))
a.errorbar(rt, rv, re, fmt="s", color="k", ms=7, capsize=3, lw=1.4,
           label="data", zorder=5)
for p, ls, nm in ((OLD, "--", "amp2609"), (NEW, "-", r"with $\tilde H$")):
    a.plot(MT, curve(p, "pi0n", "TT")/curve(p, "pi0p", "TT"), ls,
           color="#613583", lw=2.0, label=nm)
a.axhline(0.8, color="#26a269", lw=1.6, ls=":", label="GK = 0.8")
a.axhspan(0.21, 0.35, color="0.6", alpha=.30, lw=0, zorder=0)
a.text(0.40, 0.30, r"$0.28\pm0.07$", ha="right", va="center", fontsize=9,
       color="0.25")
a.set_ylim(0, 2.3)
a.set_ylabel(r"$\sigma_{TT}$  neutron / proton")
a.legend(fontsize=9, frameon=False, loc="upper left")

for a in ax:
    a.axhline(0, color="k", lw=0.8)
    a.grid(alpha=.25, lw=.5)
    a.set_xlabel(r"$|t|$   [GeV$^2$]")
    a.set_xlim(0.12, 0.42)
h = [plt.Line2D([], [], color="0.3", ls="--", lw=1.6),
     plt.Line2D([], [], color="0.3", ls="-",  lw=2.2)]
fig.legend(h, ["amp2609", r"with $\tilde H$"], ncol=2, frameon=False,
           fontsize=10, loc="upper right", bbox_to_anchor=(0.995, 0.995))
fig.suptitle(r"Hall-A $\pi^0$ off proton and neutron,  $Q^2 = 1.75$ GeV$^2$, "
             r"$x_B = 0.36$   (PRL 117 262001 / PRL 118 222002)", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig("figures/halla_Q2_1.75.png", dpi=150)
print("-> figures/halla_Q2_1.75.png")

# ---- numbers ----------------------------------------------------------------
print(f"\n{'set':10} {'obs':4} {'|t|':>6} {'data':>10} {'err':>8} "
      f"{'amp2609':>9} {'pull':>6} | {'Htilde':>9} {'pull':>6}")
tot = {"old": 0.0, "new": 0.0}
for key, ch in (("halla_y16", "pi0p"), ("halla_y16_L", "pi0p"),
                ("halla_n", "pi0n"), ("halla_n_U", "pi0n")):
    for row in SET[key]["rows"]:
        if abs(row[0] - Q2) > 1e-6: continue
        pr = SET[key]["predict"]
        po, pn = pr(OLD, row, ch), pr(NEW, row, ch)
        zo, zn = (row[4]-po)/row[5], (row[4]-pn)/row[5]
        tot["old"] += zo*zo; tot["new"] += zn*zn
        print(f"{key:10} {row[3]:4} {row[2]:6.3f} {row[4]:10.2f} {row[5]:8.2f} "
              f"{po:9.2f} {zo:6.2f} | {pn:9.2f} {zn:6.2f}")
w = 1.0/np.array(re)**2
print(f"\nsigma_TT  n/p   data  {np.sum(np.array(rv)*w)/np.sum(w):.2f} "
      f"+- {1/np.sqrt(np.sum(w)):.2f}   (published 0.28 +- 0.07)")
for p, nm in ((OLD, "amp2609"), (NEW, "Htilde ")):
    r = [A.structure(p,"pi0n",-t,XB,Q2)["TT"]/A.structure(p,"pi0p",-t,XB,Q2)["TT"]
         for t in rt]
    print(f"sigma_TT  n/p   {nm}  " + "  ".join(f"{x:.2f}" for x in r)
          + "   at |t| = " + "  ".join(f"{x:.2f}" for x in rt))
print(f"\nchi2 over this Q2 slice:  amp2609 {tot['old']:7.1f}   "
      f"with Htilde {tot['new']:7.1f}")
