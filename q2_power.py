"""The Q^2 power of sigma_L, sigma_T and R: twist counting, the two fits, GK, data.

Twist counting at fixed x_B is the same as for vector mesons.  Collinear
factorisation holds for LONGITUDINAL photons (Collins-Frankfurt-Strikman 1997),
so sigma_L is leading twist and dsigma_L/dt ~ Q^-6.  sigma_T is suppressed by one
power of 1/Q in the amplitude, dsigma_T/dt ~ Q^-8, hence R ~ Q^2 -- for vector
mesons R ~ Q^2/M_V^2.

For pi0 the transverse piece comes from the transversity GPDs with the twist-3
meson DA, whose normalisation carries mu_pi = m_pi^2/(m_u+m_d) ~ 2 GeV.  The
suppression is therefore mu_pi^2/Q^2, which is of order one at JLab, so the
asymptotic regime is far away.  That is why sigma_T dominates pi0 while sigma_L
dominates pi+, where the pion pole sits in E~.
"""
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import amplitudes as A

XB, T = 0.2, -0.3
OLD = np.load("runs/old_ctl1/fitpar.npy"); NEW = np.load("runs/Htil_free2/fitpar.npy")
gk = json.loads(open("/Users/vpk/gk_work/compass/results/"
                     f"gk_curve_xB{XB:g}_t{abs(T):g}_GK_2024.json").read())
r = [x for x in gk["rows"] if x["sigmaL"] is not None]
gq = np.array([x["Q2"] for x in r]); gL = np.array([x["sigmaL"] for x in r])
gT = np.array([x["sigmaT"] for x in r])

Q = np.linspace(1.2, 6.0, 60)
def power(v, q):          # smooth local slope from a local quadratic in log-log
    lq = np.log(q); lv = np.log(v)
    return np.gradient(lv, lq)

curves = {}
for nm, p in (("old", OLD), ("new", NEW)):
    for k in ("L", "T"):
        curves[f"{nm}{k}"] = power(np.array([A.structure(p, "pi0p", T, XB, q)[k]
                                             for q in Q]), Q)
    curves[f"{nm}R"] = curves[f"{nm}L"] - curves[f"{nm}T"]
gsm = {k: np.polyfit(np.log(gq[gq >= 2]), np.log(v[gq >= 2]), 1)[0]
       for k, v in (("L", gL), ("T", gT))}
gsm["R"] = gsm["L"] - gsm["T"]

# measured power of sigma_U, single-experiment Q^2 scans
cs = pd.read_excel("pi0_eta_database.xlsx", "cross_sections")
pts = []
for exp, xlo, xhi in (("HallA_y21", 0.44, 0.49),):
    g = cs[(cs.exp == exp) & (cs.xB >= xlo) & (cs.xB <= xhi)].copy()
    g["tpb"] = (g.tprime/0.05).round()*0.05
    for tb, h in g.groupby("tpb"):
        if h.Q2.nunique() < 3: continue
        x = np.log(h.Q2.values); y = np.log(h.s_u.values); w = h.s_u.values/h.stat_U.values
        k, b = np.polyfit(x, y, 1, w=w); res = y - (k*x + b)
        sk = np.sqrt((res**2).sum()/max(len(x)-2, 1))/np.sqrt(((x-x.mean())**2).sum())
        pts.append((np.exp(x.mean()), k, sk, tb))

fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.3))
for a, k, lab, tw in ((ax[0], "L", r"$\sigma_L$", -3.0),
                      (ax[1], "T", r"$\sigma_T$", -4.0),
                      (ax[2], "R", r"$R=\sigma_L/\sigma_T$", +1.0)):
    a.plot(Q, curves[f"old{k}"], "--", color="0.45", lw=1.8, label=r"VPK old (no $\tilde H$)")
    a.plot(Q, curves[f"new{k}"], "-", color="#1a5fb4", lw=2.0, label=r"VPK new (with $\tilde H$)")
    a.axhline(gsm[k], color="#a50f15", ls="-.", lw=1.8, label="GK 2024 (mean, $Q^2$ 2–6)")
    a.axhline(tw, color="#1b7837", ls=":", lw=2.2, label="twist counting")
    if k == "T":
        a.errorbar([p[0] for p in pts], [p[1] for p in pts], [p[2] for p in pts],
                   fmt="s", ms=6, color="k", capsize=3, zorder=5,
                   label=r"Hall A 2021 $\sigma_U$, $x_B\simeq0.46$")
    a.set_xlabel(r"$Q^2$   [GeV$^2$]")
    a.set_ylabel(rf"$d\ln {lab.strip('$')} \,/\, d\ln Q^2$")
    a.grid(alpha=.25); a.set_xlim(1.2, 6)
    a.set_title(lab, fontsize=11)
ax[0].set_ylim(-4.5, -2.0); ax[1].set_ylim(-4.5, -2.0); ax[2].set_ylim(-1.5, 1.5)
ax[1].legend(fontsize=7.5, loc="lower right", framealpha=0.9)
fig.suptitle(r"$Q^2$ power at $x_B=0.2$, $|t|=0.3$ GeV$^2$ — "
             r"twist counting wants $\sigma_L\!\sim\!Q^{-6}$, $\sigma_T\!\sim\!Q^{-8}$, "
             r"$R\!\sim\!Q^{+2}$ (as $R\!\sim\!Q^2/M_V^2$ for vector mesons)",
             fontsize=10.5)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig("figures/q2_power.png", dpi=150)
print("-> figures/q2_power.png")
for k in ("L", "T", "R"):
    print(f"  {k}: old {curves['old'+k][-20]:+.2f}  new {curves['new'+k][-20]:+.2f}"
          f"  GK {gsm[k]:+.2f}  twist {{'L':-3,'T':-4,'R':1}}[k]")
for q, k, sk, tb in pts:
    print(f"  data HallA_y21 -t'={tb:.2f}: {k:+.2f} +- {sk:.2f}")
