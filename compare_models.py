"""Data against BOTH models on the same panels, one multi-page PDF.

    python3 compare_models.py <old_tag> <new_tag> <out.pdf>

Grey dashed = the model without <Htilde> (sigma_L dies at the forward peak),
blue solid = the corrected one.  Each panel quotes both chi2 so the eye and the
number agree.  Normalisation nuisances are applied per model, since the two fits
found different ones.
"""
import json, math, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import amplitudes as amp, datasets as D
from plots import cluster, grid, OBSLAB, PLAN, ORDER

C_DAT, C_OLD, C_NEW = "#8c1d1d", "0.45", "#1a5fb4"

def page(pdf, key, obs, models, recs):
    s = D.BY[key]
    rows = [r for r in s["rows"] if r[3] == obs]
    if not rows: return None
    gs = cluster(rows)
    nr, nc = grid(len(gs))
    W = max(3.3*nc, 9.0); H = max(2.6*nr, 2.6)
    fig, ax = plt.subplots(nr, nc, figsize=(W, H + 0.95), squeeze=False)
    for a in ax.flat: a.set_visible(False)
    tot = [0.0, 0.0]; ntot = 0
    for i, g in enumerate(gs):
        a = ax.flat[i]; a.set_visible(True)
        gl = str(g[0][7])
        xv, axlab = ([r[2] for r in g], r"$-t$  [GeV$^2$]")
        if gl == "compass_Q2": xv, axlab = ([r[0] for r in g], r"$Q^2$  [GeV$^2$]")
        elif key == "bsa_demasi_phi":
            xv, axlab = ([math.degrees(r[8]) for r in g], r"$\phi$  [deg]")
        elif gl == "compass_nu":
            xv, axlab = ([r[0]/(2*0.9382720813*r[1]) for r in g], r"$\nu$  [GeV]")
        t = np.array(xv); v = np.array([r[4] for r in g])
        e = np.array([r[5] for r in g])
        o = np.argsort(t); t, v, e = t[o], v[o], e[o]
        rs = [g[j] for j in o]
        a.errorbar(t, v, e, fmt='o', ms=4, color=C_DAT, ecolor=C_DAT,
                   capsize=2, zorder=4)
        for mi, (p, rec, col, ls, lw) in enumerate(models):
            sc = rec.get("norms", {}).get(key, 1.0) or 1.0
            m = np.array([s["predict"](p, r)*sc if s["predict"](p, r) is not None
                          else np.nan for r in rs], dtype=float)
            good = np.isfinite(m)
            if good.sum() > 1 and key != "compass":
                a.plot(t[good], m[good], ls, lw=lw, color=col, zorder=2+mi)
            a.plot(t[good], m[good], 's', ms=5 if key == "compass" else 3,
                   color=col, zorder=2+mi)
            tot[mi] += float(np.nansum(((v-m)/e)**2))
            if mi == 0: ntot += int(good.sum())
        qs = [r[0] for r in g]; xs = [r[1] for r in g]; ts = [r[2] for r in g]
        qlab = (f"$Q^2$={qs[0]:.2f}" if max(qs)-min(qs) < 0.05*qs[0]
                else f"$Q^2$={min(qs):.2f}-{max(qs):.2f}")
        xlab = (f"$x_B$={xs[0]:.3f}" if max(xs)-min(xs) < 0.05*xs[0]
                else f"$x_B$={min(xs):.3f}-{max(xs):.3f}")
        tlab = (f"\n$-t$={ts[0]:.3f}" if max(ts)-min(ts) < 0.05*max(ts[0], 1e-9)
                else "")
        a.text(0.97, 0.95, f"{qlab}\n{xlab}{tlab}", transform=a.transAxes,
               ha='right', va='top', fontsize=7.5)
        a.grid(alpha=.25, lw=.5); a.tick_params(labelsize=8)
        a.margins(x=0.10, y=0.12)
        a.set_xlabel(axlab, fontsize=9)
        if i % nc == 0: a.set_ylabel(OBSLAB.get(obs, obs), fontsize=9)
    fig.suptitle(f"{s['label']}   —   {OBSLAB.get(obs,obs)}\n"
                 f"no $\\tilde H$: $\\chi^2$={tot[0]:.1f}/{ntot}={tot[0]/max(ntot,1):.2f}"
                 f"      with $\\tilde H$: $\\chi^2$={tot[1]:.1f}/{ntot}"
                 f"={tot[1]/max(ntot,1):.2f}      $\\Delta={tot[1]-tot[0]:+.1f}$",
                 fontsize=11)
    fig.legend(handles=[plt.Line2D([], [], color=C_DAT, marker='o', ls='none', label='data'),
                        plt.Line2D([], [], color=C_OLD, ls='--', label=r'no $\tilde H$ (old)'),
                        plt.Line2D([], [], color=C_NEW, ls='-', label=r'with $\tilde H$')],
               loc='lower center', ncol=3, fontsize=9, frameon=False)
    fig.tight_layout(rect=[0, 0.035, 1, 1 - 0.95/(H + 0.95)])
    pdf.savefig(fig); plt.close(fig)
    return tot[0], tot[1], ntot

if __name__ == "__main__":
    ot, nt, out = sys.argv[1], sys.argv[2], sys.argv[3]
    po = np.load(f"runs/{ot}/fitpar.npy"); ro = json.load(open(f"runs/{ot}/summary.json"))
    pn = np.load(f"runs/{nt}/fitpar.npy"); rn = json.load(open(f"runs/{nt}/summary.json"))
    models = [(po, ro, C_OLD, '--', 1.5), (pn, rn, C_NEW, '-', 1.9)]
    print(f"{'page':44} {'no Htil':>10} {'with Htil':>10} {'n':>5} {'delta':>8}")
    with PdfPages(out) as pdf:
        for key, obs in PLAN:
            for o in obs:
                r = page(pdf, key, o, models, (ro, rn))
                if r: print(f"{key+' '+o:44} {r[0]:10.1f} {r[1]:10.1f} {r[2]:5d} {r[1]-r[0]:+8.1f}")
    print("->", out)
