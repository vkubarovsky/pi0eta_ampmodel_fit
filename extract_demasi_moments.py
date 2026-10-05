"""Extract the sin(phi) moment from the De Masi phi distributions, two ways,
and compare both with the published moment values.

The paper (PRC 77 042201(R)) shows the asymmetries only as plots -- there is no
table -- which is why the fit uses the CLAS database's phi distributions
directly.  This script closes the loop: it recovers the moment the paper would
have tabulated, from the same phi points the fit sees, and puts it beside the
`asymmetries` sheet's CLAS6_demasi values and beside the model.

Two extractions, because the convention matters:
    A_sin       data fitted to  A sin(phi)                     (one parameter)
    A_full      data fitted to  A sin(phi) / (1 + b cos + c cos2),
                b and c taken from the MODEL, A the only free parameter
and the model's own
    A_model  =  sqrt(2 eps (1-eps)) sigma_LT' / (sigma_T + eps sigma_L)
"""
import math
import numpy as np
import datasets as D, amplitudes as A

LAM = 1.0630                      # the fitted normalisation nuisance
P = np.load("runs/Rosen_B/fitpar.npy")
S = {s["key"]: s for s in D.SETS}
E = 5.776

def extract(rows, use_den):
    """One-parameter weighted least squares for the sin(phi) amplitude."""
    num = den = 0.0
    for Q2, xB, mt, _, v, e, eps, _, phi in rows:
        s = A.structure(P, "pi0p", -mt, xB, Q2)
        if s is None: return None
        s0 = s["T"] + eps*s["L"]
        f = math.sin(phi)
        if use_den:
            d = (1 + math.sqrt(2*eps*(1+eps))*s["LT"]/s0*math.cos(phi)
                   + eps*s["TT"]/s0*math.cos(2*phi))
            if abs(d) < 1e-6: return None
            f /= d
        w = 1.0/e**2
        num += w*f*v; den += w*f*f
    if den <= 0: return None
    return num/den, 1.0/math.sqrt(den)

BINS = {}
for r in S["bsa_demasi_phi"]["rows"]:
    BINS.setdefault(r[7], []).append(r)

out = []
for b, rows in BINS.items():
    Q2 = np.mean([r[0] for r in rows]); xB = np.mean([r[1] for r in rows])
    mt = np.mean([r[2] for r in rows])
    a1 = extract(rows, False); a2 = extract(rows, True)
    am = A.bsa_sinphi(P, "pi0p", -mt, xB, Q2, E)
    if a1 is None or a2 is None or am is None: continue
    out.append((b, Q2, xB, mt, a1[0], a1[1], a2[0], a2[1], am, len(rows)))

PUB = [(r[0], r[1], r[2], r[4], r[5]) for r in S["bsa_demasi"]["rows"]]

print(f"{'bin':>8} {'Q2':>5} {'xB':>6} {'|t|':>6} {'n':>3} | "
      f"{'A_sin':>14} {'A_full':>14} {'A_model':>8} | {'published':>15} {'dQ2,dxB,dt':>16}")
dd, dm = [], []
for b, Q2, xB, mt, a1, e1, a2, e2, am, n in sorted(out, key=lambda z: (z[1], z[3])):
    j = min(range(len(PUB)),
            key=lambda k: ((PUB[k][0]-Q2)/0.5)**2 + ((PUB[k][1]-xB)/0.05)**2
                          + ((PUB[k][2]-mt)/0.15)**2)
    pq, px, pt, pv, pe = PUB[j]
    close = abs(pq-Q2) < 0.45 and abs(px-xB) < 0.06 and abs(pt-mt) < 0.14
    tag = f"{pv:+7.4f}+-{pe:.4f}" if close else "       no match"
    if close:
        dd.append((a1-pv)/math.hypot(e1, pe)); dm.append((am-pv)/pe)
    print(f"{b:>8} {Q2:5.2f} {xB:6.3f} {mt:6.3f} {n:3d} | "
          f"{a1:+7.4f}+-{e1:.4f} {a2:+7.4f}+-{e2:.4f} {am:+8.4f} | {tag:>15} "
          f"{(f'{pq-Q2:+.2f},{px-xB:+.3f},{pt-mt:+.3f}' if close else ''):>16}")

if dd:
    dd = np.array(dd); dm = np.array(dm)
    print(f"\nmatched bins: {len(dd)}")
    print(f"  extracted vs published : mean pull {dd.mean():+.2f}  rms {dd.std(ddof=1):.2f}  "
          f"chi2/n {np.mean(dd**2):.2f}")
    print(f"  model     vs published : mean pull {dm.mean():+.2f}  rms {dm.std(ddof=1):.2f}  "
          f"chi2/n {np.mean(dm**2):.2f}")
