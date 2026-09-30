"""One continuation chain of the chi2 profile in N_Htilde (slot 37).

    XIPOW=1 python3 scan_htilde.py <seed.npy> <N1,N2,...in walk order> <outdir>

Each node is refitted with slot 37 held fixed and seeded from the PREVIOUS
node's solution as well as from the global seed.  A one-pass scan seeded from a
single point was useless on this 29-parameter surface: nodes fell into different
local minima, R_fwd came out non-monotone in N, and one node beat the free fit
by 35 units.  Launch several overlapping chains in both directions and merge
with merge_scan.py, which keeps the lowest chi2 found at each node.
"""
import json, os, sys, numpy as np
import fitrun as R, amplitudes as A

KEYS = ("bsa_clas12,bsa_demasi_phi,bsa_zhao,clas6_eta,clas6_pi0,compass,eg1,"
        "halla_n,halla_y11,halla_y16,halla_y21").split(",")
SEED, GRID, OUT = sys.argv[1], [float(v) for v in sys.argv[2].split(",")], sys.argv[3]
os.makedirs(OUT, exist_ok=True)

R.FROZEN = R.FROZEN | {37}
R.FREE = [i for i in range(R.NPAR) if i not in set(R.TIES) | R.FROZEN]
_expand, FIXED = R.expand, [0.0]
R.expand = lambda x: (lambda p: (p.__setitem__(37, FIXED[0]), p)[1])(_expand(x))

tmin, prev = A.tmin(A.Mpi0, 2.5, 0.2), SEED
for v in GRID:
    FIXED[0] = v
    p, lam = R.fit(KEYS, seeds=tuple(s for s in (prev, SEED) if os.path.exists(s)))
    rec = R.summarise(p, KEYS, f"N{v:g}", lam)
    np.save(f"{OUT}/fitpar_N{v:g}.npy", p)
    json.dump(rec, open(f"{OUT}/summary_N{v:g}.json", "w"), indent=1)
    s0 = A.structure(p, "pi0p", tmin - 0.002, 0.2, 2.5)
    print(f"N={v:7.3f} chi2={rec['chi2_fitted']:9.2f} R_fwd={s0['L']/s0['T']:7.4f}", flush=True)
    prev = f"{OUT}/fitpar_N{v:g}.npy"
