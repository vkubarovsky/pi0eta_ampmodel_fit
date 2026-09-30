"""chi2 profile in N_Htilde (slot 37): everything else refitted at each fixed
value.  This is the only honest statement about how well the data pin down
sigma_L at the forward peak, where there is no sigma_L measurement at all and
the longitudinal sector enters only through eps*sigma_L in sigma_U and through
the interference terms sigma_LT, sigma_LT'."""
import json, os, sys, numpy as np
import fitrun as R, amplitudes as A

KEYS = ("bsa_clas12,bsa_demasi_phi,bsa_zhao,clas6_eta,clas6_pi0,compass,eg1,"
        "halla_n,halla_y11,halla_y16,halla_y21").split(",")
SEED = sys.argv[1]
GRID = [float(v) for v in sys.argv[2].split(",")]

R.FROZEN = R.FROZEN | {37}
R.FREE = [i for i in range(R.NPAR) if i not in set(R.TIES) | R.FROZEN]
_expand = R.expand
FIXED = [0.0]
def expand(x):
    p = _expand(x); p[37] = FIXED[0]; return p
R.expand = expand

tmin = A.tmin(A.Mpi0, 2.5, 0.2)
out = []
for v in GRID:
    FIXED[0] = v
    p, lam = R.fit(KEYS, seeds=(SEED,))
    rec = R.summarise(p, KEYS, f"scanHt_{v}", lam)
    s0 = A.structure(p, "pi0p", tmin - 0.002, 0.2, 2.5)
    s3 = A.structure(p, "pi0p", tmin - 0.3,   0.2, 2.5)
    row = dict(N=v, chi2=rec["chi2_fitted"], ndf=rec["ndf"],
               sigL0=s0["L"], R0=s0["L"]/s0["T"], R3=s3["L"]/s3["T"])
    out.append(row)
    print(f"N_Htil={v:7.3f}  chi2={row['chi2']:9.2f}  "
          f"sigL(-t'=0.002)={row['sigL0']:8.3f}  R_fwd={row['R0']:7.4f}  R(0.3)={row['R3']:7.4f}",
          flush=True)
json.dump(out, open("runs/scan_htilde.json", "w"), indent=1)
c = min(r["chi2"] for r in out)
print("\nchi2_min =", round(c, 2), " -> Delta chi2 = 1 window:")
for r in out:
    print(f"  N={r['N']:7.3f}  dchi2={r['chi2']-c:8.2f}  R_fwd={r['R0']:7.4f}")
