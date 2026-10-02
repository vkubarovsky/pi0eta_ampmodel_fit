"""Re-run a fit from a given seed file, to tell a local minimum from a real
difference between models.   python3 reseed.py <tag> <seed.npy> [<seed2.npy>...]
XIPOW must be exported before import, as fitrun reads it at import time."""
import json, os, sys, numpy as np
import fitrun as R, datasets as D
tag, seeds = sys.argv[1], tuple(sys.argv[2:])
keys = ("bsa_clas12,bsa_demasi_phi,bsa_zhao,clas6_eta,clas6_pi0,compass,eg1,"
        "halla_n,halla_y11,halla_y16,halla_y21").split(",")
p, lam = R.fit(keys, seeds=seeds)
out = f"runs/{tag}"; os.makedirs(out, exist_ok=True)
np.save(f"{out}/fitpar.npy", p)
rec = R.summarise(p, keys, tag, lam)
rec["xipow"] = R.XIPOW; rec["npar_slots"] = R.NPAR; rec["seeds"] = list(seeds)
json.dump(rec, open(f"{out}/summary.json", "w"), indent=1)
print(tag, "NPAR", R.NPAR, "XIPOW", R.XIPOW, "chi2", round(rec["chi2_fitted"], 1),
      "chi2/ndf", round(rec["chi2_ndf"], 4))
