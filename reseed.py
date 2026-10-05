"""Re-run a fit from a given seed file, to tell a local minimum from a real
difference between models.   python3 reseed.py <tag> <seed.npy> [<seed2.npy>...]
XIPOW must be exported before import, as fitrun reads it at import time."""
import json, os, sys, numpy as np
import fitrun as R, datasets as D
tag, seeds = sys.argv[1], tuple(sys.argv[2:])
# halla_n_U and halla_y16_L are the Rosenbluth data recovered 2026-10-04 (see
# halla_rosenbluth.py): the neutron's sigma_U at two beam energies, which is what
# was measured, and the proton's sigma_L, which never reached the workbook.
# halla_n now carries only sigma_LT and sigma_TT -- its sigma_T is derived from
# the same two sigma_U and would be counted twice.
keys = os.environ.get("FIT_KEYS", "").split(",") if os.environ.get("FIT_KEYS") else (
    "bsa_clas12,bsa_demasi_mom,bsa_zhao,clas6_eta,clas6_pi0,compass,eg1,"
    "halla_n,halla_n_U,halla_y11,halla_y16,halla_y16_L,halla_y21").split(",")
p, lam = R.fit(keys, seeds=seeds)
out = f"runs/{tag}"; os.makedirs(out, exist_ok=True)
np.save(f"{out}/fitpar.npy", p)
rec = R.summarise(p, keys, tag, lam)
rec["xipow"] = R.XIPOW; rec["npar_slots"] = R.NPAR; rec["seeds"] = list(seeds)
json.dump(rec, open(f"{out}/summary.json", "w"), indent=1)
print(tag, "NPAR", R.NPAR, "XIPOW", R.XIPOW, "chi2", round(rec["chi2_fitted"], 1),
      "chi2/ndf", round(rec["chi2_ndf"], 4))
