"""Merge the overlapping N_Htilde chains: keep the lowest chi2 at each node."""
import glob, json, os, re, shutil, numpy as np
best = {}
for d in sorted(glob.glob("runs/scan_Ht_*")):
    for f in sorted(glob.glob(f"{d}/summary_N*.json")):
        N = float(re.search(r"summary_N(-?[\d.]+)\.json", f).group(1))
        c = json.load(open(f))["chi2_fitted"]
        if N not in best or c < best[N][0]: best[N] = (c, f, d)
os.makedirs("runs/scan_Ht", exist_ok=True)
import amplitudes as A
tmin = A.tmin(A.Mpi0, 2.5, 0.2); rows = []
for N in sorted(best):
    c, f, d = best[N]
    p = np.load(f"{d}/fitpar_N{N:g}.npy")
    shutil.copy(f"{d}/fitpar_N{N:g}.npy", f"runs/scan_Ht/fitpar_N{N:g}.npy")
    shutil.copy(f, f"runs/scan_Ht/summary_N{N:g}.json")
    s0 = A.structure(p, "pi0p", tmin - 0.002, 0.2, 2.5)
    s3 = A.structure(p, "pi0p", tmin - 0.3, 0.2, 2.5)
    rows.append(dict(N=N, chi2=c, chain=d.split("_")[-1], sigL0=s0["L"],
                     R0=s0["L"]/s0["T"], R3=s3["L"]/s3["T"]))
json.dump(rows, open("runs/scan_Ht/profile.json", "w"), indent=1)
cmin = min(r["chi2"] for r in rows)
print(f"chi2_min = {cmin:.2f}")
print(f"{'N_Htil':>8} {'dchi2':>8} {'sigL_fwd':>10} {'R_fwd':>8} {'R(0.3)':>8} {'chain':>7}")
for r in rows:
    m = " <-1sig" if r["chi2"]-cmin <= 1 else (" <-2sig" if r["chi2"]-cmin <= 4 else "")
    print(f"{r['N']:8.3f} {r['chi2']-cmin:8.2f} {r['sigL0']:10.3f} {r['R0']:8.4f} "
          f"{r['R3']:8.4f} {r['chain']:>7}{m}")
