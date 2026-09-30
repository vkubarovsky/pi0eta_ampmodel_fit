"""The bracket B(phi; h, lam) must stay >= 0 over the generator box for every
helicity and target polarisation: it is the phi-differential cross section in
units of sigma_0.  Run it on any fitpar before that fit is installed as a
generator model.   python3 check_positivity.py <fitpar.npy> [<fitpar.npy> ...]"""
import sys, math, numpy as np
import amplitudes as A, polarised as P

# the box VPK specified for the next campaign, 2026-09-18
Q2G  = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]
TPG  = [0.005, 0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0]
EB   = 10.6
CHS  = ["pi0p", "etap", "pi0n", "etan"]
PHI  = np.linspace(0, 2*math.pi, 49)

def scan(p, tag):
    worst = (1e9, None); nbad = 0; ntot = 0
    for ch in CHS:
        mM = A.Meta if ch in ("etap", "etan") else A.Mpi0
        for Q2 in Q2G:
            for xB in np.arange(0.10, 0.65, 0.025):
                W2 = Q2*(1/xB - 1) + A.Mp**2
                if W2 < 1.8**2: continue
                eps = A.epsilon(xB, Q2, EB)
                if not (0.0 < eps < 1.0): continue
                tmin = A.tmin(mM, Q2, xB)
                for d in TPG:
                    s = A.structure(p, ch, tmin - d, xB, Q2)
                    if s is None: continue
                    s0 = s["T"] + eps*s["L"]
                    if s0 <= 0: nbad += 1; continue
                    for h in (-1.0, 0.0, 1.0):
                        for lam in (-1.0, 0.0, 1.0):
                            b = min(P.bracket(s, eps, f, h, lam) for f in PHI)/s0
                            ntot += 1
                            if b < worst[0]:
                                worst = (b, (ch, Q2, round(float(xB),3), d, h, lam))
                            if b < 0: nbad += 1
    print(f"{tag:22} min B/sigma0 = {worst[0]:+.4f}   violations {nbad}/{ntot}")
    print(f"{'':22} worst at ch,Q2,xB,-t',h,lam = {worst[1]}")
    return worst[0], nbad

for f in sys.argv[1:]:
    scan(np.load(f), f.split("/")[-2] if "/" in f else f)
