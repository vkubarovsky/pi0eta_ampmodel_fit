"""Does the bracket actually go negative anywhere in the campaign box?

check_positivity.py used a coarse grid and reported min B/sigma0 = +0.0009 for
the H~ fit against +0.0055 for the old one, both in pi0 n at full polarisation.
A margin that small on a coarse grid is not evidence that the model is safe: the
grid may simply have missed the minimum.  This refines it.

B(phi; h, lam) is the phi-differential cross section in units of sigma_0.  If it
goes negative the generator has no valid distribution to sample: rejection
sampling would either lose events or produce negative weights.

The phi minimum is found on a fine grid and then refined by golden section; the
kinematic minimum by a coarse pass followed by Nelder-Mead from the best cells.
"""
import sys, math, numpy as np
from scipy.optimize import minimize_scalar, minimize
import amplitudes as A, polarised as P

EB = 10.6
WMIN2 = 1.8**2
CH = sys.argv[1] if len(sys.argv) > 1 else "pi0n"
TAGS = [("old_ctl1", "no H~"), ("Htil_free2", "with H~")]

def bracket_min(p, ch, Q2, xB, dtp, h, lam):
    """min over phi of B/sigma0, or None if the point is not physical."""
    mM = A.Meta if ch in ("etap", "etan") else A.Mpi0
    if Q2*(1/xB - 1) + A.Mp**2 < WMIN2: return None
    eps = A.epsilon(xB, Q2, EB)
    if not (0.0 < eps < 1.0): return None
    s = A.structure(p, ch, A.tmin(mM, Q2, xB) - dtp, xB, Q2)
    if s is None: return None
    s0 = s["T"] + eps*s["L"]
    if s0 <= 0: return None
    f = lambda ph: P.bracket(s, eps, ph, h, lam)/s0
    grid = np.linspace(0, 2*math.pi, 721)
    v = np.array([f(x) for x in grid])
    j = int(v.argmin())
    lo, hi = grid[max(j-1, 0)], grid[min(j+1, len(grid)-1)]
    r = minimize_scalar(f, bounds=(lo, hi), method="bounded",
                        options={"xatol": 1e-10})
    return min(float(r.fun), float(v[j]))

def hunt(p, ch):
    best = (1e9, None)
    # coarse pass
    cells = []
    for Q2 in np.linspace(1.0, 6.0, 21):
        for xB in np.arange(0.10, 0.651, 0.0125):
            for dtp in [0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4,
                        0.6, 0.8, 1.0, 1.3, 1.6, 2.0]:
                for h in (-1.0, 1.0):
                    for lam in (-1.0, 1.0):
                        b = bracket_min(p, ch, Q2, xB, dtp, h, lam)
                        if b is None: continue
                        cells.append((b, Q2, xB, dtp, h, lam))
                        if b < best[0]: best = (b, (Q2, xB, dtp, h, lam))
    cells.sort()
    # local refinement from the ten worst cells
    for b0, Q2, xB, dtp, h, lam in cells[:10]:
        def neg(z):
            q, x, d = z
            if not (1.0 <= q <= 6.0 and 0.08 <= x <= 0.70 and 0.002 <= d <= 2.0):
                return 10.0
            v = bracket_min(p, ch, q, x, d, h, lam)
            return 10.0 if v is None else v
        r = minimize(neg, [Q2, xB, dtp], method="Nelder-Mead",
                     options={"xatol": 1e-6, "fatol": 1e-10, "maxiter": 600})
        if r.fun < best[0]: best = (float(r.fun), (r.x[0], r.x[1], r.x[2], h, lam))
    return best, len(cells), sum(1 for c in cells if c[0] < 0)

print(f"channel: {CH}   box Q2 1-6, W >= 1.8, -t' 0.005-2.0, h and lam = +-1")
for tag, lab in TAGS:
    p = np.load(f"runs/{tag}/fitpar.npy")
    (b, at), n, nneg = hunt(p, CH)
    Q2, xB, dtp, h, lam = at
    flag = "   *** NEGATIVE ***" if b < 0 else ""
    print(f"  {lab:9} min B/sigma0 = {b:+.6f}   at Q2={Q2:.3f} xB={xB:.4f} "
          f"-t'={dtp:.4f} h={h:+.0f} lam={lam:+.0f}   "
          f"(grid cells {n}, negative on grid {nneg}){flag}")
