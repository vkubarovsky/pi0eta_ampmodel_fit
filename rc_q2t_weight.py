"""What FRACTION of sigma_obs comes from Q2t below a threshold?

rc_q2t_exposure.py counted quadrature CALLS, which is not the same thing: the
nodes are placed by the integration scheme, not by the integrand.  Here the
model is zeroed below a Q2t threshold and sigma_obs recomputed; the drop IS the
weight.  vcut = 0 throughout (mx2_cut=None), as in production.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/vpk/exclurad_py_htilde")
from exclurad_py.models import registry
from exclurad_py.core.rc import rc_factor

class Cut:
    """Transparent proxy that zeroes the structure functions below Q2t = CUT."""
    def __init__(self, inner, cut): object.__setattr__(self, "_i", inner); object.__setattr__(self, "_c", cut)
    def __getattr__(self, k): return getattr(self._i, k)
    def __setattr__(self, k, v): setattr(self._i, k, v)
    def __call__(self, W2t, Q2t, tt):
        out = self._i(W2t, Q2t, tt)          # a dict of nine structure functions
        q = np.atleast_1d(np.asarray(Q2t, dtype=float))
        m = q < self._c
        if not m.any(): return out
        return {k: np.where(m, 0.0, np.atleast_1d(v)).reshape(np.shape(v))
                   if np.shape(v) else (0.0 if m.all() else v)
                for k, v in out.items()}

MODEL = sys.argv[1] if len(sys.argv) > 1 else "pi0.amp2614"
CUTS = [0.1, 0.25, 0.5, 0.8, 1.0]
print(f"model {MODEL}, vcut = 0, h = +1, phi = 0.7")
print(f"{'Q2':>5} {'W2':>6} {'-t':>6} {'sigma_obs':>11} | " +
      " ".join(f"{'<'+str(c):>8}" for c in CUTS))
for Q2, W2, t in [(0.9, 5.0, -0.3), (1.0, 6.0, -0.3), (1.5, 8.0, -0.3),
                  (2.0, 6.0, -0.3), (2.0, 12.0, -0.3), (4.0, 9.0, -0.5)]:
    base = registry.get(MODEL, t_nucl=-0.3)
    try:
        r0 = rc_factor(10.6, Q2, W2, t, 0.7, base, h=+1)
    except Exception as e:
        print(f"{Q2:5.1f} {W2:6.1f} {t:6.2f} | FAILED {type(e).__name__}"); continue
    s0 = r0["sigma_obs"]
    row = []
    for c in CUTS:
        try:
            r = rc_factor(10.6, Q2, W2, t, 0.7, Cut(registry.get(MODEL, t_nucl=-0.3), c), h=+1)
            row.append(f"{100*(1 - r['sigma_obs']/s0):7.2f}%")
        except Exception:
            row.append(f"{'--':>8}")
    print(f"{Q2:5.1f} {W2:6.1f} {t:6.2f} {s0:11.4g} | " + " ".join(row))
print("\n  each column: the per cent of sigma_obs contributed by Q2t below that value")

# ---- box average: the number that matters for a campaign ---------------------
if len(sys.argv) > 2 and sys.argv[2] == "box":
    rng = np.random.default_rng(7)
    tot = {c: 0.0 for c in CUTS}; norm = 0.0; n = 0
    print(f"\nbox average over the production box, model {MODEL}")
    for _ in range(int(sys.argv[3]) if len(sys.argv) > 3 else 30):
        Q2 = rng.uniform(0.8, 10.0); W2 = rng.uniform(3.24, 16.0); t = rng.uniform(-2.5, -0.001)
        try:
            s0 = rc_factor(10.6, Q2, W2, t, 0.7, registry.get(MODEL, t_nucl=-0.3), h=+1)["sigma_obs"]
            if not np.isfinite(s0) or s0 <= 0: continue
            for c in CUTS:
                r = rc_factor(10.6, Q2, W2, t, 0.7, Cut(registry.get(MODEL, t_nucl=-0.3), c), h=+1)
                tot[c] += s0 - r["sigma_obs"]
            norm += s0; n += 1
        except Exception:
            continue
    print(f"  {n} usable points")
    for c in CUTS:
        print(f"  Q2t < {c:4.2f}:  {100*tot[c]/norm:6.2f} % of the box-integrated sigma_obs")
