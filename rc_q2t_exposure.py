"""How far below the fitted range does the RC integral actually reach?

Wraps the structure-function model so every call records the shifted Q2t.  The
only question that matters for whether the Q^2 -> 0 behaviour of the amplitude
model is a practical concern: does the radiative integral sample Q2t below the
fitted Q^2 = 1, and with what weight?
"""
import sys, numpy as np
sys.path.insert(0, "/Users/vpk/exclurad_py")
from exclurad_py.models import registry
from exclurad_py.core.rc import rc_factor

SEEN = []
class Probe:
    """Transparent proxy that records Q2t on every structure-function call."""
    def __init__(self, inner): object.__setattr__(self, "_inner", inner)
    def __getattr__(self, k): return getattr(self._inner, k)
    def __setattr__(self, k, v): setattr(self._inner, k, v)
    def __call__(self, W2t, Q2t, tt):
        SEEN.append(np.atleast_1d(np.asarray(Q2t, dtype=float)).ravel().copy())
        return self._inner(W2t, Q2t, tt)

inner = registry.get("pi0.amp2609", t_nucl=-0.3)
sf = Probe(inner)
Mp = 0.938272
print(f"{'Q2':>5} {'W2':>6} {'-t':>5} | {'n calls':>9} {'Q2t min':>9} "
      f"{'frac<1':>8} {'frac<0.5':>9} {'frac<0.1':>9}")
for Q2, W2, t in [(1.0, 6.0, -0.3), (2.0, 6.0, -0.3), (2.0, 9.0, -0.3),
                  (4.0, 9.0, -0.5), (6.0, 12.0, -0.8)]:
    SEEN.clear()
    try:
        r = rc_factor(10.6, Q2, W2, t, 0.7, sf, h=+1, lam=0.0)
    except Exception as e:
        print(f"{Q2:5.1f} {W2:6.1f} {t:5.2f} | FAILED {type(e).__name__}: {e}")
        continue
    q = np.concatenate(SEEN) if SEEN else np.array([])
    if not len(q):
        print(f"{Q2:5.1f} {W2:6.1f} {t:5.2f} | no calls recorded"); continue
    print(f"{Q2:5.1f} {W2:6.1f} {t:5.2f} | {len(q):9d} {q.min():9.4f} "
          f"{(q < 1.0).mean():8.3f} {(q < 0.5).mean():9.3f} {(q < 0.1).mean():9.3f}"
          f"   eta={r.get('eta', float('nan')):.4f}")
