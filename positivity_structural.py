"""Is B(phi; h, lam) >= 0 structural, or can the model break it?

The nine structure functions are built from six channel amplitudes through the
Gram dictionary, and the bracket is the phi-differential cross section rebuilt
from them.  If the dictionary and the bracket are mutually consistent, B >= 0
should hold for ANY amplitude set, not just for the fitted one -- positivity
would then be a property of the construction, and the near-tangency found in
pi0 n would be the model saturating it rather than approaching a cliff.

Test: draw random complex amplitudes, build the structure functions exactly as
amplitudes.structure does, and minimise B over phi for every (h, lam).
"""
import math, numpy as np
import polarised as P

S2 = math.sqrt(2)
rng = np.random.default_rng(20261002)
PHI = np.linspace(0, 2*math.pi, 1441)

def sf_from_amps(T00, T01, U01):
    s2 = lambda v: sum(abs(x)**2 for x in v)
    ip = lambda x, y: sum(a.conjugate()*b for a, b in zip(x, y))
    sL, sT01, sU01 = s2(T00), s2(T01), s2(U01)
    uv, uw, vw = -S2*ip(T00, T01), -S2*ip(T00, U01), 2*ip(T01, U01)
    return dict(L=sL, T=sT01+sU01, TT=sT01-sU01, LT=uv.real, LTp=uv.imag,
                LTpLL=uw.real, LTUL=uw.imag, Tp=vw.real, TTUL=vw.imag)

def worst(s, eps):
    s0 = s["T"] + eps*s["L"]
    if s0 <= 0: return None
    return min(min(P.bracket(s, eps, f, h, lam) for f in PHI)/s0
               for h in (-1.0, 0.0, 1.0) for lam in (-1.0, 0.0, 1.0))

def draw(constrained):
    c = lambda: complex(rng.normal(), rng.normal())
    T00p, T00m, T01p, U01p = c(), c(), c(), c()
    if constrained:
        D = c()                      # the model's own degeneracy: T01m = U01m
        T01m = U01m = D/2
    else:
        T01m, U01m = c(), c()
    return (T00p, T00m), (T01p, T01m), (U01p, U01m)

for name, constrained in (("free amplitudes", False),
                          ("model degeneracy T01m = U01m", True)):
    vals = []
    for _ in range(40000):
        s = sf_from_amps(*draw(constrained))
        w = worst(s, float(rng.uniform(0.05, 0.95)))
        if w is not None: vals.append(w)
    v = np.array(vals)
    print(f"{name:32}  n={len(v):6d}  min={v.min():+.3e}  "
          f"below -1e-12: {int((v < -1e-12).sum()):5d}   "
          f"within 1e-6 of zero: {int((abs(v) < 1e-6).sum()):5d}")
