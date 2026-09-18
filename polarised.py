"""The polarised-target cross section, and reference numbers for the generator.

The generator currently samples five structure functions and therefore knows
only the beam-spin asymmetry.  The model has had all nine since the beginning,
and two of the four polarised ones are fitted -- the eg1-dvcs A_UL and A_LL
moments are in the fit.  This module writes the full bracket down once, in the
same convention the fit uses, and produces test vectors so that an independent
implementation in exclurad_py/core/born.py can be checked against it rather than
believed.

Convention.  For beam helicity h = +-1 and LONGITUDINAL target polarisation
Lam = +-1, the cross section is Gamma/(2 pi) times

    B(phi; h, Lam) =
        sigma_T + eps sigma_L
      + eps cos(2 phi)               sigma_TT
      + sqrt(2 eps (1+eps)) cos(phi) sigma_LT
      + h   sqrt(2 eps (1-eps)) sin(phi) sigma_LT'
      + Lam [ sqrt(2 eps (1+eps)) sin(phi) sigma_LT^UL
              + eps sin(2 phi)             sigma_TT^UL ]
      + h Lam [ sqrt(1 - eps^2)            sigma_T'
                + sqrt(2 eps (1-eps)) cos(phi) sigma_LT'^LL ]

The first line and a half is what born.py already computes.  The last two lines
are the new terms.  Dividing each by sigma_0 = sigma_T + eps sigma_L reproduces
the four eg1 moments exactly as datasets.py predicts them, which is the check
that this bracket and the fit agree:

    A_UL^sin phi   = sqrt(2 eps (1+eps)) sigma_LT^UL   / sigma_0
    A_UL^sin 2phi  = eps                 sigma_TT^UL   / sigma_0
    A_LL^const     = sqrt(1 - eps^2)     sigma_T'      / sigma_0
    A_LL^cos phi   = sqrt(2 eps (1-eps)) sigma_LT'^LL  / sigma_0

A transversely polarised target needs structure functions this model does not
have; only the longitudinal case is covered here.
"""
import math

import numpy as np

import amplitudes as amp

KEYS = ("T", "L", "TT", "LT", "LTp", "LTUL", "TTUL", "Tp", "LTpLL")


def bracket(s, eps, phi, h=0.0, lam=0.0):
    """B(phi; h, lam) from the nine structure functions of amp.structure."""
    a = math.sqrt(2*eps*(1 + eps))
    b = math.sqrt(2*eps*(1 - eps))
    out = (s["T"] + eps*s["L"]
           + eps*math.cos(2*phi)*s["TT"]
           + a*math.cos(phi)*s["LT"]
           + h*b*math.sin(phi)*s["LTp"])
    if lam:
        out += lam*(a*math.sin(phi)*s["LTUL"] + eps*math.sin(2*phi)*s["TTUL"])
        if h:
            out += h*lam*(math.sqrt(1 - eps*eps)*s["Tp"]
                          + b*math.cos(phi)*s["LTpLL"])
    return out


def moments(s, eps):
    """The four eg1 moments, for checking against datasets.py."""
    s0 = s["T"] + eps*s["L"]
    return dict(AULsin=math.sqrt(2*eps*(1 + eps))*s["LTUL"]/s0,
                AULsin2=eps*s["TTUL"]/s0,
                ALLc=math.sqrt(1 - eps*eps)*s["Tp"]/s0,
                ALLcos=math.sqrt(2*eps*(1 - eps))*s["LTpLL"]/s0)


def check_against_fit(p, tol=1e-12):
    """The bracket must reproduce the moments datasets.py fits.  Loudly."""
    import datasets as D
    S = D.BY["eg1"]
    worst = 0.0
    for r in S["rows"]:
        Q2, xB, mt, key = r[0], r[1], r[2], r[3]
        s = amp.structure(p, "pi0p", -mt, xB, Q2)
        mine = moments(s, r[6])[key]
        worst = max(worst, abs(mine - S["predict"](p, r)))
    if worst > tol:
        raise AssertionError(f"bracket and fit disagree by {worst:.3g}")
    return worst


def vectors(p, out="data/polarised_reference.data"):
    """Test vectors for an independent implementation to reproduce."""
    E = 10.6
    grid = [(1.75, 0.36, 0.27), (2.50, 0.25, 0.50), (3.50, 0.40, 0.80),
            (2.00, 0.20, 0.30), (5.00, 0.50, 1.20)]
    rows = []
    for Q2, xB, mt in grid:
        eps = amp.epsilon(xB, Q2, E)
        s = amp.structure(p, "pi0p", -mt, xB, Q2)
        if s is None:
            continue
        for deg in range(0, 360, 30):
            phi = math.radians(deg)
            rows.append([Q2, xB, mt, eps, deg]
                        + [s[k] for k in KEYS]
                        + [bracket(s, eps, phi, h, lam)
                           for h, lam in ((0, 0), (1, 0), (0, 1), (1, 1), (1, -1))])
    hdr = ("Reference values for the polarised-target bracket, model amp2609 =\n"
           "rc_iter1 run C1_with_compass, proton target, pi0, Ebeam 10.6 GeV.\n"
           "Structure functions in nb/GeV^2; the bracket in the same units.\n"
           "columns:\n"
           "  Q2 xB -t eps phi[deg] " + " ".join(KEYS) + "\n"
           "  B(h=0,L=0) B(h=+1,L=0) B(h=0,L=+1) B(h=+1,L=+1) B(h=+1,L=-1)\n"
           "B(0,0) is what the generator computes today; the rest is what it\n"
           "should compute once the polarised terms are in.  See polarised.py\n"
           "for the formula and for the identity that ties it to the fit.")
    np.savetxt(out, np.array(rows), fmt="%13.6g", header=hdr)
    return len(rows)


if __name__ == "__main__":
    p = np.load("runs/C1_with_compass/fitpar.npy")
    w = check_against_fit(p)
    print(f"bracket == fit moments to {w:.2g}")
    n = vectors(p)
    print(f"{n} reference rows -> data/polarised_reference.data")
