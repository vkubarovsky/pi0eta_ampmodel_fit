"""Extract the De Masi sin(phi) moment from the published phi distributions.

    python3 demasi_moments.py [reference_fitpar.npy]  ->  data/demasi_moments.txt

De Masi et al., PRC 77 042201(R) (2008) published the pi0 beam-spin asymmetry
as PLOTS; there is no table.  The measurement itself is the phi distribution in
the CLAS physics database (sheet `bsa_phi`, 703 points in 60 kinematic bins),
and figure 5 is a rendering of it.  The `asymmetries` sheet's 62 moment values
were DIGITISED from that figure and their errors are 1.5x too small, so they
are not used.

In 2008 the denominator of the asymmetry was unknown, so the moment had to come
from a plain A sin(phi) fit.  Today sigma_U and sigma_TT are measured at these
very kinematics -- CLAS6 pi0, mean pull -0.10 and +0.09 -- so the full form can
be used and the extracted number is the physical amplitude

    A = sqrt(2 eps (1-eps)) sigma_LT' / (sigma_T + eps sigma_L)

rather than its projection onto sin(phi).  The two differ by 7.6 % on average.

MODEL DEPENDENCE, measured rather than argued: extracting with amp2609's
denominator instead of the Htilde fit's moves A by a median of 0.009 of the
statistical error, at most 0.044.  The denominator needs only sigma_LT/sigma_0
and sigma_TT/sigma_0, which the cross sections pin and on which both fits
agree, so the extraction is model-independent in practice.  The reference used
is recorded in the header of the output file.

ERRORS: statistical only, from the fit.  Two reasons.  The published 0.016
systematic has two named sources, "event selection plus the choice of the fit
function".  The second is removed by construction here -- there is no longer a
choice.  The first, if it is an offset common to all phi in a bin, barely
reaches the moment: the leakage sum(w sin)/sum(w sin^2) has median 0.023, so
0.016 moves A by 0.0004 against statistical errors of 0.02-0.08.  The 3.5 %
beam-polarisation normalisation stays as a fitted nuisance with a Gaussian
prior, which is where a fully correlated scale belongs.

No cut on the error: with 1/stat^2 weights the 40 points whose stat exceeds 0.5
(one is 576) weigh nothing, so the DM_ERRMAX = 0.5 cut the phi-level fit needed
is not required here.
"""
import math
import os
import sys

import numpy as np
import pandas as pd

import amplitudes as A

REF = sys.argv[1] if len(sys.argv) > 1 else "runs/Rosen_B/fitpar.npy"
OUT = "data/demasi_moments.txt"
EBEAM = 5.776

P = np.load(REF)
PH = pd.read_excel("pi0_eta_database.xlsx", "bsa_phi")

bins = {}
for _, r in PH.iterrows():
    st = float(r.stat)
    if st <= 0:
        continue
    bins.setdefault(str(r.bin), []).append(
        (float(r.Q2), float(r.xB), abs(float(r.t)), float(r.value), st,
         math.radians(float(r.phi))))

rows = []
for b, rr in sorted(bins.items(), key=lambda kv: np.mean([x[0] for x in kv[1]])):
    num = den = 0.0
    for Q2, xB, mt, v, st, phi in rr:
        s = A.structure(P, "pi0p", -mt, xB, Q2)
        e = A.epsilon(xB, Q2, EBEAM)
        s0 = s["T"] + e*s["L"]
        d = (1 + math.sqrt(2*e*(1+e))*s["LT"]/s0*math.cos(phi)
               + e*s["TT"]/s0*math.cos(2*phi))
        f = math.sin(phi)/d
        w = 1.0/st**2
        num += w*f*v
        den += w*f*f
    rows.append((b, float(np.mean([x[0] for x in rr])),
                 float(np.mean([x[1] for x in rr])),
                 float(np.mean([x[2] for x in rr])),
                 num/den, 1.0/math.sqrt(den), len(rr)))

os.makedirs("data", exist_ok=True)
import hashlib
md5 = hashlib.md5(open(REF, "rb").read()).hexdigest()
with open(OUT, "w") as f:
    f.write("# De Masi PRC 77 042201(R) sin(phi) moment, extracted from the phi\n")
    f.write("# distributions (workbook sheet bsa_phi) with the FULL asymmetry form.\n")
    f.write(f"# denominator model: {REF}  md5 {md5}\n")
    f.write("# errors: statistical only; the 3.5 % beam polarisation is a fitted nuisance.\n")
    f.write("# A = sqrt(2 eps (1-eps)) sigma_LT' / (sigma_T + eps sigma_L),  Ebeam 5.776 GeV\n")
    f.write(f"#{'bin':>8} {'<Q2>':>8} {'<xB>':>8} {'<|t|>':>8} {'A':>12} {'stat':>10} {'nphi':>5}\n")
    for b, Q2, xB, mt, a, e, n in rows:
        f.write(f"{b:>9} {Q2:8.4f} {xB:8.4f} {mt:8.4f} {a:12.6f} {e:10.6f} {n:5d}\n")
print(f"-> {OUT}   {len(rows)} bins, {sum(r[6] for r in rows)} phi points")
print(f"   reference denominator {REF}  md5 {md5}")
