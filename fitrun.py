"""Fit the amplitude model to a chosen list of data sets.

    python3 fitrun.py <tag> <set1,set2,...|all>

Writes runs/<tag>/fitpar.npy and runs/<tag>/summary.json, the latter holding the
chi2 of EVERY data set, whether or not it was fitted, so that a run doubles as a
blind prediction for the sets left out.
"""
import json, math, os, sys
import numpy as np
from scipy.optimize import least_squares
import amplitudes as amp, fit_slope as F, datasets as D

W = 30.0
BMAX  = float(os.environ.get("BMAX", "12"))    # ceiling on the slope at xB = 1
BDMAX = os.environ.get("BDMAX")                # separate ceiling for b of H_T^d
BL = {"H_T^u": (1, 2), "H_T^d": (5, 6), "Ebar_T^u": (9, 10), "Ebar_T^d": (14, 12),
      "T00": (16, 21), "Htil": (38, 39)}
BSLOT = [1, 5, 9, 14, 16, 38]
# Slots 37-42 are the Htilde block added 2026-09-30 (see amplitudes._htilde).
# Before it, sigma_L vanished linearly in t' at the forward peak because BOTH
# longitudinal amplitudes carried sqrt(-t').  Slot 25 (rho_nf) is dead now.
# NPAR=37 drops the Htilde slots altogether and puts amplitudes.py back on the
# old code path (T00p = rho_nf*|T00m|).  It exists so that the old model can be
# refitted with exactly the same machinery and effort as the new one: this chi2
# surface has many local minima, and comparing a fresh fit against an older
# published number would not be a fair test of the Htilde block.
NPAR  = int(os.environ.get("NPAR", "43"))
XIPOW = float(os.environ.get("XIPOW", "0"))   # power of xi in front of <Etilde>
if NPAR <= 37:
    del BL["Htil"]; BSLOT.remove(38)
TIES = {7: 3, 27: 11, 6: 2, 12: 10}
# b of H_T^d tied to b of H_T^u.  Left free, the fit runs it to whatever ceiling
# it is given -- 12, 18, 25 -- and buys 2 units of chi2 over 896 points, while
# the normalisation N_d compensates by a factor sixteen.  It is a flat valley,
# not a measurement.
if os.environ.get("TIE_BD", "1") == "1": TIES[5] = 1
# TEST, off by default: the same tie in the Ebar_T sector.  b_d there comes out
# eight times steeper than b_u and is held by four neutron points alone, which
# is what H_T's flat valley looked like before b_d was tied to b_u.
if os.environ.get("TIE_BET", "0") == "1": TIES[14] = 9
# TIE_ELB: b' of <Etilde> tied to b' of <Htilde>.  Htilde and Etilde are the two
# twist-2 longitudinal GPDs and already share one u-d relative phase (p[36]);
# tying their Regge alpha' is a statement of the same kind.  The reason it is
# needed: free, b'(Etilde) came out -0.30 against -1.0 to -1.2 for every other
# block, so past the fitted |t| = 1.8 sigma_L stops falling.  That is invisible
# in the fit -- there is no data there -- but the radiative integral spends 20 %
# of its evaluations beyond |t| = 1.8 and reaches |t| = 14.6, where it made eta
# blow up to 2.6e9 against 36 for amp2609.  See exclurad_py NOTE_amp2610_sampler.
if os.environ.get("TIE_ELB", "0") == "1": TIES[21] = 39
# TIE_EQ: the WHOLE shape of <Etilde> tied to <Htilde> -- both b' (slot 21) and
# nQ (slot 17) -- leaving <Etilde> its own normalisation and d/u ratio only.
#
# Why.  Once Htilde took over the forward region, nothing held Etilde's shape:
# its nQ went +0.394 -> -0.285, i.e. the block started GROWING as Q2 -> 0, and
# its normalisation went 35.8 -> 153.9.  Neither is visible in the fit -- at
# Q2 = 2.5 the Etilde piece of T00p is 1.5 against Htilde's 3.8 -- but the
# radiative integral reaches Q2t = 0.084, where the same piece is 15990 against
# 176, ninety times the Htilde term, and eta came out 1e9 instead of 10.
# The three factors behind it all sit in this block: xi^2/(1-xi^2) x33,
# the normalisation x4.3, the nQ sign x2.3.
#
# NQ_FLOOR additionally forbids a negative nQ: a GFF that grows as Q2 -> 0 is
# not something we want to extrapolate with, whatever chi2 says about it in a
# region where it cannot be measured.
if os.environ.get("TIE_EQ", "0") == "1":
    TIES[21] = 39
    TIES[17] = 40
NQ_FLOOR = float(os.environ.get("NQ_FLOOR", "-99"))
FROZEN = {23, 28, 29, 30, 31, 32, 33} | ({25, 42} if NPAR > 37 else set())

# Sets whose papers quote an overall normalisation uncertainty.  It multiplies
# every point of the set together, so it cannot go into the per-point errors:
# it gets one free scale with a Gaussian prior of the stated width instead.
#   bsa_demasi_phi  3.5%   beam polarisation, De Masi PRC 77 042201(R) p.3
#   halla_y16       3.12%  Defurne PRL 117 262001 table III
#   halla_n         3.1%   Mazouz PRL 118 222002 p.4, the same table
NORMS = {"bsa_demasi_phi": 0.035, "halla_y16": 0.0312, "halla_n": 0.031}
LAM_LO, LAM_HI = 0.7, 1.3
FREE = [i for i in range(NPAR) if i not in set(TIES) | FROZEN]

def slopes(p, x):
    L = math.log(x)
    return {k: (p[b] + p[bp]*L) for k, (b, bp) in BL.items()}

def bounds():
    # Slots 34, 35, 36 are the d-relative phases of H_T, Ebar_T and T00; they enter
    # as exp(i p) and belong on [-pi, pi].  35 and 36 used to carry a magnitude and
    # kept its old window [0.5, 2.0] after the slots were reused, which pinned both
    # phases at 0.5 and also made the computed-phase branch (p[35] <= -9) of
    # amplitudes.py unreachable.
    LO = np.array(list(F.LO) + [-12., -5., -5., -5., -5., -math.pi, -math.pi, -math.pi, -math.pi, -math.pi]
                  + [0., 0., -5., -12., -10., 0.])          # 37-42: the Htilde block
    HI = np.array(list(F.HI) + [12., 5., 8., 5., 8., math.pi, math.pi, math.pi, math.pi, math.pi]
                  + [1e4, 12., 8., 12., 10., 1.])
    LO[3] = LO[11] = LO[17] = -12.; HI[3] = HI[11] = HI[17] = 12.
    LO[12] = -8.; HI[12] = 8.; LO[13] = -1e4; HI[13] = 1e4; LO[14] = -2.; HI[14] = 12.
    if NQ_FLOOR > -90:
        for i in (3, 11, 17, 40):      # the nQ slots
            LO[i] = max(LO[i], NQ_FLOOR)
    for i in BSLOT:
        LO[i] = 0.0
        HI[i] = BMAX          # 12 was inherited from a generic bounds array, not chosen
    if BDMAX is not None: HI[5] = float(BDMAX)
    return LO[:NPAR], HI[:NPAR]

def expand(x):
    p = np.zeros(NPAR); p[FREE] = x; p[23] = 0.0
    if NPAR > 42: p[42] = XIPOW
    for t, s in TIES.items(): p[t] = p[s]
    return p

def fit(keys, seeds=("fitpar_production_pub.npy", "fitpar_n_n_and_p.npy")):
    LO, HI = bounds()
    nn = [k for k in keys if k in NORMS]      # normalised sets actually in this fit
    nf = len(FREE)
    def resid(x):
        p = expand(x[:nf]); lam = dict(zip(nn, x[nf:])); r = []
        for k in keys:
            s = D.BY[k]; L = lam.get(k, 1.0)
            for row in s["rows"]:
                v = s["predict"](p, row)
                r.append((row[4] - L*v)/row[5] if v is not None else 5.0)
        for k in nn: r.append((lam[k] - 1.0)/NORMS[k])
        for xx in (0.05, 0.35, 1.00):
            sl = slopes(p, xx)
            for n in BL: r.append(W*min(0.0, sl[n]))
        return np.array(r)
    best = None
    for src in seeds:
        if not os.path.exists(src): continue
        z = np.load(src)
        for jitter in (0.0, 0.05):
            q = np.zeros(NPAR); q[:min(len(z), NPAR)] = z[:NPAR]
            if NPAR > 41 and q[37] == 0.0:
                # A 37-slot seed has no Htilde.  Start it with the shape of the
                # Etilde block (same slopes, same Q2 power, same d/u ratio) and
                # half its normalisation; starting at exactly zero would leave
                # least_squares with no gradient along the new directions.
                q[37], q[38], q[39], q[40], q[41] = 0.5*q[15], q[16], q[21], q[17], q[18]
            if jitter:
                rng = np.random.default_rng(len(keys))
                q = q*(1 + jitter*rng.standard_normal(NPAR))
            for t, s in TIES.items(): q[t] = q[s]
            q = np.clip(q, LO, HI)
            x0 = np.concatenate([q[FREE], np.ones(len(nn))])
            xlo = np.concatenate([LO[FREE], np.full(len(nn), LAM_LO)])
            xhi = np.concatenate([HI[FREE], np.full(len(nn), LAM_HI)])
            try:
                r = least_squares(resid, x0, bounds=(xlo, xhi),
                                  x_scale='jac', xtol=1e-13, ftol=1e-13, gtol=1e-13,
                                  max_nfev=60000)
            except Exception:
                continue
            c = float(np.sum(r.fun**2))
            if best is None or c < best[0]: best = (c, r.x)
    return expand(best[1][:nf]), dict(zip(nn, (float(v) for v in best[1][nf:])))

def summarise(p, keys, tag, lam=None):
    lam = lam or {}
    used = set(keys)
    npar = len(FREE) + len(lam)
    rec = dict(tag=tag, fitted=sorted(used), npar=npar, norms=lam, sets={})
    cf = nf = 0.0, 0
    cf, nf = 0.0, 0
    for k in D.ALL:
        c, n = D.chi2(p, k, lam.get(k, 1.0))
        rec["sets"][k] = dict(chi2=c, n=n, fitted=k in used, label=D.BY[k]["label"],
                              norm=lam.get(k))
        if k in used: cf += c; nf += n
    cf += sum(((v - 1.0)/NORMS[k])**2 for k, v in lam.items())   # the priors count
    rec["chi2_fitted"] = cf; rec["n_fitted"] = nf
    rec["ndf"] = nf - npar
    rec["chi2_ndf"] = cf/max(nf - npar, 1)
    s25 = slopes(p, 0.25)
    HT, ET, LL = amp._flavour(p, -0.3, 0.25, 2.2)
    sp = amp.structure(p, "pi0p", -0.27, 0.36, 1.75)
    sn = amp.structure(p, "pi0n", -0.27, 0.36, 1.75)
    st = amp.structure(p, "pi0p", -0.4, 0.25, 1.94)
    rec["slopes_xB025"] = {k: float(v) for k, v in s25.items()}
    rec["phases"] = dict(H_T=float(p[34]), Ebar_T=float(p[35]), T00=float(p[36]),
                         Htil=float(p[26]))
    # the forward limit of sigma_L: the whole point of the Htilde block
    fw = [amp.structure(p, "pi0p", amp.tmin(amp.Mpi0, 2.5, 0.2) - d, 0.2, 2.5)
          for d in (0.002, 0.3)]
    rec["sigL_forward"] = dict(tp0002=float(fw[0]["L"]), tp03=float(fw[1]["L"]),
                              R_tp0002=float(fw[0]["L"]/fw[0]["T"]),
                              R_tp03=float(fw[1]["L"]/fw[1]["T"]))
    rec["ratios"] = dict(du_HT=float(abs(HT[1])/HT[0]*np.sign(np.real(HT[1]))),
                         du_ET=float(abs(ET[1])/ET[0]),
                         n_over_p=float(sn["TT"]/sp["TT"]),
                         sigL_over_sigT=float(st["L"]/st["T"]))
    return rec

if __name__ == "__main__":
    tag = sys.argv[1]
    keys = D.ALL if sys.argv[2] == "all" else sys.argv[2].split(",")
    out = f"runs/{tag}"; os.makedirs(out, exist_ok=True)
    p, lam = fit(keys)
    np.save(f"{out}/fitpar.npy", p)
    rec = summarise(p, keys, tag, lam)
    json.dump(rec, open(f"{out}/summary.json", "w"), indent=1)
    print(f"=== {tag}: chi2/ndf = {rec['chi2_ndf']:.4f} "
          f"({rec['chi2_fitted']:.1f}/{rec['ndf']})")
    for k, v in rec["norms"].items():
        print(f"   norm {k:>16} = {v:.4f}  ({(v-1)/NORMS[k]:+.2f} sigma)")
    for k in D.ALL:
        s = rec["sets"][k]
        mark = "fit " if s["fitted"] else "    "
        print(f"   {mark}{k:>12}: {s['chi2']:9.1f}/{s['n']:<5d} {s['chi2']/max(s['n'],1):6.2f}")
