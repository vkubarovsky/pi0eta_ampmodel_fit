"""Every data set, loaded from the workbook, with its model prediction.

One place decides what each measurement means -- sigma_T against sigma_U, a
ratio against an asymmetry moment, which flavour combination, which beam energy
-- so that a fit and a plot can never disagree about it.

Each entry of SETS is a dict with
    rows   list of (Q2, xB, mt, obs, value, err, eps)
    label  what to print
    kind   "xs" or "asym"
and a predict(p, row) that returns the model number for that row.
"""
import json, math, os
import numpy as np
import pandas as pd
import amplitudes as amp

# Everything the fit reads comes from the workbook, exported sheet by sheet.
# Nothing is loaded from anywhere else: a second source is how the old De Masi
# amplitudes, with their errors 1.5 times too small, survived unnoticed.
XS = pd.read_csv("db_csv/cross_sections.csv")
AS = pd.read_csv("db_csv/asymmetries.csv")
PH = pd.read_csv("db_csv/bsa_phi.csv")
# Last-resort floor for a row that carries no uncertainty at all.  It is NOT a
# systematic: every set now brings its own from its paper.  Hall-A 2016 quotes a
# 2% point-to-point already inside its published errors plus a 3.12% overall
# normalisation, and Hall-A 2017 a bin-dependent band read from its figure 5;
# both normalisations are fitted as nuisances in fitrun.py, not added here.
SYST_N = 0.10

def _rows_xs(sel, meaning, ebeam=None, obs=("U", "LT", "TT"), rc=False, addsyst=0.0,
             ltscale=False):
    """ltscale converts the Drechsel-Tiator normalisation of the 2011 Hall-A paper.

    PRC 83 025201 writes the decomposition with eps_L where CLAS6 -- and the 2016,
    2017 and 2021 Hall-A papers, and we -- use eps, and defines eps_L/eps =
    4 M^2 xB^2 / Q^2 = gamma^2 just below its equation (9).  Matching the cos phi
    modulation gives sigma_LT(ours) = gamma * sigma_TL(theirs), and the same for
    sigma_LT'.  gamma is 0.43 to 0.51 at their settings, so their numbers are
    about twice ours.  sigma_U needs no conversion: the paper publishes the
    constant term itself.  sigma_TT neither: its weight is eps in both."""
    out = []
    for _, r in sel.iterrows():
        e = r.eps if np.isfinite(r.eps) else amp.epsilon(r.xB, r.Q2,
                ebeam if ebeam is not None else float(r.Ebeam))
        for o in obs:
            col = {"U": "s_u", "LT": "s_LT", "TT": "s_TT", "LTp": "s_LTp"}[o]
            sfx = "_rc" if rc else ""
            v = r[col + sfx] if o != "LTp" else r[col]
            st = r[{"U": "stat_U", "LT": "stat_LT", "TT": "stat_TT",
                    "LTp": "stat_LTp"}[o] + (sfx if o != "LTp" else "")]
            sy = r[{"U": "sys_U", "LT": "sys_LT", "TT": "sys_TT",
                    "LTp": "sys_LTp"}[o] + (sfx if o != "LTp" else "")]
            if not np.isfinite(v): continue
            err = math.hypot(st if np.isfinite(st) else 0.0,
                             sy if np.isfinite(sy) else 0.0)
            if ltscale and o in ("LT", "LTp"):
                g = 2*0.9382720813*r.xB/math.sqrt(r.Q2)
                v *= g; st = st*g if np.isfinite(st) else st
                sy = sy*g if np.isfinite(sy) else sy
                err *= g
            if addsyst: err = math.hypot(err, addsyst*abs(v))
            if err <= 0: err = SYST_N * abs(v)
            if err <= 0: continue
            out.append((float(r.Q2), float(r.xB), abs(float(r.t)),
                        meaning if o == "U" else o, float(v), err, float(e),
                        str(r.group)))
    return out

def _pred_xs(p, row, ch):
    Q2, xB, mt, obs, _, _, eps, _ = row
    s = amp.structure(p, ch, -mt, xB, Q2)
    if s is None: return None
    if obs == "U":   return s["T"] + eps*s["L"]
    if obs == "T":   return s["T"]
    if obs == "L":   return s["L"]
    if obs == "LT":  return s["LT"]
    if obs == "TT":  return s["TT"]
    if obs == "LTp": return s["LTp"]
    raise ValueError(obs)

def _mk_xs(key, label, sel, meaning, ch, ebeam=None, obs=("U","LT","TT"), rc=False,
           addsyst=0.0, ltscale=False):
    rows = _rows_xs(sel, meaning, ebeam, obs, rc, addsyst, ltscale)
    return dict(key=key, label=label, kind="xs", ch=ch, rows=rows,
                predict=lambda p, r, ch=ch: _pred_xs(p, r, ch))

def _mk_rosenbluth_p():
    """Proton sigma_L, PRL 117 262001 -- the fourth separated structure function,
    which never reached the workbook because All_experiment.xlsx has no column
    for it.  Fitted alongside sigma_T with stat and sys added in quadrature; the
    paper's warning that the two are anticorrelated cannot be honoured, the
    covariance having been lost."""
    import halla_rosenbluth as HR
    rows = []
    for Q2, tp, v, st, sy in HR.PROTON_SIGMA_L:
        mt = abs(-amp.tmin(amp.Mpi0, Q2, 0.36) + tp)
        e = math.hypot(st, sy)
        rows.append((Q2, 0.36, mt, "L", v, e, float("nan"), "ha16_L_%.2f" % Q2))
    return dict(key="halla_y16_L", kind="xs", ch="pi0p", rows=rows,
                label="Hall-A proton, $\\sigma_L$ separated (2016)",
                predict=lambda p, r, ch="pi0p": _pred_xs(p, r, ch))

def _mk_rosenbluth_n():
    """Neutron sigma_U at the two beam energies, PRL 118 222002 figure 4 -- the
    measured quantity, in place of the derived sigma_T."""
    import halla_rosenbluth as HR
    rows = []
    Q2, xB = HR.NEUTRON_Q2, HR.NEUTRON_XB
    for i, (tp, u1, s1, u2, s2) in enumerate(HR.NEUTRON_SIGMA_U):
        mt = abs(-amp.tmin(amp.Mpi0, Q2, xB) + tp)
        for u, st, sy, E in ((u1, s1, HR.NEUTRON_SYS[i], HR.NEUTRON_EBEAM[0]),
                             (u2, s2, HR.NEUTRON_SYS2[i], HR.NEUTRON_EBEAM[1])):
            eps = amp.epsilon(xB, Q2, E)
            rows.append((Q2, xB, mt, "U", u, math.hypot(st, sy), eps,
                         "ha17_U_%.3f" % E))
    return dict(key="halla_n_U", kind="xs", ch="pi0n", rows=rows,
                label="Hall-A neutron, $\\sigma_U$ at two $\\epsilon$ (2017)",
                predict=lambda p, r, ch="pi0n": _pred_xs(p, r, ch))

C6 = XS[XS.exp == "CLAS6_y12"]
SETS = [
  _mk_xs("clas6_pi0", "CLAS6 $\\pi^0$ structure functions",
         C6[C6.meson == "pi0"], "U", "pi0p", 5.75, rc=True),
  _mk_xs("clas6_eta", "CLAS6 $\\eta$ structure functions",
         C6[C6.meson == "eta"], "U", "etap", 5.75, rc=True),
  # The neutron sigma_T is REPLACED by the two beam energies' sigma_U, which is
  # what E07-007 measured; sigma_T is derived from them and is constructed to
  # have the longitudinal part removed, so it cannot see sigma_L.  sigma_LT and
  # sigma_TT keep their published values.  See halla_rosenbluth.py.
  _mk_xs("halla_n", "Hall-A neutron, $\\sigma_{LT},\\sigma_{TT}$ (2017)",
         XS[XS.exp == "HallA_y17"], "T", "pi0n", 5.55, obs=("LT","TT")),
  _mk_rosenbluth_n(),
  _mk_xs("halla_y16", "Hall-A proton, $\\sigma_T$ separated (2016)",
         XS[XS.exp == "HallA_y16"], "T", "pi0p", 5.55),
  _mk_rosenbluth_p(),
  _mk_xs("halla_y11", "Hall-A proton 6 GeV (2011)",
         XS[XS.exp == "HallA_y11"], "U", "pi0p", 5.752, obs=("U","LT","TT","LTp"),
         ltscale=True),
  _mk_xs("halla_y21", "Hall-A proton 12 GeV (2021)",
         XS[XS.exp == "HallA_y21"], "U", "pi0p", None, obs=("U","LT","TT","LTp")),
  _mk_xs("clas12_xs", "CLAS12 cross sections, preliminary",
         XS[XS.exp == "CLAS12_y25"], "U", "pi0p", 10.604),

]

# ---- COMPASS: predicted by averaging over the published (Q2,nu,|t|) grid ------
# PLB 870 139832 builds the cross section on a four-dimensional grid and combines
# the cells with the weight of the cell volume, equation (16).  One cell spans xB
# from 0.02 to 0.47, so evaluating the model at the quoted <Q2>, <xB> is not the
# measured quantity: it overshoots by a factor five to eight, while the grid
# average agrees to 1.0-1.7.
import compass_grid as _CG
_CROWS, _CBIN, _CCACHE = [], {}, {}
for _, _r in XS[XS.exp == "COMPASS_y25"].iterrows():
    _pr = str(_r.proj)
    _rng = (float(_r.bin_lo), float(_r.bin_up))
    for _o, _v, _st, _sy in (("U", _r.s_u, _r.stat_U, _r.sys_U),
                             ("TT", _r.s_TT, _r.stat_TT, _r.sys_TT)):
        if not np.isfinite(_v): continue
        _key = (_pr, _rng[0], _rng[1], _o)
        _CBIN[_key] = (_rng if _pr == "Q2" else None,
                       _rng if _pr == "nu" else None,
                       _rng if _pr in ("t", "ref27") else None)
        _CROWS.append((float(_r.Q2), float(_r.xB), abs(float(_r.t)), _o, float(_v),
                       math.hypot(_st, _sy), float(_r.eps), f"compass_{_pr}", _key))
def _pred_compass(p, row):
    ck = (id(p), row[8])
    if ck not in _CCACHE:
        q2r, nur, tr = _CBIN[row[8]]
        _CCACHE[ck] = _CG.average(p, row[8][3], q2r, nur, tr)
    return _CCACHE[ck]
SETS.append(dict(key="compass", label="COMPASS 2025, averaged over the grid",
                 kind="xs", ch="pi0p", rows=[r[:8] + (r[8],) for r in _CROWS],
                 predict=_pred_compass))

# ---- asymmetries -------------------------------------------------------------
def _asym_rows(sel, E):
    out = []
    for _, r in sel.iterrows():
        err = math.hypot(r.stat if np.isfinite(r.stat) else 0.0,
                         r.syst if np.isfinite(r.syst) else 0.0)
        if err <= 0: continue
        out.append((float(r.Q2), float(r.xB), abs(float(r.t)), r.observable,
                    float(r.value), err, amp.epsilon(r.xB, r.Q2, E), "all"))
    return out

def _pred_bsa(p, row, ch, E):
    Q2, xB, mt, obs, _, _, _, _ = row
    if obs == "A_LU^sinphi":
        return amp.bsa_sinphi(p, ch, -mt, xB, Q2, E)
    s = amp.structure(p, ch, -mt, xB, Q2)
    if s is None: return None
    e = amp.epsilon(xB, Q2, E)
    return s["LTp"]/(s["T"] + e*s["L"])

MP, MPI = 0.9382720813, 0.1349768
def above_threshold(Q2, xB, mt, m=MPI):
    W2 = MP*MP + Q2*(1-xB)/xB; W = math.sqrt(W2)
    Eg = (W2 - Q2 - MP*MP)/(2*W); pg = math.sqrt(Eg*Eg + Q2)
    Ep = (W2 + m*m - MP*MP)/(2*W); pp = math.sqrt(max(Ep*Ep - m*m, 0.0))
    return mt > -(m*m - Q2 - 2*(Eg*Ep - pg*pp))

# The digitised De Masi moments (CLAS6_demasi on the `asymmetries` sheet) were
# dropped on 2026-10-05 (VPK).  They are figure 5 read off by eye: errors 1.5x
# too small by the workbook's own note, bin labels instead of measured means,
# and nine of the sixty measured bins missing.  The moment extracted from the
# phi distributions supersedes them and reproduces them to a median 2.0 %,
# which is the accuracy of reading a plot.  Keeping them as a "cross-check"
# only invited the chart to show two bars with almost the same name.
for key, label, exp, ch, E in (
        ("bsa_zhao",   "CLAS6 $\\eta$ beam-spin asymmetry",  "CLAS6_zhao",   "etap", 5.776),
        ("bsa_clas12", "CLAS12 $\\sigma_{LT'}/\\sigma_0$",   "CLAS12_y24",   "pi0p", 10.6)):
    m = MPI if ch == "pi0p" else 0.547862
    rows = [r for r in _asym_rows(AS[AS.exp == exp], E)
            if above_threshold(r[0], r[1], r[2], m)]
    SETS.append(dict(key=key, label=label, kind="asym", ch=ch, rows=rows,
                     predict=lambda p, r, ch=ch, E=E: _pred_bsa(p, r, ch, E)))

# ---- De Masi: the sin(phi) moment, extracted from the phi distributions -------
# The paper publishes PLOTS, not a table.  The measurement is the phi
# distribution in the CLAS database (sheet bsa_phi, 703 points / 60 bins); the
# `asymmetries` sheet's 62 moment values were digitised from figure 5 and their
# errors are 1.5x too small, so they are not fitted.
#
# Until 2026-10-05 the fit used the 703 phi points directly.  That was wrong in
# one specific way: it counted 663 points (40 were cut for error > 0.5) while
# carrying about SIXTY independent numbers, so De Masi supplied 42 % of the
# fit's points and inflated every chi2/ndf we quoted.  Measured: scaling
# sigma_LT in the denominator by 0, 2 or -1 moved chi2 by 0.4 to 3.1 units out
# of 625, i.e. the phi shape beyond sin(phi) carries no information.
#
# The moment is extracted ONCE by demasi_moments.py, with the FULL asymmetry
# form rather than a plain sin(phi) fit.  In 2008 the denominator was unknown;
# today sigma_U and sigma_TT are measured at these kinematics (CLAS6 pi0, mean
# pull -0.10 and +0.09), so the extracted number is the physical amplitude
#     A = sqrt(2 eps (1-eps)) sigma_LT' / (sigma_T + eps sigma_L)
# and not its projection onto sin(phi).  The two differ by 7.6 %.
#
# Model dependence of the extraction, measured: swapping amp2609's denominator
# for the Htilde fit's moves A by a median 0.009 of the statistical error.
#
# Errors are STATISTICAL ONLY (VPK, 2026-10-05).  The published 0.016 covers
# "event selection plus the choice of the fit function"; the second is removed
# by construction, and a phi-independent first barely reaches the moment --
# the leakage sum(w sin)/sum(w sin^2) has median 0.023, so 0.016 moves A by
# 0.0004.  The 3.5 % beam polarisation stays as a fitted nuisance.
_dmrows = []
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "data", "demasi_moments.txt")) as _fh:
    for _ln in _fh:
        if _ln.startswith("#"): continue
        _b, _q, _x, _t, _a, _e, _n = _ln.split()
        _dmrows.append((float(_q), float(_x), float(_t), "A_LU^sinphi",
                        float(_a), float(_e),
                        amp.epsilon(float(_x), float(_q), 5.776), _b))
SETS.append(dict(key="bsa_demasi_mom", kind="asym", ch="pi0p", rows=_dmrows,
                 label="CLAS6 $\\pi^0$ BSA moment, from the $\\phi$ distributions",
                 predict=lambda p, r: amp.bsa_sinphi(p, "pi0p", -r[2], r[1], r[0], 5.776)))

MOM = {"A_UL^sinphi": "AULsin", "A_UL^sin2phi": "AULsin2",
       "A_LL^const": "ALLc", "A_LL^cosphi": "ALLcos"}
_eg = []
for _, _r in AS[AS.exp == "eg1dvcs"].iterrows():
    _e = math.hypot(_r.stat if np.isfinite(_r.stat) else 0.0,
                    _r.syst if np.isfinite(_r.syst) else 0.0)
    if _e <= 0: continue
    _k = MOM[_r.observable]
    # Slot 7 is the panel label.  It used to hold the moment name, which put the
    # two (Q2, xB) settings of the experiment into one panel: the model then
    # zig-zagged, because consecutive points in t belonged to different
    # kinematics.  It must be the setting.
    _eg.append((float(_r.Q2), float(_r.xB), abs(float(_r.t)), _k, float(_r.value), _e,
                amp.epsilon(float(_r.xB), float(_r.Q2), 5.9),
                f"Q2={float(_r.Q2):.2f} xB={float(_r.xB):.3f}"))
def _pred_eg1(p, row):
    Q2, xB, mt, key, _, _, e, _ = row
    s = amp.structure(p, "pi0p", -mt, xB, Q2)
    if s is None: return None
    s0 = s["T"] + e*s["L"]
    return dict(AULsin=math.sqrt(2*e*(1+e))*s["LTUL"]/s0,
                AULsin2=e*s["TTUL"]/s0,
                ALLc=math.sqrt(1-e*e)*s["Tp"]/s0,
                ALLcos=math.sqrt(2*e*(1-e))*s["LTpLL"]/s0)[key]
SETS.append(dict(key="eg1", label="eg1-dvcs, polarised target", kind="asym",
                 ch="pi0p", rows=_eg, predict=lambda p, r: _pred_eg1(p, r)))

BY = {s["key"]: s for s in SETS}
ALL = [s["key"] for s in SETS]

def chi2(p, key, scale=1.0):
    """scale is the fitted normalisation of the set, 1 when it has none."""
    s = BY[key]; c = 0.0; n = 0
    for r in s["rows"]:
        v = s["predict"](p, r)
        if v is None: continue
        c += ((r[4] - scale*v)/r[5])**2; n += 1
    return c, n

if __name__ == "__main__":
    print(f"{'set':>12} {'kind':>5} {'points':>7}   label")
    for s in SETS:
        print(f"{s['key']:>12} {s['kind']:>5} {len(s['rows']):7d}   {s['label']}")
    print(f"{'':>12} {'':>5} {sum(len(s['rows']) for s in SETS):7d}   total")
