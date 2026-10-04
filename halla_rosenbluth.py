"""The Rosenbluth-separated Hall A data that never reached the workbook.

E07-007 published four separated structure functions for the proton
(PRL 117 262001) and for the neutron (PRL 118 222002).  The workbook carries
three of them: `All_experiment.xlsx` has exactly nine cross-section columns --
s_u, s_LT, s_TT with their errors -- a layout built for UNSEPARATED sets, so
sigma_T was put in the sigma_U column and flagged through `sigma_meaning`, and
**sigma_L had nowhere to go**.  It was lost silently, in both papers, for the
one quantity a Rosenbluth separation exists to produce.

Numbers here are the authors' own (M. Defurne, private communication, 2026-10),
not digitised: they reproduce the workbook's sigma_T, sigma_LT and sigma_TT to
1e-13.

PROTON: sigma_L with its own stat and sys.  Fitted alongside sigma_T.  The paper
warns the two are "strongly anticorrelated" and the covariance is gone -- a farm
incident at Saclay destroyed it -- so this ignores the correlation.  That is the
usual approximation when a covariance is unavailable and it is stated, not
hidden.  Defurne's own suggestion was to rebuild sigma_T + eps sigma_L with a
flat 10 % error; that is NOT used, because a 10 % per-point error washes out the
difference between the two beam energies, which is the only place sigma_L lives.

NEUTRON: the two beam energies' sigma_U, which is what was MEASURED; sigma_T and
sigma_L are derived from them by solving a 2x2 system.  Checked: sigma_T +
eps sigma_L reproduces the published sigma_U to 0.02 %.  These two are genuinely
independent (different beams, different data), so no covariance is needed.  This
REPLACES the four sigma_T points rather than adding to them -- they are the same
measurement and fitting both would count it twice.

Why it matters, measured on amp2609 vs amp2612 (the fit with Htilde):
    sigma_T alone, 4 points     chi2 0.78 vs 1.26   -- cannot tell them apart
    sigma_U at two eps, 8 pts   chi2 18.9 vs 2345   -- separates them completely
sigma_T is constructed to have the longitudinal part removed, so it cannot see
sigma_L at all.  That is why the fit was content while amp2612 gave R = 1.65 for
pi0 off the neutron, and why the generator's stage-2 bound came out at 21.5
against 1.8 for amp2609.
"""
import math

# --- proton, PRL 117 262001 ------------------------------------------------
# Q2, t' = tmin - t, sigma_L [nb/GeV^2], stat, sys
PROTON_SIGMA_L = [
    (1.50, 0.0236977,  459.892, 221.428, 360.343),
    (1.50, 0.0727075,  205.479, 270.141, 339.126),
    (1.75, 0.0251812, -136.373, 150.687, 379.239),
    (1.75, 0.0742083, -259.260, 180.713, 412.568),
    (1.75, 0.1232730, -409.029, 209.363, 409.800),
    (2.00, 0.0255560,  205.509,  97.979, 152.368),
    (2.00, 0.0750411, -180.349, 130.282, 190.923),
]
PROTON_EBEAM = {1.50: (3.355, 5.55), 1.75: (4.455, 5.55), 2.00: (4.455, 5.55)}

# --- neutron, PRL 118 222002, figure 4 -------------------------------------
# t', sigma_U at E = 4.4557 and at E = 5.5489 [nb/GeV^2] with stat
NEUTRON_SIGMA_U = [
    (0.0250393, 308.148, 33.9781, 317.742, 30.7749),
    (0.0744954, 361.971, 30.2007, 331.060, 28.7380),
    (0.1240390, 334.318, 31.0850, 351.049, 28.6948),
    (0.1736350, 378.598, 33.7057, 371.694, 31.7046),
]
NEUTRON_Q2, NEUTRON_XB = 1.75, 0.36
NEUTRON_EBEAM = (4.4557, 5.5489)
# the systematics of figure 4, high/low, symmetrised the way the workbook does
NEUTRON_SYS = [(0.025949+0.028127)/2*1000, (0.055973+0.00379103)/2*1000,
               (0.029963+0.002527)/2*1000, (0.00907499+0.024052)/2*1000]
NEUTRON_SYS2 = [(0.02104+0.035251)/2*1000, (0.026804+0.001456)/2*1000,
                (0.032528+0.019711)/2*1000, (0.021793+0.008793)/2*1000]
