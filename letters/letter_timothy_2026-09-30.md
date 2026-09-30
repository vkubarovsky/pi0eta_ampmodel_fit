# Draft letter to Timothy — 2026-09-30

Dear Timothy,

I agree completely that separating σ_L and σ_T for pseudoscalar mesons is
important. Let me stress one point, though. A Rosenbluth separation rests
entirely on **absolute** cross sections: the whole signal is a small difference
between two absolute measurements at different ε, so any relative
normalisation uncertainty propagates straight into σ_L. For that reason I would
be careful with conclusions drawn from the unreleased CLAS12 RG-A and RG-K
data.

I have made an attempt to analyse all the published π⁰ and η electroproduction
data off the proton and the neutron in a single amplitude-level fit. The
database now holds 13 data sets from CLAS6, CLAS12, Hall A and COMPASS:

| # | experiment | channel | E_beam [GeV] | ε | observables |
|---|---|---|---|---|---|
| 1 | CLAS6, PRC 90 025205 (2014) | π⁰ p | 5.75 | 0.31 – 0.71 | σ_U, σ_LT, σ_TT |
| 2 | CLAS6, PRC 95 035202 (2017) | η p | 5.75 | 0.31 – 0.71 | σ_U, σ_LT, σ_TT |
| 3 | Hall A, PRL 117 262001 (2016) | π⁰ p | 3.355 + 5.55, 4.455 + 5.55 | 0.13 / 0.72 … 0.52 / 0.85 | **σ_T and σ_L separated** |
| 4 | Hall A, PRL 118 222002 (2017) | π⁰ n | E07-007 pair | as above | **σ_T and σ_L separated** |
| 5 | Hall A 6 GeV (2011) | π⁰ p | 5.752 | 0.65 – 0.77 | σ_U, σ_LT, σ_TT, σ_LT' |
| 6 | Hall A, PRL 127 152301 (2021) | π⁰ p | 4.49, 7.38, 8.52, 8.85, 10.59, 10.99 | 0.50 – 0.72 | σ_U, σ_LT, σ_TT, σ_LT' |
| 7 | COMPASS, PLB 805 135454 (2020) | π⁰ p | 160 | 0.997 | σ_U, σ_LT, σ_TT |
| 8 | COMPASS, PLB 870 139832 (2025) | π⁰ p | 160 | 0.99 – 0.999 | σ_U, σ_LT, σ_TT |
| 9 | CLAS12 (preliminary) | π⁰ p | 10.604 | 0.52 – 0.95 | σ_U, σ_LT, σ_TT — *not used in the fit* |
| 10 | CLAS6 De Masi, PRC 77 042201 | π⁰ p | 5.776 | — | A_LU(φ) |
| 11 | CLAS6 Zhao, PLB 789 426 | η p | 5.776 | — | A_LU^sinφ |
| 12 | CLAS12 Kim, PLB 849 138459 | π⁰ p | 10.6 | — | σ_LT'/σ_0 |
| 13 | eg1-dvcs | π⁰ p (pol. target) | 5.9 | — | A_UL, A_LL moments |

As you can see, the coverage in ε is already rather wide — from 0.13 at the
lowest Hall A setting up to 0.999 at COMPASS. What is *not* wide is the overlap:
the only place where two independent measurements sit at the same (Q², x_B, t)
with different ε is a single CLAS6 / Hall A 2021 cell, and it is useless because
of the CLAS6 systematics there. The genuine lever is CLAS6 at 5.75 GeV against
CLAS12 at 10.6 GeV — about 120 matched cells with Δε ≈ 0.4 — and that is exactly
where the relative normalisation, not the statistics, sets the limit.

I attach a typical result: σ_L, σ_T and R = σ_L/σ_T as functions of Q² at
x_B = 0.2 and |t| = 0.3 (GeV/c)², with the uncertainty band from the fit
(the band is the parameter covariance at fixed functional form, 16–84 %). The
grey strip marks the Q² range in which fitted data actually exist at this
(x_B, |t|); outside it the curve is the functional form extrapolating. I have
the full set of such curves for other x_B and |t| but do not want to overload
you with pictures. **The result is preliminary.**

On the same figure I show the Goloskokov–Kroll prediction for the same
structure functions, computed with their own code (libGKPi0, the 2023/2024
parameter set). The picture is interesting: there is a clear disagreement
between the model and the extraction. Three specific points:

* σ_L agrees surprisingly well — the two curves cross around Q² ≈ 3.5 and stay
  within the band above it;
* σ_T does not: GK falls much faster, by a factor ≈ 3.6 at Q² = 6;
* consequently R runs the opposite way in the two: our extraction has R
  decreasing with Q² (0.11 → 0.03 between Q² = 1 and 6), GK has it rising
  (0.02 → 0.12).

The same disagreement shows up independently in the COMPASS Q² and ν
dependences: in their Table 11 kinematics GK gives σ_U ∝ (Q²)^−1.9 against
(Q²)^−0.86 ± 0.16 measured, and ν^−1.4 against ν^−2.75 ± 0.11.

This underlines once again the problem you raised in your letter.

Best regards,
Valery
