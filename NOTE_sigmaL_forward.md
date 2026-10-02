# sigma_L died at the forward peak, because the model had no H~

**Found and fixed 2026-09-30.  Branch `htilde`.**

## The defect

`amplitudes.py` built both longitudinal amplitudes with the same kinematic
prefactor:

    kin  = -t'/(8 m^2)
    Lmag = sqrt(A*kin) * <L>              # T00m, recoil FLIP      -- correct
    T00p = rho_nf * |Lmag| * exp(i phi_nf) # recoil NON-FLIP        -- WRONG

so `sigma_L = |T00p|^2 + |T00m|^2` was proportional to `(-t')` and vanished at
the forward peak.  At `-t' = 0.002`, `R = sigma_L/sigma_T = 0.0009`.  76 % of
`sigma_L` was carried by `T00p`, the amplitude that was wrong.

This contradicts the angular-momentum table in our own amplitude note
(`~/pi0_amplitudes/notes/pi0_amplitude.tex`): with `Delta = mu - nu + nu'`, an
amplitude vanishes like `(sqrt(-t'))^|Delta|`, and `M_{0+,0+}` has `Delta = 0`.
The longitudinal non-flip amplitude is a **constant** as `t' -> 0`.

## The cause

The longitudinal sector was a *single* GFF block carrying the kinematic
prefactor of `<E~>`.  The GPD that populates the non-flip amplitude, `<H~>`, was
**absent from the model entirely**.  `T00p = rho_nf |T00m|` was a stand-in that
inherited the flip amplitude's `sqrt(-t')`.

## The fix

Slots 37-42, active only when `len(p) > 41`:

    T00p = sqrt(A (1-xi^2)) [ <H~> - xi^2/(1-xi^2) <E~> ]
    T00m = sqrt(A (-t')/(8 m^2)) xi^nx <E~>

`<H~>` has the same functional form as the other blocks
(`N exp[(b + b' ln xB) t] (Q^2)^{nQ/2}`, `H~^d = R H~^u`), its u-d relative
phase shared with `<E~>` (p[36]), and `arg<H~>` in the reused slot p[26].  The
old L block is re-read as `<E~>`; slot 25 (`rho_nf`) is dead.

`nx` (slot 42) is the power of `xi` written explicitly in front of `<E~>`, as
in GK.  It used to be absorbed in the fitted normalisation, which is not
harmless: `xi` varies by a factor 8 over the data, so absorbing it distorted the
xB dependence of the block.

**Backward compatible by vector length.**  A 37-slot fitpar takes the old code
path byte for byte -- `sigma_L(-t'=0.3) = 4.12661` before and after -- so
`amp2609` in the generator and every earlier run reproduce unchanged.

## Result

Both models refitted with identical machinery and effort (`NPAR=37` puts
`fitrun.py` back on the old configuration: 25 free + 3 norms = 28, exactly the
published one).  Four independent starts each.

| model | chi2 | npar | ndf | chi2/ndf | R at -t'=0.002 |
|---|---|---|---|---|---|
| no H~, published (`C1_with_compass`) | 2140.5 | 28 | 1527 | 1.4018 | 0.0009 |
| no H~, refitted control (`old_ctl1`) | 2140.3 | 28 | 1527 | 1.4017 | 0.0009 |
| H~, nx = 0 (`Htil_xi0_rs`)           | 2275.7 | 32 | 1523 | 1.4942 | 0.326 |
| **H~, nx = 1 (`Htil_free2`)**        | **2065.8** | 32 | 1523 | **1.3564** | **0.156** |

`Delta chi2 = -74.5` for 4 parameters, `p = 2.5e-15`.  The control reproduces
the published number to 0.2 units, so the old fit was converged and the
improvement is real.

Where it comes from -- concentrated in the sets that reach small `-t'` and
measure `sigma_LT'`, which is exactly what should constrain the forward
longitudinal amplitude:

    halla_y21  -38.8      halla_y16  +6.7
    clas6_eta  -21.7      bsa_zhao   +3.0
    halla_y11  -12.9      demasi     +1.9
    halla_n     -5.1
    compass     -5.0
    eg1         -3.1

The data clearly want the explicit `xi` as well: `nx = 0` is 210 units worse
than `nx = 1` at the same number of parameters.

## How well is sigma_L forward pinned down?

There is **no sigma_L measurement in the database at all**.  The longitudinal
sector enters only through `eps sigma_L` inside `sigma_U` and through
`sigma_LT`, `sigma_LT'`.  A profile in `N_H~` (12 nodes, six overlapping
continuation chains in both directions, `runs/scan_Ht/`) gives

    within the optimiser noise floor:  N_H~ 8-16,  R_fwd 0.096-0.214,
                                       R(-t'=0.3) 0.067-0.107
    clearly excluded:  N_H~ = 0 (+292), 3 (+39), 6 (+33), 19 (+44), 28 (+149)

The resolution is set by **optimiser noise, not by chi2 curvature**: independent
chains visiting the same node disagree by a median of 26.7 units (range 1.2 to
88) on this 29-parameter surface.  A `Delta chi2 = 1` interval would be fiction.
The true uncertainty from the data is **wider** than the interval above, not
narrower.  A defensible band needs MCMC or many more restarts per node.

## Checks

* Positivity of the bracket `B(phi; h, lam)` over the whole campaign box
  (Q^2 1-6, W >= 1.8, four channels, h and lam = -1,0,+1, 49 phi values):
  0 violations in 39600.  But the new fit sits closer to the boundary than the
  old one, `min B/sigma0 = +0.0009` (pi0 n, Q^2 = 4, xB = 0.35, -t' = 0.4,
  h = -1, lam = +1) against `+0.0055`.  Not a violation; worth watching if this
  model is installed in the generator with full target polarisation.
* Polarised-target closure (`polarised.check_against_fit`): deviation exactly 0,
  so the A_UL / A_LL machinery is unaffected.

## What H~ does NOT fix

It does not resolve the `sigma_LT` sign disagreement.  The model is negative
everywhere before and after, and slightly *more* negative after, while CLAS12
measures `+31.6 +- 6.0` at `-t' = 0.112` and `+16.0 +- 1.8` at `0.261`.  The
experiments disagree with each other on the sign of `sigma_LT`; that is
untouched by the longitudinal sector.  See the earlier commit b6f4ffa.

## A wrong turn worth recording

The first new-model fit gave 2131.7 and was reported as "the fix costs nothing,
Delta chi2 = -8.8 for 4 parameters is not significant".  That was an
unconverged minimum.  A constrained scan then found 2097 at fixed `N_H~ = 11`
with *fewer* free parameters, which is impossible at a true minimum and is what
exposed it.  The converged answer is -74.5, decisively significant.  On a
surface this rough, a single `least_squares` from two seeds is not a fit.
