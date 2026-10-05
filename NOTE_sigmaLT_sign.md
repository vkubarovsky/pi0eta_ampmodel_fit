# sigma_LT on the proton: where the sign lives, and three points nobody describes

**2026-10-04.  Joint result of this repository and the `hepgen_fit` session
(GK twist-3, `libGKPi0`), reached by each side running the other's test.**

Both sides retracted something to get here; both retractions are recorded below.

## 1. The dictionary term that GK does not have

    sigma_LT = -sqrt2 Re[ T00p* T01p  +  T00m* T01m ]
                          \________/     \________/
                           non-flip         flip

GK's `libGKPi0.cpp:1339` reads `result.Mpmp0 = result.Mppp0`.  That is the
C = E relation, and it makes the non-flip interference vanish **identically**:
GK's sigma_LT is one term, ours is two.  In this repository C = E holds at
level 2 and is unlocked at level 3 by `rho_CE = |T01p|/|U01p|`, slot 22,
fitted to **0.218**.

Measured at Q2 = 1.75, xB = 0.36, |t| = 0.232, turning `rho_CE` off:

| observable | with | without | change |
|---|---|---|---|
| sigma_T  |  802.35 |  793.62 | -1.1 % |
| sigma_L  |   94.72 |   94.72 | none |
| sigma_TT | -174.52 | -183.24 | -5 % |
| **sigma_LT** | **-17.75** | **-54.90** | **x3.1** |

**sigma_LT is the only observable that sees C != E.**  That is why nobody fits
it, and why 0.218 must not be read as a measurement: one number held by seven
points, two of which are the pair in section 3.  GK cannot test C != E even in
principle -- there is no parameter attached to that line.

## 2. The observable angle, measured

**Revised 2026-10-05 after a comparison trap caught both of us.**  The first
version of this section compared GK's INTERNAL phase,
`arg<Etilde> - arg<H_T> = 122-124 deg`, against a +-90 window.  That quantity is
not observable: the prefactors on Re and Im carry different signs, so GK's
*delivered* `sigma_LT` and `sigma_LT'` sit at `180 - delta`.  Reported as a GK
property, 122 deg would have said GK agrees with the data.  It does not.

The convention-free quantity is the angle of the line that the two structure
functions lie on,

    delta_obs = atan2(sigma_LT', sigma_LT)  mod 180

**Estimate it by profiling the likelihood in that angle, not by averaging the
ratio.**  The ratio average is badly biased: `sigma_r ~ |r|` makes the weight go
as `1/r^2`, and dividing by a `sigma_LT` with large errors is regression
dilution.  Four estimators on the same data span -0.58 to -1.64.  The profile is

    chi2(d) = sum (y cos d - x sin d)^2 / [ (sy cos d)^2 + (sx sin d)^2 ]

with x = sigma_LT, y = sigma_LT': no division, no cut, every point, and bounded
where the ratio runs away.

| | delta_obs | 68 % | chi2/ndf |
|---|---|---|---|
| **data, HallA_y11** | **121.4 deg** | [116.3, 127.2] | 0.78 |
| **data, HallA_y21** | **126.6 deg** | [123.1, 130.1] | 2.03 |
| this fit (Mom_B12), at y11 / y21 | 101.8 / 115.8 deg | rms 3.7 / 11.3 | |
| GK, from its delivered structure functions | **57 deg** | | |

The data put the ray firmly in the second quadrant -- `sigma_LT < 0`,
`sigma_LT' > 0`: at y11, 13 of 14 significant `sigma_LT` are negative and 16 of
16 `sigma_LT'` positive; at y21, 25 of 28 and 26 of 27.  GK sits in the FIRST
quadrant, both positive.

**Ranking: data 121-127, this model 102-116, GK 57.**  We are 10-20 deg out and
outside the intervals.  GK is 65-75 deg out and in the wrong quadrant.

GK's `sigma_LT'` is also a factor 4.4 small -- median 0.0197 against 0.0868 on
the 60 De Masi moments, with pulls +2.5 to +4.1 spread across the whole
kinematic range.  Combined with the rotation that is ONE statement, not two:
GK's `<Etilde>* <H_T>` product is too small and rotated by about 70 degrees.

**A caution discharged, on the GK side.**  `getCXLU` does not exist in the
original HEPGen++ library; it was written as the imaginary partner of
`getCXLT` and had never been compared to a published curve.  Against the 60
moments: model positive 60/60, data 59/60, same sign 59/60, and flipping the
model takes chi2 from 375 to 906.  The sign is validated.

## 3. Neither flip term carries sigma_LT alone

| | chi2 on the 7 HallA_y16 sigma_LT points |
|---|---|
| ours, both terms | **37.8** |
| ours, flip only (`rho_CE = 0`) | 107.9 |
| GK, own sign | 76.3 |
| GK, sign flipped by hand (a diagnostic, not a model) | 33.0 |

Our flip term alone **overshoots by three** at Q2 = 1.75 and 2.00 (-54.9
against a measured -17.2) with |T00m|/|T01m| = 0.19-0.23.  GK's, flip-only by
construction, is about right there with the same ratio at 0.012-0.025 -- **ten
times smaller**.  Two frameworks a factor ten apart in one ratio, both working
once each is allowed its own treatment of the non-flip term.

So the earlier claim "GK's flip amplitude is a factor 10-20 too small" was
**withdrawn by its author**: it presumes the flip term should carry sigma_LT
alone, and our own numbers show it should not.  The defensible form: GK's flip
amplitude is small compared with the one this fit prefers.

**sigma_LT cannot be read as a measurement of `<Etilde>` in either framework.**

## 4. Three points that neither model describes

Where the chi2 sits, over the same seven points:

| | Q2 = 1.50 pair | (1.75, 0.281) | the other four |
|---|---|---|---|
| this fit (37.8) | 21 = **56 %** | 12 = 32 % | 5 = 13 % |
| GK, sign-flipped (33.0) | 24.1 = **73 %** | 6.5 = 20 % | 2.3 = 7 % |

Same three points, same order, from models built on different assumptions with
flip amplitudes a factor ten apart.  **Four of seven points are described by
both; three are described by neither.**  That is evidence about the data.

Look at the data alone, same xB and nearly the same |t|:

| Q2 | sigma_LT (nb/GeV^2) |
|---|---|
| 1.50 | -100.1, -169.2 |
| 1.75 | -16.9, -17.2, **+13.1** |
| 2.00 | -1.6, -6.2 |

A factor six between Q2 = 1.50 and 1.75 and another three to 2.00: an
effective power near Q^-13, which no twist expansion produces.  And the point
at (1.75, 0.281) is positive where its two neighbours at the same Q2 are
negative.

The claim worth making is the weak one: **no smooth model reaches these three
points.**  Not that the measurements are wrong.

## Caveats carried forward

* The GK sign flip is a hand diagnostic with the phase fixed and nothing
  refitted.  It tests magnitude only.
* `rho_CE = 0` is the closest this model comes to GK's structure; it is not
  GK's structure.
* An earlier version of this note claimed the data require
  `arg<Etilde> - arg<H_T>` inside +-90 deg and that GK computes it outside.
  That framing is withdrawn: it applied to the SIGN of `sigma_LT` alone, where
  this model has two terms and GK one, and it was stated in a convention that is
  not observable.  Section 2 above replaces it.
* Both models' relative u-d phases remain free or frozen rather than measured,
  and the Hall-A flavour separation itself assumes zero relative phase --
  see the sigma_TT n/p discussion in `MEMORY`/`project_pi0_eta_flavor_model`.
