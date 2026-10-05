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

## 2. The sign is a phase, and GK computes it outside the window

| | arg<Etilde> - arg<H_T> |
|---|---|
| GK, **computed** from the convolutions | **122-124 deg**, flat in Q2 and \|t\| |
| this fit, **fitted** (slots 19, 20, 26) | 52 deg |
| what the data require | inside +-90 deg |

Beyond 90 deg the cosine turns negative and the explicit `-sqrt2` makes
sigma_LT **positive**; GK's is positive at all seven Hall-A points, the data
are negative at six of them.  The sign therefore sits in the phases, not in
the dictionary and not in a prefactor -- which was the alternative, and the
alarming one.

Two statements, independent, neither claiming the other's freedom:

* **what our fit establishes** -- the data want that phase difference inside
  90 deg.  Our model does not *predict* the sign; it has free phases and the
  data set them.
* **what the GK probe establishes** -- GK *computes* it outside.  Its
  `<Etilde>` sector is frozen end to end (normalisations 14 and 4, delta 0.48,
  alphastr 0.45, slope 0.9, no setter), so neither sign nor size of sigma_LT
  is reachable by any choice of GK's 16 fitted parameters.

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
* Both models' relative u-d phases remain free or frozen rather than measured,
  and the Hall-A flavour separation itself assumes zero relative phase --
  see the sigma_TT n/p discussion in `MEMORY`/`project_pi0_eta_flavor_model`.
