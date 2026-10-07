# Trying to fit photoproduction: why it cannot be done with this model

**2026-10-06.  VPK asked to try fitting something from real photoproduction, to
pin the scale of the Q2 damping by data rather than by taste.  The attempt
fails, and the reason is not the Q2 behaviour we spent the evening on.**

## The data

CLAS g12, Kunkel et al. PRC **98** 015207 (2018), figure 4: `gamma p -> pi0 p`
at six energies, `W = 2.49, 2.64, 2.79, 2.94, 3.08, 3.17 GeV`, `d sigma/dt` in
nb/GeV^2 -- the same object and the same units as our `sigma_T` at `Q2 = 0`,
since `sigma_L` is identically zero for a real photon.

HEPData is behind a bot check that neither curl nor the browser pane would
pass, so the figure was digitised at 400 dpi: red markers selected by
horizontal thickness (a marker is >= 8 px wide, an error bar 2-3), the median
marker row taken per image column, then the median over columns in bins of
0.15 GeV^2.  Calibration from the tick marks, not guessed: 90.33 px per decade
in sigma, panel-width/8 per GeV^2 in |t|.  60 points, `data/kunkel_g12_sigmaT_Q2zero.txt`.
The published uncertainties are not captured.

## Why a Q2 regulator cannot work

`sigma_T ~ (Q2)^p` at small Q2 along fixed W, and **p depends on |t|**:

| \|t\| | 0.42 | 0.72 | 1.02 | 1.47 | 2.07 |
|---|---|---|---|---|---|
| p | -0.568 | -0.347 | -0.126 | **+0.206** | **+0.648** |

independent of W.  Above |t| ~ 1.2 the model already converges; below it
diverges.  A factor `[Q2/(Q2+M^2)]^k` multiplies `sigma_T` by `(Q2)^(2k)` and
therefore cures exactly one value of p while over-suppressing every other.

## Where the t-dependence comes from, and the real defect

The amplitude carries `exp[(b + b' ln xB) t] = exp(b t) * xB^(b' t)`, and at
fixed W, `xB ∝ Q2`.  So `sigma_T` picks up `(Q2)^(2 b' t)`:
`dp/d|t| = 2b' = -2.611`, which is exactly the table above.  The Q2 divergence
IS the xB dependence of the t-slope.

And that slope is wrong at the photon point by an order of magnitude:

| W | t-slope of the DATA [GeV^-2] | model at Q2 = 1e-4 | model at Q2 = 1 |
|---|---|---|---|
| 2.49 | 1.14 | 14.21 | 2.41 |
| 2.79 | 1.06 | 14.55 | 2.70 |
| 2.94 | 1.20 | 14.70 | 2.83 |
| 3.17 | 1.90 | 14.92 | 3.03 |

`b + b' ln xB` with `b' = -1.305` grows without bound as `xB -> 0`: 15 GeV^-2
at the photon point against 1-2 measured.

**This is not a regulator problem.**  Capping the logarithm, `ln xB ->
ln(xB + x0)`, does not rescue it either: the model's slope DECREASES with
increasing xB -- 2.66 at xB = 0.13, 0.90 at xB = 0.5 -- so to reach the
measured 1.1 at `xB -> 0` the slope would have to turn around and fall, which
no choice of `x0` in a monotone form can do.  The measured photon-point slope
is the value this model gives at `xB ~ 0.45`.

## What this means for the Q2 damping

The damping of slots 43 remains correct for what it does: it makes
`sigma_L = O(Q2)` and `R = sigma_L/sigma_T -> 0`, which gauge invariance
requires and which Kroll confirmed.  But it repairs a RATIO in a limit the
model cannot reach for an independent reason.  Writing "the model now has the
correct photon-point behaviour" would be false.  The honest statement:

* `sigma_L/sigma_T -> 0` as required -- fixed by slot 43;
* `sigma_T` diverges, |t|-dependently, and is 10x too steep in |t| -- NOT fixed,
  and not fixable by any function of Q2 alone;
* so the model must not be used below roughly `Q2 = 0.5`, and the generator box
  starting at 0.8 is doing real work.

## What would be needed

A t-slope whose xB dependence saturates, fitted to photoproduction and to the
electroproduction data together.  That is a change to the xB sector, not the Q2
sector, and it would move the fitted region as well: at `xB = 0.13` the slope
is 2.66 and any saturating form changes it.  Not a night's work, and not to be
done without deciding first whether this model is meant to describe
`Q2 < 1` at all.

## The box is what makes the divergence harmless, and that is conditional

Measured across two campaigns and two models (Launch generator's amp2609
reference, 200 000 events/arm, and the amp2614 run here), the generated sample
simply does not reach the region where the slope parameter goes bad:

    born  xB < 0.128   8.1 % (amp2614)   8.28 % (amp2609)
    rad   xB < 0.128  17.5 %            16.98 %
    both  xB < 0.05     0.00 %            0.00 %

Radiation does not help it get there either.  Zeroing the model below an xB
threshold and recomputing `sigma_obs` gives 0.00 % of the weight from vertex
`xB < 0.05` at every box point tested; the one point with any weight below
0.128 is the one whose OBSERVED xB is already 0.075.  Radiation lowers `Q2t`
and `W2t` together, so vertex xB barely moves.  The exposure is entirely
through the observed kinematics, i.e. through the box.

The rad/born factor two below `xB = 0.128` is NOT smearing.  In measured-first
sampling the generator draws the observed point from `sigma_obs` directly, so
the rad arm's excess at low xB is the shape of `eta = sigma_obs/sigma_Born`,
which rises monotonically as xB falls.  Measured from the delivered
2 x 1e6-event campaign with no model evaluation at all, via
`etabar_i = [f_rad(i)/f_born(i)] * C_rad` (Launch generator session,
`docs/farm_practice.md` in exclurad_py), with
`C_rad = 1.25010 +- 0.00175` from the chunk scatter of 200 chunks per arm:

    xB bin            etabar
    0.350-1.000        0.99
    0.128-0.180        1.43
    0.100-0.128        1.91
    0.000-0.080        4.75

So below `xB = 0.08` the observed cross-section is nearly five times its own
Born value -- not an artefact an unfolding could undo, but where `sigma_obs`
genuinely sits.  Note also `etabar < 1` at high xB: the virtual plus soft
correction goes negative there and the real-photon tail does not make it up,
so the usual intuition that radiation increases the yield is built on the
low-xB end of this curve only.

A SECOND, INDEPENDENT ESTIMATE WAS WRONG AND IS WITHDRAWN.  An earlier version
of this note carried an "exact" column (0.878, 1.29, 1.63, 3.04) computed by
integrating `rc_factor` over each bin.  It built the model once with
`t_nucl = -0.30` and never updated it, and the model's t dependence lives in
`t_nucl` -- the `t` passed to `rc_factor` is the SHIFTED transfer used inside
the radiative integral (`exclurad_py/models/base.py:318`).  So it integrated a
cross section with NO t dependence over `|t| = 0.001-2.50`.  Measured cost:
`C_rad` comes out 1.143 frozen against 1.294 per point, against the campaign's
1.25010.  A frozen `t_nucl` returns a finite, plausible cross section at every
point and nothing raises, which is why it survived several rounds of checking.
Discard every number from that column; the campaign column above is untouched,
since it comes from two histograms and the generator.

**The conditional.** All of the above holds because the production box starts
at `Q2 = 0.80`, where `xB` is about 0.17 and the slope parameter is 2.33
(against 2.34 measured in the delivered structure functions -- agreement to
half a per cent).  If anyone lowers the Q2 floor below 0.80 to chase outbending
acceptance, the unbounded slope parameter stops being harmless, and nothing in
the code or the cards would say so.
