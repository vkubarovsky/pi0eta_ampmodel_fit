# Installing the H~ model in the generator — plan, and the trap

**Paused 2026-10-02 at VPK's request** (he is on the speed work).  Nothing in
`~/exclurad_py` has been touched.  `amp2609` is byte-identical to what generated
the 89.8 M September events, and those remain valid.

## Read this first: it is a code change, not a file copy

`exclurad_py` does **not** import the fit's amplitude code.  It keeps its own
copies — three of them — and none knows about the H~ block (slots 37-42, merged
into `pi0eta_ampmodel_fit` master as `2fcc412`).

| file | slots | used by |
|---|---|---|
| `_amplitude_fit.py` | 27, old slope convention | amp2021, amp2026, amp2026s |
| `_amplitude_fit2.py` | 37, `b + b' ln xB` | scalar fallback for amp2026c/n, amp2609 |
| `_amplitude_fit2_vec.py` | 37, vectorised | **fast path for amp2609** |

Dropping a 43-slot fitpar into `exclurad_py/models/` raises **no error** and
gives neither the old model nor the new one.  Measured at `pi0 p`, `t = -0.343`,
`xB = 0.2`, `Q2 = 2.5`:

| | generator would give | correct | factor |
|---|---|---|---|
| sigma_L | 28.68 | 5.75 | **x5.0** |
| sigma_LT | -22.17 | -2.09 | **x10.6** |
| sigma_T | 73.18 | 73.18 | x1.0 |
| sigma_TT | -39.06 | -39.06 | x1.0 |

Why: the new vector has `p[25] = rho_nf = 0`, so the old code computes
`T00p = rho_nf*|T00m| = 0` and sigma_L collapses to `|T00m|^2` alone; and the
refitted `<E~>` normalisation (slot 15: 35.8 -> 127.8) is calibrated to work
*alongside* `<H~>`, so by itself it is far too large.

With a 37-slot vector the copies agree with the fit repo exactly (difference 0),
which is why amp2609 is safe.

## The plan

1. **Port the H~ block into both modules amp2609 uses** — `_amplitude_fit2.py`
   (scalar) and `_amplitude_fit2_vec.py` (vectorised).  Same length gate,
   `len(p) > 41`, so every 37-slot model is untouched.

2. **Add a permanent closure test.**  Evaluate `pi0eta_ampmodel_fit/amplitudes.py`
   and both generator copies over a grid; require agreement to machine precision
   for BOTH the 37- and the 43-slot vector.  This is the test that would have
   caught this class of bug, and it is worth having regardless of H~.

3. **Install under a NEW tag** (`amp2610` or similar), leaving `amp2609` in
   place, so the September events stay reproducible.  Provenance string should
   name the merge commit `2fcc412`.

4. **Validate before generating a single event:**
   - positivity over the campaign box, four channels, every (h, lambda) —
     especially `pi0 n`, where the margin fell from 0.0055 to 0.0009;
   - Born closure (the generator's own check);
   - the RC factor eta, old model vs new, at box points.  This also closes the
     question left open on 2026-09-30: at `Q2 = 1`, 45 % of the integrand
     evaluations sit below `Q2 = 1` and reach `Q2t = 0.18`, where the two models
     differ by a factor 3 in R.  Above `Q2 = 2` the exposure is a few parts in
     a thousand (`rc_q2t_exposure.py`).
   - a pilot sample, ~100 k events, compared directly against the Born cross
     section.

5. **Only then the campaign.**

## Deliberately NOT in this plan

The duplication itself.  Three copies of one physics in two repositories is how
this happens, and the right fix is for `exclurad_py` to import the fit module or
vendor it under a checksum test.  But that is a refactor of production code and
should not ride along with a model install.  Separate task, afterwards.

## Still open on the physics side

- Positivity margin in `pi0 n` at full polarisation: 0.0055 -> 0.0009.
- The model has no photoproduction limit.  A factor `[Q2/(Q2+m^2)]^{3/2-nQ/2}`
  would impose it; the scale `m` is unconstrained by our data (nothing below
  `Q2 = 1`, and no sigma_L measurement anywhere).  Letter to Kroll asking which
  scale is defensible is drafted in `letters/` and **not sent**.
- H~ does not resolve the sigma_LT sign disagreement between CLAS12 and Hall A.
- `tex/report.tex` is from 18 September and mentions H~ zero times.
