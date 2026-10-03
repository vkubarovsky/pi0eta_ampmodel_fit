# Merging the H~ model into `exclurad_py` production — plan

Written in the shape of `exclurad_py/production/speed/VALIDATION_PLAN.md`,
because the same split applies and the same ordering rule: **the cheapest thing
that could falsify the work goes first, and nothing expensive runs until the
cheap checks pass.**

    PART A   every EXISTING model          claim: the output is UNCHANGED
    PART B   the new tag                   claim: the output is what we fitted

Part A can be proved — byte-identity.  Part B cannot be proved against anything,
because no sample with this model exists yet; it can only be checked for
internal consistency and measured against the old one.

Source: `pi0eta_ampmodel_fit` master `2fcc412` (H~ merged), fit `Htil_free2`.

---

## Stage 0 — the ground is favourable.  Checked 2026-10-03.

The speed work and the model work do not overlap:

| workstream | files |
|---|---|
| speed | `core/integrand.py`, `core/_kernel_numba.py`, `run.py` |
| model | `models/_amplitude_fit2*.py`, `models/pseudoscalar.py`, a new `.npy` |

No file under `models/` appears in the last 20 commits of `production`, and
`_kernel_numba.py` says in its own docstring that *"the structure-function model
stays outside -- it is a Python object"*.  So a merge should be conflict-free,
and the model change does not interact with the numba kernel.

**Re-check this at merge time.**  `production` is moving; if `models/` has been
touched by then, stop and re-read.

---

## Stage 1 — port H~ into the two modules amp2609 actually uses.  ~2 h.

`exclurad_py` keeps its own copies of the amplitude code.  `amp2609` uses
**`_amplitude_fit2_vec.py`** (fast path) with **`_amplitude_fit2.py`** (scalar)
as fallback.  Neither knows the H~ block.

Port `_htilde()` and the corrected longitudinal sector into both, gated on
`len(p) > 41` exactly as in the fit repo, so every 27- and 37-slot model takes
the old path untouched.  `_amplitude_fit.py` (27 slots, old convention) is not
used by amp2609 and is **not** to be touched.

---

## Stage 2 — the closure test.  Seconds.  The cheapest falsifier.

A test that evaluates, over a kinematic grid covering the campaign box and all
four channels:

    pi0eta_ampmodel_fit/amplitudes.py      (the reference)
    exclurad_py _amplitude_fit2.py         (scalar)
    exclurad_py _amplitude_fit2_vec.py     (vectorised)

and requires all three to agree to machine precision, for **both** the 37-slot
and the 43-slot vector, on all nine structure functions.

Two axes matter and both must be checked: reference-vs-generator (catches a
botched port) and scalar-vs-vectorised (catches a port done in only one of
them).  This test is worth having whatever happens to H~ — its absence is why
the duplication was invisible.

**Acceptance: exact agreement.  Anything else stops everything below.**

---

## Stage 3 — byte-identity for every existing model.  ~1 h.

Part A's whole claim is that there is nothing to discuss.  Generate short
samples before and after the Stage-1 change, same seed, and `cmp` the LUND:

| model | why |
|---|---|
| `pi0.amp2609`, `eta.amp2609`, `pi0n.amp2609`, `etan.amp2609` | all four registered channels |
| one of `amp2026c` / `amp2026n` | also 37 slots, same modules |
| one of `amp2021` / `amp2026s` | 27 slots, `_amplitude_fit.py`, must be inert |
| born and rad | different code paths |
| `ptarg != 0` | the POL_KEYS branch |

**Acceptance: byte-identical.  Anything else stops everything below** — the
87.9 M events already generated with amp2609 depend on this.

---

## Stage 4 — install the new tag.  ~1 h.

New `.npy` beside `amp2609_par.npy`, a new `_ampXXXX_sfs` and four registered
classes.  **`amp2609` is not modified and not removed**, so the September
samples stay reproducible.

Provenance string must name the fit commit `2fcc412` and the run `Htil_free2`,
and the validity ranges must be re-derived — the fitted Q2 and |t| coverage did
not change, but the string is copied by hand and has been wrong before.

Guard: a one-line clip of the phi weight at zero in `generator/events.py`.  The
bracket can legitimately touch zero at isolated points under full polarisation
(measured: +9e-09 in pi0 n), and in floating point that can land just below.

---

## Stage 5 — what the new model does, measured not proved.  ~3 h.

| check | expectation |
|---|---|
| positivity over the box, four channels, every (h, lam) | B >= 0.  Established structurally: 80 000 random amplitude draws, zero violations |
| Born closure | the generator's own existing check |
| RC factor `eta`, old tag vs new, at box points | unknown.  At `Q2 = 1`, 45 % of integrand evaluations sit below `Q2 = 1` and reach `Q2t = 0.18`, where the two models differ by a factor 3 in R.  Above `Q2 = 2` the exposure is parts in a thousand |
| generated spectra vs the Born cross section | agreement |

Expected changes, already measured on the model side (`generator_impact.py`):
`sigma_U` moves ~1 % for pi0 p, 9-16 % for eta p, up to a factor 2 for pi0 n.

---

## Stage 6 — pilot, then campaign.

~100 k events per channel first, compared against Stage 5's predictions, before
anything at production scale.

---

## Not in this plan, deliberately

**The duplication itself.**  Three copies of one physics in two repositories is
why this plan needs Stages 1-3 at all.  The fix is for `exclurad_py` to import
the fit module, or to vendor it under a checksum test.  That is a refactor of
production code and must not ride along with a model install.  Afterwards.

**The photoproduction factor** `[Q2/(Q2+m^2)]^{3/2-nQ/2}`.  The model has no
Q^2 -> 0 limit, but the RC never goes below `Q2t = 0.18`, so this does not
affect our numbers.  The scale `m` is unconstrained by our data; the letter to
Kroll asking which scale is defensible is drafted in `letters/` and **not sent**.
If it is added later it changes the model and reopens Stages 4-6.

---

## The trap, restated, because it is silent

Dropping a 43-slot fitpar into `exclurad_py/models/` **without Stage 1** raises
no error and gives neither the old model nor the new one.  Measured at
`pi0 p`, `t = -0.343`, `xB = 0.2`, `Q2 = 2.5`:

| | would give | correct | factor |
|---|---|---|---|
| sigma_L | 28.68 | 5.75 | **x5.0** |
| sigma_LT | -22.17 | -2.09 | **x10.6** |
| sigma_T, sigma_TT | correct | | x1.0 |

Because `p[25] = rho_nf` is now 0, so the old code sets `T00p = 0` and sigma_L
collapses to `|T00m|^2`; and the refitted `<E~>` normalisation (slot 15:
35.8 -> 127.8) is calibrated to work *alongside* `<H~>`.
