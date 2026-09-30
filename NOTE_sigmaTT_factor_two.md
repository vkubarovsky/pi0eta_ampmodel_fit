# The factor 2 in dσ_TT: a typo, propagated for fifteen years

**Settled 2026-09-30 by P. Kroll, private communication to V. Kubarovsky.**
Quoted in full at the bottom.  This note exists so that nobody re-derives the
question a third time.

## The statement

    |dσ_TT/dt|  ≤  dσ_T/dt

is the correct bound.  Not `≤ ½ dσ_T/dt`.

## Where the wrong factor comes from

| source | what it says | status |
|---|---|---|
| GK 2010, "An attempt…", eq. (43) | dσ_TT missing a factor 2 | **typo** |
| Kroll 4-page note 2011, eq. (7) | same expression | carries the typo |
| Kroll 4-page note 2011, eq. (11) | `dσ_T = −2 dσ_TT + (H_T term)` | the 2 is wrong |
| COMPASS PLB 870 139832 (2025), eq. (8) | `dσ_TT ∝ (t'/16 m²)|⟨Ēᵀ⟩|²` | carries the typo |
| Kroll 2312.13164 (2023), eq. (64) | | **correct** |

The 2010 and 2011 papers showed no TT cross sections, so the typo stayed
invisible in them and went unnoticed for years.  COMPASS 2025 quotes the
formalism from those papers and inherits it.

## Why it matters and why it did not

**It does not touch the fit.**  We fit measured σ_T, σ_TT, σ_LT directly; the
internal amplitudes are our own parametrisation and a factor 2 in the relation
between two of them would simply be absorbed into the fitted normalisation of
⟨Ēᵀ⟩.

**It touches two things.**  Comparing our ⟨Ēᵀ⟩ against GK's number, where a
factor 2 of pure convention would look like physics.  And any bound argument:
the wrong version forbids `|σ_TT| > σ_T/2`, which a good deal of real data
exceeds.

## Our dictionary was right

With the definitions of the amplitude note,

    T01^(ν') = ½(M_{0ν',++} − M_{0ν',−+})      natural
    U01^(ν') = ½(M_{0ν',++} + M_{0ν',−+})      unnatural

the Gram dictionary gives

    σ_T  = Σ|T01|² + Σ|U01|²
    σ_TT = Σ|T01|² − Σ|U01|²

so that

    σ_T + σ_TT = 2 Σ|T01|²  ≥ 0
    σ_T − σ_TT = 2 Σ|U01|²  ≥ 0

— the natural and unnatural cross sections, each manifestly non-negative, and
`|σ_TT| ≤ σ_T` follows identically.  That is the bound Kroll confirms.  The
correct bound is the one our dictionary produces without being told.

## A wrong turn worth recording

On 2026-09-29, reasoning from the (typo-carrying) COMPASS eq. (8), this project
concluded that the data "significantly violate" the bound `|σ_TT| ≤ σ_T/2` —
60 CLAS12 points, 7 CLAS6 points, 4 of COMPASS's own 13 — and took that as
evidence that the two-GPD leading-twist picture was incomplete.

Then, reading Kroll's note eq. (17) ("−2dσ_TT almost saturates the unseparated
cross section"), the conclusion was reversed: the median |σ_TT|/σ_U = 0.45 was
read as agreeing with him, and the violation claim retracted.

Both readings were built on the typo.  The data were never in tension with
anything; they sit comfortably inside `|σ_TT| ≤ σ_T`.  The lesson is narrow and
practical: a bound quoted from a formula is only as good as the formula, and
when data appear to violate a published inequality, suspect the inequality
before the data — especially when the same inequality follows for free from a
positivity argument that gives a different answer.

## The email

> Dear Valery,
>
> unfortunately there is a typo in Eq. (43) of the GK paper from 2010
> ("attempt…."): a factor of 2 is lacking in the expression for dsigma_TT (or
> this factor of 2 has to appear in Eq (42) in front of the dsigma_TT term).
> Since in the 2010 paper as well as in the next one from 2011 ("Transversity
> in hard exclusive electroproduction of pseudoscalar mesons") no results on
> the TT cross sections are shown this typo has no effect on these papers and
> we have not noticed this typo for quite a while.  Thus, in
> kroll_note_4_page_2011.pdf which is probably the attached note, this factor
> of 2 is lacking in Eq (7) too.  Consequently the factor of 2 in front of
> dsigma_TT in Eq (11) is wrong too.  Correct is
>
>     | d\sigma_TT/dt| ≤ d\sigma_T/dt
>
> and this can easily be shown.  Thus Eq (64) in the recent paper (2312.13164)
> is correct.
>
> Best, Peter
