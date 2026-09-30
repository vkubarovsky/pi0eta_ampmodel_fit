**Draft — not sent.**  To: P. Kroll.  From: V. Kubarovsky.  2026-09-30.
Subject: constraint on sigma_L/sigma_T as Q^2 -> 0 for pi0 electroproduction

---

Dear Peter,

thank you again for settling the factor of 2 in dsigma_TT so clearly — we have
recorded it and corrected our notes.

I have another question, at the opposite end of the Q^2 range.

We have been fitting an amplitude-level model to the CLAS, Hall A and COMPASS
pi0 and eta data.  It is built on your framework: flavour form factors for
H_T, Ebar_T, H~ and E~, combined into the channel amplitudes, with the Q^2
dependence of each block left free as a power, N (Q^2)^{n/2}, rather than taken
from a convolution.  We recently found and repaired a defect: the model had no
H~ at all, so the longitudinal non-flip amplitude carried the sqrt(-t') of the
flip amplitude and sigma_L vanished at the forward peak.  With H~ restored, the
fit improves by 74 units of chi^2 on about 1500 points.

That repair has made us look at the Q^2 dependence of R = sigma_L/sigma_T, and
we find we do not know what the theory allows at small Q^2.

**The question.**  At fixed W and t, as Q^2 -> 0, is R constrained?  Our
reasoning is the elementary one: the longitudinal polarisation vector
eps_L = (|q|, 0, 0, q_0)/sqrt(Q^2) becomes proportional to q itself at the real
photon point, so current conservation gives

    M_L = J_z sqrt(Q^2)/q_0 ,   hence   sigma_L ~ Q^2 ,   R ~ Q^2 -> 0 ,

with sigma_T going to the (finite) photoproduction cross section.  Equivalently,
sigma_U = sigma_T + eps sigma_L must reduce to the real-photon cross section,
and eps -> 0.976 rather than 0 in our kinematics, so sigma_L has to vanish.  Is
that the whole story, or is there something specific to pseudoscalar production
that modifies it?

**Why it matters to us.**  The leading-twist amplitude is ~ 1/Q, so the formula
we are using gives sigma_L ~ 1/Q^2, which diverges at the photon point.  Our
fitted model therefore has R RISING as Q^2 falls: at W = 2.2 GeV, t = -0.3
GeV^2 we get R = 0.17 at Q^2 = 1, 0.78 at Q^2 = 0.1 and 13.8 at Q^2 = 0.01 —
that is, sigma_L fourteen times sigma_T for an almost real photon.  This is of
course extrapolation far outside where the formula is meant to work, but it is
not academic for us: our radiative-correction code evaluates the structure
functions at the shifted Q^2 of the hadronic vertex, which reaches well below
the fitted range.

**What we would like to do, and where we need guidance.**  The obvious repair is
to multiply the longitudinal form factors by

    [ Q^2 / (Q^2 + m^2) ]^p ,    p = 3/2 - n/2

where n is the fitted Q^2 power of that block, so that the sqrt(Q^2) behaviour
of the amplitude is exact by construction and the factor goes to 1 at large
Q^2.  Three questions:

1. Is this the right shape, or is there a standard treatment in the GPD
   framework that we should be using instead?

2. What sets m?  Our instinct is m = m_rho, on the grounds that the turn-on of
   the longitudinal coupling is governed by the hadronic structure of the photon
   and its lowest pole — the same reasoning behind R ~ Q^2/M_V^2 for vector
   mesons.  We considered mu_pi = m_pi^2/(m_u+m_d) ~ 2 GeV, but that seems to us
   the wrong scale: it normalises the twist-3 DA and belongs in sigma_T, not in
   the L/T turn-on.  Would you agree?  The choice is not academic — at Q^2 = 1
   the factor is 0.58 with m_rho and 0.16 with mu_pi.

3. Should sigma_LT and sigma_LT' carry one power of the same factor, as the
   L x T interference structure suggests?

We are aware that our data cannot answer any of this: there is no sigma_L
measurement anywhere in our set, and nothing below Q^2 = 1.  Whatever we adopt
will be a theory prior, so we would rather adopt one you consider defensible
than one we invented.

One last observation, which we distrust precisely because it is convenient.  The
same damping factor also removes a discrepancy at the OTHER end: our fitted R
falls with Q^2 at fixed x_B (d ln R / d ln Q^2 = -0.54), whereas twist counting
and your own 2024 curves have it rising.  With m = m_rho the slope turns to
+0.51 at Q^2 = 1, and with mu_pi it lands at +1.5, essentially the twist-counting
value.  We are reluctant to read that as evidence, since the fit will simply
rescale H~ to compensate and chi^2 will barely move.  But if you think the
agreement is telling us something rather than nothing, that would change how we
proceed.

A related puzzle we cannot resolve, in case it is familiar to you: the measured
Q^2 dependence of sigma_T is much flatter than twist-3 counting.  Taking Hall A
2021 alone (x_B ~ 0.46, Q^2 from 2.67 to 6.56) we get

    d ln sigma_U / d ln Q^2 = -3.13 +- 0.20, -2.98 +- 0.17, -2.96 +- 0.21

at three t' values, i.e. sigma_T ~ Q^-6 rather than Q^-8, about six sigma away.
Your 2024 curves fall as Q^-7.2 over the same interval.  Is this understood?

With best regards, and thank you for your patience with our questions,

Valery
