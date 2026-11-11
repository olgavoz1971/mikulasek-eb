# Run-away fits: what happened, what we do, and what we do not know

This note records an experience, not a theory. The limits below are design choices of this package. The paper (Mikulášek 2015, Sect. 5.1) fits with Newton-Raphson and sets no bounds. The numbers come from synthetic spotted binaries and a few real stars, so they are orientation for the same kind of data, not a calibration.

## 1. What a run-away is

A least-squares run can return `success` with a template that is not an eclipse.

On two 1800 s light curves of a synthetic detached binary with a spot, the primary became a needle. The depth reached about 38000 mag and 5400 mag, and Gamma fell to about 1e-6 and 1e-5. The rms on the samples stayed near 24 mmag.

Two things explain it.

- As Gamma goes to zero the product of depth and Gamma stays near 0.08 to 0.11 mag. The profile tends to a cusp that is nonzero only within about 1e-6 of the eclipse centre. The depth A is a linear coefficient, so as its basis function shrinks, A grows without limit.
- The needle is narrower than the gap between samples, and narrower than the exposure. It falls between the points, so the sample rms cannot see it. The real eclipse was about 0.5 mag deep and the model had dropped it. Only a few percent of the points were in eclipse, which is why 24 mmag looked modest.

It is not simple unidentifiability. On a clean test with an eclipse narrower than the exposure, the data still rejected small Gamma clearly (chi-square 1273 at the true Gamma of 1, 2535 at 0.1, 5325 at 0.03). The run-away came with a misspecified model: on the spot curves the harmonic that carries a spot at longitude 0 or 180 degrees was missing, and the fit absorbed that wave with a primary component of width 0.25 and depth 0.1 mag.

## 2. What was tried and dropped

| Attempt | What it did | Why it was dropped |
|---|---|---|
| Cuts on depth (above 5 mag) and Gamma (below 0.02), then refit with Gamma held | Removed the needle on both curves | The numbers were a first trial. A floor on Gamma was not justified at the time. The refit still left a thin spike |
| Compare the template with itself on a fine grid, reject a step if the grid range exceeds three times the sample range | Stopped the catastrophic needle | It accepts sharp eclipses that the samples happen to land on. A narrow spike with a small comb of bumps passed. It compares the template with itself, not with the data |
| Hold a parameter on the edge of a box and refit, each step | Removed the needle on both curves and kept a second component | One result sat on the upper width edge (width 0.25), with the eclipse carried by a second component. The edge was doing the work |
| Bounded trust-region solver as the default | Gave `ok` fits on the 1800 s curve at spot longitude 180 | On XY Boo, a depth and a cosine cancelled, the depth ran up to its cap, and every candidate ended `bound-active`. The selection chose nothing |

## 3. What the package does now

### 3.1 A verdict after every fit

The default solver is Levenberg-Marquardt, with no bounds while it runs. The finished template is then judged by `assess_fit` (`mikulasek_eb.validity`):

| Check | Limit | Status when violated |
|---|---|---|
| Gamma of every profile copy | 0.3 to 30 | `bound-active` |
| Depth of every copy | between 0 and twice the robust range of the data (the second component may be negative, same size limit) | `bound-active` |
| Full width at half depth of every copy, measured on a fine grid | at most 0.5 cycles, at least half the exposure in phase | `bound-active` |
| Template against the data inside each eclipse core | mean residual below 0.2 of the folded depth, or within 5 standard errors | `unresolved` |
| Exposure integral, `nsub` against `4 nsub` | difference on a fine phase grid below the noise | `unresolved` |
| Curvature C | outside -0.6 to 0.3 | warning only |
| Jacobian rank | deficient | warning only |

Why these limits.

- Gamma from 0.3 is a design choice. The values in the paper's Table 1 are 0.62 to 1.7, and its Fig. 1 draws 1/2. The cap of 30 is high on purpose: a long flat bottom needs a large Gamma (on an ideal synthetic occultation the best Gamma was about 10.8 at a flat-to-contact width ratio of 0.58 and about 229 at 0.78).
- The width limit is on the observable (full width at half depth), not on D. D means different things for different Gamma, and a bound on D alone let a width-0.25 component through.
- The core check looks at the data, so a primary that has been removed cannot hide behind a healthy secondary. The exposure check finds the comb that a spike narrower than the sub-step leaves.
- A value within a small tolerance of an edge counts as on the edge.

Only `ok` is compared. Only `ok` is timed: `time` raises `UnreliableSolution` for any other status, and for a solution that has no status. Harmonics and baseline coefficients are not bounded. `ok` means numerically sound, not a good fit. A warning says so when rms/noise is above 1.5.

On the 27 files of the spot test (three cadences, nine spot settings) the two needles were `bound-active` and the other 25 were `ok`. The four shoulder longitudes were `ok` with an rms of about 22 mmag against a noise of about 11.5 mmag: numerically sound and poorly fitted.

### 3.2 A second solver, only as a retry

`select_model` fits a candidate with Levenberg-Marquardt. If the status is `bound-active` and `retry_bound_active` is true (the default), the same candidate is fitted once more with the trust-region solver. Depth, width and Gamma are then kept inside the boxes, with D and Gamma fitted as logarithms. The verdict of the second fit is the one compared. A retry that raises is logged and the bound-active result is kept.

What the retry did in the tests:

- On the nine 1800 s spot curves, seven candidates were retried and two became `ok`. No choice changed, because another Levenberg-Marquardt candidate was already `ok` and won.
- On a synthetic semi-detached curve (G7, mixed trend), both profile-1 candidates were `bound-active` because the secondary width was negative. With the retry off the search had no winner to carry forward. With it on, they finished at about 14.6 mmag and lost to profile 3, which was fitted under Levenberg-Marquardt.
- On a synthetic 1800 s curve (G6), one `bound-active` candidate became `ok` after the retry, at 7.03 mmag, and was chosen. This is the only chosen fit in that batch that the retry changed.

The cost is one more full staged fit for each `bound-active` candidate. A light curve on which several candidates hit a bound pays that for each of them.

A candidate that raises `FitDidNotConverge` is a different failure. It never receives a verdict, so the retry never starts. `select_model(..., retry_not_converged=True)` fits it again with the trust-region solver. It is off by default. On one synthetic file (G5) it let the (1, 2, 3) candidates finish and the evidence chose them, at a lower rms, but it is slower still because the first fit has already used its evaluation limit.

Starting with the trust-region solver is off. See the last row of the table in section 2.

### 3.3 The named three-step search turns the retry off

`tiered_search` builds every candidate with the retry off. That was copied from the synthetic tiered runs and is a choice, not a finding. With the retry off, a step in which every candidate is `bound-active` stops the search. That is what happened on a real star (SV Cen). See the notes on `start_profile` in the `mikulasek-eb` skill.

## 4. What we do not know

- No test file had every Levenberg-Marquardt candidate `bound-active` and then recovered by the retry. The retry is shown to help a search continue, not to rescue a file.
- A harmonic of high order and large amplitude, one the samples cannot resolve, is an open case. The paper's practice is a few low-order terms. No ceiling on the harmonics has been tried, and none is claimed.
- The full-width limit is judged after the fit. It is not a constraint inside either solver.
- A linear inner solution for the amplitudes, the harmonics and the baselines, so that those cannot run away, is postponed.
- Everything here is one synthetic binary with one spot, a handful of other synthetic systems, and a few real stars. It does not generalise to other stars or other surveys.

Details and the full table of checks: `mikulasek_implementation_notes.md`, section 11. The review of the first experiment note, which proposed the revised bounds and the order "never report rubbish first, repair later", is in `Runaway_solutions_in_the_Mikulasek-model_fit_diagnosis_and_bounding_advice.md`.
