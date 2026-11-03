# Runaway solutions in the Mikulášek-model fit: diagnosis and bounding advice

Notes prepared 3 November 2026 · Review of the runaway\_gamma\_experiment note (1800 s spot curves, L000 and L180)

## Verdict

Bounding the model parameters is a sound safety net, but the runaway is a symptom of two other things. The first, below, I checked in the numbers of the note and in a clean-data test. The second is a hypothesis that I could not test without the data.

## 1. What the runaway is

In all four failed solutions reported in the note, A·Γ comes out at 0.075–0.113 mag: 0.091 and 0.075 for the two runaways, 0.113 and 0.076 for the two grid-check solutions. That is the signature of the Γ → 0 limit, where the profile tends to AΓ·(−ln(1 − exp(1 − cosh x))). I confirmed this limit numerically.

- **The needle.** It exists only within about |x| ≲ 10⁻⁶ of the eclipse centre, which is why the fine grid shows a peak of 1811 mag. At the sampled phases the primary is a cusp of roughly 0.08 mag.
- **The primary has been removed.** The real eclipse is about 0.5 mag deep, so the model has effectively dropped it. This explains the “modest” 24 mmag rms: only a few percent of the points are in eclipse. Sample rms is not a safe quality measure here.
- **Why A explodes.** A is a linear coefficient. As Γ → 0 its basis function shrinks towards zero, so A grows without limit.

## 2. Why it is not simple unidentifiability

I tested a clean case: an eclipse narrower than the exposure (1800 s cadence, 21 sub-samples, no spots, 3 mmag noise, FWHM 0.013 against an exposure of 0.018 in phase). The data still reject small Γ clearly: χ² is 1273 at the true Γ = 1, 2535 at Γ = 0.1 and 5325 at Γ = 0.03, for 1296 points. So on clean data the shape is constrained even when it is narrower than an exposure.

**Hypothesis for the spot runs (not verified).** The model is misspecified in two linked ways:

- The staged model had only H2 and the O'Connell term. A spot at longitude 0° or 180° projects mainly onto cos(2πϑ), which is H1. I could not check whether H1 was free in the spot runs.
- The L000 bounded solution looks like the model absorbing that wave: its primary component sits on D = 0.25 with A = 0.1 mag. At D = 0.25 and Γ = 2.21 the FWHM is 0.74 of the orbit, which is impossible for an eclipse (the maximum FWHM is about 0.5).

**First test:** free H1 in the spot runs, and H3 as a second step.

## 3. Revised bounds

Bound observables and the regular domain of the function, not the “physical parameters”. D, Γ and C are phenomenological, and the paper says the method gives no direct geometry.

| Parameter | Interval in the note | Suggestion | Reason |
| --- | --- | --- | --- |
| Γ | 0.3–6 | Keep the floor 0.3. Raise the cap to about 30–50 and fit in log Γ. | The floor is justified by the A·Γ valley (and by Fig. 1 and Table 1 of the paper: 0.62–1.7). A cap of 6 is too low for long flat bottoms: in my ideal-geometry tests Γ was about 10.8 at flat/contact width 0.58 and about 229 at 0.79. |
| D | 0.005–0.25 | Bound the FWHM = 2D·arccosh(1 − ln(1 − 2^(−1/Γ))) to about 0.45–0.5 of the orbit or less. Lower limit tied to the exposure: FWHM ≳ 0.5 × exposure phase, otherwise flag “unresolved” and fix the shape. | D means different things for different Γ. The FWHM bound rejects the D = 0.25 solution. |
| C | −2 to 2 | About −0.6 to 0.3, default 0, freed only if significant (about 3σ). | The paper's Table 1 values are −0.34 to +0.10, and it calls C mostly negligible. At C = −2 the factor (1 + Cx²) goes negative. |
| A | 0.0001–3 mag | Upper limit about 1.5 × the folded amplitude. Drop the secondary if A₂ is below about 5σ and tie D₂, Γ₂ to eclipse 1. | At A₂ → 0, D₂ and Γ₂ are undetermined. |
| A\_b, D\_b, Γ\_b | ±3 mag, as above | A\_b ≥ 0 by default, D\_b ≤ D, A + A\_b ≤ A\_max. | The ordering breaks the exchange symmetry between the two identical components. |
| Harmonics | unbounded | k ≤ 3 (the paper's practice) plus Nyquist. | H, S and A are linear, so a runaway there means an ill-conditioned linear problem. Postponing the amplitude ceiling is fine. |

## 4. Implementation advice

- **Use bounded optimisation.** `least_squares(method="lm")` cannot take bounds, and “hold at the near edge and refit” is a crude active-set method. Use `trf` or `dogbox` with bounds, with D and Γ in logs.
- **Make A linear.** Solve A, A\_b, H, S and the baselines in the inner linear problem with box constraints, for example non-negative least squares. A then cannot run away, and the nonlinear search shrinks.
- **A bound hit means failure.** Label the run “bound-active”, mark the covariance of those parameters invalid, and do not time silently. The L000 solution, with D on its upper edge, is not a success.
- **Replace the grid check.** Compare the model with the binned folded data in eclipse, which catches “primary removed”. Add a quadrature-convergence test (results at n\_sub and 4·n\_sub), which catches the exposure comb. Neither compares the template with itself.
- **Stage C early.** The collapse appeared when C was fixed at 0 and the step refitted. Set C = 0 from the start, as the paper suggests, and restart refits from the last good solution.

## 5. Revising the experiment note

- Reduce the A1 > 5 mag and Γ < 0.02 rule and the grid check to one “tried and dropped” line each.
- State explicitly that the Γ ≥ 0.3 floor is a design choice. The note currently says a floor “is not a rule” and then adopts one.
- Add the A·Γ finding and the run status labels (ok, bound-active, unresolved).

## Open question and next step

Was H1 free in the spot runs? If not, try that before anything else. Afterwards, implement the bounded `trf` fit with the inner linear solve.
