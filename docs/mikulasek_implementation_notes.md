# Mikulášek (2015) phenomenological model for eclipsing binaries: implementation notes

**Scope.** Only eclipsing binaries (EBs) with circular orbits, single-band light curves in magnitudes, and only the method of Mikulášek (2015), A&A 584, A8. Code: the package `mikulasek_eb` (`src/mikulasek_eb/`), `examples/compare_timings.py` and the agent skill `mikulasek-eb-fitting`.

**This is an independent re-implementation from the published equations, not the author's code.**

**Legend.** ✅ implemented as in the paper · ⚠️ implemented with a deviation or simplification · ➕ added by me, not in the paper · ❌ not implemented.

**Validation status.** Synthetic tests with known truth, plus one real sector (AA And, TIC 455790537, SPOC sector 16, 2-min cadence) fitted by the user. Since 2026-10-29 the four one-star columns of Table 1 (EK Com, V477 Lyr, AV Del, AR Aur) have been compared with TESS light curves of the same stars (section 7). There is still **no** comparison with XY Boo (Tables 3 to 6), with the `ebfit` program cited in the paper (see section 4), or with MAVKA.

---

## 1. Elements of the paper that are implemented

| Paper element | Code | Status | Notes |
|---|---|---|---|
| LC model F = F<sub>e</sub> + F<sub>p</sub> + F<sub>c</sub> in magnitudes; instrumental term Y<sub>0r</sub> (Eqs. 10, 12) | `EBModel`, baselines in `fit_global` | ✅ | The intrinsic F has zero mean without eclipses; the level is absorbed by the baseline. |
| Eclipse profile, Eq. 13 = Table 2 #1: A(1 + Cφ²/D²){1 − [1 − exp(1 − cosh(φ/D))]^Γ}, eclipses at phase 0 and 0.5 | `profile=1` | ✅ | A, D, Γ, C are per eclipse in the paper. By default I share D, Γ, C between the two eclipses (the paper's usual economy: D₁≅D₂, Γ₁≅Γ₂, C₁≅C₂). Separate values are set by adding `D2, G2, C2` to the parameter dictionary. The paper itself uses different C₁, C₂ for EK Com. |
| Alternative eclipse functions, Table 2 #2–#7 | `profile=2…7` | ✅ | Formulas checked against the printed Table 2 (2026-10-30). #4 is the Andronov-type function A·Re[(1 − \|φ/D\|^Γ)^{3/2}]. #5 Gaussian, #6 triangle of half-width D, #7 exp(1 − cosh) take only A, D (no Γ, C). The staged pipeline (`run_staged_fit`) supports profiles 1 to 4 only, because it frees Γ; profiles 5 to 7 are fitted with `fit_global`. |
| Eq. 14, quartic term K (for very precise transit data) | `EBModel(profile=1, kterm=True)` | ✅ | Keys `K1` (and optional `K2`); off by default, as in the paper where K is introduced only when needed. Recovery tested on a synthetic transit. Not yet handled by `FitConfig` or the pipeline. |
| Proximity cosines and O'Connell sine, Eq. 15 | `harmonics=(…)`, `ocon=True` | ✅ | Notation differs: `H<k>` is the amplitude of cos(2πkϑ) (the paper's A<sub>k</sub>), `S1` is sin(2πϑ). Any subset of harmonics is allowed, e.g. the XY Boo model (Eq. 18) uses cos 4πϑ and cos 6πϑ only: `harmonics=(2,3)`. |
| Phase function ϑ = ϑ₁ − Ṗϑ₁²/2 and Θ(E) = M₀ + P₀E + ½P₀ṖE² (Eqs. 7, 8, 17) | `phase`, `predicted_time` | ✅ | Only the Ṗ term. The paper (footnote 6) quotes the deviation of this approximation from the exact solution for XY Boo: < 1.9×10⁻⁵ in phase and < 0.3 s in time. I have not checked it for larger Ṗ·span. |
| Instrumental term as a linear combination of functions of time (Eq. 11) | `seg`, `order` | ⚠️ | Polynomials of normalised time per segment only. Not implemented: Legendre or quasi-orthogonal polynomial-harmonic functions; weighted-variance normalisation. The coefficients are eliminated by linear least squares inside every iteration (variable projection) instead of being solved jointly. Point estimates should agree. Parameter uncertainties come from the projected Jacobian (approximately the marginal covariance). I have not checked this against a joint solution. |
| χ² minimisation and parameter uncertainties (Sect. 5.1) | `fit_global`, `correlation_matrix`, `strong_correlations`, `curve_band`, `ephemeris_at_epoch` | ⚠️ | Simultaneous weighted χ² minimisation of all template, ephemeris and baseline parameters ✅. Optimiser: Levenberg–Marquardt (SciPy) by default, instead of Newton–Raphson ⚠️. An optional trust-region solver with box bounds is `FitPolicy.bounded_solver` (section 11). After every fit a verdict decides whether the solution may be timed. Covariance (JᵀJ)⁻¹ scaled by reduced χ². **Added in this step:** correlation matrix and degeneracy warnings; propagation of the full covariance (general law of error propagation, "with correlations") to the template curve (`curve_band`); orthogonalisation of the ephemeris parameters by re-expressing M₀, P₀ at the epoch where they are uncorrelated (`ephemeris_at_epoch`). Not done: orthogonalisation including Ṗ; propagation to other derived quantities beyond the template band. |
| Weights (Sect. 5.2, principle 2) | `sigma`, `rescale`, `fit_global_iterative` | ⚠️ | User-supplied σ<sub>i</sub> plus an optional global factor √(reduced χ²), and now the paper's iterative subset-wise estimate of σ from the residual scatter. Caveat: a model misfit concentrated in a subset is absorbed into that subset's σ. |
| Parameter economy (Sect. 5.2, principles 3–4) | `free` lists, profile choice | ✅ | The rule "keep a term only if its amplitude exceeds about 3σ and χ² drops clearly" is my formulation of principle 3. |
| Virtual minima from per-interval O–C shifts (Eq. 9) | `fit_minimum`, `virtual_minima` (fixed template); `fit_global_shifts` (simultaneous) | ⚠️ | See below. |

**Timings: how close to the paper.**
- **Paper (Eq. 9, Sect. 5.2 principle 1, Sect. 6.3).** The phase is shifted per interval, ϑ = ϑ<sub>c</sub> − (O−C)<sub>r</sub>/P. Intervals are typically nights or seasons. The paper recommends processing all data simultaneously, also for individual mid-eclipse times, and propagating uncertainties with the correlations.
- **`fit_global_shifts` (new).** Template parameters and one shift per eclipse are fitted together. The ephemeris is held fixed, because M₀, P₀ are degenerate with the mean and linear trend of the shifts. The timing uncertainties come from the joint covariance, which includes the template uncertainty and its correlation with the shifts; the full matrix `cov_T` is returned.
- **Remaining differences.**
  - The interval is one eclipse (limited to a window ≤ 0.25 in phase), not a night or season.
  - The shifts relative to a fixed ephemeris are not combined with a simultaneous ephemeris fit.
  - I do not know whether the paper's second part of the solution (Sect. 6.3) refits the template; the text can be read either way.
- **Synthetic check.** Simultaneous and fixed-template uncertainties agreed within about 5% (median ratio 0.96), and the largest shift–template correlation was 0.06. The template is determined by thousands of points while each shift uses a few hundred. This was one synthetic setup; with a poorly determined template the difference could be larger. The formal σ of both versions are white-noise estimates only (see section 4).

---

## 2. What I added (not in the paper)

### 2.1 Two-component eclipse profile ➕ (modelling extension)
- **What.** Per eclipse, a second term is added to Eq. 13: A·f(φ; D, Γ) + A<sub>b</sub>·f(φ; D<sub>b</sub>, Γ<sub>b</sub>), with a shared centre, the same function f and C = 0 for the second term.
- **Keys.** Switched on by the keys `A<k>b, D<k>b, G<k>b`. `refine_two_component` tries four starting combinations.
- **Why.** A single Γ cannot make a flat centre and steep walls at once. For a short flat bottom the best single-component Γ stays near 1.2.
- **Evidence (synthetic, ideal geometric occultations).** rms 1.4–1.9 mmag against 2.7–4.4 mmag for one component, on a 0.39 mag deep eclipse.
- **Evidence (AA And S16, single-component → two-component primary).**

| Quantity | Single | Two-component |
|---|---|---|
| rms / out-of-eclipse noise | 2.35 | 1.45 |
| Primary contact ripple | ±12 mmag | ±1.4 mmag |

- **Timings.** Single and two-component timings differ by 1.0 s rms against an O–C scatter of about 16 s. It improves the *shape*, not the timing precision, in that data set.
- **Caveats.**
  - The two components are phenomenological. D and Γ of the individual components must not be interpreted physically.
  - Degeneracy risk: A and A<sub>b</sub> of opposite sign.
  - Tested only on synthetic data and one real sector.
  - The skill requires four acceptance checks before adopting it.
- **The paper's own remedy (found on re-reading).** Sect. 4.1.1: "in EBs with totalities we see that the bottoms of their occultations are flat whilst transits are convex. It can be described by introducing different parameters C1, C2 (see the case of EK Com in Table 1, Fig. 2)". So the first step for a flat bottom is separate C1 and C2, before this extension. The TESS EK Com fit with separate C1 and C2 lowers the rms from 8.5 to 6.2 mmag (section 7) but about 20 mmag residuals stay on the primary walls, so the extension is not made redundant by it. Whether the author has published a further remedy is still open.

### 2.1b Exposure-time integration ➕ (2026-10-30)
- **What.** `EBModel(exposure_days=..., nsub=..., period_days=...)` averages the model over the exposure (a boxcar in time, midpoint rule with `nsub` samples). `FitConfig(exposure_days=..., exposure_nsub=...)` passes it to the pipeline. `mikulasek_eb.exposure` has `estimate_cadence_days`, `estimate_exposure_days(t, duty_cycle)` and `suggest_nsub(exposure, fwhm)` (at least 5 samples, sub-step below a quarter of the FWHM).
- **Why.** Narrow eclipses on a long cadence are smeared; without integration D is biased wide. Synthetic test: only the integrating fit recovers D.
- **Exposure value.** There is no standard annotation for it in VOTable time series, so it is a user parameter or a cadence-based estimate (duty cycle 1 is the upper limit; TESS integrates for less than the cadence). TODO in `exposure.py`: ask for a public `volightcurve` method that returns the exposure time when the data carry it.
- **EK Com (sector 49, 10 min), pipeline fit, C1 = 0:** rms 5.89 mmag without integration, 5.90 with exposure 0.8 x cadence and 5.91 with 1.0 x cadence. The ~20 mmag wall residuals are therefore **not** explained by exposure smearing (hypothesis refuted). D moves from 0.0533 to 0.0523 and Gamma from 1.63 to 1.72 (paper 0.052 and 1.7), which may be a coincidence because the fit quality does not improve.
- **Approximation.** The exposure is converted to phase with a fixed reference period (the fitted period changes by a negligible relative amount); Pdot is ignored in that conversion.

### 2.2 Timing-mode additions ➕
- **Windowed local mode.** Fit only |ϑ − E| < `window` (window 0.25 = whole phase). The motivation was robustness to changing shape.
- **Optional per-event scales** of depth, width and kurtosis (`free=("A","D","G")`). Unstable in windows below 0.1 in my tests.
- **Initial scan over shifts** before the nonlinear fit.
- **Both-walls coverage requirement** (`wall_band`, `min_wall`). Partially covered events rest on one wall and are biased by any wall-shape mismatch.
  - Synthetic test with data gaps: the scatter about the true times fell from 16.5 s to 4.5 s (window 0.25) and from 88 s to 4.5 s (window 0.15).
  - AA And S16: the "all" and "walls" results were identical, so either there were no partial events or the filter was not applied. This is unconfirmed.
- **Cuts and clipping.** A σ cut at 5× the median, and 4σ clipping about a fitted ephemeris.

### 2.3 Fitting aids and diagnostics ➕
- **Fitting aids:**
  - `initial_guess`: M₀, A₁, A₂ and the primary FWHM from the folded curve.
  - `D_from_fwhm`: for profile 1 with C = 0, FWHM = 2D·arccosh(1 − ln(1 − 2^(−1/Γ))).
  - `multistart_global`: refits from several Γ starts.
- **Noise estimate:** `noise_floor` and `out_of_eclipse_mask`. A noise estimate that includes eclipses is inflated: 2.9 against a true 2.0 mmag in my test.
- **Residual inspection:** `residuals` and `binned_residuals`, for inspecting the eclipse and out-of-eclipse residuals.
- **Timing analysis:**
  - `window_scan`.
  - `clip_timings`.
  - `oc_diagnostics`: primary/secondary scatter, the secondary−primary offset, the primary–secondary correlation and the lag-1 autocorrelation. Tested on injected jitter only.
  - `fit_ephemeris`: a plain weighted polynomial ephemeris on the timings.
  - `compare_timings.py`: comparison with another timing list (MAVKA or other).
- **Agent skill.** All thresholds (for example "rms ≲ 1.5× noise") are my heuristics, not the paper's.

### 2.4 The paper's own extra term (K) does not replace 2.1
Eq. 14 adds a quartic factor, A(1 + Cφ²/D² + Kφ⁴/D⁴){…}. It is written for exoplanet transits (n<sub>e</sub> = 1) "in cases of very precise measurements". It is now in the package (`kterm=True`, 2026-10-30); the test below was done earlier with a separate script on ideal geometric flat-bottomed occultations:

| Flat width / contact width | rms, profile 1 (C) | rms, + K | rms, two-component |
|---|---|---|---|
| 0.19 | 4.4 mmag | 3.7 | 1.75 |
| 0.38 | 3.7 | 3.7 | 1.59 |
| 0.58 | 3.3 | 3.0 | 1.44 |
| 0.79 | 2.7 | 2.7 | 1.90 |

K gives little in this idealised test.

---

## 3. Not implemented (for EBs)

| Item in the paper | Why it might matter for you |
|---|---|
| Simultaneous solution of template, **ephemeris** and all O–C shifts together | Done for template + shifts (`fit_global_shifts`); the ephemeris is held fixed there. |
| Orthogonalisation of Ṗ together with M₀ and P₀ (Sect. 5.1) | Only the (M₀, P₀) pair is decorrelated. |
| Higher-order period terms (P̈ …, Eqs. 6–8), the exact Ṗ solution (footnote 6) | Long time spans with strong period change. |
| **Light-time effect** (third body) and **apsidal motion** phase models (Sect. 2.2) | Needed for long-term period studies. Not needed within one TESS sector. |
| **Eccentric orbits** (phase rectification; unpublished in the paper, footnote 3) | See section 3b. The code stops at "secondary not at phase 0.5". |
| **Multicolour model** (Eq. 16, Λ = λ₀/λ<sub>eff</sub> − 1; a<sub>k1</sub> + a<sub>k2</sub>Λ in Eq. 18) | For joint TESS + ground-based multiband fits. |
| Baseline alternatives: Legendre or quasi-orthogonal polynomial-harmonic functions (Eq. 11); weighted normalisation | Long segments with strong instrumental trends. |
| Joint solution with **radial velocities** (Eq. 19) and published minima times; hybrid Model II (template from a physical code) | Combined analyses as for XY Boo. |
| **Eclipse asymmetry** parameter | Not in the paper. My own gap. |
| Validation against the paper's examples | Done for the four one-star columns of Table 1 (section 7). Not done for XY Boo, `ebfit`, MAVKA. |

---

### 3b. Eccentric orbits (including high eccentricity)

**In the paper.**
- Circular orbits are assumed: Sect. 4.1.1 puts the eclipse centres at φ₀₁ = 0 and φ₀₂ = 0.5 and calls them exactly symmetric; Sect. 3.2 says the text concentrates on circular orbits and constant light curves, quoting more than 80% of EB orbits as fairly circular (CALEB).
- Footnote 3 says that with *phase rectification* (Mikulášek et al., in prep.) eccentric systems can also be solved. The technique "symmetrizes" light curves of binaries whose components move unevenly on eccentric orbits by rectifying the phase function, and it also allows apsidal-motion analysis. No formulae are given.
- Sect. 2.2 mentions apsidal motion only as one more possible phase-function model.
- I found no later publication of the rectification method in a quick web search. That is not proof that none exists.

**In the code.** Not supported.
- Eclipse 2 is fixed at phase 0.5 (`EBModel`), and the eclipses are symmetric.
- Separate D, Γ, C per eclipse are possible, which can absorb a *difference in duration* but not a displaced secondary or an asymmetric eclipse profile.
- The skill tells the agent to stop when the secondary is not near phase 0.5.
- Letting the timing shift absorb a secondary offset would hide the problem, because the template would still be fitted at 0.5.

**How others have used the model for eccentric systems.** A recent TESS study (*A Sample of Eccentric Eclipsing Binaries in TESS*, AJ) fits each eclipse of the folded light curve with the Mikulášek (2015) model to get its mid-eclipse phase and duration, then derives the eccentricity from those. I have only seen the abstract and snippets of that paper.

**Possible routes (not done; the first two are my own proposals, to be decided).**
1. *Free secondary phase.* Add φ₀₂ as a parameter (eclipse 2 at 0.5 + Δφ) with separate D₂, Γ₂, C₂. Cheap, and suitable for moderate e where each eclipse is short compared with the orbit. The standard low-e relations give e cos ω from Δφ and e sin ω from the duration ratio. They are approximations and must not be used for high e.
2. *Keplerian rectification of my own design.* Map time to the true longitude through Kepler's equation with free (e cos ω, e sin ω) and fit the usual profiles in the rectified phase. Conjunctions are then exactly 0.5 apart. For high e the eclipse shape is still not exactly symmetric and the proximity terms depend on the orbital distance, so this is an approximation, and it is not the author's method.
3. *Contact the author* about the rectification paper before investing in 2.

---

## 4. What the paper says that bears on these results

- **Accuracy floor.** The phenomenological models cannot, in the proposed form, do better than about 0.5 % of the total amplitude (Sect. 7). For AA And, removing the 1.95 mmag noise in quadrature, the two-component fit leaves about 2 mmag of misfit. That is roughly 0.4 % for a total amplitude of about 0.5 mag. I did not compute the amplitude precisely, so this is only an approximate agreement. The single-component misfit is about 4 mmag.
- **LC shape variability.** Sect. 6.4.3 and Fig. 14 show an inconstant O'Connell effect and a double wave in antiphase in XY Boo residuals. These lower the accuracy of the fit and of the eclipse timings. The author reports that the time scale of such changes is usually shorter than several days. This is consistent with the O–C scatter and the star-spot hypothesis in the AA And analysis, but the paper does not discuss primary–secondary anticorrelation.
- **Formal versus actual timing errors.** The paper states that the inner uncertainty of individual minimum times is typically about an order of magnitude better than their actual scatter about the model (6.2 min for XY Boo), caused by erratic period-rate fluctuations and LC shape changes. This supports reporting **empirical** uncertainties and not formal σ.
- **Existing timing program.** Conclusion 3 says Eq. 13 has been used for minima timings by the Czech Astronomical Society's variable-star section since 2011, citing Brát, Mikulášek & Pejcha (2012), `ebfit` (<http://var2.astro.cz/library/1350745528_ebfit.pdf>). I have not read it. It is the most natural reference for validating the timing procedure.
- **Model selection.** The paper warns that choosing the optimum model may be troublesome for beginners and that an inappropriate model may give non-convergent or bad solutions. The skill's decision tables address exactly this.

---

## 5. Suggested next steps (in order of priority)

1. Done in this step for template + shifts (`fit_global_shifts`). Next: include the ephemeris and intervals longer than one eclipse.
2. Done (section 7): Table 1 models overlaid on TESS and refitted for EK Com, AV Del, V477 Lyr, AR Aur. Still open: reproduce XY Boo (Tables 3 to 4) if the data are available.
3. Compare timings with `ebfit` and MAVKA on the same eclipses.
4. Subset-wise iterative weights.
5. Ask the author about flat-bottomed eclipses, the 2.1 extension and the status of the phase-rectification method for eccentric orbits.
6. Lower priority: profiles #5–#7, Legendre baselines, higher-order period terms, multicolour model.

---

## 6. Suggested methods paragraph (British English)

> Light curves were modelled with the phenomenological model of Mikulášek (2015), i.e. a sum of two eclipse functions (his Eq. 13), low-order harmonic terms for proximity effects, and an O'Connell term, fitted in magnitudes together with a polynomial baseline per TESS orbit and a linear ephemeris. For the flat-bottomed primary minimum of AA And the eclipse function was extended by a second component with independent width and kurtosis, a modification not included in the original method. Times of individual minima were obtained by fitting a time shift of the fixed template to the data around each eclipse, following the idea of virtual minima (his Eq. 9) but with the template held fixed. The quoted timing uncertainties are empirical, derived from the scatter about a smooth ephemeris, and not from formal fit errors.

---

## 7. Validation against Table 1 with TESS light curves (2026-10-29)

**What was done.** `working_scripts/paper_forward_check.py` builds each star's Table 1 model (AV Del evaluated at the TESS band, Lambda = 550/745.3 - 1 = -0.262, inside the BVRI range) and fits four variants to the TESS sector, each with a constant baseline per orbit segment: A paper values with only M0 fitted; B published shape fixed and amplitudes, M0, P0 free; C also D and Gamma free; D also C (C1 and C2 for EK Com) free. Results: `results/paper_forward_<star>.png`, `.txt` and the cross-star table `results/paper_comparison_summary.md`. Regression test: `tests/test_paper_comparison.py`.

**Result.**
- D and C agree with Table 1 within about 2.3 combined sigma for all four stars (EK Com also C1; EK Com C2 differs by about 2.2 sigma with the same ordering C1 > C2 as the paper).
- Gamma agrees for AV Del (1.0 sigma) and AR Aur (1.5 sigma) and is consistent for EK Com (1.2 sigma). V477 Lyr is the exception: 1.28(6) against 0.79(4), 6.7 sigma. The TESS data are noise limited there (29 mmag per point) and the band is very different (hot nucleus plus cool companion), so a band dependence of Gamma is a candidate cause; this was not tested.
- AV Del amplitudes at the TESS Lambda agree within a few percent (A1 0.8715 against 0.8465, A2 0.1740 against 0.1723, A2/A1 0.200 against 0.204). Variant A (paper values, only M0 fitted) already gives 13.6 mmag.
- EK Com: A2/A1 0.869 against 0.873. Residuals of about 20 mmag remain on the primary walls in all variants (single-Gamma ceiling for a flat occultation, 10-min exposure smearing, spots).
- V477 Lyr: TESS A1 is 0.27 of the paper depth. Third-light dilution alone does not explain it (a diluted paper model would give D about 0.015 and a sign change of C, which is not seen); a band effect is the untested hypothesis.
- AR Aur: shape agrees; amplitudes not compared (paper gives B and V only).

**Paper erratum candidate.** The Table 1 period of AR Aur, 4.13169 d, is contradicted by the data. Folding on it smears the eclipses (rms 8.6 mmag); the VOTable period 4.134651 d gives 1.7 mmag and a free fit converges to 4.134649(5) d. Probably a typo for 4.13469 d.

**Caveats.**
- Formal TESS errors are rescaled by the reduced chi2 and understate the systematic misfit (rchi2 from 1 for V477 Lyr to 2000 for AR Aur).
- The paper amplitude errors for AV Del ignore the correlation between a_k1 and a_k2.
- The comparison uses different epochs (decades apart), a different band and TESS apertures; amplitudes and spot-dependent terms (cos 2 pi theta, O'Connell) are not expected to match.
- The out-of-eclipse noise floor from `out_of_eclipse_mask` is inflated by orbital modulation in contact systems; use the rms in mmag.

---

## 8. Configuration and reuse (2026-10-30)

**Four layers** (`mikulasek_eb.config`), agreed with the developer:

| Layer | Class | Contents | Scope |
|---|---|---|---|
| Model | `ModelSpec` | profile, harmonics, O'Connell, K | the star; transferable between light curves |
| Data | `DataSetup` | baseline order, gap split, exposure, `nsub` | one light curve |
| Policy | `FitPolicy` | Gamma starts, clipping, two-component and eclipse-weight rules, thresholds | the algorithm; global defaults |
| Timing | `TimingRecipe` | window, baseline order, sigma cut, wall coverage, clipping, `template="own"` | kept separately so one recipe can serve several light curves |

`Config` bundles the four, validates them (`ConfigError`, unknown keys rejected), round-trips through JSON (`to_json`, `from_json`, `save`, `load`) and has a deterministic `config_hash()` (notes excluded) for cache keys and provenance. `FitConfig` is an internal flat view that the pipeline reads (not part of the public API; a test guards that its fields are exactly the union of the three layers). `run_staged_fit` accepts either form and stores the full `Config` in `FitPipelineResult.config`. `time_with_recipe` applies a `TimingRecipe`; `write_timing_ascii(recipe=...)` writes the recipe JSON and its hash into the ToM header.

**Rules.** Everything is plain data (no globals, no file paths, no data); the host application, a notebook or an agent can build, store and replay a config. JSON is the canonical file form. All former module constants of the pipeline and the two-component acceptance are now `FitPolicy` fields with the same defaults (a test pins them).

**Host-facing API (`mikulasek_eb.api`, done).** `fit(t_mjd, mag, sig, period_days, epoch_jd, config, data_id=None, progress=None)` returns `FitOutput` (`solution`, the arrays the fit used, the full pipeline result). `time(t_mjd, mag, sig, solution, recipe=None, which=...)` returns `TimingOutput` (timings and rejected events per eclipse type, recipe, solution hash); the recipe defaults to the one in the solution config. Both are pure: no globals, files or plots. `ModelSolution` is frozen, JSON round-trips, and has `solution_hash()` for cache keys. Errors derive from `MikulasekError`: `ConfigError`, `InvalidLightCurve` (ValueError), `FitDidNotConverge`, `NoValidSolution`, `FitCancelled` (RuntimeError). `progress(stage, fraction)` is called after each pipeline stage; raising `FitCancelled` in it cancels the fit. Module note: the engine module is `fitting.py`, since `mikulasek_eb.fit` is the function.

**Borrowing (`mikulasek_eb.borrow`, done).** `fit_borrowed(t, mag, sig, reference, data, plan)` fits a light curve with the structure and values of a reference `ModelSolution`. `ParameterPlan(free=...)` names the groups that are fitted (`shape`, `amplitudes`, `harmonics`, `zero_phase`, `period`); everything else stays fixed at the reference values. The new light curve has its own `DataSetup` (exposure, baseline, gaps); the reference supplies policy, timing recipe and the T0/period used for the epoch numbering E. Out-of-eclipse outliers are clipped as in the pipeline, then the fit is repeated. The result is a `ModelSolution` whose provenance holds the reference hash and the plan.

**AA And test (2026-10-30).** Reference: SPOC sector 16, 2-min (two-component primary). Borrowing with amplitudes, harmonics and M0 free, baseline order 3 per segment (gap 0.5 d), exposure integration at the estimated exposure: 30-min sector 16 rms 1.7 mmag (the 2-min data: 2.7 mmag with the same baseline), amplitudes differ from the reference by about 4 percent (different aperture). Timings of the 30-min light curve against the 2-min ones: 26 common primaries (mean difference 0.7 s, scatter 1.9 s, median sigma 4 s against 8 s) and 24 secondaries (0.7 s, 5.8 s, sigma 14 s against 28 s). Caveat: the two light curves share photons, so they are not independent tests. Sector 17 (30-min, neighbouring) fits to 6 mmag after clipping 16 points; its residuals show scattered-light ramps at the orbit ends and slow trends, so its timings need a closer look before use. Freeing the shape as well did not improve the fits noticeably.

**Decisions.** Timing `template="reference"` was not built; the same `TimingRecipe` plus the borrowed solution gives homogeneous timings. A separate depth scale for dilution was not needed here (free amplitudes suffice).

**Not done yet (planned).**
- Staged pipeline with some parameters fixed (not needed for borrowing, which uses `fit_global` directly).
- A joint ephemeris over several sectors.

## 9. Scope of the automatic selector (2026-10-30)

The selector (features, rules, evidence stage) applies to TESS data only; it
is not expected to work on Gaia. Calibration will use PHOEBE synthetic light
curves with TESS cadence and noise; ground truth is the fit quality, not a
visual choice of profile. See `feature_spec.md`, section 12.

## 10. TESS exposure and template smoothness (2026-10-30)

- **Exposure.** For every TESS cadence except 20 s the effective exposure is 0.8 of the cadence (onboard cosmic-ray mitigation, for example 480 s at a 600 s cadence); at 20 s it is the full cadence. This is a mission property that holds for SPOC, QLP and full-frame light curves. `tess_exposure_days` implements it; the features record the rule as `exposure_source` ("TESS: 0.8 x cadence"). Non-TESS data need an explicit exposure.
- **Sub-samples.** The exposure average is a midpoint sum over `nsub` samples. With the old minimum of 5 the sharp bottom of a profile with G < 1 produced visible corners in the template (up to 4 mmag at 30 min on NSV 807). The minimum is now 21; the template error against `nsub = 401` is 0.11 mmag. The fit rms barely changed (4.27 to 4.24 mmag), so earlier fits were not materially wrong.
- **Remaining corners.** A V-shaped minimum smeared by a boxcar becomes flat with corners at plus and minus half the exposure; that is the model, not numerical noise.
- **Check.** AA And: 30-min borrowed timings minus 2-min timings, primary mean 0.8 s with scatter 1.7 s, secondary mean 0.7 s with scatter 5.8 s (before: 0.7 and 1.9 s, 0.7 and 5.8 s).

## 11. A converged fit is not a solution (2026-11-03)

A least-squares run can return `success` with a template that is not an eclipse. On the two 1800 s light curves of a synthetic detached binary the primary became a needle: depth about 38000 mag and 5400 mag, Gamma about 10⁻⁶ and 10⁻⁵, while the rms on the samples stayed near 24 mmag, because a spike narrower than the exposure falls between the points. Reporting that as a fit is worse than saying that no solution was found. The paper (Sect. 5.1) uses Newton–Raphson and sets no bounds. The limits below are design choices of this package.

### 11.1 What a fit now carries

`assess_fit` (`mikulasek_eb.validity`) looks at the finished template and the data. `fit` stores the result in `solution.stats`:

| Status | Meaning |
|---|---|
| `ok` | No check failed. Numerically sound. Not a statement that the fit is good. |
| `bound-active` | A parameter or a derived width is outside the range allowed for an eclipse. |
| `unresolved` | The template could not be confirmed: the eclipse core is missed, the exposure integral has not converged, or a core has too few points. |

`status_reasons` says why the status is not `ok`. `status_warnings` does not change the status. `time` raises `UnreliableSolution` unless the status is `ok`, and also when a solution has no status (it was made before this check and must be refitted). `select_model` keeps every candidate in the table and in `result.fits`, and compares BIC only among those whose status is `ok`. If none are `ok`, `best` is `None`.

The checks, in `validity.py`:

| Check | Limit | Effect |
|---|---|---|
| Gamma | 0.3 to 30 | `bound-active` if on or outside the edge. Table 1 of the paper has 0.62 to 1.7. |
| Depth of one profile copy | at most twice the robust peak-to-peak range of the data (percentiles 0.5 and 99.5) | `bound-active`. The second component may be negative, with the same size limit. A is measured from the template mean, so harmonics can push it above the eclipse depth. |
| Full width at half depth of one copy | at most 0.5 cycles, and at least half the exposure in phase | `bound-active`. Measured on a fine grid. |
| Width D | 0.001 to 0.5 cycles | A solver box only, and only when the bounded solver is on. The width limit above is the one that judges the eclipse. |
| Curvature C | outside −0.6 to 0.3 | Warning only. |
| rms / out-of-eclipse noise | above 1.5 | Warning only. |
| Rank-deficient Jacobian | | Warning only. The parameter errors are unreliable. |
| Exposure integral | the model at `nsub` and at `4 nsub`, on 4001 phases, differ by more than the noise | `unresolved`. |
| Eclipse core | mean residual inside half a data-FWHM of each minimum, at least 3 points, within the larger of 0.2 of the folded depth and 5 standard errors | `unresolved`. |

A value within a small tolerance of an edge counts as on that edge. A bound violation outranks an unresolved check.

### 11.2 The solver

The default is Levenberg–Marquardt, with no bounds during the fit. The verdict is applied afterwards. That is enough to refuse a needle: both 1800 s curves come out `bound-active` and are not timed (`tests/test_validity.py`, marker `slow`).

`FitPolicy.bounded_solver` (default false) switches that one fit to the trust-region solver. Depth, width and Gamma are then kept inside the boxes above while the fit runs. D and Gamma are fitted as logarithms, and the covariance is returned in natural units. The verdict still runs afterwards, because a box on D and Gamma cannot express the full-width limit.

`FitPolicy.retry_bound_active` (default true) is the selection strategy. `select_model` fits each candidate with Levenberg-Marquardt. A candidate marked `bound-active` is fitted again with `bounded_solver=True`, and the verdict of that second fit is the one that is compared. A second fit that raises is logged, and the bound-active result is kept. Starting the fit with the trust-region solver stays off: on XY Boo that first-solver switch made every candidate `bound-active`.

The bounded solver was tried as the default and then taken back. On one fit of XY Boo with harmonics (2, 3), the same start:

| Solver | Wall time | Residual evaluations | Status | rms | A | H2 |
|---|---:|---:|---|---:|---:|---:|
| Levenberg–Marquardt | 0.40 s | 484 | `ok` | 4.23 mmag | 0.279 | 0.016 |
| Bounded trust-region | 0.95 s | 168 | `bound-active` | 4.01 mmag | 0.648 | −0.149 |

A and a cosine can cancel. The trust-region solver follows that valley until A sits on the cap (twice the robust range, 0.648 mag). The curve still looks like the star. The eclipse parameters do not describe the eclipse. With the bounded solver every XY Boo candidate was `bound-active`, so `select_model` chose none. Levenberg–Marquardt stops at the paper's set, harmonics (2, 3), and the selection test requires that. The same switch also exchanged the labels of the two primary components of AA And at an unchanged rms; with Levenberg–Marquardt the recorded primary depth matches the reference again. The one gain of the bounded solver was the 1800 s curve at spot longitude 180°, which it recovered to an `ok` fit. That is why it stays available and is not the default.

`examples/inspect_fit.py` shows the two steps on the two 1800 s curves. The first fit is the default solver. When its status is not `ok`, the same configuration is fitted again with `bounded_solver=True`, and both figures are saved (`..._bounded.png`). The comment on `refit_bounded` in that file is the short form of this section. On 2026-11-03 the 180° curve was refused by the default solver (A about 5400 mag, Gamma about 10⁻⁵) and the bounded refit was `ok` at 11.56 mmag. `time` returned 0 primary and 30 secondary: the deep minimum lies at phase 0.5 of the catalogue epoch, so the timed events are the half-integer ones. The 0° bounded refit was still `bound-active`, at 11.58 mmag, because one component had a full width of 0.80 cycles. That run uses `two_comp_max_rel_err=3` so the wide component is kept; the package default of that gate is 0.5.

### 11.3 The second component

A second copy of profile 1 on one minimum is a package extension (section 2.1), not in the 2015 paper. Both copies share the eclipse phase. The acceptance gates are `FitPolicy` fields and are unchanged: rms gain at least 0.25 unless `force_two_component` (which skips only that gain), at least two agreeing starts, D and Gamma of both copies positive, every relative error of A, D and Gamma below `two_comp_max_rel_err` (default 0.5), and the two amplitudes not a cancelling pair. `two_comp_rel_err_scope="core"` applies the relative-error gate only to the narrower copy. The default remains `"all"`.

Every staged fit records `solution.stats["two_component_trial"]`: whether the second copy was tried, whether it was adopted, which gates failed, and the rms it would have given. A disagreement of the eclipse times is stored as `timing_ok` and does not by itself reject the copy.

On the 120 s and 600 s spot curves the core-only gate changed two 600 s fits and none of the 120 s fits. Where a second copy was adopted, the usual shape was a narrow eclipse plus a wide shallow pedestal, and the half-cycle width limit is what then marks the fit `bound-active`. The comparison was made while the bounded solver was the default, so that table does not describe the present solver. The gate defaults were not changed.

### 11.4 One illustration: a spot on a detached binary

PHOEBE light curves of one detached binary, with a spot at longitudes 0° to 315° in steps of 45°, at cadences 120 s, 600 s and 1800 s. This is an illustration of the rules above, not a calibration of the selector.

The harmonic orders are measured only outside the eclipse windows (`propose_candidates`). A detached curve with a flat stretch proposes the significant cosines and no menu. A curve with continuous variation proposes (1, 2, 3) and (2,). The spot-timing runs had locked harmonics (2,) from the spot-free 120 s curve, with the O'Connell term forced on. That is a different experiment from letting each curve choose.

On the 120 s curve at 315°, `select_model` says "continuous variation, unequal depths" and the maxima differ by 39.5 sigma, so the O'Connell term is on. It chooses harmonics (1, 2, 3) and a straight baseline: one component, status `ok`, rms 15.62 mmag, and the second component is not tried (rms/noise below 1.5). The figure `examples/out/fit_from_scratch/h123_b1.png` still shows a residual wave of a few tens of mmag. The same file with only H2 adopts the second component, reaches the same rms, and is `bound-active` because the wide copy has a full width of 0.625 cycles (`h2_b1.png`). Putting H1 in by hand, on the earlier solver, lowered the shoulder rms by about 7 mmag and left the primary time shifted by about 50 s. An even term can lift both shoulders. It cannot move the minimum.

### 11.4b The instrumental term is in magnitudes (2026-11-04)

The baseline is a polynomial per segment in magnitudes. In flux that is a multiplicative trend: every flux in the segment is scaled by the same factor, so a minimum and a maximum shift by the same number of magnitudes. A trend that is additive in flux (for example a background error) does not do that. It changes a deep minimum by more magnitudes than a maximum. The magnitude baseline cannot remove such a trend, and the template then absorbs the difference. On the synthetic G4 at 600 s, where the generator adds a trend additive in flux, the primary of the two segments splits into two bands in the folded curve, and the rms falls from 20.2 to 12.1 mmag when the trend is estimated as an additive flux polynomial per segment (`examples/g4_c600_trend.py`, an experiment outside the package).

The generator adds its trend additively in flux because its author chose that, not because TESS does. A real light curve can have both kinds: an additive part (background, contamination) and a multiplicative part (pointing, aperture losses). The package does not decide between them. The real NSV 807 sector 97 file (QLP SAP, background subtracted) was tested with its own eclipses (`examples/nsv807_trend_kind.py`): the flux of the minima follows the flux of the maxima with a slope of 0.51 ± 0.06 (primary) and 0.61 ± 0.05 (secondary), which is the minimum-to-maximum ratio of the star, and with an intercept consistent with zero. A multiplicative trend predicts that, an additive one predicts a slope of 1. For that file the magnitude baseline is the right one. The same test on the synthetic G4 secondary gives a slope of 1.004 ± 0.009, as injected. One real star, one sector, a trend range of 2.4 percent. This does not generalise to other stars or sectors.

### 11.5 Not done

- The full-width limit is judged after the fit. It is not a constraint inside either solver.
- The amplitude-cosine cancellation on XY Boo is left as it is. Raising the amplitude cap, or bounding the cosines, was not chosen.
- A linear inner solution for the amplitudes, the cosines and the baselines is postponed.
- `min_wall` is still not a field of `TimingRecipe` (`todo.md`).
- The 27-file status count made before the solver choice was settled is not the result of the present default.
