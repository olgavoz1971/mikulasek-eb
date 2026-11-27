---
name: mikulasek-eb
description: Fit phenomenological (Mikulášek-type) models to eclipsing-binary light curves from any survey, tune the fit, choose a model by evidence, and measure precise times of minima (virtual minima, O-C) with the `mikulasek_eb` package. Use this whenever the user works on an eclipsing binary (EB) light curve, fitting flat-bottomed or total eclipses, times of minima, eclipse timing, O-C diagrams, period or Pdot, or automatic profile and term selection, even if they do not name the method. For TESS light curves also read the skill `mikulasek-eb-tess`.
---

# Mikulášek-type EB fitting and minima timing

Goal: (1) an adequate global fit of the phase curve of an eclipsing binary, (2) precise,
honestly-quantified times of individual minima. Work like a careful observer: every decision is based on a
number or a plot, every run is logged, and conclusions state what was and was not verified.

## 0. What you have

The package `src/mikulasek_eb` (model, fitting, timing, tuning, extensions, tom, pipeline, plotting) is an independent re-implementation of Mikulášek (2015, A&A 584, A8),
not the author's code. It is validated against the paper's Table 1 (shape parameters, four stars; see `docs/mikulasek_implementation_notes.md`, section 7).
Read `references/api.md` for signatures. Related modules: `tuning` (starting values,
multi-start, residuals, window scan, clipping), `pipeline` (staged fit), `tom` (timing files).

Layers. Stay inside the layer the task needs.
- `mikulasek_eb` is the core (numpy, scipy, astropy). The root exports only the everyday surface: `fit`, `time`, `validate_lightcurve`, `fit_borrowed`, `run_staged_fit`, the result types (`FitOutput`, `TimingOutput`, `ModelSolution`), the configuration types (`Config`, `ModelSpec`, `DataSetup`, `FitConfig`, `FitPolicy`, `TimingRecipe`), `EBModel`, `PROFILES`, `phase`, `mjd_to_jd`, `absolute_jd_to_mjd` and the errors. Everything else is public in its own module (`model`, `fitting`, `timing`, `tuning`, `tom`, `extensions`, `exposure`, `validity`, `borrow`, `descriptors`, `helpers`); import it from there, not from the root.
- `mikulasek_eb.selection` chooses a model by evidence. It holds no survey numbers (section 7b).
- `mikulasek_eb.compare` measures a minimum without the model, to compare with other methods.
- `mikulasek_eb_tess` is the TESS adapter (reader, exposure rule, calibrated settings). Its own skill is `mikulasek-eb-tess`.
The core never imports the selection layer or the adapter, and needs no survey package.

Hard limits of the code. Say so to the user when relevant, never work around them silently:
- circular orbits only (secondary must sit at phase 0.5); eccentric systems are out of scope;
- exposure integration is optional (`exposure_days` on the data layer). Without it a long cadence biases D wide; see section 10 of the implementation notes;
- no eclipse-asymmetry parameter (an asymmetric minimum biases timings);
- the default solver is Levenberg-Marquardt. Every fit then carries a status (`ok`, `bound-active`, `unresolved`); `time` refuses anything but `ok`. See `docs/mikulasek_implementation_notes.md`, section 11. `FitPolicy.bounded_solver` is an optional second attempt, not the default;
- the optional second eclipse component (section 4b) is my extension, not part of the paper's method; label it as such in every report.

## 1. Conventions (get these right first)

- y in **magnitudes** (fainter = larger), sigma in magnitudes. Eclipse depth A>0.
- Times relative to a reference: t = BJD − T_ref (many surveys already give small offsets). Float precision matters.
- Phase 0 = primary minimum, 0.5 = secondary. Eclipse number E: integer = primary, n+0.5 = secondary.
- `M0` near the **middle** of the data (decorrelates M0 and P0): `M0 += P*round((t.mean()-M0)/P)`.
- `seg`: integer id per point for baselines, one per observing segment split at long gaps (for TESS: one per orbit, gaps > ~0.5 d).
- `order` (baseline polynomial degree): 0–1 for global fits of already-detrended light curves, 1 per event in timings. Degree 2 per orbit only if PDCSAP trends remain. Too flexible baselines eat the eclipse shape.
- Before fitting remove NaNs, bad-quality points and flares (how the quality is flagged is survey specific). Clip outliers on the residuals of a preliminary model, or on the scatter about the folded curve (for example 5 sigma, and outside eclipse first when a model already exists). A cut on the raw magnitude inside an eclipse removes the eclipse itself, because the star is faintest there; check any such cut against the folded curve before accepting it.

## 2. Standard workflow

Run the stages in order. Do not free everything at once.

**Stage A – data and noise.** Estimate the white noise from OUT-OF-ECLIPSE points only:
`nf = noise_floor(t, mag, mask=out_of_eclipse_mask(t, ephem, fwhm))` (after a first ephemeris/FWHM from
`initial_guess`). The unmasked estimate is inflated by steep ingress/egress (2.9 vs true 2.0 mmag in my test).
All "rms/noise" judgements use this masked value; compare with the catalogue `sig` and report both.

**Stage B – starting values.** `g = initial_guess(t, mag, P)` returns M0, A1, A2 and the primary FWHM
(phase units). Verify the period first: if P is from a periodogram it may be P/2 (compare
odd/even eclipse depths) and the secondary must be near 0.5 in the fold. If not, stop and report (eccentric).

**Stage C – structure.** Pick `profile` and `harmonics` from the EB type (section 4). Start with profile 1.

**Stage D – staged global fit.**
1. free `A1, A2, D1, G1` + `M0, P0` with all `H*`, `S1`, `C1` = 0 (quick, locates the eclipses);
2. add `H*` (and `S1` if maxima differ), then `C1`;
3. if the two eclipses differ in shape, add `D2, G2, C2` to the parameter dict and free them.
Use `multistart_global` for stage 2 (Γ starts 1, 2, 4, 8 with D consistent with the FWHM). If the
starts disagree, take the lowest χ²; if they agree, the optimiser is not the problem.
Keep a term only if its amplitude exceeds ~3× its uncertainty AND χ² drops clearly.
Free `Pdot` only if the data span is long (many hundreds of epochs); otherwise it is degenerate with P0.

**Stage E – diagnose.** `residuals(...)` then `binned_residuals(...)`. Read section 3 to translate the pattern
into an action. Iterate D→E at most ~6 times; log each attempt (parameters, χ², rms, decision).

**Stage F – timings.** `window_scan(...)` (section 6), pick the window, clean timings, fit an ephemeris, report.

## 3. Tuning playbook: symptom → check → action

First compute: rms of residuals vs `noise_floor`; reduced χ²; binned residuals around both eclipses and
out of eclipse; parameter errors from `r.err`; eclipse depth for scale.

| Symptom | Likely cause | Action |
|---|---|---|
| rms ≲ 1.5× (out-of-eclipse) noise, no structure | adequate | stop tuning |
| Residual ripple in eclipses: centre of the eclipse offset one way, shoulders (walls) the other, alternating sign | single-component ceiling for flat-bottomed/total eclipses (smooth function vs sharp contacts) | confirm rms/noise with the masked noise floor; if still > 1.5, try the two-component extension (section 4b); otherwise accept for timing. Do not add harmonics to fix it |
| D and Γ errors ≳ 10% and strongly correlated (`r.cov`) | degenerate valley, often with C1 | fix C1=0; refit; check multi-start; use profile 2 |
| Eclipse fitted too round/wide, Γ ≈ 1 while the data show a flat bottom | start in wrong basin or shapes shared by unlike eclipses or missing out-of-eclipse terms | multi-start; D from FWHM with `D_from_fwhm`; separate D2/G2; add H3/`S1` |
| χ² drops clearly when D2,G2,C2 are freed | the two eclipses differ (occultation vs transit, limb darkening) | keep separate shapes |
| One-sided residuals at maxima / wings of eclipses | missing harmonics or O'Connell effect | `harmonics=(1,2,3)`, `ocon=True` |
| Residuals oscillate with phase everywhere | too few harmonics, or spots changing in time | add harmonics (rarely beyond 3); otherwise use local timing windows |
| Residuals differ strongly between orbits/sectors | evolving spots/shape | fit per orbit or per sector; use local windows |
| Fit returns NaN, huge errors or D<0 / Γ<0 | unbounded optimiser ran away | read `status` and `status_reasons` (section 7b); restart from the best earlier solution; free fewer parameters; take abs(D) |
| rchi2 ≫ 1 but residuals look random | `sig` too small | rely on `rescale=True`; judge by rms vs noise floor |
| Secondary not at phase 0.5 | eccentric orbit | stop, report: out of scope |
| p0 has `H3` but nothing happens | `harmonics` does not include 3 | the model silently ignores unused keys: check |

Why Γ may stay near 1.2 on a total eclipse: a single Γ cannot give a flat centre and steep walls at the same time. Best Γ of profile 1 depends on (flat half-width)/(contact half-width) of the eclipse: 0.78 → Γ≈229 (rms 2.7 mmag); 0.58 → 10.8 (3.3); 0.39 → 2.3 (3.7); 0.19 → 1.2 (4.4) on a 0.39 mag deep synthetic occultation. A fit that sits at Γ≈1.2 with all multi-starts agreeing can therefore be the correct single-component answer for a short flat bottom (large radius ratio), not an optimiser failure.

Numbers I measured on synthetic total eclipses (use as orientation, not as law):
- A geometric flat bottom is represented by profile 1 with Γ ≈ 2–3, D ≈ FWHM/3, with a residual ripple of about ±7 mmag on a 0.39 mag deep eclipse (rms floor ≈ 1.2% of depth in an ideal case).
- From any Γ start (1–30) the fit converged to the same solution in ideal data; so repeated multi-start failures point to the data/model, not the optimiser.
- FWHM / D for profile 1 (C=0): Γ=1 → 2.24; Γ=2 → 2.9; Γ=3 → 3.2; Γ=5 → 3.6; Γ=10 → 3.9. Starting D = FWHM/2.2 is only right for Γ ≈ 1.

## 4. Recipes by EB type

Always start with profile 1; use 2 if C is consistent with zero or data are noisy/sparse (one parameter fewer).
Profile 4 is the Andronov-type function. Profile 3 is an alternative tail shape. From the paper, one test star (AR Aur primary) gave scatter 3.4/3.9/5.5/6.9 mmag for profiles 1/2/3/4: indicative only.

| Type | harmonics | ocon | notes |
|---|---|---|---|
| Detached, narrow eclipses | `()` or `(2,)` | no | D ≈ 0.01–0.03; baseline order 1 |
| Detached/semi-detached with reflection/ellipsoidal | `(1,2)` | if maxima differ | |
| Contact (W UMa) | `(2,3)` (XY Boo style) | if spotted | wide eclipses D ≈ 0.05; maybe `(1,2,3)` |
| Total / flat-bottomed (e.g. AA And) | `(1,2,3)` | yes | multi-start, `D_from_fwhm`, separate shapes; accept ripple ceiling |
| Very shallow secondary | as above | | do not over-interpret its shape; shared shape is often best; timings of the secondary have larger σ |
| Eccentric | not supported | | stop |

## 4b. Optional two-component eclipse profile (extension, not in the paper)

Adds A·profile(D,Γ) + Ab·profile(Db,Gb) per eclipse (keys `A<k>b, D<k>b, G<k>b`), so a flat bottom and steep walls can coexist.
Use ONLY after the single-component fit is tuned and the residual ripple remains (rms/noise > 1.5 with the masked noise).
`best, rows = refine_two_component(t, mag, sig, model, r_single, seg=seg, order=..., eclipses=(1,2))` tries 4 starts.
Accept the second component only if ALL hold (my heuristics):
1. rms/noise improves clearly (≥ ~25%) and residual ripple shrinks;
2. the best starts agree on χ² (at least two starts at the same minimum);
3. new parameters are finite, D>0, Γ>0, relative errors < ~50%; no degenerate pair (A and Ab of opposite sign, each larger than the total depth);
4. timings from both templates (both-walls events, same window) agree within combined errors; if they do not, report the difference instead of choosing silently.
Always report both fits. In a synthetic test (short flat bottom, 2 mmag noise, data gaps) the single component gave rms 2.2× noise and the second component 1.2×; this is a simulation, not proof for real data.

## 5. Cadence rules

Points per FWHM = FWHM_days / cadence, with FWHM_days ≈ (FWHM in phase) × P.
- ≥ 20: excellent; whole-phase and local timings both fine.
- 8–20: fine for the global fit; local timings possible but check σ.
- < 8 (long cadence on narrow eclipses): use folded data for the template; individual timings mostly rely on the template (whole phase). Report that D is biased wide because exposure smearing is not modelled.

## 6. Timings

`window` is a half-width in cycles. The points fitted for eclipse number E are those with `|phase - E| < window`. The full width of that interval is twice the half-width. A half-width of 0.25 runs from quadrature to quadrature (phase from -0.25 to +0.25). Do not describe it as "±window": that phrase is read as a full width.

Default point rule of `time` / `TimingRecipe`: at least **4 points** in the window and **2 on each side of the centre**. That is the floor for a template shift (the shift plus a straight baseline). The old floor of 10 was a round number, not a requirement of the method. The wall test (3 points on each steep side, 90 to 10 percent of depth) is **not** a refusal: each timed minimum has `walls_ok` (True/False). Set `require_contact_walls=True` only if the researcher asks to drop those events.

Whole-cycle timing is a separate method, `time_cycles(...)`: one shift of the whole template per orbital cycle (primary, secondary and harmonics move together). Default: 4 points, 2 on each side of the cycle centre. Offer it alongside window timing; do not replace one with the other. On detached and semi-detached stars this is often *worse* than a local eclipse window: most points are continuum (no timing information, but they vote in the least squares with the linear baseline), and the two eclipses are forced to share one shift. On contact stars the whole cycle carries shape, so the method can match or beat a window. After `time_cycles`, offer the cycle pictures (below); do not draw them unasked.

Two modes of window timing (same function, the half-width differs):
- **whole phase** (half-width `window=0.25`): template fixed, shift + baseline per eclipse; most precise for stable shape.
- **local** (half-width `window` about 0.06 to 0.15): robust against shape changes, spots, O'Connell variability; less precise.

`free=("A","D")` lets depth/width float per eclipse. Use only when depth/width demonstrably vary, with a half-width of at least 0.1; below that it gave occasional wild solutions in my tests.

Procedure:
1. `out, rows = window_scan(...)`: reports N, median σ, O–C rms and std(z) for half-widths 0.25, 0.15, 0.10 and 0.06.
2. Choose the half-width with the **smallest O–C rms**, looking at std(z) only as a calibration check. Never choose by smallest formal σ, and never by std(z) ≈ 1 alone: narrow half-widths reach std(z) ≈ 1 because σ explodes, not because the timing is good (AA And: half-width 0.06 had std(z) 0.87 and O–C rms 115 s). A rising O–C rms towards narrow half-widths is expected (fewer points on the walls) and says nothing about template adequacy. If a narrow half-width gives a clearly SMALLER rms than the whole phase, the shape varies: prefer local.
3. Require coverage of BOTH walls: `wb = wall_band(model, r.p)` then `window_scan(..., wall_band=wb)` runs each window with all events and with both-walls events only. Partially covered events rest on one wall and are biased by any wall-shape mismatch. In a synthetic test with data gaps, removing 3 of 53 events cut the scatter about the true times from 16.5 s to 4.5 s (half-width 0.25) and from 88 s to 4.5 s (half-width 0.15). A real-data sign of the same problem: O–C rms that grows sharply for narrower half-widths while std(z) stays ≈ 1 (big-σ outliers hide in z). Report all/walls side by side. Also remove events with huge-σ events at orbit edges dominate an unweighted rms (in a test σ ranged from 3.5 to 644 s). `window_scan` also cuts σ > 5× median and clips at 4σ about the ephemeris; report how many were removed and why.
4. `fit_ephemeris(T, sigma, E, deg)` with deg=1 (linear) or 2 (quadratic); inspect O–C for structure. Genuine O–C variation (e.g. light-time effect) inflates scatter; that is not a failure of the method.
5. Error calibration: with `rescale=True`, σ is scaled by √(rchi2). When the model misfits (flat bottoms), this makes σ conservative, so std(z) < 1 is expected and acceptable (I saw 0.5–0.6 on idealised synthetic data). std(z) > 1.3 means σ underestimated: check shape variability, baselines, and timing outliers.
6. Primary and secondary timings may be systematically offset (asymmetry, spots, light-travel delay in unequal-mass close binaries). Report them separately and the offset.
7. **Understand the scatter** with `oc_diagnostics(v)`. It reports primary/secondary O–C rms, the extra scatter beyond the formal σ (`sigma_extra`), the secondary−primary offset ± s.e., the correlation of each primary with the next secondary and the lag-1 autocorrelation. Interpretation (heuristics): O–C rms ≫ formal σ and std(z) ≫ 1 with a good template means correlated noise (spots, SAP systematics) that the single-event fit cannot see; a negative primary–secondary correlation points to starspots (anticorrelated shifts have been reported for Kepler EBs), so average neighbouring primary+secondary O–C or treat it as a systematic; a positive correlation or positive lag-1 points to common-mode trends (systematics, real period variation). In that case quote the EMPIRICAL uncertainty (≈ O–C rms, or formal σ ⊕ `sigma_extra`), not the formal σ.
8. **Check that the wall filter did anything.** If the "all" and "walls" rows of `window_scan` are identical, either there are no partial events (fine, say so) or `wall_band` was not passed; confirm which before concluding.
9. A better template does not automatically give better timings. In AA And the two-component template flattened the primary (ripple ±12 → ±1.4 mmag) but timings from the single and two-component templates differed by only 1 s rms against 16 s O–C scatter, so the scatter is not a shape problem there. Compare the two timing sets (acceptance check 4 of section 4b) and say honestly which of the two improvements (shape vs timing) was achieved.

### Window scale: a hypothesis, not a rule (read before step 2)

Step 2 above picks the half-width by the smallest O-C rms. That measures precision, not accuracy: a window that cuts the minimum shifts every time by the same amount, and the rms about the mean cannot see it. Synthetic tests on semi-detached stars (G4 and G7, 120 s) showed the failure, so do not set the window from the scan alone.

Working hypothesis, preliminary and tested only on two synthetic stars:
- Measure the width of the fitted eclipse function at 5 percent of its depth (main profile plus second component, no harmonics, no O'Connell term, no other eclipse): `mikulasek_eb.descriptors.eclipse_widths(solution, level=0.05)` returns `{"primary", "secondary"}` in cycles. It exists for every profile and second-component configuration and needs no contacts on the light curve. Do not use the parameter `D`: it is not comparable between profiles or components. The level must be strictly between 0 and 1.
- A first window is 1.5 times that width (full width), so the half-width is 0.75 times it: `descriptors.window_half_widths(solution, n_widths=1.5, level=0.05, max_half=0.25)` returns `{"primary": half, "secondary": half}`; pass each to `TimingRecipe(window=half)` for that kind. Both numbers are parameters. The factor is empirical. A half-window above `max_half` raises an error; it is never capped silently, so tell the researcher and ask.
- There is **no rule yet** for enlarging it. The secondary minimum and low cadences probably need a larger window. Do not enlarge it silently.
- Use the O-C scan (rms and mean O-C against the window, in units of that width) as a **control that warns the user**, never as a silent rule that changes the window. Show the scan and say what it suggests. The researcher decides.
- Do not change the window silently between sectors or cadences of one star. The idea to test is one window **in phase units** for all data of one kind of minimum of one star. This is a hypothesis; say so when you report it.

### Timing pictures (offer them every time)

The per-eclipse timings (`time`, `virtual_minima`) are measured in separate windows, not with one template over the whole cycle. A table of times does not show whether each window was sensible, so **after every timing run, offer the researcher the timing pictures and ask which to draw. Do not draw them unasked, and do not leave them unmentioned.** Offer all of them by default (one picture per eclipse, timed or rejected) and say how many panels that is. For a long light curve (more than about 30 eclipses) offer a chosen subset: the first, one typical, the worst O-C, and every rejected one.

- `mikulasek_eb.plotting.save_timing_windows(t_mjd, mag, sig, solution, timing, kind, out_dir, eclipses=None) -> list[Path]`. **One picture file per eclipse** (`E01_primary.png`, `E00_secondary.png`), never several eclipses on one sheet; the layout follows the MAVKA pictures. The number in the name is the cycle. A secondary is the same cycle as the primary of that integer, not a half-integer eclipse number. `plot_timing_window(..., kind, E, save_path=None)` draws one. The x-axis is MJD of the observations. The title gives the found time in MJD and its sigma in seconds (or the rejection reason). The whole-cycle picture stays full width. The zoom of the window is a **separate** file in ``<kind>/zoom/``, same name, smaller figure. The template is at the fitted position (red inside the window, grey outside). A red line marks the found time. An eclipse that was not timed is drawn in grey with its reason (`too_few_points`, `one_side_uncovered`, `walls_uncovered`, `fit_failed`). `eclipses` is a list of E numbers from T0. `plot_timing_window_zoom` draws one zoom.
- `mikulasek_eb.plotting.plot_timing_overview(t_mjd, mag, sig, solution, timings, save_path=None)`. One picture: the folded light curve over two cycles with the primary and the secondary windows both shaded and the points of the timed windows coloured. `timings` is one `TimingOutput`, or a list of them when the kinds were timed with different half-windows. Offer it with the pictures: it shows how much of the minimum the window covers.
- Whole-cycle pictures, same MJD layout: `plotting.save_timing_cycles(t_mjd, mag, sig, solution, timing, out_dir, cycles=None) -> list[Path]` after `time_cycles`. One file per orbit (`E00.png`); zooms of both eclipses under the *same* shift in `zoom/E00_primary.png` and `zoom/E00_secondary.png`. The zoom band is 1.5 widths at 5 percent, so it matches the window pictures. Continuum points used in the fit are dark grey; eclipse-band points are blue/orange. `plot_timing_cycle` and `plot_timing_cycle_zoom` draw one. Cycle numbers in the names are from T0.
- Pass the same arrays and the same `solution` as to `time`. The functions raise `ValueError` if `timing` belongs to another solution.
- Limits: an eclipse removed by the outlier cut (`clip_outliers=True`) is in neither list and is not drawn. Use `clip_outliers=False` to see every window.
- Use these pictures when comparing with another method (for example the local fits of MAVKA): the same eclipses, side by side.

## 7. Acceptance and stop rules

Declare the fit "adequate for timing" when: parameters finite with D>0, Γ>0; multi-start consistent; residual rms within ~1.5× the out-of-eclipse noise, or the remaining structure is the single-component ceiling and the two-component test (4b) either was done or was judged unnecessary; timings from both-walls events are stable between half-widths 0.25 and 0.15; and `oc_diagnostics` has been run and its verdict (white vs correlated scatter) is reported together with the uncertainty actually quoted. Declare "not adequate" and say why otherwise. Heuristic thresholds (mine, not from the paper): adjust them if the user disagrees.

Budget: at most ~6 tuning iterations per star unless the user asks for more. Stop earlier when χ² improves by < 1% for a new term. Report, don't loop.

## 7b. Automatic model selection (second level)

`mikulasek_eb.selection` measures features of the folded curve, proposes candidate models and chooses by fitted evidence (BIC; a difference of `DELTA_BIC_DECISIVE` = 10 is decisive). It has no calibrated defaults on purpose:
- `extract_features(t_mjd, mag, sig, period_days, epoch_jd, settings, ...)` needs a `FeatureSettings`, and the exposure through `exposure_days` or `exposure_rule`. Without an exposure it is recorded as "not given".
- `propose_candidates(features, rules)`, `review_features(features, rules)` and `select_model(..., rules=...)` need a `RuleSettings`. A missing settings object raises `ConfigError`.
- `tiered_search(..., features=, rules=, profiles=, start_profile=)` is a named strategy on top of `select_model`: baseline, then profile, then harmonics, with a stop rule. It returns a `TieredResult` (the winner, the final comparison, every step with its BIC gain, and why it stopped early). `on_step(label, result, seconds)` lets the caller tabulate or draw each step. `make_candidate` builds one candidate by hand.
- A survey adapter supplies the settings. For TESS use `mikulasek_eb_tess` (skill `mikulasek-eb-tess`). For any other survey never reuse the TESS numbers. Stop, say that no calibration exists, and ask the user whether to calibrate or to fit by hand with sections 2 to 4.
- No candidate finishing with status `ok` raises; it is never replaced by a default model. Report such a failure as it is.

### Preliminary analysis: `extract_features`

This is the public method for the preliminary analysis of a light curve. It fits no model. Use it to look at a light curve before any fit, and never re-implement a fold, a depth or a width by hand. The result is a `Features` object (plain data; `to_json`, `from_json`, `features_hash`). Groups, each a dict:

| Attribute | Contents |
|---|---|
| `sampling` | counts, span, cadence, exposure and its source, gaps, white-noise estimate, outlier fraction |
| `geometry` | primary epoch `m0_mjd`, reference level, depths `A1`, `A2` (below the local continuum) and `A1_global`, `A2_global`, depth ratio, secondary offset, `fwhm1`, `fwhm2`, points per FWHM, exposure over FWHM |
| `shape` | bottom ratio W90/W50, wing ratio W10/W50, bottom curvature, `gamma_estimate` from profile 1, `flat_bottom` |
| `variation` | eclipse coverage, edge widths, `separated_eclipses`, `equal_depths`, `contact_candidate`, maxima difference and significance, Fourier amplitudes outside the eclipse windows |
| `ephemeris` | odd and even eclipse depths (period doubled?), drift of the centres between the two halves |
| `residuals` | fold scatter, trend amplitude, red-noise ratio, correlation time, `n_effective`, steps across gaps |
| `notes` | `"<feature>: <reason>"` for every feature that could not be measured |

Rules:
- Phases and widths are in cycles, depths and levels in magnitudes, times in MJD.
- A feature that cannot be measured is `None`, and `notes` says why. Nothing is replaced by a default. Quote the note, do not fill the value in. A missing feature can change the candidates and raises a flag in `review_features`.
- A fold with fewer than 80 percent of its phase bins filled raises `InvalidLightCurve` (a sparse light curve or a wrong period). Report it. Do not lower the threshold to get past it.
- `review_features` turns the features into reasons to ask a person (no secondary, parity pointing to twice the period, secondary offset, fewer than 8 points per FWHM). Read them before fitting.
- The measurement choices are a `FeatureSettings`. The core has no defaults. For TESS use `mikulasek_eb_tess.extract_features`, which fills in the calibrated values and the TESS exposure rule.
- The meaning of every feature and the thresholds are in `docs/feature_spec.md`.
- The `Features` object does not carry the folded curve itself. A figure of the folded light curve with the measurements drawn on it is not available yet.

### Reading the fit status

`solution.stats["status"]` is set after every fit. `status_reasons` says why. Always read both.

| Status | Meaning | What to do |
|---|---|---|
| `ok` | No check failed. The template is numerically sound. It is not a statement that the fit is good: read `status_warnings` (rms/noise above 1.5, rank-deficient Jacobian, curvature outside -0.6 to 0.3) | Compare, plot, time |
| `bound-active` | A profile copy has Gamma outside 0.3 to 30, a depth above twice the robust range of the data, or a full width above half a cycle or below half the exposure. The numbers describe a numerical artefact, even if the curve looks like the star | Do not time it and do not quote its parameters. Free fewer terms, or let the selection try other candidates (below) |
| `unresolved` | The checks could not confirm the template: it misses the data in an eclipse core, the exposure integral has not converged, or a core has too few points | Look at the fold in the eclipses. Do not time it |

`time` raises `UnreliableSolution` for anything but `ok`. A low rms on the samples proves nothing: a spike narrower than the exposure falls between the points. The history and the numbers are in `docs/runaway_fits.md`.

### Solvers and retries

- The first fit is Levenberg-Marquardt with no bounds. The status is judged afterwards.
- `select_model` refits a `bound-active` candidate once with the trust-region solver (`FitPolicy.retry_bound_active`, true). Depth, width and Gamma are then kept in boxes. The verdict of the second fit is the one compared. A failed retry keeps the bound-active result.
- A candidate that raises `FitDidNotConverge` is not retried unless `select_model(..., retry_not_converged=True)` is passed. It is an experiment switch.
- Do not start with the trust-region solver (`bounded_solver=True`) as a default. On XY Boo every candidate ended on the depth cap.
- The retry costs a second full fit per `bound-active` candidate.

### Search order and `start_profile`

`tiered_search` does not fit every combination of profile, harmonics and baseline. Step 1 chooses the baseline order with the start profile and the richest harmonic set, step 2 tries the other profiles with that order, and step 3 tries the smaller harmonic sets. The baseline is fitted together with the template in every candidate. It is not taken from profile 1 and held. It follows that the baseline order is decided under the start profile, so a start profile that does not suit the star can bias it.

`start_profile` is required. Profile 1 is the ordinary choice. If step 1 raises `NoValidSolution` because every candidate of the start profile is `bound-active`, the search never reaches the profile step. Then:

1. Say so in the report. Do not hide it.
2. Run the same search again with another start profile: 3 first, then 4 if profile 3 also fails. Profile 3 was the only profile that finished with status `ok` on SV Cen sector 64. Sectors 99 and 100 of that star needed profile 4.
3. Report which start profile was used and why, and that the baseline order was decided under it.

Do not pick a profile by hand and skip the evidence. Do not widen the search to every combination to get around the failure.

The trust-region retry is off inside `tiered_search` (`make_candidate` switches it off). With it off, a step in which every candidate is `bound-active` stops the search.

## 8. Reporting

Give: star, data and cadence, period, profile and harmonics chosen and why, final parameters ± errors, rms vs noise floor, residual plot description, the window chosen with the scan table, number of timings, O–C rms, error calibration, the comparison with another timing method if one was made, and an explicit list of caveats (the hard limits in section 0 that apply). Save: parameter JSON, timing CSV (E, T, σ, n, rchi2), figures (folded fit + residuals, O–C), and the run log. Ask the user only when the choice needs domain knowledge (period alias, eccentric orbit, which minima to include).

## 8a. Comparing fits: depth and width, not D, Gamma, C

`D`, `Gamma`, `C` and `A` are correlated inside one profile and mean different things in different Table 2 profiles. Two fits of the same star with the same profile and harmonics can return very different `D` for an almost identical curve (on V477 Lyr TESS, `D` 0.0217 against 0.0166 with `C` free against fixed at 0, while the full width at half depth is 0.0407 against 0.0409 in phase). Do not report a disagreement or an agreement between fits from their `D`, `Gamma` or `C`. Use `mikulasek_eb.descriptors.eclipse_descriptors(solution)` (or `(model, p)`): depth in mag and full width at half depth in phase of each eclipse, and the two ratios. It measures the fitted eclipse curve alone, without harmonics, baseline or exposure integration, and returns point values only (no errors). Depths still depend on band and aperture, so compare them across fits of one light curve, not across bands.

## 9. Do not

- free all parameters at once, or add harmonics beyond 3 without evidence;
- free A, D when the half-width is below 0.1;
- time partially covered eclipses without reporting them separately;
- judge rms/noise with an in-eclipse-contaminated noise estimate;
- adopt the two-component profile without the four acceptance checks of 4b;
- trust formal σ from a fit with obvious structured residuals without checking std(z);
- invent paper details: the code follows the equations as transcribed from the paper; if something is unclear, say so;
- present agreement with another timing method as validation of accuracy when it only shows consistency.
