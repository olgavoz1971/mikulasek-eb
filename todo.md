# Mikulasek EB package: work plan and progress

Scope: everything lives under the repository root. Nothing outside this folder is edited.
Reference: Mikulasek (2015), A&A 584, A8 (`docs/aa25244-14.pdf`); implementation notes in `docs/mikulasek_implementation_notes.md`.

## Goals

- (a) Successful AI-helped fit of a provided light curve.
- (b) A script that analyses a light curve by itself, chooses the model profile and the free-parameter set, and fits.
- Outputs of a run: the **model** and the **times of minima**, which may differ for primary and secondary. Timing code stays as it is for now; keep it in mind in every design decision.
- End state: a clean, tested, installable Python package.

## Decisions (developer, 29 Oct 2026)

1. TIC 450335412 is relabelled **EK Com**. No better-cadence EK Com data are available (10 min cadence stays).
2. Comparison with Table 1 uses shape parameters only (D, Gamma, C, A2/A1, harmonic ratios). Band and dilution differences accepted. Lambda extrapolation for AR Aur is a stated caveat. Every star gets a visual-inspection figure like the two existing ones.
3. All new code and changes stay inside the repository root.
4. Start with **paper defaults** (profile 1, paper parameter set). Two-component primary and eclipse weights are extensions, off by default. Intended pipeline: step 1 is light-curve analysis that results in the model profile choice for that star.
5. Silent-skip fix ordering: agent's call (see Phase 1).
6. Do not rewrite the timing part yet.
7. Rules: no code without an explicit request; discuss first. Fail fast, no silent fallbacks.

## Data

| Star | File | Cadence | Cycles | Notes |
|---|---|---|---|---|
| AA And (TIC 455790537) | `lc_discovery_TIC_455790537_SPOC_sector_16.vot` | 2 min | | Not in the paper. Already fitted (two-component adopted). |
| EK Com (TIC 450335412) | `EK_Com_sector_49.vot` (identical to `lc_discovery_TIC_450335412_...`) | 10 min | 89 | About 5.6 points per FWHM. |
| AR Aur (TIC 2236015) | `AR_Aur_sector_19.vot` | 2 min | 6 | Tmag about 6, possibly saturated. Paper has B and V only. |
| AV Del (TIC 383274804) | `AV_Del_sector_81.vot` | 3.3 min | 7 | Paper amplitudes depend on Lambda. |
| V477 Lyr (TIC 423311936) | `v477_Lyr_sector_26.vot` | 2 min | 53 | Noise about 29 mmag; dilution; strong reflection effect. |

## Phases

### Phase 0: housekeeping
- [x] Relabel EK Com: script `working_scripts/fit_ek_com_s49.py`, data `data/EK_Com_sector_49.vot`, results `results/ek_com_s49_*`. Rerun gave timings identical to the old TIC-labelled files; old script, duplicate VOTable and old results removed. Strict-cut files moved to `results/superseded_strict_cut/`.
- [x] Per-star manifest with Table 1 values: `references/stars_manifest.json` (transcribed from the PDF text; digits still to be checked against the PDF).
- [x] Record the current EK Com timing counts (primary 70, secondary 68) as the regression reference (stored in the manifest).
- [x] AR Aur period: the paper Table 1 value 4.13169 d is contradicted by the data (see the AR Aur forward-check entry); the VOTable value is right.
- [ ] Remove stray artefacts only after asking (`files (4).zip`, `.skill`).

### Phase 1: engine honesty (before new comparisons)
Reason for doing this first: the paper comparison must rest on an engine that does not drop events or accept failed fits silently. Counts may change; the reference counts above show by how much.
- [x] `virtual_minima`: no more silent `except`. Optional `rejected=` list gets one dict per skipped event (reason and detail); failed fits are logged at WARNING; coverage rejections at DEBUG. Output arrays unchanged. EK Com result: no fit failures, counts unchanged (70 and 68); the only final rejection is E=13.5 (one side uncovered).
- [x] Surface the `rejected` table in `minima_timing.py` output (header or log summary).
- [x] Convergence enforced: `FitDidNotConverge` raised by `fit_global`, `fit_minimum`, `fit_global_shifts` (inside `virtual_minima` it becomes a logged `fit_failed` rejection). Rank-deficient Jacobian logged at WARNING and stored as `FitResult.rank_deficient`. No effect on AA And or EK Com (no failures, no rank warnings).
- [x] Gamma start 0.5 added to the default multistart (`eb_tuning.multistart_global`, `FitConfig`). AA And and EK Com converge to the same minimum; EK Com timings move by at most 0.014 s. Still to test on a star that truly has Gamma < 1 (AR Aur, AV Del, V477 Lyr) and the mid-eclipse spike behaviour.
- [ ] `initial_guess`: the pipeline already stops with an explicit error when the folded secondary is not a depression (A2 < 0), which is the honest behaviour. The real fix is to measure the secondary against the local continuum, not the global median; this is a feature-extraction task, moved to Phase 4.
- [x] Docstrings (Google style), logging, units, no `print`. Units: the numerical core stays plain floats with units documented (days, cycles, mag); astropy units apply at the boundary (reader, time scales). Revisit in Phase 3 if wanted.
- [x] Tests in `tests/` (11 passing, about 6 s): synthetic parameter and timing recovery, timing error calibration, gap rejection reporting, fail-fast behaviour (non-convergence, rank deficiency, failed event fits), EK Com regression (timing counts 70/68, adopted rms). Run: `python -m pytest tests`. AA And regression still to add.

### Phase 2: paper validation
- [x] Forward check, EK Com: `working_scripts/paper_forward_check.py EK_Com` -> `results/paper_forward_ek_com.png` and `.txt`. Result: published shape (D 0.052, Gamma 1.7, C1 0.10, C2 -0.14) reproduces the TESS fold well; with paper amplitudes the depths are too large (V band, no dilution), with amplitudes refitted rms is 8.5 mmag, residuals up to about 20 mmag on the primary walls and 8 mmag out of eclipse. A2/A1 = 0.868 (paper 0.873). Refitted A1, A2, H2 are 0.93, 0.92, 0.85 of the paper values; paper H1 (0.023) and O'Connell (0.032) are not reproduced (TESS 0.0002 and 0.0037), as expected for a spotted star decades apart.
- [x] EK Com, shape freed (variants C and D in the same script): rms 8.5 (B) -> 7.5 (D, free) -> 6.2 mmag with D, Gamma, C1, C2 free. TESS D = 0.0498(3) (paper 0.052(3)), Gamma = 1.23(4) (paper 1.7(4)), C1 = 0.20(3) (paper 0.10(17)), C2 = 0.07(2) (paper -0.14(9)); formal errors are rescaled by rchi2 27 and understate the systematic misfit. Largest tension: C2 (about 2.3 sigma), same sign ordering C1 > C2 as the paper. Residual structure of about 20 mmag on the primary walls stays in all variants: single-Gamma shape ceiling plus 10-min smearing and spots, not a parameter choice.
- [x] V477 Lyr (`paper_forward_check.py V477_Lyr`, harmonics 1 to 3, no O'Connell, single shared C): TESS D = 0.0217(4) (paper 0.0219(5)), C = -0.337(16) (paper -0.342(20)), Gamma = 1.28(6) (paper 0.79(4), clear difference). Amplitudes: A1 0.42 (paper 1.56), A2 0.053 (0.078), H1 0.201 (0.321), H2 -0.024 (-0.035). Variants B, C, D all give rms 34 mmag, which is the per-point noise (sigma 29 mmag), so the data cannot separate the shape variants well. Quick check (ad hoc, not saved): diluting the paper model by third light to the TESS depth shifts D to 0.015 and flips the sign of C, which the TESS fit does not show, so plain dilution does not explain the amplitude drop. Hypothesis, untested: band effect (hot nucleus plus cool companion, so the primary is much shallower in TESS red; the paper has strong Lambda dependence of the amplitudes) and a possible band dependence of Gamma.
- [x] AV Del (`paper_forward_check.py AV_Del`; model evaluated at Lambda = -0.262: A1 = 0.8465, A2 = 0.1723, H2 = 0.0615; no cos 2 pi theta, no O'Connell, as in the paper). Variant A (all paper values, only M0 fitted) already gives rms 13.6 mmag. With the paper parameter set freed (variant D, rms 7.2 mmag): D = 0.0413(2) (paper 0.0396(7)), Gamma = 0.870(5) (0.82(5)), C = -0.2435(43) (-0.24(7)), A1 = 0.8715 (0.8465), A2 = 0.1740 (0.1723), A2/A1 = 0.200 (0.204), H2 = 0.0495 (0.0615). Best agreement of the four stars; confirms the Lambda mapping and the engine. Remaining residual: smooth out-of-eclipse asymmetry of about 15 mmag (positive near phase 0.2 to 0.3, negative near 0.6 to 0.7), which the paper model has no term for (A3 = 0, no sine); an odd harmonic or an O'Connell term would absorb it. Variants C and D coincide (C is not significant). Fitted P0 = 3.85334(4) d against VOTable 3.85343 d.
- [x] AR Aur (`paper_forward_check.py AR_Aur`, shape only; amplitudes not compared; P held at the paper value 4.13169 d as requested). Result with the paper period: D = 0.0173(1) (paper 0.0168(3)), Gamma = 0.618(4) (0.619(10)), C = -0.270(5) (-0.25(2)); A2/A1 = 0.878 (paper 0.797 in B, 0.828 in V, band dependent). But rms is 8.6 mmag with large in-eclipse residuals: the paper period smears the eclipses. Diagnostic (ad hoc, not saved): with the VOTable period rms is 1.69 mmag and D = 0.0175, Gamma = 0.604, C = -0.281; a free fit starting from the paper period converges to P = 4.134649(5) d. The paper table period is almost certainly a typo. Shape parameters agree with the paper within about 2 sigma either way; Gamma < 1 is recovered (0.60 to 0.62) from the multistart.
- [x] AR Aur reference run now uses the data period (free P from the VOTable value): rms 1.7 mmag, D = 0.0175, Gamma = 0.6035(8), C = -0.2807(10), P = 4.134649(5) d.
- [x] Paper-pure fits of all four stars (variants B to D of `paper_forward_check.py`; harmonics and O'Connell as in Table 1, shared D, Gamma, C, separate C1 and C2 for EK Com).
- [x] Comparison table: `results/paper_comparison_summary.md` (D, Gamma, C, A2/A1, rms, fitted period against Table 1; written by `paper_forward_check.py all`).
- [x] One figure per star: `results/paper_forward_<star>.png` (fold, primary, secondary, residuals; paper model and fit variants).
- [x] Implementation notes updated (`docs/mikulasek_implementation_notes.md`): paper remedy C1 != C2, new section 7 with the validation, paper period erratum candidate.

### Phase 3: package structure
- [x] Standalone installable package inside this folder: `src/mikulasek_eb/` with `model, fitting, timing, tuning, extensions, tom, helpers, pipeline, io, plotting`, `pyproject.toml` (dist name `mikulasek-eb`, extras `plot`, `tess`, `test`).
- [x] Extensions (two-component, eclipse weights and their acceptance checks) in `extensions.py`, off by default.
- [x] `sys.path` hacks removed (pytest uses `pythonpath` in `pyproject.toml`; scripts need the package installed or `PYTHONPATH=src`). [ ] One place for configuration (`FitConfig` or successor): still open, module constants are spread over `tom`, `extensions`, `pipeline`.
- [x] Timing module wrapped as is (`tom.py`); per-eclipse-type (primary, secondary) outputs preserved (EK Com 70/68 unchanged).

### Phase 3b: paper completeness (before Phase 4, agreed 2026-10-30)
- [x] Table 2 #1 to #4 verified against the printed table (page 8); #5 (Gaussian), #6 (triangle, half-width D), #7 (exp(1 - cosh)) added with `PROFILE_SHAPE` (which of D, G, C each profile uses). Tests in `tests/test_profiles.py` (Table 2 relations, fit recovery). Pipeline guard: `run_staged_fit` raises for profiles 5 to 7; generalising it belongs to Phase 4.
- [x] Eq. 14 K term: `EBModel(profile=1, kterm=True)`, keys `K1` (optional `K2`), off by default; synthetic transit test recovers K. Not yet in `FitConfig` or the pipeline (Phase 4 / exoplanet use), `wall_band` works through the model.
- [x] Exposure-time integration (not in the paper; it is model averaging over the exposure, i.e. time-binning of the model): `EBModel(exposure_days, nsub, period_days)`, `FitConfig(exposure_days, exposure_nsub)`, helpers in `exposure.py` (`estimate_cadence_days`, `estimate_exposure_days`, `suggest_nsub`), tests `tests/test_exposure.py`. EK Com: rms 5.89 / 5.90 / 5.91 mmag for none / 0.8 / 1.0 x cadence, so the ~20 mmag wall residuals are NOT exposure smearing; D and Gamma move to 0.0523 and 1.72 (paper 0.052, 1.7). TODO for volightcurve: a public method returning the exposure time when the data carry it (noted in `exposure.py`). Still open: automatic switch on cadence / FWHM (Phase 4), exposure in timing files header.

### Deferred (decided 2026-10-30, revisit when needed)
- Multicolour model, Eq. 16 (Lambda dependence a_k1 + a_k2 Lambda) with per-instrument instrumental terms (Eq. 10): when multi-band data arrive; also for a proper V477 Lyr band test.
- Legendre / quasi-orthogonal polynomial-harmonic instrumental functions (Eq. 11): same fit span for polynomials, only conditioning; harmonics only for long segments with strong trends.
- Joint solution of ephemeris, Pdot and all shifts; orthogonalisation including Pdot; higher period terms; exact Pdot solution (long time spans).
- Light-time effect, apsidal motion, eccentric orbits, radial velocities, Model II (out of scope for now; eccentric orbits stop with an explicit error).
- Eclipse asymmetry term (AV Del out-of-eclipse asymmetry of about 15 mmag).
- One place for configuration (`FitConfig` or successor): constants spread over `tom`, `extensions`, `pipeline`; settle before Phase 4 output format.
- Open validation items: V477 Lyr Gamma (1.28 against 0.79, band hypothesis untested), EK Com C2 (2.2 sigma), Table 1 AR Aur period erratum candidate (to tell Z. Mikulasek), XY Boo, `ebfit`, MAVKA comparison.
- volightcurve: request a public method that returns the exposure time of a light curve if the data have it (no VOTable standard); then read it in `mikulasek_eb.io`.
- Local-continuum measurement of the secondary depth in `initial_guess` (bright-secondary case), part of Phase 4 feature extraction.
- Write to Z. Mikulasek with the results (developer).

### Phase 3c: configuration (agreed 2026-10-30; single light curve first)
- [x] Step 1: layered `Config` (`ModelSpec`, `DataSetup`, `FitPolicy`, `TimingRecipe`), validation, JSON round trip, hash; former constants moved into `FitPolicy` / `TimingRecipe`; pipeline accepts `Config` or `FitConfig` and stores the config in its result; `time_with_recipe`; recipe and hash in the ToM header. 66 tests pass, regression numbers unchanged. Notes section 8.
- [x] Step 2: `ModelSolution` (config, values, errors, ephemeris, stats, exposure, provenance; JSON round trip, `solution_hash`, `build_model`); pure host-facing `fit(...)` and `time(...)` in `api.py` (`FitOutput`, `TimingOutput`); typed exceptions in `errors.py` (`ConfigError`, `InvalidLightCurve`, `FitDidNotConverge`, `NoValidSolution`, `FitCancelled`, all with ValueError/RuntimeError bases); progress hook `progress(stage, fraction)` (raising `FitCancelled` from it stops the fit). Module `fit.py` renamed `fitting.py` because the package attribute `fit` is now the function. 75 tests pass.
- [x] Step 3 (core): `ParameterPlan` (groups shape, amplitudes, harmonics, zero_phase, period) and `fit_borrowed(t, mag, sig, reference, data, plan)` in `borrow.py`; borrowed values stay fixed, own exposure, baseline and clipping; provenance records `borrowed_from` hash and plan. Real test on AA And (reference: SPOC S16 2-min): 30-min S16 rms 1.7 mmag; timings agree with the 2-min ones (26 primary, 24 secondary common events, mean difference 0.7 s, scatter 1.9 s and 5.8 s); S17 fits to 6 mmag after clipping 16 scattered-light outliers. 81 tests pass. Notes section 8.
- [ ] Step 3 open: timing `template="reference"` was NOT built: homogeneity is obtained by passing the same `TimingRecipe` and the borrowed solution to `time()`. Revisit only if a use case needs timing with an unrefitted reference template. Per-light-curve depth scale (dilution): not built; free amplitudes did the job on AA And (S16 30-min amplitudes differ by about 4 percent). Also open: S17 residual trends and scattered-light edges (baseline order 3 per orbit, flagged points), a joint multi-sector ephemeris.
- [x] Per-star scripts save `<stem>_config.json` and `<stem>_solution.json` next to the results (`api.solution_from_pipeline`). They still build their config as `FitConfig` (internal view); fine for now.
- [x] `FitConfig` stays as an internal resolved view only; it is not part of the public documentation (decided 2026-10-30).
- Use cases to keep in mind: host application on a remote server (stateless, cacheable, JSON), local scripts, Jupyter notebooks, analysis with an agent; not only TESS data.

### Phase 4: automatic profile selection
- [x] `features.py`: `extract_features`, `Features`, `FeatureSettings`; features in timing headers (EK Com script); `examples/features_table.py`; tests in `tests/test_features.py`. First results in `docs/feature_spec.md` section 8. Windows and contact indicator calibrated (spec section 9: edge-based windows, `separated_eclipses`, `contact_candidate`). Local-continuum depths and secondary search done (V477 Lyr secondary A2 = 0.07 mag found); `trend_order_hint` replaced by trend diagnostics (spec section 10). Open: hypothesis rules, evidence stage (baseline order with n_effective), a star near the contact boundary.
- [x] Feature spec drafted: `docs/feature_spec.md` (features A to F, draft rules, evidence stage, test stars incl. XY Boo S50). 
- [ ] Light-curve analysis features from the fold (depths, widths, flatness, out-of-eclipse curvature, reflection effect, secondary phase check).
- [ ] Hypothesis tree: detached, contact, semi-detached, reflection, eccentric flag.
- [ ] Model choice by an explicit criterion (BIC or similar on honest uncertainties), not ad hoc thresholds.
- [ ] Output: `FitConfig` plus evidence plus an "ask a human" flag.
- [ ] Validation: develop on some stars, test on the others; ground truth is Table 1 and the developer's judgement.
- [ ] AI-helped mode: the agent proposes a config from diagnostics; the script validates and records it.

## Progress log

- 2026-10-29: Review done. Found that TIC 450335412 is EK Com; paper D and Gamma agree with our fit. Decisions recorded above. This file created.
- 2026-10-29: Phase 0 done except the AR Aur period check. Phase 1: `virtual_minima` rejection reporting done; EK Com timings unchanged. Next in Phase 1: enforce `success` and singular covariance, Gamma < 1 starts, `initial_guess` for a bright secondary.
- 2026-10-29: Phase 1: convergence enforcement and Gamma 0.5 start done. AA And and EK Com fits unchanged; EK Com timings 70/68, changes below 0.014 s. Remaining in Phase 1: surface rejections in timing output, docstrings/logging/units in the engine, unit tests.
- 2026-10-29: Phase 1 tests added and passing (11). Remaining in Phase 1: surface rejections in timing output, docstrings/logging/units in the engine, AA And regression test.
- 2026-10-29: Phase 2 started. EK Com forward check done (see above). Note: rms/noise from the out-of-eclipse mask is misleading for contact binaries (noise floor 22 mmag against catalogue sigma 1 mmag); compare rms in mmag instead.
- 2026-10-29: EK Com paper-style fit done (variants A to D). Published D, Gamma, C1 agree within errors; C2 differs by about 2.3 sigma. Wall residuals of about 20 mmag remain. Next: map the other three Table 1 columns (AV Del with Lambda, V477 Lyr, AR Aur).
- 2026-10-29: V477 Lyr paper-style fit done. D and C agree with the paper, Gamma does not (1.28 against 0.79), amplitudes are 0.27 (A1) to 0.63 (H1) of the paper values. Noise-limited (29 mmag per point). Gap split for the forward-check script changed to 0.5 d (EK Com unchanged). Next: AV Del with Lambda, then AR Aur.
- 2026-10-29: AV Del paper-style fit done: shape parameters and amplitudes (at Lambda) agree with Table 1 within a few percent. Builders refactored (builder returns model, published parameters and published errors; one summary code path). Next: AR Aur (needs a decision on the B/V extrapolation and the period check).
- 2026-10-29: AR Aur shape-only check done with the paper period as requested. Shape agrees with Table 1; the paper period is a typo (data prefer 4.134649 d, rms 1.7 against 8.6 mmag). All four paper stars now have a forward check and a paper-style fit. Phase 2 remaining: summary table across stars, notes update, Phase 2 decision on next step.
- 2026-10-29: Phase 2 finished. AR Aur reference run uses the data period. Cross-star summary `results/paper_comparison_summary.md`, implementation notes section 7, regression test `tests/test_paper_comparison.py` (21 tests passing, about 11 s). Open from the validation: V477 Lyr Gamma (1.28 against 0.79, untested band hypothesis), EK Com C2 (2.2 sigma), Table 1 AR Aur period erratum candidate, amplitude comparison for AR Aur skipped by decision. Next: Phase 1 leftovers (rejection table in timing output, engine docstrings/logging/units, AA And regression test), then Phase 3.
- 2026-10-29: Phase 1 leftovers done. Rejected events (E from T0, predicted JD, reason) are written into the ToM headers via `time_both_walls(rejected=...)` and `write_timing_ascii(rejected=...)`; EK Com still 70/68, one rejected secondary (`one_side_uncovered`). Google docstrings for all functions in `mikulasek_eb.py`, `eb_tuning.py`; no `print` left in `scripts/` or `working_scripts/`; AA And regression test and rejection-header test added (23 tests). Next: Phase 3, awaiting confirmation.
- 2026-10-30: Phase 3 first pass: package `mikulasek_eb` under `src/` (old `scripts/` engine and `working_scripts` library modules split into modules and removed; `examples/` holds the two CLI examples). 23 tests pass, EK Com timings 70/68 unchanged. Open: editable install (`pip install -e .`, developer to run), one configuration place, folder rename `mikulasek_eb_2`, stray `mikulasek-eb-fitting.skill`, untracked `mikulasek_eb/` and `mikulasek_eb2.back/` at the repo root (not touched).
- 2026-10-30: Package moved to its own repository (/home/voz/projects/UPJS/mikulasek_eb, GitHub olgavoz1971/mikulasek-eb), black applied, README added, 23 tests pass. Folder rename and the editable-install note are obsolete. Still open: one configuration place; stray `.skill` file left behind in skvo_veb_2 history only.
- 2026-10-30: Agreed order before Phase 4: Phase 3b (profiles check and #5 to #7, K term, exposure-time integration); the rest is listed under Deferred.
- 2026-10-30: Phase 3b finished: Table 2 #5 to #7, K term, exposure integration. EK Com wall residuals are not exposure smearing (rms unchanged), so the cause stays open (single Gamma ceiling, spots, band). Next: Phase 4, after your decision on configuration design.
- 2026-10-30: Phase 3c step 1 done (layered configuration). Next: step 2 (`ModelSolution`, pure `fit`/`time` functions), then Phase 4.
- 2026-10-30: Phase 3c step 2 done (`ModelSolution`, `fit`, `time`, typed errors, progress hook). Next: step 3 (`ParameterPlan`, borrowing), then Phase 4.
- 2026-10-30: Phase 3c step 3 core done (borrowing tested on three AA And light curves).
- 2026-10-30: per-star scripts save config and solution JSON. Sector 17 handling and joint ephemeris postponed by the developer. Next: Phase 4 discussion.
- 2026-10-30: XY Boo (TIC 347405567, S50, 10-min) added to `data/`; feature spec drafted.
- 2026-10-30: features implemented and run on all seven light curves.
- 2026-10-30: eclipse windows now edge-based; contact indicator replaced by separated_eclipses plus contact_candidate.
- 2026-10-30: local-continuum depth, secondary peak search, trend diagnostics (red-noise ratio, correlation time, gap steps) in features.

## Phase 4: hypothesis rules and evidence stage (done)
- `candidates.py`, `selection.py`, tests in `tests/test_selection.py`, spec section 11.
- Open: n_eff for stars with no measurable correlation time (AR Aur); threshold calibration on more stars; (script `examples/select_models.py` done; accepts any .vot path).

## Calibration plan
- Selector is TESS-only (documented in README, spec section 12, notes section 9).
- Precise calibration: PHOEBE synthetic set with TESS cadence/gaps/noise (data to come).
- Meanwhile: supervised per-star analysis (run tools, judge by residual numbers and plots).

## Supervised star log
- AR Aur S19 (h12_b1): out-of-eclipse rms 0.58 mmag, in eclipse 2.9 mmag, worst phase bin 2.1 mmag. Catalogue sigma (0.16 mmag) is far below the real scatter, so rms/noise = 15 is misleading; the fit is good. The eclipse residual repeats in every cycle with the same shape at the knees (template shape defect, not timing: O-C scatter 2 s, equal to the formal error). Rule lessons: do not judge a fit by rms/catalogue sigma; n_eff not measurable on such clean data.
- Tool: `examples/diagnose_fit.py <label|path.vot>` (numbers + PNG).
- NSV 807 S04 (QLP, 30 min, P = 0.339 d, contact-like): file declared the flux zero point dimensionless; researcher confirmed 1 s-1, fixed with `working_scripts/fix_qlp_zero_point_unit.py` (assemble_volightcurve; TIMESYS TCB barycentre carried over; original files elsewhere untouched). Selected `h23_b1` (exposure integration on), rms 4.3 mmag (out of eclipse 3.3, in eclipse 5.8), phase-binned residual structure 3.2 mmag, worst bin 9.7 mmag at phase 0.95. Real noise is 2.1 mmag but the QLP catalogue sigma is 175 mmag (80 times larger): do not trust catalogue sigma of QLP files; flagged for the researcher. Only 1.9 points per FWHM, so eclipse-shape detail is limited; `h123_b1` is 3.8 BIC behind.
- Done: TESS exposure = 0.8 x cadence (20 s: full), `tess_exposure_days`; `MIN_NSUB` 5 to 21. Selections of all eight stars and the AA And borrowing/timing check rerun, no change of the chosen models.
- `docs/pipeline_explanation_for_dummies.md`: flow of `diagnose_fit.py`; profile is one Table 2 row (all seven implemented); automatic selector uses profile 1 only; a candidate is a configuration, not initial values.
- `docs/pipeline_explanation_for_dummies.md`: flow of `diagnose_fit.py`; profile is one Table 2 row (all seven implemented); automatic selector uses profile 1 only; a candidate is a configuration, not initial values.
- 2026-10-31: synthetic calibration plan agreed (`docs/synthetic_calibration_plan.md`): about 13 PHOEBE systems, TESS cadences 120/600/1800 s, 27 days, one gap, noise and trends; truth is the pure curve; times truth stored as T_conj and T_min_pure; O-C = 0 first; circular and eccentric timed and compared. Waiting for the files.
- 2026-10-31: `synthetic_check.py` + `working_scripts/check_synthetic_lc.py` (checker for generator output, 24 tests). First file (contact S05 c600): 18 pass, 5 warn (trend amplitude definition, depth value, conj_minus_min_pure_s only primary, no per-eclipse flags, smear_nsub), 0 fail.
- 2026-10-31: first synthetic batch (21 accepted rows + B22 with period_eps 1e-4) from LC_Generator. Checker fixed: truth checks use P_true (not the table P_meta) for conjunction times and the smeared truth; B22 now passes. All 22 files: 0 fail; only warn is G4 (Algol) secondary T_min_pure - T_conj = +1.47 s. Independent bisector of the pure G4 curve: secondary midpoints vary from 1 to 5 s with the level used (primary: about 0 +- 1 s), so the secondary eclipse is slightly asymmetric and the "true minimum" is ambiguous at the few-second level for that system; accepted, scoring should allow about 3 s there. Next: scoring harness (features vs truth, then candidates/regret, then timings).

## Scoring harness, step 1 (features vs generator labels), 2026-10-31
Script `examples/score_synthetic_features.py`, table `results/synthetic_scoring/features_vs_labels.md`.
- 19 of 22 files measured; G7 B16-B18 fail: P = 0.5 d is exactly 72/360/24 cadences, so the phase samples only 72 distinct values (fold sparse). Generator should use a non-commensurate period.
- Depth A1 (smeared, noisy data) vs label (pure): within 1% at 120 s, up to -10% at 1800 s (smearing, 3 points per FWHM).
- FWHM measured is systematically smaller than the label (-5 to -25%), worst for G4 primary (-23%) and contact (-12 to -18%). Check the label definition (flux or magnitude, level used) before blaming the features.
- FWHM ratio: G4 measured 1.26-1.49 vs label 1.07; G8 1.6 vs 1.48.
- Secondary offset (phase): G8 is about 0 although eccentric, so omega is probably near 90 deg (check generator).

## Scoring harness, step 2 (fits vs truth), 2026-10-31
Script `examples/score_synthetic_fits.py` (library: harmonics x baseline 0/1 x O'Connell on/off; the selector chooses among them), results in `results/synthetic_scoring/fits_vs_truth.{md,csv}`. Run took about 40 min for 19 files (`--resummarise` rebuilds the table from the CSV).
- Metric: rms of (fit - noise-free observation) in mmag; noise floor 10.8 mmag everywhere.
- The selector's choice is within 0.02 mmag of the best library member in every file (regret about 0). The ask_human flag is always yes (dBIC < 10 between near-identical candidates).
- Absolute error is not noise-limited: G1/G8 3.1 to 3.7 mmag, G2 5.3, G5 6.2, G6 7.5, G4 20 mmag. G4 (Algol, shallow secondary, narrow primary) is a real model failure of profile 1; the 3 to 7 mmag of the others includes the trend the linear baseline cannot follow, so shape and trend must be separated (next).
- Baseline order 1 always beats 0 by about 1 mmag in truth, and BIC prefers 1 by about 220: the generator trends are not constant.
- G1/G8 proposals are `h_b*` (no harmonics); the library lacked the empty set until now, so rerun G1 and G8 to include it.
- h123_* failed to converge in some detached files (known problem).

## Scoring harness, steps 2b and 3 (2026-10-31)
Scripts: `examples/score_synthetic_shape.py` (true trend removed, baseline order 0, 600 s files), `examples/score_synthetic_timing.py` (selector over baseline-1 library without O'Connell, default timing recipe, truth `T_min_pure`). Tables in `results/synthetic_scoring/`.
- Shape error alone (mmag, noise 11): G1 1.4, G8 1.2, G2 2.1, G5 2.3, G6 3.9, G4 6.8. Most of the earlier 3 to 7 mmag was the trend that the order-1 baseline cannot follow; G4 stays bad even without the trend (profile 1 limit). h123 fails to converge in several files.
- Timing, detached and contact without spots: bias within 7 s at 120 and 600 s, scatter 8 to 30 s, honesty 0.65 to 1.0 (errors slightly conservative). At 1800 s bias reaches 10 to 22 s, and for G5 B12, G6 B15 and G8 B21 (primary) no event is timed.
- G4 secondary: scatter 80 to 300 s (primary 8 to 28 s); honesty 1.1 to 1.4, errors slightly optimistic. Shallow, asymmetric eclipse.
- G6 (spot and O'Connell): bias +82 s (primary) and -98 s (secondary) against T_min_pure; T_conj - T_min_pure is +170 and -189 s. The O'Connell candidate was not in this timing library, so rerun with it on.
- G8 secondary bias -28 s (120 s file) to +98 s (1800 s file), not significant against the scatter 45 to 160 s.
- Not done: circular against eccentric comparison; G7 (period commensurate); rerun G1 and G8 with the empty harmonic set.

## Profile experiment (2026-10-31)
`examples/score_synthetic_profiles.py`: profiles 1 to 4, harmonics 23 and 123, O'Connell on/off, true trend removed, baseline 0, 600 s files. Profiles 5 to 7 are not in the staged pipeline (need `fit_global`), not tried. Shape error, best of each profile (mmag):
- G2 (total eclipse): p1 1.92, p2 2.16, p3 3.99, p4 2.44.
- G4 (Algol): p1 6.84, p2 6.84, p3 5.27, p4 7.11 (two-component).
- G6 (spot, O'Connell): p1 3.86, p2 3.86, p3 3.30, p4 5.12.
- Profile 2 equals profile 1 where C is not needed; profile 3 helps G4 and G6 but hurts G2; no profile reaches 1 to 2 mmag on G4/G6. One file each, so only a pointer for the tiered search.

## Where the G4 shape error sits (2026-10-31)
`examples/diagnose_synthetic_profiles.py G4_B08 [profile:harmonics:ocon[:2] ...]`, pictures in `results/synthetic_scoring/profile_pictures/` (`:2` forces the two-component trial with a 1% gain threshold). Error zones (mmag, 600 s, trend removed):
- profile 1, h123: primary core 19.9, primary wings 8.9, secondary 1.2, out of eclipse 0.9. The error is a ringing of +-25 mmag inside +-0.1 phase of the deep primary (1.65 mag); nothing outside.
- profile 3, h123: primary core 10.8, wings 9.7, but secondary 3.3 and out of eclipse 2.4 get worse.
- profile 1, h123, two-component forced: primary core 2.4, wings 4.1, secondary 1.2-1.5, out of eclipse 1.3.
- The two-component fit was tried by the pipeline and rejected: rms gain 14% < 25% criterion (`two_comp_min_rms_gain_frac`). The rms is dominated by the 11 mmag noise, so a 20 mmag structured error confined to 10% of the phase barely changes it. The acceptance test is blind to shape errors at this noise level.
- O'Connell on/off gives identical results for G4 (S1 is about 0).
- Candidate rule to test: judge two-component acceptance on the rms inside the eclipse windows (or on binned residuals), not on the whole-curve rms; calibrate on all synthetic files.

## Two components as a candidate (2026-10-31)
`examples/score_synthetic_twocomp.py` (600 s files, with trend, baseline 1, no O'Connell, harmonics 2/23/123; each fitted `single` and `two` = trial forced with threshold 0 and gain 1%, other acceptance checks kept; BIC chooses). Table `results/synthetic_scoring/twocomp_vs_truth.md`.
- BIC picks the two-component fit only where it was adopted: G4 (dBIC -240 against the best single). Everywhere else the forced trial was not adopted (other checks refuse it), so single and two are the same fit: no false two-component acceptance on G1, G2, G5, G6, G8.
- With the trend left in, G4 truth error moves only 20.44 to 19.57 mmag: the unmodelled trend dominates. In the trend-free run the primary core went 19.9 to 2.4, so the shape gain is real but hidden by the trend.
- G6 39.9 mmag here is because O'Connell was off in this library (S1 needed); not a two-component issue.
- Still to do: make `two` a proper candidate (config field instead of a policy hack), rerun at 120 and 1800 s, and replace the whole-curve rms gain criterion by an in-eclipse criterion and test it.

## Two components as a real candidate (2026-10-31)
Added `FitPolicy.force_two_component` (trial always runs, rms-gain criterion dropped, other acceptance checks kept; needs `allow_two_component`), `RuleSettings.two_component_candidates` (off by default: adds `<name>_2c` twins), 2 tests; 133 tests pass. Table `results/synthetic_scoring/twocomp_vs_truth.md` (all cadences, O'Connell on/off, single against forced two, BIC decides).
- BIC picks the two-component fit, and it is adopted, in 4 of 19 files: G4 at 120 and 600 s (dBIC -541, -240), G5 at 120 s (-437), G6 at 120 s (-341). Truth error gain is small with the trend in (0.2 to 0.9 mmag).
- Never adopted wrongly in G1, G8, G2 (G2 1800 s: twin is 0.1 mmag better but dBIC +2.3, BIC keeps single).
- G4 1800 s: not adopted (3 points per FWHM, no information).
- Tooling lesson: waiting with `while pgrep -f <script>` hangs because pgrep matches its own command line; use a PID or a log marker.

## Plan after the first scoring round (LC_Modeller, 2026-10-31; for discussion with the researcher)
State: selector regret about 0; trend dominates the truth error; detached shape error 1.2 to 1.4 mmag; G4 needs the second component (BIC picks it); timing good for detached and contact at 120/600 s, poor for G4 secondary, biased for G6; `ask_human` always yes (no information).
Order proposed:
1. Trend handling (largest remaining error): check which trends the generator applies (ramp, step, drift); test baseline order 2 and a split per TESS orbit on the synthetic files. Serves real stars too.
2. Two-component and O'Connell in the default proposals: switch on `two_component_candidates`; decide O'Connell by evidence instead of the maxima-significance rule (G6 needs it).
3. Calibrate `ask_human` with dBIC and selector regret: flag only when candidates differ in truth error, not for near-duplicates.
4. Timing recheck with the new selection: G4 and G6 with second component and O'Connell on; circular against eccentric on G8; T_conj offset for G6.
5. Hold-out: request part of the generator bulk (S01-S18) only after steps 1 to 3.
6. Rerun the 8 real TESS stars with the new defaults and compare with the earlier results and pictures.
7. Profiles 5 to 7 (need `fit_global`): low priority, only if something still fails.
Open, not ours to decide: G7 needs a non-commensurate period and the generator owes the FWHM label definition; pushing the ~20 local commits is the researcher's call.
Questions to the researcher: trend first or `ask_human` first? Should real stars fit the trend by default (baseline order 1 or 2) or expose it as a visible config choice?
