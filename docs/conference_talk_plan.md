# Conference talk: the `mikulasek_eb` implementation (first draft)

Status: first draft, 2026-11-02. Conference name, date and talk length: to be filled in.
Authors of the plan: Olga Vozyakova and LC_Modeller (AI assistant). To be discussed and edited.

Every completed step is reproducible from this document. A result is recorded only together with the working directory, the Python interpreter, the script and its arguments. Paths that a script fixes internally are named here too. An interactive check that was not saved as a script is marked as such and is not a source for a number in the talk.

## 1. Title and abstract (draft 2, 2026-11-02; the researcher adds her suggestions to this version)

Working titles (pick one, or mix):
1. *What is "the" time of minimum? Spots, photometric and geometric minima, and an open implementation of the Mikulášek method*
2. *Timing eclipses of spotted binaries: an open Mikulášek-method toolkit compared with WSAP/WSL*
3. *Photometric against geometric minima in TESS eclipsing binaries: a public toolkit for modelling and timing*
4. (short) *Be aware of the spot: how minima move, and a public toolkit to measure them*

Points the title should carry: (a) the toolkit is public and is an implementation of a known method;
(b) the scientific point is that "minimum" is not unique under spots. Title 1 carries both; title 4
is the shortest and fits the "warning" tone.

Abstract (draft 2, about 170 words). Square brackets are numbers to fill from the final runs; nothing in it is claimed beyond what is planned:

> The phenomenological model of Mikulášek (2015) describes an eclipsing binary with a few shape
> parameters and times the minima by shifting a fixed template. We present a public Python
> implementation (GitHub link) with model borrowing between light curves of different cadence, an
> optional second profile component for a distorted minimum, exposure integration for low-cadence
> data, a preliminary evidence-based choice of the model tuned on synthetic TESS-like light curves,
> and reproducible fit and timing files. A cool spot near an eclipse makes the time of minimum
> ambiguous. For a synthetic detached system with one spot at eight longitudes, we compare the
> geometric time of conjunction, several photometric definitions of the minimum, a Gaussian-process
> estimate, and the timings of the Mikulášek template and of the Andronov WSAP/WSL methods, at
> [120, 600 and 1800] s cadence. The offsets are [up to X s] and depend on longitude and on the
> method. We then show timings of NSV 807 and [two other stars] from several TESS sectors, with
> the caveats that spots impose on O-C interpretation, and recommend timing unflattened light curves.

Open points: keep or drop the GP in the abstract; name the two other stars or not; whether the
last clause is a sentence or only a slide; British spelling throughout; length limit of the
conference abstract (unknown).

## 1b. Message

The audience knows the Mikulášek (2015, A&A 584, A8) method well. The talk is about **this
implementation**, which is public on GitHub: what it adds, how it is built, what it does on
synthetic and real TESS light curves, and how its timings compare with the Andronov WSAP/WSL methods.
The first, synthetic part carries the scientific point: a spot can move the photometric minimum away
from the geometric one, and different methods then disagree. The aim is a **warning** ("be aware"),
not a recommendation.

Not the topic of the talk: a general study of TESS detrending. Flattening appears as one practical
recommendation with one example (slide 9).

## 2. Slide outline

| # | Slide | Content | Evidence in hand |
|---|---|---|---|
| 1 | Why another implementation | The method is known; what is new here; scope (fit, select, time) | package README, paper |
| 2 | Design | Layered `Config` (model, data, policy, timing recipe); reproducible solution file with hash and provenance; fail-fast errors, no silent fallbacks; exposure integration for low cadence | `config.py`, `solution.py`, 135 tests |
| 3 | Features | Second profile component for a minimum; profiles 1 to 7; borrowing a model from a higher cadence; timing options (window, minimum points, opt-in overlap) | sector 4 borrowing picture; G4 two-component results |
| 4 | Model selection, preliminary | Candidates (harmonics, baseline order, O'Connell) ranked by BIC; the human is asked when the evidence is weak. State clearly what was tested and what was not | scoring harness results on 18 synthetic light curves (`results/synthetic_scoring/`) |
| 5 | Spots and the time of minimum (synthetic) | One system, one cool spot, longitude sweep. Offsets of several definitions of "minimum" (bisectors at several depths, centroid, lowest point, GP estimate) against `T_conj` and `T_min_pure`; then the three timing methods on the noisy versions. Primary and secondary shown separately | to build (section 3.5, 4.2b); request sent to LC_Generator 2026-11-02 |
| 5b | Timing comparison, other synthetic cases | Mikulášek against WSAP/WSL without spots: bias, scatter, error honesty per method, cadence and window | to build (section 4.3) |
| 6 | Timing comparison, real | NSV 807, sectors 4, 30, 97: both methods on the same data | timing files in `nsv_807_oc/results/timings/` |
| 7 | NSV 807 O-C, worked example | Per-eclipse O-C, primary against secondary, changes between sectors; cautious interpretation | section 3 |
| 8 | Two more stars | Chosen to differ in cadence or period; same treatment | to choose (section 4.5) |
| 9 | Practical recommendation | Avoid flattening for timing: NSV 807 example where flattening put a wave with the TESS orbital period into O-C | researcher's example; sector 97 flattened and SAP files |
| 10 | Availability and limits | Repository link, how to run, known limits | section 4.6 |

Suggested balance: slides 5 to 8 (5 and 5b together) take about half of the time; slides 1 to 4 about a third; slides 9
and 10 the rest. Slide 4 is one slide only.

## 3. What we already have (numbers to quote, with the caveats)

### 3.1 Synthetic scoring (LC_Generator PHOEBE light curves, G1 to G8 and variants)
- The selector's regret against the best candidate is about zero; the additive trend dominates the error budget.
- Shape error is small for detached systems.
- G4 needs the second component. BIC chooses it where it helps (G4 at 120 s and 600 s, G5 at 120 s, G6 at 120 s) and never wrongly elsewhere.
- Timing is good for detached and contact systems at 120 s and 600 s. The G4 secondary is poor. G6 shows a bias to be rerun with the O'Connell term on.
- Caveat: all of it is tuned on a few TESS-like synthetic light curves. A hold-out batch (S01 to S18 from LC_Generator) has not been scored yet.

### 3.2 NSV 807 (TIC 35119266), W UMa, P = 0.3389614 d
- Light curves: sector 4 (QLP SAP, 1800 s), sector 30 (SPOC PDCSAP 120 s, and the user's FFI 600 s), sector 97 (QLP SAP 200 s, no flattening).
- Profile 3 describes the rounded contact eclipses best (S30 SPOC rms 4.68 mmag against 5.59 mmag for profile 1).
- Residual pattern, the same in all sectors and cadences: spikes of 10 to 15 mmag at minimum centres, a 5 mmag antisymmetric wiggle at both sides of each minimum, a slow 4 mmag wave. It belongs to the template, not to the noise.
- Sector 4 with the model borrowed from sector 30: rms 10.8 mmag with fixed harmonics, 4.63 mmag with H2, H3 and S1 free. The 2018 maxima are inverted against sector 30 (O'Connell term S1 = -0.022 against -0.011).
- Timing precision per eclipse: 3 to 5 s (S30, S97), about 12 to 13 s (S4, 1800 s). Rms(O-C/sigma) is 0.7 to 1.1 at window 0.25, so the quoted errors are about honest.
- Against the catalogue ephemeris: S30 about +144 s (primary) and +166 s (secondary), S97 about -135 s and -139 s; about -280 s change in 2020 to 2025. The sector 4 point is not yet on this diagram.
- Primary and secondary local periods differ within a sector (S30: 0.33896051 against 0.33896189 d), with the same sign at every window. Interpretation open (spots or a real effect).
- S97 shows a smooth wave of about +-7 s over 10 to 20 d, primary and secondary together; cause open (orbit-related systematics, edges, other).

### 3.3 Window test (S30 and S97, `nsv_807_oc/results/window_test/report.txt`)
- Wider windows (0.35, 0.5) shrink the formal sigma but increase the real O-C scatter; rms(O-C/sigma) rises to 1.4 to 2.1.
- Window 0.15 has the smallest scatter; 0.25 has about honest errors. Files for both exist.
- Consequence for the talk: the window is a tunable that must be reported with every timing.

### 3.4 Package changes made for this work
- `TimingRecipe.min_points`, `min_each_side`, `allow_overlapping_windows` (defaults reproduce the old behaviour); `virtual_minima(..., allow_overlap=False)`.
- `FitPolicy.force_two_component`; `RuleSettings.two_component_candidates`.
- These changes are not committed yet.

### 3.5 Spots and the minimum: ideas agreed on 2026-11-02 (discussion, no results yet)

Definitions
- Geometric minimum: the instant of conjunction (`T_conj` in the synthetic files).
- Photometric minimum: where the observed minimum looks symmetric or lowest; for a spotless symmetric minimum all definitions agree, for a spot-distorted one they do not. Candidate definitions: lowest point; bisector midpoint at several depths (for example 25, 50, 75 percent); centroid; KvW symmetry point; template shift (Mikulášek, window dependent); asymmetric-polynomial fit of the core (WSAP/WSL).
- `T_min_pure` (bisector of the pure curve) is one of them, not "the" ideal minimum. It stays in the comparison next to `T_conj`: **both truths are used** (two panels).
- Contact systems have a cusp at mid-eclipse (Gamma < 1), so the lowest point is poorly defined there; bisectors at several depths are more stable.

Gaussian process
- The GP implementation is in the Dash app (`skvo_veb/pages/gp_for_oc.py`, "Gaussian Process" item): scikit-learn `GaussianProcessRegressor`, kernel `ConstantKernel x [Matern(nu = 2.5) or RBF]`, noise scale and length-scale bounds as controls, optional length scale borrowed from well-covered cycles of the same light curve.
- Role in the talk for now: a shape-agnostic estimate of the curve, to apply every definition of "minimum" without a template and to give an uncertainty for each. It is **not** a fourth timing method in the comparison; that is postponed.
- To decide: whether the sklearn GP is called from the package side (a utility in the package, not importing the Dash page) or the analysis stays on the app side. The package should not depend on the Dash app.

Spot physics for the discussion
- A spot away from the eclipse gives a smooth wave that harmonics and the O'Connell term absorb; timing is hardly affected. A spot in or near an eclipse distorts ingress or egress asymmetrically and moves the photometric minimum.
- A spot fixed on one star often moves the primary and secondary minima in opposite directions; an orbital effect moves them together (the paper makes the same point). The mean of primary and secondary cancels the first-order effect only if the two shifts are equal and opposite.
- The geometric time can be recovered only by fitting a physical model with spots, which is outside this package. PHOEBE 2 has static spot features; whether it has built-in spot evolution is a question for LC_Generator (asked). One light curve per spot configuration is enough for the warning; evolution matters later for the interpretation of real O-C.

Design of the synthetic experiment (updated after the reply of LC_Generator, 2026-11-02)
- One system type only, a detached one, so the shift comes from the spot longitude and not from the shape of the eclipse. A contact binary such as G6, and NSV 807 itself, has a cusp at mid-eclipse, so the lowest point is poorly defined: a tiny change at the bottom moves it, and a figure comparing the lowest point, the bisectors and the centroid would mix the cusp with the spot. The detached system has a rounded bottom and clear wings, so a spot on one wing moves the photometric minimum by a measurable amount and the definitions can be compared. This does not mean spots matter only for detached stars. The experiment cannot be quoted for a contact binary; NSV 807 stays a separate, real-data example. LC_Generator therefore did not use G6. The base is a new **short-period detached** system, circular: P = 0.8731 d (not commensurate), q = 0.95, i = 87 deg, fill 0.32 / 0.30, Teff 5800 / 5600 K, blackbody atmospheres. Eclipse depth on the spot-free twin: 0.541 mag.
- Nine parents with prefix `SP_`: `SP_none` (no spot) and `SP_L000` to `SP_L315` in steps of 45 degrees (one cool spot on the primary, `relteff` 0.82, radius 32, colatitude 50). Observed at 120, 600 and 1800 s with the same seed for every file, so differences come from the spot. 27 accepted files; manifest `results/manifest_spotsweep.csv` in the generator repository. The weak-spot variant (`relteff` 0.92) was skipped.
- PHOEBE spots are fixed (no drift or size change within a light curve), one file per configuration. Longitude 0 faces the companion, 180 faces away; at primary eclipse a spot at 0 is partly hidden, at 180 it is on the visible outer face; 90 and 270 are the O'Connell (quadrature) cases.
- `T_min_pure` is now defined in the metadata (`t_min_pure_definition`): median of the bisector midpoints at 20, 35, 50, 65 and 80 percent of the way from the local flux minimum to the 90th percentile continuum, on the native-step pure curve in flux. So it is already a multi-depth definition. The spread between the five levels is itself a measure of the asymmetry, and should be kept in the analysis, not only the median.
- The FWHM label is defined too (`eclipse_fwhm_definition`; flux, 90th percentile continuum); our FWHM will come out 5 to 25 percent smaller for proximity-effect systems, as seen on G4, G5, G6.
- Output of our analysis: offset against `T_conj` and against `T_min_pure` per definition and per method, as a function of longitude, primary and secondary separately.

### 3.6 Spot sweep, first look at the metadata (2026-11-02, no timing run yet)
Values of `T_conj - T_min_pure` (seconds, from the file metadata of the 120 s files; P = 0.8731 d, circular orbit; the checker passes: 22 for the spot-free, 20 for the spotted files, offsets reported as info).

Where the two times come from (generator code `experiments/truth_eclipses.py`; neither is a PHOEBE input parameter):

- `T_conj` is the geometric conjunction time, computed from the orbital elements that were given to PHOEBE: the epoch of superior conjunction (`t0_supconj`), the period, the eccentricity and the argument of periastron. For this circular system the primary is at `t0_supconj + E P` and the secondary half a period later. The checker confirms these times match `epoch + E P` to 0.0000 s. A spot does not move `T_conj`.
- `T_min_pure` is measured on the noise-free one-orbit flux curve that PHOEBE produced. The measuring code does not call PHOEBE. Around each eclipse it takes the bisector midpoints at 20, 35, 50, 65 and 80 percent of the way from the local flux minimum up to the 90th percentile of the flux in that window, and `T_min_pure` is the median of those five midpoints. The spot changes that curve, so it moves `T_min_pure`.
- The table is the mean of `T_conj - T_min_pure` in seconds (file key `conj_minus_min_pure_*`). A positive value means the photometric minimum falls later than the conjunction. The whole offset, up to about 74 s, is the movement of the measured photometric minimum. The spread between the five bisector levels is not in the metadata and has to be computed from the pure curve.

| Longitude | primary (s) | secondary (s) | eclipse depth (mag) |
|---|---|---|---|
| none | +0.1 | -0.2 | 0.547 |
| 0 | -2.4 | -0.2 | 0.621 |
| 45 | +71.4 | -0.7 | 0.616 |
| 90 | +21.3 | -19.1 | 0.563 |
| 135 | +0.4 | -41.1 | 0.547 |
| 180 | +0.1 | -5.8 | 0.590 |
| 225 | -0.2 | +44.4 | 0.547 |
| 270 | -22.2 | +19.2 | 0.564 |
| 315 | -73.9 | +0.4 | 0.619 |

- The photometric minimum moves away from the geometric one by up to about 74 s (primary) and 44 s (secondary), against a timing precision of a few seconds at 120 s cadence in real data. The effect is large enough to matter.
- Pattern: at quadrature (90, 270) the primary and secondary shifts have opposite signs and similar size (about 20 s); at 45, 135, 225, 315 mostly one eclipse is shifted. This is what the paper and the spot discussion predict: a spot-driven shift moves primary and secondary in opposite directions, an orbital effect moves them together.
- The eclipse depth also changes with longitude (0.547 to 0.621 mag), so amplitude and timing are both affected.

### 3.7 Pure curves, measured (2026-11-02)

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
/home/voz/projects/UPJS/skvo_veb_2/.venv/bin/python examples/spot_pure_minima.py
```

No arguments. The script itself reads every one-orbit parent `*SP_*.vot` in `/home/voz/projects/UPJS/EB_synthetic_lc/results/detached` (names containing `_c120`, `_c600` or `_c1800` are the degraded curves and are skipped), measures the primary at phase 0 and the secondary at phase 0.5 with `measure_curve_minimum`, `kind="flux"`, window `0.12 * period`, and compares the bisector median with `conj_minus_min_pure_*_s` from the matching 120 s file. It writes `results/synthetic_scoring/spot_pure/pure_minima.csv` and `pure_minima.png`. The column `generator_conj_minus_min_s` is the table in section 3.6; that table was first read interactively and is not a separate script.

- The median of our five bisectors reproduces the generator's `T_min_pure` to 0.00 s for every longitude and both eclipses. The two codes measure the same thing.
- The five bisector levels do not agree with each other. Their spread is 0.2 s with no spot and up to 74 s (secondary at longitude 225). So "the photometric minimum" moves by tens of seconds depending on the depth used.
- The lowest sample barely moves. The primary stays at the conjunction except at longitude 45 (-156 s) and 315 (+120 s); the secondary stays at +2.1 s at every longitude. A spot shifts the shape long before it shifts the single lowest sample.
- The centroid over this window is not a time of minimum. It reaches about +-1200 s, because the window includes the out-of-eclipse spot wave and the centroid follows that wave. It needs a narrower window, or a continuum taken outside the spot, before it can go on a slide.

## 4a. Updated plan v2 (2026-11-02, after the second abstract draft)

Source of truth for the claims: `Abstract_Second_draft.md` (researcher's version). Every sentence of the abstract needs evidence by the conference; the table says what we have and what is missing. The GP is **not** in the abstract, so the separate GP package (4.2c) is off the critical path: priority C unless the spot figure needs it.

| Abstract claim | Evidence we have | Missing | Phase |
|---|---|---|---|
| Public Python implementation of the Mikulášek model, timing of minima | package, 135 tests | public repository ready (README, tutorial, licence, clean clone, CI) | P4 |
| Semi-automatic selection: candidates, BIC, flag for review; preliminary tests on synthetic TESS-like curves | scoring harness on 18 PHOEBE files, regret about 0; flag raised in every case so far | one honest summary table (what was tested, hold-out status); the flag wording matches the code (done in the abstract) | P2, P4 |
| Reuse of a high-cadence model for lower cadence | S4 borrowed from S30: rms 10.8 to 4.63 mmag | a short, simple-to-use example in the package (public API) | P4 |
| Second profile component | G4 and G6 results | one clean figure and the statement of limits | P2 |
| Finite exposure of long-cadence data | `exposure.py`, tests | one demonstration of the effect (smeared against unsmeared fit, 1800 s) | P2 |
| Comparison with WSAP/WSL on synthetic TESS-like curves (PHOEBE, noise, instrumental trends) | our timing against truth on G1 to G8 (`timing_vs_truth`) | WSAP/WSL results on the same files; common protocol | P1, P2 |
| Spotted system, geometric time of conjunction against photometric minimum | 27 `SP_*` files from LC_Generator (received) | all analysis | P1 |
| Real TESS light curves of NSV 807 and two other stars, preliminary O-C with uncertainties | NSV 807: timings for S4, S30 (SPOC, FFI), S97 with ASCII files, window test | WSAP/WSL on the same data; sector 4 on the O-C diagram; flattened pair; two more stars | P3 |
| Cases where each method performs better | window test (ours) | the common protocol results | P1 to P3 |

### Phases (order of work)
- **P1, now: spot sweep (synthetic, first part of the talk).** Steps: (1) first look at the 27 `SP_*` files with the checker and plots; (2) code for the definitions of the minimum on a pure curve (bisectors at several depths, centroid, lowest point) in the package utilities, with tests; (3) our timings of the noisy files at the three cadences, per eclipse type, against `T_conj` and `T_min_pure`; (4) the list of files and the output format for the WSAP/WSL runs of the researcher; (5) the figure: offset against spot longitude, per definition and per method, primary and secondary separately.
- **P2: common protocol on the wider synthetic set.** Fix protocol (truth, sigma convention, event selection, statistics: bias, scatter, rms(error / actual), rejected fraction); run ours and WSAP/WSL on G1 to G8 and `SP_none`; one table for the selection results; the second component and exposure demonstrations; ask LC_Generator for the hold-out batch only after the protocol is frozen and score it once.
- **P3: real stars.** NSV 807: sector 4 into the O-C (TDB), cross-sector O-C with the cycle count, WSAP/WSL on the same files (researcher), flattened against unflattened sector 97 as one practical example, the sector 97 wave and edge events; choose and process two more stars with the same recipe.
- **P4: repository and tools.** Commit and push decisions; paths, README, licence, CITATION, tutorial; a simple borrowing example; calibrate or document the review flag; time-scale handling between sectors.
- **P5 (priority C): separate GP package.** Only if the spot figure needs a shape-agnostic minimum, or after the conference. Read `docs/gp_extrema_modeller.md` in the host application first (it already designs the GP extremum tool).

### Next step
Mikulášek timings of the degraded `SP_*` curves (120, 600 and 1800 s). The script for that run will be named here, with its arguments, before the numbers are quoted.

## 4. Detailed checklists (by topic)

Each item has an owner (O = Olga, M = LC_Modeller) and a priority (A must, B should, C if time).

### 4.1 Repository ready for the public (A)
- [ ] Commit the pending package changes; decide what to push (O).
- [ ] No personal paths in `examples/` and `working_scripts/`; scripts run from a clean clone (M).
- [ ] README: install, a 20-line example from light curve to timing file, paper reference, scope and limits (M, O reviews).
- [ ] Licence, CITATION file, version number (O decides the licence).
- [ ] Decide about the real-star data: include, or give a download script (O).
- [ ] CI that runs the tests (M); run the tests in a fresh environment.
- [ ] One runnable tutorial notebook or script, light curve to timing file (M).

### 4.2 Package items that the talk exposes (B)
- [ ] Save the timing recipe and the window in every output (already in the `#` header of the timing files).
- [ ] Time scale handling: record and convert the time scale between sectors (S4 file says TCB, the researcher says TDB; fix the data or the loader's check).
- [ ] Decide whether `fit_global_shifts` should support overlapping windows (probably not for the talk).
- [ ] Rerun G6 timing with the O'Connell term on; rerun G1 and G8 with the empty harmonic set.

### 4.2c Separate GP package (decision of the researcher, 2026-11-02)
- Decision: a **separate package** for Gaussian-process analysis, used both by the host application (Dash) and by offline timing scripts. It does not live inside `mikulasek_eb` and does not import the Dash page.
- Working name: to be chosen (placeholder `gp_extrema`). Own repository, installable, with tests, in the same style as `mikulasek_eb` (public API, typed errors, no silent fallbacks, Google docstrings, `logging` and no `print`).
- Starting point: ideas that already work in the host application (`skvo_veb/pages/gp_for_oc.py`, "Gaussian Process" item, and the code behind it): scikit-learn `GaussianProcessRegressor`; kernel `ConstantKernel x [Matern(nu = 2.5) or RBF]`; noise scale as a multiplier of the given or MAD-estimated errors; length-scale min, initial and max; length scale borrowed from well-covered cycles of the same light curve for sparsely covered extrema.
- First scope for the talk: given a (pure or noisy, folded or per-cycle) light curve, return the GP curve with its uncertainty, and the time of the minimum under several definitions (lowest point, bisector at chosen depths, centroid), each with a posterior spread (samples of the GP). This is the "shape-agnostic true minimum" for slide 5.
- Later scope: per-cycle timing for the comparison, once it is wanted (postponed).
- Decisions needed from the researcher: the name; whether the host application is changed to use the package now or later; the licence; the rule for the kernel family (Matern against RBF) and how the length-scale bounds are set for the synthetic study (they must not be tuned per curve).
- Needs a look at the host-application code in detail before the package API is fixed (M, after the researcher's go-ahead).
- Open question under the project rules: the dependency order (Astropy, SciPy, NumPy) has no entry for a GP library. Using scikit-learn's GP is a researcher's decision, recorded here.

### 4.2b Spot sweep (A, first part of the talk)
- [x] Request to LC_Generator written (`agent_mailroom/LC_Modeller_to_LC_Generator_spot_sweep.md`, 2026-11-02) (M).
- [x] Reply received 2026-11-02: 27 files `SP_*`, detached base, definitions added (section 3.5).
- [ ] Acknowledge the reply and confirm the files to LC_Generator; ask nothing new unless the first look finds a problem (M).
- [ ] First look at the sweep files: run the checker, plot the pure curves, check `T_conj` against `T_min_pure` per longitude (M).
- [x] Code for the definitions of the minimum on a pure curve (`mikulasek_eb.measure_curve_minimum`: bisectors at 0.20, 0.35, 0.50, 0.65, 0.80, their median, the centroid and the lowest sample; flux or magnitude; tests in `tests/test_curve_minimum.py`) (M, 2026-11-02).
- [ ] GP-based minimum on the pure and on the noisy curve: decide where the code lives (O + M); reuse the settings of the Dash GP item.
- [ ] Figure: offset against longitude for each definition, with `T_conj` and `T_min_pure` panels (M).
- [ ] Run the three timing methods on the noisy versions; same figure with the methods (O supplies the WSAP/WSL results, M the Mikulášek ones).
- [ ] Wording of the warning slide and of the discussion (photometric against geometric minima; spot evolution research against orbital research) (O).

### 4.3 Timing comparison with WSAP/WSL on synthetic data (A)
- [ ] Agree the protocol: same light curves, truth `T_min_pure`, same sigma convention, same selection of eclipses (O + M).
- [ ] Run both methods for each cadence (120 s, 600 s, 1800 s where available) and for the contact, detached and spotted cases.
- [ ] Report bias, scatter, rms(error / actual) and the fraction of rejected events.
- [ ] Hold-out batch S01 to S18 from LC_Generator: score once, after the protocol is fixed, no tuning afterwards (M asks LC_Generator).
- [ ] Postponed: Gaussian-process template followed by template shifting.

### 4.4 NSV 807 (A)
- [ ] Add the sector 4 point to the O-C diagram; handle TCB/TDB (S4 is TDB according to the researcher).
- [ ] Build the cross-sector O-C with the cycle count shown, with primary and secondary separately.
- [ ] Show the flattened against the unflattened O-C of sector 97 (the researcher's own example; M prepares the pair from the existing files).
- [ ] Investigate the S97 wave and the segment-end events; decide how edges are treated (O + M).
- [ ] Run both methods (ours and WSAP/WSL) on the three sectors with the same data files.
- [ ] Interpretation text: state the facts (change of about -280 s, different primary and secondary periods), keep the explanations (spots, light-time effect, other) as hypotheses. Do not apply a light-time fit mechanically.

### 4.5 Two more real stars (B)
- [ ] Choose two stars that differ from NSV 807 in cadence or period range, with an unflattened light curve (O).
- [ ] Process them with the same recipe, no tuning by star; report the same table as for NSV 807 (M).

### 4.6 Slides and wording (A)
- [ ] One-page results table: method, light curve, window, n events, median sigma, O-C scatter, rms(O-C/sigma).
- [ ] The model-borrowing figure (S4 before and after) and the profile-fit figures (`nsv_807_oc/results/step5/`).
- [ ] Slide 4 wording: "preliminary", with the test set and the hold-out status stated.
- [ ] Known limits slide: profiles 5 to 7 only through `fit_global`; the staged pipeline does not free the K term; no light-time or apsidal fit; selection tuned on TESS-like synthetic data; the ask-the-human flag is always raised so far (to calibrate).

## 5. Risks and open questions
- "Mikulášek sometimes wins": it must come from the common protocol (4.3) and not from selected cases.
- The over-learning risk: the selector and the scoring thresholds were tuned on few synthetic light curves; the real stars are the independent test.
- The cause of the S97 wave and of the primary-secondary period difference is open; do not present either as a discovery before the checks in 4.4.
- Sector 4 uses a constant sigma (3.64 mmag) chosen from the residuals, because the file errors (175 mmag) are unusable; say so on the slide.
- The spot sweep has one spot size and contrast only; the dependence on contrast is not tested (weak-spot variant skipped). Say so on the slide.
- Conference details (date, talk length, abstract deadline) decide how much of section 4 is feasible; item priorities should be revisited when they are known.

## 6. Recipes already run

Interpreter for every command below:

```bash
PY=/home/voz/projects/UPJS/skvo_veb_2/.venv/bin/python
```

### 6.1 Synthetic scoring (section 3.1)

From `/home/voz/projects/UPJS/mikulasek_eb`. With no arguments each script reads `/home/voz/projects/UPJS/EB_synthetic_lc/results/manifest.csv` and writes `results/synthetic_scoring/`.

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
$PY examples/score_synthetic_features.py
$PY examples/score_synthetic_fits.py
$PY examples/score_synthetic_shape.py
$PY examples/score_synthetic_timing.py
$PY examples/score_synthetic_profiles.py
$PY examples/score_synthetic_twocomp.py
$PY examples/diagnose_synthetic_profiles.py G4_B08
```

### 6.2 Spot sweep, pure curves (sections 3.6 and 3.7)

The command is in section 3.7. Tests of the measurement itself:

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
$PY -m pytest -q tests/test_curve_minimum.py -p no:cacheprovider
```

### 6.3 NSV 807 (sections 3.2 and 3.3)

From `/home/voz/projects/UPJS/nsv_807_oc`. These commands regenerate the files now in `results/`. `scripts/prepare_files.py` as it stands writes only the FFI copy `data/lc_tess_FFI_sector_30_user_prepared.vot`.

```bash
cd /home/voz/projects/UPJS/nsv_807_oc
$PY scripts/prepare_files.py
$PY scripts/prepare_qlp_sigma.py
$PY scripts/step1_select.py data/lc_discovery_TIC_35119266_SPOC_sector_30_prepared.vot
$PY scripts/step1_select.py data/lc_discovery_TIC_35119266_QLP_sector_97_sigma.vot
$PY scripts/step1_select.py data/lc_tess_FFI_sector_30_user_prepared.vot
$PY scripts/step1_select.py data/NSV_807_sector_30_SPOC_pdcsap.vot
$PY scripts/step1_select.py data/NSV_807_sector_97_flatten.vot
$PY scripts/step1_plot.py data/lc_discovery_TIC_35119266_SPOC_sector_30_prepared.vot data/lc_discovery_TIC_35119266_QLP_sector_97_sigma.vot data/lc_tess_FFI_sector_30_user_prepared.vot
$PY scripts/step2_shape.py data/lc_discovery_TIC_35119266_SPOC_sector_30_prepared.vot
$PY scripts/step3_timing.py data/lc_discovery_TIC_35119266_SPOC_sector_30_prepared.vot
$PY scripts/step3_timing.py data/lc_discovery_TIC_35119266_QLP_sector_97_sigma.vot
$PY scripts/step4_sector4.py
$PY scripts/step5_fit_pictures.py
$PY scripts/step6_timing_files.py
$PY scripts/step7_window_test.py
$PY scripts/step8_timing_files_w015.py
```

`step4_sector4.py` does not produce timings: at 1800 s every event is rejected. The sector 4 timing files come from `step6_timing_files.py`.

### 6.4 Degraded spot curves, Mikulášek timings, trend and profile pictures

Both pictures are public functions in `mikulasek_eb.plotting` (matplotlib is optional, so that module is not imported from the package root). `plot_instrumental_trend` is the instrumental trend. `plot_global_fit` is the four-panel profile picture: time series, folded curve, primary zoom, residuals. This script calls both and writes one of each per light curve, under `trends/` and `fits/`. The 120 s run below is the one that wrote both.

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
$PY examples/spot_degraded_timing.py --cadence 120 --window 0.12
```

The same command with `--cadence 600 1800` is the rest. Trend pictures go to `results/synthetic_scoring/spot_degraded/trends/` and profile pictures to `fits/`. The table is `timing_<cadences>_w<window>.csv`. The model structure is chosen on the spot-free 120 s curve and the O'Connell term is kept on; each file is fitted on its own, with its own baseline.

## 7. Timeline (to fill in)
- Abstract deadline: __
- Slides ready for the first review: __
- Repository public and tested from a clean clone: __
- Conference: __
