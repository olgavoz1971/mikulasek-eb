# Talk tickets

Living list of work and of decisions, including the ones we rejected.
A result that can be shown is written in `docs/talk_outline.md`, in the paragraph it belongs to.
This file is where a choice is recorded. The archive `docs/conference_talk_plan.md` is not updated.

Owners: O = Olga Vozyakova, M = LC_Modeller. Priority A must, B should, C if time.

## Open

### A. Spot figure at 120 s

- [x] (A, M) Fit the nine 120 s curves under the current default. Done. Every chosen fit is `ok` and is timed. Numbers, figure paths and the command are in the outline, section "Mikulášek timings of the nine 120 s curves". Script: `examples/spot_120_select.py`. All 30 candidate profile figures and all 30 trend figures are kept.
- [x] (A, Mavka_user) MAVKA timings of the same nine 120 s files arrived 2026-11-03. Reply: `agent_mailroom/Mavka_user_to_LC_Modeller.md`. Tables and figures: `/home/voz/projects/UPJS/Mavka_user_folder/output/detached_120/`. WSAP is the table that enters the spot figure. WSL is kept on disk and is not used for this test.
- [ ] (A, M) Add the WSAP offsets to the longitude figure. Same sign: seconds from conjunction, positive when later. No WSL points.
- [x] (A, M) Figure of offset against longitude. Done. Primary and secondary panels. Bisector median and Mikulášek time, both as seconds from conjunction, positive when later. The centroid is not drawn. Figure `results/synthetic_scoring/spot_120_select/offset_vs_longitude.png`, from `examples/spot_120_offset_figure.py`.
- [ ] (A, O) Wording of the warning: photometric against geometric minima. One spot size and one contrast. The figure is 120 s only.

### A. Comparison that is not the longitude figure

- [ ] (A, O+M) Freeze the protocol: same light curves, truth `T_min_pure`, same sigma convention, same eclipses, statistics (bias, scatter, rms of error over actual error, fraction rejected). A rejected event includes a fit the verdict did not mark `ok`.
- [ ] (A, O+M) Run both methods on G1 to G8 and on `SP_none`, at 120 s, 600 s and 1800 s where the files exist.
- [ ] (A, M) Ask LC_Generator for the hold-out batch S01 to S18 only after the protocol is frozen. Score it once. No tuning afterwards.
- [ ] (B, M) Rerun G6 timing with the O'Connell term on. Rerun G1 and G8 with the empty harmonic set.
- [ ] (B, M) One clean figure of the second component, and one demonstration of exposure integration on an 1800 s curve (smeared against unsmeared). These are general tests.

### A. NSV 807 and the other real stars

- [ ] (A, M) Repeat the NSV 807 timings after the verdict. The numbers in the outline were made before it existed.
- [ ] (A, O+M) Put the sector 4 point on the O-C. The researcher says the sector 4 times are TDB. The file says TCB. Fix the data or the loader check before the point is plotted.
- [ ] (A, M) Cross-sector O-C with the cycle count, primary and secondary separately.
- [ ] (A, M) Flattened against unflattened O-C of sector 97, from the files that already exist. One example for the slide.
- [ ] (A, O+M) Sector 97 wave and the events at the ends of segments. Decide how edges are treated before calling the wave anything.
- [ ] (A, O) WSAP/WSL on sectors 4, 30 and 97, same files as ours.
- [ ] (A, O) Interpretation text: the change of about -280 s, and the different primary and secondary periods, stated as facts. Spots, a light-time effect, and anything else stay hypotheses. No mechanical light-time fit.
- [ ] (B, O) Choose two more stars, different from NSV 807 in sampling or period, with an unflattened light curve.
- [ ] (B, M) Process those two with the same recipe. No tuning by star.

### A. Slides and the public repository

- [ ] (A, O) Pick the title. Fill the abstract blanks: references, GitHub link, the two other stars. The spot sentence can name the 120 s figure when she next edits `docs/Abstract_Second_draft.md`. This outline does not edit that file.
- [ ] (A, M+O) One-page table: method, light curve, window, number of events, median sigma, O-C scatter, rms(O-C/sigma).
- [ ] (A, O) Known limits on the last slide: profiles 5 to 7 only through `fit_global`; the staged pipeline does not free K; no light-time or apsidal fit; selection tuned on these TESS-like curves; the ask-the-human flag has been raised in every case so far; a converged fit is timed only when the verdict is `ok`.
- [ ] (A, O) Decide what to commit and what to push. The verdict, the solver switch and the outline work are uncommitted.
- [ ] (A, M) No personal paths in `examples/` and `working_scripts/`. Scripts run from a clean clone.
- [ ] (A, M, O reviews) README: install, a short example from light curve to timing file, the paper, the scope, and the verdict.
- [ ] (A, O) Licence, CITATION, version number.
- [ ] (A, O) Real-star data: include it, or ship a download script.
- [ ] (A, M) CI that runs the tests, then a run in a fresh environment.
- [ ] (B, M) Calibrate or document the review flag. It has been raised on every scored curve so far.
- [ ] (C) Conference date, talk length and abstract deadline. Priorities move when those are known.

### Held, do not start unless asked

- [ ] (C) A separate Gaussian-process package. Off this talk. Read `docs/gp_extrema_modeller.md` in the host application before any API. The package does not live inside `mikulasek_eb` and does not import the Dash page. Name, licence, and whether the host application switches now or later, are open. scikit-learn as the GP library needs an explicit decision, because the dependency order has no entry for it.
- [ ] (C, M) `min_wall` as a field of `TimingRecipe`. It is still only a ticket in `todo.md`.
- [ ] (C) Full width at half depth inside the solver. Today it is judged after the fit.
- [ ] (C) A linear inner solution for the amplitudes, the cosines and the baselines.
- [ ] (M) Acknowledge the spot-sweep reply to LC_Generator. Nothing new is asked. The mail has not been sent.
- [ ] Overlapping timing windows in `fit_global_shifts`. Left out of the talk.

## Decided

### The talk

- The talk is a warning. A ranking of methods comes only from the common protocol, on the whole set.
- Flattening is one practical example (NSV 807, sector 97), not a study of TESS detrending.
- `docs/Abstract_Second_draft.md` is the researcher's abstract. It is not rewritten from here.
- The Gaussian process is not in the abstract and is not a timing method in this talk.
- Every timing table that enters the outline is stored with the fit figure of each row.
- An interactive look that was not saved as a script is not a number in the outline.

### The spot experiment

- One detached system, not a contact binary. NSV 807 stays a separate real-star example. The base is circular, P = 0.8731 d, q = 0.95, i = 87 deg, fill 0.32 / 0.30, Teff 5800 / 5600 K. Eclipse depth of the spot-free twin is 0.541 mag. Parents `SP_none` and `SP_L000` to `SP_L315`. Spot on the primary, `relteff` 0.82, radius 32, colatitude 50. Manifest `results/manifest_spotsweep.csv` in the generator repository.
- The longitude figure is the nine curves at 120 s.
- The 600 s and 1800 s copies of the same parents are general tests and the comparison that is not about longitude. The two 1800 s curves are the fixture the verdict must refuse.
- Both truths are kept: `T_conj` and `T_min_pure`. `T_min_pure` is the median of the five bisectors, not "the" minimum. The spread of the five levels is part of the result.
- Harmonics for a spotted curve are chosen on that curve. The O'Connell term is on when that curve's maxima require it.
- The timing window of the detached spot curves is a half-width of 0.12 cycles for the Mikulášek times. The NSV 807 window does not set it. The MAVKA WSAP run used its own window, a full phase width of 0.058, frozen on the spot-free file.
- On this spot test the comparison method is WSAP. WSL is for flat-bottomed eclipses. This minimum is not flat, so the WSL table is not a result of the spot test.
- `min_points` 10, `min_each_side` 3, and the contact walls stay as they are. They are not lowered to recover an 1800 s time.
- One spot size and one contrast. The weak-spot variant (`relteff` 0.92) was not generated.
- PHOEBE spots are fixed inside one light curve. Spot evolution is outside this package.

### The fitter

- The default solver is Levenberg-Marquardt, with no bounds while it runs. Written in section 11 of `docs/mikulasek_implementation_notes.md`.
- After every fit the verdict is `ok`, `bound-active` or `unresolved`. `time` raises `UnreliableSolution` unless the status is `ok`. `ok` means the checks passed. It does not mean the fit is good.
- `FitPolicy.bounded_solver` (default false) is the optional second look. It is how `examples/inspect_fit.py` retries the two 1800 s curves. The verdict still runs afterwards.
- The amplitude cap stays at twice the robust peak-to-peak range. It was not raised.
- The second-component gates stay: rms gain at least 0.25 unless `force_two_component` is set, at least two agreeing starts, positive D and Gamma, every relative error of the six parameters below 0.5, and the two amplitudes not a cancelling pair. `two_comp_rel_err_scope` stays `"all"`.
- Extra weight on the minima stays opt-in (`FitPolicy.eclipse_weight`, default 1). It is not part of the spot recipe. The agent does not turn it on unless the researcher asks.
- On XY Boo the paper's harmonic set (2, 3) is the selection the test requires. With Levenberg-Marquardt that test passes.

### NSV 807, already fixed

- Sector 4 sigma is the constant 3.64 mmag taken from the residuals. The file errors of 175 mmag are not used. This is said on the slide.
- Window 0.25 is the one with about honest errors. Window 0.15 is the one with the smallest scatter. Both files exist. Wider windows are not the quoted timing.

## Rejected

### WSL on this spot test

WSL is the flat-bottom method. The spotted system in this test is detached, and the minima are not flat. The WSL times were computed, and they stay in `/home/voz/projects/UPJS/Mavka_user_folder/output/detached_120/tables/WSL.csv` with their figures. They are not plotted and not quoted for this test. WSAP is the method that goes on the longitude figure.

### Bounded solver as the default

Tried as the default, then removed. On one fit of XY Boo, harmonics (2, 3), the same start: Levenberg-Marquardt 0.40 s, 484 residual evaluations, status `ok`, rms 4.23 mmag, A = 0.279, H2 = 0.016. The bounded trust-region solver 0.95 s, 168 evaluations, status `bound-active`, rms 4.01 mmag, A on the cap (0.648), H2 = -0.149. A cosine can cancel a large A, and the curve still looks like the star. With the bounded solver every XY Boo candidate was `bound-active`, so `select_model` chose none. The same switch also exchanged the two primary components of AA And at an unchanged rms. The one gain was the 1800 s curve at 180 degrees, which the bounded solver recovered to `ok`. That is why the switch remains and is not the default.

Numbers: section 11.2 of `docs/mikulasek_implementation_notes.md`. There is no separate script for that one-curve timing. The current `examples/inspect_fit.py` is the Levenberg-Marquardt path plus the optional retry.

### Loosened second-component gates

On XY Boo the chosen single component is V-shaped (Gamma 0.63) and the points at the bottom are rounder. Forcing a second component and leaving the gates at their defaults rejects the trial: one agreeing start, and the relative errors of the wider copy are above 0.5. Raising the relative-error tolerance from 0.5 to 50 and accepting one agreeing start draws the rounder line (rms about 3.03 mmag, status `ok`). The wider copy is not measured (primary only: depth 0.020 ± 0.021, Gamma 0.80 ± 1.04; both minima: the wider D is 0.16 ± 0.26). H2 can jump from 0.016 to 0.10. The verdict does not check whether a parameter is measured, so those two gates are the only check. They stay at 0.5 and at two starts.

Figures, from `examples/inspect_fit.py`:

- `examples/out/inspect_fit/xy_boo_forced_primary.png` (rejected trial; the drawn line is the single component)
- `examples/out/inspect_fit/xy_boo_forced_primary_shown.png` (gates loosened, primary only)
- `examples/out/inspect_fit/xy_boo_forced_both.png`
- `examples/out/inspect_fit/xy_boo_forced_both_shown.png`

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg $PY examples/inspect_fit.py
```

`$PY` is `/home/voz/projects/UPJS/skvo_veb_2/.venv/bin/python`. The same command writes the other figures in that directory.

### Extra weight on the minima, as a way to round XY Boo

Weight 3 inside 0.12 cycles of each minimum, ephemeris and cosines frozen, on the chosen XY Boo model. Gamma moved from 0.63 to 0.69. The unweighted binned residuals in the minima still swing by about ±4 mmag in the same pattern. The rms of the whole curve went from 4.23 mmag to 4.93 mmag. Status `ok`, timed 55 primary and 53 secondary. 1395 of 2932 points already lie inside the window, so the bright parts were not drowning the eclipses. A single Gamma did not become round at the bottom and right on the wings.

Figure: `examples/out/inspect_fit/xy_boo_weighted_minima.png`

Same command as the previous ticket. The subtitle of an earlier copy of this figure said "Two-component primary" because the two curves differed. The fit is one component. The script now passes the subtitle "One component. Minima weight 3; cosines and ephemeris frozen."

### Relative-error gate on the narrower component only

Compared on the 120 s and 600 s spot curves while the bounded solver was the default. At 120 s the two scopes agreed. The shoulder longitudes that adopted a second component then failed the half-cycle width limit. At 600 s the core-only scope newly kept 45 degrees (then the width limit refused it) and 315 degrees (status `ok`). The gate default was not changed. The run took about 28 minutes, wrote no figures and no parameter files, and does not describe the present solver. It is not repeated.

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg $PY experiments/two_component_gate/compare.py
```

Table: `experiments/two_component_gate/out/compare_120_600.csv`

### Harmonics locked to the spot-free curve, as the spot timing

The structure was taken from the spot-free 120 s curve: harmonics (2,), O'Connell term on, window 0.12. Each file was fitted on its own. Primary offsets against conjunction were about -49 s at 45 degrees and +51 s at 315 degrees, at 120 s, and the 600 s run repeated that pattern. At 1800 s every eclipse was rejected (31 primary and 30 secondary on each file) and the timing limits were not lowered. Those times predate the verdict and are not quoted. The nine-curve run that lets each curve choose (outline) still finds -49.23 s at 45 degrees and +51.52 s against conjunction on the chosen models. The offset is not an artefact of locking the harmonics. The pictures in `examples/out/fit_from_scratch/` are an earlier look at the 315 degree file and are not the table in the outline.

120 s table: `results/synthetic_scoring/spot_degraded/timing_120_w0.12.csv`
600 s and 1800 s table: `results/synthetic_scoring/spot_degraded/timing_600_1800_w0.12.csv`
Profile figures: `results/synthetic_scoring/spot_degraded/fits/`
Trend figures: `results/synthetic_scoring/spot_degraded/trends/`
The recorded command line is `results/synthetic_scoring/spot_degraded/command.txt`.

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
$PY examples/spot_degraded_timing.py --cadence 120 --window 0.12
MPLBACKEND=Agg $PY examples/spot_degraded_timing.py --cadence 600 1800 --window 0.12
```

### Adding H1 by hand, and then stopping

On the 120 s curves, harmonics (1, 2) instead of (2,), on the solver in use at that time. The shoulder rms fell (about 22.7 mmag to about 15.9 mmag at 45 degrees). The primary time stayed shifted by about -49 s at 45 degrees and about +51 s at 315 degrees. An even term can lift both shoulders. It cannot move the minimum. H3 and a sine term were not added.

Table: `results/synthetic_scoring/spot_degraded/h1h2/timing_120_w0.12.csv`
Command record: `results/synthetic_scoring/spot_degraded/h1h2/command.txt`

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
$PY examples/spot_degraded_timing.py --cadence 120 --window 0.12 --harmonics 1 2 --arm h1h2
```

### The centroid over a 0.12-cycle window

On the pure curves it follows the out-of-eclipse spot wave and reaches about +-1200 s. It is not a time of minimum at that window. A narrower window, or a continuum taken outside the spot, would be a new measurement.

### A count of how many of the 27 files were `ok`

The count was made before the solver choice was settled, and once again while the bounded solver was the default. Neither count is the present default. It is not quoted.
