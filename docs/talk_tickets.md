# Talk tickets

Living list of work and of decisions, including the ones we rejected.
A result that can be shown is written in `docs/talk_outline.md`, in the paragraph it belongs to.
This file is where a choice is recorded. The archive `docs/conference_talk_plan.md` is not updated.

Owners: O = Olga Vozyakova, M = LC_Modeller. Priority A must, B should, C if time.

## Open

### A. Spot figure at 120 s

- [x] (A, M) Fit the nine 120 s curves under the current default. Done. Every chosen fit is `ok` and is timed. Numbers, figure paths and the command are in the outline, section "Mikulášek timings of the nine 120 s curves". Script: `examples/spot_120_select.py`. All 30 candidate profile figures and all 30 trend figures are kept.
- [x] (A, Mavka_user) MAVKA timings replaced 2026-11-04. Reply: `agent_mailroom/Mavka_user_to_LC_Modeller.md`. Tables and figures: `/home/voz/projects/UPJS/Mavka_user_folder/output/detached_120/`. AP and WSAP enter the spot figure. WSL is on disk and is not used for this test. The earlier scan window (full width 0.058) was overwritten and is not quoted.
- [x] (A, M) AP and WSAP are on the longitude figure. Offsets are in `results/synthetic_scoring/spot_120_select/offset_vs_longitude.csv`. WSL is not plotted.
- [x] (A, M) Figure of offset against longitude. Done. Primary and secondary panels. Bisector median, Mikulášek, AP and WSAP, all as seconds from conjunction, positive when later. The centroid is not drawn. WSL is not drawn. Figure `results/synthetic_scoring/spot_120_select/offset_vs_longitude.png`, from `examples/spot_120_offset_figure.py`.
- [ ] (A, O) Wording of the warning: photometric against geometric minima. One spot size and one contrast. The figure is 120 s only.

### A. Noise of the synthetic curves

- [x] (A, M) Measure the real noise against cadence and magnitude. Done. Nine files, seven stars. Numbers, figure and command are in the outline, section "Noise: the synthetic curves and the real TESS curves". Script `examples/real_noise_vs_cadence.py`.
- [x] (A, O) Choose the noise. Tmag 11, one magnitude, per-cadence `rel_sigma` 0.0024 / 0.0014 / 0.00125 at 120 / 600 / 1800 s. The 1 percent set is kept and labelled.
- [ ] (A, LC_Generator) Degradation only, no PHOEBE, spot sweep (27 files) and first batch (B01 to B22) at the three cadences, new names, new manifests, old files untouched. Request `agent_mailroom/LC_Modeller_to_LC_Generator_tmag11.md`. Waiting for the reply and the checker output.
- [x] (A, M) Check the m11 files against the request. Done, 49 files. Injected noise and Tmag 11 scale are as requested. Details and the command are in the outline. One finding: the fold-residual noise estimator is inflated by narrow eclipses (2 times at 120 s on the spot files), so it must not be used as a noise check here.
- [x] (A, M) Folders moved: scripts and the outline now read `results/clean/` (parents), `results/n01/` (1 percent timings in this document) and `results/m11/` (Tmag 11, not yet timed). Mavka_user's tables were made on the old 1 percent files; those files now live in `n01/` and were not overwritten.
- [ ] (B, M) Look at what the inflated noise estimator does to the package's rms-over-noise warning and to the two-component trigger on these narrow eclipses.
- [x] (A, M) Nine 120 s Mikulášek timings at Tmag 11. Done. The shift against conjunction is unchanged (45 degrees -49.8 s, 315 degrees +50.8 s). The scatter falls to 1.5 s and 2.1 s. Table and command in the outline, section "Mikulášek timings of the same nine curves at Tmag 11".
- [x] (A, Mavka_user) AP and WSAP on the nine m11 120 s files. Reply `agent_mailroom/Mavka_user_to_LC_Modeller_m11.md`. Tables `/home/voz/projects/UPJS/Mavka_user_folder/output/detached_120_m11/tables/`. Same windows. WSL not run. The 1 percent directory was not overwritten.
- [x] (A, M) Tmag 11 longitude figure. Offsets agree with the 1 percent figure to about 1 s. Scatter is about 3 to 4 s (AP), 2 to 3 s (WSAP), 1.5 and 2.1 s (Mikulášek). Figure `results/synthetic_scoring/spot_120_select_m11/offset_vs_longitude.png`. The 1 percent figure is a separate file.
- [x] (A, M) From-scratch fit of the 27 Tmag 11 spot-sweep files (spot-free and eight longitudes, at 120, 600 and 1800 s). Done. Every chosen fit is `ok`. Twelve other candidates failed and were not chosen, including the 1800 s needles. Table and command in the outline, section "Fit solution on the Tmag 11 spotted curves". Script `examples/spot_m11_fit.py`.
- [x] (A, M) Trust-region fallback for a `bound-active` candidate, as the default of `select_model` (`FitPolicy.retry_bound_active`). Tested on the nine 1800 s Tmag 11 files only. The chosen model is unchanged on every file. Two retried candidates became `ok` and were not chosen. The needles are gone; five retries stay `bound-active` (Gamma on 0.3, or a component wider than half a cycle). Outline, same section.
- [x] (B, M) Timing precision across cadences is not measured. The cadence run is a fit-solution test. 2026-11-04.

### A. Conclusions of the spot chapter

- [x] (A, M) Conclusions written in the outline, section "Conclusions: timing a spotted detached binary". Seven statements and a list of what is not concluded. No ranking of methods.
- [ ] (A, O) Accept or change the slide sentence at the end of that section, and the warning-slide wording.
- [x] (A, O) 2026-11-04: do not extend the timing comparison. No common-window rerun, no second noise seed, and no timing of the 600 s and 1800 s files in this chapter. The cadence work is the fit solution, recorded under "Noise of the synthetic curves".

### A. Comparison that is not the longitude figure

Order, decided 2026-11-04. The spot chapter stays as it is. The next tests use the other synthetic stars and leave the spotted longitudes out.

- [x] (A, M) Next, and first. Fit solution only, no timing, on the Tmag 11 first batch (`results/m11/manifest.csv`, every cadence). Done. 22 files. G6 is included with the spot LC_Generator generated. Every chosen fit is `ok`. Outline, section "Fit solution on the other Tmag 11 stars". Script `examples/batch_m11_fit.py`.
- [x] (A, M) Save the chosen solution of each of those files with `ModelSolution.save` (JSON, one file per light curve), next to the candidate table. Done. 22 files in `results/synthetic_scoring/batch_m11_fit/solutions/`. Chosen profile figures copied to `fits/chosen/`. A file with no chosen solution did not occur.
- [ ] (A, O+M) Freeze the protocol before any timing is scored: same light curves, truth `T_min_pure`, same sigma convention, same eclipses, statistics (bias, scatter, rms of error over actual error, fraction rejected). A rejected event includes a fit the verdict did not mark `ok`. Spotted longitudes (`SP_L000` to `SP_L315`) are not in this protocol. `SP_none` is spot-free and may be included. Decide that when this test starts.
- [ ] (A, O+M) After the protocol is frozen: time the stored solutions. Mikulášek from the saved solution. The other method is a separate run. Cadences 120 s, 600 s and 1800 s where the files exist. Not started.
- [ ] (A, M) Borrowing, after the fit-solution step. Reference is the saved 120 s solution of the same star. Targets are the 600 s and 1800 s files of that star. Stars: G1, G2, G4, G5, G6, G7, G8. G6 is included, with the spot LC_Generator generated. The detached spot-longitude sweep (`SP_L000` to `SP_L315`) stays out. Questions, to be answered from the run: whether borrowing stays `ok` where a from-scratch 1800 s fit was refused, and the times against the known truth. Not started.
- [ ] (A, M) Ask LC_Generator for the hold-out batch S01 to S18 only after the protocol is frozen. Score it once. No tuning afterwards.
- [ ] (B, M) Rerun G6 timing with the O'Connell term on. Rerun G1 and G8 with the empty harmonic set.
- [ ] (B, M) One clean figure of the second component, and one demonstration of exposure integration on an 1800 s curve (smeared against unsmeared). These are general tests.

Questions parked here, to be answered when that test starts, not before:

- Mavka at 600 s and 1800 s on the non-spot files. Expectation, not a result: the expert windows (± 0.0176 and ± 0.0178 cycles) hold about one sample at 1800 s, so rejection or a poor time is likely, and 600 s is marginal. Undecided: one run at those expert windows only, or a second run with windows widened to the cadence. A widened window would be labelled as a different test. WSL is not assumed.
- Whether `SP_none` enters the non-spot timing and borrowing, or only the first-batch stars. It is spot-free. The spotted longitudes do not.

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

- One detached system, not a contact binary. NSV 807 stays a separate real-star example. The base is circular, P = 0.8731 d, q = 0.95, i = 87 deg, fill 0.32 / 0.30, Teff 5800 / 5600 K. Eclipse depth of the spot-free twin is 0.541 mag. Parents `SP_none` and `SP_L000` to `SP_L315`. Spot on the primary, `relteff` 0.82, radius 32, colatitude 50. Parents `results/clean/detached/`. 1 percent copies `results/n01/`. Tmag 11 copies `results/m11/`. Manifests sit in each set folder (`manifest_spotsweep.csv`).
- The longitude figure is the nine curves at 120 s.
- The 600 s and 1800 s copies of the same parents are general tests and the comparison that is not about longitude. The two 1800 s curves are the fixture the verdict must refuse. That fixture is the earlier locked-harmonic test on the 1 percent files.
- Fit solution on the Tmag 11 spot set (2026-11-04): from-scratch candidate selection on all 27 files (spot-free and eight longitudes, at 120, 600 and 1800 s). The solver sees only the degraded light curve. This run does not time the eclipses and does not read truth columns.
- Next tests (2026-11-04): the spot chapter is closed. The fit-solution step is the Tmag 11 first batch, no timing, with the chosen `ModelSolution` saved. G6 is in that fit, with the spot LC_Generator generated, and G6 is in the borrowing plan. The detached spot-longitude sweep stays out of borrowing. The Mavka window question (expert windows only, or also a cadence-widened run) is parked until that test starts.
- Both truths are kept: `T_conj` and `T_min_pure`. `T_min_pure` is the median of the five bisectors, not "the" minimum. The spread of the five levels is part of the result.
- Harmonics for a spotted curve are chosen on that curve. The O'Connell term is on when that curve's maxima require it.
- A timing window in this test is the centre ± a half-width, in cycles. Mikulášek uses ± 0.12 (full width 0.24). AP and WSAP use the expert windows measured on the spot-free one-orbit curve: primary ± 0.0176, secondary ± 0.0178. The NSV 807 window does not set any of these.
- On this spot test the comparison methods are AP and WSAP. WSL is for flat-bottomed eclipses. This minimum is not flat, so the WSL table is not a result of the spot test.
- `min_points` 10, `min_each_side` 3, and the contact walls stay as they are. They are not lowered to recover an 1800 s time.
- One spot size and one contrast. The weak-spot variant (`relteff` 0.92) was not generated.
- PHOEBE spots are fixed inside one light curve. Spot evolution is outside this package.
- Noise (2026-11-04): the synthetic curves are degraded again at Tmag 11, with `rel_sigma` 0.0024 (120 s), 0.0014 (600 s), 0.00125 (1800 s). The published TESS model gives 0.0024 / 0.00107 / 0.00062. The factors 1 / 1.3 / 2 on top of it are our own measurement on real files and are judgements from few files. Tmag 11 is the magnitude of our real test stars (10.6 to 11.2), not a claim of a typical TESS binary. PHOEBE is not rerun. The 1 percent set is kept, labelled "1 percent noise", and is not plotted on one figure with the new set. The 1 percent value was a round working default and is not called TESS noise.

### The fitter

- The default solver is Levenberg-Marquardt, with no bounds while it runs. Written in section 11 of `docs/mikulasek_implementation_notes.md`.
- After every fit the verdict is `ok`, `bound-active` or `unresolved`. `time` raises `UnreliableSolution` unless the status is `ok`. `ok` means the checks passed. It does not mean the fit is good.
- `select_model` retries a candidate marked `bound-active` with the trust-region solver (`FitPolicy.retry_bound_active`, default true). The verdict of the second fit is the one that is compared. Starting with the trust-region solver stays off.
- `FitPolicy.bounded_solver` (default false) is that second fit, and the optional first-solver switch. `examples/inspect_fit.py` still shows the two steps on a locked configuration. The verdict still runs afterwards.
- The amplitude cap stays at twice the robust peak-to-peak range. It was not raised.
- The second-component gates stay: rms gain at least 0.25 unless `force_two_component` is set, at least two agreeing starts, positive D and Gamma, every relative error of the six parameters below 0.5, and the two amplitudes not a cancelling pair. `two_comp_rel_err_scope` stays `"all"`.
- Extra weight on the minima stays opt-in (`FitPolicy.eclipse_weight`, default 1). It is not part of the spot recipe. The agent does not turn it on unless the researcher asks.
- On XY Boo the paper's harmonic set (2, 3) is the selection the test requires. With Levenberg-Marquardt that test passes.

### NSV 807, already fixed

- Sector 4 sigma is the constant 3.64 mmag taken from the residuals. The file errors of 175 mmag are not used. This is said on the slide.
- Window 0.25 is the one with about honest errors. Window 0.15 is the one with the smallest scatter. Both files exist. Wider windows are not the quoted timing.

## Rejected

### WSL on this spot test

WSL is the flat-bottom method. The spotted system in this test is detached, and the minima are not flat. The WSL times were computed again with the expert windows, and they stay in `/home/voz/projects/UPJS/Mavka_user_folder/output/detached_120/tables/WSL.csv` with their figures. They are not plotted and not quoted for this test. AP and WSAP are the methods on the longitude figure.

### The scanned WSAP window

The first MAVKA run chose a full phase width of 0.058 by minimising the formal uncertainty on the spot-free file. That window was judged too wide for these eclipses and the times were rubbish. The run was repeated and those files were overwritten. The quoted windows are the expert half-widths, centre ± 0.0176 cycles (primary) and ± 0.0178 cycles (secondary).

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
