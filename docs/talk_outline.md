# Talk outline

Living document for the talk. Title, abstract, slides, and the results that can still be shown.
A new measurement replaces the paragraph it belongs to. A number that can no longer be shown is removed from this file.
The archive of the mixed plan is `docs/conference_talk_plan.md`. It is not updated.
Decisions, including the ones we stopped, are in `docs/talk_tickets.md`.

Conference name, date, length and abstract deadline are not set.

The audience knows Mikulášek (2015, A&A 584, A8). The talk is this public implementation: what it adds, what it does on synthetic and real TESS light curves, and how its timings compare with WSAP/WSL. The synthetic part carries the warning: a spot can move the photometric minimum away from the geometric conjunction, and two methods can then disagree.

Interpreter for every command below:

```bash
PY=/home/voz/projects/UPJS/skvo_veb_2/.venv/bin/python
```

## Title

Not chosen. The current candidates are in `docs/Abstract_Second_draft.md`:

1. The Mikulášek Model for Eclipsing-Binary Timing: Implementation and Comparison with WSAP/WSL
2. Mikulášek and WSAP/WSL eclipse timing compared: an open implementation, synthetic and TESS data
3. The Mikulášek Model for Timing Eclipsing Binaries: Implementation and Comparison

Earlier alternatives, still unused:

1. What is "the" time of minimum? Spots, photometric and geometric minima, and an open implementation of the Mikulášek method
2. Timing eclipses of spotted binaries: an open Mikulášek-method toolkit compared with WSAP/WSL
3. Photometric against geometric minima in TESS eclipsing binaries: a public toolkit for modelling and timing
4. Be aware of the spot: how minima move, and a public toolkit to measure them

## Abstract

The text below is the researcher's draft in `docs/Abstract_Second_draft.md`. This outline does not edit that file. Square brackets are still empty. The spot figure in the slide list is the nine curves at 120 s. The draft does not name that sampling.

> Several well-established methods are available for measuring times of minima of eclipsing binaries, including the Kwee-van Woerden method, the wall-supported asymptotic-parabola (WSAP) and flat-bottom (WSL) methods of Andronov [refs], and the phenomenological model of Mikulášek (2015). We present a public Python implementation (GitHub link) of the Mikulášek model and its application to timing of minima.
>
> One of the features of the implementation is a semi-automatic model-selection procedure. It fits a set of candidate models, ranks them using the Bayesian information criterion (BIC), and flags the result for the user's review when no candidate is clearly preferred or when the data raise other warnings. This provides an automated alternative to the manual model-selection procedure described in the original paper. We present preliminary tests of this feature on synthetic TESS-like light curves.
>
> The implementation also includes several practical extensions and features intended to simplify its application to survey data. These include reusing a model fitted to a high-cadence light curve as a template for fitting lower-cadence observations of the same system, an optional second profile component to describe distorted minima, and accounting for the finite exposure time of long-cadence observations.
>
> We compare the timing results obtained with the Mikulášek model and the WSAP/WSL methods using synthetic light curves generated with PHOEBE and degraded to resemble TESS observations by adding noise and instrumental trends. The test set includes a spotted system, for which we compare the results with both the geometric time of conjunction and the photometric minimum.
>
> We also apply the methods to real TESS light curves of NSV 807 and [two other stars] from multiple sectors. These examples demonstrate the practical application of the methods to real data and provide preliminary O-C diagrams with their associated uncertainties. We discuss cases where each method performs better.

## Slides

About half the time is slides 5 to 8. About a third is slides 1 to 4. Slides 9 and 10 take the rest. Slide 4 is one slide.

| # | Slide | What it says |
|---|---|---|
| 1 | Why another implementation | The method is known. This package fits, selects and times, and is public. |
| 2 | Design | Layered configuration (model, data, policy, timing recipe). A solution file carries a hash and its provenance. Errors are explicit. Exposure is integrated when the caller passes it. |
| 3 | Features | A second profile component on a minimum. Profiles 1 to 7. A model fitted at short sampling reused at long sampling. Timing window, minimum points, walls. |
| 4 | Model selection, preliminary | Candidates in harmonics, baseline order and the O'Connell term, ranked by BIC. The person is asked when the evidence is weak. The test set and the missing hold-out are stated on the slide. |
| 5 | Spots and the time of minimum | One detached system, one cool spot, eight longitudes plus the spot-free parent, 120 s only. Bisectors at several depths and the lowest sample, against conjunction and against `T_min_pure`. Then Mikulášek and WSAP on the same nine noisy curves. Primary and secondary apart. WSL is not on this figure: it is for flat-bottomed eclipses, and this minimum is not flat. |
| 5b | Other synthetic timings | Mikulášek and WSAP/WSL on the spot-free and the G1 to G8 set, at 120 s, 600 s and 1800 s where those files exist. Bias, scatter, and whether the quoted error matches the scatter. |
| 6 | Real timings | NSV 807, sectors 4, 30 and 97. Both methods on the same files. |
| 7 | NSV 807 O-C | Per eclipse, primary and secondary, and the change between sectors. The explanation stays a hypothesis. |
| 8 | Two more stars | Same recipe. The stars are not chosen yet. |
| 9 | Flattening | One example: flattening put a wave at the TESS orbital period into the O-C of NSV 807 sector 97. |
| 10 | Availability and limits | Repository, how to run, and the limits. A converged fit is timed only when the verdict is `ok`. |

## Results

### Photometric minimum against conjunction, pure curves

One detached system, circular, P = 0.8731 d. One cool spot on the primary at longitudes 0 to 315 degrees in steps of 45, plus the spot-free parent. `T_conj` is the geometric conjunction and does not move when the spot moves. `T_min_pure` is the median of five bisector midpoints (20, 35, 50, 65 and 80 percent of the way from the local flux minimum up to the 90th percentile continuum), measured on the noise-free one-orbit curve.

`T_conj - T_min_pure` in seconds. A positive value means the bisector median is earlier than the conjunction.

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

The photometric minimum moves by up to about 74 s on the primary and 44 s on the secondary. At 90 and 270 degrees the two eclipses move by about 20 s with opposite signs. At 45, 135, 225 and 315 degrees mostly one eclipse moves. The eclipse depth changes from 0.547 to 0.621 mag.

Our five bisectors reproduce the generator's `T_min_pure` to 0.00 s at every longitude and both eclipses. The five levels do not agree with each other: the spread is 0.2 s with no spot and up to 74 s (secondary at 225 degrees). The lowest sample barely moves, except the primary at 45 degrees (-156 s) and at 315 degrees (+120 s). The secondary lowest sample stays at +2.1 s at every longitude. The centroid over the 0.12-cycle window reaches about +-1200 s, because the window contains the out-of-eclipse spot wave. It is not used as a time of minimum.

One spot size and one contrast only. The figure is this system at the sampling of the pure curve. It is not a contact binary, and it is not NSV 807.

Figure: `results/synthetic_scoring/spot_pure/pure_minima.png`
Table: `results/synthetic_scoring/spot_pure/pure_minima.csv` (column `generator_conj_minus_min_s` is the table above)

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
$PY examples/spot_pure_minima.py
```

No arguments. The script reads every one-orbit parent `*SP_*.vot` in `/home/voz/projects/UPJS/EB_synthetic_lc/results/detached`, skips names containing `_c120`, `_c600` or `_c1800`, and measures both eclipses with `measure_curve_minimum`, `kind="flux"`, window `0.12 * period`.

The measurement itself:

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
$PY -m pytest -q tests/test_curve_minimum.py -p no:cacheprovider
```

### Mikulášek timings of the nine 120 s curves

Each curve chose its own harmonics, O'Connell term and baseline. The solver is Levenberg-Marquardt. The second-component gates are the package defaults. There is no extra weight on the minima. The timing half-width is 0.12 cycles. Every chosen fit has status `ok`, so every curve is timed. The offset is the Mikulášek time minus the truth, in seconds. A positive value means the template time is later than that truth.

| Longitude | Chosen | Harmonics | O'Connell | Second component | rms (mmag) | Primary vs conj (s) | Primary vs pure (s) | Secondary vs conj (s) | Secondary vs pure (s) |
|---|---|---|---|---|---|---|---|---|---|
| none | h2_b1 | 2 | no | no | 11.37 | +0.79 | +0.89 | +1.07 | +0.84 |
| 0 | h2_b1 | 2 | no | yes | 11.59 | +2.75 | +0.38 | +1.07 | +0.84 |
| 45 | h123_b1 | 1, 2, 3 | yes | no | 15.84 | -49.23 | +22.19 | -1.60 | -2.28 |
| 90 | h2_b1 | 2 | yes | no | 11.50 | -4.29 | +16.96 | +3.17 | -15.91 |
| 135 | h123_b1 | 1, 2, 3 | yes | no | 15.36 | +2.94 | +3.33 | +11.28 | -29.83 |
| 180 | h123_b1 | 1, 2, 3 | no | no | 11.52 | +0.79 | +0.88 | +1.14 | -4.68 |
| 225 | h123_b1 | 1, 2, 3 | yes | no | 15.52 | -1.38 | -1.59 | -8.71 | +35.65 |
| 270 | h2_b1 | 2 | yes | no | 11.50 | +6.49 | -15.75 | -1.08 | +18.11 |
| 315 | h123_b1 | 1, 2, 3 | yes | no | 15.62 | +51.52 | -22.39 | +3.76 | +4.15 |

The spot-free primary and secondary sit within 1 s of both truths. Where the pure-curve bisector moves by tens of seconds, the template moves in the same direction and by less. At 45 degrees the primary is 49 s earlier than conjunction and 22 s later than the bisector median. At 315 degrees it is 52 s later than conjunction and 22 s earlier than the bisector median. The secondary shows the same partial shift at 135 degrees (11 s later than conjunction, 30 s earlier than the bisector median) and at 225 degrees (9 s earlier than conjunction, 36 s later than the bisector median).

Each curve keeps 28 primary times and rejects 4, and keeps 30 secondary times and rejects none. The scatter about the pure-curve minimum is about 6 s on the primary and about 9 s on the secondary. The median formal sigma is about 8 s and 9 s. rms(residual / sigma) is about 0.8 to 0.9.

The shoulder longitudes 45, 135, 225 and 315 degrees choose harmonics 1, 2 and 3, with rms 15 to 16 mmag. The spot-free curve, longitude 0, and the quadrature longitudes choose harmonic 2, with rms about 11.5 mmag. Longitude 0 is the only chosen fit that keeps a second component. Its times stay within 3 s of conjunction.

Profile figures, one per candidate, 30 files: `results/synthetic_scoring/spot_120_select/fits/`
Trend figures, the same names, 30 files: `results/synthetic_scoring/spot_120_select/trends/`
A file is `{light-curve stem}__{candidate}.png`. The chosen figure for 315 degrees is `fits/detached_none_nooconnell_SP_L315_SP_L315_120_c120__h123_b1.png`. The candidate table names every other file.

The longitude figure puts the bisector median and the Mikulášek time on one pair of axes, primary and secondary. Both are seconds from conjunction, positive when later. The bar on the Mikulášek points is the scatter of the eclipses. The centroid is not on this figure. WSAP will be added from Mavka_user's table. WSL will not: WSL is for flat-bottomed eclipses, and this spotted detached minimum is not flat.

Figure: `results/synthetic_scoring/spot_120_select/offset_vs_longitude.png`
Plotted values: `results/synthetic_scoring/spot_120_select/offset_vs_longitude.csv`

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg $PY examples/spot_120_offset_figure.py
```

No arguments. It reads `pure_minima.csv` and `timing_120_w0.12.csv` and does not refit.

Tables:

- `results/synthetic_scoring/spot_120_select/timing_120_w0.12.csv`
- `results/synthetic_scoring/spot_120_select/candidates_120.csv`

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg $PY examples/spot_120_select.py
```

No arguments. The nine paths are the accepted 120 s rows of `/home/voz/projects/UPJS/EB_synthetic_lc/results/manifest_spotsweep.csv`. The command line is also stored in `results/synthetic_scoring/spot_120_select/command.txt`.

### A converged fit is not yet a time

On the two 1800 s curves of the same detached system the default solver can return a needle: a huge depth and a tiny Gamma, with the rms still near 24 mmag, because the spike falls between the samples. The verdict marks that fit `bound-active`. `time` does not return those eclipses.

On 2026-11-03, default solver, harmonics locked to (2,), O'Connell term on, relative-error tolerance 3 (the package default of that gate is 0.5; this example sets 3):

- Longitude 0: `bound-active`. Depth about 38000 mag, Gamma about 10^-6, full width 2.5e-5 cycles. rms 24.08 mmag. Not timed.
- Longitude 180: `bound-active`. Depth about 5400 mag, Gamma about 10^-5, full width 2.5e-5 cycles. The exposure integral at 21 and at 84 samples differs by 256 mag. Not timed.

The optional bounded solver, same configuration:

- Longitude 180: status `ok`, rms 11.56 mmag. `time` returned 0 primary and 30 secondary. The deep minimum lies at phase 0.5 of the catalogue epoch.
- Longitude 0: still `bound-active`, rms 11.58 mmag, because one copy has a full width of 0.80 cycles.

Figures:

- `examples/out/inspect_fit/detached_none_nooconnell_SP_L000_SP_L000_1800_c1800.png`
- `examples/out/inspect_fit/detached_none_nooconnell_SP_L000_SP_L000_1800_c1800_bounded.png`
- `examples/out/inspect_fit/detached_none_nooconnell_SP_L180_SP_L180_1800_c1800.png`
- `examples/out/inspect_fit/detached_none_nooconnell_SP_L180_SP_L180_1800_c1800_bounded.png`

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg $PY examples/inspect_fit.py
```

No arguments. The same run also writes the XY Boo figures under `examples/out/inspect_fit/`. Those XY Boo figures are not a talk result. The check that the two 1800 s default fits are not `ok`, and that `time` refuses them:

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg $PY -m pytest -q tests/test_validity.py::test_runaway_1800_curves_are_not_ok -p no:cacheprovider -m slow
```

The limits and the solver switch are written in section 11 of `docs/mikulasek_implementation_notes.md`.

### Selection on the synthetic TESS-like set

Eighteen PHOEBE light curves, G1 to G8 and variants, from `/home/voz/projects/UPJS/EB_synthetic_lc/results/manifest.csv`.

The selector's regret against the best candidate is about zero. The additive trend dominates the error budget. Shape error is small for detached systems. The second component is chosen where it helps (G4 at 120 s and 600 s, G5 at 120 s, G6 at 120 s) and is not chosen on the other files. Timing is good for detached and contact systems at 120 s and 600 s. The G4 secondary is poor. G6 has a bias that has not been rerun with the O'Connell term on.

The set is the one the rules were tuned on. The hold-out batch S01 to S18 has not been scored. These scripts write tables and no figures.

Tables in `results/synthetic_scoring/`:

- `features_vs_labels.csv`
- `fits_vs_truth.csv`
- `shape_vs_truth.csv`
- `timing_vs_truth.csv`
- `profiles_vs_truth.csv`
- `twocomp_vs_truth.csv`

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

### NSV 807

TIC 35119266, W UMa, P = 0.3389614 d.

Light curves: sector 4 (QLP SAP, 1800 s), sector 30 (SPOC PDCSAP 120 s, and the FFI 600 s), sector 97 (QLP SAP 200 s, no flattening).

Profile 3 describes the rounded contact eclipses best. On sector 30 SPOC the rms is 4.68 mmag for profile 3 and 5.59 mmag for profile 1. The residual pattern is the same in every sector: spikes of 10 to 15 mmag at minimum centres, a 5 mmag antisymmetric wiggle on both sides of each minimum, and a slow 4 mmag wave. It belongs to the template.

Sector 4 with the model borrowed from sector 30: rms 10.8 mmag with the harmonics fixed, 4.63 mmag with H2, H3 and S1 free. The 2018 maxima are inverted against sector 30 (S1 = -0.022 against -0.011). Sector 4 uses a constant sigma of 3.64 mmag taken from the residuals, because the file errors (175 mmag) are unusable.

Timing precision per eclipse: 3 to 5 s (sectors 30 and 97), about 12 to 13 s (sector 4). At window 0.25, rms(O-C/sigma) is 0.7 to 1.1. Against the catalogue ephemeris: sector 30 about +144 s (primary) and +166 s (secondary), sector 97 about -135 s and -139 s. The change from 2020 to 2025 is about -280 s. The sector 4 point is not on this diagram. Inside sector 30 the local periods differ: 0.33896051 d (primary) and 0.33896189 d (secondary), with the same sign at every window. Sector 97 shows a smooth wave of about +-7 s over 10 to 20 d, primary and secondary together. Neither the period difference nor the wave is presented as a discovered orbit.

These timings were made before the verdict existed. They have not been repeated.

Figures, in `/home/voz/projects/UPJS/nsv_807_oc/results/step5/`:

- `S30_SPOC_120s_fit.png`
- `S30_FFI_600s_fit.png`
- `S97_QLP_200s_fit.png`
- `S04_QLP_1800s_borrowed_ampl_epoch_fit.png`
- `S04_QLP_1800s_free_harmonics_fit.png`

Timing files, in `/home/voz/projects/UPJS/nsv_807_oc/results/timings/`:

- `S30_SPOC_120s_primary.txt`, `S30_SPOC_120s_secondary.txt`
- `S30_FFI_600s_primary.txt`, `S30_FFI_600s_secondary.txt`
- `S97_QLP_200s_primary.txt`, `S97_QLP_200s_secondary.txt`
- `S04_QLP_1800s_free_harmonics_primary.txt`, `S04_QLP_1800s_free_harmonics_secondary.txt`

Window 0.15 copies of the same names carry `_w015` in the filename.

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

`prepare_files.py` writes only the FFI copy `data/lc_tess_FFI_sector_30_user_prepared.vot`. `step4_sector4.py` produces no timings: at 1800 s every event is rejected. The sector 4 timing files come from `step6_timing_files.py`.

### Timing window on NSV 807

Sectors 30 and 97. Wider windows (0.35 and 0.5) shrink the formal sigma and increase the real O-C scatter. rms(O-C/sigma) rises to 1.4 to 2.1. Window 0.15 has the smallest scatter. Window 0.25 has about honest errors. The window is reported with every timing. The NSV 807 window is not the window of the detached spot experiment. That experiment uses a half-width of 0.12 cycles.

Report: `/home/voz/projects/UPJS/nsv_807_oc/results/window_test/report.txt`

The report is written by `step7_window_test.py` in the command block above. No figure was saved.
