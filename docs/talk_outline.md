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
| 5 | Spots and the time of minimum | One detached system, one cool spot, eight longitudes plus the spot-free parent, 120 s only. Bisectors at several depths and the lowest sample, against conjunction and against `T_min_pure`. Then Mikulášek, AP and WSAP on the same nine noisy curves, at two noise levels (1 percent, and Tmag 11). Primary and secondary apart. WSL is not on this figure: it is for flat-bottomed eclipses, and this minimum is not flat. The slide states the conclusions of the section "Conclusions: timing a spotted detached binary", including what is not concluded (no ranking, one system, one spot). |
| 5b | Other synthetic timings | Mikulášek and WSAP/WSL on the spot-free and the G1 to G8 set, at 120 s, 600 s and 1800 s where those files exist. Bias, scatter, and whether the quoted error matches the scatter. |
| 6 | Real timings | NSV 807, sectors 4, 30 and 97. Both methods on the same files. |
| 7 | NSV 807 O-C | Per eclipse, primary and secondary, and the change between sectors. The explanation stays a hypothesis. |
| 8 | Two more stars | Same recipe. The stars are not chosen yet. |
| 9 | Flattening | One example: flattening put a wave at the TESS orbital period into the O-C of NSV 807 sector 97. |
| 10 | Availability and limits | Repository, how to run, and the limits. A converged fit is timed only when the verdict is `ok`. |

## Results

### Photometric minimum against conjunction, pure curves

The PHOEBE model of this experiment is a new short-period detached binary, not G6 and not G1 to G8. LC_Generator named the family `SP_`. Circular orbit, P = 0.8731 d (not commensurate with a day), q = 0.95, i = 87 deg, fill 0.32 / 0.30, Teff 5800 / 5600 K, blackbody atmospheres. Eclipse depth of the spot-free twin is 0.541 mag. One cool spot on the primary when present: `relteff` 0.82, radius 32 deg, colatitude 50 deg, longitude 0 to 315 deg in steps of 45. Longitude 0 faces the companion.

The undegraded PHOEBE files are the one-orbit parents. They have no TESS sampling, no noise and no trend. Directory `/home/voz/projects/UPJS/EB_synthetic_lc/results/clean/detached/`. The spot-free parent is:

`/home/voz/projects/UPJS/EB_synthetic_lc/results/clean/detached/detached_none_nooconnell_SP_none.vot`

The eight spotted parents are the same stem with `SP_L000` to `SP_L315` in place of `SP_none`. The 1 percent degraded copies (120 s, 600 s, 1800 s) are in `/home/voz/projects/UPJS/EB_synthetic_lc/results/n01/detached/`. Manifest of that set: `/home/voz/projects/UPJS/EB_synthetic_lc/results/n01/manifest_spotsweep.csv`. The Tmag 11 copies are in `results/m11/` (see the noise section). Every synthetic timing in this section is the 1 percent set.

`T_conj` is the geometric conjunction and does not move when the spot moves. `T_min_pure` is the median of five bisector midpoints (20, 35, 50, 65 and 80 percent of the way from the local flux minimum up to the 90th percentile continuum), measured on the noise-free one-orbit curve.

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

No arguments. The script reads every one-orbit parent `*SP_*.vot` in `/home/voz/projects/UPJS/EB_synthetic_lc/results/clean/detached` and takes `T_conj - T_min_pure` from the matching 120 s file in `results/n01/detached`. It measures both eclipses with `measure_curve_minimum`, `kind="flux"`, window `0.12 * period`.

The measurement itself:

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
$PY -m pytest -q tests/test_curve_minimum.py -p no:cacheprovider
```

### Mikulášek timings of the nine 120 s curves

The curves carry white noise of 1 percent of the flux (not matched to TESS; see the noise section). Each curve chose its own harmonics, O'Connell term and baseline. The solver is Levenberg-Marquardt. The second-component gates are the package defaults. There is no extra weight on the minima. The timing window is the centre ± 0.12 cycles (full width 0.24). Every chosen fit has status `ok`, so every curve is timed. The offset is the Mikulášek time minus the truth, in seconds. A positive value means the template time is later than that truth.

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

The longitude figure of this 1 percent set puts the bisector median, the Mikulášek time, AP and WSAP on one pair of axes, primary and secondary. All four are seconds from conjunction, positive when later. The bars are the scatter of the eclipses. The centroid is not on this figure. WSL is not on it: WSL is for flat-bottomed eclipses, and this spotted detached minimum is not flat.

A window here is the centre ± a half-width. AP and WSAP use the same windows, set on the spot-free one-orbit curve as the half-width at half the full depth: primary ± 0.0176 cycles, secondary ± 0.0178 cycles. They were not copied from one eclipse to the other, and they were not chosen by a scan. Mikulášek uses ± 0.12 cycles (full width 0.24). Both MAVKA methods fitted the flux and used no weights.

With those windows, AP and WSAP follow the primary bisector. At 45 degrees the primary bisector is -71.4 s, AP -72.0 s, WSAP -70.9 s, and Mikulášek -49.2 s. At 315 degrees the primary bisector is +73.9 s, AP +70.9 s, WSAP +70.3 s, and Mikulášek +51.5 s. On the secondary they move the same way and stop short of the bisector: at 135 degrees the bisector is +41.1 s, WSAP +25.1 s, AP +16.3 s, Mikulášek +11.3 s; at 225 degrees the bisector is -44.4 s, WSAP -23.2 s, AP -19.4 s, Mikulášek -8.7 s. The spot-free point is within about 2 s of conjunction for AP, WSAP and Mikulášek. The eclipse scatter is about 15 s for AP, about 8 to 11 s for WSAP, and about 6 to 9 s for Mikulášek.

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

No arguments. The nine paths are the accepted 120 s rows of `/home/voz/projects/UPJS/EB_synthetic_lc/results/n01/manifest_spotsweep.csv`. The command line is also stored in `results/synthetic_scoring/spot_120_select/command.txt`.

### Mikulášek timings of the same nine curves at Tmag 11

Same selection, same solver, same window (centre ± 0.12 cycles), no extra weight on the minima. The light curves are the Tmag 11 copies (`rel_sigma` 0.0024, about 2.6 mmag). Every chosen fit has status `ok` and is timed. A positive offset means later than that truth.

| Longitude | Chosen | Harmonics | O'Connell | Second component | rms (mmag) | Primary vs conj (s) | Primary vs pure (s) | Secondary vs conj (s) | Secondary vs pure (s) |
|---|---|---|---|---|---|---|---|---|---|
| none | h23_b1 | 2, 3 | yes | no | 4.27 | +0.09 | +0.18 | +0.24 | +0.01 |
| 0 | h2_b1 | 2 | yes | yes | 5.07 | +1.98 | -0.39 | +0.24 | +0.00 |
| 45 | h123_b1 | 1, 2, 3 | yes | no | 11.91 | -49.78 | +21.63 | -2.42 | -3.10 |
| 90 | h23_b1 | 2, 3 | yes | no | 4.65 | -4.95 | +16.30 | +2.35 | -16.74 |
| 135 | h123_b1 | 1, 2, 3 | yes | no | 11.28 | +2.24 | +2.63 | +10.48 | -30.63 |
| 180 | h23_b1 | 2, 3 | yes | yes | 4.77 | +0.09 | +0.19 | +0.39 | -5.43 |
| 225 | h123_b1 | 1, 2, 3 | yes | no | 11.48 | -2.08 | -2.29 | -9.47 | +34.88 |
| 270 | h23_b1 | 2, 3 | yes | no | 4.60 | +5.77 | -16.47 | -1.91 | +17.28 |
| 315 | h123_b1 | 1, 2, 3 | yes | no | 11.67 | +50.81 | -23.11 | +2.87 | +3.25 |

The shift against conjunction is the same as at 1 percent noise. At 45 degrees the primary is 50 s earlier than conjunction and 22 s later than the bisector. At 315 degrees it is 51 s later than conjunction and 23 s earlier than the bisector. The secondary at 135 and 225 degrees is again about 10 s from conjunction and about 30 to 35 s from the bisector. Lowering the noise did not move the template.

What did change is the scatter. Each curve still keeps 28 primary times (4 rejected) and 30 secondary times. The scatter about the pure-curve minimum is 1.5 s on the primary and 2.1 s on the secondary, against about 6 s and 9 s at 1 percent noise. The median formal sigma is about 3.1 to 3.7 s. rms(residual / sigma) is about 0.4 to 0.6, so the formal error is larger than the scatter.

The spot-free rms is 4.3 mmag, against 11.4 mmag at 1 percent. The shoulder longitudes stay at 11 to 12 mmag: that residual is the spot wave the template does not follow, and it is now several times the 2.6 mmag noise. Every chosen fit has status `ok` with the rms/noise warning.

The lower noise also changes the chosen model. The O'Connell term is on for all nine curves, including the spot-free curve, longitude 0 and longitude 180, where it was off at 1 percent. Longitude 180 now keeps a second component; its times stay within 0.4 s of conjunction. Longitude 0 still keeps one, and its times stay within 2 s.

Profile and trend figures, 30 of each: `results/synthetic_scoring/spot_120_select_m11/fits/` and `trends/`. The 1 percent figures in `spot_120_select/` are untouched.

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg $PY examples/spot_120_select.py --set m11
```

`--set n01` repeats the 1 percent run. The command line is stored in `results/synthetic_scoring/spot_120_select_m11/command.txt`.

The Tmag 11 longitude figure uses the same bisector, the Mikulášek times above, and AP and WSAP from the repeat run. The windows are unchanged: AP and WSAP primary ± 0.0176 and secondary ± 0.0178, Mikulášek ± 0.12. WSL was not run. Each method accepted 522 eclipses and rejected 36, the same four unobserved events per file as before. No observed eclipse was rejected.

The offsets are the same as at 1 percent noise, to about 1 s. At 45 degrees the primary bisector is -71.4 s, AP -71.0 s, WSAP -70.4 s, and Mikulášek -49.8 s. At 315 degrees the bisector is +73.9 s, AP +72.0 s, WSAP +70.4 s, and Mikulášek +50.8 s. On the secondary they again stop short of the bisector: at 135 degrees the bisector is +41.1 s, WSAP +25.1 s, AP +17.5 s, Mikulášek +10.5 s; at 225 degrees the bisector is -44.4 s, WSAP -23.1 s, AP -18.2 s, Mikulášek -9.5 s. The spot-free point is within 1 s of conjunction for all three.

The scatter fell with the noise. It is about 3.1 to 3.9 s for AP, about 2.0 to 3.0 s for WSAP, and 1.5 s (primary) and 2.1 s (secondary) for Mikulášek. At 1 percent noise those were about 15 s, 8 to 11 s, and 6 to 9 s.

Figure: `results/synthetic_scoring/spot_120_select_m11/offset_vs_longitude.png`
Plotted values: `results/synthetic_scoring/spot_120_select_m11/offset_vs_longitude.csv`
MAVKA tables: `/home/voz/projects/UPJS/Mavka_user_folder/output/detached_120_m11/tables/`. Reply: `agent_mailroom/Mavka_user_to_LC_Modeller_m11.md`.

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg $PY examples/spot_120_offset_figure.py --set m11
```

### Conclusions: timing a spotted detached binary

Basis: one detached system, one cool spot (`relteff` 0.82, radius 32 deg, colatitude 50 deg) at eight longitudes and the spot-free twin, 120 s sampling, two noise levels (1 percent of the flux, and Tmag 11), one noise realisation per file. Truths: the geometric conjunction, and the pure-curve bisector median (the median of five bisector levels). Methods: the Mikulášek template (window centre ± 0.12 cycles), AP and WSAP (primary ± 0.0176, secondary ± 0.0178 cycles). WSL was not used. Figures: `results/synthetic_scoring/spot_120_select/offset_vs_longitude.png` (1 percent) and `results/synthetic_scoring/spot_120_select_m11/offset_vs_longitude.png` (Tmag 11).

What the data support:

1. **Without a spot all three methods return conjunction.** On the spot-free curve the offsets are within 1 s of conjunction at Tmag 11 (within 2 s at 1 percent). The offsets below are therefore not a fault of the fits on a clean curve.

2. **A spot moves the photometric minimum away from conjunction.** On the noise-free curve the bisector median is up to 74 s from conjunction on the primary and 44 s on the secondary. The shift depends on the longitude. At 45 and 315 degrees mostly the primary moves, in opposite directions. At 135 and 225 degrees mostly the secondary moves, in opposite directions. At 90 and 270 degrees both move by about 20 s with opposite signs, so the secondary is displaced from phase 0.5 relative to the primary. The five bisector levels do not agree with each other (up to 74 s), so "the photometric minimum" is itself not one number.

3. **Every method follows the photometric minimum in the same direction, and none returns conjunction on the spotted curves.** Largest offset from conjunction: AP 72.0 s, WSAP 70.4 s, Mikulášek 50.8 s (all Tmag 11). The sign is always that of the bisector shift.

4. **The methods differ in how far they follow.** Fraction of the bisector shift, Tmag 11:
   - primary at 45 and 315 degrees: AP 0.97 to 1.00, WSAP 0.95 to 0.99, Mikulášek 0.69 to 0.70;
   - secondary at 135 and 225 degrees: AP 0.41 to 0.43, WSAP 0.52 to 0.61, Mikulášek 0.21 to 0.26;
   - quadrature (90 and 270 degrees, both eclipses): AP 0.48 to 0.58, WSAP 0.70 to 0.75, Mikulášek 0.10 to 0.26.
   The methods disagree with each other by up to 21 s on the primary (AP against Mikulášek at 315 degrees) and about 14 s on the secondary (WSAP against Mikulášek at 135 and 225 degrees). In no case does a method return the bisector exactly on the secondary.

5. **The offset is a bias, not noise.** Lowering the noise from 1 percent to Tmag 11 left every offset unchanged to about 1.3 s, while the eclipse scatter fell by a factor of about 4 (Mikulášek 6 and 9 s to 1.5 and 2.1 s; WSAP 8 to 11 s to 2 to 3 s; AP 15 s to 3 to 4 s). The Mikulášek median formal sigma is 3 to 4 s at Tmag 11, so a shift of 50 s is about 14 times it. A quoted timing error does not contain this shift.

6. **Residuals and the verdict do not reveal it.** At the shoulder longitudes the Mikulášek fit has rms 11 to 12 mmag against 2.6 mmag noise, status `ok`, and a warning that the rms exceeds the noise. The warning says the template leaves structure. It does not say the time is displaced. A fit that passes the verdict and the quoted error can still be about 50 s from conjunction.

7. **Which truth is the right one depends on the use, and the data do not decide it.** A timing that sits on the photometric minimum is close to the bisector, and a timing that stays nearer to conjunction is closer to the geometric event. Neither is a merit by itself. For an O-C study of an ephemeris or of a period change, a spot that migrates in longitude produces a shift of the same size as the table above in the times, and the primary and the secondary can move in opposite directions.

What is not concluded:

- **No ranking of methods.** The windows differ by a factor of about 7 (Mikulášek ± 0.12, AP and WSAP ± 0.0176 to 0.0178). The difference between the Mikulášek template and the other two is the method and the window together, and they cannot be separated from this data. The Mikulášek template also has per-curve harmonics, which AP and WSAP do not. The comparison has not been repeated with a common window.
- **No statement for other spots or other stars.** One size, one contrast, one orbital period (0.8731 d), one inclination, circular orbit. Smaller or weaker spots, different periods and eccentric orbits are not tested. The size of the shift (up to about 74 s) is for this system and is not transferred to NSV 807 or to any real star.
- **No statement for longer sampling.** The 600 s and 1800 s copies have not been timed for this question.
- **No statement on noise realisations.** There is one realisation per file. The scatter is the scatter of the eclipses within it. A second seed has not been tried.
- **No physical explanation.** That the shift comes from one wing of the minimum being made shallower by the spot is a likely reading and has not been tested.
- **No WSL, and no Gaussian-process estimate.** WSL is for flat-bottomed eclipses. The Gaussian process is not a timing method in this talk.

Sentence for the slide, for the researcher to accept or change: on a detached system with one cool spot, three timing methods, applied to the same noisy curves, all move in the direction of the photometric minimum by up to 50 to 72 s relative to conjunction. They do not move by the same amount, the shift does not shrink with lower noise, and the quoted error and the fit verdict do not show it.

Verdicts of the researcher, 2026-11-04 (these replace the open items above for this chapter; the experiment is not extended):

- **Short conclusion.** The timing methods tested (Mikulášek, AP, WSAP) all give times that lie closer to the photometric minimum than to the geometric conjunction, to different degrees. We do not claim that two methods agree: the fractions of the bisector shift differ (see conclusion 4).
- **Open question.** How to extract the geometric time of minimum from spotted curves is open. It needs more research and literature reading. The common-window test and the second noise seed are dropped from this chapter for now.
- **Fit reliability.** The from-scratch fits are in the next section. A converged fit of the curve is not a geometric time. Both statements stay in the chapter side by side.

### Fit solution on the Tmag 11 spotted curves

This is not a timing test. Each of the 27 Tmag 11 spot-sweep files (spot-free and eight longitudes, at 120, 600 and 1800 s) was fitted from scratch. The harmonic orders, the O'Connell term and the baseline were chosen on that file. The only inputs were the degraded light curve: time, magnitude, error, and the period and epoch stored in the file. No truth flux, no conjunction, no bisector, and no timing. One system, one spot, one noise realisation.

Every chosen fit has status `ok`. None is `bound-active` or `unresolved`. Twelve other candidates failed. None of those twelve was chosen.

Chosen model and rms (mmag). The O'Connell term is on unless the cell says otherwise. `b1` is baseline order 1, `b0` is a constant baseline. Two-component is marked where the gates adopted it.

| Longitude | 120 s | rms | 600 s | rms | 1800 s | rms |
| --- | --- | --- | --- | --- | --- | --- |
| spot-free | (2, 3), b1 | 4.27 | no harmonics, no O'Connell, b1 | 3.40 | no harmonics, no O'Connell, b1 | 2.96 |
| 0 | (2,), b1, two-component | 5.07 | (1, 2, 3), no O'Connell, b1 | 3.91 | (1, 2, 3), no O'Connell, b1 | 3.39 |
| 45 | (1, 2, 3), b1 | 11.91 | (1, 2, 3), b1 | 11.61 | (1, 2, 3), b1 | 11.52 |
| 90 | (2, 3), b1 | 4.65 | (2,), b1 | 3.79 | (2,), b1 | 3.30 |
| 135 | (1, 2, 3), b1 | 11.28 | (1, 2, 3), b1 | 11.11 | (1, 2, 3), b1 | 10.90 |
| 180 | (2, 3), b1, two-component | 4.77 | (1, 2, 3), no O'Connell, b1 | 3.75 | (1, 2, 3), no O'Connell, b0 | 4.81 |
| 225 | (1, 2, 3), b1 | 11.48 | (1, 2, 3), b1 | 11.19 | (1, 2, 3), b1 | 11.08 |
| 270 | (2, 3), b1 | 4.60 | (2,), b1 | 3.77 | (2,), b1 | 3.35 |
| 315 | (1, 2, 3), b1 | 11.67 | (1, 2, 3), b1 | 11.50 | (1, 2, 3), b1 | 11.21 |

What this supports, on this one system:

- The chosen solution converges at all three cadences, including the curves the spot distorts most (45, 135, 225 and 315 degrees). On those four the rms stays near 11 mmag at every cadence. The noise put into the files is 2.6, 1.5 and 1.4 mmag, so that residual is the spot, not the noise. On the spot-free curve and at 0, 90 and 270 degrees the rms falls from about 4.3 to 5.1 mmag at 120 s to about 3.0 to 3.9 mmag at 1800 s.
- The model is chosen per file. The same longitude does not keep the same harmonics when the cadence changes. A second component is adopted only at 120 s, and only at 0 and 180 degrees.
- At 1800 s the needle still appears, and selection does not adopt it. On longitude 0 the candidate with harmonics (2,) and baseline order 1 is `bound-active`: depth about 26000 mag, Gamma about 3.5e-6, full width 2.5e-5 cycles, rms 21 mmag. Two further candidates on that file did not converge. The chosen fit is harmonics (1, 2, 3), status `ok`, rms 3.39 mmag. On longitude 180 the candidate with harmonics (2, 3) and no baseline is the same kind of needle (depth about 30 mag, rms 16.9 mmag), and two candidates did not converge. The chosen fit is `ok`, rms 4.81 mmag, with baseline order 0 because the order-1 copy of that harmonic set did not converge. At 135, 225 and 315 degrees the (2, 3) candidates are `bound-active` (Gamma on its low bound, and a needle at 135 degrees). The chosen (1, 2, 3) fits are `ok`. This is the same needle as in "A converged fit is not yet a time", which locked the harmonics to (2,). Here the lock is not applied, and the needle loses the selection.

- One 120 s candidate is `bound-active` and not chosen: longitude 180, harmonics (2, 3), baseline order 0, because a second component is wider than half a cycle.

What this does not support:

- It does not say the time of minimum is the geometric conjunction. These fits were not timed.
- `ok` means the checks passed. On every 120 s and 600 s chosen fit the residual warning says the rms is more than 1.5 times the noise estimate, so the template leaves structure. That estimate is the fold-residual estimator, which is inflated on these narrow eclipses, so the ratio is not a measure of the injected noise. At 600 s the selection also says the eclipses are poorly resolved (about 4 points per FWHM), and at 1800 s about 1.5 to 1.7 points per FWHM.
- On every file the selection used the number of points as the effective sample size, because the correlation time was not measured. The note on the file says the evidence is overstated.
- One spot, one system, one noise realisation. The 1 percent set was not refitted here.

Table: `results/synthetic_scoring/spot_m11_fit/chosen.csv`. Every candidate: `results/synthetic_scoring/spot_m11_fit/candidates.csv`. Profile figures of every fitted candidate: `results/synthetic_scoring/spot_m11_fit/fits/`.

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg $PY examples/spot_m11_fit.py
```

### Test of the trust-region fallback on the 1800 s spotted curves

Strategy under test. Levenberg-Marquardt runs first, as in every fit above. A candidate whose fit is marked `bound-active` is fitted again with the trust-region solver (depth, width and Gamma kept inside their boxes while it runs), and the verdict of that second fit is the one compared. A candidate that did not converge is not retried. Starting with the trust-region solver stays off, because on XY Boo that made every candidate `bound-active`. The strategy is now the default of `select_model` (`FitPolicy.retry_bound_active`, true).

Test. The nine Tmag 11 files at 1800 s only (spot-free and eight longitudes), the same input as before: the degraded light curve, nothing else. Everything else is the same as in the run above. Compared with the 1800 s column of the table above.

Result.

- The chosen model is the same on all nine files: the same harmonics, O'Connell term, baseline, solver (Levenberg-Marquardt) and rms.
- Seven candidates were retried (each was `bound-active` at 1800 s, on longitudes 0, 135, 180, 225 and 315). None of the nine files needed the retry to have a chosen fit.
- Two of the seven became `ok`: longitude 135, harmonics (2, 3), baseline order 0, at 17.71 mmag (the chosen fit is 10.90 mmag); and longitude 225, harmonics (2, 3), baseline order 1, at 11.08 mmag with a second component (the chosen fit, harmonics (1, 2, 3), is 11.08 mmag). Neither was chosen.
- Five stayed `bound-active`. Longitude 0, harmonics (2,): one component is 0.755 cycles wide (rms 4.14 mmag). Longitude 135 with baseline order 1, 180 with baseline order 0, and 315 with both baseline orders: Gamma sits on the box floor 0.3 (rms 17 to 20 mmag).
- The needles are gone from the retried fits. Depth is no longer 26000 mag, and the full width is no longer 2.5e-5 cycles.
- Candidates that did not converge stayed as they were (four on the 1800 s files: two on longitude 0 and two on longitude 180).

Conclusion.

1. On this test the fallback does what it is meant to do on the candidates it touches: it removes the run-away depth. It does not make those candidates good. Five of seven still fail the verdict, now on the box floor of Gamma or on a component wider than half a cycle.
2. It did not change a single choice. On all nine files another candidate was already `ok` under Levenberg-Marquardt, and it won. So this test does not show that the fallback rescues a file. A file on which every Levenberg-Marquardt candidate is `bound-active` was not in this set, and the strategy is not tested on one. The earlier case where the bounded solver recovered a curve (longitude 180, harmonics locked to (2,), 1800 s, 1 percent noise) is a different test.
3. The recovered `ok` candidates lose on evidence (the 135 degree one also on rms), so the fallback does not push selection towards worse fits here. One of them adopted a second component at 1800 s, which no chosen fit does.
4. The cost is one more fit per `bound-active` candidate. The fallback is safe to leave on, because the verdict still runs on the second fit and `time` still refuses a fit that is not `ok`. This statement is for this system only.
5. The suggestion to users is therefore: keep it on, and report in the selection table which solver produced each candidate. The table `candidates.csv` has a `solver` column for that.

What this test does not support: any statement for other stars, for 120 s and 600 s (no `bound-active` candidate appeared there, except the one at 120 s, longitude 180, which was not refitted in this test), for the 1 percent set, or for real light curves. The selection tests of the package pass with the fallback on (`tests/test_selection.py`, 7 passed, includes XY Boo). The XY Boo selection was not fitted again with the fallback and compared in detail.

Table: `results/synthetic_scoring/spot_m11_fit_c1800/chosen.csv`. Candidates with the solver of each: `results/synthetic_scoring/spot_m11_fit_c1800/candidates.csv`.

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg $PY examples/spot_m11_fit.py --cadence 1800
```

### Fit solution on the other Tmag 11 stars

Fit only, no timing. The 22 accepted files of the Tmag 11 first batch (`results/m11/manifest.csv`): G1, G2, G4, G5, G6, G7 and G8, at 120, 600 and 1800 s. G1 also has a second 600 s file (B22). G6 carries the spot LC_Generator generated, and it is fitted as produced. Each file was fitted from scratch. The only inputs were the degraded light curve. The trust-region fallback was on. No truth column was read.

Every chosen fit has status `ok`. The chosen solution of each file is saved as JSON for a later timing or borrowing run (`ModelSolution.save`). No second component was adopted. The profile figure of each chosen candidate is copied to `results/synthetic_scoring/batch_m11_fit/fits/chosen/` (22 files). The other candidate figures stay in `fits/`.

O'Connell is off unless the cell says otherwise. `b1` is baseline order 1. A blank harmonic set is written "none".

| Star | 120 s | rms | 600 s | rms | 1800 s | rms |
| --- | --- | --- | --- | --- | --- | --- |
| G1, detached | (1,), b1 | 4.13 | none, b1 | 3.60 | none, b1 | 3.47 |
| G2, detached | (2, 3), O'Connell, b1 | 5.36 | (2,), O'Connell, b1 | 4.80 | none, b1 | 4.72 |
| G4, semi-detached | (1, 2, 3), b1 | 20.43 | (1, 2, 3), b1 | 20.50 | (1, 2, 3), b1 | 20.17 |
| G5, contact | (2, 3), b1 | 6.62 | (2, 3), b1 | 6.18 | (1, 2, 3), b1 | 5.69 |
| G6, contact, one spot | (2, 3), O'Connell, b1 | 7.86 | (2, 3), O'Connell, b1 | 7.39 | (1, 2, 3), O'Connell, b1 | 7.03 |
| G7, semi-detached | (1, 2, 3), O'Connell, b1 | 17.25 | (1, 2, 3), b1 | 16.54 | (1, 2, 3), b1 | 14.80 |
| G8, detached, eccentric | none, b1 | 3.98 | none, b1 | 3.40 | none, b1 | 3.31 |

The second G1 file at 600 s (B22) chose the same model as B02, at the same rms of 3.60 mmag.

What this supports, on these files:

- A chosen solution exists at every cadence, including 1800 s, for these seven stars. The rms of G1, G2, G5, G6 and G8 falls as the cadence lengthens. G6 goes from 7.86 mmag at 120 s to 7.03 mmag at 1800 s. The noise put into the files is 2.6, 1.5 and 1.4 mmag.
- G4 stays near 20 mmag at every cadence, and G7 stays near 15 to 17 mmag. That residual does not follow the noise. The same pattern appeared on the spotted detached curves. G6, which does have a spot, does not show it: its rms stays near 7 to 8 mmag.
- On G6 at 1800 s the harmonics (1, 2, 3) with baseline order 1 were `bound-active` under Levenberg-Marquardt (depth 15.9 mag, on the cap). The trust-region retry of that same candidate is `ok` at 7.03 mmag, and it is the chosen fit. The other three candidates on that file were already `ok` under Levenberg-Marquardt, the nearest at 7.11 mmag, so the file had a solution without the fallback. This is the one chosen fit in the batch that the fallback changed.
- Harmonics (1, 2, 3) did not converge on G5 and on G6 at 120 s and at 600 s. Harmonics (2, 3) were chosen there. No other candidate failed. On G5 the (1, 2, 3) fit dies in the next stage after the harmonics are freed: the depth runs to about 9 mag and is cancelled by H1 of about -4 mag, and the solver then exceeds its evaluation limit. There is no finished (1, 2, 3) fit, so it never enters the comparison. The (2, 3) fit has no H1. Its primary Gamma is 0.73 and its secondary Gamma is 0.97, and H2 is 0.15 mag. That is the pointed primary. The finished (1, 2, 3) fit at 1800 s has H1 of only 0.016 mag, so the rounded bottom there is not a large H1. The trust-region fallback does not run on this failure: the candidate raises "did not converge", it is not marked `bound-active`. The batch selection does not skip that fallback because another candidate has already finished. A candidate that raises never receives a verdict, so the retry never starts. Among candidates that do finish with status ok, the choice is the evidence, not the order of the list.

G5 was refitted on its own with that retry turned on for every "did not converge" (`examples/g5_converge_fallback.py`, results in `results/synthetic_scoring/g5_converge_fallback/`). The switch is off in ordinary selection. At 120 s and at 600 s both (1, 2, 3) candidates then finished under the trust-region solver, status ok. The evidence chose `h123_b1`: rms 6.48 mmag at 120 s and 6.04 mmag at 600 s, against 6.62 and 6.18 for the (2, 3) fit that the batch had chosen. At 1800 s nothing needed the retry, and the choice stayed the Levenberg-Marquardt `h123_b1` at 5.69 mmag. The adopted (1, 2, 3) solution has H1 of about 0.017 mag and a primary depth of about 0.30 mag, the same family as the 1800 s fit, not the cancellation at 9 mag and -4 mag. Primary Gamma is 0.77 at 120 s and 0.76 at 600 s; secondary Gamma is 0.91 at both. The curved bottom of the primary is smooth. A residual spike near phase 0 remains. Other profiles were not tried.

The same 600 s file was then fitted again with the trust-region fallback off and Table 2 profiles 2, 3 and 4 added beside profile 1 (`examples/g5_c600_profiles.py`, results in `results/synthetic_scoring/g5_c600_profiles/`). Profiles 5 to 7 stay outside this loop: the staged fit frees Gamma, and those three profiles have none. Every finished candidate keeps its profile figure and a picture of the fitted baseline. The two profile-1 candidates with harmonics (1, 2, 3) again did not converge, and profile 2 with those harmonics finished on the depth bound (about 40 mag) and was left there. The evidence chose profile 3, harmonics (1, 2, 3), linear baseline, rms 5.91 mmag. Profile 3 with harmonics (2, 3) is close, at 5.98 mmag. Profile 4 is worse, about 6.9 mmag. The fitted baseline is a straight line on each side of the gap, with a step of about 0.03 mag at the gap. That is the whole instrumental term the fit removed.

The same file was then searched in three steps instead of all at once (`examples/g5_c600_tiered.py`, results in `results/synthetic_scoring/g5_c600_tiered/`): baseline order with profile 1 and harmonics (2, 3); then profiles 1 to 4 with that baseline; then the larger harmonic set (1, 2, 3) with the winning profile. Each step compares the previous winner with the new candidates. The trust-region fallback was off. The result is the same as the exhaustive run: profile 3, harmonics (1, 2, 3), linear baseline, rms 5.91 mmag. Eight fits were run (six different candidates, the previous winner being fitted again in each later step) and the search took 35 s, against 16 fits and about 5.5 min. The two profile-1 (1, 2, 3) fits that do not converge were never started. The BIC gain of each step was 1332, 247 and 70, so the stop rule (two steps in a row below 10) did not fire. The BIC uses the number of points as the number of independent points, so these gains are overstated; the rms gain of the last step is 0.07 mmag.

The same search was then run on all eight 600 s files of the first batch (`examples/batch_m11_tiered_c600.py`, results in `results/synthetic_scoring/batch_m11_tiered_c600/`). It took about 4 min. The trust-region fallback stayed off, and the early-stop rule did not fire on any file. G1, G2 and G8 keep profile 1 and the same rms as the batch (3.60, 4.80 and 3.40 mmag); the profile step gained nothing there. G5, G6, G4 and G7 move to profile 3 with harmonics (1, 2, 3). The rms becomes 5.91, 7.18, 20.17 and 16.28 mmag, against 6.18, 7.39, 20.50 and 16.54 in the profile-1 batch. G4 and G7 stay far above the noise. A different Table 2 profile does not repair those two light curves.

G4 at 600 s was then used to test how the instrumental trend is removed (`examples/g4_c600_trend.py`, results in `results/synthetic_scoring/g4_c600_trend/`). The package removes a per-segment polynomial in magnitudes, and the light curve is already split at the gap. In the fit the deep primary minima of the first segment lie about 0.15 mag below those of the second, and the folded primary shows two bands. A trend that is additive in flux does not shift a deep minimum and a maximum by the same number of magnitudes, so a magnitude baseline cannot remove it. The test used the observed curve only: an additive polynomial in flux per segment, estimated from the residual of the current model, subtracted from the observed flux, with the package fit on the result, repeated until the rms stopped changing. The profile, the harmonics and the gap split were the same in all four runs (profile 3, harmonics (1, 2, 3)). The rms was 21.4 mmag with a constant per segment, 20.2 mmag with a straight line per segment (the batch value), 12.1 mmag with a straight additive flux trend, and 9.9 mmag with a quadratic one. The primary core rms fell from 72 to 39 and 27 mmag, and the two bands in the folded primary are gone. The remaining large residuals are at the first and last points of each segment, where the steep ingress or egress sits at the edge of the data. The package has no flux-domain baseline, and the back-fitting loop lives in the example only. If this is adopted, it belongs in the package. LC_Generator confirmed afterwards that the instrumental trend he injects is additive in flux. The test therefore matches the generator. It does not show that the trend of a real TESS light curve is additive in flux, and the choice between additive and multiplicative trend must not be fixed from the synthetic files alone.

The real file NSV 807, sector 97 (QLP SAP, background subtracted, 200 s) was tested with its own eclipses (`examples/nsv807_trend_kind.py`, results in `results/nsv807_s97_trend/`). For each of 134 orbital cycles the mean flux was taken at the primary minimum, the secondary minimum and the two quadratures. The flux of a minimum follows the flux of the maxima with a slope of 0.51 ± 0.06 at the primary and 0.61 ± 0.05 at the secondary, and the intercepts are 0.08 ± 0.07 and 0.01 ± 0.06. A multiplicative trend gives a line through the origin whose slope is the minimum-to-maximum ratio (0.58 and 0.62 here). An additive trend gives a slope of 1, which is excluded at about 8 sigma for both minima. The trend of this file is multiplicative to within the errors. The fraction of the variation of the maxima that is additive comes out at -0.17 ± 0.15 (primary) and -0.03 ± 0.13 (secondary), that is, zero within the errors. The test on the synthetic G4 secondary returns a slope of 1.004 ± 0.009, as injected, so the method does distinguish the two cases. G4 primary gives 0.72 ± 0.26 and is not used, because the eclipse is narrow and the window is wide. Limits: one star, one sector, and the maxima vary by 2.4 percent over the sector. Spots and the O'Connell effect change the ratio of the minima to the maxima as well and are not separated from the trend. The flux errors in this file are all NaN, so `load_tess_votable` refuses it; the test reads time and flux directly with `volightcurve`. A request to LC_Generator (2026-11-04, `agent_mailroom/LC_Modeller_to_LC_Generator_trend.md`) asks for the degradation to be repeated with the trend amplitude set to 3 percent of the median out-of-eclipse flux, 80 percent of it multiplicative and 20 percent additive, and the same ramp, gap-step and drift shapes. The present files are not overwritten.

The tiered search was repeated on the eight 600 s files of that new set (`examples/batch_m11_tiered_c600.py --manifest results/mix/manifest.csv`, results in `results/synthetic_scoring/batch_mix_tiered_c600/`). The trust-region fallback stayed off. G4 falls from 20.17 mmag to 5.71 mmag, still profile 3 with harmonics (1, 2, 3), and the two bands in the folded primary are gone. G5 falls from 5.91 to 2.09 mmag and G6 from 7.18 to 3.92 mmag, the same profile and harmonics. G1, G2 and G8 stay on profile 1, at 2.34, 2.67 and 2.32 mmag against 3.60, 4.80 and 3.40 before; the harmonic set is not the same on every file. G7 produced no fit: both baseline candidates, profile 1 with harmonic 2, are `bound-active` because the secondary width is negative, and with the fallback off the search stops. Injected noise at 600 s is 1.5 mmag, so G4 at 5.71 mmag still leaves structure. G7 alone was then repeated with the trust-region fallback on (`results/synthetic_scoring/g7_mix_c600_fallback/`). The two profile-1 candidates finish after the retry, at about 14.6 mmag, and are not chosen. The evidence chooses profile 3, harmonic 2, linear baseline, Levenberg-Marquardt, rms 5.34 mmag. Harmonics (1, 2, 3) on that profile are worse, at 5.79 mmag. On the old additive trend the same star was 16.28 mmag.
- The harmonics are chosen per file and are not the same at every cadence (G1, G2, G5 and G6 change). G4 and G7 keep (1, 2, 3).

What this does not support:

- No timing, and no comparison with another method. The saved solutions are the input for that later step.
- On G7 at 600 s the selection says the odd and even secondary eclipses differ by 3.7 sigma, so the period may be 2P. That is a remark, not a tested period.
- The same caveats as the spot run: `ok` means the checks passed, the rms-over-noise warning uses the fold-residual estimator, and on every file the selection used the number of points as the effective sample size because the correlation time was not measured. At 1800 s the selection says the eclipses are poorly resolved (1.7 to 5.6 points per FWHM).
- G6 is in this result. Borrowing still does not use a spotted light curve. One noise realisation. The 1 percent copies were not fitted.

Table: `results/synthetic_scoring/batch_m11_fit/chosen.csv`. Candidates: `results/synthetic_scoring/batch_m11_fit/candidates.csv`. Solutions: `results/synthetic_scoring/batch_m11_fit/solutions/` (22 JSON files). Chosen profile figures: `results/synthetic_scoring/batch_m11_fit/fits/chosen/`. All candidate figures: `results/synthetic_scoring/batch_m11_fit/fits/`.

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg $PY examples/batch_m11_fit.py
```

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

### Noise: the synthetic curves and the real TESS curves

**What the synthetic curves carry now.** Every degraded file (G1 to G8 and the spot sweep) has independent white Gaussian noise added after the exposure smear and the additive trend. The standard deviation is 1 percent of the flux of that point (`rel_sigma = 0.01`, about 10.9 mmag), the same at 120 s, 600 s and 1800 s. The error column equals that sigma. Seed 1 is shared on the sweep, so differences between longitudes come from the spot. The generator's own note calls 0.01 a round working default set early in the project: it has not been matched to TESS and is not a measured TESS precision. Do not call these files "calibrated TESS noise". Real white noise falls as 1/sqrt(exposure), so at 600 s and 1800 s the synthetic files are far noisier than real data. The nine 120 s spot curves, the G1 to G8 scoring and the selection table were all made with this noise. Label them "1 percent noise". Source: `/home/voz/projects/UPJS/EB_synthetic_lc/config/spoil_grid.yaml` and the reply `agent_mailroom/LC_Generator_to_LC_Modeller_noise.md`.

**What real TESS curves show.** The white noise of nine real light curves of seven stars was measured as the robust sigma of consecutive differences of the residuals about the fold, so the eclipses are removed. The magnitude is the 10th percentile of the magnitude column built from the file's own photometric calibration. The files carry no Tmag, so this is not the TIC magnitude. The reference is the published noise model of 2-min targets, sigma(1 hr) = 50.2 + 97.4 10^(0.2 (T - 10)) + 92.9 10^(0.4 (T - 10)) ppm for the median star (arXiv:2202.03656, quoted in `docs/synthetic_noise.md`), scaled to one hour by the square root of the cadence.

| Star and file | Cadence | Magnitude | Sigma per point (mmag) | Sigma of 1 hr (ppm) | Model, 1 hr (ppm) | Measured / model |
|---|---|---|---|---|---|---|
| AA And S16 | 120 s | 10.63 | 1.99 | 335 | 346 | 0.97 |
| AV Del S81 | 200 s | 10.76 | 1.44 | 312 | 376 | 0.83 |
| XY Boo S50 | 600 s | 9.86 | 0.86 | 325 | 224 | 1.45 |
| EK Com S49 | 600 s | 11.24 | 1.56 | 586 | 513 | 1.14 |
| AA And S16 | 1800 s | 10.66 | 0.96 | 623 | 353 | 1.77 |
| AA And S17 | 1800 s | 10.66 | 1.06 | 693 | 353 | 1.96 |
| NSV 807 S04 (QLP) | 1800 s | 11.08 | 2.13 | 1387 | 460 | 3.01 |
| AR Aur S19 | 120 s | 6.11 | 0.45 | 75 | 69 | 1.09 |
| V477 Lyr S26 | 120 s | 14.56 | 27.2 | 4568 | 7056 | 0.65 |

At 120 s and 600 s the model holds within a factor of 1.5, mostly within 20 percent. AA And at 2 min is on the model. At 1800 s both AA And sectors lie about 1.8 to 2 times above it, and NSV 807 about 3 times. AA And sector 16 has both a 2-min and a 30-min file of the same star: 335 ppm and 623 ppm of 1 hr, so the excess is in the 30-min photometry, not in the star. The catalogue errors are not a noise estimate: measured noise is 1.8 to 2.8 times the file error on three files, and NSV 807 has a file error of 175 mmag against 2.1 mmag measured. AR Aur (6th magnitude, saturated) and V477 Lyr (faint, extrapolated model) are shown for completeness and are not used for the conclusion.

Limits: nine files and seven stars, one star near magnitude 11 in each cadence class. This is a check of a published model, not a calibration. White noise only: red noise and the instrumental trend are not measured here.

**Decision (2026-11-04): one star magnitude, Tmag 11.** The synthetic curves are degraded again with the noise of a Tmag 11 star at each cadence. PHOEBE is not rerun. Only the degradation step is repeated on the existing one-orbit parents. The 1 percent files stay as a labelled "1 percent noise" set.

| Cadence | `rel_sigma` | mmag | Basis |
|---|---|---|---|
| 120 s | 0.00240 | 2.6 | published model at Tmag 11 (438 ppm in 1 hr) times sqrt(3600 / 120) |
| 600 s | 0.00140 | 1.5 | model 1070 ppm, times 1.3 (real 10-min data) |
| 1800 s | 0.00125 | 1.4 | model 620 ppm, times 2 (real 30-min FFI data) |

Why Tmag 11:
- Four of the seven real stars in `data/` lie between 10.6 and 11.2 (AA And 10.6, AV Del 10.8, NSV 807 11.1, EK Com 11.2). The synthetic test is meant to look like the stars the package will time, so it takes their magnitude. It is not claimed to be the typical magnitude of a TESS eclipsing binary.
- At Tmag 11 the measured noise at 120 s agrees with the published model (AA And, ratio 0.97), so the 120 s value rests on a direct check.
- One magnitude for all cadences separates the effect of sampling from the effect of noise. Tmag 11 is also far enough from the bright-star floor that the noise, not a systematic floor, limits the timing.
- The factors 1.3 (600 s) and 2 (1800 s) come from two stars at 600 s and from one star in two sectors at 1800 s, all measured on the real files above. They are judgements from few files and are stated as such. Without them the model would be too optimistic at the longer cadences.

Not covered: red noise and the instrumental trend are not matched to real data. The timing precision at different cadences is not assumed to be constant. The nine 120 s curves have been timed at Tmag 11 (the section above). The 600 s and 1800 s m11 files, and the G1 to G8 m11 files, have not.

Request to LC_Generator: `agent_mailroom/LC_Modeller_to_LC_Generator_tmag11.md`. Spot sweep and first batch, 120, 600 and 1800 s, degradation step only. Replies: `LC_Generator_to_LC_Modeller_tmag11.md` and `LC_Generator_to_LC_Modeller_dirs_cal.md`.

**Where the files are now.** Under `/home/voz/projects/UPJS/EB_synthetic_lc/results/`: `clean/` one-orbit parents; `n01/` the 1 percent set (unscaled flux, bright state about 1.8 ct/s, which is mag 9.4 on the project PhotCal, not 10 as its README says); `m11/` the Tmag 11 set (27 spot files and 22 first-batch files, parents scaled to the same flux). The empty class folders `results/detached/`, `results/contact/` and `results/semidetached/` are leftovers of the move. Manifests are inside each set folder (`n01/manifest.csv`, `n01/manifest_spotsweep.csv`, and the same names in `m11/`). Scripts that reproduce the 1 percent numbers in this document read `n01/`. `examples/spot_pure_minima.py` reads the parents from `clean/detached/`.

**What was checked in the m11 files** (all 49, `examples/check_m11_files.py`, independent of the generator's checker). Every file reads with `volightcurve` and has the PhotCal Pogson zp_mag 10, zp_flux 1 ct/s. The bright state of the noise-free smeared curve is mag 11.000 to 11.002. The 10th percentile of the observed magnitudes is 10.95 to 10.99 (trend and noise make it slightly brighter). The injected noise, `(flux - truth - trend) / trended truth`, has a standard deviation of 0.00238 at 120 s, 0.00140 at 600 s and 0.00124 at 1800 s, against the requested 0.0024, 0.0014 and 0.00125. The error column over the trended truth equals the requested value. In mmag at the bright state it is 2.61, 1.52 and 1.36. The lag-1 autocorrelation is between -0.05 and +0.002, so the noise is white. Time step 120 s, 600 s and 1800 s, 18720 points at 120 s (27 d minus the 1 d gap). The noise realisation is the same in all nine 120 s spot files (correlation 1.0000 between the normalised noise of the eight longitudes and the first file), so differences between longitudes come from the spot.

The noise estimator that was used on the real files (robust sigma of consecutive differences of the residuals about the fold) is NOT a valid check on these files. On the nine 120 s spot files it gives 5.1 to 5.9 mmag against 2.6 injected, and on the G-batch up to 28 mmag at 1800 s. The eclipses are narrow (FWHM about 0.03 d, 1.4 samples at 1800 s) and the fold with its bins cannot follow them, so eclipse structure enters the residuals. The injected noise was therefore verified against the generator's own truth columns, not with that estimator. Consequence still to examine: the noise-over-catalogue and rms-over-noise numbers that the package derives from the same estimator are inflated on this system.

Figure: `results/real_noise/m11_check.png` (left: the synthetic files on the real-noise figure, as the estimator sees them; right: one 120 s curve against truth plus trend). Table: `results/real_noise/m11_check.csv`. Shared-noise test: `results/real_noise/m11_same_noise.txt`.

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg $PY examples/check_m11_files.py
```

Figure: `results/real_noise/noise_vs_cadence.png`
Table: `results/real_noise/noise_vs_cadence.csv`

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg $PY examples/real_noise_vs_cadence.py
```

### Selection on the synthetic TESS-like set

Eighteen PHOEBE light curves, G1 to G8 and variants, with the 1 percent noise described above, from `/home/voz/projects/UPJS/EB_synthetic_lc/results/n01/manifest.csv`.

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
