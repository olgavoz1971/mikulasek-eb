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
| 5 | Spots and the time of minimum | One detached system, one cool spot, eight longitudes plus the spot-free parent, 120 s only. The spot is fixed within each curve and Mikulášek is fitted from scratch on each curve, which is the best case for its adaptation, so the slide says "a fixed spot" and not "spots". Bisectors at several depths and the lowest sample, against conjunction and against `T_min_pure`. Then Mikulášek, AP and WSAP on the same nine noisy curves, at two noise levels (1 percent, and Tmag 11). Primary and secondary apart. WSL is not on this figure: it is for flat-bottomed eclipses, and this minimum is not flat. The slide states the conclusions of the section "Conclusions: timing a spotted detached binary", including what is not concluded (no ranking, one system, one spot). |
| 5b | Other synthetic timings | Mikulášek and WSAP/WSL on the spot-free and the G1 to G8 set, at 120 s, 600 s and 1800 s where those files exist. Bias, scatter, and whether the quoted error matches the scatter. |
| 6 | Real timings | NSV 807, sectors 4, 30 and 97. Both methods on the same files. |
| 7 | NSV 807 O-C | Per eclipse, primary and secondary, and the change between sectors. The explanation stays a hypothesis. |
| 8 | Two more stars | Same recipe. The stars are not chosen yet. |
| 9 | Flattening | One example: flattening put a wave at the TESS orbital period into the O-C of NSV 807 sector 97. |
| 10 | Availability and limits | Repository, how to run, and the limits. A converged fit is timed only when the verdict is `ok`. |

## Results

### How a Mikulášek time is measured

This is the measurement behind every Mikulášek time in this outline. The numbers below use the G4 solution at 120 s (`p4_h123_b1`), because that is the case the questions were asked on. The same steps apply to the other stars. Only the profile, the harmonics and the second component change.

**What the stored model contains.** The function \(F\) is the whole cycle. For an eclipse it is a depth times a profile of the phase from that eclipse's centre. G4 uses profile 4. With \(x = \varphi/D\) the profile is \((1 - |x|^\Gamma)^{3/2}\) while \(|\varphi| < D\), and it is exactly 0 for \(|\varphi| \ge D\). There is no wing. \(D\) is the half-width of that component: the distance from the centre to the phase where the component has reached zero. The primary of G4 also has a second component, with its own depth and its own \(D_{1b} = 0.077\), wider than \(D_1 = 0.041\). The secondary is one component, \(D_2 = 0.075\).

Outside both eclipses every profile term is zero. The curve that remains between the minima is the harmonics, and only those:

\[
H_1\cos(2\pi\theta) + H_2\cos(4\pi\theta) + H_3\cos(6\pi\theta)
\]

with \(H_1 = 0.017\), \(H_2 = 0.031\), \(H_3 = -0.0025\) mag for this G4 solution. The O'Connell term is off. These cosines are part of \(F\). They are not an instrumental trend.

**The instrumental term is not kept.** The global fit, the step that produces the solution, gives every segment of the light curve its own polynomial. For this file the degree is 1 (the `b1` in the solution name) and a new segment starts at a gap longer than 0.5 d. That polynomial is the instrumental term of the fit. Its coefficients are not written into the solution. Timing does not reuse them.

**What timing does with the points.** For eclipse number \(E\), the points with \(|\theta - E|\) inside the half-window are kept and the rest are left out of that eclipse. On the kept points the model is still the whole \(F\), both eclipses and the harmonics. The other eclipse is evaluated too, and it is zero unless the window reaches it.

**What is free, and along which axis.** Nothing in \(F\) is scaled and nothing in \(F\) is refitted. The depth, \(D\), Gamma and the harmonic coefficients stay at the stored values. One shift \(\delta\), in days, moves the whole of \(F\) along the time axis, and the time of minimum is the ephemeris prediction plus \(\delta\). There is no slide along the vertical axis of the eclipse itself: the depth of the dip is the stored depth.

A straight line is fitted on that same window and removed:

\[
y(t) = a + b\,(t - \bar t) + F\!\left(\theta(t) - \delta/P\right)
\]

\(a\) is a vertical offset and \(b\) is a slope across the window. This is the local stand-in for a slow trend, instrumental or otherwise, and it is estimated again for every eclipse. It is not only a vertical shift. A straight line cannot make the cosine curvature, so on a wide window that curvature is still the harmonics inside \(F\), shifted by the same \(\delta\). The lower-level fitter can also scale the depth, \(D\) or Gamma of the eclipse being timed. `time` does not ask for that, and none of the timings in this outline used it.

On the overview light curves the same half-widths are drawn as chords, and the legend names the drawn length. The red chord is centred on phase 0 and runs from \(-D_1\) to \(+D_1\), so its length is \(2 D_1\). The green chord is centred on phase 0.5 and its length is \(2 D_2\). The blue chord is also centred on phase 0. The MAVKA curves stop at a half-window of 3 \(D\), and for the primary that \(D\) is \(D_1\), so the full width of that window is \(2 \times 3 D_1\). The blue chord has that length. One legend, on the first light curve, applies to every row. The phase axis is \(-1\) to \(1\). The figure is described under "MAVKA on the same eclipses".

### Goal and design of the spotted-curve experiment

**Question.** One question only. For a spot that does not move, how far from the geometric conjunction does each timing method land, and how closely does it follow the asymmetric minimum that the spot produces? The aim is to show examples where the methods differ. It is not an error budget for real stars, and it is not a ranking.

**Design.**

- One detached system. The spot-free twin and the same cool spot at eight longitudes, nine curves in all. The spot is fixed in longitude in the model, and the parent is one orbit, so every cycle of a curve has the same shape. The only difference between cycles is the noise.
- Truths: the geometric conjunction, and the pure-curve bisector median (`T_min_pure`).
- AP, WSAP and WSL are local methods. Each eclipse is measured in its own window, and a window knows nothing about the other cycles. WSL is not used here, because the minima are not flat.
- Mikulášek: each curve is fitted from scratch, with its own selection of profile, harmonics, O'Connell term and baseline. The curve of the spot-free twin is not used as a reference for the others. The fitted template is then held fixed, and each eclipse is timed as a shift of that template, in a window of the centre ± 0.12 cycles. Because the cycles are identical, the template takes its shape from that curve, so the fit behaves as a local fit: the shape belongs to this curve.
- A shared template made on the spot-free curve and refitted on the spotted curves was not used. It would have one fixed shape for every cycle, which imitates a moving spot smoothed into an average shape. That is a different question, and it is not answered here.

**What this cannot answer.**

- **A moving spot.** In a real star the spot migrates. The Mikulášek template would then be an average over cycles, and the shifts would also carry the motion of the spot. The result here is therefore the best case for the Mikulášek method, and it is not a typical error. For AP and WSAP each eclipse is measured separately, so a fixed spot should not favour them. That is an expectation, and it has not been tested.
- Other spot sizes, contrasts, periods, eccentric orbits, other stars.
- Which method is better. The windows differ by a factor of about 7, and the Mikulášek template has harmonics that AP and WSAP do not.

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

### Mikulášek timings of the mixed-trend 120 s curves (2026-11-05)

The same nine curves, now from `results/mix/`: Tmag 11 noise and the mixed trend. No fit was run. Each curve was timed from the solution already chosen for it by the full comparison (`results/synthetic_scoring/spot_mix_full/solutions/`). Every chosen template is profile 1 and status `ok`. The window is again the centre ± 0.12 cycles. The stored solution does not keep the fit's clip mask. Four fits had dropped 1 to 9 points out of 18720, so those points are included here. The 1 percent timings and the additive Tmag 11 timings are untouched.

| Longitude | Chosen | Harmonics | O'Connell | Second component | rms (mmag) | Primary vs conj (s) | Primary vs pure (s) | Secondary vs conj (s) | Secondary vs pure (s) |
|---|---|---|---|---|---|---|---|---|---|
| none | p1_h23_b1 | 2, 3 | no | no | 3.30 | +0.10 | +0.19 | +0.36 | +0.12 |
| 0 | p1_h2_b1 | 2 | yes | yes | 3.42 | +2.00 | -0.38 | +0.36 | +0.13 |
| 45 | p1_h123_b1 | 1, 2, 3 | yes | no | 11.37 | -49.86 | +21.55 | -2.28 | -2.96 |
| 90 | p1_h2_b1 | 2 | yes | no | 3.58 | -4.97 | +16.28 | +2.49 | -16.59 |
| 135 | p1_h123_b1 | 1, 2, 3 | yes | no | 10.91 | +2.25 | +2.64 | +10.61 | -30.50 |
| 180 | p1_h23_b1 | 2, 3 | no | yes | 3.33 | +0.10 | +0.19 | +0.53 | -5.29 |
| 225 | p1_h123_b1 | 1, 2, 3 | yes | no | 11.06 | -2.07 | -2.28 | -9.39 | +34.96 |
| 270 | p1_h2_b1 | 2 | yes | no | 3.59 | +5.81 | -16.43 | -1.82 | +17.37 |
| 315 | p1_h123_b1 | 1, 2, 3 | yes | no | 11.24 | +50.91 | -23.01 | +2.98 | +3.37 |

Each curve keeps 28 primary times and rejects 4, and keeps 30 secondary times and rejects none, as before. The offsets against conjunction agree with the additive Tmag 11 run to 0.2 s. The scatter about the pure-curve minimum is 1.4 to 1.5 s on the primary and 1.9 to 2.1 s on the secondary. The median formal sigma is 2.4 to 2.9 s.

AP and WSAP were then run on the same nine files, with the same windows (primary ± 0.0176, secondary ± 0.0178). WSL was not run. Each method accepted 522 eclipses and rejected 36, the same four unobserved events per file. No observed eclipse was rejected. The offsets agree with the additive Tmag 11 run to about 0.2 s. At 45 degrees the primary bisector is -71.4 s, AP -71.0 s, WSAP -70.3 s, and Mikulášek -49.9 s. At 315 degrees the bisector is +73.9 s, AP +72.0 s, WSAP +70.4 s, and Mikulášek +50.9 s. On the secondary they again stop short of the bisector: at 135 degrees the bisector is +41.1 s, WSAP +25.3 s, AP +17.7 s, Mikulášek +10.6 s; at 225 degrees the bisector is -44.4 s, WSAP -23.0 s, AP -18.1 s, Mikulášek -9.4 s. The spot-free point is within 0.5 s of conjunction for all three. The eclipse scatter is about 3.1 to 3.7 s for AP, about 1.9 to 2.8 s for WSAP, and 1.4 to 1.5 s (primary) and 1.9 to 2.1 s (secondary) for Mikulášek. The mixed trend does not move the comparison.

Figure: `results/synthetic_scoring/spot_120_mix_time/offset_vs_longitude.png`
Plotted values: `results/synthetic_scoring/spot_120_mix_time/offset_vs_longitude.csv`
MAVKA tables: `/home/voz/projects/UPJS/Mavka_user_folder/output/detached_120_mix/tables/`. Reply: `agent_mailroom/Mavka_user_to_LC_Modeller_mix.md`.
Mikulášek table: `results/synthetic_scoring/spot_120_mix_time/timing_120_w0.12.csv`.

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg $PY examples/spot_120_mix_time.py
MPLBACKEND=Agg $PY examples/spot_120_offset_figure.py --set mix
```

### Conclusions: timing a spotted detached binary

Basis: one detached system, one cool spot (`relteff` 0.82, radius 32 deg, colatitude 50 deg) at eight longitudes and the spot-free twin, 120 s sampling, two noise levels (1 percent of the flux, and Tmag 11), one noise realisation per file. The spot is fixed within each curve. Each Mikulášek fit was made from scratch on its own curve, so the template could adapt to the spot. This is the best case for the Mikulášek method (see "Goal and design of the spotted-curve experiment"). Truths: the geometric conjunction, and the pure-curve bisector median (the median of five bisector levels). Methods: the Mikulášek template (window centre ± 0.12 cycles), AP and WSAP (primary ± 0.0176, secondary ± 0.0178 cycles). WSL was not used. Figures: `results/synthetic_scoring/spot_120_select/offset_vs_longitude.png` (1 percent), `results/synthetic_scoring/spot_120_select_m11/offset_vs_longitude.png` (Tmag 11, additive trend) and `results/synthetic_scoring/spot_120_mix_time/offset_vs_longitude.png` (Tmag 11, mixed trend). The mixed-trend offsets agree with the additive Tmag 11 offsets to about 0.2 s.

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

7. **Which truth is the right one depends on the use, and the data do not decide it.** A timing that sits on the photometric minimum is close to the bisector, and a timing that stays nearer to conjunction is closer to the geometric event. Neither is a merit by itself. For an O-C study of an ephemeris or of a period change, the table above suggests that a spot which migrates in longitude could shift the times by amounts of this size, with the primary and the secondary moving in opposite directions. This was not tested, because no spot moves in these curves.

What is not concluded:

- **No ranking of methods.** The windows differ by a factor of about 7 (Mikulášek ± 0.12, AP and WSAP ± 0.0176 to 0.0178). The difference between the Mikulášek template and the other two is the method and the window together, and they cannot be separated from this data. The Mikulášek template also has per-curve harmonics, which AP and WSAP do not. The comparison has not been repeated with a common window.
- **No moving spot.** The spot does not move within a curve, and each Mikulášek template is fitted on that curve. The Mikulášek offsets are therefore the best case for adaptation, not a typical error. With a spot that migrates, the template would be an average and the shifts would also carry the motion. That is an expectation, and it has not been tested. For AP and WSAP each eclipse is measured separately, and a fixed spot is not expected to favour them, which is also untested.
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

### The same sweep with the mixed trend (2026-11-04)

The 27 files were fitted again from `results/mix/`, with the full comparison: profiles 1 to 4, the harmonic sets the rules propose, both baseline orders, a trust-region retry after `bound-active` and after "did not converge". The reported solution is the best BIC among status `ok`. No truth column, no timing. Every file has a chosen fit, and every chosen fit is Levenberg-Marquardt. Results: `results/synthetic_scoring/spot_mix_full/`.

| Longitude | 120 s | rms | 600 s | rms | 1800 s | rms |
| --- | --- | --- | --- | --- | --- | --- |
| spot-free | profile 1, (2, 3) | 3.30 | profile 1, (2,) | 2.39 | profile 1, no harmonics | 2.01 |
| 0 | profile 1, (2,) | 3.42 | profile 4, (1, 2, 3) | 2.55 | profile 4, (1, 2, 3) | 1.89 |
| 45 | profile 1, (1, 2, 3) | 11.37 | profile 1, (1, 2, 3) | 11.05 | profile 4, (1, 2, 3) | 10.97 |
| 90 | profile 1, (2,) | 3.58 | profile 1, (2,) | 2.81 | profile 4, (2,) | 2.33 |
| 135 | profile 1, (1, 2, 3) | 10.91 | profile 1, (1, 2, 3) | 10.75 | profile 4, (1, 2, 3) | 10.59 |
| 180 | profile 1, (2, 3) | 3.33 | profile 1, (2, 3) | 2.46 | profile 4, (1, 2, 3) | 1.95 |
| 225 | profile 1, (1, 2, 3) | 11.06 | profile 1, (1, 2, 3) | 10.75 | profile 4, (1, 2, 3) | 10.70 |
| 270 | profile 1, (2,) | 3.59 | profile 1, (2,) | 2.80 | profile 4, (2,) | 2.38 |
| 315 | profile 1, (1, 2, 3) | 11.24 | profile 1, (1, 2, 3) | 11.06 | profile 1, (1, 2, 3) | 10.92 |

O'Connell is on except for the spot-free curve at 600 s and 1800 s, longitude 0 at 600 s and 1800 s, and longitude 180 at 120 s and 600 s. The baseline is order 1 on every chosen fit. At 45, 135, 225 and 315 degrees the rms stays near 11 mmag at every cadence, as before: that residual is the spot. On the other five curves the rms is lower than on the additive-only trend, from about 3.3 to 3.6 mmag at 120 s down to about 1.9 to 2.4 mmag at 1800 s. Profile 4 is the choice on seven of the nine 1800 s files. At 120 s and 600 s the choice is profile 1, except longitude 0 at 600 s.

Postponed. This run stores the chosen solution only. The pictures of the finished candidates are kept. Repeating it so that every candidate solution and the selection protocol are stored, as in the first-batch section below, is postponed. The selection rule would be the same, so the chosen fits in the table would not change.

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

### Tiered search on the mixed trend, 600 s

The same three steps, trust-region fallback off, on the eight 600 s files of `results/mix/` (`results/synthetic_scoring/batch_mix_tiered_c600/`). The old additive-trend run is the comparison.

| Star | Mixed trend | rms (mmag) | Additive trend, same search |
| --- | --- | --- | --- |
| G1 B02 | profile 1, no harmonics | 2.34 | 3.60 |
| G1 B22 | profile 1, harmonic 1 | 2.42 | 3.60, no harmonics |
| G2 | profile 1, harmonics (2, 3) | 2.67 | 4.80, harmonic 2 |
| G4 | profile 3, harmonics (1, 2, 3) | 5.71 | 20.17, same model |
| G5 | profile 3, harmonics (1, 2, 3) | 2.09 | 5.91, same model |
| G6 | profile 3, harmonics (1, 2, 3) | 3.92 | 7.18, same model |
| G8 | profile 1, harmonic 3 | 2.32 | 3.40, no harmonics |

The two bands in the folded primary of G4 are gone. The residual still spikes by about 30 mmag at the bottom of that primary. Injected noise at 600 s is 1.5 mmag, so 5.71 mmag still leaves structure. The early-stop rule did not fire.

G7 produced no fit with the fallback off. Both first-step candidates, profile 1 with harmonic 2, are `bound-active` because the secondary width is negative, and the later steps never run. Repeated with the fallback on (`results/synthetic_scoring/g7_mix_c600_fallback/`, 111 s): those two candidates finish after the trust-region retry at 14.73 and 14.62 mmag and are not chosen. The evidence chooses profile 3, harmonic 2, linear baseline, Levenberg-Marquardt, rms 5.34 mmag. Harmonics (1, 2, 3) on that profile are worse, at 5.79 mmag. On the additive trend the same star was 16.28 mmag. The chosen fit did not itself use the trust-region solver. The retry only let the search continue. The first step alone took 43 s, because both candidates were fitted twice.

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

### Full comparison on the mixed-trend first batch (2026-11-05)

The 22 files of `results/mix/manifest.csv` (G1, G2, G4, G5, G6, G7 and G8, at 120, 600 and 1800 s, plus the second G1 file at 600 s) were fitted with the full comparison. Profiles 1 to 4, the harmonic sets the rules propose, and both baseline orders. Levenberg-Marquardt first, then a trust-region retry after `bound-active` and after "did not converge". The reported solution is the best BIC among status `ok`. Profiles 5 to 7 were not tried. No truth column, no timing. This is the comparison to quote for these files. The tiered 600 s numbers above are the faster experiment, and the 20 mmag and 16 mmag figures are the additive-trend test.

Every file has a chosen fit, and every chosen fit is Levenberg-Marquardt. All 272 candidates finished with status `ok`. Twenty-two of them finished only after the trust-region retry: harmonics (1, 2, 3) on profiles 1 and 2 for G5 and G6, and harmonic 2 on profiles 1 and 2 for G7 at 600 s. None of those twenty-two was chosen. Every file already had an `ok` Levenberg-Marquardt candidate. The baseline is order 1 on every chosen fit. O'Connell is off unless the cell says otherwise.

| Star | 120 s | rms | 600 s | rms | 1800 s | rms |
| --- | --- | --- | --- | --- | --- | --- |
| G1, detached | profile 1, (1, 2, 4) | 3.17 | profile 1, no harmonics | 2.34 | profile 1, no harmonics | 2.18 |
| G2, detached | profile 1, (2, 3, 4), O'Connell | 3.42 | profile 1, (2, 3), O'Connell | 2.67 | profile 1, (2,) | 2.16 |
| G4, semi-detached | profile 4, (1, 2, 3), second component | 5.74 | profile 4, (1, 2, 3), second component | 5.38 | profile 1, (1, 2, 3), second component | 2.79 |
| G5, contact | profile 3, (1, 2, 3) | 3.08 | profile 3, (1, 2, 3) | 2.09 | profile 3, (1, 2, 3) | 1.91 |
| G6, contact, one spot | profile 3, (1, 2, 3), O'Connell | 4.46 | profile 3, (1, 2, 3), O'Connell | 3.92 | profile 3, (1, 2, 3), O'Connell | 3.87 |
| G7, semi-detached | profile 1, (1, 2, 3), second component | 4.44 | profile 3, (2,), second component | 5.34 | profile 1, (1, 2, 3) | 3.55 |
| G8, detached, eccentric | profile 1, (1,) | 3.08 | profile 1, (3,) | 2.32 | profile 1, no harmonics | 2.17 |

The second G1 file at 600 s (B22) chose profile 1, harmonic 1, linear baseline, at 2.42 mmag.

A second component was adopted on G4 at all three cadences and on G7 at 120 s and 600 s. It was not adopted on G7 at 1800 s, nor on G1, G2, G5, G6 or G8. The second component is the package extension that is tried when the rms exceeds 1.5 times the noise estimate. It is not part of Mikulášek 2015.

On G4 at 1800 s, profile 2 with the same harmonics is within 0.0 BIC of profile 1, so the profile is not decided. On G7 at 1800 s, profile 2 is within 7.1 BIC, and at 120 s within 9.8 BIC. The decisive difference used elsewhere is 10. On G7 at 120 s the selection also says the odd and even primary eclipses differ by 4.9 sigma, so the period may be 2P. That is a remark, not a tested period.

At 600 s the full comparison agrees with the tiered search on seven of the eight files, to the reported precision: G1 B02 2.34, G1 B22 2.42, G2 2.67, G5 2.09, G6 3.92, G7 5.34, G8 2.32 mmag, and the same profile and harmonics. G4 moves from profile 3 at 5.71 mmag to profile 4 at 5.38 mmag. G7 reaches profile 3, harmonic 2, under Levenberg-Marquardt inside the full loop. The bound-active profile 1 and 2 candidates are retried and then lose. They do not stop the search.

Against the additive-trend profile-1 batch, the rms falls on every star. G4 goes from about 20 mmag at every cadence to 5.74, 5.38 and 2.79 mmag. G7 goes from 17.25, 16.54 and 14.80 mmag to 4.44, 5.34 and 3.55 mmag. G5 goes from 6.62, 6.18 and 5.69 mmag to 3.08, 2.09 and 1.91 mmag. G6 goes from 7.86, 7.39 and 7.03 mmag to 4.46, 3.92 and 3.87 mmag.

The noise put into the files is 2.6, 1.5 and 1.4 mmag. Every chosen rms is above that. The nearest is G5 at 1800 s, 1.91 mmag against 1.4 mmag. G1, G2 and G8 at 1800 s sit near 2.2 mmag. G4 at 120 s and 600 s (5.74 and 5.38 mmag) and G7 at 600 s (5.34 mmag) still leave several times the noise. G6 stays near 4 mmag at 600 s and 1800 s; that file carries a spot. The package warning compares the rms with its own fold-residual estimate, which is larger than the injected noise on the narrow eclipses, so a file can sit above the injected noise and still not carry that warning.

What this does not support:

- No Mikulášek timing. The saved solutions are the input for that later step, which waits on the window discussion in `docs/talk_tickets.md`. The MAVKA times of these files are in the next section.
- Profiles 5 to 7 were not in the loop. One noise realisation. The 1 percent copies were not fitted. G6 is in this result, and borrowing still does not use a spotted light curve.
- The selection used the number of points as the effective sample size, because the correlation time was not measured. At 1800 s it says the eclipses are poorly resolved (1.7 to 5.5 points per FWHM).
- The retry after "did not converge" was on for this run. It is not the ordinary default of `select_model`.

Table: `results/synthetic_scoring/batch_mix_full/chosen.csv`. Candidates: `results/synthetic_scoring/batch_mix_full/candidates.csv`. In `candidates.csv`, `rms_fit_mmag` is the rms on the points that candidate kept, and `rms_common_mmag` is the rms on the points common to every status `ok` candidate, which is the rms that entered the BIC. Every finished candidate has a solution JSON beside its pictures, `solutions/{light curve}__{candidate}.json`. `selection/{light curve}.json` is the decision record: the best name, the comparison table, the reasons, the remarks and `n_eff`. The `solution` column of `chosen.csv` points at the winner's file. A candidate that raises has a protocol line and no solution. The light curve is not copied.

This archive is kept so the choice can be repeated without another fit. It can be used to tune the tiered search, by replaying a shorter sequence of candidates against the winner already stored here. It can also be used to try another selection rule, for example an rms measured around the minima, or a limit on which harmonics are allowed. The BIC choice reported in the table stays the record of this run.

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg $PY examples/spot_mix_full.py \
  --manifest /home/voz/projects/UPJS/EB_synthetic_lc/results/mix/manifest.csv \
  --out results/synthetic_scoring/batch_mix_full \
  --expected 22
```

### MAVKA timings of the mixed-trend first batch (2026-11-05)

AP, WSAP and WSL were run on the same 22 files. Flux, no weights. The period and the epoch are the file parameters, held fixed. Primary and secondary follow the kind stored on the file. The pure parents were used only to choose the window. Mikulášek has not been timed on this batch.

The window is one primary half-width and one secondary half-width per star, the same for all three methods and for every cadence of that star. The half-width is half the phase span at half the depth of that valley on the pure curve. The top of the depth is the fainter of the two maxima that separate the valley from the other one. G7 was then set to 0.0400 and 0.0400. Where an observed eclipse's full window held fewer than 6 points, that eclipse's half-width was increased by one cadence in phase and counted again, until it held 6 points. The other kind kept its adopted half-width in that test. Growth stopped before the two half-widths would have added to more than 0.5. No eclipse reached that stop. This widening was requested, so a time from a widened window is part of the result and is marked by `n_steps` in the tables. At 120 s only the adopted width was used. At 600 s the only increases are 11 G7 primaries and 11 G7 secondaries, each by one step (AP widened 7 of those primaries and all 11 secondaries; the other accepted G7 times at 600 s kept the adopted width). At 1800 s, G1, G2, G5, G6, G7 and G8 needed 1 to 4 steps on some or all of their observed eclipses. G4 needed none. The largest increase is 4 steps, on one G1 primary.

The PHOEBE parameters below are the ones stored on the light curves. Fill is the stellar radius divided by its Roche lobe, except the contact stars, where the fill-out factor is the PHOEBE value. The eclipse depth is the deeper eclipse. `conj - pure` is the geometric conjunction minus the pure-curve photometric minimum, in seconds, primary then secondary.

| star | class | P (d) | q | i (deg) | e | fill | depth (mag) | FWHM pri / sec (d) | conj - pure (s) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| G1 | detached | 2.8 | 0.95 | 87 | 0 | 0.28 / 0.26 | 0.507 | 0.0866 / 0.0875 | +0.59 / +0.55 |
| G2 | detached, total | 2.1 | 1.00 | 89.5 | 0 | 0.36 / 0.36 | 0.743 | 0.0786 / 0.0786 | +0.23 / +0.23 |
| G4 | semi-detached | 1.85 | 0.35 | 85 | 0 | 0.45 / 1 | 1.648 | 0.155 / 0.167 | -0.05 / -1.47 |
| G5 | contact | 0.34 | 0.40 | 80 | 0 | fill-out 0.40 | 0.635 | 0.0597 / 0.0626 | -0.04 / -0.11 |
| G6 | contact, spot | 0.36 | 0.45 | 79 | 0 | fill-out 0.38 | 0.672 | 0.0665 / 0.0697 | +170 / -189 |
| G7 | semi-detached | 0.5037 | 0.30 | 84 | 0 | 0.22 / 1 | 1.291 | 0.0415 / 0.0542 | +0.04 / +0.99 |
| G8 | detached | 3.2 | 0.80 | 86 | 0.2 | 0.30 / 0.28 | 0.523 | 0.0725 / 0.107 | +0.37 / -0.04 |

Temperatures are 5800/5600 K (G1), 6200/6150 (G2), 8500/4800 (G4), 5800/5550 (G5), 5850/5600 (G6), 9000/4300 (G7) and 6500/5500 (G8). G6 has one cool spot on the primary: longitude 90 deg, colatitude 50 deg, radius 32 deg, relative temperature 0.82, and the O'Connell flag is set. G8 has argument of periastron 90 deg. On these G8 files the stored secondary conjunction still falls half a period after the primary, to 2e-15 d, so the eccentricity does not move the secondary off phase 0.5 in the metadata that was timed against.

Half of the file FWHM, expressed in phase, sits next to the adopted half-width:

| star | adopted half-width, pri / sec | half the FWHM, in phase |
| --- | --- | --- |
| G1 | 0.0155 / 0.0156 | 0.0155 / 0.0156 |
| G2 | 0.0187 / 0.0187 | 0.0187 / 0.0187 |
| G4 | 0.0421 / 0.0454 | 0.0420 / 0.0451 |
| G5 | 0.0890 / 0.0934 | 0.0878 / 0.0920 |
| G6 | 0.0790 / 0.0814 | 0.0923 / 0.0969 |
| G7 | 0.0400 / 0.0400 | 0.0412 / 0.0538 |
| G8 | 0.0113 / 0.0168 | 0.0113 / 0.0167 |

On G1, G2 and G8 the adopted half-width and half the FWHM agree to 0.0001 in phase. G4 is within 0.0003 and G5 within 0.0014. G6 is narrower than half its FWHM, on both eclipses. G7 was set by hand: 0.0400 is near half the primary FWHM and narrower than half the secondary FWHM.

An offset is the accepted time minus the matched conjunction of the same kind, on the file clock. The same match against `eclipse_T_min_pure_days` gives the offset from the pure-curve minimum. On G1, G2, G4, G5, G7 and G8 those two references differ by less than 1.5 s, so a median quoted from conjunction is the same, to that tolerance, as a median from the pure minimum. On G6 they differ by 170 s and 189 s, and both are quoted. The window centre is not the reference. On B22 the file period is long by enough that the centre walks off the conjunction by up to 0.00266 d (230 s) at secondary cycle 9. That difference is `ephemeris_minus_conj_d`. The offsets below are against the conjunction, so the catalogue period error is not inside them. B22 at 600 s follows the other G1 file at 600 s: the same primary cycle 5, and medians within a few seconds on the eclipses that method accepted.

At 120 s, median of the accepted times, in seconds from conjunction. The standard deviation is in parentheses. Primary, then secondary.

| star | AP | WSAP | WSL |
| --- | --- | --- | --- |
| G1 primary | -1.5 (7.1), 8 kept | -0.1 (52), 9 kept | -0.6 (3.6), 8 kept |
| G1 secondary | -2.7 (5.3) | +0.1 (4.5) | +0.3 (4.4) |
| G2 primary | -3.6 (4.9) | -1.1 (2.0) | -1.2 (1.9) |
| G2 secondary | -0.1 (3.7) | +1.1 (2.7) | +1.0 (2.7) |
| G4 primary | 0.0 (4.2) | -0.4 (1.7) | -0.4 (1.5) |
| G4 secondary | -1.0 (31) | -1.0 (29) | -0.6 (29) |
| G5 primary | +0.7 (4.2) | +0.0 (2.9) | +0.1 (2.9) |
| G5 secondary | -0.6 (4.5) | -0.2 (3.8) | -0.2 (3.8) |
| G6 primary | -97 (4.7) | -148 (3.2) | -148 (3.2) |
| G6 secondary | +124 (5.6) | +170 (4.1) | +170 (4.1) |
| G7 primary | +2.2 (8.4), 38 kept | -0.4 (1.0) | -0.3 (0.8) |
| G7 secondary | +0.1 (16) | -0.9 (14) | -0.5 (13) |
| G8 primary | -1.5 (5.9) | -0.4 (3.2) | -0.3 (3.2) |
| G8 secondary | -3.0 (19) | -2.4 (14) | -2.3 (13) |

Counts are the accepted eclipses. Where the count is omitted, every observed eclipse of that kind was accepted. G1 has 9 observed primaries and 10 observed secondaries; G7 has 51 and 52. The standard deviation on G1 WSAP primary is one eclipse, cycle 5, at +164 s. The same cycle is rejected by AP (extremum outside the parabolic part) and by WSL (the formal uncertainty is not finite). It is rejected again on both 600 s copies of G1, and WSAP accepts it there at -888 s, which is the whole of the 279 s standard deviation on that primary. At 1800 s the same cycle was widened by 4 steps and accepted by AP at +7900 s and by WSAP at +8700 s. WSL still rejects it. The median of the accepted G1 primaries stays within 10 s at every cadence, because one wild time does not move a median. On G4 the secondary standard deviation of about 30 s is mostly cycle 14, which lies about 90 s from conjunction at 120 s and at 600 s. The median stays near 1 s.

On G6 the pure-curve minimum lies 170 s before conjunction on the primary and 189 s after it on the secondary. At 120 s and at 600 s, where the adopted window was used, WSAP and WSL land about 20 s from that pure minimum (primary median +21 s and +23 s from it; secondary -18 s and -21 s) and 146 to 170 s from conjunction. AP lands between the two references: about 97 s and 124 s from conjunction, and about 73 s and 65 s from the pure minimum. The scatter on G6 at those cadences is 3 to 10 s.

At 600 s the other medians stay within a few seconds of the 120 s table. G7 is the star whose window grew, by one cadence, on 11 primaries and 11 secondaries. WSAP medians there stay at -0.2 s and -2.8 s. WSL primary median is -0.7 s, with the standard deviation up to 15 s. AP accepts 40 primaries, median +4.4 s, standard deviation 108 s.

At 1800 s the widened window is the one that was used, except on G4, which stayed at the adopted width, and on 2 of the 8 G8 secondaries. Medians of the unspotted primaries stay within about 10 s of conjunction, and G7 moves to -12 s (WSAP) and -17 s (WSL). The standard deviations grow into tens of seconds, and into about 100 s where many eclipses were rejected (AP on G7: 20 primaries and 15 secondaries). On G6, widened by 2 steps, the WSAP medians are -346 s and +382 s from conjunction. That is 176 s and 193 s beyond the pure-curve minimum, further from conjunction than the minimum itself, with a standard deviation of 32 s. WSL medians are -326 s and +391 s, with standard deviations 134 s and 190 s. AP medians are -81 s and +93 s, with a standard deviation of about 63 s, and one primary rejected.

Each method has 1553 rows. Accepted / rejected: AP 1414 / 139, WSAP 1477 / 76, WSL 1463 / 90. Of the rejections, 76 on every method are `eclipse_observed is false`, the same cycles at every cadence of a star. The remaining AP rejections are 36 that reached the iteration limit, 24 with a non-finite formal uncertainty, and 3 with the extremum outside the parabolic part (G1 primary cycle 5 on the three files that were not the 1800 s file). The remaining WSL rejections are 8 at the iteration limit and 6 with a non-finite formal uncertainty, including that same G1 cycle on all four G1 files. WSAP returned a finite time for every observed eclipse. An empty time is a rejection. It is not a time of zero.

What this section supports is the location of these three methods relative to the conjunction and the pure-curve minimum, on this one noise realisation, with the window rule above. The per-file medians are `results/synthetic_scoring/batch_mix_mavka/offsets.csv`. The parameters are `results/synthetic_scoring/batch_mix_mavka/physics.csv`. The times and the failure reasons are `Mavka_user_folder/output/batch_mix/tables/{AP,WSAP,WSL}.csv`.

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
$PY examples/batch_mix_mavka_offsets.py
```

### Mikulášek half-window scan (2026-11-05)

The batch has not been timed with Mikulášek. This scan asks how the half-window changes the scatter of the times, before a width is chosen.

**Setup.**

- One noise realisation of the mixed-trend first batch. G1, G2, G4, G5 and G7 at 120 s, 600 s and 1800 s. G6 and G8 were added at 120 s only, on 2026-11-06, from the solutions already chosen for them. No fit was repeated.
- The template is the chosen solution of that file. D1 and D2 are the widths of that solution, so D is not the same at every cadence and not the same for the two eclipses.
- No truth column is read. O is the accepted time. C is the file ephemeris, epoch plus E periods, on the file clock. E is an integer for a primary and a half-integer for a secondary. The rms is about the mean of the timed values.
- The grid step is D/10, from D/10 up to 0.25 cycles. That is the widest half-window `time` allows without overlapping the neighbouring minimum.
- The rms at every width uses one fixed set of eclipses, chosen once at the widest window: samples reach both ends, and no gap inside the window is longer than 2 cadences. A cycle cut by the data gap, or by the start or end of the series, stays out of every width. Where a narrower window times only part of that set, the rms is of that part.
- Corrected 2026-11-06. The coverage test had added a further half period to every secondary, so those flags were taken on the following primary. The stored O-C already used epoch plus E periods and was not refitted. The secondary flags, and the rms that depends on them, were recomputed. The primary set did not change. G1 did not change at all. The numbers below are the corrected ones.
- `clip_outliers` is off, so the outlier cut cannot remove an eclipse. The 10-point minimum and the contact-wall check stay at the package defaults. They are what reject G5 and G7 at 1800 s.
- The hypothesis, stated before the run: the widest window should do best, because the template follows the eclipse past the half-depth points, where AP, WSAP and WSL do not. The scan runs out to the 0.25 cap so that this can be seen.

**Why the Mikulášek curves have different lengths.** Every eclipse is scanned over the same phase interval, from D/10 to 0.25 cycles. The step is a tenth of that eclipse's own D, so the samples are evenly spaced in units of D. The right-hand end is not a chosen number of D. It is 0.25 cycles divided by D, and D is different for each eclipse, so the same phase cap falls at a different place on an axis drawn in units of D. A narrow eclipse runs to many D. A wide eclipse stops after a few. This cap was fixed before the rms was known: it is the widest half-window `time` accepts without one eclipse's window meeting the next, and the scan was taken out to it to test the widest window. The 4 D proposal was read off the result afterwards. On the comparison figures the cap is a dash-dot vertical line, labelled "0.25 cycles", in the same way the MAVKA cutoff is a dotted line labelled "3 D". At 120 s the cap falls at:

| Star | Primary | Secondary |
| --- | --- | --- |
| G1 | 18.7 D | 18.7 D |
| G2 | 13.4 D | 13.3 D |
| G4 | 6.0 D | 3.3 D |
| G5 | 5.3 D | 5.2 D |
| G6 | 5.4 D | 5.3 D |
| G7 | 17.0 D | 7.4 D |
| G8 | 26.2 D | 14.4 D |

A window that timed nothing is not drawn, so a curve can also start later than 0.1 D. The G7 primary starts at about 2.1 D. The different lengths across the figure are the right-hand ends.

**MAVKA's window on the same axis.** The 120 s figures were redrawn from the saved tables. No timing was repeated. The dashed line is the adopted MAVKA half-width of that star and that eclipse, the half-depth window used for the batch, divided by the Mikulášek D of the same eclipse. Primary and secondary are separate, because both the adopted half-width and D differ.

| Star | MAVKA half-width, pri / sec | In units of D |
| --- | --- | --- |
| G1 | 0.0155 / 0.0156 | 1.16 / 1.17 |
| G2 | 0.0187 / 0.0187 | 1.00 / 1.00 |
| G4 | 0.0421 / 0.0454 | 1.02 / 0.60 |
| G5 | 0.0890 / 0.0934 | 1.91 / 1.96 |
| G6 | 0.0790 / 0.0814 | 1.71 / 1.73 |
| G7 | 0.0400 / 0.0400 | 2.72 / 1.19 |
| G8 | 0.0113 / 0.0168 | 1.19 / 0.97 |

On G1, G2 and G8 the MAVKA window sits at about 1 D, on the steep part. On the contact stars it is near 2 D, where the rms is still falling. G4's secondary is 0.60 D, because D2 is 0.075 and the adopted half-width is 0.0454. G7's primary is 2.72 D, because D1 is 0.015. The marked figures are `results/synthetic_scoring/window_scan_mik/{star}_c120/rms_vs_half_width_mavka.png`.

At 120 s the rms falls steeply up to about 3 D or 4 D and then changes little on the detached stars. G1 goes from about 65 s (primary) and 200 s (secondary) at 0.5 D to about 2.2 s at 0.25. G2 levels near 4 D (primary 1.4 s, secondary 2.6 s) and the widest window is no better. G4's primary is already about 1.2 s from 2 D; its secondary is still falling at the cap, 13 s, because 0.25 cycles is only 3.3 D. G5 keeps falling all the way: primary 27 s at 1 D and 2.9 s at 0.25, secondary 70 s and 3.3 s. G7's primary falls from 1.7 s at 3 D to 0.63 s at 0.25; the secondary is still at 11 s, with the mean near -2.4 s. G8's primary falls from 12 s at 1 D to 2.1 s at the cap. Its secondary falls from 107 s at 1 D to 11 s at the cap and is still near 11 s, with the mean near -8 s. G6, the spotted contact star, keeps falling: primary 25 s at 1 D, 4.1 s at 3 D, 2.9 s at the cap; secondary 44 s, 5.0 s and 3.3 s. On G6 the mean also moves as the window widens, from -25 s to -47 s on the primary and from +25 s to +39 s on the secondary, so a wider window changes the bias and not only the scatter. Beyond a few D, rms(O-C)/sigma is about 0.5 on G1, so the formal sigma is about twice the scatter, and about 1 on G5 and G6, so there the formal error matches the scatter.

At 600 s the scatter is larger and is still falling, or has just levelled, at 0.25. G1 primary goes from 22 s at 1 D to 3.7 s. G2 from 20 s to 4.7 s (primary) and from 21 s to 5.6 s (secondary). G4's primary reaches 1.1 s by 3 D; the secondary is still 18 s at the cap. G5 times the whole set only from about 0.14 (3 D), and the rms there is 5.2 s (primary) and 5.3 s (secondary), falling to 3.6 s and 3.9 s at 0.25. G7's primary is never timed: from 0.10 onward the reason is `walls_uncovered` on 51 of 53 eclipses. The secondary falls from 24 s at 2 D to 16 s at the cap.

At 1800 s, G1's primary is about 6 s from 6 D to the cap, and the secondary times 6 of 8 and stays near 9 to 11 s. G2 is about 6.4 s and 5.6 s at the cap. G4's primary is about 2.9 s from 3 D; the secondary is 48 s at 2 D and then stays near 35 to 38 s, and does not improve further. G5 times nothing at any width, because even 0.25 holds fewer than 10 points (`too_few_points` on every eclipse). G7 times nothing either (`walls_uncovered` on 51 eclipses of each kind). The 0.25 cap, the 10-point check and the contact-wall check are what stop the short-period stars at 1800 s.

On G1, G5, G6 and G7 at 120 s, and on the 600 s curves that do time, the widest window is the best of the grid or tied with it. On G2 at 120 s it does not beat about 4 D. On G4 at 1800 s the secondary does not improve as the window widens. The unmarked figures are `results/synthetic_scoring/window_scan_mik/{star}_c{cadence}/rms_vs_half_width.png`. The numbers are `scan.csv` in the same folders.

**Preliminary conclusions.** <u>All of them come from the synthetic set, one noise realisation each, no real light curve. The 120 s set is G1, G2, G4, G5, G6, G7 and G8. The 600 s and 1800 s sets are G1, G2, G4, G5 and G7.</u>

1. A window that is too narrow is costly, and the lower limit scales with D and not with a fixed phase. Below about 2 D the rms is several times its best value. From 3 D to 4 D it is within about 1.5 to 2 times, or better. In phase that is about 0.04 to 0.05 cycles on G1 and G2 but about 0.14 to 0.19 on G5. The MAVKA half-depth windows, converted above, run from 0.60 D to 2.72 D. On G1, G2 and G8 they sit on the steep part. A fixed ±0.12 is about 9 D on G1 and about 2.6 D on G5, where it costs roughly 30 to 40 percent over the best window.
2. The hypothesis that the widest window wins holds weakly. No star got worse as the window widened. G1, G2 and the G4 primary are flat beyond about 4 D, and what remains is within what 8 to 13 eclipses can resolve. G5, G7 and the wide secondaries (G4, G7) are still improving at the cap. The scan never shows a minimum, so it gives a lower bound on the width and no optimum. "Widest" means "up to 0.25 cycles", which is a limit of `time` and not of the data.
3. The formal sigma depends on the window. At 120 s it is about twice the scatter on G1, G4 and G7 and about right on G5. At 600 s and 1800 s, on G2 and the G4 secondary, rms(O-C)/sigma rises to 1.0 to 1.3 at wide windows, so the errors become optimistic. A width rule has to say which sigma is quoted.
4. At 1800 s the windowed method runs out of room on the short-period stars. G5 and G7 time nothing at any width. Part of this is the 10-point minimum: at 0.25 cycles G5 holds about 8 samples, so the gate fails by construction. G7 fails the contact-wall check, and it also fails at 600 s. Both are recipe gates, not necessarily the data. A rescan of these O-C curves with the gates relaxed was a ticket. It is postponed, and it may be cancelled: the timing experiment used the 120 s scans and did not ask for a new scan at 600 s or 1800 s.

What the scan cannot say:

- It measures precision, not accuracy. The rms is against the file ephemeris and shows how tightly the times cluster, not how close they sit to the truth.
- G6 at 120 s is the only spotted star in the scan. There the mean O-C moves by about 20 s as the window widens, so the scatter alone does not describe it. The nine-longitude spot experiment used ±0.12 and is a different test.
- The G4 and G7 secondaries stay at about 11 s to 40 s once the window is wide. That looks like a shallow, low-S/N eclipse and not a window effect, but it has not been shown.
- One realisation per file. Differences under about 10 percent are not resolved.

> **Proposed rule for the batch.** One rule only, not a pair of windows. The half-window of each eclipse is 4 D, with D1 for the primary and D2 for the secondary of the chosen solution of that file, and it does not exceed the largest window that does not overlap the neighbouring minimum (0.25 cycles for a circular orbit; G8 needs the actual separation of its two minima). Where 4 D is larger than that limit, the limit is used (the G4 secondary at 120 s and 600 s). The window is reported with every timing. <u>This rule is drawn from the synthetic set only. G6 at 120 s is in the scan, and there a wider window also moves the mean. The rule has not been used to time the batch, and it has not been tried on a real light curve.</u>

No width has been chosen for the batch yet. The rule above is a proposal for the researcher to confirm.

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg $PY examples/window_scan_mik.py --physics G1 --cadence 120
MPLBACKEND=Agg $PY examples/window_scan_mik.py --physics G1 --mark-only
```

**Recommendations before the real stars.** <u>They rest on the synthetic set only.</u>

1. Compare each method at its own sensible window and not at one common window. The expert (MAVKA) half-windows, 0.6 D to 2.7 D, are too narrow for Mikulášek. They sit on the steep part of its curve.
2. On the same fixed eclipse set, at 3 D, Mikulášek and WSAP agree to a few seconds. WSAP is a little tighter on the contact stars and on the G8 primary. The wide secondaries (G4, G7, G8) stay near 12 s to 15 s for both. The earlier batch standard deviations were a different quantity: every accepted eclipse, not this fixed set.
3. Most of the gain comes by 3 D to 4 D. Beyond that the rms falls by 10 to 30 percent, except on the S/N-limited secondaries.
4. The result may not transfer to real stars. G6 already shows a mean that moves by hundreds of seconds as the WSAP window widens, while the rms stays flat. Real spots, period changes and O'Connell variability can all act on the outer parts of a minimum. The window rule used for the batch is the section "Timing experiment: how the windows were chosen", and it does not use the mean.
5. Run the same scan on every real star. No truth is needed. Use the rms of O-C on a fixed eclipse set. The mean O-C is plotted against the window and reported. It does not choose the window, because on a real star it is set by the epoch of the ephemeris.
6. One eclipse set and one ephemeris for both methods, and each method reports its own window. Compare the methods by the paired difference per eclipse. The mean of the difference is the relative bias and its scatter is the disagreement. It does not say which method is right.
7. Expect windowed timing to fail on short-period stars at 1800 s. The low-cadence O-C rescan that would have tested the gates again is postponed, and it may be cancelled. Whole-cycle timing is the section "Whole-cycle timing (2026-11-07)".

**Whole-cycle timing (2026-11-07).** One shift of the stored template per orbital cycle, and no eclipse window. The template is not refitted. The shift moves the primary, the secondary and the harmonics by the same amount. The points of cycle n are those with phase in (n - 0.25, n + 0.75), so both minima sit inside the cycle. The per-cycle baseline is an offset and a slope. A cycle with fewer than 6 points, or fewer than 3 on either side of its centre, is not timed. Script `examples/batch_whole_cycle_time.py`. Tables `results/synthetic_scoring/batch_whole_cycle/summary.csv` and `eclipses.csv`. Figure `scatter_vs_pure.png`, against the relaxed run at a half-window of 0.25 cycles. Every file returned times. The rejected cycles are the ones at the ends of the series.

The two eclipses of one cycle are not timed separately. Their scatter about the pure minimum is the scatter of one shift, so a shallow secondary can look precise because it inherits the time of the deep primary. That is what happens on G4 and G7. The G7 primary matches the windowed run (0.6 s, 3.0 s and 7.9 s at 120, 600 and 1800 s). The G7 secondary then has that same scatter, where the windowed secondary was 11 s, 16 s and 44 s. The G4 secondary falls from 21 s to 40 s in the windowed run to 7 s to 10 s here, and the primary is worse (7 s to 10 s against 1 s to 3 s).

Where the two minima want different shifts, one number cannot be both. On G6 the windowed times sit near −47 s and +39 s from conjunction. The shared shift sits near −12 s for both, at every cadence, with a scatter of 2.3 s, 2.7 s and 4.8 s. The small scatter is the repeatability of that compromise. It is not a time of either minimum.

On the detached stars the window is the better measurement. G1 scatter is 12 s, 16 s and 19 s here, against 2 s to 9 s in the windowed run, and the shift is 8 s to 11 s early. G8 is 18 s to 21 s against 3 s to 8 s on the primary. G2 is close to the windowed run (5 s to 7 s). G5, where the eclipse fills the cycle, is slightly better here: 2.2 s, 2.3 s and 3.8 s against 2.8 s, 3.6 s and 5.4 s. The result is one noise realisation. It says the whole-cycle shift is a check on a shared ephemeris, and the eclipse window remains the time of one minimum.

**Higher baseline on the same whole-cycle shift (2026-11-07).** The degree-1 run is kept. Degree 2 and degree 3 are new directories, `results/synthetic_scoring/batch_whole_cycle_b2/` and `batch_whole_cycle_b3/`. The comparison figure is `batch_whole_cycle_baseline/scatter_vs_pure.png`. Nothing else was changed. A straight line over a cycle of 2 to 3 days was the suspected cause of the early shift on the long-period stars. Raising the degree tests that and nothing else.

It is the cause on the long-period stars, and it is not the whole of it. G1 primary scatter goes from 12.3 s to 6.1 s to 2.4 s at 120 s as the degree goes from 1 to 3, against 2.2 s for the eclipse window. The early shift goes from -7.8 s to -4.5 s to -2.0 s, against -0.5 s. At 1800 s the cubic still sits at 8.7 s and -5.1 s, against 6.2 s and -4.0 s for the window, so a cubic over 2.8 days has not caught the window. G8 moves the same way: the 120 s primary goes from 17.7 s and -7.4 s to 5.0 s and -2.9 s, against 3.2 s and -0.9 s. G2 at degree 3 is as tight as the window or tighter. G5 and G7 do not move, which is what a short cycle predicts, because a straight line over 0.3 to 0.5 days was already enough. G6 gets slightly worse (2.3 s to 5.2 s at 120 s) and the shared shift stays near -12 s. A polynomial cannot stand in for the spot, and it does not give the two minima different shifts. So the baseline explains the detached long-period stars. The shallow secondary still inherits the primary, and G6 still has one shift for two minima.

**Timing experiment: how the windows were chosen (2026-11-07).** Mikulášek has been timed. The MAVKA run at the windows below is still the request in the mail. The batch is timed twice for MAVKA, once at the expert half-depth windows already on disk and once at the windows below, and once for Mikulášek at its own windows below. The point of the second run is that the two methods get their windows from the same rule, not from the same number.

**Why not one shared window.** The expert half-depth windows sit at 0.6 D to 2.7 D. On the Mikulášek scan that range is still on the steep drop, so a comparison at those windows punishes Mikulášek for a window it was not built for. A single wide window does the opposite to AP, WSAP and WSL, which are built for the bottom of the minimum. Each method is therefore read off its own 120 s scan.

**Why the rms, and why not the smallest one.** The quantity is the rms of O-C about the mean, on the fixed eclipse set, against the file ephemeris. No truth column enters the choice, so the same rule can be used on a real star. The smallest rms is not taken. On Mikulášek it usually sits at the last point of the scan, and the change from the start of the flat region to that last point is smaller than the uncertainty of the rms. The rms from N eclipses is uncertain by `rms / sqrt(2(N-1))`, with N the number of eclipses in that rms. The adopted window is the narrowest one whose rms is no higher than the best rms by more than this uncertainty. No percentage was chosen by eye. Where N is small, about 8, the uncertainty is about a quarter of the rms and the window stays near 3 D. Where N is about 70, the uncertainty is a few tenths of a second and the window follows a real continued fall, which is why G5, G6 and the G7 primary land at 4 D to 8 D for Mikulášek.

**Why the mean O-C is not in the rule.** A mean far from zero is set by the epoch of the ephemeris. On a real star an honest window can sit tens or hundreds of seconds from the published epoch, and choosing the window that pulls the mean toward zero would be fitting the ephemeris rather than the eclipse. The researcher rejected that. What the mean does show, and what is kept as a result rather than a selector, is a mean that changes with the window. On G6 the WSAP primary mean runs from about -180 s to about -390 s across the region where the rms is flat. The rule still takes 2.2 D there. The timing reports that mean with the rms.

**Where on the existing scans.** No new scan was asked for. Mikulášek is searched out to its 0.25-cycle cap. MAVKA is searched over the curves as drawn: to 3 D, except G5 and G6, which have no continuum and are searched to 0.25 cycles, and the G7 primary, which is searched to 4.2 D1, the shoulder where the steep wall has become the slow outer rise. MAVKA's window is read from the WSAP curve. WSL falls within 0.1 D of it. AP, WSAP and WSL then share that one half-width, so the three methods are compared at the same window. AP's own scan is incomplete on G1, G2 and G8, which is a second reason not to give AP a window of its own.

**Why 120 s, and what happens at the other cadences.** The 120 s curve has the shape of the eclipse without the smearing of a long exposure, so the multiple of D is chosen there. D is the Mikulášek width of that eclipse at 120 s. Mikulášek carries the multiple to 600 s and 1800 s and multiplies by the D of the solution at that cadence, because a longer exposure changes D. MAVKA is not scanned again at 600 s and 1800 s. It keeps the 120 s half-width in cycles. Where that window holds fewer than 6 points, it is widened by one cadence in phase and counted again, until it holds 6. The widening stops if the two half-widths of the star would add to more than 0.5. The number of widenings is part of the result.

**What is compared, and what is left untouched.** The new MAVKA times are compared with the expert-window times, with Mikulášek at the windows in the table, and with the PHOEBE conjunction and the pure-curve minimum. The expert-window run stays in `output/batch_mix/`. The new run is requested in `agent_mailroom/LC_Modeller_to_Mavka_user_plateau_windows.md` and is to be written to `output/batch_mix_plateau/`. The window choice uses no truth. The score uses the truth. The scatter is still measured on the same noise realisation that chose the window, so that part is partly circular and the slide says so. These windows replace the 4 D proposal for this comparison. The 4 D proposal is not deleted: it was the earlier reading of the same scans, and this rule is what replaced it. Mikulášek was timed on 2026-11-07 by `examples/batch_plateau_time.py`. The tables are `results/synthetic_scoring/batch_plateau_mik/summary.csv` and `eclipses.csv`. No window had to be cut at 0.25 cycles. The outlier cut was off. The package gates were left as they are, so several 1800 s eclipses, and the G7 primary at 600 s, returned no time. The reason is `too_few_points` except the G7 primary at 600 s, which is `walls_uncovered`. On G1 B22 the file period is 2.80028 d, so the O-C against the file ephemeris walks by about 120 s, while the offset from conjunction stays under 2 s.

**Relaxed-gate run, maximal window (2026-11-07).** The aim was a Mikulášek time for every eclipse of every light curve. This is a second run, labelled as such, and the plateau run above is unchanged. Script `examples/batch_relaxed_time.py`, tables `results/synthetic_scoring/batch_relaxed_mik/summary.csv` and `eclipses.csv`. Three things differ from the plateau run. The half-window is 0.25 cycles for every eclipse, the widest `time` accepts without overlap. The minimum is 6 points per window (the MAVKA floor), with 3 on each side as in the package. The contact-wall check is off. It is not dropped silently: a second pass with the wall check on marks each time that would have failed it (`walls_ok` False, counted in `n_walls_failed`). The outlier cut is off. The solutions are the stored ones. Nothing is refitted or filled in. An eclipse that still cannot be timed keeps a row with its reason. The rule is no longer the plateau rule, so these windows are not to be read as a recommendation for the method. They test whether the strict gates are what removes the low-cadence times.

Result. Every one of the 44 kind-and-file rows now has times. The wall check would have failed in 10 rows. It fails for every timed eclipse of G5, G6 and G7 at 1800 s (6 rows), G7 primary at 600 s, and G8 primary at 1800 s, and for 2 G1 secondaries and 1 G2 secondary at 1800 s. The remaining rejections are `one_side_uncovered` (1 to 2 per row, the first or last eclipse at the edge of the data) and a few `too_few_points`. Where the plateau run had times, the relaxed times agree with them within the scatter, and the rms is mostly lower or equal because the window is wider. The mean offset from conjunction stays within a few seconds: G5 at 1800 s +0.3 s primary and +0.1 s secondary, G7 primary at 600 s and 1800 s +0.1 s and +0.6 s, G6 primary and secondary about -50 s and +40 s at every cadence as at 120 s. The scatter grows with cadence, for example G5 primary 2.8, 3.6 and 5.4 s, and G7 secondary 11, 16 and 44 s. G4 secondary and G8 secondary stay the poorest: G4 secondary scatter is 21 to 40 s and sits about -9 s from conjunction at 1800 s, and G8 secondary is 14 s from conjunction at 1800 s with a 22 s scatter. On this noise realisation the wall check therefore did not protect against a bias, since the walled-off times of the contact stars are as close to conjunction as the others. One realisation does not show that the check is unnecessary. The formal sigma is not trusted at 6 to 8 points, and the scatter of the times is the number to quote.

**Comparison of all methods on the batch (2026-11-07).** Mavka_user delivered AP, WSAP and WSL at the plateau windows (`agent_mailroom/Mavka_user_to_LC_Modeller_plateau_windows.md`, tables in `Mavka_user_folder/output/batch_mix_plateau/`, expert-window tables untouched). Script `examples/batch_compare_methods.py` reads the tables and times nothing. Every accepted time is matched to the nearest pure-curve minimum of the same kind and compared with it and with the PHOEBE conjunction. Results in `results/synthetic_scoring/batch_compare/`: `summary.csv`, `scatter_vs_pure.png`, `bias_vs_conj.png`, `bias_vs_pure.png`, and the same scored on one common set of eclipses (`summary_common.csv`, `scatter_vs_pure_common.png`, `bias_vs_conj_common.png`). The second G1 600 s file (B22, wrong period) is left out. Colours as before, Mikulášek red (solid: plateau windows, dashed: relaxed gates with 0.25 cycle windows), MAVKA at plateau windows solid and at the expert windows dotted.

Why the common set. The methods did not time the same eclipses, and one eclipse can dominate a scatter. On the G1 primary the plain WSAP scatter is 124 s at 120 s, 224 s at 600 s and 2700 s at 1800 s. This is one eclipse, E = 5, which has 38 of the 76 points of the window (half of it is in a gap). WSAP accepted it with a formal sigma of 163 s and a time about 400 s off. AP and WSL rejected it, and Mikulášek rejected it as `one_side_uncovered`. The scatter in `scatter_vs_pure_common.png` is on the eclipses that Mikulášek relaxed, AP, WSAP and WSL all timed. The conclusions below use that figure. Note that WSAP's own formal sigma flagged this eclipse, so the failure is detectable, but it is a time the method returned.

What the numbers say (scatter about the pure minimum, seconds, one noise realisation, common set).

- 120 s. All methods are within a factor of about two and mostly within 1 s to 3 s of each other. Primary: Mikulášek plateau 0.7 to 3.3 s against WSAP 0.8 to 3.1 s and WSL 1.3 to 3.1 s (G1 2.6 against 1.7 and 1.8, G2 1.8 against 1.0 and 1.3, G4 2.9 against 1.6 and 2.2, G5 3.2 against 2.9 and 2.9, G7 0.7 against 0.8 and 1.8, G8 3.3 against 3.1 and 3.1). The relaxed Mikulášek run is lower than the plateau run on the primaries (G4 1.8 s). AP is the weakest on most stars (6.5 s on G1, 6.3 s on G2, 6.5 s on G7). The G4 secondary is poor for all (15 to 32 s) and Mikulášek plateau is the lowest at 14.8 s.
- 600 s. On the contact star G6 Mikulášek holds 4.6 s (relaxed 3.9 s) while AP, WSAP and WSL go to 7.1, 10.2 and 9.9 s. On G5 Mikulášek and WSAP and WSL are together at 3.7 to 3.9 s. On G2 and G4 the MAVKA methods are better than Mikulášek plateau (G2 primary WSAP 2.3 s against 5.8 s, G4 primary 2.4 s against 3.8 s), and the relaxed Mikulášek run recovers some of it (G4 1.1 s). The Mikulášek plateau G2 secondary is 14 s against 5.6 s for the relaxed run, so the 1.5 D window chosen at 120 s is too narrow for the 600 s file.
- 1800 s. The strict Mikulášek gates time only the G1 and G2 primaries and 2 eclipses of the G4 secondary. The relaxed run times everything. There it is the best or close to the best on G6 (6.4 s against 29.8, 34.1 and 137 s), G7 primary (9.5 s against 178, 61 and 151 s), G2 and G4 primaries (6.6 s and 2.9 s against 12.7 s and 3.9 s for WSAP). WSL is as good on G1 (5.8 s), G5 (5.0 s against 5.4 s) and G8 (4.8 s against 8.4 s). AP is the worst at 1800 s on almost every star, with 14 to 178 s. The MAVKA windows at 1800 s are the 120 s half-widths widened by up to three cadence steps, so they are narrow, while the relaxed Mikulášek window is 0.25 cycles. This part of the comparison is partly a comparison of windows, not only of methods.
- Offsets from conjunction. At 120 s and 600 s all methods are within a few seconds of conjunction on every star except G6 and the secondary of G8. G6 (O'Connell spot): Mikulášek -45 s primary and +39 s secondary at every cadence, AP -93 and +123 s, WSAP and WSL -183 and +184 s at 120 s and -350 and +397 s for WSAP at 1800 s. The pure-curve minimum is +125 and -150 s from Mikulášek, so no method sits on both references. The Mikulášek times are the closest to conjunction. The G8 secondary sits at +5 to +7 s at 600 s and +14 to +20 s at 1800 s for every method including Mikulášek, so that shift is shared and is not a property of one method.
- MAVKA plateau windows against the expert windows. They are close. The plateau windows are better on some stars at 120 s (G1 primary WSAP 1.7 s against 3.5 s) and worse on others (G2 primary at 600 s, WSAP 2.3 s against 1.0 s). The rms-plateau rule did not give a general improvement for MAVKA on the common set.

Conclusions the data support. At 120 s, Mikulášek, WSAP and WSL give the same precision, and AP is somewhat worse. The difference between Mikulášek and the other methods grows with the exposure and the width of the minimum. On the contact stars and the shallow-secondary star at 600 s and 1800 s the template over a wide window keeps the scatter at 4 to 9 s where the narrow-window methods reach 30 s to 180 s. This needs the relaxed gates, which are not the package defaults, and the wall check would have rejected those files. The offset from conjunction on G6 shows that agreement with the pure minimum, with conjunction and with the other methods are three different questions, and none of them says which time is physically right. What is not shown: whether the result holds with another noise realisation, with a spotted real star, or with MAVKA windows as wide as the template; the windows were chosen on this same realisation, so the scatter is partly circular; the 1800 s MAVKA windows were not chosen for that cadence.

**Borrowing the 120 s template for the 1800 s files (2026-11-07).** A test that the scripts work, on synthetic files, before the real-star use (a 120 s sector as the base for a 1800 s one). Script `examples/borrow_batch_1800.py`, results in `results/synthetic_scoring/borrow_batch_1800/` (`fits.csv`, `timing.csv`, `timing_common.csv`, `residuals_vs_phase.png`, `timing_scatter.png`, solutions as JSON). Stars G1, G2, G4, G5, G7, G8, one 1800 s file each. G6 is left out by decision, since borrowing from a spotted template is not what this test is for. The reference is the stored chosen 120 s solution of the same star.

What is frozen. The profile, Gamma, C, K, the harmonics and the period are frozen in every borrowed fit, because they describe the shape of the star's eclipse and should not follow the cadence. The free groups are the choice. A new group `widths` (D alone) was added to `ParameterPlan` so that D can be freed while Gamma stays fixed. The name `borrow_phase` means one thing only: the epoch of the primary minimum is the template parameter that is fitted. It does not mean a phase window. The baseline of the 1800 s file (its own polynomial, split at gaps) is fitted in every variant, borrowed or not, because that baseline belongs to the light curve and not to the star. Frozen in every borrowed fit: the profile number, Gamma, C, K, the harmonics and the period. The four variants then differ by what else is fitted. `borrow_phase`: epoch only, so the depths and the widths stay at the 120 s values. `borrow_amp`: epoch and the eclipse depths. `borrow_width`: epoch and the eclipse widths D. `borrow_amp_width`: epoch, depths and widths. The profile figure and the trend figure of every solution are in `fits/` and `trends/` next to the tables. Amplitudes and widths are left as options because depth and width can change for physical reasons, not only through exposure smearing. The 1800 s smearing is calculated: the borrowed model is integrated over the exposure written in that file, 1440 s (0.8 of the cadence, 21 sub-samples). The 120 s solutions were fitted with the exposure integration off, because 96 s is under a tenth of the eclipse width, so that small smear stays inside the borrowed D and is not removed before the 1440 s integral. The baseline and gap settings are those of that file's own stored solution. The timing is the relaxed recipe (0.25 cycles, 6 points, wall check off), on the whole file, scored against the pure minimum and the conjunction, on the eclipses timed by every solution.

Correction to the plan. I had said the stored 1800 s from-scratch fits were made without the exposure. They were not: the stored solutions already carry `exposure_days` of 1440 s with 21 sub-samples, because the candidate generator switches integration on. The extra from-scratch refit with the exposure (`scratch_exp`) therefore reproduces the stored solution to every printed digit and adds nothing. It stays in the tables as a check that the refit is repeatable. The real comparison is the stored fit against the four borrowed fits.

Result (rms in mmag, scatter of the times about the pure minimum in seconds).

| Star | from scratch rms | borrowed rms, phase / amp / width / amp+width | scratch scatter pri / sec | borrowed (amp+width) pri / sec |
| --- | --- | --- | --- | --- |
| G1 | 2.18 | 2.17 / 2.16 / 2.16 / 2.16 | 6.2 / 8.6 | 6.0 / 8.4 |
| G2 | 2.16 | 2.26 / 2.22 / 2.24 / 2.22 | 6.6 / 5.4 | 4.8 / 4.5 |
| G4 | 2.79 | 4.66 / 4.61 / 4.60 / 4.60 | 2.9 / 39.8 | 3.0 / 40.2 |
| G5 | 1.91 | 2.10 / 2.05 / 1.97 / 1.96 | 5.4 / 7.2 | 5.6 / 7.5 |
| G7 | 3.55 | 10.8 / 3.66 / 3.74 / 3.24 | 8.0 / 44.2 | 7.8 / 48.1 |
| G8 | 2.17 | 2.13 / 2.12 / 2.12 / 2.12 | 8.4 / 22.1 | 9.3 / 22.3 |

- The borrowed template follows the 1800 s data as well as the 1800 s fit from scratch on G1, G2, G5 and G8, with the rms within about 3 percent and, on G1, G8 and G2 (timing), at or below the from-scratch values. Freeing the amplitudes and widths gains little there.
- G7 needs freedom. With the phase only, the rms is 10.8 mmag and the primary scatter 31 s. Freeing the amplitudes brings the rms to 3.7 and the scatter to 12 s, and amplitudes with widths give 3.2 mmag and 7.8 s, equal to the from-scratch fit. The 1800 s from-scratch fit has a primary amplitude of 1.20 against 0.47 at 120 s, so the depth is not carried over by the 120 s template.
- G4 is worse in rms (4.6 against 2.8 mmag). The residual plot shows spikes of up to 10 mmag at the eclipses. The 1800 s from-scratch fit chose a different profile (profile 1) from the 120 s reference (profile 4), and the frozen profile 4 with its Gamma does not reproduce the 1800 s eclipse. The timing is unaffected (2.9 s against 3.0 s on the primary), so the rms on the eclipse wall does not move the time of the minimum.
- The timing biases against the conjunction are the same to within a second for every solution. The G4 and G7 secondaries and the G8 secondary stay as poor as before, for every method.
- The result is on one noise realisation, six stars, one cadence, and a template with the same noise-free shape. It shows that the code runs and that a frozen shape is enough on four stars of six, which is what the test was for. It does not show how well the method works on a real star, where the 120 s and the 1800 s sectors differ by more than exposure and noise. That is the next step.

| Star | Mikulášek, pri / sec | MAVKA, pri / sec |
| --- | --- | --- |
| G1 | 2.8 D, 0.0373 / 2.6 D, 0.0347 | 1.4 D, 0.0187 / 1.5 D, 0.0200 |
| G2 | 2.9 D, 0.0540 / 1.5 D, 0.0281 | 0.8 D, 0.0149 / 0.6 D, 0.0112 |
| G4 | 1.4 D, 0.0579 / 1.5 D, 0.1130 | 0.7 D, 0.0289 / 0.7 D, 0.0527 |
| G5 | 4.7 D, 0.2185 / 4.7 D, 0.2238 | 1.7 D, 0.0790 / 1.5 D, 0.0714 |
| G6 | 4.3 D, 0.1983 / 4.7 D, 0.2214 | 2.2 D, 0.1014 / 1.9 D, 0.0895 |
| G7 | 7.8 D, 0.1147 / 6.8 D, 0.2280 | 2.4 D, 0.0353 / 1.4 D, 0.0469 |
| G8 | 2.6 D, 0.0248 / 1.7 D, 0.0294 | 1.3 D, 0.0124 / 1.0 D, 0.0173 |

The second number of each pair is the half-width in cycles at 120 s.

**Which y-axis: rms or rms/sigma.** The rms (in seconds) is the precision. It is the quantity to optimise. The ratio rms/sigma compares the scatter with the formal error bar: about 1 is honest, below 1 is conservative, above 1 is optimistic. It is a calibration and not a quality. Optimising it would pick the narrowest windows, where the rms is worst, because the formal sigma is largest there. Use it as a constraint: the smallest rms such that rms/sigma stays near or below 1. With 8 to 13 eclipses the ratio itself has an uncertainty of about 25 to 30 percent. On NSV 807 the scatter was smallest at 0.15 and the errors were honest at 0.25. On real stars the rms also contains real O-C variation.

**MAVKA on the same eclipses (2026-11-07).** AP, WSAP and WSL were run on the same 120 s files, the same fixed eclipse set and the same file ephemeris as the Mikulášek scan. The rms is about the mean of that set. No outlier cut. The step is D/10. G4, G5, G6 and G7 were run out to 0.25 cycles, for all three methods. G1, G2 and G8 have AP only as far as the stopped job had reached, and WSAP and WSL out to 3 D. G8 has no AP secondary. The reply is `agent_mailroom/Mavka_user_to_LC_Modeller_window_scan_wsap_wsl.md`. Tables: `/home/voz/projects/UPJS/Mavka_user_folder/output/window_scan_mavka_120/`.

**Why the MAVKA curves stop at 3 D.** AP, WSAP and WSL are built to fit the bottom of the minimum, not the continuum outside it. D is the scale of the Mikulášek template of that eclipse. On a normal eclipse the template is still a large fraction of its central depth at one D, and by about two D it is down to a few percent, which is the wing. Three D is past the wing. G4 is fitted with profile 4, which is exactly zero once the phase from the centre exceeds D, so there its D is the edge of the model eclipse. The G7 primary is still deep at one D. The rms of WSAP and WSL has already stopped falling there, so a wider window adds continuum these methods are not meant to fit. On the long-period stars that continuum was also the slow part of the AP run. The default on the figure is therefore 3 D. This is a limit of the drawing, not of the tables: G4, G5, G6 and G7 still hold the grid out to 0.25. Mikulášek is drawn in full, in red, because its template is the one being tested past the well. Three eclipses are drawn past 3 D, for the reasons in the next paragraph.

rms of O-C at 3 D, in seconds. One noise realisation. G8 secondary has no AP time.

| Star | Kind | Mikulášek | AP | WSAP | WSL |
| --- | --- | --- | --- | --- | --- |
| G1 | primary | 2.7 | 4.8 | 2.9 | 2.8 |
| G1 | secondary | 2.6 | 4.5 | 2.4 | 2.1 |
| G2 | primary | 1.8 | 3.3 | 2.1 | 2.1 |
| G2 | secondary | 2.3 | 3.9 | 2.4 | 2.2 |
| G4 | primary | 1.3 | 3.0 | 1.9 | 1.6 |
| G4 | secondary | 13.5 | 16.5 | 13.6 | 14.8 |
| G5 | primary | 3.9 | 3.9 | 3.0 | 3.1 |
| G5 | secondary | 4.9 | 5.7 | 3.5 | 3.5 |
| G6 | primary | 4.1 | 4.5 | 3.6 | 4.0 |
| G6 | secondary | 5.0 | 5.8 | 4.1 | 5.0 |
| G7 | primary | 1.7 | 4.5 | 1.6 | 0.8 |
| G7 | secondary | 16.2 | 21.3 | 12.9 | 13.0 |
| G8 | primary | 3.2 | 5.0 | 1.6 | 2.1 |
| G8 | secondary | 15.2 |  | 12.0 | 12.9 |

WSAP and WSL have reached this floor by about 1 D to 2 D, except the G7 primary, which is still unusable at 1 D (WSAP 108 s) and has fallen by its adopted window at 2.7 D. Widening MAVKA past 3 D does not tighten the rms. On G6 the rms is already flat while the mean keeps moving: at 3 D the WSAP mean is -244 s (primary) and +281 s (secondary). That is the spotted star, and the rms alone does not describe it. Mikulášek at the adopted MAVKA window is still on the steep part of its own curve (G1 about 20 s, G8 secondary 107 s) and only meets these values near 3 D.

**Where the MAVKA curves are drawn past 3 D.** No further run was asked for. The points were already in the tables, and only the figure was cutting them off.

G5 and G6 are contact stars. The minimum and the variation between the maxima fill the orbit, so there is no flat continuum at which a window of 3 D can be said to have left the eclipse. Three D is about 0.14 cycles on these stars. The saved scans already reach 0.25 cycles, the same cap as the Mikulášek scan, and the figure now draws AP, WSAP and WSL out to that cap. The rms does not keep falling there. At 0.25 cycles, in seconds:

| Star | Kind | Mikulášek | AP | WSAP | WSL |
| --- | --- | --- | --- | --- | --- |
| G5 | primary | 2.9 | 6.6 | 3.2 | 3.3 |
| G5 | secondary | 3.3 | 7.1 | 4.1 | 4.2 |
| G6 | primary | 2.9 | 10.7 | 3.6 | 4.8 |
| G6 | secondary | 3.3 | 6.7 | 5.6 | 7.0 |

AP is noisier at the cap than at 3 D. WSAP and WSL stay near their 3 D values, or a little worse. On G6 the mean is still moving at 3 D (WSAP -244 s on the primary and +281 s on the secondary), so a flat rms along this longer stretch is not a settled time.

The G7 primary is a different case. Its solution is profile 1 twice: the first component has D1 = 0.0147 cycles and Gamma 7.30, the second has D = 0.0230 cycles (1.56 D1) and Gamma 5.37. Profile 1 never reaches exactly zero, so it has no phase at which the eclipse ends. The inflection of each component is the middle of its wall, at 1.88 D1 and 2.79 D1, and about half that component's depth is still left there. The sum of the two has its inflection at 1.98 D1, with 72 percent of the combined depth still left. That is the steep part. The wing after it is short, because both Gammas are large: 23 percent of the combined depth is left at 3 D1, 1.3 percent at 4 D1, and none by 5 D1. On the pure curve the steep slope has given way to the slow outer rise at 4.2 D1. The figure draws the G7 primary to that shoulder. The saved scan continues to 17 D1 and is not all drawn. At 4.2 D1 the rms is 1.1 s (Mikulášek), 1.3 s (AP), 4.2 s (WSAP) and 0.7 s (WSL). WSAP is worse there than at 3 D, where it was 1.6 s. The G7 secondary stays at 3 D. Its Gamma is 1.00, and its wing is already down to a few percent by about 2 D.

**Figures.** Script `examples/window_scan_compare.py`. It reads the saved scan tables and does not fit anything. One figure per star, `{star}_rms_and_mean.png`: primary and secondary, each with the rms of O-C and the mean O-C against half-window over D. Mikulášek is the red line, drawn for the whole grid. AP is blue, WSAP orange, WSL green. They stop at 3 D, except where the previous paragraph draws them further: G5 and G6 to 0.25 cycles, and the G7 primary to 4.2 D1. The dashed line is the adopted MAVKA half-width, the dotted line is 3 D, and the dash-dot line is the Mikulášek cap at 0.25 cycles. The overview is `overview_rms.png`: the same rms curves (thick, left axis, logarithmic), the mean O-C of the same methods in the same colours (thin, right axis, in seconds, symmetric logarithmic: linear inside ±10 s), and a third column with the noise-free one-orbit light curve of that star over phase -1 to 1. The period is written in the bottom-right corner. Each light curve has three horizontal chords, at the flux of the pure curve at the two ends of the chord. Red, centred on phase 0: length 2 D1, from -D1 to +D1. Green, centred on phase 0.5: length 2 D2. Blue, centred on phase 0: length 2 x 3 D1, the full width of the MAVKA half-window of 3 primary D. The legend on the first light curve says "2 D1", "2 D2" and "2 x 3 D1", and it applies to every row. What D is, and what the timing fits, is in "How a Mikulášek time is measured". Files: `results/synthetic_scoring/window_scan_compare/`.

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg $PY examples/window_scan_compare.py
```

What this comparison cannot say: which time is closer to the true minimum. The rms is the scatter about the mean, not the offset from conjunction. One realisation. The formal sigma of the two methods is not defined in the same way, so rms/sigma is not used to rank them.

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

### Sector 4 error column, checked on the MAST file (2026-11-08)

The redownloaded sector 4 VOTable reports a SAP error of 0.47 mmag. That is the MAST column, not a conversion mistake. In `hlsp_tess-spoc_tess_phot_0000000035119266-s0004_tess_v1_lc.fits` the column `SAP_FLUX_ERR` has median 2.45 e-/s on a median `SAP_FLUX` of 5720 e-/s, which is 0.47 mmag, and `PDCSAP_FLUX_ERR` has median 2.98 e-/s. The times in that file are BJD minus 2457000. The keyword `EXPOSURE` is 20.53 days, the sector livetime. The cadence is `TIMEDEL` = 1800 s, and `INT_TIME` is 1.98 s. The 175 mmag error and the constant residual sigma of 3.64 mmag belong to the earlier QLP file. No constant sigma is applied to this sector 4 file. Sector 4 is not fitted in this step.

### Tiered search on NSV 807, sector 30 SAP (2026-11-08)

Both sector 30 files were fitted from the SAP flux only. `pdcsap_flux` was not read. Magnitudes use the photcal of `sap_flux`. The exposure is 0.8 of the cadence, 96 s and 480 s, with 21 subsamples. It is taken from the times, not from the file `exptime`. At 120 s the package gate would have left the integration off, because 96 s is about 0.02 of the measured FWHM and the gate asks for 0.1. It was switched on for this star. At 600 s the gate itself asked for it (0.12 of the FWHM).

The search is the three-step one (`examples/nsv807_s30_tiered.py`, results in `results/nsv807/s30_tiered/`). Trust-region fallback off. The features call the star contact-like: equal depths, no flat stretch, so the harmonic sets are (2, 3) and (1, 2, 3). The O'Connell term is on. Each step compares the previous winner with the new candidates. A step that gains less than 10 in BIC does not, by itself, stop the search. The search stops only after two such steps in a row. That stop did not fire. Every finished candidate keeps its profile figure and its trend figure (12 of each). Nothing is timed.

| File | Step 1, baseline | Step 2, profile | Step 3, harmonics | Chosen | rms |
| --- | --- | --- | --- | --- | --- |
| 120 s SPOC | constant per segment beats a line by 8.7; gain 0 | profile 3, gain 28.1; profiles 1 and 2 stay at 5.08 mmag; profile 4 is 6.74 | H1 does not pay: rms 4.41 against 4.48, BIC worse by 10.5 | profile 3, harmonics (2, 3), constant per segment | 4.50 mmag |
| 600 s TESS-SPOC | same, a line is worse by 8.9 | profile 3, gain 72.9 | H1 wins by 4.6, under the decisive difference of 10 | profile 3, harmonics (1, 2, 3), constant per segment | 5.90 mmag |

The two cadences agree on profile 3 and on a constant baseline in each of the two segments. They disagree on the first harmonic, and the 600 s gain that brings it in is smaller than the difference the search treats as decisive. Both chosen solutions have status ok, and both report a rank-deficient Jacobian, so the formal parameter errors are not used. The chosen profile figure and trend figure are copied to `fits/chosen/` and `trends/chosen/`. The eclipse set for the later timing is the one fixed by coverage at the widest window, the same set for both cadences.

The primary is flat at the bottom and none of these single-component templates follow it. The binned primary residual of the 120 s choice swings from +4.9 mmag to -6.9 mmag. The search does not propose a two-component candidate. Inside a fit the second component is tried only when rms/noise exceeds 1.5. Profile 3 is at 1.35, so it was skipped. Profiles 1 and 4 did try it and rejected it: the rms gain was 13 percent, below the gate of 25 percent, and only one start agreed. Forcing the second component on the chosen profile 3, primary only, lowers the 120 s rms from 4.50 mmag to 4.16 mmag and the 600 s rms from 5.90 mmag to 5.79 mmag. The starts still do not agree, so the package does not adopt it. The drawn trial, with that one gate loosened so the curve is visible, is `fits/TIC35119266_s30_SPOC_120s__2c_primary_shown.png` and the same name for the 600 s file. The primary bins are still about 4 mmag. The second component flattens the bottom a little. It does not follow the flat bottom.

### Protocol of the tiered fit for a real light curve (locked, 2026-11-08)

Status: locked on 2026-11-08. The two files named `..._2c_primary_shown` in that section are exploratory trials, not results: a condition of the package was loosened to draw them, and their titles say "adopted" although the package rejected them. They are to be removed from the results folder.

Code: `examples/nsv807_s30_tiered.py` runs the protocol, steps 0 to 7. By default it reads the NSV 807 sector 30 file and writes `results/nsv807/s30_protocol/`. Another light curve is run with `--file`, `--out`, `--flux` and `--flux-err`, for example for AA And: `python nsv807_s30_tiered.py --file ../data/lc_discovery_TIC_455790537_SPOC_sector_16.vot --out ../results/tic455790537_s16_protocol --flux flux --flux-err flux_err`. It uses `select_model` and the candidate builder of the package (`src/mikulasek_eb/selection.py`, `candidates.py`), and the option `keep_curvature` added to the fit policy (`config.py`, `pipeline.py`). `examples/nsv807_s30_forced.py` is not part of the protocol. It forces the second component on a saved candidate, as a look.

Marks: **Decided** is a choice you made. **Package** is how the package works now. **Suggestion** is not decided.

**Step 0. Input.** Decided: SAP flux only, converted with the photcal of the file. Times are TDB, barycentric. Period and epoch come from the file. The HJD against TDB question of the epoch is ignored. Decided: the exposure is 0.8 of the cadence (96 s at 120 s, 480 s at 600 s), because the `exptime` of the file is the cadence. Package: exposure integration is on only if the exposure divided by the eclipse width exceeds 0.1. Decided: we follow that rule and do not force it on.

**Step 1. Outlier control.** Two options were tried on sector 30, without refitting. Decided: develop option A first. Option B stays as a comparison.

*Option A. 5 sigma from the folded light curve.* Fold by the file period and epoch. Take the median magnitude in 400 phase bins. The residual of each point is its magnitude minus the median of its bin. The local sigma is 1.4826 times the median absolute deviation of the residuals in that bin. A point is flagged when its residual is more than 5 local sigma. It needs only time and magnitude, so it works on any file. Result on sector 30:

| | 120 s | 600 s |
| --- | --- | --- |
| Points flagged | 106 of 18354 | 50 of 3666 |
| Within 0.08 cycles of an eclipse | 17 | 5 |
| Brighter / fainter than the fold | 98 / 8 | 24 / 26 |
| Scatter of the whole curve | 4.45 mmag | 6.30 mmag |

The flagged points lie at all phases, none above the maximum or below the minimum. They fall on days 11 to 13 and 26 to 27 after the start, where `sap_bkg` rises from a median of about 1800 e-/s to 70 000 (120 s) and 340 000 (600 s) e-/s. Some flagged points sit almost on the curve, on the steep branches, where the local sigma of a bin is small. Nothing was tuned. The threshold of 5 and the 400 bins were first choices. At 600 s 400 bins hold about 9 points each, so the local sigma is poorly measured there. Open: whether flagged points are removed once or in rounds, and the number of bins. Pictures, in `results/nsv807/s30_outliers/`: `TIC35119266_s30_SPOC_120s_outliers.png` and `TIC35119266_s30_TESS-SPOC_600s_outliers.png` (residual against phase, against time, and `sap_bkg`), and `..._folded_outliers.png` for both files (flagged points on the folded curve).

*Option B. Steep background.* A point is removed when the rate of change of `sap_bkg` is above the median of that rate plus 5 robust sigma. The threshold comes from the data: 9850 e-/s per day at 120 s and 2341 at 600 s. It needs the `sap_bkg` column, which another file may not have. Result:

| | 120 s | 600 s |
| --- | --- | --- |
| Points removed | 1009 (5.5%) | 371 (10.1%) |
| Scatter of kept points | 3.77 mmag | 3.35 mmag |
| Robust scatter (MAD sigma) | 3.43 mmag | 2.84 mmag |
| Option A flags left among kept points | 7 (was 106) | 14 (was 50) |

It removes the whole rising branch of each orbit end, so good points are removed with the bad ones, evenly at all phases. A rule on the background level (median plus 5 robust sigma) removes 19.5% of the points and was not used. Pictures: `..._bkg_rule.png` for both files, the folded curve and `sap_bkg` against time with the removed points in red.

Rule for either option: removed points are listed and are not deleted from the file. The rms is reported for the kept points, and the number removed is reported next to it.

**Step 2. Candidates.** Package: from the folded curve the package measures features and proposes the harmonic sets, the baseline orders, the O'Connell term and whether to integrate over the exposure. For sector 30 the star is called contact-like (equal depths, no flat stretch), so the harmonic sets are (2, 3) and (1, 2, 3). The baseline orders are 0 and 1. The O'Connell term is on. Package: the baseline is one level or one line on each piece of the light curve, and a new piece starts wherever two points are more than 0.5 d apart. Sector 30 has one such gap, of about 1.07 d. Decided: the baseline is split there, one polynomial per orbit, with one shared eclipse template. This is what the package does and what the base paper does, with one instrumental term for each observational subset (Sect. 3.1). The 0.5 d value is the package's, not the paper's. Decided: the search starts from the richest proposed harmonic set, (1, 2, 3), and step 3 tries the smaller ones. Why: a set that is too small leaves a large rms/noise and sends every candidate into the second-component trial.

**Step 3. Tiered search.** Each step compares the previous winner with new candidates on the same points.

1. Baseline. Profile 1, the smaller harmonic set, baseline order 0 against 1. Why: the baseline is cheap to decide and does not change the shape.
2. Profile. The winning baseline and harmonics, profiles 1 to 4 of Table 2 of Mikulášek (2015), and profile 1 again with the curvature C1 and C2 kept free and never set to zero (`keep_curvature`, name suffix `_C`). Why: the profile is the main choice of shape, and the staged fit otherwise sets C to zero for flat-bottomed minima. Profiles 5 to 7 are not tried. Package change: `keep_curvature` is a new option of the fit policy.
3. Harmonics. The winning profile and baseline, every other proposed harmonic set.

The selection criterion is BIC, and a difference of 10 counts as decisive. The search stops when two steps in a row each gain less than 10. A candidate counts only if its fit status is `ok`. The trust-region retry is off: it was copied from the synthetic tiered experiment and is our choice. Known weaknesses: the baseline is decided with profile 1, which is wrong for this star. BIC treats every point as independent, so its gains are overstated. BIC judges the whole curve, and timing needs the eclipse shape.

**Step 4. Second profile component.** Package: it is not a candidate of the search. Inside the fit of each candidate it is tried only when the root mean square residual is more than 1.5 times the noise measured outside the eclipses. It is kept only if the rms falls by at least 25 percent and at least two of the four starting values reach the same fit. Suggestion: the package has an option that adds a twin of every candidate with the second component forced on. The twins would be judged by the same BIC, with no rms condition and no later look.

**Step 5. Baseline.** One polynomial per orbit, as in step 2, with one shared eclipse template. The order is chosen in step 3 between the orders the package proposes. The baseline is not cut at scattered-light rises, only at gaps over 0.5 d.

**Step 6. Pay more attention to the rms in the minima.** Suggestion, not decided. Timing uses the eclipse shape, and the root mean square residual of the whole curve can look good while the flat bottom is missed. Report the rms in a window around each minimum, and the binned residual at the minima, next to the whole-curve rms. It would be a reported number. It would not replace BIC.

**Step 7. Output.** One folder `chosen` with the final solution and its two figures per light curve, nothing else. The title of every figure has the light curve name, a line break, then the profile number with the harmonics and baseline, the fit status (`ok` or `bound-active`), `adopted` or `not adopted` for the second component, and the rms in mmag. Every candidate keeps its figures in a separate folder. Figures are redrawn from the saved solutions. A run is not repeated without being asked.

### Results of the locked protocol on NSV 807 and AA And (2026-11-08)

*NSV 807, sector 30, 120 s (18354 rows, 106 removed by step 1; the whole run takes 43 s).*

| Step | Winner | rms (mmag) | BIC gain |
|---|---|---|---|
| 1 baseline | p1_h123_b0 | 4.83 | 6.9 over order 1 |
| 2 profile | p3_h123_b0 | 4.28 | 1.4 over profile 1 |

Profile 2 is 4.84 mmag, profile 4 is 6.51 mmag and the curvature candidate `p1_h123_b0_C` is 4.72 mmag. The search stops after step 2. The earlier run with the smallest harmonic set gave the same winner, p3 with rms 4.35 mmag. C1 and C2 of that `_C` candidate were -0.078 and -0.191, so the curvature made the bottom more pointed, not flatter.

*AA And (TIC 455790537, sector 16, 120 s, 16717 rows).* The whole run takes about 45 s. Baseline order 1 wins (rms 2.81 mmag, BIC better by 9.0), profile 1 and profile 2 tie, and the package adopts a second component by its own rules. The adopted fit follows the flat bottom of the primary and the residuals are about 5 mmag with no spike at the minima. This was read from the figure. The trial record was not opened.

**Conclusion: the flat-bottomed, asymmetric minima of NSV 807 are not followed.** Under the locked protocol no candidate of the search reproduces the minima. The fit is pointed where the data are flat, and the residuals have spikes of 15 to 35 mmag at phases 0 and 0.5 with an antisymmetric wiggle beside them. What was tried and what it gave:

- Curvature C free (`_C`): rms 4.72 mmag, BIC 6.8 worse than the winner. C went negative.
- Second component forced on the chosen profile 3: the rms gain is 6 percent on the primary only and 10 percent on both minima. The package rejects it because only one of four starts agrees.
- Second component forced on profile 1 with harmonics (1, 2, 3): gain 12 percent (primary) and 19 percent (both). The same single failed condition, one start of four.
- Weight 3 inside 0.12 cycles of the minima: the trial numbers are identical to weight 1 to every digit. It seems the weight does not reach the second-component refit. This was not confirmed in the code.
- Centres of the minima, measured by mirror symmetry against the fitted ephemeris: primary -0.001 phase (about -33 s), secondary 0.000. The centres are not the problem.

Asymmetry of the minima is our explanation, **not checked**: the flat bottom sits to one side of the fitted centre and the residual wiggle is antisymmetric, as a symmetric profile leaves on an asymmetric minimum. The odd part of the minimum around its centre was not measured. The package warns "A1 is not much deeper than A2; period may be P/2" for this star. That was not examined either.

Decided: we stop here. NSV 807 is not timed from a Mikulášek fit that follows the minima. We go on to the other steps of the plan. Open for later: the per-start spread of the second-component trial, the odd part of the minima, the period question, and whether the weight reaches the second-component refit.

### Protocol of the timing for a real light curve (locked, 2026-11-09)

Status: locked on 2026-11-09. Marks as in the tiered-fit protocol: **Decided**, **Package**, **Suggestion**.

Code: `examples/nsv807_timing.py`, with two modes, `scan` and `time`. It reuses `scan_kind` and `window_is_covered` of `examples/window_scan_mik.py` and `eb.time` of the package. The template is the saved `chosen/` solution of the tiered fit of the same light curve, and the points are the same cleaned points.

Two rule sets are kept apart. The scan chooses a window and is a research step. The timing gives the result and keeps as many minima as the data allow.

**Scan rules.**

1. Decided: one fixed set of eclipses, chosen once at the widest window (the cap of 0.25 cycles). An eclipse belongs to the set when samples reach both ends of the window and no gap inside it is longer than 2 cadences. The same set is used at every window and for both cadences. Package: the rule is the one of `window_scan_mik.py`.
2. Decided: the half-window runs from D/10 in steps of D/10 up to 0.25 cycles. D is the width of that eclipse kind in the solution (`D1` primary, `D2` secondary). The outlier cut of the timing is off, so that a missing time can only come from the window.
3. Decided: the quantity is the rms of O-C on that fixed set, taken about a straight line fitted to the same times against the eclipse number. It is not taken about the file ephemeris, because the file period may have changed. The mean O-C is reported and does not choose the window.
4. Decided: the adopted window is the narrowest one whose rms is within rms/sqrt(2(N-1)) of the best rms, N being the number of eclipses in the set.
5. Decided: the window is scanned on the 120 s file and carried to the 600 s file as the same multiple of D, with the D of the 600 s solution.
6. Decided: the scan uses the same point rule as the timing (below), so that the two agree.

**Timing rules.**

1. Decided: the window is the adopted multiple of D, capped at 0.25 cycles.
2. Decided: a time is returned only when the fit can be identified. The timing fit has three free parameters, the shift and a baseline offset and slope. It needs more points than parameters and points on both sides of the predicted centre, so that both walls constrain the shift. The minimum is 4 points with at least 2 on each side. If every point lies on one side the shift cannot be told from the baseline, and no time is returned. The value 2 is our proposal. The package default is 3 per side.
3. Decided: the number of points is not a gate beyond that. The MAVKA floor of 6 points is not used for the Mikulášek fit.
4. Decided: nothing else rejects a time. The contact-wall check, the sigma cut and the clip of the package are off. Flags are kept instead: `walls_ok` (False when the wall check would have failed), `gap_in_window` (a gap in the window longer than 2 cadences), and the number of points and the degrees of freedom of each time.
5. Decided: edge eclipses are timed when they meet rule 2. Nothing is rejected for lying at the start or the end of the series.
6. Decided: there is no ladder of narrower windows. An eclipse that fails rule 2 at the adopted window gets no time and keeps its reason.
7. Decided: the two methods do not share conditions. Each method uses its own window and point rule. Where they are compared it is on the eclipses that both timed, as paired differences.
8. Package: the file epoch and period are used only to number the eclipses and to label the output. The predicted centres come from the fitted ephemeris of the solution.

**Whole-cycle timing of a real light curve.** Code `examples/nsv807_whole_cycle.py`, which calls `time_cycles` of `examples/batch_whole_cycle_time.py` (the method of the section "Whole-cycle timing (2026-11-07)"). Decided: one shift of the whole saved template per orbital cycle, points with phase in (n - 0.25, n + 0.75) of the fitted ephemeris, a baseline of degree 1 in each cycle. Decided: the point rule is the one of the window timing, at least 4 points with at least 2 on each side of the cycle centre, not the MAVKA floor of 6 and 3. The flag is `gap` when a gap longer than 2 cadences lies inside the cycle, and `ok` otherwise. The primary and the secondary of a cycle share one time shift, so their times are not independent. The outputs are as for the window timing: `eclipses.csv`, the two `.dat` files, a summary figure per kind and one figure per timed cycle in `whole_cycle/fits/`.

### Results of the timing of NSV 807 under the locked protocol (2026-11-09)

Window multiples 3.0 D (primary) and 3.2 D (secondary), from the scan of sector 30 at 120 s (adopted 3.0 D and 3.2 D, rms about a straight line 2.98 s and 3.54 s on 69 and 65 eclipses of the fixed set; figure `results/nsv807/s30_protocol/timing_scan/rms_vs_half_width.png`). Each folder `timing/` holds `eclipses.csv`, the two `.dat` files, a summary figure per kind and one figure per timed eclipse.

| Light curve | Primary: timed / no time | Secondary: timed / no time | Median sigma (s) primary, secondary | Gap flag | Wall check would fail |
|---|---|---|---|---|---|
| Sector 30, 120 s | 77 / 0 | 76 / 0 | 4.8, 4.5 | 5, 8 | 0, 0 |
| Sector 30, 600 s | 77 / 0 | 76 / 0 | 7.5, 7.8 | 5, 5 | 4, 3 |
| Sector 4, 1800 s | 62 / 3 | 60 / 6 | 6.3, 7.7 | 0, 2 | 62, 60 |
| Sector 97, 200 s (constant sigma 8.73 mmag) | 145 / 4 | 143 / 4 | 6.9, 5.5 | 4, 3 | 1, 2 |

Eclipses without a time are `too_few_points` or `one_side_uncovered`, that is, fewer than 4 points or fewer than 2 on a side. On sector 30 at 120 s the template follows the walls and the binned residual at the minimum is +6 mmag beside the centre and -7 mmag at the centre, the flat bottom missed. On the primary one time is about -96 s from the line, next to the data gap, and it carries the gap flag. The secondary figure and the flagged eclipses were not examined.

**Doubts about the median sigma on sector 4.**

- A median sigma of 6 to 8 s is small for 1800 s data. The windows hold as few as 4 points, and an eclipse of 4 points has 1 degree of freedom, so its formal sigma is not reliable.
- The sector 4 flux errors are tiny (0.47 mmag, checked on the MAST file), and the formal sigma of the shift scales with them. It is not known whether the real scatter of the times matches.
- The scatter of the times about a straight line was not compared with the formal sigma on this sector, so rms(O-C)/sigma is not known.
- Every timed eclipse of sector 4 would fail the contact-wall check. The times rest on a few points on the walls and are flagged only.
- The window is the 120 s multiple carried over. The 1800 s data were not scanned.
- The column `n_points` of `eclipses.csv` is counted by `nsv807_timing.py` around the fitted time. The package selects the points around the predicted centre, so the two counts can differ by a point. On sector 97 the smallest `n_points` is 3 for a secondary, below the minimum of 4 that the rule requires, which shows the difference.

### Results of the whole-cycle timing of NSV 807 (2026-11-09)

Median formal sigma and rms of O-C about a straight line, for the window timing of the section above and for the whole cycle. The rms is taken on all timed eclipses, with no clip, so a single bad time raises it (the window rms of the primary at sector 30 includes one time about -96 s from the line next to the data gap, and the sector 97 secondary one outlier). The paired difference is whole cycle minus window, on the eclipses that both timed.

| Light curve | Kind | Window: n, rms (s), median sigma (s) | Whole cycle: n, rms (s), median sigma (s) | Paired mean difference (s), scatter (s) |
|---|---|---|---|---|
| Sector 30, 120 s | primary | 77, 11.3, 4.8 | 76, 3.0, 2.7 | +25.9, 4.2 |
| Sector 30, 120 s | secondary | 76, 4.2, 4.5 | 76, 3.1, 2.7 | +7.0, 5.4 |
| Sector 30, 600 s | primary | 77, 12.1, 7.5 | 76, 4.2, 5.0 | +25.8, 4.7 |
| Sector 30, 600 s | secondary | 76, 5.2, 7.8 | 76, 4.2, 5.0 | +7.9, 6.2 |
| Sector 4, 1800 s | primary | 62, 5.3, 6.3 | 64, 3.7, 8.1 | +31.4, 5.6 |
| Sector 4, 1800 s | secondary | 60, 5.9, 7.7 | 64, 3.7, 8.1 | -7.3, 5.8 |
| Sector 97, 200 s | primary | 145, 4.7, 6.9 | 145, 5.9, 3.2 | +7.6, 6.0 |
| Sector 97, 200 s | secondary | 143, 52.8, 5.5 | 145, 5.9, 3.2 | +8.9, 53.3 |

The whole-cycle times of the primary sit later than the window times by about 26 s on sector 30 at both cadences and 31 s on sector 4, and the secondary by 7 s and -7 s. The offset is the same at 120 s and 600 s. A single shift for both minima is a compromise between them, so a difference of about 19 s between the primary and the secondary window times on sector 30 would show as a shared offset. That reading was not tested. The figures of the whole-cycle summary and of the per-cycle fits were looked at only for the sector 30 primary at 120 s, where the O-C is within about ±10 s with a slow wave of about 5 s.

### Protocol of the borrowed template for a real light curve (2026-11-09)

Code: `examples/nsv807_borrow.py`. It calls `eb.fit_borrowed` (`src/mikulasek_eb/borrow.py`) and then `run_time` of `examples/nsv807_timing.py`.

1. Decided: the reference is the chosen solution of the tiered fit of the 120 s light curve of sector 30. The target is sector 4 (1800 s), with all points and its own file errors.
2. Package: the borrowed values are fixed. The target keeps its own data layer: the baseline order and the gap split of its own tiered fit, and its own exposure (1440 s).
3. Decided: the epoch is always fitted. The variants differ by the parameter groups that are also free: none, depths, widths, depths and widths, and depths, widths and harmonics.
4. Decided: no variant is selected in the protocol. All variants are timed with the window timing of this section "Protocol of the timing", with the multiples 3.0 D (primary) and 3.2 D (secondary), D being the D of the variant.

### Results of the borrowed template on NSV 807, sector 4 (2026-11-09)

Reference sector 30 at 120 s, target sector 4 at 1800 s. Tables: `results/nsv807/s4_protocol/borrow/fits.csv` and `borrow/<variant>/timing/eclipses.csv`.

| Variant (free besides the epoch) | rms (mmag) | Timed primary / secondary | Median sigma primary / secondary (s) |
|---|---|---|---|
| none | 28.7 | 62 / 60 | 209 / 212 |
| widths | 24.7 | 62 / 60 | 221 / 182 |
| depths | 13.9 | 62 / 60 | 30 / 59 |
| depths, widths | 12.6 | 62 / 60 | 13 / 26 |
| depths, widths, harmonics | 8.5 | 62 / 60 | 7.1 / 7.9 |
| tiered fit of sector 4 itself | 8.5 | 62 / 60 | 6.3 / 7.7 |

Conclusions the data support. With the harmonics fixed to the 2020 values the template does not fit the 2018 light curve, and the formal sigma is 4 to 30 times larger. With the harmonics free the borrowed template reaches the rms, the timed counts and the formal sigma of the tiered fit of sector 4 itself. On this star and this sector borrowing gave no gain over the direct fit. What is not shown: whether the borrowed times agree with the direct times eclipse by eclipse, and their real scatter. The harmonics variant used 993 points, not 1025, through the clipping of the stored data layer. That was not examined.

### Protocol of the MAVKA timing of NSV 807 (locked, 2026-11-09)

Status: rules for the request in `agent_mailroom/LC_Modeller_to_Mavka_user_nsv807.md`. The methods are AP, WSAP and WSL.

1. Decided: the four light curves of NSV 807 (sector 30 at 120 s and 600 s, sector 4 at 1800 s, sector 97 at 200 s), SAP flux, as delivered. The period of the file is held fixed. The epoch of the file is not used to find minima.
2. Decided: the minima are found by the MAVKA fit itself, from the data. Numbering and the primary or secondary label are assigned afterwards by matching to the times of the other method.
3. Decided: the scan is run on the 120 s file by the same rules as in the Mikulášek protocol: one fixed set of eclipses chosen at the widest window, the rms of O-C about a straight line fitted to the same times, and the narrowest window within rms/sqrt(2(N-1)) of the best. The grid is D-dependent: the step is D/10 of the Mikulášek D of that kind (0.048450 primary, 0.049041 secondary, in cycles). W UMa is a contact star with no continuum, so the scan runs to 0.25 cycles, as for G5 and G6. The half-width is in cycles, one value per method and per kind. At the other cadences it is the same number in cycles, with no widening.
6. Decided: MAVKA times each minimum separately in its window. There is no whole-cycle timing for MAVKA. One primary and one secondary anchor time per sector are given for labelling, and the period of the file places the rest.
4. Decided: the point rule, the clean-up of outliers and the rejections are those of MAVKA, and each is written in the reply. They are not harmonised with the Mikulášek rules. Where the methods are compared it is on the eclipses that both timed.
5. Decided: a failed fit keeps its row with an empty time and the reason. Output as `.dat` files, one per method and kind (`jd  err_days  flag`), a scan figure, a figure of every fit, and the table.

### MAVKA timing of NSV 807 as it was performed (2026-11-09)

Reply: `agent_mailroom/Mavka_user_to_LC_Modeller_nsv807_cleaned.md`. Products: `/home/voz/projects/UPJS/Mavka_user_folder/output/nsv807_cleaned/`. A first run on the uncleaned files (`Mavka_user_to_LC_Modeller_nsv807.md`, `output/nsv807/`) was discarded, because sector 97 holds many outliers. It was not rewritten.

**Data.** Sectors 30 (120 s and 600 s) and 97 are the cleaned files of `results/nsv807/cleaned/` (the outlier step of the tiered-fit protocol: 97, 45 and 892 points removed). Sector 4 is the original file, 1025 points, because the outlier step cannot run on it. The fit is on the SAP flux with no weights. In sector 97, 1163 rows with non-finite SAP are left out, so 20909 values were fitted. The period of the file is held fixed and the epoch of the file is not used. Cycle 0 is the anchor of each kind that was given in the request.

**Departures from the request.**
1. The scan was run for AP and WSAP only, on the cleaned 120 s curve. WSL was not scanned. Its half-width is the WSAP half-width, and on the 120 s file its centres are the WSAP centres. This is written in `scan/chosen.txt`.
2. A window with fewer than 6 points was widened by one cadence in phase until it held 6, or until the two half-widths would add to more than 0.5. This is the choice of MAVKA, and such times are flagged `extended`. We did not ask for the widening and said the half-width is not widened.

**Windows.** Chosen on the 120 s file, in cycles, carried to the other files unchanged. They equal the numbers of the uncleaned curve.

| Method | Primary (D) | Secondary (D) | Set, primary / secondary |
|---|---|---|---|
| AP | 0.067830 (1.4) | 0.058848 (1.2) | 66 / 62 |
| WSAP | 0.053295 (1.1) | 0.044136 (0.9) | 65 / 62 |
| WSL | 0.053295 (borrowed) | 0.044136 (borrowed) | none |

**Times accepted** (primary / secondary):

| Light curve | AP | WSAP | WSL |
|---|---|---|---|
| Sector 30, 120 s | 77 / 76 | 77 / 78 | none |
| Sector 30, 600 s | 78 / 77 | 78 / 77 | none |
| Sector 4, 1800 s | 52 / 49 | 65 / 65 | 8 / 13 |
| Sector 97, 200 s | 149 / 145 | 151 / 150 | 146 / 146 |

WSL gives no time on sector 30, because the formal sigma is not finite at every window that holds 6 points. Windows were extended in 13 fits per method on sector 30 at 120 s, and in 36 (AP), 150 (WSAP and WSL) at 600 s. At 1800 s they were extended in 157 fits per method, and at 200 s in 40 per method. On sector 4 some extended AP times have a formal sigma of days to years. The largest such case is secondary E=57, with a sigma of about 8e4 days. The labels agree with depth: the median SAP at the minimum is lower for the primary. The timings were not compared with the Mikulášek times.

### Collected timings of NSV 807 (2026-11-09)

All `.dat` files of all methods are in `results/nsv807/timing_dat/` (50 files and `index.csv`), written by `examples/nsv807_collect_dat.py`, which only copies. The file names are `NSV807_<light curve>_<method>_<kind>.dat`. The methods are the Mikulášek window timing, the Mikulášek whole-cycle timing, the borrowed sector 4 variants (five) and MAVKA AP, WSAP and WSL. Every header line starts with `#` and names the star, the file, the minimum, the method, the window, the cleaning and the columns. A row without a time (MAVKA rejections) is written as a `#` line with its reason, so every remaining line is data: `jd  err_days  flag`.

### Protocol for V1383 Tau (locked, 2026-11-10)

Status: locked on 2026-11-10. Method only. Results go in separate sections. The star is V1383 Tau (TIC 388995726), a contact binary with a strong O'Connell effect. The protocols of NSV 807 apply, with the changes below.

Files, in `lc_shipper/data/V1383_Tau/`: sectors 42, 43, 44 (TESS-SPOC, 600 s) and sectors 70, 71 (TESS-SPOC, 200 s). Period of the file 0.35064 d.

**Mikulášek, changes against the NSV 807 protocols.**
1. Decided: SAP flux. The files have a real SAP error column, so there is no constant-sigma pilot.
2. Decided: every sector is fitted on its own. Sectors are not combined, because the level of the maxima and the depth of the primary differ between sectors.
3. Decided: the outlier step is the locked rule (400 phase bins, 5 sigma, local MAD, removed once, file period only, phase zero at `t.min()`), with the picture of the flagged points. If a sector has fewer than 5 points per bin the step stops and the case is reported.
4. Decided: the baseline is one polynomial per segment with a shared template, split at gaps longer than 0.3 d. (For NSV 807 the split was 0.5 d.)
5. Decided: the tiered search and the exposure rule are as locked. No second component is forced.
6. Decided: window timing only, with the timing rules of NSV 807 (at least 4 points, at least 2 on each side, flags kept, no ladder). There is no whole-cycle timing, because of the O'Connell effect. There is no template borrowing, because the cadences are sufficient.
7. Decided: the O-C scan is run once, on sector 70 (200 s), under the scan rules of NSV 807. The adopted window is carried to the other sectors as the same multiple of D, with the D of that sector's own solution.
8. Package: O-C between sectors uses the fitted ephemeris, not the five-digit file period.
9. Decided (2026-11-10, added after the sector 44 scan): a window is adopted only if it times the whole fixed set. Windows that time fewer eclipses of the set are not compared, because an rms over a few eclipses is not comparable with one over the whole set. This amends rule 4 of the scan rules of the timing protocol for this star. `nsv807_timing.py adopt` applies it to a saved scan.
10. Decided (2026-11-10): on a light curve of lower cadence only windows of at least the multiple adopted on the higher-cadence light curve are considered, and the adoption rule is applied among them (the narrowest within rms/sqrt(2(N-1)) of the best of those). A scan at lower cadence can therefore only widen the window. `nsv807_timing.py adopt --start-over-d PRIMARY SECONDARY` applies it.

**MAVKA (AP, WSAP, WSL).**
1. Decided: MAVKA receives the cleaned files (the points removed by item 3 above, with the lists of removed points), from the start.
2. Decided: the main scan is run on sector 70 (200 s) for all three methods, under the scan rules of the NSV 807 request: the step is D/10 of the Mikulášek D of that kind, up to 0.25 cycles, one fixed set of eclipses, the rms of O-C about a straight line, and the narrowest window within rms/sqrt(2(N-1)) of the best. The half-width in cycles is carried to the other sectors unchanged.
3. Decided: a second scan experiment, WSAP only, on sector 44 (600 s), with the same rules. It is labelled as a second experiment. Its purpose is to test the window carried from 200 s against a window chosen at 600 s.
4. Decided: MAVKA may widen a window that holds fewer than his minimum number of points. Every time obtained on a widened window is flagged `extended`, and the half-width actually used is written with it.
5. Decided: one primary and one secondary anchor time per sector are given for labelling. The period of the file places the rest. The file epoch is not used.
6. Decided: outputs are `.dat` files with `#` header lines and `#` lines for rows without a time, the tables, the scan figures and a figure of every fit.

### Protocol of the cleaning of V1383 Tau (2026-11-10)

Two cleaning rules exist for this star. Neither replaces the other. Both write a removed-point list, a folded figure, a time figure with the SAP background below it on the same time axis, and the original VOTable without the removed points (`examples/nsv807_cleaned_files.py`, `clean_one`).

1. Folded-curve rule, `examples/v1383_clean.py`. It is the outlier step of the locked tiered-fit protocol (5 sigma from the median of the phase bin, local MAD, removed once, file period only). The number of bins is 400 for the 200 s files and 200 for the 600 s files, because the 600 s files do not fill 400 bins evenly (see the results). Output in `results/v1383_tau/cleaned/`.
2. Background rule, `examples/v1383_clean_background.py`. A point is removed when `sap_bkg / median(sap_bkg)` of its file exceeds a ratio, default 10. It does not look at the magnitude. It is an option for the scattered-light episodes, not a rule for every case. Output in `results/v1383_tau/cleaned_bkg/`.
3. Diagnostic, `examples/v1383_bkg_diagnostic.py`: the residual from the bin median of the folded curve against `sap_bkg / median`, with the robust scatter in ratio bins. Nothing is removed. Figure `results/v1383_tau/bkg_ratio_diagnostic.png`.

### Results of the cleaning of V1383 Tau (2026-11-10)

Files: sectors 42, 43, 44 (600 s) and 70, 71 (200 s), SAP flux.

**Bins of the folded-curve rule.** At 400 bins the three 600 s files have bins with 1 to 3 points, because the cadence (0.0198 cycles) is much longer than the bin (0.0025 cycles). The step stopped. At 200 bins the smallest bin holds 9 to 11 points. The 200 s files have at least 23 points per bin at 400 bins.

**Folded-curve rule, removed points.**

| Sector | Cadence | Removed | Left |
|---|---|---|---|
| 42 | 600 s | 13 of 3440 | 3427 |
| 43 | 600 s | 8 of 3327 | 3319 |
| 44 | 600 s | 18 of 3248 | 3230 |
| 70 | 200 s | 68 of 10539 | 10471 |
| 71 | 200 s | 107 of 10596 | 10489 |

In time (sectors 70, 71 and 44 were looked at, sectors 42 and 43 were not), almost all removed points of the 200 s files lie in short stretches at the starts and ends of orbits, where the SAP background is 100000 to 750000 e/s against about 1100 to 1600 in quiet stretches. Sector 44 has isolated removed points at quiet background, a few on the steep wall. In sector 71 the points kept around days 13 to 14.4, where the background stays at tens of thousands, scatter widely, so the rule left part of that episode in the file. Whether any removed point is a real point of the minimum was not decided from these pictures.

**Background in sigma, rejected.** Removing points more than 5 robust sigma above the median background removed 6 to 23 percent of each file (778, 246, 210, 1545 and 2240 points), including whole stretches where the background is only some tens to hundreds of e/s above quiet. The quiet background is stable (sigma 39 to 59 e/s in four of five files), so any multiple of it is far below a level that harms the photometry. This version is not kept.

**Residual against the background ratio** (robust sigma of the residual in mmag):

| Sector | Quiet (ratio near 1) | Behaviour at larger ratio |
|---|---|---|
| 42, 43, 44 (600 s) | 4 to 5 | 4 to 7 up to ratio 30 |
| 70 (200 s) | 6 to 7 | 9 to 14 from ratio 1.5 to 2; 19 above 30; 95 above 100 |
| 71 (200 s) | 4 to 6 | 8 at 5 to 10; 10 at 10 to 30; 12 at 30 to 100; 48 above 100 |

The residual is taken from the bin median, so it includes real changes of the curve between cycles and the bad points themselves. The rise of the scatter is a hint, not proof, that the points are bad. The 600 s files show no rise within the range they cover. Sector 71 has a stripe of about 600 points near a ratio of 17.

**Background rule at ratio 10.**

| Sector | Cadence | Removed | Cut (e/s) |
|---|---|---|---|
| 42 | 600 s | 161 of 3440 | 16198 |
| 43 | 600 s | 40 of 3327 | 14406 |
| 44 | 600 s | 67 of 3248 | 14162 |
| 70 | 200 s | 335 of 10539 | 10996 |
| 71 | 200 s | 793 of 10596 | 12317 |

The ratio of 10 removes more than the folded-curve rule on every file, 3 to 8 percent of the 200 s files. On the 600 s files the diagnostic shows no need for a cut, yet the cut removes 40 to 161 points, among them the sector 44 rise to 47000 e/s where the magnitudes looked normal. The new pictures at ratio 10 and the stripe in sector 71 were not looked at.

### Tiered fit of V1383 Tau on the background-cleaned files (2026-11-10)

Method. The files are those of `results/v1383_tau/cleaned_bkg/` (background rule, ratio 10), chosen by the researcher. The protocol of "Protocol for V1383 Tau (locked, 2026-11-10)" applies: SAP flux, real errors, each sector on its own, tiered search of `examples/nsv807_s30_tiered.py` with `--skip-outliers` (the files are already cleaned) and `--gap-days 0.3`. The option `--gap-days` was added to that script. It passes the value to `FeatureSettings`, and the default is the package value of 0.5 d. The runs are in `results/v1383_tau/<sector>_fit_bkg/`, one folder per sector.

Results.

| Sector | Cadence | Points | Chosen | rms (mmag) | Baseline pieces | Exposure integrated (s) |
|---|---|---|---|---|---|---|
| 42 | 600 s | 3279 | p1_h123_b0 | 5.26 | 4 | 480 |
| 43 | 600 s | 3287 | p1_h123_b1 | 5.30 | 3 | 480 |
| 44 | 600 s | 3181 | p1_h123_b1 | 5.94 | 3 | 480 |
| 70 | 200 s | 10204 | p1_h123_b0_C | 7.77 | 2 | 160 |
| 71 | 200 s | 9803 | p1_h123_b0 | 7.74 | 2 | 160 |

All fits have status ok. Every search stopped after two steps that gained less than the decisive BIC difference, with profile 1 and the harmonic set 1, 2, 3. On sector 70 the free-C variant was chosen. The two-component trial was tried on sectors 70 and 71 and rejected (rms gain below 0.25). Sectors 42, 43 and 44 report a rank-deficient Jacobian (rank 12 of 14 on 42 and 43, 10 of 11 on 44), so their parameter errors are not reliable. The rms is 5 to 8 mmag against a formal photometric error of about 1.7 mmag at 600 s and 2.6 mmag at 200 s. The folded picture of sector 70 follows both minima and the O'Connell difference of the maxima. The residuals are structured, with an amplitude of about 10 mmag, repeating over the phase. That structure was not analysed. The figures of the other sectors were not looked at.

### O-C scan of V1383 Tau, sector 70 (2026-11-10)

Method: the scan rules of "Protocol of the timing for a real light curve" (`examples/nsv807_timing.py scan`), on the background-cleaned file of sector 70 (200 s) and its saved solution `p1_h123_b0_C`. Output in `results/v1383_tau/s70_fit_bkg/timing_scan/`: the scan tables, `adopted.csv` and `rms_vs_half_width.png`.

Results.

| Kind | D (cycles) | Adopted half-window | Set N | rms of O-C at the adopted window (s) | Best rms (s) |
|---|---|---|---|---|---|
| primary | 0.04615 | 0.1523 cycles = 3.3 D | 64 | 6.96 | 6.40 |
| secondary | 0.05589 | 0.2068 cycles = 3.7 D | 67 | 17.50 | 16.13 |

The primary rms falls steeply up to about 2 D and then slowly, from 15 s at 1.5 D to 6.4 s at 4 D. The secondary rms is larger at every window and falls from 40 s at 1.4 D to 16 s at 4 D. The secondary scan ends at 4.4 D because the cap of 0.25 cycles is reached. Below about 1 D (secondary) and 0.8 D (primary) no eclipse could be timed. The rms of O-C over the formal sigma is between 0.5 and 1.2 for the primary and between 0.5 and 0.8 for the secondary beyond 1 D, so the formal errors are not smaller than the scatter. The reason for the larger scatter of the secondary was not investigated.

### Window timing of V1383 Tau at the adopted windows (2026-11-10)

Method: `examples/nsv807_timing.py time` with the multiples of the sector 70 scan, 3.3 D (primary) and 3.7 D (secondary), each sector with its own D, on the background-cleaned files and the saved solutions of `<sector>_fit_bkg`. The timing rules are those of "Protocol of the timing for a real light curve". Output in `results/v1383_tau/<sector>_fit_bkg/timing/`: `eclipses.csv`, the two `.dat` files, a summary figure per kind and a figure of every timed eclipse.

Results (timed eclipses, median formal sigma, median number of points in the window, eclipses without a time and their reason):

| Sector | Cadence | Primary: timed, sigma (s), points | Secondary: timed, sigma (s), points | Not timed |
|---|---|---|---|---|
| 42 | 600 s | 68, 9.2, 16 | 65, 38.9, 17 | 1 secondary (one side uncovered) |
| 43 | 600 s | 66, 10.6, 15 | 64, 26.0, 17 | 1 primary (too few points), 1 secondary (one side uncovered) |
| 44 | 600 s | 64, 9.1, 15 | 63, 27.8, 17 | 2 primary (one side uncovered) |
| 70 | 200 s | 68, 6.6, 46 | 67, 22.1, 63 | 1 primary, 1 secondary (one side uncovered) |
| 71 | 200 s | 65, 6.5, 47 | 64, 24.3, 48 | 2 primary, 1 secondary (one side uncovered) |

The half-window is 0.152 to 0.154 cycles for the primary in all sectors. For the secondary it is 0.207 cycles on sector 70, 0.159 on sector 71 and 0.167 to 0.168 on the 600 s sectors, because D is the D of each sector's own solution. The gap flag is set on 0 to 5 timed eclipses of each kind and the wall check is false on at most 3. The secondary sigma is three to four times the primary sigma. The O-C figure of the sector 70 secondary shows no trend, with scatter of about 15 to 20 s and one time about 70 s late near JD 2460222. The residuals of the template in that window show a repeating pattern of about ±5 mmag. Only the sector 70 secondary figure was looked at. The 600 s files have about 16 points in a window, so their sigma should be read with the same doubts as for sector 4 of NSV 807. The scatter of the real times was not compared with the formal sigma, and the timings of the sectors were not compared with each other.

### O-C scan of V1383 Tau, sector 44 (600 s) with Mikulášek (2026-11-10)

Method: as for sector 70, on the 600 s file of sector 44 and its saved solution `p1_h123_b1`. It is a second scan experiment. Its windows were not used for the timing above. Output in `results/v1383_tau/s44_fit_bkg/timing_scan/`.

Results. The fixed set holds 57 primary and 61 secondary eclipses. The scan took about 12 minutes.

| Kind | D (cycles) | Adopted by the rule | rms (s) | N used by the rule |
|---|---|---|---|---|
| primary | 0.04663 | 0.0699 cycles = 1.5 D | 9.45 | 4 |
| secondary | 0.04526 | 0.1901 cycles = 4.2 D | 19.5 | 61 |

The secondary rms falls from 59 s at 1.9 D to 18 s at 5.5 D, and rms/sigma is 0.65 to 0.95 beyond 3 D. The primary adoption is an artefact of the rule. At 1.5 D only 4 of the 57 eclipses of the set are timed (the others have fewer than 4 points or points on one side only), and the rms of 9.45 s is taken over those 4, so it is the smallest value on the scan. With 57 eclipses timed (from 1.9 D) the primary rms is 22 s at 1.9 D and 26 to 33 s above 2 D, falling slowly to 26 s at 4.9 D. The rule compares rms values taken over different numbers of eclipses. For sector 70 this did not matter, because the adopted window timed 64 eclipses. The rule was then amended (item 9 of the V1383 Tau protocol): only windows that time the whole set are compared. Applied to the saved scans without repeating them, it gives for sector 44 a primary window of 0.0886 cycles = 1.9 D (rms 22.2 s, N 57) and an unchanged secondary of 4.2 D. The sector 70 windows are unchanged (3.3 D and 3.7 D). The sector 44 primary rms of 22.2 s is the first window of the scan that times all 57 eclipses, and it is the lowest value of the full-set windows. Beyond it the rms is higher (26 to 33 s). The adopted 1.9 D is therefore at the edge where timing becomes possible, so the primary scan has no interior minimum at 600 s. This was not investigated further. The 1.9 D window was rejected as too narrow (the primary window at 200 s is 3.3 D). With rule 10 (windows from 3.3 D for the primary and 3.7 D for the secondary) the sector 44 windows are 3.6 D (0.1679 cycles, rms 27.1 s, best 25.4 s, N 57) for the primary and 4.2 D (0.1901 cycles, rms 19.5 s, N 61) for the secondary. Above 3.3 D the primary rms is flat, 25 to 29 s. The 600 s timings made before this rule used 3.3 D and 3.7 D. They are not yet repeated at these windows. The scan figure is `rms_vs_half_width.png` in that folder.

### Window timing of the 600 s sectors of V1383 Tau at the windows of rule 10 (2026-11-10)

Method: `nsv807_timing.py time` with 3.6 D (primary) and 4.2 D (secondary), the windows of sector 44 under rule 10, each sector with its own D, on sectors 42, 43 and 44. The output is in `timing_w36_42/` of each run folder. The earlier set in `timing/` (3.3 D and 3.7 D) was not rewritten and its files were checked unchanged.

Results. The counts of timed eclipses are the same as in the first set (68 / 65, 66 / 64 and 64 / 63 for primary / secondary). The half-windows are 0.168 cycles (primary) and 0.189 to 0.190 cycles (secondary). The paired difference is new minus first set, on the eclipses timed in both:

| Sector | Kind | Median sigma first / new (s) | Paired mean difference (s) | Scatter (s) | Largest absolute (s) |
|---|---|---|---|---|---|
| 42 | primary | 9.2 / 8.8 | +1.0 | 3.3 | 19 |
| 42 | secondary | 38.9 / 37.8 | +8.9 | 10.5 | 32 |
| 43 | primary | 10.6 / 10.1 | +3.2 | 2.6 | 9 |
| 43 | secondary | 26.0 / 25.9 | +12.5 | 8.4 | 32 |
| 44 | primary | 9.1 / 9.3 | +2.8 | 3.1 | 13 |
| 44 | secondary | 27.8 / 27.0 | +12.0 | 6.9 | 29 |

The wider windows change the formal sigma by less than 5 percent. The primary times move by 1 to 3 s on average. The secondary times move later by 9 to 12 s, with a scatter of 7 to 10 s, which is comparable to the secondary sigma. Which set is closer to the true minima is not known. The pictures of the new set were not looked at, and the `.dat` files of the new set are not yet collected.

### Collected timings of V1383 Tau (2026-11-10)

`results/v1383_tau/timing_dat/` holds the ten window-timing `.dat` files (five sectors, two kinds) and `index.csv`, written by `examples/v1383_collect_dat.py`. Names: `V1383Tau_<sector>_mikulasek_window_<kind>.dat`. Every header line starts with `#` and gives the star, file, minimum, method, window, cleaning, template, flags and columns. All other lines are data. MAVKA files are not yet included.

### Request to MAVKA for V1383 Tau (2026-11-10)

The request is `agent_mailroom/LC_Modeller_to_Mavka_user_v1383_tau.md`, with the reply expected in `Mavka_user_to_LC_Modeller_v1383_tau.md` and the products in `Mavka_user_folder/output/v1383_tau/`. It follows "Protocol for V1383 Tau" (MAVKA part) with these additions, which were decided after the protocol was locked:

1. The files are the background-cleaned files (`results/v1383_tau/cleaned_bkg/`), given from the start.
2. Scan 1 is on sector 70 for AP, WSAP and WSL: step D/10 of the Mikulášek D of that kind (0.04615 primary, 0.05589 secondary), up to 0.25 cycles. Only half-widths that time the whole fixed set take part in the choice (item 9 of the Mikulášek protocol, applied to MAVKA too).
3. Scan 2 is on sector 44 for WSAP only (D 0.04663 primary, 0.04526 secondary), with the whole grid written. Only half-widths from the WSAP half-width of scan 1 upward are considered (item 10).
4. WSAP is timed twice on the 600 s sectors, with the window of scan 1 and with the window of scan 2, and the two sets are labelled apart (`WSAP2`). AP and WSL use the window of scan 1.
5. Widening is allowed and flagged `extended`. Anchors are one primary and one secondary JD per sector. The `.dat` files start with `#` comment lines.

### First MAVKA run on V1383 Tau rejected (2026-11-10)

The first reply (`Mavka_user_to_LC_Modeller_v1383_tau.md`, `output/v1383_tau/`) is not used. MAVKA placed the next window centre at the fitted time of the previous eclipse plus one period whenever that fit fell inside the window, and the scan was run at a half-width of 0.25 cycles. On sector 70 a WSAP fit that was 6828 s early (0.19 periods) moved the next centre, and WSAP stood on the maximum for cycles 54 to 62, with model fluxes of 2351 to 2374 ct/s against 1380 to 1390 ct/s at the minimum. AP and WSL did not follow it. The centres of the three methods were therefore not the same and the run does not compare the methods. The anchor plus n periods, with the period of the file, agrees with the Mikulášek fitted times to about 30 s typically, and to 226 s at most, over a sector. The window half-widths are 4000 to 7500 s. A new request (`LC_Modeller_to_Mavka_user_v1383_tau_centres.md`) asks for fixed centres `anchor + n x period` in every window, no centre taken from a fitted time, a rejection of a fit more than 0.1 period from its centre (`far_from_centre`), and a check table for sector 70 before the full run. The new products go to `output/v1383_tau_2/`.

### WSL on V1383 Tau: all times rejected, request to keep them (2026-11-10)

On the second run (`output/v1383_tau_2/`, note `Mavka_user_to_LC_Modeller_v1383_tau_2_note.md`) WSL gave no accepted time on any of the five files, at the windows of its own scan (0.152295 and 0.190026 cycles) and after a refit at the WSAP half-widths (0.050765 and 0.067068 cycles). Of 145 to 148 rows per file, 130 to 140 were rejected because the formal sigma was not finite and 8 to 13 because they held fewer than 6 points (one `maxfev` in sector 71). On NSV 807 sector 30 the same reason removed every WSL time. The request `LC_Modeller_to_Mavka_user_v1383_tau_wsl_sigma.md` asks for WSL again on all five files at the WSAP windows, with one new rule for V1383 Tau only: a converged fit whose formal sigma is not finite keeps its time, the error is written as `None` (the agreed way to say that the formal error is infinite), and the time is flagged `sigma_none`. The other rejections stay. The new run goes to `output/v1383_tau_2/WSL_sigma_none/`, and the earlier tables are kept. Such times have no usable formal error, so they cannot be weighted in a comparison and enter it only as times.

### WSL on V1383 Tau with times kept at an infinite sigma (2026-11-10)

Reply: `Mavka_user_to_LC_Modeller_v1383_tau_wsl_sigma.md`, products in `Mavka_user_folder/output/v1383_tau_2/WSL_sigma_none/`. WSL was run again on the five files at the WSAP half-widths (0.050765 primary, 0.067068 secondary), centres at anchor + n x period, with the rule that a converged fit with a non-finite formal sigma keeps its time, the error is `None` and the flag is `sigma_none` (`extended+sigma_none` on a widened window). The limit for `far_from_centre` was 0.1 period (3029.5 s).

Counts (kept with error `None`, rejected): s42: 68, 6 (primary), 65, 9 (secondary); s43: 66, 7 and 64, 8; s44: 65, 7 and 63, 8; s70: 69, 5 and 67, 7; s71: 66, 7 and 64, 9. No WSL time has a finite sigma. Rejections are `far_from_centre` (1 to 5 per kind), fewer than 6 points (3 to 7), and one `maxfev` on s71. Time minus centre, computed from the saved tables with the anchors and the period of the file: for 53 to 69 of the kept times of each kind in sectors 42, 70 and 71 and for the primary of 43 and 44 the difference is under 100 s (median 19 to 50 s). The secondary of sector 43 has a median of 158 s with 62 of 64 times between 100 and 300 s, and the secondary of sector 44 has a median of 109 s, so the WSL secondary minima lie systematically 100 to 200 s from the anchor prediction in those sectors. The largest differences are +2836 s (s70 primary, cycle 36, `extended+sigma_none`, at the last point of a rising stretch), +432 s (s42 primary, cycle 30) and -330 s (s44 primary, cycle 32). The 0.1 period limit let the 2836 s time through, so that limit is too loose for a single time that is plainly not a minimum. These WSL times were not compared with the Mikulášek times.

### Results of the MAVKA run on V1383 Tau and comparison of all methods (2026-11-10)

Reply: `Mavka_user_to_LC_Modeller_v1383_tau_2.md`, products in `Mavka_user_folder/output/v1383_tau_2/`. All collected `.dat` files (52, Mikulášek first set, Mikulášek second set for 600 s, MAVKA AP, WSAP, WSAP2 and WSL with the `None` rule) are in `results/v1383_tau/timing_dat/`. The comparison, `examples/v1383_compare.py` (reads only those files), is in `results/v1383_tau/method_comparison.csv`.

**What MAVKA did.** Centres at anchor + n x period, rejection of a time more than 3029.5 s from the centre (39 times, differences of 4250 to 282334 s, most on AP and WSAP at 600 s), widening to 6 points flagged `extended`. His scan 1 on sector 70 chose the same half-width for AP and WSAP, 0.050765 (1.1 D) and 0.067068 (1.2 D), and for WSL 0.152295 and 0.190026 (the WSL windows were replaced by the WSAP ones). His scan 2 on sector 44 (WSAP only, from the scan 1 half-width upward) chose 0.083934 (1.8 D, primary) and 0.085994 (1.9 D, secondary). The window of the Mikulášek timing is 3.3 D and 3.7 D (first set), and 3.6 D and 4.2 D (second set), so the MAVKA windows are about a third to a half of them. The sector 70 centre check gave at most 107 s between centre and the Mikulášek time. One accepted WSAP time (s70, primary, cycle 36) lies +2671 s from its centre on six points of a rising stretch, and the same cycle gives +2836 s in WSL. No WSL formal sigma is finite.

**Numbers of times (primary / secondary), accepted.**

| Sector | Mikulášek | AP | WSAP | WSL | WSAP2 |
|---|---|---|---|---|---|
| 42 | 68 / 65 | 68 / 64 | 65 / 64 | 68 / 65 | 68 / 64 |
| 43 | 66 / 64 | 66 / 64 | 66 / 64 | 66 / 64 | 66 / 64 |
| 44 | 64 / 63 | 65 / 63 | 64 / 62 | 65 / 63 | 65 / 63 |
| 70 | 68 / 67 | 67 / 67 | 69 / 67 | 69 / 67 | |
| 71 | 65 / 64 | 65 / 64 | 66 / 64 | 66 / 64 | |

**Rms of O-C about an unweighted straight line (s), and the median of the time minus the Mikulášek first-set time on the common eclipses (s, MAD scatter in brackets).** The rms is not clipped, except where noted.

| Sector | Kind | Mik. rms | AP rms, diff | WSAP rms, diff | WSL rms, diff |
|---|---|---|---|---|---|
| 42 | primary | 8.6 | 66, +7 (8) | 50, +15 (6) | 50, +17 (9) |
| 42 | secondary | 35.6 | 35, -30 (23) | 39, -31 (28) | 40, -33 (28) |
| 43 | primary | 6.5 | 8.2, +26 (9) | 7.7, +40 (12) | 8.5, +41 (11) |
| 43 | secondary | 22.8 | 29.6, -116 (22) | 25.8, -117 (17) | 26.1, -117 (18) |
| 44 | primary | 27.4 | 35.8, +34 (7) | 21.3, +45 (12) | 50.8, +46 (10) |
| 44 | secondary | 23.3 | 42.2, -116 (34) | 27.4, -128 (26) | 26.0, -129 (24) |
| 70 | primary | 6.9 | 8.1, +12 (8) | 316 (14.4 clipped), +13 (8) | 335 (7.6 clipped), +14 (8) |
| 70 | secondary | 17.5 | 16.9, -25 (13) | 17.2, -22 (12) | 17.2, -22 (11) |
| 71 | primary | 5.6 | 8.1, +6 (7) | 16.5, +7 (6) | 30.9, +7 (6) |
| 71 | secondary | 22.7 | 30.3, +61 (25) | 25.0, +66 (21) | 24.8, +66 (24) |

Median formal sigma of the secondary (s): Mikulášek 22 to 39, AP 23 to 41, WSAP 15.5 to 20, WSL none. For the primary: Mikulášek 6.5 to 10.6, AP 7.3 to 9.0, WSAP 5.1 to 9.8.

**Conclusions the data support.**
1. AP, WSAP and WSL agree with one another. Their times differ from each other by a few seconds on the common eclipses, and each differs from the Mikulášek time by the same amount. The difference is not random: the MAVKA primary times are later than the Mikulášek ones by 6 to 66 s, and the secondary times earlier by 22 to 144 s on sectors 42 to 70, but later by 61 to 66 s on sector 71. It is largest on the 600 s sectors 43 and 44 (secondary about -116 to -129 s, primary +26 to +46 s). The second WSAP set (wider window) moves the MAVKA times a further 8 to 16 s away from the Mikulášek ones. The three MAVKA methods are local fits of the bottom of the minimum, and the Mikulášek times follow a template over a wide window. The asymmetric minimum is the probable cause (this star has a strong O'Connell effect), but this was not tested. This difference of 100 s or more at 600 s is larger than any formal sigma, so absolute minimum times of the two kinds of method cannot be mixed in one O-C diagram. Which time is the true minimum is not decided by this comparison: there is no truth for this star.
2. The scatter about a straight line is about the same for all methods in the secondary at 200 s (17 s) and the primary at 200 s (7 to 8 s for Mikulášek and AP, after removing one bad time for WSAP and WSL). At 600 s the Mikulášek times are smoother on sectors 42 (primary 8.6 s against 50 to 66 s) and 44 (secondary 23 s against 26 to 42 s), and MAVKA is smoother in the primary of sector 44 (WSAP 21 s against 27 s). On sector 42 the MAVKA primary scatter comes from times close to the data gap (cycle 30, +422 to +507 s from the centre, flagged `gap` and `edge`), where the Mikulášek times do not show it. The rms is not clipped, and a few bad times dominate several entries, so these numbers are rough.
3. The MAVKA formal sigma of the secondary (WSAP 15 to 20 s) is below the scatter about the line (17 to 39 s) on every sector, so it is too small. For the Mikulášek secondary the sigma (22 to 39 s) is about the scatter (17.5 to 35.6 s). The WSL sigma is not finite for every time. It is a failure of the fit to provide an error.
4. The windows are different: MAVKA at 1.1 D to 1.9 D, Mikulášek at 3.3 D to 4.2 D. The methods were run under their own conditions, as agreed. This comparison does not separate the effect of the window from the effect of the method.
5. Both methods time the minimum in nearly every cycle. The numbers of times are within 1 to 4 of one another on every sector. A few bad times are accepted: s70 WSAP and WSL primary cycle 36 (+2671 s and +2836 s from the centre, on a gap) and several times on sector 42 near its gap. The rejection limit of 0.1 period is too loose for them.

What is not shown: the cause of the systematic differences, the effect of the wider MAVKA windows on the offset beyond the WSAP2 test, the agreement with an independent time such as a spectroscopic ephemeris, the 200 s sectors 70 and 71 at the Mikulášek second-set windows, and the pictures of the individual fits (only the sector 70 secondary Mikulášek figure and the MAVKA folded figures were mentioned in the reports and none of the MAVKA figures was looked at).

### Protocol for SV Cen (locked, 2026-11-10)

Star SV Cen, TIC 324101484. Aim: see how the Mikulášek method deals with a fast period change, not to compare it with MAVKA (MAVKA is postponed). Expected rate about -2.1 s/yr, period about 1.658 d. Data: `lc_shipper/data/SV_Cen`, QLP SAP flux, sectors 10, 11 (1800 s), 37, 38 (600 s), 64, 100 (200 s), 99 (600 s), 2531 d in all. The SPOC file of sector 64 is ignored. About 16 cycles per sector, about 105 eclipses of each kind in all.

1. Decided: no template borrowing between sectors and no whole-cycle timing. The shape changes between sectors (maxima differ by about 0.25 mag, a variable O'Connell effect), and the cadences differ.
2. Decided: the cleaning is per sector. The background ratio rule only (ratio 10). There is no folded rule. Pictures with the outliers marked, folded and unfolded with the background, before anything is fitted.
3. Decided: flux errors. The files have no `sap_flux_err`, and `kspsap_flux_err` is judged unrealistic. The sigma is taken from the residual of the fit, per sector.
4. Decided: two kinds of fit.
   - A. Sector-specific: the tiered fit of the locked protocol on each sector, with the timing of each eclipse. This gives the main results and the O-C.
   - B. One template for all sectors: one shape, no tier search, with the period and a free `Pdot` fitted to all sectors together. The shape is fixed once, from fit A of sector 64. This fit runs second.
5. Decided: timing is window timing with the rules of V1383 Tau (at least 4 points, at least 2 on each side, flags kept). The window scan is on sector 64 (200 s) with a straight-line reference per sector, and the adopted window is applied to the other sectors, with the rule for lower cadences (item 10 of V1383 Tau).
6. Decided (2026-11-10): the window is locked in phase (cycles), not as a multiple of D. D is a width of the fitted profile, so the same multiple of D is a different fraction of the orbit for each minimum and each sector. On sector 64 the fit gives D1 = 0.0385 and D2 = 0.1239 cycles, and 1.6 D is 0.062 cycles for the primary but 0.198 cycles for the secondary. The locked half-width is that wider one: 0.198 cycles, taken from the sector 64 secondary at 1.6 D, and it is used for the primary and the secondary of every sector. The timings made at 1.6 D of each sector's own profile are not kept.
6. Decided: the O-C of all sectors is first drawn about the straight line, then a quadratic is fitted to all the times. Primary and secondary are kept apart. The quadratic gives Pdot from the O-C.
7. Decided: the local period of each sector from fit A is extracted and plotted against time, as an extra view.
8. Decided: Pdot from the O-C (item 6), Pdot from the one-template fit (item 4B) and the slope of the local periods (item 7) are shown side by side. The cause of a difference is not concluded without a test. A variable O'Connell effect can shift the times and mimic or hide a period change.
9. Result, not decided. The tiered search stopped at its first step on sector 64. That step fits only profile 1, and every profile 1 candidate ended with status ``bound-active`` (the secondary amplitude A2 at 8.8 against a bound of 4.0), so the search raised an error and never reached the profile step, which is where profiles 2, 3 and 4 are tried. The search was trapped by its own order of steps. A diagnostic run of profiles 1 to 4 on sector 64, outside the search, showed profile 3 as the only one where the fit finishes with status ``ok`` (rms 13.7 mmag). The option ``--start-profile`` was added to ``nsv807_s30_tiered.py`` and the search was rerun on sector 64 with profile 3 as the start. It chose ``p3_h123_b1``, status ``ok``, rms 13.72 mmag, and profiles 1 and 2 failed again in the profile step. The residuals still reach about 75 mmag at the walls of the primary and the template leaves structure. The constant sigma of that fit is 9.23 mmag, from the scatter about the phase-bin medians, because the pilot fit of the protocol also ended ``bound-active``. The ephemeris of this fit has ``Pdot`` fixed at 0. The other sectors are not fitted yet.
10. Before the fits, the expected size is: quadratic O-C about 3 h across the span, inside one sector about 1 s, and a cycle count that is unambiguous across the gaps. These are estimates from the quoted rate, to be checked against the data.

### Tiered fits of SV Cen, per sector (2026-11-10)

Method. The files are `results/sv_cen/cleaned/` (folded-curve rule, 5 sigma, the three flagged points of sector 64 restored). The search is `examples/nsv807_s30_tiered.py` with `--skip-outliers`, `--no-flux-err`, `--sigma-from-bins` and `--start-profile`. The constant sigma is the robust scatter about the phase-bin medians, because the pilot fit of the protocol ends `bound-active` on these files. The number of bins is 400 at 200 s, 200 at 600 s and 70 at 1800 s. The start profile is 3, except sectors 99 and 100, where profiles 1, 2 and 3 all end `bound-active` at the first step and the search was rerun with start profile 4. `Pdot` is fixed at 0 in every fit. The figures were looked at for sectors 64, 99 and 100 only.

| Sector | Cadence | Chosen | rms (mmag) | Sigma (mmag) | Baseline pieces | BIC gains | Stopped early |
|---|---|---|---|---|---|---|---|
| 10 | 1800 s | p3_h123_b0 | 17.69 | 22.68 | 1 |  | yes |
| 11 | 1800 s | p3_h123_b1 | 15.46 | 30.73 | 1 |  | yes |
| 37 | 600 s | p3_h2_b1 | 13.95 | 20.31 | 1 |  | no |
| 38 | 600 s | p3_h123_b0 | 15.18 | 10.55 | 1 |  | yes |
| 64 | 200 s | p3_h123_b1 | 13.72 | 9.23 | 1 | 5.7, 0.0 | yes |
| 99 | 600 s | p4_h123_b1 | 10.45 | 16.03 | 2 | 632.7, 0.0, 0.0 | no |
| 100 | 200 s | p4_h123_b0 | 13.06 | 10.95 | 1 | 0.0, 0.0 | yes |

All chosen fits have status ok. On sectors 99 and 100 the other three profiles end `bound-active` (an amplitude of 3.4 to 17.9 against a bound of about 2.7). Sector 99 has a gap of about 2 days and two baseline pieces, and its baseline order 1 beats order 0 by 633 in BIC. Sector 100 keeps baseline order 0: order 1 has the lower rms (12.51 mmag) but a BIC gain of 6.2, below the decisive difference of 10. Sector 37 drops to the single harmonic 2. The residuals of sectors 64, 99 and 100 are waves of 20 to 75 mmag that repeat over the phase, so the template leaves structure. The figures of sectors 99 and 100 show a second component on the primary. The runs are in `results/sv_cen/s<sector>_fit/`, and `s99_fit_p4/` and `s100_fit_p4/` for the two that needed profile 4.

### Window timing of SV Cen at the locked phase window (2026-11-10)

The window is item 6 of "Protocol for SV Cen": one half-width of 0.19826 cycles for both minima of every sector. It is 1.6 times D2 of the sector 64 solution (`D2 = 0.12391` cycles). The earlier timing, made at 1.6 D of each sector's own profile, was deleted, including its figures and its `timing_dat`.

Run from `mikulasek_eb`, with `MPLBACKEND=Agg`. `examples/nsv807_timing.py` was not changed. `examples/sv_cen_timing.py time` calls it. The script takes the window as a multiple of D, so the multiple passed for each minimum is `0.19826 / D` of that minimum, which makes the half-width 0.19826 cycles. `--bins` is the number used for that sector's fit. The cleaned file is `--file`.

```bash
python examples/sv_cen_timing.py time --run results/sv_cen/s10_fit \
    --file results/sv_cen/cleaned/TIC324101484_s10_QLP_1800s.vot --bins 70 --half-over-d 5.300 1.621
python examples/sv_cen_timing.py time --run results/sv_cen/s11_fit \
    --file results/sv_cen/cleaned/TIC324101484_s11_QLP_1800s.vot --bins 70 --half-over-d 5.212 1.562
python examples/sv_cen_timing.py time --run results/sv_cen/s37_fit \
    --file results/sv_cen/cleaned/TIC324101484_s37_QLP_600s.vot --bins 200 --half-over-d 6.027 3.163
python examples/sv_cen_timing.py time --run results/sv_cen/s38_fit \
    --file results/sv_cen/cleaned/TIC324101484_s38_QLP_600s.vot --bins 200 --half-over-d 5.177 1.635
python examples/sv_cen_timing.py time --run results/sv_cen/s64_fit \
    --file results/sv_cen/cleaned/TIC324101484_s64_QLP_200s.vot --bins 400 --half-over-d 5.151 1.600
python examples/sv_cen_timing.py time --run results/sv_cen/s99_fit_p4 \
    --file results/sv_cen/cleaned/TIC324101484_s99_QLP_600s.vot --bins 200 --half-over-d 1.444 1.754
python examples/sv_cen_timing.py time --run results/sv_cen/s100_fit_p4 \
    --file results/sv_cen/cleaned/TIC324101484_s100_QLP_200s.vot --bins 400 --half-over-d 1.436 1.752
```

Every run reported a half-window of 0.19826 cycles for both minima. The output is `results/sv_cen/<run>/timing/`.

`nsv807_timing.py` draws each eclipse with the phase measured from the fitted time. The cycle numbers of these files are counted from the file epoch, so that phase falls half a cycle off and the template is drawn on the other minimum. Those figures were deleted. `examples/sv_cen_redraw_s64.py` redrew all of them, evaluating the template with `mikulasek_eb.phase` at the saved M0 and P0, plus one constant offset. On the sector 64 primary the residual of that drawing stays within about 20 mmag. The other sectors' figures were redrawn and not looked at one by one.

`examples/sv_cen_window_mark.py` draws the locked window on the folded curve of every sector and writes `windows_on_fold.png` in the fit folder. Both minima are shaded with the same phase width.

`examples/sv_cen_collect_dat.py` writes `results/sv_cen/timing_dat/`: 14 `.dat` files and `index.csv`. Each header names the star, the sector, the minimum, the method, the locked window, the template, the cleaning and the sigma. A rejected eclipse is a `#` line. Every other line is `jd  err_days  flag`.

| Sector | Primary timed | Secondary timed | Not timed | Median sigma, primary (s) | Median sigma, secondary (s) |
|---|---|---|---|---|---|
| 10 | 14 | 15 | 4 | 40.5 | 33.1 |
| 11 | 15 | 15 | 4 | 42.7 | 33.5 |
| 37 | 14 | 14 | 3 | 25.0 | 15.3 |
| 38 | 15 | 16 | 1 | 21.9 | 18.9 |
| 64 | 16 | 16 | 1 | 12.5 | 11.4 |
| 99 | 13 | 12 | 1 | 18.8 | 12.3 |
| 100 | 16 | 16 | 0 | 17.8 | 14.3 |

The eclipses without a time are at the ends of a sector, where one side of the window is not covered. The O-C over all sectors and the quadratic fit are not done yet.

### Timing window on NSV 807

Sectors 30 and 97. Wider windows (0.35 and 0.5) shrink the formal sigma and increase the real O-C scatter. rms(O-C/sigma) rises to 1.4 to 2.1. Window 0.15 has the smallest scatter. Window 0.25 has about honest errors. The window is reported with every timing. The NSV 807 window is not the window of the detached spot experiment. That experiment uses a half-width of 0.12 cycles.

Report: `/home/voz/projects/UPJS/nsv_807_oc/results/window_test/report.txt`

The report is written by `step7_window_test.py` in the command block above. No figure was saved.
