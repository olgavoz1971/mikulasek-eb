# Runaway primary on the 1800 s spot curves

The staged fit of the 1800 s light curves at spot longitude 0° and 180° returned `success` with an rms near 24 mmag, while the adopted primary was not an eclipse. Longitude 0° finished at amplitude about 38000 mag and Gamma about 0.000002. Longitude 180° finished at about 5400 mag and Gamma about 0.000014. Drawn on a fine phase grid, that profile is a needle of hundreds to thousands of magnitudes. The samples sit 1800 s apart, so the rms on the data points stays modest. The collapse happens when the curvature is fixed at zero and the step is refitted: amplitude runs up and Gamma runs down, and the Jacobian is rank-deficient.

This note describes an experiment. The rule is not in the package.

## The rule

After every call of `fit_global` in the staged pipeline, the primary amplitude A1 and Gamma G1 are inspected. The step is a runaway when either of these is true:

- A1 is above 5 mag
- G1 is below 0.02

Either cut is enough. The numbers are the cuts used in this experiment.

A flagged step is solved a second time. Gamma is taken out of the free list and kept at the value that step started from. The other free parameters of that step are fitted again. The first solution is discarded. The next stage starts from the held result, and Gamma may be freed again there. The same check is applied to that later stage.

If Gamma was not free in the step, the step is left as it is.

## What was run

The subjects are the two failed files only:

- `detached_none_nooconnell_SP_L000_SP_L000_1800_c1800.vot`
- `detached_none_nooconnell_SP_L180_SP_L180_1800_c1800.vot`

The script is `experiments/runaway_1800/run.py`. It fits each file twice: once with the current solver, once with the wrapper above. The model structure is the one used for the degraded spot run: profile 1, harmonic H2, O'Connell term on, straight baseline, and the exposure and gaps of that file.

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg /home/voz/projects/UPJS/skvo_veb_2/.venv/bin/python experiments/runaway_1800/run.py
```

The table is `experiments/runaway_1800/out/comparison.csv`. Pictures are `out/*_current.png` and `out/*_held_gamma.png`. The peak is the largest value of the template on a fine grid across the primary, which is where the needle appears.

| Curve | Mode | A1 (mag) | Gamma | rms (mmag) | Peak (mag) |
| --- | --- | ---: | ---: | ---: | ---: |
| 0° | current | 38015 | 0.000002 | 24.1 | 1811 |
| 0° | Gamma held | 0.75 | 0.31 | 24.1 | 0.55 |
| 180° | current | 5382 | 0.000014 | 23.7 | 257 |
| 180° | Gamma held | 0.50 | 0.80 | 11.6 | 0.53 |

Holding Gamma removes the needle on both curves. At 180° the second component is then accepted and the rms falls to 11.6 mmag. At 0° the profile is finite and the rms on the samples stays 24 mmag: the out-of-eclipse spot wave is still there, and the second component was not kept.

The cuts A1 > 5 mag and Gamma < 0.02 were only this first trial. Gamma = 1/2 is a normal eclipse in Figure 1 of the paper, so a floor on Gamma is not a rule. The script no longer uses those cuts.

## The check used now

The sample rms cannot see a needle, a collapsed width, a wild curvature term, a second-component spike, or a high harmonic sitting in a gap. After every `fit_global` the template is evaluated on a fine phase grid covering one cycle and at the phases of the data in that step. The step is rejected when the range on the grid is more than three times the range at the samples. The factor of three is the cut in this experiment: a real eclipse has nearly the same range in both, and the needles above are hundreds of times larger on the grid. A rank-deficient Jacobian is recorded with the step. It is not, by itself, a rejection: a weak parameter can be degenerate while the curve at the samples and on the grid still agree.

A rejected step is not solved again with one parameter held. The parameters of the previous step that passed are kept, and the log names the grid range, the sample range, and the rank. The next stage starts from that kept solution and is checked in the same way. If the first step fails, there is nothing to keep: the solution stays, and the run is marked unresolved. The picture is still written. This experiment does not time the curves.

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg /home/voz/projects/UPJS/skvo_veb_2/.venv/bin/python experiments/runaway_1800/run.py
```

Pictures for this check are `out/*_current.png` and `out/*_grid_check.png`. The table `out/comparison.csv` has one row per file and mode. This command was run on 2026-11-03 and exited 0.

| Curve | Mode | A1 (mag) | Gamma | Grid (mag) | Samples (mag) | Steps reverted | Flagged |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| 0° | current | 38015 | 0.000002 | 1811 | 0.56 | 0 | yes |
| 0° | grid check | 2.57 | 0.044 | 0.64 | 0.57 | 1 | no |
| 180° | current | 5382 | 0.000014 | 257 | 0.52 | 0 | yes |
| 180° | grid check | 16.9 | 0.0045 | 1.28 | 0.54 | 1 | no |

At 0° the rejected step is the one that fixes the curvature at zero. The kept solution has a grid range and a sample range that agree, and the rms stays 24 mmag. At 180° that same step was itself already narrow (amplitude 17 mag, Gamma 0.0045). Its grid range is 2.4 times the sample range, which is under the factor of three, so it was kept. The following step, the 257 mag needle, was rejected. The factor of three therefore stops the catastrophe and still lets a narrow primary through.

## The grid check is dropped

The cost of the grid is negligible beside the least-squares fit. The check is dropped because the pictures it accepts are still failures. It compares the template with itself. A sharp eclipse that the 1800 s samples happen to land on has nearly the same range on the grid and at the samples, so it is kept. Both `grid_check` pictures show that. At 180° the folded primary still runs to about 1.2 mag while the points reach about half a magnitude, and a high-frequency wave of small amplitude remains inside the minimum. At 0° the ranges agree (0.64 mag against 0.57 mag) and the residuals are still a wave of about 70 mmag.

That wave is not a harmonic. Gamma of the kept step is 0.044 at 0° and 0.0045 at 180°, so the eclipse function is still a thin spike, only no longer thousands of magnitudes tall. The exposure average is 21 samples across 1440 s, and the spike is narrower than the gap between those samples. The average is then a comb of small bumps. The time-series panel, drawn only at the data times, hides the comb.

## Bounds tried on the same two curves

The grid check stays above as the record of that trial. The script now does something else. After every solver step, a shape parameter outside the interval below is held at the near edge and the step is fitted again. The intervals are wide on purpose. Gamma from 0.3 includes the cusp in Figure 1 (the drawn value is 1/2) and excludes the needle. Harmonics are not bounded. This is not in the package. The same numbers are the dictionary `SHAPE_BOUNDS` in `experiments/runaway_1800/run.py`.

| Parameter | What it is | Interval |
|---|---|---|
| A1, A2 | Depth of the primary, of the secondary | 0.0001 to 3 mag |
| A1b, A2b | Amplitude of a second component on that minimum | −3 to 3 mag |
| D1, D2, D1b, D2b | Width | 0.005 to 0.25 of the orbit |
| G1, G2, G1b, G2b | Shape, Gamma | 0.3 to 6 |
| C1, C2 | Curvature | −2 to 2 |
| H2, S1, epoch, period | Harmonic, O'Connell term, ephemeris | not bounded |

```bash
cd /home/voz/projects/UPJS/mikulasek_eb
MPLBACKEND=Agg /home/voz/projects/UPJS/skvo_veb_2/.venv/bin/python experiments/runaway_1800/run.py
```

Pictures are `out/*_current.png` and `out/*_bounded.png`, drawn with the larger type and with both minima in the zoom panel. The adopted template is the red line in every panel. The grey line is the single-component template, and it is drawn only when a second component was kept. One component has no grey line. In the zoom, blue circles are the primary and orange squares are the secondary, both measured from their own centre, and a red or grey line is the same function in both places.

The run on 2026-11-03 finished with exit code 0.

| File | Mode | A1 | Gamma1 | D1 | rms (mmag) | Outside at the end |
|---|---|---|---|---|---|---|
| L000 | current | 38015 | 2.4e-6 | 0.055 | 24.08 | A1, Gamma1 |
| L000 | bounded | 0.099 | 2.21 | 0.250 | 11.58 | D1 on the upper edge |
| L180 | current | 5382 | 1.4e-5 | 0.070 | 23.67 | A1, Gamma1 |
| L180 | bounded | 0.501 | 0.795 | 0.016 | 11.56 | none |

The needle is gone in both bounded pictures, and both kept the second component. L180 is an interior solution: the primary is a normal eclipse and the residuals are a cloud of about ±40 mmag. L000 is not interior. Its first primary component sits on D1 = 0.25, a quarter of the orbit, with A1 only 0.1 mag, and the eclipse itself is carried by the second component. The upper edge of the width interval was doing the work there. These intervals are still the trial values, not a package rule.

## Postponed: a harmonic of high order and large amplitude

Not the comb above. This is a free cosine of order k, amplitude H_k, that grows large and oscillates between the observations.

What is already known. Section 3.2 of the paper says a Fourier series of high degree fills phase gaps with unreal artefacts: the samples can look fitted while the curve between them cannot. The paper's practice is a few low-order terms, not a bound and not a test. In this package the selector offers only harmonics it has measured as significant. On the spot-free 120 s curve that set was H2 alone. H2 and the O'Connell term stayed of order 0.01 mag in the runaway fits. No high order was free, so these two curves do not show the problem.

What is not known, and is not to be invented here. A rule that would refuse an H_k the sampling cannot resolve (fewer than two samples per cycle of that cosine) has not been tried. A numerical ceiling on |H_k| has not been tried, and no value for it is claimed. A harmonic that the samples do resolve, and that is still too large because a spot or a gap is pushing it, is an open case. The same is true of a high-order baseline polynomial on a short segment. Leave this until it appears in a fit, and do not treat the exposure comb in the pictures above as evidence for or against a harmonic rule.

## Implemented in the package: a verdict on every fit (2026-11-03)

The review of this note (`Runaway_solutions_in_the_Mikulasek-model_fit_diagnosis_and_bounding_advice.md`) agreed on the order: first make sure rubbish is never reported as a result, then try to repair the fit. The first part is in `src/mikulasek_eb/validity.py`. Every `fit` now stores `status`, `status_reasons` and `status_warnings` in `solution.stats`, and `time` raises `UnreliableSolution` for any status but `ok`, or for a missing status (a solution made before this check must be refitted). The fit itself still runs to the end, so the picture can be drawn and the reason read.

| Status | Meaning |
|---|---|
| `ok` | No check failed. This is numerically sound, not a good fit: a warning says so when rms/noise is above 1.5. |
| `bound-active` | A profile copy is outside the range of a real eclipse. |
| `unresolved` | The checks could not confirm the template. |

The limits are design choices of this package. The paper sets no bounds.

| Check | Limit | Status when violated |
|---|---|---|
| Gamma of every profile copy | 0.3 to 30 | `bound-active` |
| Depth of every copy | between 0 and 1.5 times the folded depth (the second component may be negative, same size limit) | `bound-active` |
| FWHM of every copy (measured on a fine grid) | at most 0.5 cycles, at least half the exposure in phase | `bound-active` |
| Template against the data in each eclipse core, per minimum | mean residual in the core below 0.2 of the folded depth, or within 5 standard errors | `unresolved` |
| Exposure integral, `nsub` against `4 nsub` | difference on a fine phase grid below the noise floor | `unresolved` |
| Curvature C | outside −0.6 to 0.3 | warning only |
| Jacobian rank | deficient | warning only |

Test on all 27 spot files (`experiments/runaway_1800/status_sweep.py --cadence 120 600 1800`, table `experiments/runaway_1800/out/status_sweep.csv`):

- The two curves that ran away, 1800 s L000 and L180, are `bound-active`. Reasons: Gamma 2.4e-6 and 1.4e-5, depth 38000 mag and 5382 mag, FWHM 2.5e-5 cycles, and an exposure integral that differs between `nsub` 21 and 84 by 1810 mag and 256 mag.
- The other 25 files are `ok`, including the spot-free curves at all three cadences and EK Com. Their timings are unchanged.
- The four shoulder longitudes (45°, 135°, 225°, 315°) are `ok` at 120 s and 600 s with rms about 22 mmag against a noise of about 11.5 mmag. They carry the rms/noise warning. They are numerically sound and poorly fitted, and the status does not say they are good.

Tests: `tests/test_validity.py` (12 fast tests; 3 slow tests on the generated curves, marker `slow`). The full suite passes (157 tests before these were added).

Not done yet: the bounded `trf` solver, the linear inner solve for A, H, S and the baselines, and a rule for harmonics. A fit that is `bound-active` is reported as such and not repaired.

## H1 arm at 120 s (2026-11-03)

`examples/spot_degraded_timing.py --cadence 120 --window 0.12 --harmonics 1 2 --arm h1h2` frees H1 as well. The structure then differs from the spot-free choice, so the results are in `results/synthetic_scoring/spot_degraded/h1h2/`. rms falls at 45° from 22.7 to 15.9 mmag and at 225° and 315° from about 22 to about 15.5 mmag. The primary shifts do not change (−49.3 s at 45°, +51.5 s at 315°). A wave of one cycle per orbit remains at 45°.
