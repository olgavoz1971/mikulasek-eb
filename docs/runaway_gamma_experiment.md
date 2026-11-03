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

Pictures for this check are `out/*_current.png` and `out/*_grid_check.png`. The table `out/comparison.csv` has one row per file and mode.
