# mikulasek-eb

Phenomenological light-curve model and eclipse timing for eclipsing binaries, implementing the method of

> Mikulášek, Z. 2015, A&A 584, A8, [doi:10.1051/0004-6361/201425244](https://doi.org/10.1051/0004-6361/201425244)

**This is an independent re-implementation written from the published description. It is not the author's code and is not endorsed by him.** Differences and open questions against the paper's Table 1 are recorded in `docs/mikulasek_implementation_notes.md` (section 7).

## For AI agents

If you were asked to install this package from GitHub, do these steps in order. Stop before you write any fit.

1. Clone `https://github.com/olgavoz1971/mikulasek-eb` into the working directory and leave the clone there. The clone is the documentation: `AGENTS.md`, `llms.txt`, `skills/`, `examples/`, `docs/` and `data/` stay on disk for the rest of the task.
2. Read `AGENTS.md`, then `skills/mikulasek-eb/SKILL.md`. When the light curve is from TESS, also read `skills/mikulasek-eb-tess/SKILL.md`.
3. From that clone run `pip install ".[plot,tess]"`.

Do not run `pip install git+https://github.com/olgavoz1971/mikulasek-eb`. That command keeps `src/` and discards the documentation.

## What it does

- Global fit of the phase curve (Eqs. 13, 15, 17) with per-segment baselines and propagated uncertainties.
- Times of individual minima ("virtual minima", Eq. 9) with an explicit record of every eclipse that could not be timed.
- A staged fitting pipeline for TESS light curves, plus optional extensions that are not in the paper (second profile component, eclipse weighting), off by default.
- Fail-fast behaviour: non-converged or degenerate fits raise or are flagged, nothing is silently replaced. A finished fit also carries a status (`ok`, `bound-active`, `unresolved`); only `ok` may be timed. What the statuses mean and how the solver and the limits came about are in `docs/runaway_fits.md`; the full table of checks is in `docs/mikulasek_implementation_notes.md`, section 11. How the model is chosen is in "How a model is chosen" below.

Circular orbits and single-band magnitudes only.

A plain-language description of the automatic script (`diagnose_fit.py`) is in
`docs/pipeline_explanation_for_dummies.md` (profiles versus candidates).

## Layout

```
src/mikulasek_eb/            core: model, fitting, timing, tuning, pipeline. Needs numpy, scipy, astropy only.
src/mikulasek_eb/selection/  survey-neutral selection machinery. No calibrated defaults: you pass the settings.
src/mikulasek_eb/compare/    model-free minimum measurement, to compare with other methods.
src/mikulasek_eb/plotting.py optional, needs matplotlib (extra "plot").
src/mikulasek_eb_tess/       TESS adapter: reader, exposure rule, calibrated settings, cleaning. Needs volightcurve (extra "tess").
AGENTS.md  llms.txt         how an agent should use this repository, and an index of the docs
examples/                    scripts that read the files in data/
skills/                      agent skills: mikulasek-eb (method, core, selection) and mikulasek-eb-tess (adapter)
working_scripts/             per-star fit scripts used by the tests, and the Table 1 comparison
tests/core/                  tests that build their own arrays and need no data files
tests/tess/                  tests that open the TESS files in data/ (regression, paper comparison)
docs/  references/           implementation notes, star manifest (Table 1 values)
data/                        TESS light curves used by the tests and the paper comparison
```

The core never imports the selection layer or the adapter. A test (`tests/core/test_core_independence.py`) fails if it does.
For another survey, write a small adapter that reads the data and builds a `FeatureSettings` and a `RuleSettings` calibrated for that survey; `mikulasek_eb_tess/rules.py` is the worked example.

`examples/` reads the light curves in `data/`. `clean_folded.py` and `scan_windows.py` show cleaning and the timing-window scan. `fit_ek_com.py`, `fit_aa_and.py` and `fit_xy_boo.py` are the fits from the regression and selection tests. `fit_v477_lyr.py`, `fit_av_del.py` and `fit_ar_aur.py` lay the published Table 1 model on the TESS file, as the paper test does. `time_ek_com.py` times the eclipses of EK Com. Figures go to `examples/out/`.

## Install and test

```bash
pip install -e ".[plot,tess,test]"   # tess extra needs the volightcurve package
pytest
```

Without installing, `PYTHONPATH=src` works for the scripts; pytest uses `pythonpath` from `pyproject.toml`.

## How a model is chosen

Every fit has to settle four questions: the shape of the eclipse (one of the profiles of Table 2 of the paper), which harmonics describe the out-of-eclipse variation, whether to add the O'Connell term, and how to remove the instrumental trend. The selection does **not** try every combination of these. That would be many full fits, and some combinations never converge.

### The ingredients

- **The instrumental trend** is a polynomial per stretch of the light curve. A new stretch starts where two points are more than half a day apart. It is fitted together with the eclipse template, in magnitudes, so it acts as a multiplicative trend. A trend that is additive in flux (background, contamination) is not removed by it.
- **Features** are measured on the folded curve: depths, flat stretches, the out-of-eclipse variation, the noise and the correlation time of the residuals. From these, rules propose a *short* list of candidates: a few harmonic sets, the baseline orders and the O'Connell term.
- **Fit status.** Every finished fit is judged `ok`, `bound-active` or `unresolved`. A candidate that is not `ok` is kept in the table but is never chosen. See "Run-away fits" below.
- **BIC** decides between candidates that are `ok`. It uses a Gaussian likelihood with a free variance, evaluated on the points that every candidate kept, and the effective number of independent points, because the residuals of space photometry are correlated. A difference of 10 is treated as decisive. When the choice is not clear, the result says that a person should look.

### Two ways of searching

`select_model` fits the proposed candidates and compares them.

`tiered_search` is a shorter named strategy on top of it. It works in three steps, and each step compares the previous winner with a few new candidates on the same points:

1. **Baseline.** The start profile, the richest harmonic set the rules propose, and each proposed baseline order.
2. **Profile.** The winning baseline and harmonics, with the other profiles of Table 2 (1 to 4). If the start profile is 1, profile 1 with the curvature kept free is also tried, for flat-bottomed minima.
3. **Harmonics.** The winning profile and baseline, with the smaller harmonic sets.

The search stops before step 3 when steps 1 and 2 each gained less than 10 in BIC. Profiles 5 to 7 are not part of it, because the staged fit frees Gamma and those profiles have none.

On one test file (a synthetic contact binary, 600 s cadence) the three steps chose the same model as fitting all 16 combinations, in 35 s instead of about 5.5 min, and two fits that do not converge were never started.

### Which profile to start from

The start profile is an argument of `tiered_search`. Profile 1 is the ordinary start. Because the baseline order is decided in step 1, it is decided under the start profile. When the start profile does not suit the star, every candidate of step 1 can end `bound-active`, and the search stops with an error before it reaches the profile step. On the real eclipsing binary SV Cen that happened on several sectors. Rerunning with another start profile (3, and on two sectors 4) gave fits with status `ok`. The package does not do this on its own: the caller says which start profile to use and reports it.

### Limits of the evidence

- BIC counts every point as independent when the correlation time could not be measured, so its gains are then overstated. The result says when that happens.
- BIC judges the whole curve. A timing depends on the eclipse shape, and a good whole-curve rms can hide a flat bottom that the template misses.

### Run-away fits

A least-squares run can finish with `success` and a template that is not an eclipse, for example a needle of enormous depth narrower than the exposure, while the rms on the samples looks fine. The package judges every finished fit, refuses to time one that is not `ok`, and refits a `bound-active` candidate once with a bounded solver during selection. What was tried, what was dropped, and what is still unknown is in `docs/runaway_fits.md`.

## Scope of the automatic selector

The selection machinery (`mikulasek_eb.selection`) has no calibrated numbers and cannot run without
explicit settings. The only calibrated set is `mikulasek_eb_tess.FEATURES` and `mikulasek_eb_tess.RULES`,
designed for **TESS light curves only** (cadences of 2, 10 and 30 minutes, TESS noise and systematics).
They are not expected to work on other surveys such as Gaia, and must not be used there without a new
calibration. The fitting engine and the timing code do not have this restriction.
Change one calibrated value with `dataclasses.replace(mikulasek_eb_tess.FEATURES, gap_days=...)`.

## Status

Early development. See `todo.md` for the plan and progress log.

## Licence

MIT, see `LICENSE`.
