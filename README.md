# mikulasek-eb

Phenomenological light-curve model and eclipse timing for eclipsing binaries, implementing the method of

> Mikulášek, Z. 2015, A&A 584, A8, [doi:10.1051/0004-6361/201425244](https://doi.org/10.1051/0004-6361/201425244)

**This is an independent re-implementation written from the published description. It is not the author's code and is not endorsed by him.** Differences and open questions against the paper's Table 1 are recorded in `docs/mikulasek_implementation_notes.md` (section 7).

## What it does

- Global fit of the phase curve (Eqs. 13, 15, 17) with per-segment baselines and propagated uncertainties.
- Times of individual minima ("virtual minima", Eq. 9) with an explicit record of every eclipse that could not be timed.
- A staged fitting pipeline for TESS light curves, plus optional extensions that are not in the paper (second profile component, eclipse weighting), off by default.
- Fail-fast behaviour: non-converged or degenerate fits raise or are flagged, nothing is silently replaced.

Circular orbits and single-band magnitudes only.

## Layout

```
src/mikulasek_eb/   the package (model, fit, timing, tuning, extensions, tom, pipeline, io, plotting)
working_scripts/    per-star fit scripts and the Table 1 comparison
examples/           synthetic workflow and timing comparison
tests/              synthetic, fail-fast, regression on real TESS data, paper comparison
docs/  references/  implementation notes, star manifest (Table 1 values)
data/  results/     TESS light curves (VOTable) and outputs
```

## Install and test

```bash
pip install -e ".[plot,tess,test]"   # tess extra needs the volightcurve package
pytest
```

Without installing, `PYTHONPATH=src` works for the scripts; pytest uses `pythonpath` from `pyproject.toml`.

## Status

Early development. See `todo.md` for the plan and progress log.
