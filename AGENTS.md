# Using mikulasek-eb

Instructions for an agent working in this repository. The scientific method is not repeated here.

## Read first

1. `skills/mikulasek-eb/SKILL.md` for any eclipsing-binary fit, model choice or timing.
2. `skills/mikulasek-eb-tess/SKILL.md` as well when the light curve is from TESS. Read it after the first skill. Do not use its calibrated settings on another survey.
3. `README.md`, section "How a model is chosen", for the human account of selection.
4. The docstring of the example you are about to copy.

Do not guess a function signature, a parameter name or a default. Check `references/api.md` and the docstring in `src/`. If they disagree, the source wins. If the paper is unclear, say so. Do not invent a step the code does not have.

You are in the clone of this repository. Leave it in the working directory. Install from here with `pip install ".[plot,tess]"`. Do not run `pip install git+https://github.com/olgavoz1971/mikulasek-eb`: that command keeps `src/` and discards this file, the skills, the examples and `data/`.

## Where the code lives

- `src/mikulasek_eb` is the core: model, fit, timing. It needs numpy, scipy and astropy. It must not import `volightcurve`, `lightkurve`, `lc_discovery` or `mikulasek_eb_tess`.
- `src/mikulasek_eb/selection` chooses a model. It has no calibrated numbers. Missing settings raise `ConfigError`.
- `src/mikulasek_eb/compare` measures a minimum without the model.
- `src/mikulasek_eb_tess` reads a TESS VOTable, applies the TESS exposure rule and the calibrated `FEATURES` and `RULES`, and cleans light curves. It needs the `tess` extra (`volightcurve`).

Import the everyday names from `mikulasek_eb`. Import everything else from the module that defines it (`fitting`, `timing`, `tuning`, `selection`, and so on). The list is the top of `references/api.md`.

## Examples and light curves

`examples/` holds the scripts to copy. Run them from the repository root. They read the TESS files in `data/`. Those files are the light curves for the examples and the TESS tests. Use them. Do not point an example at another tree.

`examples/fit_paper_star.py` (and the three thin scripts that call it) also imports `working_scripts/paper_forward_check.py`. That is the Table 1 comparison. Prefer `examples/` for a new script.

Figures go to `examples/out/`, which is gitignored.

## What not to treat as the interface

- `working_scripts/` is how the tests call the engine, including the Table 1 comparison. It is not a second public API.
- `references/` holds the star manifest and `api.md`. Use `api.md` to check a name. The manifest is test data.
- `todo.md` is the development log. It is not a procedure for a fit.
- `docs/mikulasek_implementation_notes.md` records differences from the paper. `docs/runaway_fits.md` records why a finished fit is marked `ok`, `bound-active` or `unresolved`. `docs/feature_spec.md` records how features are measured. Read them when the skill points there. Do not promote a note in `todo.md` over the skill.
