# Pipeline explanation for Dummies

This note describes what `examples/diagnose_fit.py` actually does today.
It is a wrapper around `select_model`. The numbered steps below are that
function, plus the plots at the end.

British English. The selector is for TESS light curves only.

## The flow in five steps

1. **Load the light curve.** The VOTable is read. We take times (as MJD),
   magnitudes, errors, the period and the epoch. Every later step uses these
   raw points. They are not binned.

2. **Extract features (no fit).** The curve is folded on the given period.
   A median curve in phase is used only to measure things: eclipse depths and
   widths, whether the bottoms are flat, whether the eclipses are separated or
   the variation is continuous (contact-like), the O'Connell effect (unequal
   maxima), the TESS exposure (0.8 times the cadence, except at 20 s), gaps
   and noise. The output is a table of numbers. It is not a model.

3. **Propose candidates.** The features become a short list of
   configurations to try. Each has a short name such as `h123_b0`. That name
   is **not** a profile number. See the next two sections.

4. **Fit every candidate and pick one.** Each candidate is fitted to the
   unfolded light curve (all points, in time). They are compared with BIC.
   The winner is kept. If the best two are close, a candidate failed, or a
   review rule fired, the result is marked "ask human". There is **no**
   second round that looks at the residuals and then tries another Table 2
   profile. That loop is not built.

5. **Show the winner.** Residuals (data minus model, after removing the
   baselines) are printed in millimagnitudes. A six-panel figure is saved.
   With `--show` the figure also opens in a matplotlib window.

Inside one fit the pipeline may still add a **second component** to the
primary eclipse if the single-component residuals are large. That is an extra
bump of the **same** profile, not a switch to another row of Table 2.

## What a profile is

The paper (Mikulasek 2015, Table 2) lists **seven alternative eclipse
functions**. Each row is a different mathematical shape for one eclipse:
how the brightness drops and recovers. That row is what we call a
**profile**.

| Profile | In plain words | Shape parameters besides depth A |
|---|---|---|
| 1 | The paper's default (Eq. 13). Round or slightly flattened bottom, with optional extra flattening C | D, Gamma, C |
| 2 | Profile 1 with C fixed at 0 | D, Gamma |
| 3 | A Gaussian raised to Gamma | D, Gamma |
| 4 | Andronov-type, finite width | D, Gamma |
| 5 | A plain Gaussian | D only |
| 6 | A triangle | D only |
| 7 | exp(1 - cosh), no Gamma | D only |

**Yes: a profile is one row of that table.**

**Yes: all seven rows are implemented** in `EBModel`. You can ask the engine
to fit any of them by setting `ModelSpec(profile=...)`.

**No: the automatic script does not use three of them.** Today
`propose_candidates` always sets **profile 1**. Profiles 2 to 7 are in the
package, but `diagnose_fit.py` / `select_model` never offer them. The staged
pipeline can fit profiles 1 to 4 (they have Gamma). Profiles 5 to 7 exist in
the model and can be fitted with `fit_global`, but the automatic selector
does not call them.

If you saw three things in the candidate list, those were almost certainly
**harmonic sets**, not profiles. Typical names are `h123`, `h23` and `h2`.

## What a candidate is

A **candidate is not** "a profile plus initial values of the parameters".

A candidate is a **full configuration of the fit**, written as a `Config`:

- **Profile:** always 1 in this script.
- **Harmonics:** extra cosine waves added on top of the two eclipses
  (`cos 2 pi k theta`). They describe the out-of-eclipse variation
  (ellipsoidal shape, reflection). They are **not** another Table 2 row.
  `h123` means orders 1, 2 and 3; `h23` means 2 and 3; `h2` means 2 only.
- **O'Connell term** (`ocon`): a sine wave that lets the two maxima differ.
  On or off for the whole list, according to the features.
- **Baseline order:** `b0` is a constant level per stretch of data; `b1` is
  a straight line per stretch. This absorbs instrumental drift. It is not
  part of the star's eclipse shape.
- **Exposure integration:** on when the exposure is a large fraction of the
  eclipse width. The model is averaged over the TESS exposure at each point.
- **Gap split:** if there are gaps, each stretch of data gets its own
  baseline.

**Initial values come later.** When a candidate is fitted, `initial_guess`
measures starting depths and widths from the fold. Those starts are not stored
in the candidate. Two candidates can share the same profile and still start
from the same guess; they differ by which extra terms they are allowed to use.

The short name encodes only the parts that vary in the list. Example:
`h23_b1` means profile 1 (silent), harmonics 2 and 3, linear baseline.

How many candidates? Usually **four**: two harmonic sets times two baseline
orders. For a clearly detached star there may be only two (one harmonic set
times two baselines). The harmonic sets are chosen from the features:

- contact-like (continuous variation, similar depths): try `(1, 2, 3)` and
  `(2, 3)`;
- continuous variation with unequal depths: try `(1, 2, 3)` and `(2,)`;
- separated eclipses: try only the harmonics that already look significant.

## What this script does **not** do yet

- It does not choose among Table 2 profiles 1 to 7.
- It does not look at the residual plot and then retry another profile.
- It does not change the period or the epoch (those come from the file).
- It is not calibrated for Gaia or ground-based data.

Those belong to later work (AI-helped mode and the PHOEBE synthetic
calibration). Until then, a human judges the winner by the residual numbers
and the plots.
