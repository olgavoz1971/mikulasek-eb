---
name: mikulasek-eb-tess
description: Fit and select models for TESS eclipsing-binary light curves with the `mikulasek_eb_tess` adapter of the `mikulasek_eb` package. It reads a TESS VOTable through volightcurve, uses the SAP flux, applies the TESS exposure rule and the calibrated settings, cleans outliers, and runs the automatic selection and the O-C window scan. Use this whenever the user has a TESS light curve of an EB (SPOC, QLP or full-frame, 20 s to 30 min cadence). Read the skill `mikulasek-eb` first for the method itself. Do not use the calibrated settings for other surveys.
---

# TESS adapter of mikulasek_eb

The method (model, fit, timing, tuning rules) is in the skill `mikulasek-eb`. This skill covers only what is specific to TESS. The adapter is the package `mikulasek_eb_tess`. The core `mikulasek_eb` never imports it.

Worked examples that read the files shipped in `data/` are `examples/clean_folded.py`, `examples/clean_ek_com.py` and `examples/scan_windows.py`. The two cleaning scripts draw the flagged points and do not modify the file.

## 0. Install and scope

- `pip install mikulasek-eb[tess]` pulls `volightcurve`. Without that extra the adapter does not import; the core still does.
- Calibrated for TESS cadences of 2, 10 and 30 minutes, with TESS noise and systematics. These numbers are our own heuristics, set on the test stars (`docs/feature_spec.md`). They are not validated beyond that.
- For any other survey do not use this adapter's settings. Use `mikulasek_eb.selection` with settings that the user calibrates, as described in `mikulasek-eb`, section 7b.

## 1. Use the SAP flux

A SPOC file usually carries two fluxes. **SAP** (simple aperture photometry) is the one to fit. **PDCSAP** is the mission's own detrended product, and that detrending can change the eclipse shape. This is easy to miss, because both columns are real fluxes of the same star.

`load_tess_votable(path)` loads a file that has exactly one flux column, which is what the files in `data/` are. A file with both `sap_flux` and `pdcsap_flux` is refused. Name the SAP column:

```python
load_tess_votable(path, flux="sap_flux", flux_err="sap_flux_err")
```

The names must be columns that `volightcurve` reports as flux and flux error. The function does not choose SAP on its own.

## 2. Reading the data

`load_tess_votable` reads a light-curve VOTable with `volightcurve` and returns plain arrays and metadata (time as absolute JD, magnitudes, errors, the catalogue period and epoch).
- Time system, flux columns and photometric calibration are decided by `volightcurve`. Read its skill before changing anything about them. Do not re-implement the reading.
- Convert times for the engine with `mikulasek_eb.absolute_jd_to_mjd`.
- If the file lacks what the reader needs, the reader raises. Report it. Do not patch the data.

## 3. The exposure rule

`tess_exposure_days(t)` returns the effective exposure from the median cadence: 0.8 of the cadence at every cadence except 20 s, where it is the full cadence. The 20 percent loss is the onboard cosmic-ray mitigation, a mission property for SPOC, QLP and full-frame light curves alike. An irregular cadence raises `ValueError`.
The exposure matters when exposure/FWHM is large. Without it a long cadence biases D wide. The `exptime` column of the file is not used, because VOTable has no standard for it.

## 4. Calibrated settings

`mikulasek_eb_tess.FEATURES` and `mikulasek_eb_tess.RULES` are in `src/mikulasek_eb_tess/rules.py`.
- Change one value with `dataclasses.replace(FEATURES, gap_days=0.3)`. Say in the report which value you changed and why.
- Do not edit `rules.py` to make one star work.

## 5. Selection functions

The wrappers fill in the TESS values when the settings argument is None:
- `extract_features`, `review_features`, `propose_candidates`, `select_model`.
- `tiered_search(..., features=, profiles=, start_profile=)`: baseline, then profile, then harmonics. A decisive BIC difference is `DELTA_BIC_DECISIVE` = 10. Profile 1 as the start also tries profile 1 with the curvature kept free.
- If every candidate of the start profile ends `bound-active`, `tiered_search` raises `NoValidSolution` at step 1. Report it, then rerun with `start_profile=3`, and with 4 if that also fails. The rules are in section 7b of the `mikulasek-eb` skill. The `RULES` numbers are not changed for this.
- Candidates that fail are recorded and block a confident answer. They are not replaced by a default model.

## 6. Cleaning, and the window for the timings

`flag_folded_outliers(t, mag, period_days)` compares each point with the median of its phase bin. The scatter of the bin is 1.4826 times the median absolute deviation. The defaults are 400 bins and 5 sigma, with at least 5 points in a bin. A bin with fewer points, or with no scatter, raises. The catalogue epoch is not used: the phase zero is the first time, and any zero gives the same bins. Load the file with `load_tess_votable` only. Do not import a private helper, and do not write a cleaned copy into `data/`. The original file stays as it is. Draw the fold and the time series with the flagged points marked, and write the figure under `examples/out/`. `examples/clean_folded.py` does this for `data/v477_Lyr_sector_26.vot` at the default 400 bins. `examples/clean_ek_com.py` does it for `data/EK_Com_sector_49.vot`: 400 bins raises there because phase bin 11 has 4 points, so that example uses 200 bins and 5 sigma (15 of 2651 points). With 100 bins the same file flags none. Use the stated number. Do not try a list of bin counts and keep the one that flags a convenient number.

A cut on the raw magnitude is a poor way to reject points that sit inside an eclipse: the eclipse is the faint part of the curve, so the cut takes the eclipse with the outliers. Prefer the residual from the folded bin, and look at the flagged points before accepting them. On a short or gapped light curve, pass a smaller `n_bins`.

`read_sap_background(path)` reads the SPOC `sap_bkg` column and raises if it is absent. The files in `data/` do not have it. `flag_sap_background(t_jd, t_background, background, ratio=...)` flags points whose background is more than `ratio` times the median background of those points. There is no default ratio. The V1383 Tau light curves were cleaned at 10. A cut in units of the background scatter was tried and rejected: the quiet background is so stable that any multiple of its scatter removes ordinary points. The rule does not look at the magnitude, so a point on an eclipse wall stays unless the background itself is high.

`window_scan` tries several half-widths and reports the O-C scatter of each, so the half-width is chosen from that scatter rather than from the formal error. `window` is a half-width: the points fitted are those with `|phase - E| < window`. The function lives in `mikulasek_eb.tuning`, because it has no TESS number in it, and the adapter re-exports it. `wall_band` is re-exported with it. `examples/scan_windows.py` fits one candidate of `data/EK_Com_sector_49.vot` and runs the scan. The choice of window is described in `mikulasek-eb`, section 6.

## 7. Order of work on a TESS light curve

1. SAP flux, through `load_tess_votable` (section 1).
2. Outliers, if the light curve needs them: `flag_folded_outliers`, and `flag_sap_background` when the file has `sap_bkg`. Write down what was removed.
3. Features with `extract_features`, then `tiered_search` or `select_model`.
4. Report the chosen name, its status, rms, the gains of the steps and whether the search stopped early. Time the minima only from a fit whose status is `ok`, then use `window_scan`.

The fits the tests run on the files in `data/` are examples: `fit_ek_com.py`, `fit_aa_and.py`, `fit_xy_boo.py`, `fit_v477_lyr.py`, `fit_av_del.py`, `fit_ar_aur.py`. Each writes `fit.png` and `trend.png` under `examples/out/`. `time_ek_com.py` times EK Com and writes the fit, the trend, the O-C, the folded half-width and one figure per timed eclipse. `fit_v477_lyr.py`, `fit_av_del.py` and `fit_ar_aur.py` also write the Table 1 comparison figure.

## 8. TESS specifics of the core rules

- Baselines: one polynomial per orbit, split at gaps longer than about 0.5 d (`FEATURES.gap_days`). Order 0 or 1 for an already detrended light curve.
- A 10 to 30 minute cadence on a narrow eclipse gives few points per FWHM. Use the folded curve for the template and read `mikulasek-eb`, section 5.
- The 20 s cadence has no cosmic-ray mitigation, so the exposure is the whole cadence.

## 9. Do not

- use `FEATURES` or `RULES` for a non-TESS light curve;
- fit PDCSAP when SAP is in the file, without saying why;
- lower a threshold in `rules.py` to rescue one star;
- time or report a fit whose status is not `ok`;
- replace a failed candidate by a default model;
- read the VOTable with anything but `load_tess_votable` and `volightcurve`;
- write a cleaned light curve into `data/`, or modify the original file;
- import a private name (a function whose name starts with `_`) to build that file.
