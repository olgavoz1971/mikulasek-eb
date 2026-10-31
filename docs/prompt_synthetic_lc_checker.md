# Prompt for the agent generating synthetic light curves: the acceptance checker

You generate synthetic eclipsing-binary light curves (PHOEBE, spoiled to look like
TESS data). Each file must meet the metadata requirements in
`prompt_synthetic_lc_metadata.md`. A checker verifies that. This text tells you how
to use it and what you may and may not do with it.

British English in all text you write. Do not use the long dash character in file
names or inside files.

## What the checker is

`mikulasek_eb/working_scripts/check_synthetic_lc.py`, with its logic in
`mikulasek_eb/src/mikulasek_eb/synthetic_check.py`. It belongs to the developer. It
only reads files; it never changes them.

Run it on every file you generate, before you report the file as finished:

```bash
cd <path to mikulasek_eb>
PYTHONPATH=src python working_scripts/check_synthetic_lc.py <file1.vot> [<file2.vot> ...]
```

It reads the pure high-resolution curve named by the `source_file` parameter from
the same folder, so keep the pure curve next to the spoiled file. The exit code is 1
if any check failed.

## What it checks

- The file reads through `volightcurve` with exactly one time column, one flux
  column and one flux-error column, and no calibration or time problems. The truth
  columns must not have a UCD that contains `phot.flux`.
- All required parameters are present (see the metadata prompt), the labels are
  valid, and the booleans are written as `true` or `false`.
- The declared cadence, exposure, duration, gap, period and time system agree with
  the data and with each other.
- The noise and trends: the rms of (flux minus truth minus trend) divided by the
  error equals `1 / sigma_factor`; the relative noise equals `rel_sigma`; the
  peak-to-peak trend divided by the flux depth equals `trend_amp_frac`;
  `eclipse_depth_mag` equals the depth of the smeared pure curve.
- The eclipse truth: consistent array lengths, valid kinds, conjunction times equal
  to epoch plus E times P (circular orbits), the offset between `T_min_pure` and
  `T_conj` per kind, `conj_minus_min_pure_s` covering both kinds, and per-eclipse
  flags for eclipses in the gap or outside the data.
- The smeared truth: `truth_flux_pure_smeared` must agree with an independent
  average of the pure curve over a window centred on the time stamp. A window that
  starts or ends at the stamp must fail.

## Outcomes

| Mark | Meaning | What you do |
|---|---|---|
| PASS | the check agrees | nothing |
| info | information only (for example an offset that is expected for an eccentric or spotted system) | nothing |
| skip | the check could not run (for example no pure curve) | find out why and tell the developer |
| WARN | an inconsistency or an item that is missing | fix the cause if it is a mistake; if you think it is intended, explain it in your report |
| FAIL | the file is not acceptable | fix it and run the checker again |

A file is accepted only with **no FAIL**. The developer decides about every WARN.

## Rules

1. **Do not edit the checker or its tolerances.** If you think a check is wrong,
   say so in your report and leave it as it is.
2. **Do not hide a warning.** Do not change a number only to make a warning
   disappear. Fix the quantity it describes. If you change a parameter, make sure
   the data still agree with it (the checker tests that).
3. **Report the checker output as it is**, with the counts of PASS, WARN and FAIL
   for every file, and a sentence for every WARN or skip: what it says, what you
   did, or why you left it.
4. **Run the checker on every file you generate**, including every cadence, noise
   and trend variant, not only on a sample. A good place is the acceptance step of
   your own tests.
5. **Do not invent values.** If you do not know something, say so in the file and
   tell the developer.
6. **Do not change the time system, the zero point or the units** of the existing
   columns.

## What the checker cannot tell you

It checks consistency. It does not tell you whether PHOEBE was set up correctly,
whether your method for `T_min_pure` is good, or whether the system really is what
`eb_class` says. Those remain your responsibility; describe what you did in the file
description and in your report.

## Notes for the first batch

The checker was written and tested on one file only (a circular contact system,
cadence 600 s). Its tolerances for detached, eccentric and spotted systems have not
been tried on real output, so expect a few warnings that come from the checker and
not from your files. Report them as they are. The developer will adjust the
thresholds, not you.

## What to hand back

For every file:

- the file and its pure curve;
- the checker output (verbatim);
- a short list of the WARN and skip items with your explanation;
- anything you could not provide.
