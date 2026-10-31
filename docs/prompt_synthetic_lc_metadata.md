# Prompt for the agent generating the synthetic light curves (EB_synthetic_lc)

You generate synthetic eclipsing-binary light curves with PHOEBE and spoil them
to look like TESS data. The files are used to calibrate and test an automatic
model-selection and eclipse-timing package (`mikulasek_eb`). Read this whole
text before changing anything.

Reviewed example: `results/contact/contact_overcontact_nooconnell_G5_S05.vot`.
It already loads through `volightcurve`, its time system is declared (TDB,
barycentric, origin 2458325.0), and its zero point has a unit. Keep all of that.

Use British English in all text you write. Do not use the long dash character
in file names or inside the files.

## What is wrong or missing

1. **Everything important is in one free-text description.** Cadence, exposure,
   seed, noise level, trend parameters and the physics labels can only be read
   by parsing prose. Astropy also wraps that text at about 72 characters and
   splits words, so values get broken (for example `fillout_factor@contac` on
   one line and `t_envelope@envelope@component=0.4` on the next).
2. **The truth is not in the file** (true times of minima, noise-free flux,
   trend).
3. **No machine-readable class or flags** (contact, detached, semi-detached;
   O'Connell effect; spots).
4. **The gap and the trends are not described exactly** (no gap start and end,
   no functional form of ramps, step and drift).
5. **The true exposure is not stated.** Only `cadence_s=600` and
   `smear_frac=0.8` are in the prose.
6. **The file name does not contain the cadence.**

## What to do

### 1. Write the metadata as named PARAM elements

Keep the human-readable description, but also write every item below as a
`PARAM` (name, datatype, unit, value) in the table. Do not require anybody to
parse the description.

Spoiling record (numbers are in the units shown):

| PARAM name | Unit | Meaning |
|---|---|---|
| `product_kind` | none | `tess_observe` for the spoiled curve, `highres_one_orbit` for the pure curve |
| `recipe_id` | none | for example `S05` |
| `physics_id` | none | for example `G5` |
| `source_file` | none | name of the pure-curve file |
| `cadence_s` | s | median time step |
| `exposure_s` | s | true exposure, equal to `smear_frac` times the cadence |
| `smear_frac` | none | 0.8 for 120, 600 and 1800 s |
| `duration_days` | d | 27 |
| `gap_start_time`, `gap_end_time` | d | in the time column origin (same as the time column) |
| `noise_type` | none | `white` or `red` |
| `rel_sigma` | none | relative noise per point |
| `sigma_factor` | none | factor between the error column and the noise actually added (1.0 means honest errors) |
| `seed` | none | random seed |
| `P_true_days` | d | period used to generate the data |
| `P_meta_days` | d | period written in the file (`period`) |
| `period_eps` | none | relative difference between them |
| `trend_types` | none | for example `ramps,gap_step,slow_drift` |
| `trend_amp_frac` | none | peak-to-peak trend divided by eclipse depth |
| `eclipse_depth_mag` | mag | depth used for that scaling |
| `tess_spacecraft_orbit_days` | d | 13.7 |
| `ramp_days` | d | ramp length |
| `gap_step_amplitude` | same as flux, or relative | size of the step across the gap |
| `drift_form` | none | for example `linear` or `quadratic`, with its coefficients as further PARAMs |

If a recipe has no trend or no gap, write the PARAM with the value `none` or 0,
do not leave it out. A missing PARAM must always mean "not recorded", never
"not applied".

Physics labels for the selector (new):

| PARAM name | Meaning |
|---|---|
| `eb_class` | one of `detached`, `semi_detached`, `contact` |
| `fillout_factor` | for contact systems, from PHOEBE; empty for others |
| `has_oconnell` | true or false (is the O'Connell effect, from spots or otherwise, present) |
| `has_spots` | true or false |
| `eccentric` | true or false |
| `total_eclipse` | true or false (flat bottom in the pure curve) |

Also keep every PHOEBE parameter as a separate `PARAM`, with its full name as
in PHOEBE (`requiv@primary@star@component`, and so on), PHOEBE version and
passband. The file already lists them in the description; make them `PARAM`s
too, one per element, so nothing can be split by line wrapping.

State clearly, in a `PARAM` and in the description, the time convention: times
are in the same system as the `TIMESYS` of the file (scale, reference position,
origin). PHOEBE itself does not declare a time scale, so you must.

### 2. Store the truth

Add truth columns to the table of the spoiled curve. **Do not give them a UCD
that contains `phot.flux`**, because `volightcurve` discovers columns by UCD and
a second flux column would make the loader refuse the file. Use for example
`meta.modelled` or `meta.ref`, and give each column a name that starts with
`truth_`.

- `truth_flux_pure_smeared`: the pure PHOEBE flux averaged over the exposure at
  the observed time, with no noise and no trend. This lets the package score a
  fit directly at the data times.
- `truth_trend`: the trend alone (multiplicative or additive, say which in a
  `PARAM`), so that the recovery of the baseline can be scored.

Times of minima, for every eclipse in the unfolded light curve. Put them in a
second table of the same VOTable (one row per eclipse), not in the photometry
table. Columns:

| Column | Meaning |
|---|---|
| `E` | cycle number (integer) |
| `kind` | `primary` or `secondary` |
| `T_conj` | geometric conjunction from the PHOEBE parameters (for a circular orbit: `t0_supconj` plus E times P for the primary, half a period later for the secondary; for an eccentric orbit use e and omega through Kepler's equation) |
| `T_min_pure` | time of the minimum of the noise-free high-resolution curve, found by a method independent of ours (for example a bisector at several depths of the high-resolution eclipse) |
| `O_C_injected` | injected time shift. For now it is 0.0 for every row |

PHOEBE's `t0_supconj` is the time of superior conjunction of the primary, which
is the primary eclipse. `t0@system` is only the reference for the period change
and omega. For circular, unspotted, symmetric systems `T_conj` and `T_min_pure`
must agree to a fraction of a second; report the difference in a `PARAM`
(`conj_minus_min_pure_s`) and flag any system where it is large.

### 3. Names of the files

Include the cadence in the file name, for example
`contact_overcontact_nooconnell_G5_S05_c600.vot`. Keep the same stem for the
pure curve.

## Things to check before you hand the files back

Run the checker on every file (it reads the pure curve named by `source_file`
from the same folder). It must report no FAIL, and every WARN must be fixed or
explained to the developer:

```bash
PYTHONPATH=<mikulasek_eb>/src python <mikulasek_eb>/working_scripts/check_synthetic_lc.py file1.vot file2.vot
```

The checks it makes:

- Every file loads with `volightcurve.read_lightcurve` with no error and with an
  empty `conversion_problems`.
- The photometry table still has exactly one flux column and one flux-error
  column found by UCD.
- `cadence_s` equals the median time step in the data, and `exposure_s` equals
  `smear_frac` times `cadence_s`.
- The gap in the data is where `gap_start_time` and `gap_end_time` say.
- `P_true_days` and the `period` PARAM agree up to `period_eps`.
- There are no truncated or wrapped values in any `PARAM`.

## Do not

- Do not leave any value only in the free text.
- Do not change the time system, the zero point or the units of the existing
  columns.
- Do not invent values. If something is not known, say so in the file and tell
  the developer.

Report to the developer what you changed and anything that you could not
provide.
