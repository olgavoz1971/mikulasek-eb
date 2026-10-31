# Synthetic calibration plan (PHOEBE, TESS-like)

Status: agreed in discussion, 2026-10-31. No code yet.
Scope: the selector and the timing are calibrated for TESS data only.

## 1. The synthetic set

- About 13 eclipsing binaries computed with PHOEBE. For each, a folded
  high-resolution "pure" light curve is stored.
- The pure curve is turned into unfolded TESS-like light curves covering 27 days:
  - cadence 120, 600 and 1800 s;
  - the exposure is 0.8 times the cadence (the model is averaged over it, not
    only sampled);
  - one gap of about 1 day in the middle;
  - noise (white, plus one red-noise variant, plus one variant with a wrong
    catalogue sigma, as seen in QLP files);
  - instrument-like trends (ramps at the orbit ends, a step across the gap, a
    slow drift).
- All PHOEBE parameters used for generation are written into the VOTable
  description. The spoiling recipe (seed, cadence, gap, noise, trends) is
  recorded in the same place, so every file is reproducible.
- Coverage wanted: detached, semi-detached, contact; a flat-bottomed total
  eclipse; a very shallow secondary; a spotted system (O'Connell effect); at
  least one eccentric system; periods from about 0.3 d to a few days; one case
  with a slightly wrong given period.
- Use a fractional design, not the full product. Several noise seeds for a few
  systems.
- Tune on about 8 systems, hold back about 5 for checking.

## 2. Truth

The truth for model selection is not "the right profile" (a physical model has
no true Mikulasek profile). It is the **pure light curve**.

- Fit quality: rms of the fitted template against the pure curve (mmag).
- Regret of a choice: its rms against the pure curve minus the best rms among
  all candidates fitted to the same light curve.

### Times of minima

PHOEBE `t0_supconj` is the time of superior conjunction of the primary (the
primary eclipse). `t0@system` is only the reference for the period change and
omega. PHOEBE does not declare a time scale: the description of each file must
say that its times are in the same system as our data (MJD, TCB, barycentric)
or whatever was actually used.

Three times are stored for each eclipse of the unfolded curve:

| Name | Meaning | Use |
|---|---|---|
| `T_conj` | geometric conjunction from the PHOEBE parameters (circular: `t0_supconj + E P`, secondary half a period later; eccentric: from e and omega) | physical orbit |
| `T_min_pure` | minimum of the noise-free high-resolution curve, found by a method independent of ours (for example a bisector at several depths) | best any photometric method can do |
| `T_fit` | our virtual minimum | what is tested |

- **Method error** = `T_fit - T_min_pure`, after removing the mean. It tests
  noise, cadence, gap and trends. The mean offset is reported separately,
  because the epoch is a convention.
- **Astrophysical offset** = `T_conj - T_min_pure`. This is the effect of spots
  or eccentricity, not a method error. For a circular, unspotted, symmetric
  system it must be near zero; if not, the chain has a mistake.
- Timings are compared for both circular and eccentric systems.
- Honesty of errors: check that `(T_fit - T_min_pure) / sigma` has unit scatter.

### O-C is zero by construction (first stage)

The curves are unfolded from one folded model, so every eclipse is identical and
the true O-C is exactly zero. This tests scatter and bias, **not** the ability
to follow a changing period. Injecting known time shifts (a sinusoid, a step, a
quadratic drift) is a possible second stage.

## 3. What the set is used for

1. **Features:** compare the measured depths, widths, A2/A1, contact or
   separated, maxima difference with the known values; look for bias against
   cadence, noise and trends.
2. **Rules:** fit all candidates and see which features predict the best one;
   change thresholds only where the evidence agrees across systems.
3. **Evidence stage:** does BIC choose the candidate with the lowest regret?
   Calibrate the BIC difference that triggers "ask human"; count false alarms
   and missed alarms.
4. **Timings:** method error, error honesty, and the astrophysical offset.
5. **Ablation:** change one spoiling factor at a time to see which one breaks
   which step.

The real stars (AR Aur, EK Com, NSV 807, ...) remain the final check, because
PHOEBE has no real instrument.

## 4. Other Table 2 profiles

Order driven by evidence:

1. Offline: fit all seven profiles to the pure and the degraded curves with the
   same harmonics and baseline; compare regret.
2. Add to the candidate rules only the profiles that win for some class of
   star, each with a reason in terms of features.
3. Search in tiers: choose harmonics and baseline with profile 1, then vary the
   profile with those fixed.
4. Engineering needed: the staged pipeline supports profiles 1 to 4; the start
   values for the width are derived for profile 1 only and need a version for
   each other profile; profiles 5 to 7 have no Gamma.
5. The residual-driven loop (propose a different profile after looking at the
   residuals) comes after that.

## 5. Order of work

1. Load the synthetic files and build a scoring harness (truth metrics above).
2. Validate the features.
3. Run the seven-profile experiment.
4. Tune thresholds on the tuning set; check on the held-back set.
5. Validate on the real stars.
