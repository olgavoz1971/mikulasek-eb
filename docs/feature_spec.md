# Light-curve features for automatic model selection (draft for review)

Status: first implementation 2026-10-30 (`src/mikulasek_eb/features.py`, `FeatureSettings` holds every window and threshold). Nothing here is validated; thresholds are placeholders to be set from the test stars and labelled as our heuristics. See section 8 for the first results and for where the implementation differs from this text.

## 1. Purpose and rules

The features describe a light curve before any model is fitted. They replace the implicit judgement ("flat or round bottom", "O'Connell effect", "contact or detached") with numbers that can be logged, tested and shown to a human or an agent.

Rules:
- Input: times (MJD), magnitudes, errors, and the catalogue period and epoch. No model fit is used, only the folded curve and simple statistics.
- Every feature is a number with a unit, and an uncertainty where one can be given honestly. A feature that cannot be measured is reported as `None` with a reason, never replaced by a default.
- Features are deterministic and stored as plain JSON (they go into the evidence report of Phase 4).
- Features only propose candidate configurations. The choice between candidates is made by fitted evidence (section 5).

Common preparation (all features use it):
1. Fold with the catalogue period and epoch, then locate the primary minimum with `initial_guess` (deepest bin of a 200-bin median fold). Phase 0 is that minimum.
2. Binned fold: 200 phase bins, median per bin, error of the median from the bin scatter (at least 3 points per bin, otherwise the bin is empty).
3. Noise sigma: `noise_floor` (robust, from consecutive differences), also reported next to the median catalogue sigma.

## 2. Feature list

Units: phase in cycles, depths in mag, widths in cycles unless stated.

### A. Sampling and data quality
| Id | Feature | Definition | Used for |
|---|---|---|---|
| A1 | `n_points`, `span_days`, `n_cycles` | Count, last minus first time, span over period | sanity, drift test power |
| A2 | `cadence_s` | Median time step (`estimate_cadence_days`) | exposure choice |
| A3 | `exposure_s` | `estimate_exposure_days` (duty cycle 1) unless the user supplies it | exposure integration |
| A4 | `points_per_cycle` | Points divided by cycles | whether the shape is resolved |
| A5 | `gaps` | Number and total length of gaps longer than 0.5 d | `gap_split_days`, baseline segments |
| A6 | `sigma_noise_mag`, `sigma_catalogue_mag`, `noise_over_catalogue` | See section 1 | weights, honesty of chi-square |
| A7 | `outlier_fraction` | Fraction of points beyond 5 MAD from the fold (all phases) | clipping, bad data warning |

### B. Eclipse geometry (from the binned fold)
| Id | Feature | Definition | Used for |
|---|---|---|---|
| B1 | `A1`, `A2` | Depth of the primary and of the secondary minimum below a local reference level (median of the bins in a window 0.10 to 0.20 cycles from the minimum, on the side away from the other eclipse) | amplitude starts |
| B2 | `depth_ratio` | A2 / A1 | P/2 warning (close to 1), detachment hint |
| B3 | `secondary_offset` | Phase of the secondary minimum minus 0.5, from a parabola fitted to the bottom core (bins deeper than 0.8 of A2), with error | eccentricity flag |
| B4 | `fwhm1`, `fwhm2` | Full width at half depth (linear interpolation on the binned fold) | D starts, resolved check |
| B5 | `fwhm_ratio` | fwhm2 / fwhm1 | separate D2, G2 |
| B6 | `points_per_fwhm` | Points per cycle times fwhm1 | resolved check (Phase 3 rule: 8) |
| B7 | `exposure_over_fwhm` | A3 in days, divided by fwhm1 times period | exposure integration on or off, `nsub` |

### C. Eclipse shape (what the bottom and the wings look like)
Width at level q means the width of the region deeper than q times the eclipse depth.
| Id | Feature | Definition | Used for |
|---|---|---|---|
| C1 | `bottom_ratio` | W90 / W50 (width at 90 percent of depth over FWHM) | flat bottom (total eclipse) versus round |
| C2 | `wing_ratio` | W10 / W50 | extended wings versus compact |
| C3 | `bottom_curvature` | Quadratic coefficient of the bottom core (bins deeper than 0.9 A1), in mag per cycle squared, with error | sign tells flat, concave or pointed |
| C4 | `gamma_estimate` | Gamma of profile 1 (C = 0) that reproduces the pair (C1, C2): a lookup computed once from the profile, not a fit | Gamma starts (replaces a blind multistart where possible) |
| C5 | `flat_bottom` | True if C1 exceeds a threshold and C3 is consistent with zero | enable two-component trial early |

### D. Out-of-eclipse variation and contact character
| Id | Feature | Definition | Used for |
|---|---|---|---|
| D1 | `ooe_amplitude_over_depth` | Peak-to-peak of the binned fold outside both eclipse windows (centre plus or minus 1.5 FWHM), divided by A1; `None` if the windows cover the cycle | contact (large) versus detached (small) |
| D2 | `eclipse_coverage` | Fraction of the cycle inside the two eclipse windows | whether D1 is defined |
| D3 | `continuous_variation` | True if D2 is at least 0.9 or fwhm1 is above 0.2 cycles | contact branch of the hypothesis tree |
| D4 | `maxima_difference` | Binned magnitude near phase 0.25 minus near 0.75 (bins within 0.05 cycles), in mag, with error | O'Connell effect |
| D5 | `maxima_significance` | D4 divided by its error | `ocon` on or off |
| D6 | `harmonic_amplitudes` | Amplitudes and errors of cos and sin terms k = 1 to 4, from a linear fit to the fold with the eclipse windows masked. Defined only when D2 allows it; for contact curves the amplitudes come from the evidence stage instead | harmonics set |

### E. Period and ephemeris sanity
| Id | Feature | Definition | Used for |
|---|---|---|---|
| E1 | `half_period_test` | Fold at 2P; difference of the two primary depths in units of their error. A small difference with A2 close to A1 suggests the true period is 2P | period warning |
| E2 | `drift_shift_s` | Fold the first and the second half of the data separately with a fixed period, cross-correlate the folds, report the shift in seconds with its error | whether M0, P0 are stable; early hint of period change |
| E3 | `drift_power` | The smallest shift that E2 can detect (about the error of E2) | says whether a Pdot test is meaningful for this span |

### F. Residual structure of the fold (no model)
| Id | Feature | Definition | Used for |
|---|---|---|---|
| F1 | `fold_scatter_ratio` | Robust rms of points about the binned fold, divided by `sigma_noise_mag` | above 1: not strictly periodic or trends present |
| F2 | `trend_amplitude_mag` | Peak-to-peak of the residuals about the fold after smoothing in time (window 2 periods), per gap segment | `baseline_order`, spots |
| F3 | `trend_order_hint` | Lowest polynomial degree per segment (0 to 3) for which the smoothed residuals fall below the noise, checked by BIC | `baseline_order` |

## 3. Output

A `Features` object (frozen dataclass, JSON round trip like `Config`) with the groups above, a `notes` list for measurements that could not be made, and the package version. Same input always gives the same output, so it can sit in the cache next to the solution.

## 4. From features to candidate configurations (draft rules)

The rules only decide which candidates are tried. All thresholds are placeholders and will be written down with the data that set them.

| Decision | Rule (draft) |
|---|---|
| Branch | contact if D3, otherwise detached |
| Profile | profile 1 for both branches at first; profile 2 to 4 only when the evidence stage fits better (they have no C) |
| Gamma start | C4 if defined, else the current multistart set |
| Second component | try if C5 or if the first fit leaves rms/noise above the existing threshold |
| Harmonics | contact: (1, 2, 3) as candidate and (2, 3) as a smaller one; detached: those k with significant D6 |
| O'Connell | candidate with `ocon` on if D5 is above 3 |
| C free | candidate only if the eclipse is wide (fwhm1 above 0.1) |
| Separate D2, G2 | candidate if B5 differs from 1 by more than its error |
| Exposure | on if B7 above 0.1, `nsub` from the existing rule |
| Baseline | F3 |
| Warnings | E1 (P/2), B3 (eccentric), A4 or B6 (unresolved), A7 (bad data), E2 versus E3 (drifting ephemeris) |

## 5. Evidence stage (after the features)

Each candidate is fitted with the existing pipeline. The comparison uses BIC with a residual-based sigma (the iterative subset weighting of the paper, already implemented), because the catalogue sigma is much smaller than the model misfit for contact stars (EK Com: 1 mmag against about 6 mmag rms). Best candidate wins unless the BIC gap to the runner-up is small; then the result carries `ask_a_human = True` with both candidates listed.

## 6. Test stars and what to check

| Star | Character | What the features should show (hypothesis, to be measured) |
|---|---|---|
| AA And (S16, 2-min and 30-min) | detached, total primary, flat bottom | C5 true, D3 false, second component candidate |
| EK Com S49 | contact-like wide eclipses | D3 true, large fwhm1, C free |
| V477 Lyr, AV Del, AR Aur | other Table 1 columns | to be measured; used to set thresholds |
| XY Boo S50 (TIC 347405567, 10-min, P = 0.3706 d, 26 d) | contact, about 70 cycles, paper model uses harmonics (2, 3) | D3 true; used to test E2 and E3 for a changing period |

For XY Boo the ephemeris question is whether a single 26-day sector can detect the period change at all; E3 answers that before any Pdot fit is attempted. No expectation is stated here because I have not measured it.

Development and test split: thresholds are set on AA And, EK Com and two of the three Table 1 stars; the remaining Table 1 star and XY Boo are held out. With so few stars this is a first calibration only, not a proof that the rules generalise.

## 7. Open questions for the developer

1. Window sizes (1.5 FWHM for the eclipse windows, 0.10 to 0.20 for the local reference) are my choice. Are they acceptable?
2. Is 200 phase bins right for 10-min data (2932 points, about 70 cycles, so about 15 points per bin)? For 30-min data about 6 per bin.
3. Do you want the features also to be written into the header of the timing files?

## 8. First implementation (2026-10-30)

Differences from the draft above:
- **Reference level.** Depths are measured from the mean of the median fold at the two quadratures (0.25 and 0.75, window plus or minus 0.05), not from a window 0.10 to 0.20 cycles away. Harmonic terms therefore count in the depth (tested on a synthetic star).
- **Noise.** `sigma_noise_mag` is the robust consecutive-difference noise of the residuals about the fold, not of the raw curve, so eclipses do not inflate it.
- **Bins.** `n_bins=None` means `min(1000, max(100, n_points // 10))`. With 200 bins the narrow AR Aur eclipse had 3 bins in its core.
- **Drift (E2).** Eclipse centre = vertex of a parabola on the bottom core in each half; reported with the implied period correction. For flat-bottom eclipses the error is large (honest). Parity (E1) is the even/odd median magnitude within 0.25 FWHM of each eclipse.
- **Missing values.** Unmeasurable features are `None` plus a note in `Features.notes`. An irregular cadence gives `None` for cadence and exposure instead of an exception.
- **Timing headers.** `write_timing_ascii(features=...)` writes `# features:` (one-line JSON) and `# features_hash:`. The EK Com script does this and also saves `results/ek_com_s49_features.json`.
- **Not implemented yet:** the lookup `gamma_estimate` assumes profile 1 with C = 0 (the match distance is reported); the hypothesis rules (section 4) and the evidence stage (section 5).

First results (`examples/features_table.py`, `results/features_summary.md`), to be read as observations, not as validation:
- The contact indicator (`eclipse coverage` above 0.9 or FWHM above 0.2) says contact for AA And, whose eclipses are wide (FWHM 0.14, 1.5 FWHM windows cover 93 percent of the cycle) but which has real flat maxima (out-of-eclipse amplitude 3 percent of the depth). The windows are too generous or the indicator should use `ooe_amplitude_over_depth`. To be settled on the data.
- `trend_order_hint` is 3 (the maximum) for every star, so it carries no information yet; the residuals about the fold contain structure that the BIC always wants to absorb.
- `gamma_estimate` for the two-component AA And profile does not match profile 1 (distance 0.26), which is the expected warning.
- XY Boo S50: both eclipses are measured; primary centre moves by -51 +- 12 s between the two halves of the sector while the secondary moves by +9 +- 12 s. Opposite behaviour of the two minima is the signature of spots, not of a period change; a period change would move both the same way. This is an observation, not a conclusion. The 10-min data can detect shifts of about 10 s at one sigma.
- Questions 1 and 2 of section 7 stay open and are now `FeatureSettings` fields.
