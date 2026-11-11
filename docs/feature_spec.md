# Light-curve features for automatic model selection (draft for review)

Status: first implementation 2026-10-30 (`src/mikulasek_eb/selection/features.py`; calibrated values in `src/mikulasek_eb_tess/rules.py`, `FeatureSettings` holds every window and threshold). Nothing here is validated; thresholds are placeholders to be set from the test stars and labelled as our heuristics. See section 8 for the first results and for where the implementation differs from this text.

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

## 9. Window and contact calibration (2026-10-30)

Problem found in section 8: windows of 1.5 FWHM were an arbitrary size and the contact indicator mislabelled the stars. Looking at the folds of the seven light curves showed what actually separates them.

Folds (see `results/features_summary.md`): only AR Aur has a flat stretch between narrow eclipses. AA And, AV Del, EK Com and XY Boo vary continuously; V477 Lyr is one eclipse on a strong reflection-like variation with no real secondary. Earlier I called AA And "detached with flat maxima"; its fold does not show that, so the first indicator was not wrong about it.

Changes (`FeatureSettings`):
- **Eclipse edge.** `edge_level` (0.05): the eclipse edge is where the eclipse has fallen to 5 percent of its depth. `edge_width1`, `edge_width2` are the full widths there, in cycles. The window of each eclipse is `window_margin` (1.0) times that width. An eclipse whose edge is not reached within half a cycle fills the cycle (a note says so).
- **Separated eclipses.** `eclipse_coverage` (the fraction of the cycle inside the windows) below `separated_coverage` (0.3). `continuous_variation` is its negation.
- **Contact candidate.** Not separated and `equal_depths` (A2 / A1 at least 0.8). Reasoning: contact components have nearly equal temperatures, so equal depths; AA And and AV Del have depths in ratio 0.26 and 0.28, EK Com and XY Boo 0.90 and 0.93. `contact_candidate` is None if there is no secondary.
- Removed `window_fwhm`, `contact_fwhm`, `contact_coverage`.

Results at `edge_level` 0.05: coverage AR Aur 0.11; AV Del 0.62; AA And 0.71; EK Com 0.76; XY Boo 0.80. Contact candidate: EK Com and XY Boo only. The ordering is the same for edge levels 0.02, 0.05, 0.1 (tested); at 0.2 AV Del drops to 0.33, close to the 0.3 threshold.

Limits (honest): one star on the separated side, so the 0.3 threshold is only placed in the gap (0.11 to 0.62) and the 0.8 depth ratio in the gap (0.28 to 0.90). Another star (a detached twin with equal depths, e.g. AR Aur itself has A2/A1 = 0.86 but is separated) is handled by the coverage test first. A semi-detached system with a flat stretch and unequal depths would be "separated, not contact", which is what we want. I have no star yet that tests the boundary. The labels are my reading of the folds, not a catalogue classification.

## 10. Local depth and trend diagnostics (2026-10-30)

### Local-continuum depth (replaces the quadrature level for depths)
The quadrature reference level fails when the out-of-eclipse variation is strong: the shallow secondary of V477 Lyr lies above that level (global depth -0.18 mag). Now:
- The **local level** of an eclipse is the mean of the brightest fold values on its two sides, between `local_inner` (0.03) and `local_outer` (0.25) cycles from its centre. `A1`, `A2` are depths below that local level; `A1_global`, `A2_global` (below the quadrature level) are kept for reference. `level_mag` is the global level.
- The **secondary** is the most prominent local flux minimum within `secondary_search` (0.2) of phase 0.5, accepted only if its prominence is at least `min_prominence_sigma` (5) median bin errors **and** it is at least `min_secondary_width_bins` (3) bins wide (a first version without the width test accepted a noise peak in a synthetic star). Otherwise `A2` is None with a note.
- For the six calibrated stars the local and global depths differ by 1 to 3 percent, and `separated_eclipses` and `contact_candidate` are unchanged. V477 Lyr: A2 = 0.072 mag (A2/A1 = 0.12), FWHM 0.020, centre 0.0035 cycles before phase 0.5; the secondary-centre error is not available (the bottom parabola is not convex). Its contact candidate is False.
- Edge-level sensitivity of the coverage (tested for 0.05 and 0.1): at 0.02 the AV Del eclipses fill the cycle (coverage 1.0, edge not reached), so the ordering AV Del below EK Com does not hold at that level.

### Trend diagnostics (replace `trend_order_hint`)
In TESS data the slow structure is mostly instrumental, so the features only describe it and make no model choice. The polynomial BIC under a white-noise assumption always chose the maximum order and is removed (`trend_max_order`, `trend_order_hint`, `trend_bic`).
- `trend_amplitude_mag`, `trend_over_depth`: largest peak-to-peak of the residual medians in time bins of `trend_bin_days` (1 d, rounded to whole cycles, at least one) per gap segment.
- `red_noise_ratio`: rms of those medians over the white-noise expectation. Values of 7 to 60 on all stars mean the residuals are far from white; it is a comparison number, not a test.
- `correlation_time_days`, `n_effective`: first lag where the autocorrelation of finely binned residuals drops below 0.5, and the span over that time. None if it never does within half the longest segment (AR Aur, with a note).
- `gap_steps`: for every gap longer than `gap_days`, the difference of the residual medians over `step_window_days` (whole cycles) after and before the gap, with its error. The unit tests recover an injected 6 mmag step.
- Observation: V477 Lyr has a step of 9.8 +- 2.0 mmag at its gap (instrumental); the other stars show steps of 0 to 2.5 mmag.
- The baseline order itself will be decided in the evidence stage with `n_effective` (section 5), not here.

## 11. Hypothesis rules and the evidence stage

Rules (`candidates.py`) turn the features into a short list of candidate
configurations, each with written reasons. They never fit anything.

- Harmonics: a contact candidate gives (1, 2, 3) and (2, 3); continuous
  variation with unequal depths gives (1, 2, 3) and (2,); separated eclipses
  give the significant harmonics only.
- Baseline order 0 and 1 are both tried. `ocon` is added when the maxima
  differ significantly (> 3 sigma).
- Exposure integration is switched on when exposure / FWHM > 0.1.
- A gap split is added when the features found gaps.
- `review_features` lists reasons to ask a human (for example a possible
  double period from odd/even differences); these are never silently resolved.

Evidence (`selection.py`): every candidate is fitted by the staged pipeline and
compared on the common set of points with a residual-based sigma,
BIC = n_eff ln(rss/n) + k ln(n_eff), where n_eff = span / correlation time.
A failed candidate stays in the table with its error. The answer is flagged
`ask_human` when the best two differ by BIC < 10, when any candidate failed,
when the review raised reasons, or when n_eff could not be measured (then the
evidence is overstated). Thresholds are placeholders until more stars exist.

## 12. Scope and calibration plan

**Scope.** Sections 1 to 11 are valid for TESS data only. The thresholds
(edge level, contact indicator, equal-depth ratio, ΔBIC, exposure switch) were
set on TESS light curves and depend on TESS cadence, noise and systematics.
The selector is not expected to work on Gaia or ground-based data without a new
calibration.

**Ground truth.** The ground truth for calibration is not a human choice of
profile and harmonics but how closely the fit follows the observations. Precise
calibration will use synthetic light curves computed with PHOEBE and sampled
with the cadence, gaps and noise of real TESS data, so that the true geometry
is known. Until then the thresholds are placeholders.

**Interim work.** Stars are analysed one by one in a supervised mode: the
existing tools are run, the chosen fits are judged by the residual numbers and
by plots, and rules are changed only when the evidence supports it.
