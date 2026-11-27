# API reference

Import locations. The root `mikulasek_eb` exports only the everyday surface (`fit`, `time`, `validate_lightcurve`, `fit_borrowed`, `run_staged_fit`, `FitOutput`, `TimingOutput`, `ModelSolution`, `FitPipelineResult`, `ParameterPlan`, `Config`, `ModelSpec`, `DataSetup`, `FitConfig`, `FitPolicy`, `TimingRecipe`, `EBModel`, `PROFILES`, `phase`, `mjd_to_jd`, `absolute_jd_to_mjd` and the errors). The functions below come from their own modules: `model`, `fitting`, `timing`, `tuning`, `tom`, `extensions`, `exposure`, `validity`, `borrow`, `descriptors`, `helpers`.
Selection (survey-neutral, settings are arguments): `mikulasek_eb.selection` has `extract_features`, `review_features`, `propose_candidates`, `select_model`, `tiered_search`, `make_candidate`, `FeatureSettings`, `RuleSettings`, `Candidate`, `SelectionResult`, `TieredResult`. TESS adapter with calibrated settings: `mikulasek_eb_tess`.

The signatures below are for the engine.

## EBModel(profile=1, harmonics=(2,), ocon=False)
Callable `model(theta, p) -> mag`. Parameter dict keys: `A1 A2 D1 G1 C1` (+ optional `D2 G2 C2`; absent means equal to D1,G1,C1), `H<k>` for every k in `harmonics` (term cos(2πkϑ)), `S1` if `ocon` (sin 2πϑ). All keys needed by the model must be present in the dict, free or not. `model.names(shared_shape=True)` lists them.
Profiles: 1 = A(1+Cx²){1−[1−exp(1−cosh x)]^Γ}; 2 = same without C; 3 = A{1−[1−exp(−x²/2)]^Γ}; 4 = A·Re[(1−|x|^Γ)^{3/2}] (x = φ/D). Γ<1 cusp, Γ=1 rounded, Γ>1 flatter bottom.

## Phase and times
`phase(t, M0, P0, Pdot=0)`; `predicted_time(E, M0, P0, Pdot=0) = M0 + P0 E + 0.5 P0 Pdot E²`.

## fit_global(t, y, sigma, model, p0, free, ephem0, free_ephem=("M0","P0"), seg=None, order=0, rescale=True) -> FitResult
`free`: template keys to vary; `free_ephem`: subset of M0, P0, Pdot. Per-segment polynomial baselines (degree `order`) are solved linearly. FitResult: `p`, `ephem`, `err` (dict incl. ephemeris), `chi2`, `dof`, `rchi2`, `rms` (mag, after baselines), `names`, `cov` (order = `names`), `success`. `rescale` multiplies the covariance by rchi2.

## descriptors.eclipse_descriptors(source, p=None) -> dict
Depth and full width at half depth of both fitted eclipses, as properties of the fitted curve and not of its parameters. `source` is a `ModelSolution` (then `p` is ignored) or an `EBModel` with `p`. Returns `{"eclipse_1": {"depth_mag", "fwhm_cycles"}, "eclipse_2": {...}, "depth_ratio", "fwhm_ratio"}` (ratios are eclipse 2 over eclipse 1). Each eclipse is the sum of its profile copies (main profile plus the second component, if any). Harmonics, O'Connell, the other eclipse, the baseline and the exposure integration are excluded. Depth is the deepest point (mag, fainter positive), width is the distance between the half-depth crossings on both sides of that point (linear interpolation on a 2.5e-5 cycle grid). Raises `ValueError` if `p` is missing, if an eclipse has no positive depth, or if it does not fall below half depth within 0.25 cycles. No uncertainties: a solution stores parameter errors but not their covariance.

## descriptors.eclipse_widths(source, p=None, level=0.05) -> dict
Full width of each fitted eclipse at `level` of its depth, in cycles: `{"primary": w, "secondary": w}`. Same eclipse function as `eclipse_descriptors` (main profile plus second component; no harmonics, O'Connell, other eclipse or exposure), so it needs no contacts and does not depend on the fitted harmonics. `level` must lie strictly between 0 and 1 (profiles such as profile 1 only tend to zero, so a width at 0 is not defined); 0.5 gives the FWHM. Raises `ValueError` for a bad level, a missing `p`, no positive depth, or an eclipse that does not fall below the level within 0.25 cycles.

## descriptors.window_half_widths(source, p=None, n_widths=1.5, level=0.05, max_half=0.25) -> dict
Half-windows for timing: `half = n_widths * eclipse_widths(...)[kind] / 2`, returned as `{"primary": half, "secondary": half}`. Pass each to `TimingRecipe(window=half)` for that kind. The default of 1.5 widths at 5 percent is preliminary, tested on two synthetic stars only. Nothing is capped or enlarged silently: a half-window above `max_half` raises `ValueError`. The O-C scan is a warning for the researcher, not an input to this function.

## fit_minimum(t, y, sigma, model, p, ephem, E, window=0.1, order=1, free=(), scan=0.1, rescale=True) -> MinResult
Points with |ϑ−E| < window. Free: delta [d] (+ multiplicative scales of A, D, G of that eclipse named in `free`). `scan` = initial shift grid half-width in phase. MinResult: `E, T, sigma, delta, n, rchi2, scales`.

## time_cycles(t_mjd, mag, sig, solution, min_points=4, min_each_side=2, order=1) -> CycleTimingOutput
One shift of the whole template per orbital cycle, a separate method from `time`. Points of cycle n have phase in (n − 0.25, n + 0.75). Default floor: 4 points, 2 on each side of the cycle centre. `cycles` holds E (integer), delta (days), T_primary and T_secondary (absolute JD), sigma, n, rchi2. The same shift is used for both eclipses of the cycle.

## virtual_minima(t, y, sigma, model, p, ephem, window=0.25, order=1, free=(), min_points=4, min_each_side=2) -> dict of arrays (E, T, sigma, delta, n, rchi2, walls_ok)
All eclipses with enough points on both sides. Default floor 4 points, 2 per side. `walls_ok` is True when both steep walls are sampled; uncovered walls skip the event only if `reject_uncovered_walls=True`. Call per instrument/filter. window is capped at 0.25.

## TimingRecipe defaults
`min_points=4`, `min_each_side=2`, `require_contact_walls=False`. The wall test is stored as `walls_ok` on each timed minimum of `time`.

## plotting.save_timing_windows(t_mjd, mag, sig, solution, timing, kind, out_dir, eclipses=None) -> list[Path]
Writes one picture per eclipse of `kind` timed by `time` (`E01_primary.png`, `E00_secondary.png`): whole-cycle points, fitted window, template at the fitted position, found time and prediction; rejected eclipses with their reason. `plotting.plot_timing_window(..., kind, E)` draws one. Needs matplotlib. Same arrays and solution as passed to `time`. Uses `timing.event_curve(t, y, sigma, model, p, ephem, E, delta, window, order)`, which repeats no fit.

## plotting.plot_timing_overview(t_mjd, mag, sig, solution, timings, save_path=None) -> Figure
Folded light curve over two cycles with the primary and secondary timing windows in one picture.

## plotting.save_timing_cycles(t_mjd, mag, sig, solution, timing, out_dir, cycles=None) -> list[Path]
Writes one picture per orbit timed by `time_cycles` (`E00.png`): all points used in that cycle, template at the fitted shift, primary and secondary times (the same shift). Zooms of each eclipse under that shift go in `out_dir/zoom/` as `E00_primary.png` and `E00_secondary.png`; the zoom band is 1.5 widths at 5 percent of depth. `plot_timing_cycle(..., E)` and `plot_timing_cycle_zoom(..., E, kind)` draw one. Cycle `E` is counted from T0. Needs matplotlib. Same arrays and solution as passed to `time_cycles`.

## fit_ephemeris(T, sigma, E, deg=2) -> (coefficients, residuals, weighted rms)
Weighted polynomial T(E).

## eb_tuning.py
`noise_floor(t, y)`, `D_from_fwhm(fwhm, G)`, `initial_guess(t, mag, P)`, `multistart_global(...)` -> (best, rows), `residuals(t, mag, sig, model, r, seg, order)` -> (res, theta), `binned_residuals(res, theta)`, `clip_timings(vm, deg, nsig)`, `window_scan(...)` -> (timings per window, summary rows).

## compare_timings.py
`python compare_timings.py ours.csv other.csv --P P --M0 M0 [--tol 0.02] [--deg 2] [--png out.png]`
