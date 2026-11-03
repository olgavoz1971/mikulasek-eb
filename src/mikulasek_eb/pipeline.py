"""Staged global fit: locate, multistart, shared C, separate D2/G2, clip, optional extensions."""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

from .tuning import (
    D_from_fwhm,
    binned_residuals,
    initial_guess,
    multistart_global,
    noise_floor,
    out_of_eclipse_mask,
    residuals,
)
from .model import EBModel
from .extensions import (
    _fold_phase_from_t0,
    _mjd_to_jd,
    _two_component_timing_ok,
    eclipse_reweighted_sigma,
    refine_two_component,
    two_component_degenerate,
)
from .helpers import complete_template, out_of_eclipse_phase_mask, report_fit
from .config import Config, FitConfig
from .errors import ConfigError, InvalidLightCurve, NoValidSolution
from .exposure import suggest_nsub
from .fitting import fit_global
from .validity import FitVerdict, amplitude_scale, assess_fit, shape_bounds

logger = logging.getLogger(__name__)


@dataclass
class FitPipelineResult:
    """Outputs of the staged fit pipeline.

    Attributes:
        model (EBModel): Phenomenological model instance.
        adopted: Final ``FitResult`` (single- or two-component).
        single_component: Last valid single-component ``FitResult``.
        seg (ndarray | None): Baseline segment indices (None = one segment).
        order (int): Baseline polynomial degree.
        t_mjd (ndarray): Times used in the fit (MJD).
        mag (ndarray): Magnitudes (possibly clipped).
        sig (ndarray): Uncertainties (catalogue).
        noise_floor (float): Out-of-eclipse noise (mag).
        initial_guess (dict): Folded start (M0, A1, A2, fwhm).
        config (Config): The full layered configuration used (provenance).
        verdict (FitVerdict | None): Status of the adopted template (``ok``,
            ``bound-active`` or ``unresolved``). Always set by ``run_staged_fit``.
    """

    model: EBModel
    adopted: object
    single_component: object
    seg: np.ndarray | None
    order: int
    t_mjd: np.ndarray
    mag: np.ndarray
    sig: np.ndarray
    noise_floor: float
    initial_guess: dict
    config: Config
    verdict: FitVerdict | None = None


def _shape_free_names(r, model: EBModel) -> list[str]:
    """Eclipse-shape keys, excluding frozen harmonics and ephemeris.

    Args:
        r: Fit result whose ``names`` list is the last free set.
        model (EBModel): Template (harmonics / O'Connell).

    Returns:
        list[str]: Names to free in the weighted shape refit.
    """
    frozen = set(_harmonic_free_names(model)) | {"S1", "M0", "P0", "Pdot"}
    return [n for n in r.names if n not in frozen]


def _harmonic_free_names(model: EBModel) -> list[str]:
    """Return H* and S1 names to free for this model.

    Args:
        model (EBModel): Configured model.

    Returns:
        list[str]: Parameter names.
    """
    names = [f"H{k}" for k in model.harmonics]
    if model.ocon:
        names.append("S1")
    return names


def run_staged_fit(
    t_mjd: np.ndarray,
    mag: np.ndarray,
    sig: np.ndarray,
    period_days: float,
    epoch_jd: float,
    config: Config | FitConfig,
    progress=None,
) -> FitPipelineResult:
    """Run staged global fit and optional two-component primary extension.

    Catalogue sigma only. Adoption of two-component follows ``aa_and_timing_oc``
    shape checks 1-3 and virtual-minima check 4 when ``config.allow_two_component``.

    Args:
        t_mjd (ndarray): Observation times in MJD.
        mag (ndarray): Magnitudes.
        sig (ndarray): Catalogue magnitude uncertainties.
        period_days (float): VOTable PARAM period (days).
        epoch_jd (float): VOTable PARAM epoch as absolute JD (T0 for E).
        config (Config | FitConfig): Layered configuration, or its flat view (the default
            timing recipe is then recorded in the result).
        progress (callable | None): ``progress(stage, fraction)`` called after each stage,
            fraction in [0, 1]. The hook may raise ``FitCancelled`` to stop the fit.

    Returns:
        FitPipelineResult: Adopted and diagnostic fit products.

    Raises:
        NoValidSolution: If multi-start finds no valid solution.
        InvalidLightCurve: If the folded secondary is not an eclipse.
        ConfigError: If the profile is not one of the Table 2 profiles 1 to 4 (the staged
            fit frees Gamma) or the K term is requested.
    """

    def notify(stage: str, fraction: float) -> None:
        """Report a finished stage to the optional progress hook."""
        if progress is not None:
            progress(stage, fraction)

    if isinstance(config, Config):
        layered = config.validate()
        config = layered.to_fit_config()
    else:
        layered = Config.from_fit_config(config)
    if config.kterm:
        raise ConfigError(
            "The staged pipeline does not free the K term (Eq. 14); fit it with fit_global."
        )
    if config.profile not in (1, 2, 3, 4):
        raise ConfigError(
            f"run_staged_fit supports Table 2 profiles 1 to 4 (they have Gamma); "
            f"got profile {config.profile}. Fit profiles 5 to 7 with fit_global."
        )
    order = config.baseline_order
    if config.gap_split_days is None:
        seg = None
    else:
        seg = np.cumsum(np.r_[0, np.diff(t_mjd) > config.gap_split_days])
    g = initial_guess(t_mjd, mag, period_days)
    m0_jd = _mjd_to_jd(float(g["M0"]))
    phi_pri = _fold_phase_from_t0(m0_jd, epoch_jd, period_days)
    logger.info(
        "T0=%.10f JD; deepest eclipse M0=%.10f JD; phase vs T0=%.4f",
        epoch_jd,
        m0_jd,
        phi_pri,
    )
    e0 = dict(M0=g["M0"], P0=period_days, Pdot=0.0)
    oo_mask = out_of_eclipse_mask(t_mjd, e0, g["fwhm"])
    nf = noise_floor(t_mjd, mag, mask=oo_mask)
    nf_raw = noise_floor(t_mjd, mag)
    med_sig = float(np.median(sig))
    logger.info(
        "noise floor OOE=%.5f mag unmasked=%.5f median catalogue sig=%.5f",
        nf,
        nf_raw,
        med_sig,
    )
    cadence_d = float(np.median(np.diff(np.sort(t_mjd))))
    fwhm_d = g["fwhm"] * period_days
    nseg = 1 if seg is None else int(np.max(seg) + 1)
    logger.info(
        "start M0=%.6f MJD A1=%.4f A2=%.4f FWHM=%.4f phase (%.4f d) "
        "points/FWHM=%.1f cadence=%.0f s  baseline nseg=%s order=%s "
        "profile=%s harmonics=%s ocon=%s",
        g["M0"],
        g["A1"],
        g["A2"],
        g["fwhm"],
        fwhm_d,
        fwhm_d / cadence_d,
        86400.0 * cadence_d,
        nseg,
        order,
        config.profile,
        config.harmonics,
        config.ocon,
    )
    if fwhm_d / cadence_d < config.min_points_per_fwhm:
        logger.warning(
            "Fewer than %.0f points per FWHM; contacts may be exposure-smeared "
            "and D biased wide (set FitConfig.exposure_days to integrate over the exposure).",
            config.min_points_per_fwhm,
        )
    if g["A2"] < 0:
        raise InvalidLightCurve(
            f"Folded secondary depth A2={g['A2']:.4f} mag is not an eclipse."
        )
    if g["A1"] < config.primary_a2_over_a1_min_ratio * g["A2"]:
        logger.warning(
            "A1 is not much deeper than A2; period may be P/2. Continuing with given P."
        )

    if config.exposure_days is None:
        model = EBModel(
            profile=config.profile, harmonics=config.harmonics, ocon=config.ocon
        )
    else:
        nsub = config.exposure_nsub or suggest_nsub(config.exposure_days, fwhm_d)
        logger.info(
            "exposure integration: %.1f s in %d sub-samples",
            86400.0 * config.exposure_days,
            nsub,
        )
        model = EBModel(
            profile=config.profile,
            harmonics=config.harmonics,
            ocon=config.ocon,
            exposure_days=config.exposure_days,
            nsub=nsub,
            period_days=period_days,
        )
    hnames = _harmonic_free_names(model)
    amp_scale = amplitude_scale(mag)
    bounds = shape_bounds(amp_scale)
    p0 = complete_template(
        dict(A1=g["A1"], A2=g["A2"], D1=D_from_fwhm(g["fwhm"], 1.0), G1=1.0),
        model,
    )
    free_ephem = ("M0", "P0")

    r1 = fit_global(
        t_mjd,
        mag,
        sig,
        model,
        p0,
        ["A1", "A2", "D1", "G1"],
        e0,
        free_ephem,
        seg=seg,
        order=order,
        bounds=bounds,
    )
    report_fit("1 locate eclipses", r1, nf)
    notify("locate eclipses", 0.15)

    r2, rows = multistart_global(
        t_mjd,
        mag,
        sig,
        model,
        complete_template(r1.p, model),
        ["A1", "A2", "D1", "G1"] + hnames,
        r1.ephem,
        g["fwhm"],
        starts=config.gamma_starts,
        free_ephem=free_ephem,
        seg=seg,
        order=order,
        bounds=bounds,
    )
    for g0, chi2, d1, g1, st in rows:
        logger.info("  G0=%s chi2=%s D1=%s G1=%s %s", g0, chi2, d1, g1, st)
    if r2 is None:
        raise NoValidSolution("Multi-start returned no valid (D>0, G>0) solution.")
    report_fit("2 best multi-start", r2, nf)
    notify("multi-start", 0.35)

    r3 = fit_global(
        t_mjd,
        mag,
        sig,
        model,
        complete_template(r2.p, model),
        ["A1", "A2", "D1", "G1", "C1"] + hnames,
        r2.ephem,
        free_ephem,
        seg=seg,
        order=order,
        bounds=bounds,
    )
    report_fit("3 add C1", r3, nf)
    notify("curvature C1", 0.5)
    d_rel = r3.err["D1"] / abs(r3.p["D1"])
    g_rel = r3.err["G1"] / abs(r3.p["G1"])
    c1_weak = abs(r3.p["C1"]) < 3.0 * r3.err["C1"]
    if c1_weak or d_rel > 0.10 or g_rel > 0.10:
        logger.info("Fix C1=C2=0 (weak C1 or large D/G errors).")
        p_fix = complete_template(r2.p, model)
        p_fix["C1"] = 0.0
        p_fix["C2"] = 0.0
        r = fit_global(
            t_mjd,
            mag,
            sig,
            model,
            p_fix,
            ["A1", "A2", "D1", "G1"] + hnames,
            r2.ephem,
            free_ephem,
            seg=seg,
            order=order,
            bounds=bounds,
        )
        report_fit("3b C1=C2=0", r, nf)
        free_shape = ["A1", "A2", "D1", "G1", "D2", "G2"] + hnames
    else:
        r = r3
        free_shape = ["A1", "A2", "D1", "G1", "C1", "D2", "G2", "C2"] + hnames

    r4 = fit_global(
        t_mjd,
        mag,
        sig,
        model,
        complete_template(r.p, model),
        free_shape,
        r.ephem,
        free_ephem,
        seg=seg,
        order=order,
        bounds=bounds,
    )
    report_fit("4 separate eclipse shapes", r4, nf)
    notify("separate eclipse shapes", 0.6)
    if r.chi2 > 0 and (r.chi2 - r4.chi2) / r.chi2 < config.separate_shape_chi2_frac:
        logger.info(
            "Keeping shared eclipse shape (< %.0f%% chi2 drop).",
            100 * config.separate_shape_chi2_frac,
        )
    else:
        r = r4
        logger.info("Keeping separate D2, G2.")

    res, th = residuals(t_mjd, mag, sig, model, r, seg=seg, order=order)
    oo = out_of_eclipse_phase_mask(th)
    mad = (
        1.4826 * np.median(np.abs(res[oo] - np.median(res[oo])))
        if np.any(oo)
        else float(np.median(np.abs(res)))
    )
    keep = np.ones(len(t_mjd), dtype=bool)
    keep[oo] = np.abs(res[oo]) < config.clip_nsig * max(mad, 1e-8)
    n_clip = int((~keep).sum())
    logger.info(
        "%s-sigma clip out of eclipse: dropped %s points mad=%.2f mmag",
        config.clip_nsig,
        n_clip,
        1e3 * mad,
    )
    t_work, mag_work, sig_work = t_mjd, mag, sig
    seg_work = seg
    if n_clip:
        free_tpl = [n for n in r.names if n not in ("M0", "P0", "Pdot")]
        seg_clip = None if seg is None else seg[keep]
        r5 = fit_global(
            t_mjd[keep],
            mag[keep],
            sig[keep],
            model,
            complete_template(r.p, model),
            free_tpl,
            r.ephem,
            free_ephem,
            seg=seg_clip,
            order=order,
            bounds=bounds,
        )
        report_fit("5 after clip", r5, nf)
        notify("outlier clip", 0.7)
        if r5.chi2 < r.chi2:
            r = r5
            t_work, mag_work, sig_work = t_mjd[keep], mag[keep], sig[keep]
            seg_work = seg_clip

    res, th = residuals(t_work, mag_work, sig_work, model, r, seg=seg_work, order=order)
    b = binned_residuals(res, th)
    logger.info(
        "single-component binned primary (mmag): %s",
        np.array2string(b["primary"], precision=2),
    )
    logger.info(
        "single-component binned secondary (mmag): %s",
        np.array2string(b["secondary"], precision=2),
    )
    oo_mask = out_of_eclipse_mask(t_work, r.ephem, g["fwhm"])
    nf = noise_floor(t_work, mag_work, mask=oo_mask)
    ratio = r.rms / nf
    logger.info("single-component rms/noise=%.2f G1=%.3f", ratio, r.p["G1"])
    r_single = r

    if config.allow_two_component and (
        config.force_two_component or ratio > config.rms_noise_two_comp_threshold
    ):
        logger.info(
            "Two-component primary (extension, not Mikulasek 2015); rms/noise > %.1f",
            config.rms_noise_two_comp_threshold,
        )
        r_two, tc_rows = refine_two_component(
            t_work,
            mag_work,
            sig_work,
            model,
            r_single,
            seg=seg_work,
            order=order,
            eclipses=config.two_comp_eclipses,
            bounds=bounds,
        )
        for wf, g2, chi2, rms_mmag, st in tc_rows:
            logger.info(
                "  start wf=%s Gb=%s chi2=%s rms=%s mmag %s",
                wf,
                g2,
                chi2,
                rms_mmag,
                st,
            )
        if r_two is None:
            logger.info("No valid two-component start; keeping single-component.")
        else:
            report_fit("6 two-component primary", r_two, nf)
            ok_chi = [
                row[2] for row in tc_rows if row[4] == "ok" and np.isfinite(row[2])
            ]
            n_agree = int(
                np.sum(np.abs(np.asarray(ok_chi) - r_two.chi2) / r_two.chi2 < 0.01)
            )
            rms_gain = 1.0 - r_two.rms / r_single.rms
            keys = ("A1", "D1", "G1", "A1b", "D1b", "G1b")
            rel = {k: r_two.err[k] / abs(r_two.p[k]) for k in keys if k in r_two.err}
            deg = two_component_degenerate(r_two.p, 1, r_single.p["A1"])
            pos = all(r_two.p[k] > 0 for k in ("D1", "G1", "D1b", "G1b"))
            rel_ok = (
                all(v < config.two_comp_max_rel_err for v in rel.values())
                if rel
                else False
            )
            check123 = (
                (config.force_two_component or rms_gain >= config.two_comp_min_rms_gain_frac)
                and n_agree >= config.two_comp_min_agreeing_starts
                and pos
                and rel_ok
                and not deg
            )
            timing_ok, timing_msg = _two_component_timing_ok(
                t_work,
                mag_work,
                sig_work,
                model,
                r_single,
                r_two,
                config.two_comp_timing_agreement_sigma_max,
            )
            logger.info("  %s", timing_msg)
            res2, th2 = residuals(
                t_work, mag_work, sig_work, model, r_two, seg=seg_work, order=order
            )
            b2 = binned_residuals(res2, th2)
            logger.info(
                "  binned primary [two-comp] (mmag): %s",
                np.array2string(b2["primary"], precision=2),
            )
            logger.info(
                "  rms gain=%.1f%% agree_starts=%s check123=%s timing_ok=%s",
                100 * rms_gain,
                n_agree,
                check123,
                timing_ok,
            )
            if check123:
                r = r_two
                logger.info("Two-component adopted for template.")
                if not timing_ok:
                    logger.warning(
                        "Shape checks passed but virtual-minima timings disagree; "
                        "see check4 line above."
                    )
            else:
                logger.info("Two-component rejected; keeping single-component.")
    elif not config.allow_two_component:
        logger.info("Skipping two-component: not allowed in this star's FitConfig.")
    else:
        logger.info(
            "Skipping two-component: rms/noise %.2f <= %.1f.",
            ratio,
            config.rms_noise_two_comp_threshold,
        )

    r_catalogue = r
    if config.eclipse_weight > 1.0:
        logger.info(
            "Eclipse-weighted shape refit: freeze M0, P0 and harmonics; "
            "weight=%g (chi2 is not catalogue goodness-of-fit).",
            config.eclipse_weight,
        )
        sig_w = eclipse_reweighted_sigma(
            t_work,
            sig_work,
            r.ephem,
            config.eclipse_weight,
            config.eclipse_phase_half,
        )
        free_ecl = _shape_free_names(r, model)
        if not free_ecl:
            raise NoValidSolution(
                "No eclipse parameters left to free in the weighted refit."
            )
        r_w = fit_global(
            t_work,
            mag_work,
            sig_w,
            model,
            complete_template(r.p, model),
            free_ecl,
            r.ephem,
            (),
            seg=seg_work,
            order=order,
            bounds=bounds,
        )
        report_fit("weighted eclipse shape (H and ephemeris frozen)", r_w, nf)
        res_w, th_w = residuals(
            t_work, mag_work, sig_work, model, r_w, seg=seg_work, order=order
        )
        b_w = binned_residuals(res_w, th_w)
        logger.info(
            "unweighted binned primary after weight (mmag): %s",
            np.array2string(b_w["primary"], precision=2),
        )
        logger.info(
            "unweighted binned secondary after weight (mmag): %s",
            np.array2string(b_w["secondary"], precision=2),
        )
        r = r_w

    verdict = assess_fit(
        t_work,
        mag_work,
        sig_work,
        period_days,
        model,
        r,
        seg=seg_work,
        order=order,
        start=g,
        noise=nf,
        scale=amp_scale,
    )
    notify("done", 1.0)
    return FitPipelineResult(
        model=model,
        adopted=r,
        single_component=r_catalogue if config.eclipse_weight > 1.0 else r_single,
        seg=seg_work,
        order=order,
        t_mjd=t_work,
        mag=mag_work,
        sig=sig_work,
        noise_floor=nf,
        initial_guess=g,
        config=layered,
        verdict=verdict,
    )
