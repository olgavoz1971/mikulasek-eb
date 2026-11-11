"""Regression tests on the real TESS light curves (numbers recorded 2026-10-29).

These run the full staged pipeline and take a few seconds each. Reference values
live in ``references/stars_manifest.json``.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from mikulasek_eb.helpers import absolute_jd_to_mjd
from mikulasek_eb_tess.io import load_tess_votable
from mikulasek_eb.pipeline import run_staged_fit
from mikulasek_eb.tom import HALF_ORBIT_WINDOW, time_both_walls

ROOT = Path(__file__).resolve().parent.parent.parent
MANIFEST = json.loads((ROOT / "references" / "stars_manifest.json").read_text())
ADOPTED_RMS_REL_TOL = 0.01


@pytest.fixture(scope="module")
def ek_com_run():
    """Fit EK Com with the settings of ``fit_ek_com_s49.py`` and return the pieces."""
    import fit_ek_com_s49 as star

    payload = load_tess_votable(star.VOTABLE_PATH)
    pipeline = run_staged_fit(
        absolute_jd_to_mjd(payload["t"]),
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        star.FIT_CONFIG,
    )
    return payload, pipeline


@pytest.mark.parametrize(
    "which,key", [("primary", "timing_primary_n"), ("secondary", "timing_secondary_n")]
)
def test_ek_com_timing_counts(ek_com_run, which, key):
    """Both-sides timing counts for EK Com match the recorded reference."""
    payload, pl = ek_com_run
    vm = time_both_walls(
        pl.t_mjd,
        pl.mag,
        pl.sig,
        pl.model,
        pl.adopted,
        HALF_ORBIT_WINDOW,
        which,
        payload["epoch_jd"],
        payload["period"],
        require_contact_walls=False,
        clip_outliers=False,
    )
    expected = MANIFEST["stars"]["EK_Com"]["regression_reference"][key]
    assert len(vm["T"]) == expected


def test_ek_com_adopted_rms(ek_com_run):
    """The adopted EK Com fit keeps its recorded rms."""
    _, pl = ek_com_run
    expected = MANIFEST["stars"]["EK_Com"]["regression_reference"]["adopted_rms_mag"]
    assert pl.adopted.rms == pytest.approx(expected, rel=ADOPTED_RMS_REL_TOL)


@pytest.fixture(scope="module")
def aa_and_run():
    """Fit AA And (TIC 455790537, sector 16) with the settings of its fit script."""
    import fit_tic455790537_spoc_s16 as star

    payload = load_tess_votable(star.VOTABLE_PATH)
    return run_staged_fit(
        absolute_jd_to_mjd(payload["t"]),
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        star.FIT_CONFIG,
    )


def test_aa_and_adopted_fit(aa_and_run):
    """AA And keeps its recorded rms, noise ratio, two-component adoption and period."""
    ref = MANIFEST["stars"]["AA_And"]["regression_reference"]
    adopted = aa_and_run.adopted
    assert adopted.rms == pytest.approx(ref["adopted_rms_mag"], rel=ADOPTED_RMS_REL_TOL)
    assert adopted.rms / aa_and_run.noise_floor == pytest.approx(
        ref["rms_over_noise"], rel=0.02
    )
    assert ("A1b" in adopted.p) == ref["two_component_adopted"]
    assert adopted.ephem["P0"] == pytest.approx(ref["P0_d"], rel=1e-8)
    assert adopted.p["A1"] == pytest.approx(ref["A1"], rel=0.02)
    assert adopted.p["A2"] == pytest.approx(ref["A2"], rel=0.02)


def test_ek_com_rejections_reach_the_output(ek_com_run, tmp_path):
    """Rejected eclipses are listed in the ToM file header, never dropped silently."""
    from mikulasek_eb.tom import write_timing_ascii

    payload, pl = ek_com_run
    rejected: list = []
    vm = time_both_walls(
        pl.t_mjd,
        pl.mag,
        pl.sig,
        pl.model,
        pl.adopted,
        HALF_ORBIT_WINDOW,
        "secondary",
        payload["epoch_jd"],
        payload["period"],
        require_contact_walls=False,
        clip_outliers=False,
        rejected=rejected,
    )
    assert [row["reason"] for row in rejected] == ["one_side_uncovered"]
    out = tmp_path / "sec.dat"
    write_timing_ascii(
        out,
        vm,
        lightcurve=Path("x.vot"),
        period=payload["period"],
        t0_jd=payload["epoch_jd"],
        method="test",
        window=HALF_ORBIT_WINDOW,
        events="secondary",
        extra_comments=[],
        rejected=rejected,
    )
    header = out.read_text().splitlines()
    assert "# rejected_events: 1 (one_side_uncovered=1)" in header
    assert any(line.startswith("# rejected E=") for line in header)


def test_layered_config_gives_the_same_fit_as_the_flat_one(ek_com_run):
    """Config and FitConfig inputs give identical results; the result records the config."""
    import fit_ek_com_s49 as star
    from mikulasek_eb import Config

    payload, flat_result = ek_com_run
    layered = Config.from_fit_config(star.FIT_CONFIG)
    result = run_staged_fit(
        absolute_jd_to_mjd(payload["t"]),
        payload["mag"],
        payload["sig"],
        payload["period"],
        payload["epoch_jd"],
        layered,
    )
    assert result.adopted.rms == flat_result.adopted.rms
    assert result.config == layered
    assert flat_result.config.config_hash() == layered.config_hash()


def test_timing_recipe_is_recorded_in_the_header(ek_com_run, tmp_path):
    """time_with_recipe equals time_both_walls, and the header carries recipe and hash."""
    from mikulasek_eb import TimingRecipe
    from mikulasek_eb.tom import time_with_recipe, write_timing_ascii

    payload, pl = ek_com_run
    recipe = TimingRecipe(require_contact_walls=False, clip_outliers=False)
    vm = time_with_recipe(
        pl.t_mjd,
        pl.mag,
        pl.sig,
        pl.model,
        pl.adopted,
        "primary",
        payload["epoch_jd"],
        payload["period"],
        recipe,
    )
    expected = MANIFEST["stars"]["EK_Com"]["regression_reference"]["timing_primary_n"]
    assert len(vm["T"]) == expected
    out = tmp_path / "pri.dat"
    write_timing_ascii(
        out,
        vm,
        lightcurve=Path("x.vot"),
        period=payload["period"],
        t0_jd=payload["epoch_jd"],
        method="test",
        window=recipe.window,
        events="primary",
        extra_comments=[],
        recipe=recipe,
    )
    header = out.read_text().splitlines()
    assert any(line.startswith("# timing_recipe: {") for line in header)
    assert any(line.startswith("# timing_recipe_hash: ") for line in header)
