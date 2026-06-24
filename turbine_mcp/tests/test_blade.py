"""Regression tests for the rotor-blade design step.

Reproduces what `turbine_design_stage` does: feed the mean-line rotor angles and
annulus into design_blade() with the default BladeChoices and pin the geometry.
"""
import pytest

from turbine_design import (StageInputs, DesignCoeffs, run_meanline,
                            BladeChoices, design_blade)


@pytest.fixture(scope="module")
def blade():
    ml = run_meanline(StageInputs(), DesignCoeffs())
    s2, s3, a = ml["station2"], ml["station3"], ml["annulus"]
    return design_blade(s2["beta_deg"], s3["beta_deg"], s2["alpha_deg"],
                        s3["alpha_deg"], a["h2_m"], a["D_mean_m"],
                        BladeChoices())


def test_geometry_headline(blade):
    assert blade["chord_m"] == pytest.approx(0.192013, rel=1e-4)
    assert blade["axial_chord_m"] == pytest.approx(0.149222, rel=1e-4)
    assert blade["pitch_m"] == pytest.approx(0.169463, rel=1e-4)
    assert blade["n_blades"] == 42


def test_axial_chord_is_chord_times_cos_stagger(blade):
    import math
    assert blade["axial_chord_m"] == pytest.approx(
        blade["chord_m"] * math.cos(math.radians(blade["stagger_deg"])), rel=1e-9)


def test_parablade_block_signs(blade):
    """ParaBlade wants stagger and exit angle negated, inlet angle as-is."""
    pb = blade["parablade"]
    assert pb["stagger"] == pytest.approx(-blade["stagger_deg"], rel=1e-9)
    assert pb["theta_in"] == pytest.approx(blade["beta2_deg"], rel=1e-9)
    assert pb["theta_out"] == pytest.approx(-blade["beta3_deg"], rel=1e-9)
    # edge radii are normalised by axial chord
    assert pb["radius_in"] == pytest.approx(blade["LE_radius_m"] / blade["axial_chord_m"], rel=1e-9)


def test_blade_count_scales_with_pitch():
    """Halving the Zweifel coefficient tightens the pitch and adds blades."""
    ml = run_meanline(StageInputs(), DesignCoeffs())
    s2, s3, a = ml["station2"], ml["station3"], ml["annulus"]
    base = design_blade(s2["beta_deg"], s3["beta_deg"], s2["alpha_deg"],
                        s3["alpha_deg"], a["h2_m"], a["D_mean_m"], BladeChoices())
    tight = design_blade(s2["beta_deg"], s3["beta_deg"], s2["alpha_deg"],
                         s3["alpha_deg"], a["h2_m"], a["D_mean_m"],
                         BladeChoices(zweifel=0.5))
    assert tight["pitch_m"] < base["pitch_m"]
    assert tight["n_blades"] > base["n_blades"]
