"""Regression tests pinning the 1D mean-line design to the reference HP stage.

The reference numbers are the module's own output for the default inputs
(ψ=1.8, φ=0.5, DOR=0.35), which reproduce the project's LMECA2323 brief. Any
change to the physics that moves these values trips the test.
"""
import math

import pytest

from turbine_design import StageInputs, DesignCoeffs, run_meanline


@pytest.fixture(scope="module")
def ref():
    """Default-input mean-line result (the reference design)."""
    return run_meanline(StageInputs(), DesignCoeffs())


def test_derived_headline_numbers(ref):
    d = ref["derived"]
    assert d["U_mps"] == pytest.approx(349.0716, rel=1e-4)
    assert d["specific_work_Jkg"] == pytest.approx(219331.82, rel=1e-4)
    assert d["P3_Pa"] == pytest.approx(1106801.13, rel=1e-4)
    # the P3 secant loop must actually converge on the work residual
    assert abs(d["iteration_residual_Jkg"]) < 1.0


def test_station_angles_and_mach(ref):
    s2, s3 = ref["station2"], ref["station3"]
    assert s2["M"] == pytest.approx(0.69904, rel=1e-4)
    assert s2["beta_deg"] == pytest.approx(42.6474, rel=1e-4)
    assert s3["beta_deg"] == pytest.approx(69.5302, rel=1e-4)
    assert s3["alpha_deg"] == pytest.approx(34.1735, rel=1e-4)
    assert s3["turning_deg"] == pytest.approx(112.1776, rel=1e-4)
    # turning must equal beta2 + beta3 by construction
    assert s3["turning_deg"] == pytest.approx(s2["beta_deg"] + s3["beta_deg"], rel=1e-9)


def test_annulus(ref):
    a = ref["annulus"]
    assert a["R_mean_m"] == pytest.approx(1.11113, rel=1e-4)
    assert a["h3_m"] > a["h2_m"]  # diverging annulus
    assert a["span_ratio_h3_h2"] == pytest.approx(a["h3_m"] / a["h2_m"], rel=1e-9)


def test_default_constraints_flag_M2_and_h3h2(ref):
    """The default design sits on the edge: M2 just below 0.70, h3/h2 above 1.20."""
    c = ref["constraints"]
    assert c["all_pass"] is False
    assert c["M2"]["pass"] is False
    assert c["h3_h2"]["pass"] is False
    # everything else should pass at the reference point
    passing = {k for k, v in c.items()
               if isinstance(v, dict) and v["pass"]}
    assert {"beta2_deg", "Mw2", "Mw3", "beta3_deg",
            "turning_deg", "alpha3_deg", "AN2_m2rpm2"} <= passing


def test_offdesign_low_turning_flags_turning():
    """Low loading + high flow coefficient drops turning below the 110° floor."""
    r = run_meanline(StageInputs(), DesignCoeffs(psi=1.6, phi=0.6, DOR=0.4))
    c = r["constraints"]
    assert r["station3"]["turning_deg"] < 110.0
    assert c["turning_deg"]["pass"] is False


def test_offdesign_high_reaction_flags_alpha3():
    """A 50%-reaction, high-loading stage pushes exit swirl α3 past the 35° cap."""
    r = run_meanline(StageInputs(), DesignCoeffs(psi=2.0, phi=0.45, DOR=0.5))
    c = r["constraints"]
    assert r["station3"]["alpha_deg"] > 35.0
    assert c["alpha3_deg"]["pass"] is False


def test_blade_speed_capped_at_mechanical_limit():
    """Very low loading would demand U > U_limit; it must be clamped to the cap."""
    inp = StageInputs()
    r = run_meanline(inp, DesignCoeffs(psi=1.0))  # sqrt(dH/1.0) > 350
    assert r["derived"]["U_mps"] == pytest.approx(inp.U_limit, rel=1e-9)


def test_reference_design_is_feasible_and_unwarned(ref):
    """The reference case converges, doesn't hit the U cap, and raises no warning."""
    f = ref["feasibility"]
    assert f["work_converged"] is True
    assert f["U_capped"] is False
    assert abs(f["residual_Jkg"]) < 1.0
    assert ref["warnings"] == []
    # when converged, delivered work matches the required work
    assert ref["derived"]["delivered_work_Jkg"] == pytest.approx(
        ref["derived"]["specific_work_Jkg"], abs=1.0)


def test_infeasible_high_power_case_warns():
    """120 MW / 3 stages at psi=1.5 demands U>350: must clamp U AND fail work-match."""
    r = run_meanline(StageInputs(P_shaft=120e6, mdot=300.0, pi_C=18.0,
                                 TIT_C=1300.0, n_stages=3),
                     DesignCoeffs(psi=1.5, phi=0.6, DOR=0.5))
    f = r["feasibility"]
    assert f["U_capped"] is True
    assert f["work_converged"] is False
    # the headline 'required' work must exceed what the stage actually delivers
    assert r["derived"]["delivered_work_Jkg"] < r["derived"]["specific_work_Jkg"]
    # both conditions must each produce a warning string
    assert len(r["warnings"]) == 2


def test_constraint_record_shape(ref):
    """Every constraint carries value/min/max/pass and pass is internally consistent."""
    for name, rec in ref["constraints"].items():
        if name == "all_pass":
            continue
        assert set(rec) == {"value", "min", "max", "pass"}
        lo, hi = rec["min"], rec["max"]
        expected = (lo is None or rec["value"] >= lo) and (hi is None or rec["value"] <= hi)
        assert rec["pass"] is expected
