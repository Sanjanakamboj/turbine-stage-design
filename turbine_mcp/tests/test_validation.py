"""Tests for the MCP input boundary — the pydantic models on the tool signatures.

These guard the validation layer that rejects bad LLM-supplied arguments before
any physics runs. Importing turbine_mcp_server constructs the FastMCP instance
but does not start it (mcp.run() is under __main__), so import is cheap and safe.
"""
import pytest
from pydantic import ValidationError

from turbine_mcp_server import (MeanlineInput, StageDesignInput,
                                AirfoilInput, CfdInput, ResponseFormat)


def test_meanline_defaults_are_the_reference_point():
    m = MeanlineInput()
    assert m.psi == 1.8 and m.phi == 0.5 and m.DOR == 0.35
    assert m.response_format == ResponseFormat.MARKDOWN


@pytest.mark.parametrize("field,value", [
    ("psi", 0.0),       # gt=0
    ("psi", -1.0),
    ("phi", 0.0),       # gt=0
    ("M3", 1.0),        # lt=1
    ("M3", 0.0),        # gt=0
    ("DOR", 1.5),       # le=1
    ("DOR", -0.1),      # ge=0
    ("pi_C", 1.0),      # gt=1
    ("eta_TT", 1.1),    # le=1
    ("n_stages", 0),    # ge=1
    ("gamma", 1.0),     # gt=1
])
def test_meanline_rejects_out_of_bounds(field, value):
    with pytest.raises(ValidationError):
        MeanlineInput(**{field: value})


def test_meanline_forbids_unknown_fields():
    """extra='forbid' means a typo'd argument is an error, not silently ignored."""
    with pytest.raises(ValidationError):
        MeanlineInput(psii=1.8)


def test_stage_input_inherits_meanline_bounds_and_adds_blade_fields():
    s = StageDesignInput()
    assert s.psi == 1.8                 # inherited
    assert s.AR == 1.0 and s.zweifel == 1.0 and s.stagger_deg == 39.0
    with pytest.raises(ValidationError):
        StageDesignInput(zweifel=0.0)   # gt=0
    with pytest.raises(ValidationError):
        StageDesignInput(AR=-1.0)       # gt=0


def test_airfoil_required_fields():
    """beta angles, axial chord and pitch are required (no defaults)."""
    with pytest.raises(ValidationError):
        AirfoilInput()
    ok = AirfoilInput(beta2_deg=42.6, beta3_deg=69.5,
                      axial_chord_m=0.149, pitch_m=0.169)
    assert ok.stagger_deg == 39.0       # default
    assert ok.out_dat == "blade_coordinates.dat"
    with pytest.raises(ValidationError):
        AirfoilInput(beta2_deg=42.6, beta3_deg=69.5,
                     axial_chord_m=-0.1, pitch_m=0.169)  # axial_chord gt=0


def test_cfd_iteration_bounds():
    base = dict(pitch_m=0.169, axial_chord_m=0.149, beta2_deg=42.6, beta3_deg=69.5,
                T0rel_K=1577.0, P0rel_Pa=1.5e6, P3_Pa=1.1e6, mach_in_rel=0.31)
    assert CfdInput(**base).inner_iter == 4000          # default
    with pytest.raises(ValidationError):
        CfdInput(**base, inner_iter=50)                 # ge=100
    with pytest.raises(ValidationError):
        CfdInput(**base, inner_iter=50000)              # le=20000
    with pytest.raises(ValidationError):
        CfdInput(**{**base, "mach_in_rel": 0.0})        # gt=0
