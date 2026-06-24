"""Tests for the figure orchestration layer (report.py).

report.py renders each owner module's plot_* onto a figure, SAVES it as a PNG and
returns a {key: filename} manifest. The mean-line figures are pure-matplotlib and
fast; the CFD figures need pyvista + a real solution dir, so they stay opt-in
(RUN_CFD_TESTS=1).
"""
import os

import pytest

from turbine_design import StageInputs, DesignCoeffs, run_meanline
from turbine_design import report


def test_meanline_figures_write_pngs(tmp_path):
    r = run_meanline(StageInputs(), DesignCoeffs())
    figs = report.meanline_figures(r, str(tmp_path))
    assert set(figs) == {"hs", "tri", "ann", "con"}
    for name in figs.values():
        p = tmp_path / name
        assert name.endswith(".png")
        assert p.exists() and p.stat().st_size > 3000


def test_meanline_figures_survive_infeasible(tmp_path):
    """Rendering must not crash on a non-converged / warning-bearing design."""
    r = run_meanline(StageInputs(P_shaft=120e6, mdot=300.0, pi_C=18.0,
                                 TIT_C=1300.0, n_stages=3),
                     DesignCoeffs(psi=1.5, phi=0.6, DOR=0.5))
    figs = report.meanline_figures(r, str(tmp_path))
    assert (tmp_path / figs["con"]).exists()


# --- heavy CFD figures: need pyvista + a real solution dir --------------------
_PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_VTU = os.path.join(_PROJ, "cfd_run_mcptest3", "flow.vtu")

@pytest.mark.skipif(os.environ.get("RUN_CFD_TESTS") != "1" or not os.path.exists(_VTU),
                    reason="set RUN_CFD_TESTS=1 and provide cfd_run_mcptest3/ outputs")
def test_cfd_figures_render(tmp_path):
    cfd = dict(
        blade_dat=os.path.join(_PROJ, "blade_coordinates_mcptest3.dat"),
        surface_csv=os.path.join(_PROJ, "cfd_run_mcptest3", "surface_flow.csv"),
        flow_vtu=_VTU, history_csv=os.path.join(_PROJ, "cfd_run_mcptest3", "history.csv"),
        axial_chord_m=0.14922, pitch_m=0.16946, beta2_deg=42.647, beta3_deg=69.53,
        P0rel_Pa=1521102.8, P3_Pa=1106801.1, mach_in_rel=0.30782, gamma=4 / 3,
    )
    figs = report.cfd_figures(cfd, str(tmp_path))
    assert {"cfdgeom", "mesh", "load", "surfp", "colormaps", "pw", "wake", "conv"} <= set(figs)
    for name in figs.values():
        assert (tmp_path / name).stat().st_size > 1500
