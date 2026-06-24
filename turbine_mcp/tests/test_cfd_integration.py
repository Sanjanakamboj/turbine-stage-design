"""Heavy end-to-end CFD test — skipped by default.

This is the only test that touches GMSH, SU2 and pyvista, and it takes minutes.
It is opt-in so the default suite stays fast and pure-Python. Enable with:

    RUN_CFD_TESTS=1 pytest tests/test_cfd_integration.py

It needs the SU2_CFD binary (env SU2_BIN or the project default) and a blade
.dat. It verifies the pipeline runs and — importantly — that a non-converged
solve is reported honestly (converged flag + warning) rather than passed off as
a trustworthy result.
"""
import os

import pytest

RUN = os.environ.get("RUN_CFD_TESTS") == "1"
SU2_BIN = os.environ.get("SU2_BIN", "/Users/sanju/Downloads/bin/SU2_CFD")
PROJECT_DIR = os.environ.get(
    "TURBINE_PROJECT_DIR",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
BLADE_DAT = os.path.join(PROJECT_DIR, "blade_coordinates.dat")

pytestmark = pytest.mark.skipif(
    not (RUN and os.path.exists(SU2_BIN) and os.path.exists(BLADE_DAT)),
    reason="set RUN_CFD_TESTS=1 and ensure SU2_BIN + blade_coordinates.dat exist")


def test_cfd_pipeline_runs_and_reports_convergence_honestly():
    from turbine_design import cfd  # noqa: F401  (import deferred to keep default suite light)
    import asyncio
    from turbine_mcp_server import turbine_run_cascade_cfd, CfdInput
    import json

    class _Ctx:  # minimal stand-in for the MCP Context (only report_progress is used)
        async def report_progress(self, *a, **k):
            return None

    params = CfdInput(
        blade_dat="blade_coordinates.dat",
        pitch_m=0.16946, axial_chord_m=0.14922,
        beta2_deg=42.647, beta3_deg=69.53,
        T0rel_K=1577.47, P0rel_Pa=1521102.8, P3_Pa=1106801.1,
        mach_in_rel=0.30782, inner_iter=300, work_dir="cfd_run_pytest")

    out = json.loads(asyncio.run(turbine_run_cascade_cfd(params, _Ctx())))

    assert out["mesh"]["n_nodes"] > 1000
    assert "converged" in out["convergence"]
    # A short 300-iter run will not hit the -9 target; the tool must say so.
    if not out["convergence"]["converged"]:
        assert "warning" in out
    assert "exit_massaveraged" in out["results"]
