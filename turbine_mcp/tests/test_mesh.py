"""Geometry regression test for the cascade mesh.

Guards the β3 sign bug: the rotor exit flow leaves at -β3 (downward), so the
outlet wall extensions must slant DOWN. A previous port used +tan(β3), extending
the outlet UP against the flow — which produced a kinked domain, a wrong Mach
field, and the SU2 limit-cycle. Needs gmsh; skipped if a blade .dat is absent.
"""
import os

import pytest

_PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _blade():
    for f in ("blade_coordinates_mcptest3.dat", "blade_coordinates.dat"):
        p = os.path.join(_PROJ, f)
        if os.path.exists(p):
            return p
    return None


def _read_su2_ys(path):
    ys = []
    with open(path) as f:
        it = iter(f)
        for line in it:
            if line.startswith("NPOIN="):
                npoin = int(line.split("=")[1].split()[0])
                for _ in range(npoin):
                    ys.append(float(next(it).split()[1]))
                break
    return ys


@pytest.mark.skipif(_blade() is None, reason="no blade .dat available to mesh")
def test_outlet_extends_downstream_with_exit_flow(tmp_path):
    pytest.importorskip("gmsh")
    from turbine_design.cfd import build_cascade_mesh
    pitch, Cx, b2, b3 = 0.16946, 0.14922, 42.647, 69.53
    out = str(tmp_path / "m.su2")
    build_cascade_mesh(_blade(), pitch, Cx, b2, b3, out)
    ymin = min(_read_su2_ys(out))
    # the lower wall must extend well below the blade (exit flow at -β3, downward).
    # With the old +tan(β3) sign the domain stayed above ≈ -0.30.
    assert ymin < -0.40, f"outlet not extending downward (ymin={ymin:.3f}) — β3 sign regressed?"
