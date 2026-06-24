"""2D CFD of the rotor cascade: build the three-blade mesh (GMSH), write and run
the SU2 RANS config, and post-process the converged flow.

Domain: three-blade linear cascade in the relative frame. The central blade is a
viscous (no-slip) wall; the upper boundary is the pressure side of the upper
neighbour and the lower boundary the suction side of the lower neighbour, both
treated as inviscid (Euler) walls. No periodicity.
"""
from __future__ import annotations
from typing import Dict, Any, Optional, Callable
import os
import subprocess
import numpy as np


# ----------------------------------------------------------------------------
# Mesh
# ----------------------------------------------------------------------------
def build_cascade_mesh(blade_dat: str, pitch_m: float, axial_chord_m: float,
                       beta2_deg: float, beta3_deg: float, out_su2: str,
                       out_msh: Optional[str] = None, size_min_frac: float = 1 / 110,
                       size_max_frac: float = 1 / 7) -> Dict[str, Any]:
    """Build the three-blade cascade mesh and write a .su2 file."""
    import gmsh
    c = np.loadtxt(blade_dat)
    n = len(c) // 2
    pres = c[:n]                                   # pressure side, LE->TE
    suc = c[n:]                                    # suction side, LE->TE
    Cx = axial_chord_m
    xs_us, ys_us = suc[::-1, 0], suc[::-1, 1]      # middle suction, TE->LE
    xs_ls, ys_ls = pres[:, 0], pres[:, 1]          # middle pressure, LE->TE
    uw = np.column_stack([pres[:, 0], pres[:, 1] + pitch_m])   # upper wall (LE->TE)
    lw = np.column_stack([suc[:, 0], suc[:, 1] - pitch_m])     # lower wall (LE->TE)

    LE_x = 0.0
    TE_x = Cx
    ext = 0.5 * Cx
    x_in, x_out = LE_x - ext, TE_x + ext
    si = np.tan(np.radians(beta2_deg))       # inlet flow ENTERS at +beta2
    so = -np.tan(np.radians(beta3_deg))      # outlet flow LEAVES at -beta3 (downward)
    up = lambda p: np.array([x_in, p[1] + si * (x_in - p[0])])
    dn = lambda p: np.array([x_out, p[1] + so * (x_out - p[0])])
    uw_f = np.vstack([up(uw[0]), uw, dn(uw[-1])])
    lw_f = np.vstack([up(lw[0]), lw, dn(lw[-1])])

    try:
        gmsh.finalize()
    except Exception:
        pass
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("cascade3")
    geo = gmsh.model.geo
    lc = Cx / 8
    S = 4
    addpts = lambda a: [geo.addPoint(float(p[0]), float(p[1]), 0, lc) for p in a]

    pLi = geo.addPoint(float(lw_f[0, 0]), float(lw_f[0, 1]), 0, lc)
    pLo = geo.addPoint(float(lw_f[-1, 0]), float(lw_f[-1, 1]), 0, lc)
    pUi = geo.addPoint(float(uw_f[0, 0]), float(uw_f[0, 1]), 0, lc)
    pUo = geo.addPoint(float(uw_f[-1, 0]), float(uw_f[-1, 1]), 0, lc)
    pLE = geo.addPoint(LE_x, 0.0, 0, lc)
    pTE = geo.addPoint(TE_x, float(ys_ls[-1]), 0, lc)
    sp_lw = geo.addBSpline([pLi] + addpts(lw_f[1:-1:S]) + [pLo])
    sp_uw = geo.addBSpline([pUi] + addpts(uw_f[1:-1:S]) + [pUo])
    sp_suc = geo.addBSpline([pTE] + addpts(np.column_stack([xs_us, ys_us])[S:n - S:S]) + [pLE])
    sp_pre = geo.addBSpline([pLE] + addpts(np.column_stack([xs_ls, ys_ls])[S:n - S:S]) + [pTE])
    l_inlet = geo.addLine(pLi, pUi)
    l_outlet = geo.addLine(pUo, pLo)
    outer = geo.addCurveLoop([sp_lw, -l_outlet, -sp_uw, -l_inlet])
    blade = geo.addCurveLoop([-sp_suc, -sp_pre])
    surf = geo.addPlaneSurface([outer, blade])
    geo.synchronize()
    for nm, t in [("INLET", [l_inlet]), ("OUTLET", [l_outlet]), ("WALL_UP", [sp_uw]),
                  ("WALL_LO", [sp_lw]), ("WALL_MID", [sp_suc, sp_pre])]:
        gmsh.model.addPhysicalGroup(1, t, name=nm)
    gmsh.model.addPhysicalGroup(2, [surf], name="FLUID")

    df = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(df, "CurvesList", [sp_suc, sp_pre, sp_uw, sp_lw])
    gmsh.model.mesh.field.setNumber(df, "Sampling", 400)
    th = gmsh.model.mesh.field.add("Threshold")
    for kk, vv in [("InField", df), ("SizeMin", Cx * size_min_frac),
                   ("SizeMax", Cx * size_max_frac), ("DistMin", 0.004), ("DistMax", 0.05)]:
        gmsh.model.mesh.field.setNumber(th, kk, vv)
    gmsh.model.mesh.field.setAsBackgroundMesh(th)
    gmsh.model.mesh.generate(2)
    gmsh.write(out_su2)
    if out_msh:
        gmsh.write(out_msh)
    gmsh.finalize()
    # report the counts SU2 actually meshes (NPOIN/NELEM in the written file),
    # not gmsh's raw node tags, which include curve nodes SU2 does not emit.
    n_nodes, n_elements = _su2_counts(out_su2)
    return {"mesh_file": out_su2, "n_nodes": n_nodes,
            "n_elements": n_elements, "size_min_frac": float(size_min_frac)}


def _su2_counts(su2_path: str) -> tuple:
    """Read NPOIN/NELEM from a written .su2 mesh (authoritative element/point counts)."""
    n_points = n_elems = None
    with open(su2_path) as f:
        for line in f:
            if line.startswith("NPOIN="):
                n_points = int(line.split("=")[1].split()[0])
            elif line.startswith("NELEM="):
                n_elems = int(line.split("=")[1].split()[0])
            if n_points is not None and n_elems is not None:
                break
    return n_points, n_elems


# ----------------------------------------------------------------------------
# SU2 config + run
# ----------------------------------------------------------------------------
def inlet_reynolds(T0rel_K: float, P0rel_Pa: float, mach_in_rel: float,
                   axial_chord_m: float, gamma: float = 4.0 / 3.0,
                   R: float = 287.0) -> float:
    """Reynolds number on the axial chord from the rotor-inlet relative state.

    Uses an air-based Sutherland viscosity, extrapolated to the (hot) combustion
    gas temperature — approximate, but consistent and design-independent.

    Note: SU2 ignores ``REYNOLDS_NUMBER`` under ``REF_DIMENSIONALIZATION=
    DIMENSIONAL`` (viscosity comes from its own model), so this value is written
    for the record only and does not drive the solution.
    """
    g = gamma
    T_static = T0rel_K / (1 + (g - 1) / 2 * mach_in_rel ** 2)
    P_static = P0rel_Pa * (T_static / T0rel_K) ** (g / (g - 1))
    rho_in = P_static / (R * T_static)
    V_in = mach_in_rel * np.sqrt(g * R * T_static)
    mu = 1.716e-5 * (T_static / 273.15) ** 1.5 * (273.15 + 110.4) / (T_static + 110.4)
    return rho_in * V_in * axial_chord_m / mu


def write_su2_config(cfg_path: str, mesh_file: str, T0rel_K: float, P0rel_Pa: float,
                     beta2_deg: float, P3_Pa: float, axial_chord_m: float,
                     mach_in_rel: float, gamma: float = 4.0 / 3.0, R: float = 287.0,
                     inner_iter: int = 4000, conv_minval: float = -9.0,
                     reynolds_number: Optional[float] = None,
                     flow_scheme: str = "JST") -> str:
    """Write the SU2 config for the 3-blade cascade (Euler neighbour walls).

    The freestream/initial state and gas model are derived from the rotor-inlet
    relative Mach number and the gas properties passed in.

    ``flow_scheme`` selects the convective scheme for the mean flow:
      * ``"JST"`` (default) — 2nd-order central + scalar dissipation; matches the
        original reference setup that converges cleanly on a well-formed mesh.
      * ``"ROE"`` — 2nd-order upwind with a Venkatakrishnan limiter; available as
        an alternative, though on the current unstructured mesh it converged no
        better (limiter limit-cycle).
    """
    g = gamma
    dirx = np.cos(np.radians(beta2_deg))
    diry = np.sin(np.radians(beta2_deg))
    # static inlet state from the relative total state and the inlet relative Mach
    T_static = T0rel_K / (1 + (g - 1) / 2 * mach_in_rel ** 2)
    P_static = P0rel_Pa * (T_static / T0rel_K) ** (g / (g - 1))
    Re = (reynolds_number if reynolds_number is not None
          else inlet_reynolds(T0rel_K, P0rel_Pa, mach_in_rel, axial_chord_m, g, R))
    if flow_scheme.upper() == "JST":
        scheme_block = ("CONV_NUM_METHOD_FLOW= JST\n"
                        "MUSCL_FLOW= NO\n"
                        "JST_SENSOR_COEFF= (0.5, 0.02)")
    else:
        scheme_block = ("CONV_NUM_METHOD_FLOW= ROE\n"
                        "MUSCL_FLOW= YES\n"
                        "SLOPE_LIMITER_FLOW= VENKATAKRISHNAN\n"
                        "VENKAT_LIMITER_COEFF= 0.05")
    cfg = f"""SOLVER= RANS
KIND_TURB_MODEL= SA
MATH_PROBLEM= DIRECT
RESTART_SOL= NO
FLUID_MODEL= IDEAL_GAS
GAMMA_VALUE= {g:.6f}
GAS_CONSTANT= {R:.4f}
MACH_NUMBER= {mach_in_rel:.6f}
FREESTREAM_TEMPERATURE= {T_static:.4f}
FREESTREAM_PRESSURE= {P_static:.4f}
AOA= 0.0
INIT_OPTION= TD_CONDITIONS
REYNOLDS_NUMBER= {Re:.6e}
REYNOLDS_LENGTH= {axial_chord_m:.8f}
KIND_TRANS_MODEL= NONE
REF_LENGTH= {axial_chord_m:.8f}
REF_AREA= 0.0
REF_DIMENSIONALIZATION= DIMENSIONAL
INLET_TYPE= TOTAL_CONDITIONS
MARKER_INLET= (INLET, {T0rel_K:.4f}, {P0rel_Pa:.4f}, {dirx:.10f}, {diry:.10f}, 0.0)
MARKER_OUTLET= (OUTLET, {P3_Pa:.4f})
MARKER_HEATFLUX= (WALL_MID, 0.0)
MARKER_EULER= (WALL_UP, WALL_LO)
MARKER_PLOTTING= (WALL_MID)
MARKER_MONITORING= (WALL_MID)
NUM_METHOD_GRAD= GREEN_GAUSS
CFL_NUMBER= 2.0
CFL_ADAPT= YES
CFL_ADAPT_PARAM= (0.1, 1.2, 0.5, 50.0)
{scheme_block}
TIME_DISCRE_FLOW= EULER_IMPLICIT
CONV_NUM_METHOD_TURB= SCALAR_UPWIND
MUSCL_TURB= NO
TIME_DISCRE_TURB= EULER_IMPLICIT
CFL_REDUCTION_TURB= 0.5
LINEAR_SOLVER= FGMRES
LINEAR_SOLVER_PREC= ILU
LINEAR_SOLVER_ERROR= 1E-4
LINEAR_SOLVER_ITER= 15
INNER_ITER= {inner_iter}
CONV_RESIDUAL_MINVAL= {conv_minval}
CONV_STARTITER= 50
CONV_FIELD= RMS_DENSITY
MESH_FILENAME= {os.path.basename(mesh_file)}
MESH_FORMAT= SU2
OUTPUT_FILES= (PARAVIEW, SURFACE_CSV, RESTART)
RESTART_FILENAME= restart_flow.dat
VOLUME_FILENAME= flow
SURFACE_FILENAME= surface_flow
CONV_FILENAME= history
SCREEN_WRT_FREQ_INNER= 250
OUTPUT_WRT_FREQ= 1000
"""
    with open(cfg_path, "w") as f:
        f.write(cfg)
    return cfg_path


def run_su2(su2_bin: str, cfg_path: str, work_dir: str,
            on_line: Optional[Callable[[str], None]] = None,
            conv_minval: float = -9.0) -> Dict[str, Any]:
    """Run SU2_CFD and return the convergence status.

    ``converged`` is True only when the final RMS density residual actually
    reached the target ``conv_minval``; the solver hitting its iteration cap
    without meeting the target leaves it False so callers can flag the result.
    """
    if not os.path.exists(su2_bin):
        return {"success": False, "error": f"SU2 binary not found at {su2_bin}"}
    proc = subprocess.Popen([su2_bin, os.path.basename(cfg_path)], cwd=work_dir,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, bufsize=1)
    for line in proc.stdout:
        if on_line:
            on_line(line.rstrip())
    proc.wait()
    hist = os.path.join(work_dir, "history.csv")
    final_res, n_iter = None, None
    if os.path.exists(hist):
        rows = open(hist).read().strip().splitlines()
        n_iter = len(rows) - 1
        try:
            final_res = float(rows[-1].split(",")[3])
        except Exception:
            pass
    converged = final_res is not None and final_res <= conv_minval
    return {"success": proc.returncode == 0, "returncode": proc.returncode,
            "iterations": n_iter, "final_rms_density": final_res,
            "converged": converged, "conv_target": conv_minval,
            "flow_vtu": os.path.join(work_dir, "flow.vtu"),
            "surface_csv": os.path.join(work_dir, "surface_flow.csv")}


# ----------------------------------------------------------------------------
# Post-processing
# ----------------------------------------------------------------------------
def postprocess(flow_vtu: str, surface_csv: str, blade_dat: str,
                axial_chord_m: float, pitch_m: float, Pt0_Pa: float,
                beta3_deg: float, gamma: float = 4.0 / 3.0) -> Dict[str, Any]:
    """Extract blade loading and mass-averaged exit conditions from the solution."""
    import pandas as pd
    import pyvista as pv
    from scipy.spatial import cKDTree

    c = np.loadtxt(blade_dat)
    n = len(c) // 2
    pres, suc = c[:n], c[n:]
    xs_us, ys_us = suc[::-1, 0], suc[::-1, 1]
    xs_ls, ys_ls = pres[:, 0], pres[:, 1]
    Cx, TE_x, TE_y = axial_chord_m, axial_chord_m, float(ys_ls[-1])

    # --- surface isentropic Mach (loading) ---
    df = pd.read_csv(surface_csv, skipinitialspace=True)
    df.columns = df.columns.str.strip().str.replace('"', '')
    x, y = df["x"].values, df["y"].values
    rho = df["Density"].values
    mx, my = df["Momentum_x"].values, df["Momentum_y"].values
    E = df["Energy"].values
    P = (gamma - 1.0) * (E - 0.5 * (mx ** 2 + my ** 2) / rho)
    Mis = np.sqrt((2 / (gamma - 1)) * np.maximum((Pt0_Pa / P) ** ((gamma - 1) / gamma) - 1, 0))
    ds, _ = cKDTree(np.column_stack([xs_us, ys_us])).query(np.column_stack([x, y]))
    dp, _ = cKDTree(np.column_stack([xs_ls, ys_ls])).query(np.column_stack([x, y]))
    iss = ds <= dp

    # --- mass-averaged exit (one pitch around the central wake) ---
    m = pv.read(flow_vtu)
    pts = m.points
    V = np.array(m["Velocity"])
    Mach = np.array(m["Mach"])
    rho_v = np.array(m["Density"])
    exit = {}
    for frac in (0.15, 0.30):
        xq = TE_x + frac * Cx
        # the wake convects along the (downward) relative exit angle -beta3
        yc = TE_y + np.tan(np.radians(-beta3_deg)) * (xq - TE_x)
        sel = (np.abs(pts[:, 0] - xq) < 0.004) & (np.abs(pts[:, 1] - yc) < pitch_m / 2)
        if sel.sum() >= 8:
            w = rho_v[sel] * np.abs(V[sel, 0])
            ang = np.degrees(np.arctan2(V[sel, 1], V[sel, 0]))
            exit[f"x_TE+{frac:.2f}Cx"] = {
                "n_points": int(sel.sum()),
                "exit_flow_angle_deg": float(np.average(ang, weights=w)),
                "exit_mach": float(np.average(Mach[sel], weights=w)),
            }

    return {
        "loading": {
            "suction_peak_Mis": float(np.nanmax(Mis[iss])),
            "pressure_peak_Mis": float(np.nanmax(Mis[~iss])),
            "n_surface_points": int(len(x)),
        },
        "field": {
            "mach_min": float(Mach.min()), "mach_max": float(Mach.max()),
        },
        "exit_massaveraged": exit,
    }


# =========================================================================== #
#  CFD figures (step 4: cfd.py owns its own plots)
#
#  Each draws onto a caller-supplied ``ax`` so ``report.py`` can compose and save
#  them as PNGs. matplotlib / pandas / scipy / pyvista are imported lazily so the
#  mesh+solve pipeline above still runs without the plotting stack.
# =========================================================================== #
def _load_blade(blade_dat: str):
    c = np.loadtxt(blade_dat)
    n = len(c) // 2
    return c[:n], c[n:], n                      # pressure side, suction side, per-side count


def _blade_outline(pres, suc):
    x = np.concatenate([pres[:, 0], suc[::-1, 0]])
    y = np.concatenate([pres[:, 1], suc[::-1, 1]])
    return x, y


def isentropic_mach(ax, surface_csv: str, blade_dat: str, axial_chord_m: float,
                    Pt0_Pa: float, gamma: float = 4.0 / 3.0) -> None:
    """Surface isentropic-Mach loading (suction + pressure side vs x/Cx)."""
    import pandas as pd
    from scipy.spatial import cKDTree
    pres, suc, n = _load_blade(blade_dat)
    xs_us, ys_us = suc[::-1, 0], suc[::-1, 1]    # TE->LE
    xs_ls, ys_ls = pres[:, 0], pres[:, 1]        # LE->TE
    df = pd.read_csv(surface_csv, skipinitialspace=True)
    df.columns = df.columns.str.strip().str.replace('"', "")
    x, y = df["x"].values, df["y"].values
    P = (gamma - 1.0) * df["Energy"].values      # no-slip wall: V≈0 → P=(γ-1)ρE
    Mis = np.sqrt((2.0 / (gamma - 1.0)) * np.maximum((Pt0_Pa / P) ** ((gamma - 1) / gamma) - 1, 0))
    xnorm = x / axial_chord_m
    d_s, i_s = cKDTree(np.column_stack([xs_us, ys_us])).query(np.column_stack([x, y]))
    d_p, i_p = cKDTree(np.column_stack([xs_ls, ys_ls])).query(np.column_stack([x, y]))
    is_s = d_s <= d_p
    arc_s = (n - 1) - i_s[is_s]; arc_p = i_p[~is_s]

    def _sm(a):
        from scipy.signal import savgol_filter
        nn = len(a)
        if nn < 5:
            return a
        w = min(15, nn); w -= (w % 2 == 0)
        return savgol_filter(a, w, 2) if w > 2 else a

    for mask, arc, col, ls, lab in ((is_s, arc_s, "#1f77b4", "-", "suction side"),
                                    (~is_s, arc_p, "#d62728", "--", "pressure side")):
        o = np.argsort(arc)
        xv, mv = xnorm[mask][o], Mis[mask][o]
        v = np.isfinite(mv) & (mv > 0)
        ax.plot(xv[v], _sm(mv[v]), color=col, ls=ls, lw=2.2, label=lab)
    ax.set_xlabel(r"$x\,/\,C_x$")
    ax.set_ylabel(r"isentropic Mach number $M_\mathrm{is}$")
    ax.set_title("Blade Surface Isentropic Mach Number Distribution")
    ax.set_xlim(0, 1); ax.set_ylim(bottom=0)
    ax.legend(fontsize=8); ax.grid(alpha=0.25)


def surface_pressure(ax, surface_csv: str, blade_dat: str, axial_chord_m: float,
                     Pt0_Pa: float, gamma: float = 4.0 / 3.0) -> None:
    """Surface static-pressure distribution along PS and SS, in [kPa] (notebook
    style): p = (γ−1)(ρE − ½|ρV|²/ρ) from the SU2 conservative variables."""
    import pandas as pd
    from scipy.spatial import cKDTree
    pres, suc, n = _load_blade(blade_dat)
    xs_us, ys_us = suc[::-1, 0], suc[::-1, 1]
    xs_ls, ys_ls = pres[:, 0], pres[:, 1]
    df = pd.read_csv(surface_csv, skipinitialspace=True)
    df.columns = df.columns.str.strip().str.replace('"', "")
    x, y = df["x"].values, df["y"].values
    rho = df["Density"].values
    mx, my = df["Momentum_x"].values, df["Momentum_y"].values
    E = df["Energy"].values
    P = (gamma - 1.0) * (E - 0.5 * (mx ** 2 + my ** 2) / rho)   # static pressure [Pa]
    d_s, i_s = cKDTree(np.column_stack([xs_us, ys_us])).query(np.column_stack([x, y]))
    d_p, i_p = cKDTree(np.column_stack([xs_ls, ys_ls])).query(np.column_stack([x, y]))
    is_s = d_s <= d_p
    arc_s = (n - 1) - i_s[is_s]; arc_p = i_p[~is_s]
    xnorm = x / axial_chord_m
    ax.plot(xnorm[~is_s][np.argsort(arc_p)], P[~is_s][np.argsort(arc_p)] / 1e3,
            color="#1f77b4", lw=2, label="pressure side")
    ax.plot(xnorm[is_s][np.argsort(arc_s)], P[is_s][np.argsort(arc_s)] / 1e3,
            color="#ff7f0e", lw=2, label="suction side")
    ax.set_xlabel("axial position  x / Cx")
    ax.set_ylabel("static pressure [kPa]")
    ax.set_title("Static pressure along PS and SS")
    ax.legend(fontsize=8); ax.grid(alpha=0.25); ax.set_xlim(0, 1)


def _overlay_domain_walls(ax, blade_dat: str, pitch_m: float, pts_xy=None) -> None:
    """Overlay the middle blade (filled white) and the correct domain walls.

    WALL_UP = pressure side of upper neighbour blade (pres coords + pitch).
    WALL_LO = suction  side of lower  neighbour blade (suc  coords − pitch).
    Regions outside the domain (above WALL_UP, below WALL_LO) are masked white.
    """
    pres, suc, _ = _load_blade(blade_dat)
    # middle blade interior
    bx = np.r_[pres[:, 0], suc[::-1, 0]]
    by = np.r_[pres[:, 1], suc[::-1, 1]]
    ax.fill(bx, by, color="white", ec="0.2", lw=0.8, zorder=3)
    if pitch_m:
        uw_x, uw_y = pres[:, 0], pres[:, 1] + pitch_m   # WALL_UP (LE→TE)
        lw_x, lw_y = suc[:, 0],  suc[:, 1]  - pitch_m   # WALL_LO (LE→TE)
        y_top = float(np.max(uw_y)) + 0.5 * pitch_m
        y_bot = float(np.min(lw_y)) - 0.5 * pitch_m
        # mask above WALL_UP
        ax.fill(np.r_[[uw_x[0]], uw_x, [uw_x[-1]]],
                np.r_[[y_top],   uw_y, [y_top]],
                color="white", lw=0, zorder=3)
        # mask below WALL_LO
        ax.fill(np.r_[[lw_x[0]], lw_x, [lw_x[-1]]],
                np.r_[[y_bot],   lw_y, [y_bot]],
                color="white", lw=0, zorder=3)
        # draw wall lines on top
        ax.plot(uw_x, uw_y, "-", color="0.25", lw=0.9, zorder=4)
        ax.plot(lw_x, lw_y, "-", color="0.25", lw=0.9, zorder=4)


def field_colormap(ax, flow_vtu, field="Mach", label=None, blade_dat=None,
                   pitch_m=None, cmap="turbo", scale=1.0):
    """Filled-contour colormap of any VTU point field (Mach, Pressure, ...)."""
    import matplotlib.pyplot as plt
    import pyvista as pv
    m = pv.read(flow_vtu); pts = m.points
    v = np.array(m[field]) * scale
    if v.ndim > 1:
        v = np.linalg.norm(v, axis=1)
    tcf = ax.tricontourf(pts[:, 0], pts[:, 1], v, levels=40, cmap=cmap)
    plt.colorbar(tcf, ax=ax, label=(label or field), fraction=0.046, pad=0.02)
    if blade_dat:
        _overlay_domain_walls(ax, blade_dat, pitch_m or 0.0)
    ax.set_aspect("equal"); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
    ax.set_title(f"{label or field} field")


def _exit_probe(flow_vtu, blade_dat, axial_chord_m, pitch_m, beta3_deg, frac, gamma):
    """Sample one pitch across the wake at x = x_TE + frac·Cx with a dense line
    probe interpolated onto the mesh (notebook method: ``pv.Line`` + ``sample``).

    The line spans exactly one pitch and is centred on the central-blade wake
    streamline — i.e. the TE point convected downstream along the exit flow
    direction (−β3), so the wake stays inside the window. ``sample`` interpolates
    onto 400 evenly-spaced points (smooth profiles regardless of mesh density);
    points falling outside the fluid domain are flagged invalid.
    """
    import pyvista as pv
    Cx = axial_chord_m
    m = pv.read(flow_vtu)
    pres, _, _ = _load_blade(blade_dat)
    TE_x, TE_y = float(pres[-1, 0]), float(pres[-1, 1])
    x_max = m.bounds[1]
    used = frac
    for fr in (frac, 0.20, 0.12):                    # back off if past the outlet
        if TE_x + fr * Cx < x_max:
            used = fr; break
    xq = TE_x + used * Cx
    # wake centre = TE streamline convected to xq along the exit flow angle (−β3)
    y_c = TE_y + (xq - TE_x) * np.tan(np.radians(-float(beta3_deg)))
    y0, y1 = y_c - pitch_m / 2, y_c + pitch_m / 2
    line = pv.Line(pointa=(xq, y0, 0.0), pointb=(xq, y1, 0.0), resolution=400)
    res = line.sample(m)
    yp = res.points[:, 1]
    if "vtkValidPointMask" in res.array_names:
        valid = np.array(res["vtkValidPointMask"]).astype(bool)
    else:
        valid = np.ones(len(yp), bool)
    order = np.argsort(yp)
    yn = (yp[order] - y0) / pitch_m
    return dict(yn=yn, valid=valid[order], xq=xq, used_frac=used,
                P=np.array(res["Pressure"])[order], T=np.array(res["Temperature"])[order],
                vel=np.array(res["Velocity"])[order], Mach=np.array(res["Mach"])[order])


def _bin_avg(vals, yn, nbins=60):
    edges = np.linspace(0, 1, nbins + 1)
    yc = 0.5 * (edges[:-1] + edges[1:])
    out = np.full(nbins, np.nan)
    for i in range(nbins):
        m = (yn >= edges[i]) & (yn < edges[i + 1])
        if m.sum():
            out[i] = np.nanmean(vals[m])
    return yc, out


def pitchwise(axes, flow_vtu, blade_dat, axial_chord_m, pitch_m, beta3_deg,
              Pt0_Pa, gamma=4.0 / 3.0, frac=0.30) -> None:
    """Pitchwise Pt/Pt0, Mach and flow angle across the wake (3 axes).

    Raw scatter (faint) + pitch-wise bin average (bold line). The flow-angle
    panel includes a dashed reference at the design exit angle −β3.
    """
    R = 287.0; cp = gamma * R / (gamma - 1)
    d = _exit_probe(flow_vtu, blade_dat, axial_chord_m, pitch_m, beta3_deg, frac, gamma)
    v = d["valid"]
    if v.sum() < 8:
        raise RuntimeError(f"only {int(v.sum())} valid probe points")
    yn, P, T, vel, Mach = d["yn"][v], d["P"][v], d["T"][v], d["vel"][v], d["Mach"][v]
    V2 = vel[:, 0] ** 2 + vel[:, 1] ** 2
    Tt = T + V2 / (2 * cp); Pt = P * (Tt / T) ** (gamma / (gamma - 1))
    ang = np.degrees(np.arctan2(vel[:, 1], vel[:, 0]))
    pt_ratio = Pt / Pt0_Pa
    b3 = float(beta3_deg)                   # SU2 relative frame: exit leaves at −β3

    panels = [
        (axes[0], pt_ratio, r"$P_t\,/\,P_{t,\mathrm{inlet}}$", "#1f77b4", "pt"),
        (axes[1], Mach,     "Mach number",                     "#d62728", "m"),
        (axes[2], ang,      "flow angle [°]",                  "#2ca02c", "ang"),
    ]
    for ax, raw, lab, col, kind in panels:
        yb, vb = _bin_avg(raw, yn, nbins=60)
        mask = np.isfinite(vb)
        ax.plot(vb[mask], yb[mask], "-o", lw=2.0, ms=3, color=col, zorder=2)
        if kind == "pt":
            ax.axvline(1.0, color="grey", ls="--", lw=1, label=r"inlet $P_t$")
            ax.legend(fontsize=7)
        if kind == "ang":
            ax.axvline(-b3, color="grey", ls="--", lw=1.2, label=f"design $-\\beta_3$={-b3:.1f}°")
            ax.legend(fontsize=7)
        ax.set_xlabel(lab, fontsize=8); ax.grid(alpha=0.25)
    axes[0].set_ylabel("pitchwise position  y / pitch")
    fr = d.get("used_frac", frac)
    axes[1].set_title(f"Pitchwise distributions at $x_{{TE}} + {fr:.2f}\\,C_x$", fontsize=9)


def wake_profile(ax, flow_vtu, blade_dat, axial_chord_m, pitch_m, beta3_deg,
                 gamma=4.0 / 3.0, frac=0.30, beta2_deg=None, mach_in_rel=None,
                 T0rel_K=None, R=287.0) -> None:
    """Wake axial-velocity profile Vx/Vx,in with the momentum-deficit region
    shaded (notebook style). Normalises by the design inlet axial velocity when
    the inlet conditions are supplied, else by the freestream Vx of the slice."""
    d = _exit_probe(flow_vtu, blade_dat, axial_chord_m, pitch_m, beta3_deg, frac, gamma)
    v = d["valid"]
    if v.sum() < 8:
        raise RuntimeError(f"only {int(v.sum())} valid probe points")
    yn, vel = d["yn"][v], d["vel"][v]
    Vx = vel[:, 0]
    yb, vb = _bin_avg(Vx, yn, nbins=25)
    if beta2_deg is not None and mach_in_rel is not None and T0rel_K is not None:
        Tin = T0rel_K / (1 + (gamma - 1) / 2 * mach_in_rel ** 2)
        Vx_in = mach_in_rel * np.sqrt(gamma * R * Tin) * np.cos(np.radians(beta2_deg))
        xlabel = r"$V_x\,/\,V_{x,\mathrm{in}}$"
    else:
        Vx_in = float(np.nanmax(np.abs(vb[np.isfinite(vb)]))) or 1.0
        xlabel = r"$V_x\,/\,V_{x,\max}$"
    Vn = vb / Vx_in
    ok = np.isfinite(Vn)
    ax.plot(Vn[ok], yb[ok], "-o", lw=1.8, ms=4, color="#1f77b4")
    ax.axvline(1.0, color="k", ls="--", lw=1.2, label=r"inlet $V_x$")
    ax.fill_betweenx(yb[ok], Vn[ok], 1.0, where=(Vn[ok] < 1.0),
                     alpha=0.25, color="red", label="momentum deficit")
    ax.set_xlabel(xlabel); ax.set_ylabel("pitchwise position  y / pitch")
    fr = d.get("used_frac", frac)
    ax.set_title(f"Wake axial-velocity profile at x = x$_{{TE}}$ + {fr:.2f}·Cx")
    ax.legend(fontsize=8); ax.grid(alpha=0.25); ax.set_ylim(0, 1)


def convergence(ax, history_csv: str) -> None:
    """SU2 residual convergence history (rms columns vs inner iteration)."""
    import pandas as pd
    df = pd.read_csv(history_csv, skipinitialspace=True)
    df.columns = df.columns.str.strip().str.replace('"', "")
    xcol = "Inner_Iter" if "Inner_Iter" in df.columns else df.columns[0]
    for col in df.columns:
        if col.startswith("rms"):
            ax.plot(df[xcol], df[col], lw=1.2, label=col)
    ax.set_xlabel("inner iteration"); ax.set_ylabel("RMS residual (log10)")
    ax.set_title("Convergence history")
    ax.legend(fontsize=7.5, ncol=2); ax.grid(alpha=0.25)


def _parse_su2_mesh(su2_path: str):
    """Parse a SU2 text mesh; return (points[N,2], triangles[M,3]).

    SU2 format (2-D):
      NELEM= M   then M lines: ``5 n1 n2 n3 id``  (VTK_TRIANGLE = 5)
      NPOIN= N   then N lines: ``x y id``
    """
    points = []; triangles = []
    with open(su2_path) as fh:
        lines = fh.readlines()
    i = 0
    while i < len(lines):
        ln = lines[i].strip()
        if ln.startswith("NELEM="):
            nelem = int(ln.split("=")[1].split()[0])
            i += 1
            for _ in range(nelem):
                p = lines[i].split()
                if int(p[0]) == 5:          # VTK_TRIANGLE
                    triangles.append([int(p[1]), int(p[2]), int(p[3])])
                i += 1
        elif ln.startswith("NPOIN="):
            npoin = int(ln.split("=")[1].split()[0])
            i += 1
            for _ in range(npoin):
                p = lines[i].split()
                points.append([float(p[0]), float(p[1])])
                i += 1
        else:
            i += 1
    return np.array(points), np.array(triangles, dtype=int)


def plot_mesh(ax, su2_path: str, blade_dat: Optional[str] = None) -> None:
    """Cascade mesh triangulation from the SU2 mesh file.

    Plots each triangle as a wire-frame with fine lines so the mesh topology
    and refinement near the blade surface are visible (matching the GMSH output
    cells shown in the notebook).
    """
    pts, tris = _parse_su2_mesh(su2_path)
    if len(pts) == 0 or len(tris) == 0:
        ax.text(0.5, 0.5, "mesh data unavailable", ha="center", va="center",
                color="#b22222", fontsize=10, transform=ax.transAxes)
        ax.set_title("Cascade mesh"); return
    ax.triplot(pts[:, 0], pts[:, 1], tris, lw=0.25, color="#3399cc", alpha=0.65)
    if blade_dat:
        pres, suc, _ = _load_blade(blade_dat)
        bx = np.r_[pres[:, 0], suc[::-1, 0]]
        by = np.r_[pres[:, 1], suc[::-1, 1]]
        ax.fill(bx, by, color="white", ec="0.15", lw=1.2, zorder=3)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
    ax.set_title(f"GMSH cascade mesh  "
                 f"({len(pts):,} nodes,  {len(tris):,} triangles)")
    ax.grid(False)


def colormaps_4panel(axes, flow_vtu: str, blade_dat: Optional[str] = None,
                     pitch_m: Optional[float] = None, Pt0_Pa: Optional[float] = None,
                     gamma: float = 4.0 / 3.0, T0_K: Optional[float] = None,
                     R: float = 287.0) -> None:
    """2 × 2 flow-field colormap panel matching the notebook: Mach (jet), specific
    total enthalpy (plasma), entropy rise Δs (RdYlBu_r) and static pressure in bar
    (coolwarm). Uses the true mesh connectivity so the blade interiors stay as
    holes (no Delaunay fill-in), 60 filled contour levels, 2–98 percentile clip.
    """
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    import pyvista as pv

    m = pv.read(flow_vtu).triangulate()
    cells = m.cells.reshape(-1, 4)
    tris = cells[:, 1:]
    xm, ym = m.points[:, 0], m.points[:, 1]
    triang = mtri.Triangulation(xm, ym, tris)

    cp = gamma * R / (gamma - 1)
    rho = np.array(m["Density"]); P = np.array(m["Pressure"])
    T = np.array(m["Temperature"]); rhoE = np.array(m["Energy"])
    Mach = np.array(m["Mach"])
    vel = np.array(m["Velocity"])
    Tt0 = T0_K if T0_K else float(np.nanmax(T + (vel[:, 0] ** 2 + vel[:, 1] ** 2) / (2 * cp)))
    Pt0 = Pt0_Pa if Pt0_Pa else float(np.nanmax(P))
    h_t = (rhoE + P) / rho                              # specific total enthalpy [J/kg]
    ds = cp * np.log(T / Tt0) - R * np.log(P / Pt0)     # entropy rise [J/kg/K]

    panels = [
        (axes[0, 0], Mach,      "Mach Number",                                       "jet"),
        (axes[0, 1], h_t / 1e6, "Specific Total Enthalpy [MJ/kg]",                   "plasma"),
        (axes[1, 0], ds,        r"Entropy Rise $\Delta s$ [J kg$^{-1}$ K$^{-1}$]",   "RdYlBu_r"),
        (axes[1, 1], P / 1e5,   "Static Pressure [bar]",                             "coolwarm"),
    ]
    for ax, field, title, cmap in panels:
        lo, hi = np.percentile(field, [2, 98])
        tcf = ax.tricontourf(triang, field, levels=60, cmap=cmap, vmin=lo, vmax=hi)
        plt.colorbar(tcf, ax=ax, shrink=0.85, pad=0.02)
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.set_xlabel("x [m]", fontsize=8); ax.set_ylabel("y [m]", fontsize=8)
        ax.set_aspect("equal"); ax.tick_params(labelsize=7)
