#!/usr/bin/env python3
"""MCP server for axial-turbine stage design.

Exposes the turbine design chain as tools an LLM can call:
  * turbine_design_meanline   - 1D mean-line design (pure Python, instant)
  * turbine_design_stage      - mean-line + rotor-blade design (workflow)
  * turbine_generate_airfoil  - build the airfoil with ParaBlade
  * turbine_run_cascade_cfd   - mesh (GMSH) + RANS solve (SU2) + post-process

Run locally over stdio:  python turbine_mcp_server.py
"""
from __future__ import annotations
from typing import Optional, List
from enum import Enum
import os
import re
import json
import datetime

from pydantic import BaseModel, Field, ConfigDict
from mcp.server.fastmcp import FastMCP, Context

from turbine_design.meanline import StageInputs, DesignCoeffs, run_meanline
from turbine_design.bladedesign import BladeChoices, design_blade

mcp = FastMCP("turbine_mcp")

# --- environment defaults (override via env vars) ---
# Repo root = the directory that holds this turbine_mcp/ package, derived from the
# script location so a fresh clone works anywhere without editing any paths.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJECT_DIR = os.environ.get("TURBINE_PROJECT_DIR", _REPO_ROOT)
# SU2 is not bundled (platform-specific binary). Default to the name on PATH;
# override with SU2_BIN=/abs/path/to/SU2_CFD if it lives elsewhere.
SU2_BIN = os.environ.get("SU2_BIN", "SU2_CFD")
PARABLADE = os.environ.get("PARABLADE_PATH",
                           os.path.join(PROJECT_DIR, "parablade-master", "parablade"))


class ResponseFormat(str, Enum):
    MARKDOWN = "markdown"
    JSON = "json"


# ---------------------------------------------------------------------------
def _coeffs(p) -> DesignCoeffs:
    return DesignCoeffs(psi=p.psi, phi=p.phi, DOR=p.DOR, eta_TT=p.eta_TT, M3=p.M3)


def _inputs(p) -> StageInputs:
    return StageInputs(P_shaft=p.P_shaft_MW * 1e6, pi_C=p.pi_C, TIT_C=p.TIT_C,
                       mdot=p.mdot, M_comb=p.M_comb, eta_isC=p.eta_isC, eta_m=p.eta_m,
                       eta_row=p.eta_row, n_stages=p.n_stages, gamma=p.gamma, R=p.R)


def _meanline_markdown(r: dict) -> str:
    d, s2, s3, a, cc = (r["derived"], r["station2"], r["station3"],
                        r["annulus"], r["constraints"])
    L = ["# Mean-line design result", ""]
    for w in r.get("warnings", []):
        L.append(f"> ⚠️ {w}")
        L.append("")
    if not r.get("feasibility", {}).get("work_converged", True):
        work_line = (f"required work **{d['specific_work_Jkg']/1e3:.1f} kJ/kg** "
                     f"(delivered only **{d['delivered_work_Jkg']/1e3:.1f} kJ/kg**)")
    else:
        work_line = f"specific work **{d['specific_work_Jkg']/1e3:.1f} kJ/kg**"
    L.append(f"- Blade speed **U = {d['U_mps']:.1f} m/s**, {work_line}, "
             f"exit static pressure **P3 = {d['P3_Pa']/1e3:.1f} kPa**")
    L.append("")
    L.append("| Quantity | Station 2 (rotor in) | Station 3 (rotor out) |")
    L.append("|---|---|---|")
    L.append(f"| Mach M | {s2['M']:.3f} | {s3['M']:.3f} |")
    L.append(f"| Rel. Mach Mw | {s2['Mrel']:.3f} | {s3['Mrel']:.3f} |")
    L.append(f"| Flow angle α [°] | {s2['alpha_deg']:.2f} | {s3['alpha_deg']:.2f} |")
    L.append(f"| Metal angle β [°] | {s2['beta_deg']:.2f} | {s3['beta_deg']:.2f} |")
    L.append(f"| Span h [m] | {a['h2_m']:.3f} | {a['h3_m']:.3f} |")
    L.append("")
    L.append(f"Turning Δβ = **{s3['turning_deg']:.1f}°**, "
             f"AN² = {a['AN2_m2rpm2']:.2e}, "
             f"mean radius {a['R_mean_m']:.3f} m.")
    L.append("")
    status = "✅ all constraints satisfied" if cc["all_pass"] else "⚠️ some constraints violated"
    L.append(f"**Design check: {status}**")
    if not cc["all_pass"]:
        bad = [k for k, v in cc.items() if isinstance(v, dict) and not v["pass"]]
        L.append("Violations: " + ", ".join(bad))
    return "\n".join(L)


LOG_PATH = os.path.join(PROJECT_DIR, "design_cases_log.xlsx")


def _slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", s)[:60] or "case"


def _case_name(params, default: str) -> str:
    return (getattr(params, "case_name", "") or default)


def _emit_design(r: dict, params, *, blade=None, blade_dat=None, title: str) -> None:
    """Append a row to the Excel case-log and (if make_report) write a LaTeX report."""
    name = _case_name(params,
                      f"psi{params.psi}_phi{params.phi}_R{params.DOR}_{params.n_stages}stg")
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    try:
        from turbine_design.caselog import append_case
        append_case(LOG_PATH, r, blade=blade, case_name=name, timestamp=ts)
        r["case_log"] = LOG_PATH
    except Exception as e:  # noqa: BLE001
        r["log_error"] = f"{type(e).__name__}: {e}"
    if params.make_report:
        try:
            from turbine_design.latexreport import write_latex_report
            out = write_latex_report(r, os.path.join(PROJECT_DIR, "reports", _slug(name)),
                                     case_name=name, blade=blade, blade_dat=blade_dat, title=title)
            r["report_pdf"] = out.get("pdf"); r["report_tex"] = out.get("tex")
            if out.get("error"):
                r["report_error"] = out["error"]
        except Exception as e:  # noqa: BLE001
            r["report_error"] = f"{type(e).__name__}: {e}"


# ---------------------------------------------------------------------------
class MeanlineInput(BaseModel):
    """Inputs for the 1D mean-line design."""
    model_config = ConfigDict(validate_assignment=True, extra="forbid")
    psi: float = Field(1.8, description="Stage loading ψ = ΔH_real/U² (HP range 1.5–2.0)", gt=0)
    phi: float = Field(0.5, description="Flow coefficient φ = Vx/U", gt=0)
    DOR: float = Field(0.35, description="Degree of reaction (0 impulse … 0.5 reaction)", ge=0, le=1)
    eta_TT: float = Field(0.87, description="Assumed total-to-total efficiency", gt=0, le=1)
    M3: float = Field(0.30, description="Target absolute exit Mach number", gt=0, lt=1)
    P_shaft_MW: float = Field(330.0, description="Net shaft (extracted) power [MW]", gt=0)
    pi_C: float = Field(20.0, description="Compressor pressure ratio", gt=1)
    TIT_C: float = Field(1400.0, description="Turbine inlet static temperature [°C]", gt=0)
    mdot: float = Field(750.0, description="Turbine mass-flow rate [kg/s]", gt=0)
    M_comb: float = Field(0.15, description="Combustor-exit Mach number", gt=0, lt=1)
    # --- design assumptions you trial (efficiencies, stage count, gas model) ---
    eta_isC: float = Field(0.92, description="Compressor isentropic efficiency", gt=0, le=1)
    eta_m: float = Field(0.95, description="Mechanical efficiency (shaft)", gt=0, le=1)
    eta_row: float = Field(0.90, description="Stator/rotor blade-row efficiency", gt=0, le=1)
    n_stages: int = Field(4, description="Number of turbine stages sharing the total power", ge=1)
    gamma: float = Field(4.0 / 3.0, description="Ratio of specific heats of the combustion gas", gt=1)
    R: float = Field(287.0, description="Gas constant [J/kg/K]", gt=0)
    response_format: ResponseFormat = Field(ResponseFormat.MARKDOWN,
                                            description="'markdown' or 'json'")
    make_report: bool = Field(True, description="Also write a sectioned LaTeX design "
                              "report (tables + h-s, triangles, annulus, constraints) "
                              "and return its PDF path")
    case_name: str = Field("", description="Optional name for this test case; used for "
                           "the report folder and as the key in the Excel case-log "
                           "(re-running the same name updates that row)")


@mcp.tool(name="turbine_design_meanline", annotations={
    "title": "1D Mean-Line Turbine Design", "readOnlyHint": True,
    "destructiveHint": False, "idempotentHint": True, "openWorldHint": False})
async def turbine_design_meanline(params: MeanlineInput) -> str:
    """Run a 1D mean-line design of an axial turbine stage.

    Computes the thermodynamic and kinematic state at the three stations
    (stator inlet, rotor inlet, rotor exit), the velocity-triangle angles, the
    annulus geometry and a design-constraint check, from the imposed boundary
    conditions and the chosen design coefficients (ψ, φ, degree of reaction).

    Returns (markdown or JSON):
      derived (U, specific work, P3), station1/2/3 (P, T, V, angles, Mach),
      annulus (areas, spans, radii, AN²) and constraints (pass/fail per limit).

    Use when: "design a turbine stage", "what are the velocity triangles for
    ψ=1.8, φ=0.5", "size the annulus". This tool is read-only and instant.
    """
    r = run_meanline(_inputs(params), _coeffs(params))
    _emit_design(r, params, title="Axial turbine mean-line design")
    if params.response_format == ResponseFormat.JSON:
        return json.dumps(r, indent=2)
    md = _meanline_markdown(r)
    if r.get("report_pdf"):
        md += f"\n\n📄 LaTeX report (sections + tables): `{r['report_pdf']}`"
    elif r.get("report_error"):
        md += f"\n\n⚠️ Report not generated: {r['report_error']}"
    if r.get("case_log"):
        md += f"\n📊 Logged to case sheet: `{r['case_log']}`"
    return md


# ---------------------------------------------------------------------------
class StageDesignInput(MeanlineInput):
    """Mean-line + rotor-blade design in one call."""
    AR: float = Field(1.0, description="Blade aspect ratio span/chord", gt=0)
    stagger_deg: float = Field(39.0, description="Stagger angle [°]")
    zweifel: float = Field(1.0, description="Zweifel loading coefficient", gt=0)


@mcp.tool(name="turbine_design_stage", annotations={
    "title": "Full Stage + Rotor-Blade Design", "readOnlyHint": True,
    "destructiveHint": False, "idempotentHint": True, "openWorldHint": False})
async def turbine_design_stage(params: StageDesignInput) -> str:
    """Design the stage end-to-end: mean-line plus rotor-blade geometry.

    Runs the 1D mean-line design, then converts the rotor velocity triangles into
    a blade definition (chord, axial chord, stagger, pitch, blade count, edge
    radii, thickness). Returns a JSON object with 'meanline' and 'blade' sections.

    Use when: "design the HP stage and the rotor blade", "how many blades and
    what pitch for this stage". Read-only and instant (no CFD).
    """
    ml = run_meanline(_inputs(params), _coeffs(params))
    s2, s3, a = ml["station2"], ml["station3"], ml["annulus"]
    bl = design_blade(s2["beta_deg"], s3["beta_deg"], s2["alpha_deg"],
                      s3["alpha_deg"], a["h2_m"], a["D_mean_m"],
                      BladeChoices(AR=params.AR, stagger_deg=params.stagger_deg,
                                   zweifel=params.zweifel),
                      span_h3_m=a["h3_m"])
    out = {"meanline": ml, "blade": bl}
    # build the rotor airfoil so the report's blade-geometry figures (section,
    # cascade, thickness, curvature) are populated for this test case
    name = _case_name(params, f"psi{params.psi}_phi{params.phi}_R{params.DOR}_{params.n_stages}stg")
    blade_dat = None
    try:
        from turbine_design.geometry import generate_airfoil
        d = os.path.join(PROJECT_DIR, "reports", _slug(name))
        os.makedirs(d, exist_ok=True)
        dat = os.path.join(d, "blade_coordinates.dat")
        generate_airfoil(bl["parablade"], bl["axial_chord_m"], out_dat=dat,
                         pitch_m=bl["pitch_m"], parablade_path=PARABLADE)
        blade_dat = dat
    except Exception as e:  # noqa: BLE001
        out["airfoil_error"] = f"{type(e).__name__}: {e}"
    _emit_design(ml, params, blade=bl, blade_dat=blade_dat,
                 title="Axial turbine stage + rotor-blade design")
    for k in ("report_pdf", "report_tex", "report_error", "case_log", "log_error"):
        if ml.get(k):
            out[k] = ml[k]
    return json.dumps(out, indent=2)


# ---------------------------------------------------------------------------
class AirfoilInput(BaseModel):
    """Build the rotor airfoil with ParaBlade."""
    model_config = ConfigDict(validate_assignment=True, extra="forbid")
    beta2_deg: float = Field(..., description="Rotor inlet relative (metal) angle [°]")
    beta3_deg: float = Field(..., description="Rotor exit relative (metal) angle [°]")
    axial_chord_m: float = Field(..., description="Axial chord [m]", gt=0)
    pitch_m: float = Field(..., description="Blade pitch [m]", gt=0)
    stagger_deg: float = Field(39.0, description="Stagger angle [°]")
    out_dat: str = Field("blade_coordinates.dat",
                         description="Output file for dimensional blade coordinates")


@mcp.tool(name="turbine_generate_airfoil", annotations={
    "title": "Generate Rotor Airfoil (ParaBlade)", "readOnlyHint": False,
    "destructiveHint": False, "idempotentHint": True, "openWorldHint": False})
async def turbine_generate_airfoil(params: AirfoilInput) -> str:
    """Construct the rotor airfoil geometry with ParaBlade and save coordinates.

    Builds a camber-line + thickness airfoil from the rotor metal angles and
    stagger, scales it to the axial chord, and writes the dimensional blade
    surface coordinates (suction + pressure) to a .dat file used by the CFD step.

    Returns JSON: LE/TE points, point count, and the output .dat path.
    Use when: "generate the blade profile", "make the airfoil for CFD".
    """
    from turbine_design.geometry import generate_airfoil
    from turbine_design.bladedesign import design_blade, BladeChoices
    bl = design_blade(params.beta2_deg, params.beta3_deg, 0.0, 0.0,
                      params.axial_chord_m * 1.0, params.pitch_m,
                      BladeChoices(stagger_deg=params.stagger_deg))
    out = os.path.join(PROJECT_DIR, params.out_dat) if not os.path.isabs(params.out_dat) else params.out_dat
    res = generate_airfoil(bl["parablade"], params.axial_chord_m, out_dat=out,
                           parablade_path=PARABLADE)
    return json.dumps({"LE": res["LE"], "TE": res["TE"],
                       "n_points": res["n_points"], "out_dat": res["out_dat"]}, indent=2)


# ---------------------------------------------------------------------------
class CfdInput(BaseModel):
    """Mesh + SU2 RANS solve + post-process the rotor cascade."""
    model_config = ConfigDict(validate_assignment=True, extra="forbid")
    blade_dat: str = Field("blade_coordinates.dat",
                           description="Dimensional blade coordinates (.dat) from generate_airfoil")
    pitch_m: float = Field(..., description="Blade pitch [m]", gt=0)
    axial_chord_m: float = Field(..., description="Axial chord [m]", gt=0)
    beta2_deg: float = Field(..., description="Rotor inlet relative angle [°]")
    beta3_deg: float = Field(..., description="Rotor exit relative angle [°]")
    T0rel_K: float = Field(..., description="Relative total temperature at rotor inlet [K]")
    P0rel_Pa: float = Field(..., description="Relative total pressure at rotor inlet [Pa]")
    P3_Pa: float = Field(..., description="Rotor exit static pressure [Pa]")
    mach_in_rel: float = Field(..., description="Rotor inlet relative Mach number Mw2 "
                               "(from the mean-line) — seeds the freestream/init state", gt=0)
    gamma: float = Field(4.0 / 3.0, description="Ratio of specific heats of the gas", gt=1)
    R: float = Field(287.0, description="Gas constant [J/kg/K]", gt=0)
    inner_iter: int = Field(4000, description="Max SU2 inner iterations", ge=100, le=20000)
    work_dir: str = Field("cfd_run", description="Working directory for mesh/solver files")
    make_report: bool = Field(True, description="Also write a sectioned LaTeX CFD "
                              "report (geometry, mesh, SU2 config, loading, Mach & "
                              "pressure fields, pitchwise & wake, convergence) and "
                              "return its path")
    case_name: str = Field("", description="Optional test-case name; matches the design "
                           "report folder and the Excel case-log row so the CFD results "
                           "are added to the same case entry")


@mcp.tool(name="turbine_run_cascade_cfd", annotations={
    "title": "Run Rotor-Cascade CFD (GMSH + SU2)", "readOnlyHint": False,
    "destructiveHint": False, "idempotentHint": True, "openWorldHint": False})
async def turbine_run_cascade_cfd(params: CfdInput, ctx: Context) -> str:
    """Mesh the three-blade rotor cascade, solve the 2D RANS flow with SU2, and
    post-process the blade loading and exit conditions.

    LONG-RUNNING (typically a few minutes): builds a GMSH mesh of the three-blade
    cascade (central viscous blade, inviscid neighbour walls), writes and runs the
    SU2 RANS config in the relative frame, then extracts the surface isentropic
    Mach loading and the mass-averaged exit flow angle / Mach number.

    Returns JSON: mesh stats, convergence (final density residual, iterations),
    and post-processed loading + exit conditions. Requires the SU2 binary and the
    gmsh/pyvista Python packages to be installed.
    """
    from turbine_design import cfd
    wd = os.path.join(PROJECT_DIR, params.work_dir) if not os.path.isabs(params.work_dir) else params.work_dir
    os.makedirs(wd, exist_ok=True)
    blade = os.path.join(PROJECT_DIR, params.blade_dat) if not os.path.isabs(params.blade_dat) else params.blade_dat

    await ctx.report_progress(0.1, 1.0, "Building GMSH mesh…")
    mesh_su2 = os.path.join(wd, "blade_mesh.su2")
    mstat = cfd.build_cascade_mesh(blade, params.pitch_m, params.axial_chord_m,
                                   params.beta2_deg, params.beta3_deg, mesh_su2)

    await ctx.report_progress(0.3, 1.0, f"Mesh: {mstat['n_nodes']} nodes. Writing SU2 config…")
    cfg = os.path.join(wd, "blade_su2.cfg")
    conv_target = -9.0
    cfd.write_su2_config(cfg, mesh_su2, params.T0rel_K, params.P0rel_Pa,
                         params.beta2_deg, params.P3_Pa, params.axial_chord_m,
                         mach_in_rel=params.mach_in_rel, gamma=params.gamma, R=params.R,
                         inner_iter=params.inner_iter, conv_minval=conv_target)

    await ctx.report_progress(0.4, 1.0, "Running SU2 (this can take a few minutes)…")
    run = cfd.run_su2(SU2_BIN, cfg, wd, conv_minval=conv_target)
    if not run.get("success"):
        return json.dumps({"error": "SU2 run failed", "details": run}, indent=2)

    await ctx.report_progress(0.9, 1.0, "Post-processing…")
    post = cfd.postprocess(run["flow_vtu"], run["surface_csv"], blade,
                           params.axial_chord_m, params.pitch_m, params.P0rel_Pa,
                           beta3_deg=params.beta3_deg, gamma=params.gamma)
    await ctx.report_progress(1.0, 1.0, "Done.")
    out = {"mesh": mstat,
           "convergence": {"final_rms_density": run["final_rms_density"],
                           "iterations": run["iterations"],
                           "converged": run["converged"],
                           "target_rms_density": conv_target},
           "results": post}
    if not run["converged"]:
        out["warning"] = (
            f"SU2 stopped at RMS density {run['final_rms_density']:.2f} after "
            f"{run['iterations']} iterations without reaching the {conv_target} target. "
            "The post-processed loading and exit conditions below are from a "
            "non-converged field and should be treated as unreliable — raise inner_iter "
            "or refine the mesh before trusting them.")

    cfd_data = {
        "blade_dat": blade, "surface_csv": run["surface_csv"],
        "flow_vtu": run["flow_vtu"], "history_csv": os.path.join(wd, "history.csv"),
        "axial_chord_m": params.axial_chord_m, "pitch_m": params.pitch_m,
        "beta2_deg": params.beta2_deg, "beta3_deg": params.beta3_deg,
        "T0rel_K": params.T0rel_K, "P0rel_Pa": params.P0rel_Pa, "P3_Pa": params.P3_Pa,
        "mach_in_rel": params.mach_in_rel, "gamma": params.gamma, "R": params.R,
        "mesh": mstat, "convergence": out["convergence"], "results": post,
    }
    name = params.case_name or f"cfd_b2{params.beta2_deg:.0f}_b3{params.beta3_deg:.0f}"
    # update the Excel case-log with the CFD results (matches the design row by name)
    try:
        from turbine_design.caselog import append_case
        ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
        append_case(LOG_PATH, cfd=cfd_data, case_name=name, timestamp=ts)
        out["case_log"] = LOG_PATH
    except Exception as e:  # noqa: BLE001
        out["log_error"] = f"{type(e).__name__}: {e}"
    if params.make_report:
        try:
            from turbine_design.latexreport import write_latex_report
            rep = write_latex_report(None, os.path.join(PROJECT_DIR, "reports", _slug(name)),
                                     case_name=name, cfd=cfd_data, title="Rotor Cascade CFD",
                                     basename="cfd_report")
            out["report_pdf"] = rep.get("pdf"); out["report_tex"] = rep.get("tex")
            if rep.get("error"):
                out["report_error"] = rep["error"]
        except Exception as e:  # noqa: BLE001 - report is a bonus, never fail the run
            out["report_error"] = f"{type(e).__name__}: {e}"
    return json.dumps(out, indent=2)


if __name__ == "__main__":
    mcp.run()
