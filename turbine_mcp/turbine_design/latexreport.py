"""Generate a sectioned LaTeX design report (compiled to PDF) for one test case.

Structure (tables-first, no intro/method prose):
  1. Mean-line design   — input parameters, design coefficients, station
       quantities, h-s diagram, annulus sizing, constraint check
  2. Rotor blade design — blade parameters, section, cascade, passage geometry,
       thickness & curvature distributions
  3. CFD evaluation     — geometry check, mesh, SU2 configuration, Mach & loading,
       field colormaps, pitchwise distributions, wake

Figures come from turbine_design.report (data-driven, per test case). Requires a
LaTeX toolchain (latexmk or pdflatex). Falls back to leaving the .tex if no
compiler is present.
"""
from __future__ import annotations
import os
import shutil
import subprocess
from typing import Dict, Any, Optional

from . import report as R
from . import meanline as ML


# --------------------------------------------------------------------------- #
_UNI = [("°", r"$^\circ$"), ("ψ", r"$\psi$"), ("φ", r"$\phi$"), ("α", r"$\alpha$"),
        ("β", r"$\beta$"), ("Δ", r"$\Delta$"), ("²", r"$^2$"), ("³", r"$^3$"),
        ("→", r"$\to$"), ("·", r"$\cdot$"), ("≤", r"$\le$"), ("≥", r"$\ge$"),
        ("—", "---"), ("–", "--")]


def _u2l(s) -> str:
    """Convert unicode symbols to LaTeX; safe on strings that already hold LaTeX
    (author-controlled table cells use intentional math like $\\pi_C$)."""
    s = str(s)
    for a, b in _UNI:
        s = s.replace(a, b)
    return s


def _textsafe(s) -> str:
    """Full escape for FREE text (title, case name, warnings) — no intentional LaTeX."""
    s = str(s)
    for a, b in (("\\", r"\textbackslash{}"), ("{", r"\{"), ("}", r"\}"),
                 ("_", r"\_"), ("%", r"\%"), ("#", r"\#"), ("&", r"\&"),
                 ("$", r"\$"), ("~", r"\textasciitilde{}"), ("^", r"\textasciicircum{}")):
        s = s.replace(a, b)
    return _u2l(s)


def _table(rows, header=None, caption="", label="", colspec=None):
    n = len(rows[0])
    colspec = colspec or ("l" + "r" * (n - 1))
    out = [r"\begin{table}[H]", r"\centering", r"\begin{tabular}{%s}" % colspec, r"\toprule"]
    if header:
        out.append(" & ".join(r"\textbf{%s}" % _u2l(h) for h in header) + r" \\")
        out.append(r"\midrule")
    for row in rows:
        out.append(" & ".join(_u2l(c) for c in row) + r" \\")
    out += [r"\bottomrule", r"\end{tabular}"]
    if caption:
        out.append(r"\caption{%s}" % _u2l(caption))
    if label:
        out.append(r"\label{%s}" % label)
    out.append(r"\end{table}")
    return "\n".join(out)


def _figure(rel, caption="", label="", width=0.82):
    out = [r"\begin{figure}[H]", r"\centering",
           r"\includegraphics[width=%s\textwidth]{%s}" % (width, rel)]
    if caption:
        out.append(r"\caption{%s}" % _u2l(caption))
    if label:
        out.append(r"\label{%s}" % label)
    out.append(r"\end{figure}")
    return "\n".join(out)


def _subfig_row(items, caption="", label="", widths=None):
    """Lay two (or more) figures side by side as subfigures.

    ``items`` is a list of ``(rel_path, subcaption)``; ``widths`` an optional list
    of \\textwidth fractions (defaults to an even split).
    """
    n = len(items)
    widths = widths or [round(0.96 / n, 3)] * n
    out = [r"\begin{figure}[H]", r"\centering"]
    for i, (rel, sub) in enumerate(items):
        out.append(r"\begin{subfigure}[t]{%.3f\textwidth}" % widths[i])
        out.append(r"\centering")
        out.append(r"\includegraphics[width=\textwidth]{%s}" % rel)
        if sub:
            out.append(r"\caption{%s}" % _u2l(sub))
        out.append(r"\end{subfigure}")
        if i < n - 1:
            out.append(r"\hfill")
    if caption:
        out.append(r"\caption{%s}" % _u2l(caption))
    if label:
        out.append(r"\label{%s}" % label)
    out.append(r"\end{figure}")
    return "\n".join(out)


def _f(v, p="{:.3f}"):
    try:
        return p.format(v)
    except Exception:
        return str(v)


# --------------------------------------------------------------------------- #
def _meanline_section(result, figdir, relfig):
    d, a, c = result["derived"], result["annulus"], result["constraints"]
    co, inp = result["coefficients"], result["inputs"]
    s1, s2, s3 = result["station1"], result["station2"], result["station3"]
    L = [r"\section{Mean-line design}"]

    for w in result.get("warnings", []):
        L.append(r"\noindent\fbox{\parbox{\linewidth}{\textbf{Warning:} %s}}\medskip" % _textsafe(w))

    L.append(r"\subsection{Input parameters}")
    L.append(_table(
        [["Shaft power $P$ [MW]", _f(inp["P_shaft"] / 1e6, "{:.1f}")],
         ["Mass flow $\\dot m$ [kg/s]", _f(inp["mdot"], "{:.1f}")],
         ["Compressor PR $\\pi_C$", _f(inp["pi_C"], "{:.1f}")],
         ["TIT [$^\\circ$C]", _f(inp["TIT_C"], "{:.0f}")],
         ["Rotational speed $N$ [rpm]", _f(inp["N_rpm"], "{:.0f}")],
         ["Stages", str(inp["n_stages"])],
         ["$\\gamma$, $R$ [J/kg/K]", "%s, %s" % (_f(inp["gamma"], "{:.3f}"), _f(inp["R"], "{:.0f}"))],
         ["Blade-speed limit [m/s]", _f(inp["U_limit"], "{:.0f}")]],
        header=["Quantity", "Value"], caption="Imposed boundary conditions.",
        label="tab:inputs"))

    L.append(r"\subsection{Design coefficients}")
    L.append(_table(
        [["Loading $\\psi$", _f(co["psi"], "{:.3f}")],
         ["Flow coefficient $\\phi$", _f(co["phi"], "{:.3f}")],
         ["Reaction $R$", _f(co["DOR"], "{:.3f}")],
         ["Blade speed $U$ [m/s]", _f(d["U_mps"], "{:.1f}")],
         ["Specific work [kJ/kg]", _f(d["specific_work_Jkg"] / 1e3, "{:.1f}")]],
        header=["Coefficient", "Value"], caption="Chosen design coefficients and derived loading.",
        label="tab:coeffs"))

    L.append(r"\subsection{Thermodynamic station quantities}")

    def srow(label, key, fmt="{:.1f}"):
        def cell(s):
            v = s.get(key)
            return "--" if v is None else _f(v, fmt)
        return [label, cell(s1), cell(s2), cell(s3)]

    body = [
        srow("Static pressure $P$ [Pa]", "P_Pa", "{:.0f}"),
        srow("Total pressure (abs.) $P_0$ [Pa]", "P0_Pa", "{:.0f}"),
        srow("Rel. total pressure $P_{0R}$ [Pa]", "P0rel_Pa", "{:.0f}"),
        srow("Static temperature $T$ [K]", "T_K", "{:.1f}"),
        srow("Total temperature (abs.) $T_0$ [K]", "T0_K", "{:.1f}"),
        srow("Rel. total temperature $T_{0R}$ [K]", "T0rel_K", "{:.1f}"),
        srow("Absolute velocity $V$ [m/s]", "V_mps", "{:.1f}"),
        srow("Axial velocity $V_x$ [m/s]", "Vx_mps", "{:.1f}"),
        srow("Tangential velocity $V_u$ [m/s]", "Vu_mps", "{:.1f}"),
        srow("Relative velocity $W$ [m/s]", "W_mps", "{:.1f}"),
        srow("Yaw angle $\\alpha$ [$^\\circ$]", "alpha_deg", "{:.2f}"),
        srow("Relative angle $\\beta$ [$^\\circ$]", "beta_deg", "{:.2f}"),
        srow("Absolute Mach $M$", "M", "{:.3f}"),
        srow("Relative Mach $M_w$", "Mrel", "{:.3f}"),
    ]
    L.append(_table(body, header=["Quantity", "Station 1", "Station 2", "Station 3"],
                    caption="Thermodynamic and kinematic quantities at the three stations (converged design).",
                    label="tab:stations"))

    figs = R.meanline_figures(result, figdir)

    L.append(r"\subsection{Velocity triangles}")
    L.append(_figure(relfig(figs["tri"]), "Velocity triangles at rotor inlet (2) and exit (3).", "fig:tri", 0.7))

    L.append(r"\subsection{h--s diagram}")
    lad = ML.enthalpy_ladder(result)

    def _J(x):
        return _f(x, "{:.0f}")

    erows = [
        [r"$H_{01}=H_{02}$", _J(lad["H01"]), r"$h_3$", _J(lad["h3"])],
        [r"$h_1$", _J(lad["h1"]), r"$h_{3s}$", _J(lad["h3s"])],
        [r"$h_2$", _J(lad["h2"]), r"$h_{3ss}$", _J(lad["h3ss"])],
        [r"$h_{2s}$", _J(lad["h2s"]), r"$H_{03}$", _J(lad["H03"])],
        [r"$H_{0R}$", _J(lad["H0R"]), r"$H_{03ss}$", _J(lad["H03ss"])],
    ]
    L.append(_table(erows, header=["Quantity", "Value", "Quantity", "Value"],
                    caption="Enthalpy ladder of the stage (J/kg).", label="tab:hladder",
                    colspec="lrlr"))
    L.append(_figure(relfig(figs["hs"]), "h--s diagram of the first turbine stage.", "fig:hs", 0.78))

    L.append(r"\subsection{Annulus sizing}")
    L.append(_table(
        [["Mean radius $R_m$ [m]", _f(a["R_mean_m"])],
         ["Span $h_2$ [mm]", _f(a["h2_m"] * 1e3, "{:.1f}")],
         ["Span $h_3$ [mm]", _f(a["h3_m"] * 1e3, "{:.1f}")],
         ["Tip radius st.3 [m]", _f(a["Rtip3_m"])],
         ["Hub radius st.3 [m]", _f(a["Rhub3_m"])],
         ["$AN^2$ [m$^2$ rpm$^2$]", _f(a["AN2_m2rpm2"], "{:.3e}")]],
        header=["Quantity", "Value"], caption="Annulus geometry (constant mean radius).",
        label="tab:annulus"))
    L.append(_figure(relfig(figs["ann"]), "Meridional annulus from station 2 to 3.", "fig:ann", 0.58))

    L.append(r"\subsection{Design constraint check}")
    crows = []
    for k, v in c.items():
        if not isinstance(v, dict):
            continue
        lim = ("[%g, %g]" % (v["min"], v["max"]) if v["min"] is not None and v["max"] is not None
               else "<= %g" % v["max"] if v["max"] is not None
               else ">= %g" % v["min"] if v["min"] is not None else "-")
        crows.append([R._CONSTRAINT_LABELS.get(k, k), _f(v["value"], "{:.4g}"), lim,
                      "PASS" if v["pass"] else "FAIL"])
    L.append(_table(crows, header=["Constraint", "Value", "Allowed", "Status"],
                    caption="Design-constraint verification (%s)." %
                    ("all satisfied" if c["all_pass"] else "violations present"),
                    label="tab:constraints"))
    L.append(_figure(relfig(figs["con"]), "Each constraint vs its allowed band.", "fig:con", 0.8))
    return "\n\n".join(L)


def _blade_section(blade, blade_dat, figdir, relfig):
    L = [r"\section{Rotor blade design}"]
    pitch = blade["pitch_m"]; Cx = blade["axial_chord_m"]
    L.append(r"\subsection{Main blade parameters}")
    brows = [
        ["Inlet flow angle $\\alpha_2$ [$^\\circ$]", _f(blade["alpha2_deg"], "{:.2f}")],
        ["Inlet metal (rel.) angle $\\beta_2$ [$^\\circ$]", _f(blade["beta2_deg"], "{:.2f}")],
        ["Exit flow angle $\\alpha_3$ [$^\\circ$]", _f(blade["alpha3_deg"], "{:.2f}")],
        ["Exit metal (rel.) angle $\\beta_3$ [$^\\circ$]", _f(blade["beta3_deg"], "{:.2f}")],
        ["Aspect ratio $AR$", _f(blade["aspect_ratio"], "{:.1f}")],
        ["Stagger angle $\\xi$ [$^\\circ$]", _f(blade["stagger_deg"], "{:.0f}")],
        ["Zweifel coefficient $Z$", _f(blade["zweifel"], "{:.1f}")],
        ["Chord $c$ [m]", _f(blade["chord_m"], "{:.3f}")],
        ["Axial chord $c_x$ [m]", _f(Cx, "{:.3f}")],
        ["Pitch $s$ [m]", _f(pitch, "{:.3f}")],
        ["Pitch-to-axial-chord $s/c_x$", _f(blade["pitch_to_axial_chord"], "{:.2f}")],
        ["Number of blades $N_b$", str(blade["n_blades"])],
        ["Max. thickness-to-chord $t_{max}/c$", _f(blade["t_max_over_chord"], "{:.2f}")],
        ["Leading-edge radius $R_{LE}$ [m]", _f(blade["LE_radius_m"], "{:.4f}")],
        ["Trailing-edge radius $R_{TE}$ [m]", _f(blade["TE_radius_m"], "{:.4f}")],
    ]
    if blade.get("flare_angle_deg") is not None:
        brows.append(["End-wall flare angle $\\varepsilon$ [$^\\circ$]",
                      _f(blade["flare_angle_deg"], "{:.2f}")])
    L.append(_table(brows, header=["Parameter", "Value"],
                    caption="Main rotor-blade parameters.", label="tab:blade"))
    figs = R.blade_figures(blade, blade_dat, figdir)
    if figs.get("sec"):
        L.append(_figure(relfig(figs["sec"]),
                         "Constructed rotor-blade section (camber line and thickness distribution).",
                         "fig:sec", 0.52))
    if figs.get("cas") and figs.get("pass"):
        L.append(_subfig_row(
            [(relfig(figs["cas"]), "Cascade of rotor blades."),
             (relfig(figs["pass"]), "Blade-to-blade passage and throat.")],
            "Rotor-blade cascade and passage geometry generated by ParaBlade.",
            "fig:cascade", widths=[0.40, 0.56]))
    elif figs.get("cas"):
        L.append(_figure(relfig(figs["cas"]), "Cascade of rotor blades.", "fig:cas", 0.5))
    if figs.get("thk") and figs.get("cur"):
        L.append(_subfig_row(
            [(relfig(figs["thk"]), "Thickness distribution."),
             (relfig(figs["cur"]), "Curvature $\\kappa$ distribution.")],
            "Profile thickness and curvature along the normalised surface coordinate.",
            "fig:thkcur"))
    return "\n\n".join(L)


def _cfd_section(cfd, figdir, relfig):
    L = [r"\section{CFD evaluation}"]
    g, Pt0 = cfd.get("gamma", 4 / 3), cfd["P0rel_Pa"]
    conv = cfd.get("convergence", {}); mesh = cfd.get("mesh", {})
    figs = R.cfd_figures(cfd, figdir)

    L.append(r"\subsection{Geometry check and mesh}")
    if figs.get("mesh"):
        L.append(_subfig_row(
            [(relfig(figs["cfdgeom"]), "Three-blade cascade layout."),
             (relfig(figs["mesh"]), "GMSH mesh.")],
            "Computational domain (central blade = viscous wall; neighbour surfaces = "
            "inviscid walls) and mesh.", "fig:domain", widths=[0.42, 0.54]))
    else:
        L.append(_figure(relfig(figs["cfdgeom"]), "Three-blade cascade domain (CFD).", "fig:cfdgeom", 0.5))
    L.append(_table([["Nodes", str(mesh.get("n_nodes", "--"))],
                     ["Elements", str(mesh.get("n_elements", "--"))]],
                    header=["Mesh", "Count"], caption="GMSH cascade mesh.", label="tab:mesh"))

    L.append(r"\subsection{SU2 configuration}")
    L.append(_table(
        [["Solver", "RANS (Spalart--Allmaras)"],
         ["Convective scheme", "JST (central + scalar diss.)"],
         ["Gas model", "ideal, $\\gamma=%s$, $R=%s$" % (_f(g, "{:.3f}"), _f(cfd.get("R", 287), "{:.0f}"))],
         ["Inlet (total, rel.)", "$T_0=%s$ K, $P_0=%s$ kPa" % (_f(cfd.get("T0rel_K", float('nan')), "{:.0f}"), _f(Pt0 / 1e3, "{:.1f}"))],
         ["Inlet flow angle $\\beta_2$ [$^\\circ$]", _f(cfd["beta2_deg"], "{:.2f}")],
         ["Outlet static $P_3$ [kPa]", _f(cfd["P3_Pa"] / 1e3, "{:.1f}")],
         ["Inlet rel. Mach $M_{w2}$", _f(cfd["mach_in_rel"], "{:.3f}")],
         ["Neighbour walls", "Euler (inviscid)"],
         ["Convergence target", "RMS $\\rho$ = %s" % str(conv.get("target_rms_density", -9.0))],
         ["Final RMS $\\rho$ / iters", "%s / %s" % (_f(conv.get("final_rms_density", float('nan')), "{:.2f}"), str(conv.get("iterations", "--")))],
         ["Converged", "yes" if conv.get("converged") else "NO"]],
        header=["Setting", "Value"], caption="SU2 RANS configuration and convergence.",
        label="tab:su2"))

    L.append(r"\subsection{Blade loading and field}")
    L.append(_subfig_row(
        [(relfig(figs["load"]), "Isentropic-Mach loading."),
         (relfig(figs["surfp"]), "Surface static pressure.")],
        "Blade-surface loading: isentropic Mach number and static-pressure distribution.",
        "fig:loading"))
    if figs.get("colormaps"):
        L.append(_figure(relfig(figs["colormaps"]),
                         r"Flow-field colormaps: Mach number, specific total enthalpy, "
                         r"entropy rise $\Delta s$, and static pressure.",
                         "fig:colormaps", 0.9))

    L.append(r"\subsection{Wake and pitchwise distributions}")
    L.append(_figure(relfig(figs["pw"]),
                     r"Pitchwise $P_t/P_{0,rel}$, Mach and flow angle at $x_{TE}+0.3\,C_x$.",
                     "fig:pwcfd", 0.88))
    L.append(_subfig_row(
        [(relfig(figs["wake"]), "Wake axial-velocity profile."),
         (relfig(figs["conv"]), "SU2 residual convergence.")],
        "Downstream wake profile and solver convergence history.",
        "fig:wakeconv", widths=[0.40, 0.56]))
    return "\n\n".join(L)


_PREAMBLE = r"""\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{graphicx}
\usepackage{amsmath,amssymb}
\usepackage[margin=1in]{geometry}
\usepackage{booktabs}
\usepackage{float}
\usepackage{caption}
\usepackage{subcaption}
\usepackage{hyperref}
\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=blue}
\captionsetup{font=small,labelfont=bf}
\captionsetup[sub]{font=small,labelfont=bf}
\setlength{\parindent}{0pt}
\title{%(title)s}
\author{Turbine MCP --- automated design report}
\date{Case: %(case)s}
\begin{document}
\maketitle
\tableofcontents
\newpage
"""


def write_latex_report(result: Optional[Dict[str, Any]], out_dir: str, case_name: str = "case",
                       *, blade: Optional[Dict] = None, blade_dat: Optional[str] = None,
                       cfd: Optional[Dict] = None,
                       title: str = "Axial Turbine Stage Design Report",
                       basename: str = "report", compile_pdf: bool = True) -> Dict[str, Any]:
    """Write (and optionally compile) a sectioned LaTeX report. Returns paths/status."""
    os.makedirs(out_dir, exist_ok=True)
    figdir = os.path.join(out_dir, "figures")
    os.makedirs(figdir, exist_ok=True)

    def relfig(name):
        return "figures/" + name

    body = [_PREAMBLE % {"title": _textsafe(title), "case": _textsafe(case_name)}]
    if result is not None:
        body.append(_meanline_section(result, figdir, relfig))
    if blade is not None:
        body.append(_blade_section(blade, blade_dat, figdir, relfig))
    if cfd is not None:
        body.append(_cfd_section(cfd, figdir, relfig))
    body.append(r"\end{document}")

    tex_path = os.path.join(out_dir, basename + ".tex")
    with open(tex_path, "w") as f:
        f.write("\n\n".join(body))

    out = {"tex": tex_path, "compiled": False, "pdf": None}
    if compile_pdf:
        engine = shutil.which("latexmk") or shutil.which("pdflatex")
        if engine:
            try:
                texname = basename + ".tex"
                if engine.endswith("latexmk"):
                    cmd = [engine, "-pdf", "-interaction=nonstopmode", "-halt-on-error", texname]
                else:
                    cmd = [engine, "-interaction=nonstopmode", "-halt-on-error", texname]
                for _ in range(1 if engine.endswith("latexmk") else 2):
                    subprocess.run(cmd, cwd=out_dir, stdout=subprocess.DEVNULL,
                                   stderr=subprocess.DEVNULL, timeout=180)
                pdf = os.path.join(out_dir, basename + ".pdf")
                if os.path.exists(pdf):
                    out["compiled"] = True; out["pdf"] = pdf
                else:
                    out["error"] = "compiler ran but report.pdf not produced (see report.log)"
            except Exception as e:  # noqa: BLE001
                out["error"] = f"{type(e).__name__}: {e}"
        else:
            out["error"] = "no latexmk/pdflatex found; .tex written but not compiled"
    return out
