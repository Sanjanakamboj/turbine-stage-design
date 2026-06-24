"""Report orchestration layer (step 5: report.py is a pure aggregator).

This module owns *no* plot logic of its own. The figures live in the modules that
own the data:

    meanline.py  -> hs_diagram, velocity_triangles, annulus_sketch, constraint_panel
    geometry.py  -> blade section / cascade / thickness / curvature / passage / verification
    cfd.py       -> isentropic_mach, surface_pressure, field_colormap, pitchwise,
                    wake_profile, convergence

``report.py`` imports those ``plot_*`` functions, renders each onto a figure and
SAVES it as a PNG, returning a ``{key: filename}`` manifest. ``latexreport.py``
consumes the manifest and ``\\includegraphics`` the PNGs — so this is the single
rendering path (the old matplotlib ``PdfPages`` writers were retired here).
"""
from __future__ import annotations
from typing import Dict, Any, Optional, Tuple

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from . import meanline as ML
from . import geometry as G
from . import cfd as CF
# re-exported so latexreport's constraint table can label rows
from .meanline import _CONSTRAINT_LABELS  # noqa: F401

DPI = 150


# --------------------------------------------------------------------------- #
#  Rendering helpers
# --------------------------------------------------------------------------- #
def _save(draw, path: str, figsize: Tuple[float, float] = (6.2, 4.4)) -> str:
    """Render a single-axes drawing to a PNG; on failure draw an error box rather
    than crash the whole report."""
    fig, ax = plt.subplots(figsize=figsize)
    try:
        draw(ax)
    except Exception as e:  # noqa: BLE001 - the report must complete regardless
        ax.clear(); ax.axis("off")
        ax.text(0.5, 0.5, f"plot unavailable:\n{type(e).__name__}: {e}",
                ha="center", va="center", color="#b22222", fontsize=8, wrap=True)
    fig.tight_layout(); fig.savefig(path, dpi=DPI); plt.close(fig)
    return path


def _save_multi(draw, path: str, figsize: Tuple[float, float], ncols: int) -> str:
    """Render a multi-axes drawing (e.g. the 3-panel pitchwise figure) to one PNG."""
    fig, axes = plt.subplots(1, ncols, figsize=figsize, sharey=True)
    try:
        draw(axes)
    except Exception as e:  # noqa: BLE001
        for a in axes:
            a.clear(); a.axis("off")
        axes[ncols // 2].text(0.5, 0.5, f"unavailable:\n{e}", ha="center", va="center",
                              color="#b22222", fontsize=8, wrap=True)
    fig.tight_layout(); fig.savefig(path, dpi=DPI); plt.close(fig)
    return path


def _save_2x2(draw, path: str, figsize: Tuple[float, float] = (10.0, 8.5)) -> str:
    """Render a 2×2 multi-axes drawing to one PNG (e.g. the 4-panel colormap)."""
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    try:
        draw(axes)
    except Exception as e:  # noqa: BLE001
        for a in axes.flat:
            a.clear(); a.axis("off")
        axes[0, 0].text(0.5, 0.5, f"unavailable:\n{e}", ha="center", va="center",
                        color="#b22222", fontsize=8, wrap=True)
    fig.tight_layout(); fig.savefig(path, dpi=DPI); plt.close(fig)
    return path


# --------------------------------------------------------------------------- #
#  Mean-line figures (owned by meanline.py)
# --------------------------------------------------------------------------- #
def meanline_figures(result: Dict[str, Any], figdir: str) -> Dict[str, str]:
    """Render the four mean-line figures into ``figdir``; return {key: filename}."""
    os.makedirs(figdir, exist_ok=True)
    figs: Dict[str, str] = {}
    spec = [
        ("hs",  lambda ax: ML.hs_diagram(ax, result),         (7.0, 4.6)),
        ("tri", lambda ax: ML.velocity_triangles(ax, result), (6.2, 4.4)),
        ("ann", lambda ax: ML.annulus_sketch(ax, result),     (6.2, 4.4)),
        ("con", lambda ax: ML.constraint_panel(ax, result),   (7.0, 4.6)),
    ]
    for key, draw, size in spec:
        name = key + ".png"
        _save(draw, os.path.join(figdir, name), size)
        figs[key] = name
    return figs


# --------------------------------------------------------------------------- #
#  Blade-geometry figures (owned by geometry.py, via ParaBlade analytic geometry)
# --------------------------------------------------------------------------- #
def _resolve_blade(blade: Optional[Dict], blade_dat: Optional[str]):
    """Return (blade_obj, s_over_Cx, axial_chord_m) from either the design dict
    (``blade['parablade']``) or the ``<dat>.parablade.json`` sidecar; None if the
    ParaBlade design isn't available (so the caller can skip/fall back)."""
    pdict = None; Cx = None; s_over_Cx = None
    if blade and blade.get("parablade"):
        pdict = blade["parablade"]
        Cx = blade.get("axial_chord_m")
        if blade.get("pitch_m") and Cx:
            s_over_Cx = blade["pitch_m"] / Cx
    if pdict is None and blade_dat:
        side = G.load_parablade(blade_dat)
        if side:
            pdict = side.get("parablade")
            Cx = side.get("axial_chord_m", Cx)
            s_over_Cx = side.get("s_over_Cx", s_over_Cx)
    if pdict is None:
        return None
    try:
        obj = G.build_blade(pdict, s_over_Cx=s_over_Cx)
    except Exception:  # noqa: BLE001 - ParaBlade unavailable / bad design
        return None
    return obj, s_over_Cx, Cx


def blade_figures(blade: Optional[Dict], blade_dat: Optional[str], figdir: str) -> Dict[str, str]:
    """Render the rotor-blade geometry figures (section, cascade, thickness,
    curvature, passage geometry, passage width). Returns {} if the ParaBlade
    design is unavailable."""
    os.makedirs(figdir, exist_ok=True)
    resolved = _resolve_blade(blade, blade_dat)
    if resolved is None:
        return {}
    obj, s_over_Cx, _Cx = resolved
    figs: Dict[str, str] = {}

    def add(key, draw, size=(6.2, 4.4)):
        name = key + ".png"
        _save(draw, os.path.join(figdir, name), size)
        figs[key] = name

    add("sec", lambda ax: G.plot_blade_section(ax, obj))
    add("cas", lambda ax: G.plot_blade_cascade(ax, obj), (5.5, 7.0))
    add("thk", lambda ax: G.plot_thickness_distribution(ax, obj))
    add("cur", lambda ax: G.plot_curvature_distribution(ax, obj))
    if s_over_Cx:
        add("pass", lambda ax: G.plot_passage_geometry(ax, obj, s_over_Cx), (6.5, 5.5))
    return figs


# --------------------------------------------------------------------------- #
#  CFD figures (owned by cfd.py)
# --------------------------------------------------------------------------- #
def _cfd_geometry_draw(ax, cfd: Dict[str, Any]):
    """Geometry-check figure for the CFD section: the 3-blade verification view
    if the ParaBlade design is available (sidecar), else a plain .dat cascade."""
    bd = cfd["blade_dat"]; Cx = cfd["axial_chord_m"]; pitch = cfd["pitch_m"]
    resolved = _resolve_blade(None, bd)
    if resolved is not None:
        obj, _s, _Cx = resolved
        G.plot_geometry_verification(ax, obj, Cx, pitch, cfd["beta2_deg"], cfd["beta3_deg"])
        return
    # fallback: cascade outline straight from the .dat coordinates
    pres, suc, _ = CF._load_blade(bd)
    bx, by = CF._blade_outline(pres, suc)
    for k, col in ((-1, "0.7"), (1, "0.7"), (0, "#1f77b4")):
        ax.fill(bx, by + k * pitch, color=col, alpha=(0.9 if k == 0 else 0.5), ec="0.2", lw=1.0)
    ax.set_aspect("equal"); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
    ax.set_title(f"Rotor cascade (pitch = {pitch*1e3:.0f} mm, Cx = {Cx*1e3:.0f} mm)")
    ax.grid(alpha=0.2)


def cfd_figures(cfd: Dict[str, Any], figdir: str) -> Dict[str, str]:
    """Render the CFD-section figures into ``figdir``; return {key: filename}."""
    os.makedirs(figdir, exist_ok=True)
    bd, sc, vtu, hist = cfd["blade_dat"], cfd["surface_csv"], cfd["flow_vtu"], cfd["history_csv"]
    Cx, pitch = cfd["axial_chord_m"], cfd["pitch_m"]
    b3, Pt0, g = cfd["beta3_deg"], cfd["P0rel_Pa"], cfd.get("gamma", 4.0 / 3.0)
    b2 = cfd.get("beta2_deg"); Mw2 = cfd.get("mach_in_rel")
    T0r = cfd.get("T0rel_K"); Rg = cfd.get("R", 287.0)
    figs: Dict[str, str] = {}

    def add(key, draw, size=(6.2, 4.4)):
        name = key + ".png"
        _save(draw, os.path.join(figdir, name), size)
        figs[key] = name

    add("cfdgeom", lambda ax: _cfd_geometry_draw(ax, cfd), (5.5, 7.0))
    # GMSH mesh figure — inferred from the flow.vtu directory (blade_mesh.su2)
    mesh_su2 = cfd.get("mesh_su2") or os.path.join(os.path.dirname(vtu), "blade_mesh.su2")
    if os.path.exists(mesh_su2):
        add("mesh", lambda ax: CF.plot_mesh(ax, mesh_su2, bd), (7.0, 6.5))
    add("load", lambda ax: CF.isentropic_mach(ax, sc, bd, Cx, Pt0, g))
    add("surfp", lambda ax: CF.surface_pressure(ax, sc, bd, Cx, Pt0, g))
    # 4-panel colormap (Mach, Pt/Pt0, entropy, P/Pt0) replaces separate mach+pf figures
    figs["colormaps"] = "colormaps.png"
    _save_2x2(lambda axes: CF.colormaps_4panel(axes, vtu, bd, pitch, Pt0, g, T0_K=T0r, R=Rg),
              os.path.join(figdir, "colormaps.png"), (10.0, 8.5))
    figs["pw"] = "pw.png"
    _save_multi(lambda axes: CF.pitchwise(axes, vtu, bd, Cx, pitch, b3, Pt0, g),
                os.path.join(figdir, "pw.png"), (8.0, 4.6), 3)
    add("wake", lambda ax: CF.wake_profile(ax, vtu, bd, Cx, pitch, b3, g,
                                           beta2_deg=b2, mach_in_rel=Mw2, T0rel_K=T0r, R=Rg),
        (5.5, 5.5))
    add("conv", lambda ax: CF.convergence(ax, hist))
    return figs
