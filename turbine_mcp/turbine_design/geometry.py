"""ParaBlade geometry layer: build the rotor airfoil from the blade-design
parameters and own every *geometry* figure (step 3 of the report chain).

Responsibilities
----------------
* ``generate_airfoil`` — build the 2D ParaBlade section, scale to the axial
  chord and write the dimensional ``.dat`` (used by the CFD mesher). It also
  drops a tiny ``<dat>.parablade.json`` sidecar holding the ParaBlade design
  variables so any later consumer (notably the CFD report, which only knows the
  ``.dat`` path) can rebuild the analytic blade.
* ``build_blade`` — construct the ParaBlade ``Blade2DCamberThickness`` object
  once, set the cascade spacing, and return it for plotting.
* ``passage_metrics`` — analytic passage-gap / geometric-throat distribution.
* ``plot_*`` — the geometry figures, each drawing onto a caller-supplied ``ax``
  so ``report.py`` can compose and save them as PNGs. Blade section, cascade,
  thickness and curvature wrap ParaBlade's *native* analytic plotters (matching
  the notebook); passage geometry, passage-width and the 3-blade
  geometry-verification reproduce the notebook's custom figures.

matplotlib is imported lazily inside the plot functions, so the pure-geometry
math (``generate_airfoil`` / ``build_blade`` / ``passage_metrics``) still runs on
a machine with no plotting stack.
"""
from __future__ import annotations
from typing import Dict, Any, Optional
import os
import sys
import json
import numpy as np

# Vendored ParaBlade lives at <repo>/parablade-master/parablade; derive that from
# this file's location (turbine_mcp/turbine_design/geometry.py) so clones work.
DEFAULT_PARABLADE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "parablade-master", "parablade")


# --------------------------------------------------------------------------- #
# ParaBlade construction
# --------------------------------------------------------------------------- #
def _design_variables(parablade_design: Dict[str, Any]) -> Dict[str, np.ndarray]:
    """Map a bladedesign ``parablade`` sub-dict into ParaBlade's design-variable
    dict (singleton numpy arrays, with the thickness control points expanded)."""
    dv: Dict[str, np.ndarray] = {}
    for key in ("stagger", "theta_in", "theta_out", "radius_in", "radius_out",
                "dist_in", "dist_out"):
        dv[key] = np.asarray(parablade_design[key])
    for i, t in enumerate(parablade_design["thickness_upper"], start=1):
        dv[f"thickness_upper_{i}"] = np.asarray(t)
    for i, t in enumerate(parablade_design["thickness_lower"], start=1):
        dv[f"thickness_lower_{i}"] = np.asarray(t)
    return dv


def build_blade(parablade_design: Dict[str, Any], n_section: int = 1000,
                s_over_Cx: Optional[float] = None,
                parablade_path: str = DEFAULT_PARABLADE):
    """Construct and initialise a ParaBlade ``Blade2DCamberThickness``.

    Args:
        parablade_design: the 'parablade' sub-dict from ``design_blade``.
        n_section: number of section sample points used to initialise the blade.
        s_over_Cx: non-dimensional pitch (pitch / axial chord); when given it is
            set as ``blade.spacing`` so the native cascade plot uses the real
            pitch instead of ParaBlade's 0.75*chord default.
        parablade_path: path to the parablade package.

    Returns:
        the initialised ParaBlade blade object.
    """
    if parablade_path not in sys.path:
        sys.path.append(parablade_path)
    from parablade.blade_2D_camber_thickness import Blade2DCamberThickness

    blade = Blade2DCamberThickness(_design_variables(parablade_design))
    blade.get_section_coordinates(np.linspace(0.0, 1.0, n_section))
    if s_over_Cx is not None:
        blade.spacing = float(s_over_Cx)
    return blade


# --------------------------------------------------------------------------- #
# .dat + design sidecar
# --------------------------------------------------------------------------- #
def _sidecar_path(dat_path: str) -> str:
    base, _ = os.path.splitext(dat_path)
    return base + ".parablade.json"


def save_parablade(parablade_design: Dict[str, Any], dat_path: str,
                   extra: Optional[Dict[str, Any]] = None) -> str:
    """Write the ParaBlade design variables next to ``dat_path`` as JSON so the
    analytic blade can be rebuilt later from just the ``.dat`` path."""
    def _coerce(v):
        if isinstance(v, np.ndarray):
            return v.tolist()
        if isinstance(v, (list, tuple)):
            return [_coerce(x) for x in v]
        if isinstance(v, (np.floating, np.integer)):
            return float(v)
        return v
    payload = {"parablade": {k: _coerce(v) for k, v in parablade_design.items()}}
    if extra:
        payload.update({k: _coerce(v) for k, v in extra.items()})
    out = _sidecar_path(dat_path)
    with open(out, "w") as f:
        json.dump(payload, f, indent=2)
    return out


def load_parablade(dat_path: str) -> Optional[Dict[str, Any]]:
    """Load the design sidecar written by ``save_parablade`` (or None if absent).
    Returns the full payload; the design dict is under key ``'parablade'``."""
    path = _sidecar_path(dat_path)
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def generate_airfoil(parablade_design: Dict[str, Any], axial_chord_m: float,
                     out_dat: Optional[str] = None, n_points: int = 600,
                     pitch_m: Optional[float] = None,
                     parablade_path: str = DEFAULT_PARABLADE) -> Dict[str, Any]:
    """Build the 2D airfoil with ParaBlade and return its surface coordinates.

    Args:
        parablade_design: the 'parablade' sub-dict from design_blade().
        axial_chord_m: dimensional axial chord, used to scale the coordinates.
        out_dat: optional path to write the dimensional coordinates (.dat). When
            given, a ``<dat>.parablade.json`` design sidecar is written alongside.
        n_points: number of points per surface.
        pitch_m: optional blade pitch [m]; stored in the sidecar (s/Cx) for the
            geometry/cascade figures.
        parablade_path: path to the parablade package.

    Returns:
        dict with suction/pressure surface coordinates, LE/TE points, and the
        output paths (.dat and design sidecar).
    """
    blade = build_blade(parablade_design, n_section=1000, parablade_path=parablade_path)

    u = np.linspace(0.0, 1.0, n_points)
    us = np.real(blade.get_upper_side_coordinates(u))   # suction, u=0->TE,1->LE
    ls = np.real(blade.get_lower_side_coordinates(u))   # pressure, u=0->LE,1->TE
    Cx = axial_chord_m
    suction = np.column_stack([us[0] * Cx, us[1] * Cx])
    pressure = np.column_stack([ls[0] * Cx, ls[1] * Cx])

    LE = (float(suction[-1, 0]), float(suction[-1, 1]))
    TE = (float(suction[0, 0]), float(suction[0, 1]))

    parablade_json = None
    if out_dat:
        os.makedirs(os.path.dirname(os.path.abspath(out_dat)), exist_ok=True)
        # pressure (LE->TE) then suction (LE->TE), matching the notebook layout
        pres_le2te = pressure
        suc_le2te = suction[::-1]
        arr = np.vstack([pres_le2te, suc_le2te])
        np.savetxt(out_dat, arr, header="x[m]  y[m]", fmt="%.10e")
        extra = {"axial_chord_m": float(Cx)}
        if pitch_m is not None:
            extra["pitch_m"] = float(pitch_m)
            extra["s_over_Cx"] = float(pitch_m / Cx)
        parablade_json = save_parablade(parablade_design, out_dat, extra=extra)

    return {
        "axial_chord_m": Cx,
        "LE": LE, "TE": TE,
        "n_points": n_points,
        "suction_xy": suction.tolist(),
        "pressure_xy": pressure.tolist(),
        "out_dat": out_dat,
        "parablade_json": parablade_json,
    }


# --------------------------------------------------------------------------- #
# Analytic passage / throat
# --------------------------------------------------------------------------- #
def passage_metrics(blade, s_over_Cx: float, n: int = 1200) -> Dict[str, Any]:
    """Non-dimensional passage-gap distribution and geometric throat.

    Builds the channel between one blade's suction surface and the neighbour's
    pressure surface (shifted by the pitch ``s_over_Cx``) and finds the minimum
    gap (throat). All lengths are normalised by the axial chord ``Cx``.
    """
    u = np.linspace(0.0, 1.0, 2000)
    us = np.real(blade.get_upper_side_coordinates(u))   # suction  (2, N)
    ls = np.real(blade.get_lower_side_coordinates(u))   # pressure (2, N)
    cl = np.real(blade.get_camberline_coordinates(u))   # camber   (2, N)

    us_s = us[:, np.argsort(us[0, :])]
    ls_s = ls[:, np.argsort(ls[0, :])]
    x_lo = max(us_s[0, 0], ls_s[0, 0])
    x_hi = min(us_s[0, -1], ls_s[0, -1])
    xg = np.linspace(x_lo, x_hi, n)
    y_suc = np.interp(xg, us_s[0, :], us_s[1, :])
    y_pres = np.interp(xg, ls_s[0, :], ls_s[1, :]) + s_over_Cx
    gap = y_pres - y_suc

    i_thr = int(np.argmin(gap))
    o_over_Cx = float(gap[i_thr])
    return {
        "x": xg, "gap": gap, "y_suc": y_suc, "y_pres": y_pres,
        "camber_xy": cl,
        "x_throat": float(xg[i_thr]),
        "o_over_Cx": o_over_Cx,
        "o_over_s": float(o_over_Cx / s_over_Cx) if s_over_Cx else float("nan"),
        "s_over_Cx": float(s_over_Cx),
    }


# --------------------------------------------------------------------------- #
# Plot functions (lazy matplotlib; each draws onto a supplied ax)
# --------------------------------------------------------------------------- #
def _ensure_tick_label_compat() -> None:
    """ParaBlade's ``plot_*`` methods call ``tick.label.set_fontsize(...)``, but
    matplotlib removed ``Tick.label`` in 3.8 (split into ``label1``/``label2``).
    Restore the attribute so the unmodified upstream ParaBlade (pulled as a git
    submodule) works on modern matplotlib without patching its source."""
    import matplotlib.axis as _maxis
    if not hasattr(_maxis.Tick, "label"):
        _maxis.Tick.label = property(lambda self: self.label1)


def plot_blade_section(ax, blade) -> None:
    """Single rotor blade section (surfaces + camber line) via ParaBlade."""
    _ensure_tick_label_compat()
    blade.plot_blade_section(fig=ax.figure, ax=ax,
                             upper_side='yes', lower_side='yes',
                             upper_side_control_points='no', lower_side_control_points='no',
                             camberline='yes')
    ax.set_title("Rotor blade section")


def plot_blade_cascade(ax, blade) -> None:
    """Three-blade cascade via ParaBlade (uses blade.spacing = s/Cx)."""
    _ensure_tick_label_compat()
    blade.plot_blade_cascade(fig=ax.figure, ax=ax)
    ax.set_title("Rotor blade cascade")


def plot_thickness_distribution(ax, blade) -> None:
    """Prescribed thickness distribution (upper/lower B-splines) via ParaBlade."""
    _ensure_tick_label_compat()
    blade.plot_thickness_distribution(fig=ax.figure, ax=ax)
    ax.set_title("Thickness distribution")


def plot_curvature_distribution(ax, blade) -> None:
    """Analytic section curvature distribution via ParaBlade."""
    _ensure_tick_label_compat()
    blade.plot_curvature_distribution(fig=ax.figure, ax=ax)
    ax.set_title("Curvature distribution")


def plot_passage_geometry(ax, blade, s_over_Cx: float) -> Dict[str, Any]:
    """Rotor blade passage geometry: suction + neighbour pressure surface, the
    shaded channel, camber line at mid-passage, and the geometric throat
    (reproduces the notebook 'Rotor Blade Passage Geometry' figure)."""
    m = passage_metrics(blade, s_over_Cx)
    u = np.linspace(0.0, 1.0, 2000)
    us = np.real(blade.get_upper_side_coordinates(u))
    ls = np.real(blade.get_lower_side_coordinates(u))
    cl = m["camber_xy"]
    xg, y_suc, y_pres = m["x"], m["y_suc"], m["y_pres"]
    x_thr, g_thr = m["x_throat"], m["o_over_Cx"]
    i_thr = int(np.argmin(m["gap"]))
    y_tmid = 0.5 * (y_suc[i_thr] + y_pres[i_thr])

    ax.plot(us[0, :], us[1, :], color="steelblue", lw=2.2, label="Blade 1 — suction surface")
    ax.plot(ls[0, :], ls[1, :] + s_over_Cx, color="crimson", lw=2.2, label="Blade 2 — pressure surface")
    ax.fill_between(xg, y_suc, y_pres, alpha=0.14, color="limegreen", label="Passage channel")
    ax.plot(cl[0, :], cl[1, :] + s_over_Cx / 2.0, "--", color="dimgray", lw=1.2,
            label="Camberline (midpassage)")
    ax.annotate("", xy=(x_thr, y_pres[i_thr]), xytext=(x_thr, y_suc[i_thr]),
                arrowprops=dict(arrowstyle="<|-|>", color="darkgreen", lw=1.8, mutation_scale=14))
    ax.text(x_thr + 0.03, y_tmid, f"Throat\n$o/C_x$ = {g_thr:.3f}", color="darkgreen",
            fontsize=9, va="center")
    ax.set_xlabel(r"$x / C_x$  (axial)")
    ax.set_ylabel(r"$y / C_x$  (pitchwise)")
    ax.set_title("Rotor blade passage geometry")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    return m


def plot_passage_width(ax, blade, s_over_Cx: float) -> Dict[str, Any]:
    """Passage gap g/Cx along the axial chord with the geometric throat marked
    (reproduces the notebook 'Passage Width Variation' figure)."""
    m = passage_metrics(blade, s_over_Cx)
    xg, gap = m["x"], m["gap"]
    x_thr, g_thr = m["x_throat"], m["o_over_Cx"]
    ax.plot(xg, gap, color="royalblue", lw=2.2, label=r"Gap  $g/C_x$")
    ax.axhline(g_thr, color="red", ls="--", lw=1.4, label=f"Throat  $o/C_x$ = {g_thr:.4f}")
    ax.axvline(x_thr, color="darkgreen", ls=":", lw=1.2, label=f"Throat at $x/C_x$ = {x_thr:.3f}")
    ax.fill_between(xg, 0, gap, alpha=0.10, color="royalblue")
    ax.set_xlabel(r"$x / C_x$")
    ax.set_ylabel(r"Passage gap  $g / C_x$")
    ax.set_title("Passage width variation along axial chord")
    ax.legend(loc="lower center", fontsize=8)
    ax.grid(True, alpha=0.3)
    return m


def plot_geometry_verification(ax, blade, axial_chord_m: float, pitch_m: float,
                               beta2_deg: float, beta3_deg: float) -> None:
    """Three-blade cascade geometry verification (dimensional): middle blade,
    the two neighbour walls actually used as CFD boundaries, and the
    flow-aligned inlet/outlet planes (reproduces the notebook pre-mesh check)."""
    Cx, pitch = axial_chord_m, pitch_m
    u = np.linspace(0.0, 1.0, 600)
    us = np.real(blade.get_upper_side_coordinates(u))   # suction
    ls = np.real(blade.get_lower_side_coordinates(u))   # pressure
    xsu, ysu = us[0] * Cx, us[1] * Cx
    xsl, ysl = ls[0] * Cx, ls[1] * Cx
    LEx, TEx = float(xsu[-1]), float(xsu[0])            # u=1 -> LE, u=0 -> TE

    bx = np.r_[xsu, xsl]
    by = np.r_[ysu, ysl]

    uw = np.column_stack([xsl, ysl + pitch])            # upper wall = pressure(+pitch)
    lw = np.column_stack([xsu[::-1], ysu[::-1] - pitch])  # lower wall = suction(-pitch)

    xin, xout = LEx - 0.5 * Cx, TEx + 0.5 * Cx
    si = np.tan(np.radians(beta2_deg))                  # inlet flow angle (+beta2)
    so = -np.tan(np.radians(beta3_deg))                 # outlet flow angle (-beta3): downward
    ext_in = lambda p: np.array([xin, p[1] + si * (xin - p[0])])
    ext_out = lambda p: np.array([xout, p[1] + so * (xout - p[0])])
    uw_f = np.vstack([ext_in(uw[0]), uw, ext_out(uw[-1])])
    lw_f = np.vstack([ext_in(lw[0]), lw, ext_out(lw[-1])])

    ax.fill(bx, by, color="0.75", ec="k", lw=1.5, label="Middle blade (measured)", zorder=4)
    ax.fill(bx, by + pitch, color="0.92", ec="0.6", lw=0.8, zorder=2)
    ax.fill(bx, by - pitch, color="0.92", ec="0.6", lw=0.8, zorder=2)
    ax.plot(uw_f[:, 0], uw_f[:, 1], "r-", lw=2.5, label="WALL_UP = pressure side, upper blade")
    ax.plot(lw_f[:, 0], lw_f[:, 1], "b-", lw=2.5, label="WALL_LO = suction side, lower blade")
    ax.plot([xin, xin], [lw_f[0, 1], uw_f[0, 1]], "g-", lw=2, label="INLET")
    ax.plot([xout, xout], [uw_f[-1, 1], lw_f[-1, 1]], "m-", lw=2, label="OUTLET")
    ax.set_aspect("equal")
    ax.grid(alpha=0.3)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("Three-blade cascade: geometry verification")
    ax.legend(loc="upper right", fontsize=8)
