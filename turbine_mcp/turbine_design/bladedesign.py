"""Preliminary rotor-blade design from the mean-line velocity triangles.

Converts the rotor inlet/exit angles and the annulus geometry into a blade-row
definition: chord, stagger, pitch, blade count, edge radii and a thickness
distribution suitable for ParaBlade.
"""
from __future__ import annotations
from dataclasses import dataclass, field, asdict
from typing import Dict, Any, List
import math


@dataclass
class BladeChoices:
    """Designer choices for the rotor blade."""
    AR: float = 1.0                   # aspect ratio  span / chord
    stagger_deg: float = 39.0         # stagger angle (Kacker-Okapuu correlation)
    zweifel: float = 1.0              # Zweifel loading coefficient
    LE_frac_pitch: float = 0.05       # leading-edge radius as fraction of pitch
    TE_radius_m: float = 0.0007       # trailing-edge radius [m]
    thickness_upper: List[float] = field(
        default_factory=lambda: [0.15, 0.20, 0.13, 0.07, 0.03, 0.02])
    thickness_lower: List[float] = field(
        default_factory=lambda: [0.15, 0.20, 0.13, 0.07, 0.03, 0.02])


def design_blade(beta2_deg: float, beta3_deg: float, alpha2_deg: float,
                 alpha3_deg: float, span_h2_m: float, D_mean_m: float,
                 ch: BladeChoices, span_h3_m: float = None) -> Dict[str, Any]:
    """Compute the rotor-blade geometric parameters.

    Args:
        beta2_deg, beta3_deg: rotor inlet/exit relative (metal) angles.
        alpha2_deg, alpha3_deg: rotor inlet/exit absolute flow angles.
        span_h2_m: blade span at rotor inlet (= annulus h2).
        D_mean_m: mean diameter.
        ch: blade design choices.
        span_h3_m: blade span at rotor exit (= annulus h3); if given, the end-wall
            flare angle ε = atan((h3 - h2) / cx) is reported.
    """
    b2 = math.radians(beta2_deg)
    b3 = math.radians(beta3_deg)
    xi = math.radians(ch.stagger_deg)

    chord = span_h2_m / ch.AR
    cx = chord * math.cos(xi)
    # optimum pitch from the (incompressible) Zweifel criterion
    pitch = ch.zweifel * cx / (2.0 * math.cos(b3) ** 2 * (math.tan(b2) + math.tan(b3)))
    n_blades = math.ceil(math.pi * D_mean_m / pitch)
    R_LE = ch.LE_frac_pitch * pitch
    t_max_c = max(max(ch.thickness_upper), max(ch.thickness_lower))

    out = {
        "choices": asdict(ch),
        "alpha2_deg": alpha2_deg, "beta2_deg": beta2_deg,
        "alpha3_deg": alpha3_deg, "beta3_deg": beta3_deg,
        "aspect_ratio": ch.AR,
        "stagger_deg": ch.stagger_deg,
        "zweifel": ch.zweifel,
        "chord_m": chord,
        "axial_chord_m": cx,
        "pitch_m": pitch,
        "pitch_to_axial_chord": pitch / cx,
        "n_blades": int(n_blades),
        "t_max_over_chord": t_max_c,
        "LE_radius_m": R_LE,
        "TE_radius_m": ch.TE_radius_m,
        # ParaBlade-ready (non-dimensional, normalised by axial chord)
        "parablade": {
            "stagger": -ch.stagger_deg,
            "theta_in": beta2_deg,
            "theta_out": -beta3_deg,
            "radius_in": R_LE / cx,
            "radius_out": ch.TE_radius_m / cx,
            "dist_in": math.cos(math.radians(beta2_deg)),
            "dist_out": math.sin(math.radians(beta2_deg)),
            "thickness_upper": ch.thickness_upper,
            "thickness_lower": ch.thickness_lower,
        },
    }

    if span_h3_m is not None:
        out["flare_angle_deg"] = math.degrees(math.atan((span_h3_m - span_h2_m) / cx))
    return out
