"""1D mean-line design of an axial turbine stage.

Pure-Python, no external binaries. Given the imposed boundary conditions and a
set of design coefficients (stage loading, flow coefficient, degree of reaction,
total-to-total efficiency), the module returns the full thermodynamic and
kinematic state at the three stations, the velocity triangles, the annulus
geometry and a design-constraint check.

All angles are in degrees, pressures in Pa, temperatures in K, velocities in m/s.
"""
from __future__ import annotations
from dataclasses import dataclass, field, asdict
from typing import Dict, Any
import math


@dataclass
class StageInputs:
    """Imposed boundary conditions (defaults = the LMECA2323 HP-stage brief)."""
    gamma: float = 4.0 / 3.0          # ratio of specific heats (combustion gas)
    R: float = 287.0                  # gas constant [J/kg/K]
    mdot: float = 750.0               # turbine mass flow [kg/s]
    far: float = 0.025                # fuel-air ratio
    P_shaft: float = 330e6            # net shaft (extracted) power [W]
    pi_C: float = 20.0                # compressor pressure ratio
    p_adm: float = 101325.0           # compressor inlet static pressure [Pa]
    t_adm_C: float = 15.0             # compressor inlet static temperature [C]
    eta_isC: float = 0.92             # compressor isentropic efficiency
    eta_m: float = 0.95               # mechanical efficiency
    M_comb: float = 0.15              # combustor-exit Mach number
    TIT_C: float = 1400.0             # turbine inlet static temperature [C]
    N_rpm: float = 3000.0             # shaft speed [rpm]
    U_limit: float = 350.0            # blade peripheral-speed limit [m/s]
    eta_row: float = 0.90             # fixed stator/rotor blade-row efficiency
    n_stages: int = 4                 # stages sharing the total power
    gamma_air: float = 1.4            # compressor (air) ratio of specific heats


@dataclass
class DesignCoeffs:
    """Free design coefficients (the variables of the iteration)."""
    psi: float = 1.8                  # stage loading  dH_real / U^2
    phi: float = 0.5                  # flow coefficient  Vx / U
    DOR: float = 0.35                 # degree of reaction
    eta_TT: float = 0.87              # assumed total-to-total efficiency
    M3: float = 0.30                  # target absolute exit Mach number


def _isen_p_from_T(p0: float, T0: float, T: float, g: float) -> float:
    return p0 * (T / T0) ** (g / (g - 1.0))


def run_meanline(inp: StageInputs, cf: DesignCoeffs,
                 tol: float = 1.0, max_iter: int = 100) -> Dict[str, Any]:
    """Run the full 1D mean-line design and return all results as a dict.

    The rotor exit static pressure P3 is the unknown that makes the actual stage
    work equal to the required specific work; it is found by a bisection/secant
    loop on the enthalpy residual.
    """
    g, R = inp.gamma, inp.R
    cp = g * R / (g - 1.0)
    cv = R / (g - 1.0)
    T_adm = inp.t_adm_C + 273.15
    T1 = inp.TIT_C + 273.15
    omega = 2.0 * math.pi * inp.N_rpm / 60.0

    # --- compressor power & inlet conditions ---
    k = inp.gamma_air
    mdot_C = inp.mdot / (1.0 + inp.far)
    P_C = (mdot_C / (inp.eta_isC * inp.eta_m)) * (k / (k - 1.0)) * R * T_adm \
        * (inp.pi_C ** ((k - 1.0) / k) - 1.0)
    P1 = inp.pi_C * inp.p_adm
    P01 = P1 * (1.0 + 0.5 * (g - 1.0) * inp.M_comb ** 2) ** (g / (g - 1.0))
    T01 = T1 * (1.0 + 0.5 * (g - 1.0) * inp.M_comb ** 2)

    # --- required specific work for one stage ---
    dH_real = (P_C + inp.P_shaft) / inp.n_stages / inp.mdot

    # --- blade speed from stage loading (capped at mechanical limit) ---
    U_uncapped = math.sqrt(dH_real / cf.psi)
    U = min(U_uncapped, inp.U_limit)
    U_capped = U_uncapped > inp.U_limit
    Vx = cf.phi * U

    # --- exit-state first guess (to bracket P3) ---
    dH_is = dH_real / inp.eta_row
    T03is = T01 - dH_is / cp
    P03_guess = _isen_p_from_T(P01, T01, T03is, g)

    def stage_from_P3(P3: float):
        # ----- station 2 (stator exit / rotor inlet) -----
        P2 = cf.DOR * (P1 - P3) + P3
        V2is = math.sqrt(max(2.0 * cp * T01 * (1.0 - (P2 / P1) ** ((g - 1.0) / g)), 0.0))
        V2 = math.sqrt(inp.eta_row * V2is ** 2)
        T2 = T01 - V2 ** 2 / (2.0 * cp)
        rho2 = P2 / (R * T2)
        a2 = math.sqrt(g * R * T2)
        M2 = V2 / a2
        T02 = T01
        P02 = _isen_p_from_T(P2, T2, T02, g)
        Vx2 = Vx
        Vu2 = math.sqrt(max(V2 ** 2 - Vx2 ** 2, 0.0))
        alpha2 = math.degrees(math.acos(min(Vx2 / V2, 1.0)))
        Wu2 = Vu2 - U
        W2 = math.hypot(Wu2, Vx2)
        beta2 = math.degrees(math.atan2(abs(Wu2), Vx2))
        Mw2 = W2 / a2
        T02R = T2 + W2 ** 2 / (2.0 * cp)
        P02R = _isen_p_from_T(P2, T2, T02R, g)

        # ----- station 3 (rotor exit) -----
        W3is = math.sqrt(max(2.0 * cp * T02R * (1.0 - (P3 / P02R) ** ((g - 1.0) / g)), 0.0))
        W3 = math.sqrt(inp.eta_row * W3is ** 2)
        T3 = T02R - W3 ** 2 / (2.0 * cp)
        rho3 = P3 / (R * T3)
        a3 = math.sqrt(g * R * T3)
        Mw3 = W3 / a3
        Vx3 = Vx2                                  # constant axial velocity
        Wu3 = math.sqrt(max(W3 ** 2 - Vx3 ** 2, 0.0))
        beta3 = math.degrees(math.atan2(Wu3, Vx3))
        Vu3 = Wu3 - U
        V3 = math.hypot(Vu3, Vx3)
        alpha3 = math.degrees(math.atan2(abs(Vu3), Vx3))
        M3 = V3 / a3
        T03 = T3 + V3 ** 2 / (2.0 * cp)
        P03 = _isen_p_from_T(P3, T3, T03, g)
        dH0 = cp * (T01 - T03)
        return locals()

    # --- secant iteration on P3 so that dH0 == dH_real ---
    lo, hi = 0.40 * P01, 0.98 * P01
    f_lo = stage_from_P3(lo)["dH0"] - dH_real
    f_hi = stage_from_P3(hi)["dH0"] - dH_real
    P3 = P03_guess
    work_converged = False
    for _ in range(max_iter):
        s = stage_from_P3(P3)
        f = s["dH0"] - dH_real
        if abs(f) < tol:
            work_converged = True
            break
        # keep a valid bracket and bisect/secant
        if f_lo * f < 0:
            hi, f_hi = P3, f
        else:
            lo, f_lo = P3, f
        denom = (f_hi - f_lo)
        P3 = (lo + hi) / 2.0 if denom == 0 else hi - f_hi * (hi - lo) / denom
        P3 = min(max(P3, lo), hi)
    s = stage_from_P3(P3)

    # --- annulus sizing (constant mean radius) ---
    R_mean = U / omega
    D_m = 2.0 * R_mean
    A2 = inp.mdot / (s["rho2"] * s["Vx2"])
    h2 = A2 / (math.pi * D_m)
    A3 = inp.mdot / (s["rho3"] * s["Vx3"])
    h3 = A3 / (math.pi * D_m)
    AN2 = 0.5 * (A2 + A3) * inp.N_rpm ** 2
    turning = s["beta2"] + s["beta3"]

    # --- feasibility / solver health (so an infeasible design can't pass silently) ---
    residual = s["dH0"] - dH_real
    psi_eff = dH_real / U ** 2
    warnings: list = []
    if U_capped:
        warnings.append(
            f"Blade speed clamped to the {inp.U_limit:.0f} m/s limit: loading "
            f"psi={cf.psi:.2f} would need U={U_uncapped:.0f} m/s, so the EFFECTIVE "
            f"loading is psi={psi_eff:.2f}, not {cf.psi:.2f}.")
    if not work_converged:
        warnings.append(
            f"Work-match did NOT converge: the stage delivers ~{s['dH0']/1e3:.1f} "
            f"kJ/kg but {dH_real/1e3:.1f} kJ/kg is required "
            f"(residual {residual/1e3:.1f} kJ/kg). The design is infeasible at these "
            "inputs - add stages / lower psi, or raise U_limit. The stations, angles "
            "and constraints below are from this unconverged state and are unreliable.")

    result: Dict[str, Any] = {
        "inputs": asdict(inp),
        "coefficients": asdict(cf),
        "warnings": warnings,
        "feasibility": {
            "work_converged": work_converged,
            "U_capped": U_capped,
            "U_uncapped_mps": U_uncapped,
            "psi_effective": psi_eff,
            "residual_Jkg": residual,
        },
        "derived": {
            "cp": cp, "cv": cv, "compressor_power_W": P_C,
            "specific_work_Jkg": dH_real, "delivered_work_Jkg": s["dH0"],
            "U_mps": U, "Vx_mps": Vx,
            "P3_Pa": P3, "iteration_residual_Jkg": residual,
        },
        "station1": {
            "P_Pa": P1, "P0_Pa": P01, "T_K": T1, "T0_K": T01,
            "V_mps": Vx, "Vx_mps": Vx, "alpha_deg": inp.__dict__.get("alpha1", 0.0),
            "M": Vx / math.sqrt(g * R * T1),
        },
        "station2": {
            "P_Pa": s["P2"], "P0_Pa": s["P02"], "P0rel_Pa": s["P02R"],
            "T_K": s["T2"], "T0_K": s["T02"], "T0rel_K": s["T02R"],
            "V_mps": s["V2"], "Vx_mps": s["Vx2"], "Vu_mps": s["Vu2"],
            "W_mps": s["W2"], "alpha_deg": s["alpha2"], "beta_deg": s["beta2"],
            "M": s["M2"], "Mrel": s["Mw2"], "rho": s["rho2"],
        },
        "station3": {
            "P_Pa": P3, "P0_Pa": s["P03"], "P0rel_Pa": s["P02R"],
            "T_K": s["T3"], "T0_K": s["T03"], "T0rel_K": s["T02R"],
            "V_mps": s["V3"], "Vx_mps": s["Vx3"], "Vu_mps": s["Vu3"],
            "W_mps": s["W3"], "alpha_deg": s["alpha3"], "beta_deg": s["beta3"],
            "M": s["M3"], "Mrel": s["Mw3"], "rho": s["rho3"], "turning_deg": turning,
        },
        "annulus": {
            "R_mean_m": R_mean, "D_mean_m": D_m,
            "A2_m2": A2, "h2_m": h2, "Rtip2_m": R_mean + h2 / 2, "Rhub2_m": R_mean - h2 / 2,
            "A3_m2": A3, "h3_m": h3, "Rtip3_m": R_mean + h3 / 2, "Rhub3_m": R_mean - h3 / 2,
            "AN2_m2rpm2": AN2, "span_ratio_h3_h2": h3 / h2,
        },
    }
    result["constraints"] = _check_constraints(result, inp)
    return result


def _check_constraints(res: Dict[str, Any], inp: StageInputs) -> Dict[str, Any]:
    s2, s3, ann = res["station2"], res["station3"], res["annulus"]
    checks = {
        "M2": (s2["M"], 0.70, 0.85),
        "beta2_deg": (s2["beta_deg"], None, 47.5),
        "Mw2": (s2["Mrel"], None, 0.50),
        "Mw3": (s3["Mrel"], 0.65, 0.80),
        "beta3_deg": (s3["beta_deg"], None, 75.0),
        "turning_deg": (s3["turning_deg"], 110.0, 120.0),
        "alpha3_deg": (s3["alpha_deg"], None, 35.0),
        "h3_h2": (ann["span_ratio_h3_h2"], None, 1.20),
        "AN2_m2rpm2": (ann["AN2_m2rpm2"], None, 3e7),
    }
    out = {}
    all_ok = True
    for name, (val, lo, hi) in checks.items():
        ok = (lo is None or val >= lo) and (hi is None or val <= hi)
        all_ok = all_ok and ok
        out[name] = {"value": val, "min": lo, "max": hi, "pass": ok}
    out["all_pass"] = all_ok
    return out


# =========================================================================== #
#  Mean-line figures (step 1: meanline.py owns its own plots)
#
#  Each draws onto a caller-supplied ``ax`` so ``report.py`` can compose and save
#  them as PNGs. matplotlib is imported lazily so the design math above still runs
#  on a machine with no plotting stack. All consume the dict from run_meanline.
# =========================================================================== #
_CONSTRAINT_LABELS = {
    "M2": "M2 (stator-exit Mach)", "beta2_deg": "β2 [°]", "Mw2": "Mw2 (rel.)",
    "Mw3": "Mw3 (rel.)", "beta3_deg": "β3 [°]", "turning_deg": "turning Δβ [°]",
    "alpha3_deg": "α3 [°]", "h3_h2": "span ratio h3/h2", "AN2_m2rpm2": "AN²",
}


def enthalpy_ladder(result: Dict[str, Any]) -> Dict[str, float]:
    """Enthalpy-entropy ladder of the stage (faithful to the notebook h-s cell).

    Returns the static / total / isentropic enthalpy levels [J/kg] and the static
    and total entropies [J/kg/K] referenced to the station-1 static state. Used by
    both ``hs_diagram`` and the LaTeX report's enthalpy-ladder table.
    """
    import numpy as np
    g = result["inputs"]["gamma"]; R = result["inputs"]["R"]; cp = result["derived"]["cp"]
    s1c, s2c, s3c = result["station1"], result["station2"], result["station3"]
    T1, P1, P01 = s1c["T_K"], s1c["P_Pa"], s1c["P0_Pa"]
    T2, P2, P02 = s2c["T_K"], s2c["P_Pa"], s2c["P0_Pa"]
    T02R, P02R = s2c["T0rel_K"], s2c["P0rel_Pa"]
    T3, P3, T03, P03 = s3c["T_K"], s3c["P_Pa"], s3c["T0_K"], s3c["P0_Pa"]
    V3, W3, T01 = s3c["V_mps"], s3c["W_mps"], s1c["T0_K"]

    h1 = cp * T1; H01 = cp * T01; H02 = H01
    h2 = cp * T2; H02R = cp * T02R
    h2s = cp * (T01 * (P2 / P1) ** ((g - 1) / g))
    # rotor isentropic exit (relative frame), then back to static
    T3s_rel = T02R * (P3 / P02R) ** ((g - 1) / g)
    W3s = float(np.sqrt(max(2 * cp * (T02R - T3s_rel), 0.0)))
    h3s = cp * (T02R - W3s ** 2 / (2 * cp))
    h3 = cp * T3; H03 = cp * T03
    h3ss = cp * (T01 * (P3 / P01) ** ((g - 1) / g))
    H03ss = h3ss + V3 ** 2 / 2.0
    H03R = h3 + 0.5 * W3 ** 2

    def ent(T, P):
        return cp * np.log(T / T1) - R * np.log(P / P1)

    return {
        "H01": H01, "H02": H02, "h1": h1, "h2": h2, "h2s": h2s, "H0R": H02R,
        "h3": h3, "h3s": h3s, "h3ss": h3ss, "H03": H03, "H03ss": H03ss, "H03R": H03R,
        "s1": 0.0, "s2": float(ent(T2, P2)), "s3": float(ent(T3, P3)),
        "s01": float(ent(T01, P01)), "s02": float(ent(T01, P02)), "s03": float(ent(T03, P03)),
    }


def velocity_triangles(ax, result: Dict[str, Any]) -> None:
    """Rotor inlet (2) and exit (3) velocity triangles — filled triangles with
    vector arrows and α/β angle arcs (notebook rendering). x = tangential,
    y = axial direction."""
    import numpy as np
    from matplotlib.patches import Arc, Polygon, FancyArrowPatch
    d, s2, s3 = result["derived"], result["station2"], result["station3"]
    U = d["U_mps"]
    W2, b2, V2 = s2["W_mps"], s2["beta_deg"], s2["V_mps"]
    W3, b3, V3 = s3["W_mps"], s3["beta_deg"], s3["V_mps"]
    Vx = s2["Vx_mps"]

    A = np.array([0.0, 0.0]); C = np.array([U, 0.0])
    B = np.array([-W2 * np.sin(np.radians(b2)), W2 * np.cos(np.radians(b2))])  # inlet apex
    D = np.array([W3 * np.sin(np.radians(b3)), W3 * np.cos(np.radians(b3))])   # exit apex

    cU, cV2, cW2, cV3, cW3 = "#2b2b2b", "#0571b0", "#c1272d", "#2f9e44", "#e8893a"
    ax.add_patch(Polygon([A, B, C], closed=True, facecolor=cV2, alpha=0.07, edgecolor="none"))
    ax.add_patch(Polygon([A, D, C], closed=True, facecolor=cW3, alpha=0.07, edgecolor="none"))
    ax.plot([B[0], D[0]], [Vx, Vx], ls=(0, (4, 4)), color="#9aa0a6", lw=1.1, zorder=1)
    ax.annotate(f"$V_x = {Vx:.0f}$ m/s", xy=((B[0] + D[0]) / 2, Vx), xytext=(0, 7),
                textcoords="offset points", ha="center", color="#5f6368", fontsize=9)

    def vec(p0, p1, color, label, mag, loff):
        ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=18, lw=2.2,
                                     color=color, shrinkA=0, shrinkB=0, zorder=4))
        m = 0.52 * np.array(p0) + 0.48 * np.array(p1)
        ax.annotate(f"{label}\n{mag:.0f} m/s", xy=m, xytext=loff, textcoords="offset points",
                    color=color, fontsize=10, ha="center", va="center", fontweight="bold")

    vec(A, C, cU, r"$\vec{U}$", U, (0, -18))
    vec(B, A, cW2, r"$\vec{W}_2$", W2, (-24, 6))
    vec(B, C, cV2, r"$\vec{V}_2$", V2, (18, 12))
    vec(D, A, cW3, r"$\vec{W}_3$", W3, (26, 6))
    vec(D, C, cV3, r"$\vec{V}_3$", V3, (24, 6))

    def arc_to(apex, target, color, label, Rr, lr=1.32):
        dd = np.array(target) - np.array(apex)
        th = np.degrees(np.arctan2(dd[1], dd[0])) % 360
        t1, t2 = min(th, 270.0), max(th, 270.0)
        ax.add_patch(Arc(apex, 2 * Rr, 2 * Rr, theta1=t1, theta2=t2, color=color, lw=1.6, zorder=5))
        mid = np.radians((t1 + t2) / 2)
        ax.annotate(label, xy=(apex[0] + lr * Rr * np.cos(mid), apex[1] + lr * Rr * np.sin(mid)),
                    color=color, fontsize=11, ha="center", va="center", fontweight="bold")

    arc_to(B, C, cV2, r"$\alpha_2$", 58); arc_to(B, A, cW2, r"$\beta_2$", 34)
    arc_to(D, C, cV3, r"$\alpha_3$", 58); arc_to(D, A, cW3, r"$\beta_3$", 34)
    ax.plot(*B, "o", color="#444", ms=4, zorder=6); ax.plot(*D, "o", color="#444", ms=4, zorder=6)
    ax.annotate("Rotor inlet (2)", xy=B, xytext=(8, 12), textcoords="offset points",
                ha="left", fontsize=9, color="#444", style="italic")
    ax.annotate("Rotor exit (3)", xy=D, xytext=(-8, 12), textcoords="offset points",
                ha="right", fontsize=9, color="#444", style="italic")
    xs = [A[0], B[0], C[0], D[0]]; ys = [A[1], B[1], C[1], D[1]]
    ax.set_xlim(min(xs) - 90, max(xs) + 120); ax.set_ylim(min(ys) - 60, max(ys) + 70)
    ax.set_aspect("equal")
    ax.set_xlabel("Tangential direction  [m/s]"); ax.set_ylabel("Axial direction  [m/s]")
    ax.set_title("Velocity Triangles — Rotor Inlet (2) and Exit (3)")
    ax.grid(True, color="#e8eaed", lw=0.8)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def annulus_sketch(ax, result: Dict[str, Any]) -> None:
    """Meridional annulus from station 2 to 3 (constant mean radius)."""
    a = result["annulus"]
    x = [0.0, 1.0]
    ax.fill_between(x, [a["Rhub2_m"], a["Rhub3_m"]], [a["Rtip2_m"], a["Rtip3_m"]],
                    color="#cfe0f3", alpha=0.7)
    ax.plot(x, [a["Rtip2_m"], a["Rtip3_m"]], "k-", lw=1.5)
    ax.plot(x, [a["Rhub2_m"], a["Rhub3_m"]], "k-", lw=1.5)
    ax.axhline(a["R_mean_m"], color="0.5", ls="--", lw=1, label="mean radius")
    ax.set_xticks(x)
    ax.set_xticklabels(["station 2", "station 3"])
    ax.set_ylabel("radius [m]")
    ax.set_title(f"Annulus  (h2={a['h2_m']*1e3:.0f} mm → h3={a['h3_m']*1e3:.0f} mm)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.2)


def constraint_panel(ax, result: Dict[str, Any]) -> None:
    """Each design constraint plotted against its allowed band (green=pass)."""
    cons = {k: v for k, v in result["constraints"].items() if isinstance(v, dict)}
    names = list(cons.keys())
    for i, n in enumerate(names):
        rec = cons[n]
        val, lo, hi, ok = rec["value"], rec["min"], rec["max"], rec["pass"]
        if lo is not None and hi is not None:
            x = (val - lo) / (hi - lo); lim_txt = f"[{lo:g}, {hi:g}]"
        elif hi is not None:
            x = val / hi; lim_txt = f"≤ {hi:g}"
        elif lo is not None:
            x = lo / val if val else 0.0; lim_txt = f"≥ {lo:g}"
        else:
            x, lim_txt = 0.5, ""
        ax.barh(i, 1.0, left=0.0, color="#cfe8cf", height=0.55, zorder=1)
        ax.axvline(0.0, color="0.6", lw=0.8); ax.axvline(1.0, color="0.6", lw=0.8)
        ax.plot(min(max(x, -0.12), 1.12), i, "o", ms=9,
                color=("#2ca02c" if ok else "#d62728"), zorder=3)
        ax.text(1.18, i, f"{val:.3g}  ({lim_txt})", va="center", fontsize=7.5,
                color=("0.25" if ok else "#b22222"))
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels([_CONSTRAINT_LABELS.get(n, n) for n in names], fontsize=8)
    ax.invert_yaxis(); ax.set_xlim(-0.3, 1.7)
    ax.set_title("Constraint check  (green = pass, red = violated)")
    ax.set_xlabel("position within allowed band  (0 = lower, 1 = upper limit)")
    ax.grid(axis="x", alpha=0.2)


def hs_diagram(ax, result: Dict[str, Any]) -> None:
    """h–s ladder of the stage: static / total / isentropic states, left-labelled
    reference lines and the V²/2 and W²/2 kinetic-energy ladders (notebook style)."""
    from matplotlib.patches import FancyArrowPatch
    L = enthalpy_ladder(result)
    s1, s2, s3 = L["s1"], L["s2"], L["s3"]
    s01, s02, s03 = L["s01"], L["s02"], L["s03"]
    H01, H02, H0R = L["H01"], L["H02"], L["H0R"]
    h1, h2, h2s = L["h1"], L["h2"], L["h2s"]
    h3, h3s, h3ss = L["h3"], L["h3s"], L["h3ss"]
    H03, H03ss, H03R = L["H03"], L["H03ss"], L["H03R"]

    # process lines: total (01->02), static (1->2->3) and the per-station verticals
    ax.plot([s1, s2], [H01, H02], lw=2, color="#16a34a", zorder=3)
    ax.plot([s1, s2, s3], [h1, h2, h3], lw=2, color="#ea580c", zorder=3)
    ax.plot([s1, s1, s1, s1, s1], [H01, h1, h2s, H03ss, h3ss], lw=1.5, color="#16a34a", zorder=2)
    ax.plot([s2, s2, s2, s2], [H02, H0R, h2, h3s], lw=1.5, color="#ea580c", zorder=2)
    ax.plot([s3, s3, s3], [H0R, H03, h3], lw=1.5, color="#7c3aed", zorder=2)
    for s, h, t in [(s1, H01, "01"), (s2, H02, "02"), (s1, h1, "1"), (s2, h2, "2"),
                    (s3, h3, "3"), (s1, h2s, "2s"), (s2, h3s, "3s"),
                    (s1, h3ss, "3ss"), (s3, H03, "03")]:
        ax.scatter([s], [h], s=26, zorder=5, color="#333")
        ax.annotate(f" {t}", (s, h), fontsize=7.5)

    # left-labelled horizontal reference lines (enthalpy levels)
    s_label = min(s1, s2, s3) - 6.0
    x_right = max(s01, s02, s03, s3)
    for h, lab, ls in [(H01, r"$H_{01}=H_{02}$", "dashed"), (h1, r"$h_1$", "dotted"),
                       (h2, r"$h_2$", "dotted"), (h2s, r"$h_{2s}$", "dotted"),
                       (H0R, r"$H_{0R}$", "dashdot"), (H03, r"$H_{03}$", "dashed"),
                       (h3, r"$h_3$", "dotted"), (H03ss, r"$H_{03ss}$", "dotted"),
                       (h3s, r"$h_{3s}$", "dotted")]:
        ax.hlines(h, s_label, x_right, linestyles=ls, linewidth=0.9, color="0.45")
        ax.text(s_label, h, lab, va="center", ha="right", fontsize=8, style="italic")

    # kinetic-energy ladders (double-headed arrows)
    ds = 2.0
    astyle = dict(arrowstyle="<|-|>", lw=1.0, color="black", shrinkA=0, shrinkB=0)

    def ladder(sx, lo, hi, label, side):
        ax.add_patch(FancyArrowPatch((sx, lo), (sx, hi), **astyle))
        ax.annotate(label, xy=(sx + (0.5 if side == "r" else -0.5), 0.5 * (lo + hi)),
                    ha="left" if side == "r" else "right", va="center", fontsize=11)

    ladder(s1 + ds, h1, H01, r"$\frac{V_1^2}{2}$", "r")
    ladder(s2 + ds, h2, H02, r"$\frac{V_2^2}{2}$", "r")
    ladder(s2 - ds, h2, H0R, r"$\frac{W_2^2}{2}$", "l")
    ladder(s3 + ds, h3, H03, r"$\frac{V_3^2}{2}$", "r")
    ladder(s3 - ds, h3, H03R, r"$\frac{W_3^2}{2}$", "l")

    ax.set_xlim(s_label - 6, x_right + 6)
    ax.set_xlabel("Entropy  s − s$_1$  [J/kg·K]")
    ax.set_ylabel("Enthalpy  h, H  [J/kg]")
    ax.set_title("h–s Diagram for Turbine Stage")
    ax.grid(True, axis="x", alpha=0.2)
