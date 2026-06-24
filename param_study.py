#!/usr/bin/env python3
"""Parametric study + case-log plots for the axial-turbine mean-line design.

Two things, re-run any time (they always reflect the *current* model + case-log):

  1. Coefficient sweeps -- hold the reference engine fixed and vary one design
     coefficient (degree of reaction DOR, flow coefficient phi, stage loading
     psi) to show how it ripples through the velocity triangles, Mach numbers,
     turning and annulus. Points are coloured by feasibility.

  2. Case-log overview -- read design_cases_log.xlsx (every turbine MCP design
     call appends a row) and scatter the logged cases, so the picture fills in
     as you run more test cases.

    python3 param_study.py

Figures are written to param_study/.
"""
import os
import sys

PROJ = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(PROJ, "turbine_mcp"))
os.environ.setdefault("PARABLADE_PATH", os.path.join(PROJ, "parablade-master", "parablade"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from turbine_design.meanline import StageInputs, DesignCoeffs, run_meanline

OUT = os.path.join(PROJ, "param_study")
os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(PROJ, "design_cases_log.xlsx")


# --------------------------------------------------------------------------- #
# Sweep machinery
# --------------------------------------------------------------------------- #
def _metrics(r):
    s2, s3, a, d, f = (r["station2"], r["station3"], r["annulus"],
                       r["derived"], r["feasibility"])
    return dict(
        U=d["U_mps"], Vx=d["Vx_mps"], work=d["specific_work_Jkg"] / 1e3,
        P3=d["P3_Pa"] / 1e3,
        M2=s2["M"], Mw2=s2["Mrel"], M3=s3["M"], Mw3=s3["Mrel"],
        beta2=s2["beta_deg"], beta3=s3["beta_deg"],
        alpha2=s2["alpha_deg"], alpha3=s3["alpha_deg"],
        V2=s2["V_mps"], V3=s3["V_mps"], W2=s2["W_mps"], W3=s3["W_mps"],
        P2=s2["P_Pa"] / 1e3, turning=s3["turning_deg"],
        h2=a["h2_m"] * 1e3, h3=a["h3_m"] * 1e3,
        AN2=a["AN2_m2rpm2"] / 1e6, span_ratio=a["span_ratio_h3_h2"],
        # physically realizable: the power balance converged and the blade
        # speed was not clamped to the mechanical limit
        feasible=bool(f["work_converged"] and not f["U_capped"]),
        constraints_pass=bool(r["constraints"]["all_pass"]),
    )


def sweep(param, values, **base):
    out = []
    for v in values:
        co = dict(psi=base.get("psi", 1.8), phi=base.get("phi", 0.5),
                  DOR=base.get("DOR", 0.35), eta_TT=base.get("eta_TT", 0.87),
                  M3=base.get("M3", 0.3))
        co[param] = float(v)
        m = _metrics(run_meanline(StageInputs(), DesignCoeffs(**co)))
        m[param] = float(v)
        out.append(m)
    return out


def plot_sweep(rows, xkey, xlabel, panels, title, fname):
    n = len(panels)
    ncol, nrow = 3, int(np.ceil(n / 3))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.4 * ncol, 3.2 * nrow))
    axes = np.atleast_1d(axes).flatten()
    x = [r[xkey] for r in rows]
    feas = [r["feasible"] for r in rows]
    colours = ["tab:green" if fe else "tab:red" for fe in feas]
    for ax, (key, lab) in zip(axes, panels):
        y = [r[key] for r in rows]
        ax.plot(x, y, "-", color="0.65", lw=1.3, zorder=1)
        ax.scatter(x, y, c=colours, s=42, zorder=2, edgecolors="k", linewidths=0.4)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(lab)
        ax.grid(alpha=0.3)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle(title, fontsize=14, fontweight="bold")
    fig.legend(handles=[Line2D([], [], marker="o", ls="", color="tab:green",
                               mec="k", label="physically feasible"),
                        Line2D([], [], marker="o", ls="", color="tab:red",
                               mec="k", label="infeasible (work not matched / U capped)")],
               loc="lower center", ncol=2, frameon=False)
    plt.tight_layout(rect=[0, 0.045, 1, 0.96])
    plt.savefig(os.path.join(OUT, fname), dpi=150)
    plt.close()
    print("wrote param_study/" + fname)


# --------------------------------------------------------------------------- #
# Case-log overview
# --------------------------------------------------------------------------- #
def plot_caselog():
    try:
        import openpyxl
    except ImportError:
        print("openpyxl not installed; skipping case-log plot")
        return
    if not os.path.exists(LOG):
        print("no case-log found; skipping")
        return
    ws = openpyxl.load_workbook(LOG, read_only=True).active
    rows = list(ws.iter_rows(values_only=True))
    hdr, data = rows[0], rows[1:]
    if not data:
        print("case-log empty; skipping")
        return
    idx = {h: i for i, h in enumerate(hdr)}

    def col(name):
        return [(row[idx[name]] if name in idx else None) for row in data]

    names = col("Case")
    # (x-column, y-column, x-label, y-label) cross-plots over the logged cases
    pairs = [("psi", "U [m/s]", r"loading $\psi$", "blade speed U [m/s]"),
             ("phi", "M2", r"flow coeff. $\phi$", "stator-exit $M_2$"),
             ("R", "turn [deg]", "reaction R", r"turning $\Delta\beta$ [deg]"),
             ("AN2", "h3/h2", r"$AN^2$", "span ratio $h_3/h_2$")]
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    for ax, (xc, yc, xl, yl) in zip(axes.flat, pairs):
        xs, ys, ls = [], [], []
        for nm, xv, yv in zip(names, col(xc), col(yc)):
            if isinstance(xv, (int, float)) and isinstance(yv, (int, float)):
                xs.append(xv); ys.append(yv); ls.append(nm)
        ax.scatter(xs, ys, s=70, c="tab:blue", edgecolors="k", zorder=3)
        for xv, yv, nm in zip(xs, ys, ls):
            ax.annotate(str(nm)[:14], (xv, yv), fontsize=7,
                        xytext=(4, 4), textcoords="offset points")
        ax.set_xlabel(xl); ax.set_ylabel(yl); ax.grid(alpha=0.3)
    fig.suptitle(f"Logged design cases  (n = {len(data)})  -- design_cases_log.xlsx",
                 fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(os.path.join(OUT, "caselog_overview.png"), dpi=150)
    plt.close()
    print(f"wrote param_study/caselog_overview.png  ({len(data)} cases)")


# --------------------------------------------------------------------------- #
def main():
    print("Sweeping coefficients about the reference engine...")
    plot_sweep(
        sweep("DOR", np.linspace(0.20, 0.50, 13)), "DOR", "degree of reaction $r$",
        [("M2", "stator-exit $M_2$"), ("Mw3", "rotor-exit $M_{w3}$"),
         ("P2", "stator-exit $P_2$ [kPa]"), ("beta2", r"$\beta_2$ [deg]"),
         ("V2", "$V_2$ [m/s]"), ("W3", "$W_3$ [m/s]")],
        "Effect of degree of reaction (DOR)", "effect_of_DOR.png")

    plot_sweep(
        sweep("phi", np.linspace(0.30, 0.70, 13)), "phi", r"flow coefficient $\phi$",
        [("Vx", "axial velocity $V_x$ [m/s]"), ("M2", "stator-exit $M_2$"),
         ("M3", "exit $M_3$"), ("turning", r"turning $\Delta\beta$ [deg]"),
         ("h3", "exit span $h_3$ [mm]"), ("AN2", r"$AN^2$ [$10^6$]")],
        "Effect of flow coefficient ($\\phi$)", "effect_of_phi.png")

    plot_sweep(
        sweep("psi", np.linspace(1.4, 2.2, 13)), "psi", r"stage loading $\psi$",
        [("U", "blade speed $U$ [m/s]"), ("turning", r"turning $\Delta\beta$ [deg]"),
         ("M2", "stator-exit $M_2$"), ("M3", "exit $M_3$"),
         ("h3", "exit span $h_3$ [mm]"), ("AN2", r"$AN^2$ [$10^6$]")],
        "Effect of stage loading ($\\psi$)", "effect_of_psi.png")

    print("Reading the case-log...")
    plot_caselog()
    print("Done -> param_study/")


if __name__ == "__main__":
    main()
