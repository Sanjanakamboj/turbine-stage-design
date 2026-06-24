#!/usr/bin/env python3
"""Run the full axial-turbine design chain end-to-end, no LLM / MCP client needed.

    mean-line  ->  rotor blade  ->  ParaBlade airfoil  ->  [GMSH + SU2 CFD]  ->  LaTeX report

This is the headless/automation entry point: it calls the same turbine_design
functions the MCP tools wrap, so a cron job, CI step, or shell loop can drive the
whole design without an interactive Claude session.

Examples
--------
    # mean-line + blade + airfoil + report (instant)
    python3 run_pipeline.py --case base

    # change the design coefficients
    python3 run_pipeline.py --psi 1.8 --phi 0.5 --dor 0.35 --case hp_design

    # also mesh + solve the cascade CFD (takes minutes; needs SU2_BIN)
    python3 run_pipeline.py --case hp_cfd --cfd

    # parameter sweep
    for p in 1.6 1.8 2.0; do python3 run_pipeline.py --psi $p --case psi_$p; done
"""
import argparse
import os
import sys

PROJ = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(PROJ, "turbine_mcp"))
os.environ.setdefault("PARABLADE_PATH", os.path.join(PROJ, "parablade-master", "parablade"))

from turbine_design.meanline import StageInputs, DesignCoeffs, run_meanline
from turbine_design.bladedesign import BladeChoices, design_blade
from turbine_design.geometry import generate_airfoil
from turbine_design.latexreport import write_latex_report


def main() -> None:
    ap = argparse.ArgumentParser(description="Axial-turbine design pipeline (headless).")
    ap.add_argument("--case", default="run", help="case name -> reports/<case>/")
    ap.add_argument("--psi", type=float, default=1.8, help="stage loading")
    ap.add_argument("--phi", type=float, default=0.5, help="flow coefficient")
    ap.add_argument("--dor", type=float, default=0.35, help="degree of reaction")
    ap.add_argument("--gamma", type=float, default=1.33, help="ratio of specific heats")
    ap.add_argument("--stagger", type=float, default=39.0, help="blade stagger [deg]")
    ap.add_argument("--zweifel", type=float, default=1.0, help="Zweifel coefficient")
    ap.add_argument("--cfd", action="store_true", help="also run GMSH+SU2 CFD (minutes)")
    ap.add_argument("--no-report", action="store_true", help="skip the LaTeX report")
    args = ap.parse_args()

    outdir = os.path.join(PROJ, "reports", args.case)
    os.makedirs(outdir, exist_ok=True)

    # 1) mean-line ---------------------------------------------------------------
    ml = run_meanline(StageInputs(gamma=args.gamma),
                      DesignCoeffs(psi=args.psi, phi=args.phi, DOR=args.dor))
    s2, s3, a = ml["station2"], ml["station3"], ml["annulus"]
    print(f"[mean-line] U={ml['derived']['U_mps']:.1f} m/s  "
          f"beta2={s2['beta_deg']:.2f}  beta3={s3['beta_deg']:.2f}  "
          f"constraints_all_pass={ml['constraints']['all_pass']}")

    # 2) rotor blade -------------------------------------------------------------
    bl = design_blade(s2["beta_deg"], s3["beta_deg"], s2["alpha_deg"], s3["alpha_deg"],
                      a["h2_m"], a["D_mean_m"],
                      BladeChoices(AR=1.0, stagger_deg=args.stagger, zweifel=args.zweifel),
                      span_h3_m=a["h3_m"])
    print(f"[blade] chord={bl['chord_m']:.3f} m  pitch={bl['pitch_m']:.3f} m  "
          f"Nb={bl['n_blades']}  eps={bl.get('flare_angle_deg', float('nan')):.2f} deg")

    # 3) ParaBlade airfoil -------------------------------------------------------
    blade_dat = os.path.join(outdir, "blade_coordinates.dat")
    generate_airfoil(bl["parablade"], bl["axial_chord_m"], out_dat=blade_dat,
                     pitch_m=bl["pitch_m"], parablade_path=os.environ["PARABLADE_PATH"])
    print(f"[airfoil] {blade_dat}")

    # 4) optional CFD ------------------------------------------------------------
    cfd_data = None
    if args.cfd:
        from turbine_design import cfd as CFD
        su2 = os.environ.get("SU2_BIN", "/Users/sanju/Downloads/bin/SU2_CFD")
        wd = os.path.join(outdir, "cfd_run")
        os.makedirs(wd, exist_ok=True)
        mesh_su2 = os.path.join(wd, "blade_mesh.su2")
        mstat = CFD.build_cascade_mesh(blade_dat, bl["pitch_m"], bl["axial_chord_m"],
                                       s2["beta_deg"], s3["beta_deg"], mesh_su2)
        cfg = os.path.join(wd, "blade_su2.cfg")
        CFD.write_su2_config(cfg, mesh_su2, s2["T0rel_K"], s2["P0rel_Pa"],
                             s2["beta_deg"], s3["P_Pa"], bl["axial_chord_m"],
                             mach_in_rel=s2["Mrel"], gamma=args.gamma, R=287.0,
                             inner_iter=4000, conv_minval=-9.0)
        print(f"[cfd] mesh {mstat['n_nodes']} nodes; running SU2 (minutes)...")
        run = CFD.run_su2(su2, cfg, wd, conv_minval=-9.0)
        if not run.get("success"):
            print("[cfd] SU2 FAILED:", run)
        else:
            post = CFD.postprocess(run["flow_vtu"], run["surface_csv"], blade_dat,
                                   bl["axial_chord_m"], bl["pitch_m"], s2["P0rel_Pa"],
                                   beta3_deg=s3["beta_deg"], gamma=args.gamma)
            cfd_data = {
                "blade_dat": blade_dat, "surface_csv": run["surface_csv"],
                "flow_vtu": run["flow_vtu"], "history_csv": os.path.join(wd, "history.csv"),
                "axial_chord_m": bl["axial_chord_m"], "pitch_m": bl["pitch_m"],
                "beta2_deg": s2["beta_deg"], "beta3_deg": s3["beta_deg"],
                "T0rel_K": s2["T0rel_K"], "P0rel_Pa": s2["P0rel_Pa"], "P3_Pa": s3["P_Pa"],
                "mach_in_rel": s2["Mrel"], "gamma": args.gamma, "R": 287.0, "mesh": mstat,
                "convergence": {"final_rms_density": run["final_rms_density"],
                                "iterations": run["iterations"], "converged": run["converged"],
                                "target_rms_density": -9.0},
                "results": post,
            }
            print(f"[cfd] converged={run['converged']}  rms={run['final_rms_density']:.2f}  "
                  f"iters={run['iterations']}")

    # 5) report ------------------------------------------------------------------
    if not args.no_report:
        out = write_latex_report(ml, outdir, case_name=args.case, blade=bl,
                                 blade_dat=blade_dat, cfd=cfd_data,
                                 title="Axial turbine stage design")
        print(f"[report] {out.get('pdf') or out.get('tex')}  compiled={out.get('compiled')}")


if __name__ == "__main__":
    main()
