# SKILLS — engineering knowledge codified for the turbine design chain

This file is the *tribal knowledge* behind the `turbine` MCP tools: the input
ranges the models are actually valid over, the constraints that decide whether a
design is sensible, where the 1D mean-line method stops being trustworthy, and
how to read the CFD output. It is written so that an LLM (or a new engineer)
calling the tools knows **what a number means and when not to believe it** —
not just how to call the function.

The numbers below are the ones enforced or checked in code, with file pointers so
they can't drift from the implementation:
- input bounds → `turbine_mcp_server.py` (pydantic `Field(...)`)
- design defaults → `turbine_design/meanline.py` (`StageInputs`, `DesignCoeffs`)
- constraint bands → `turbine_design/meanline.py::_check_constraints`
- CFD config / convergence → `turbine_design/cfd.py`

---

## 1. The four tools and the pure-Python ↔ solver boundary

| Tool | Layer | Cost | Hard dependency |
|------|-------|------|-----------------|
| `turbine_design_meanline` | 1D thermo/kinematics | instant | none (pure Python) |
| `turbine_design_stage` | mean-line **+** rotor-blade sizing | instant | none (pure Python) |
| `turbine_generate_airfoil` | ParaBlade camber/thickness → `.dat` | ~1 s | ParaBlade (vendored) |
| `turbine_run_cascade_cfd` | GMSH mesh + SU2 RANS + post | minutes | `SU2_CFD` binary, gmsh, pyvista |

The boundary is deliberate: the first three need **no native toolchain**, so the
design space can be swept cheaply and deterministically. Only the last tool
crosses into the solver world. Treat the mean-line as the source of truth for the
*design intent* and CFD as the *verification* — when they disagree by more than a
few percent, suspect the CFD setup (mesh, BCs, convergence) before the design.

---

## 2. Valid input ranges

Hard bounds are enforced by pydantic and a bad value is rejected before any
physics runs. "Designed-for" is the regime the correlations and assumptions were
built around (the LMECA2323 HP-stage brief); going outside it is allowed but the
results need a second look.

### Design coefficients (the variables you actually sweep)

| Input | Hard bound (rejected outside) | Designed-for | Notes |
|-------|-------------------------------|--------------|-------|
| `psi` (stage loading ΔH/U²) | `> 0` | **1.5 – 2.0** (HP) | Above ~2.5 the stage is over-loaded → high turning, thick blades, loss. Below ~1.0 it demands U beyond the speed limit (see §4). |
| `phi` (flow coefficient Vx/U) | `> 0` | **0.4 – 0.7** | Low φ → small annulus, high swirl; high φ → low turning, weak work. |
| `DOR` (degree of reaction) | `0 … 1` | **0.3 – 0.5** | 0 = impulse, 0.5 = 50 % reaction. >0.5 pushes exit swirl α₃ past its cap (§4). This is a *pressure-based* reaction, P₂ = DOR·(P₁−P₃)+P₃. |
| `M3` (target exit Mach) | `0 < M3 < 1` | 0.3 | **Dead input** — declared but not consumed; exit Mach floats out of the work match. See §4. |
| `eta_TT` | `0 < … ≤ 1` | 0.87 | **Dead input** — not used by the solver; `eta_row` (below) drives the losses instead. |

### Boundary conditions / machine inputs

| Input | Hard bound | Default (reference) | Notes |
|-------|------------|---------------------|-------|
| `TIT_C` (turbine inlet temp) | `> 0 °C` | 1400 °C | Real metal/cooling limits sit ~1500–1700 °C; the model has **no cooling or material limit** — it will happily accept 2000 °C and give nonsense. Stay in the 1100–1600 °C band for credible gas-turbine designs. |
| `pi_C` (compressor PR) | `> 1` | 20 | Sets P₁ = π_C·p_adm and the compressor power that the turbine must supply. |
| `mdot` (mass flow) | `> 0 kg/s` | 750 | Sizes the annulus (A = ṁ/ρVx); does not change the per-kg thermodynamics. |
| `P_shaft_MW` (net shaft power) | `> 0 MW` | 330 | **Total** for the machine; split across `n_stages`. |
| `n_stages` | `≥ 1` | 4 | Power is divided **equally** per stage (see assumption in §4). |
| `M_comb` (combustor-exit Mach) | `0 < … < 1` | 0.15 | Sets the inlet total/static split. |
| `gamma` | `> 1` | 1.333 | Combustion-gas value; **don't** pass air's 1.4 for the hot section. |
| `R` | `> 0` | 287 | |
| `eta_isC`, `eta_m`, `eta_row` | `0 < … ≤ 1` | 0.92, 0.95, 0.90 | `eta_row` is the **fixed** stator/rotor efficiency that produces all the entropy in the mean-line — it is assumed, **not** estimated from a loss correlation. |
| `U_limit` | (in `StageInputs`) | 350 m/s | Mechanical peripheral-speed cap; the design is clamped to it (§4). |

### Blade-design inputs (`turbine_design_stage`, `turbine_generate_airfoil`)

| Input | Hard bound | Default | Notes |
|-------|------------|---------|-------|
| `AR` (aspect ratio span/chord) | `> 0` | 1.0 | |
| `stagger_deg` | (free) | 39° | Drives axial chord = chord·cos(stagger). |
| `zweifel` | `> 0` | 1.0 | Loading coefficient → pitch → blade count. Lower Zweifel ⇒ tighter pitch ⇒ more blades. |

> ⚠️ **Sign convention trap:** the Zweifel/pitch formula needs the mean-line's
> **positive** β₃ (~+69.5°). ParaBlade, however, wants stagger and exit angle
> **negated** and the inlet angle as-is. The stage tool handles this; if you call
> `design_blade` / `generate_airfoil` by hand, feeding a negative β₃ gives a
> negative pitch and blade count. (`bladedesign.py`, guarded by `test_blade.py`.)

---

## 3. Design-constraint bands (`_check_constraints`)

Every mean-line result carries a `constraints` block. Each entry is
`{value, min, max, pass}` and `all_pass` is the AND. These bands encode the
"is this a sane HP stage?" judgement that a designer would otherwise apply by eye.

| Constraint | Band | What it protects against |
|------------|------|--------------------------|
| `M2` (stator-exit abs. Mach) | **0.70 – 0.85** | Too low → under-expanded stator; too high → transonic stator losses / shocks. |
| `beta2_deg` (rotor inlet rel. angle) | **≤ 47.5°** | Excessive incidence / rotor inlet blockage. |
| `Mw2` (rotor inlet rel. Mach) | **≤ 0.50** | Keeps the rotor inlet subsonic with margin. |
| `Mw3` (rotor exit rel. Mach) | **0.65 – 0.80** | The rotor should accelerate the relative flow but stay below transonic exit. |
| `beta3_deg` (rotor exit rel. angle) | **≤ 75°** | Very high turning ⇒ thin trailing flow, separation risk. |
| `turning_deg` (β₂+β₃) | **110 – 120°** | The HP-turning sweet spot; below → weak work, above → loss-heavy. |
| `alpha3_deg` (exit swirl) | **≤ 35°** | Residual swirl the next stage / diffuser must absorb. |
| `h3_h2` (span ratio) | **≤ 1.20** | Annulus flare; steep divergence ⇒ end-wall separation. |
| `AN2_m2rpm2` (AN²) | **≤ 3×10⁷** | Disk/blade-root stress proxy (area × speed²). |

**Reading a failing constraint:** the default reference design intentionally sits
on the edge — `M2 ≈ 0.699` (just under 0.70) and `h3/h2 > 1.20` both "fail." That
is expected and documented; it reflects the real brief, not a bug. A failure is a
flag for *engineering judgement*, not an automatic reject — but two or more
failures, or a large margin of failure, means the (ψ, φ, DOR) point is poor.

Off-design illustrations (pinned in `test_meanline.py`):
- ψ=1.6, φ=0.6, DOR=0.4 → turning drops **below 110°** (under-worked).
- ψ=2.0, φ=0.45, DOR=0.5 → exit swirl **α₃ > 35°** (high-reaction penalty).

---

## 4. When the mean-line breaks down

The 1D method is a **design-point** tool with explicit assumptions. It stops being
trustworthy here:

1. **Blade-speed cap hit (`feasibility.U_capped`).** Loading sets
   U = √(ΔH/ψ). If that exceeds `U_limit` (350 m/s), U is clamped and the
   *effective* loading ψ_eff = ΔH/U² is reported in a warning. Low ψ (≈1.0) or
   high per-stage power triggers this. When capped, the requested ψ is **not** the
   delivered ψ — read `feasibility.psi_effective`.

2. **Work match did not converge (`feasibility.work_converged == False`).** P₃ is
   solved by a bracketed secant loop so that delivered ΔH₀ equals the required
   specific work. If the design is infeasible (e.g. 120 MW / 3 stages at ψ=1.5),
   the loop can't match and the residual is large. **Every station, angle and
   constraint downstream is then unreliable** — the tool says so in `warnings`.
   Fix by adding stages, lowering ψ, or raising `U_limit`.

3. **High reaction (DOR → 0.5 and up).** Pushes exit swirl α₃ past 35° and erodes
   the rotor turning split. The constraint check catches it, but the 1D loss model
   (fixed `eta_row`) does **not** capture the extra secondary loss — CFD will be
   worse than the mean-line suggests.

4. **Transonic regimes.** The model is isentropic-with-fixed-efficiency. As M₂ → 0.85+
   or Mw₃ → 0.80+ the real flow develops shocks and shock/boundary-layer
   interaction that a constant `eta_row` cannot represent. Inside the bands (§3)
   you're fine; outside, trust CFD, not the mean-line.

5. **Multi-stage interaction is not modelled.** Total power is split **equally**
   across `n_stages` and each stage is designed in isolation at the same (ψ, φ,
   DOR). Real machines load stages unequally and have inter-stage swirl/leakage
   coupling. Use this for a first cut, not a final multi-stage matching.

6. **Losses are assumed, not predicted.** `eta_row` is a fixed input; there is **no
   profile/secondary/tip-leakage loss correlation** (the reference workflow's
   loss-estimation step is intentionally not implemented). Efficiency out ≈
   efficiency in.

7. **Known dead inputs.** `M3` and `eta_TT` are declared design coefficients but
   are **never consumed** — exit Mach floats out of the work match rather than
   being targeted. Don't tune them expecting an effect.

8. **Small representational quirks (not bugs):** T_ref = 273.15 K (vs 273.0 in the
   original notebook → ≤0.03 % drift); station-1 absolute velocity is reported as
   Vx = φU rather than from M_comb (a documented inlet-block inconsistency carried
   over from the brief). These are pinned, expected, and ≤0.1 %.

---

## 5. Interpreting the SU2 / CFD output

`turbine_run_cascade_cfd` returns a `converged` flag and post-processed loading and
exit flow. Read them in this order:

1. **`converged` first — gate everything on it.** It is `True` *only* when the
   final RMS-density residual reaches the **−9** target (`cfd.py::run_su2`). Hitting
   the iteration cap is **not** convergence. A non-converged field gives unreliable
   loading and exit numbers — the tool attaches a warning; do not quote those
   results.

2. **Cross-check against the mean-line.** On the reference case the converged CFD
   exit angle (≈ −69.5°) and exit Mach (≈ 0.69) match the mean-line design (β₃,
   Mw₂) to ~1 %. If CFD and mean-line disagree by ≫ a few percent, suspect the CFD
   setup (BCs, mesh, not-actually-converged) before rewriting the design.

3. **Entropy rise = loss, and it should localize in the wake.** The entropy-rise
   colormap (Δs = cₚ·ln(T/T_t0) − R·ln(P/P_t0)) should be ~zero in the free stream
   and concentrated in the blade wake and boundary layers. Loss spread across the
   whole passage means a mesh/BC problem, not real physics.

4. **Total enthalpy should be ~uniform** across the passage (rothalpy is conserved
   in the relative frame). A non-uniform total-enthalpy field is a red flag.

5. **Mach limits.** A peak suction-side Mach pushing past ~1 means the cascade has
   gone transonic — expect a shock, a kink in the isentropic-Mach loading, and
   loss the mean-line never predicted.

6. **Reynolds caveat.** SU2 runs in `DIMENSIONAL` mode and **ignores
   `REYNOLDS_NUMBER`** — viscosity comes from its own (Sutherland) model. The Re we
   compute is informational; don't expect changing it to move the solution.

---

## 6. Known failure modes (and what they actually mean)

| Symptom | Real cause | Fix / status |
|---------|------------|--------------|
| SU2 limit-cycles at residual ≈ **−5.6**, never reaches −9 | **β₃ sign error** in the outlet-wall extension: rotor exit leaves at −β₃ (downward), so the outlet domain must slant *down* (`so = −tan β₃`). The original used `+tan β₃`, slanting the domain against the flow → malformed field. | **Fixed** (`cfd.py::build_cascade_mesh`, guarded by `test_mesh.py`). Was previously misattributed to GMSH versions. Now converges to −9 in ~1300 iters. |
| GMSH meshing fails / degenerate cells | Bad/insufficient `blade_coordinates.dat` (self-intersecting or too few points), or pitch/chord inconsistent with the blade. | Regenerate the airfoil with `turbine_generate_airfoil`; ensure the `.dat`, `pitch_m`, `axial_chord_m` and β angles all come from the **same** design. |
| SU2 reports `converged=False` at the iteration cap | Genuinely hard case (transonic / off-design), too few `inner_iter`, or wrong BCs. | Raise `inner_iter` (bounded 100–20000), check the design is feasible (§4), confirm β angles/pressures match the mean-line. **Don't trust the post-processed numbers** until it's true. |
| Pitchwise / wake profiles look flat or jagged | The exit probe slice was too sparse (historical) — now a dense 401-pt line probe centred on the wake streamline. | Resolved (`cfd.py::_exit_probe`); a too-short probe now draws a "data unavailable" note rather than lying. |
| Report figure shows a "data unavailable" box | The underlying field/probe was missing or the run didn't converge. | Expected guard behaviour — fix the upstream run, then the panel fills in. |
| `SU2 binary not found` | `SU2_CFD` not on PATH and `SU2_BIN` unset. | Install SU2 and set `SU2_BIN` (see README). The other three tools don't need it. |

---

## 7. Quick heuristics

- Start from the **reference point** (ψ=1.8, φ=0.5, DOR=0.35) and move one
  coefficient at a time; watch which constraint flips.
- If U pins at 350 m/s, you're **speed-limited** — lower power-per-stage (more
  stages) or accept lower loading before touching anything else.
- Want more work per stage? Raise ψ — but watch turning (→120° cap) and Mw₃.
- Want a smaller, faster machine? Raise φ — but turning drops toward the 110° floor.
- Treat any design with `work_converged=False` or ≥2 failed constraints as **not a
  design** — it's a diagnostic telling you to move in the feasible direction.
- Mean-line and CFD should agree to ~1 % on exit angle/Mach. When they don't, the
  CFD setup is guilty until proven innocent.
