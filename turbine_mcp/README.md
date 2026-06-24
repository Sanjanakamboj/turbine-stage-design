# Turbine MCP server

An [MCP](https://modelcontextprotocol.io) server that exposes the axial-turbine
stage-design chain as tools an LLM (Claude, etc.) can call.

## Tools

| Tool | What it does | Cost | Needs |
|------|--------------|------|-------|
| `turbine_design_meanline` | 1D mean-line design → stations, velocity triangles, annulus, constraint check | instant | pure Python |
| `turbine_design_stage` | mean-line **+** rotor-blade geometry (chord, stagger, pitch, blade count, edge radii) | instant | pure Python |
| `turbine_generate_airfoil` | build the airfoil with ParaBlade, save coordinates | ~1 s | parablade |
| `turbine_run_cascade_cfd` | GMSH 3-blade mesh + SU2 RANS solve + post-process loading & exit flow | minutes | SU2 binary, gmsh, pyvista |

The design logic lives in the importable `turbine_design/` package
(`meanline.py`, `bladedesign.py`, `geometry.py`, `cfd.py`) and is reused by the
server. The mean-line module reproduces the project's reference design to <0.1 %.

### Reports + case-log (returned by every tool)

Every tool writes a **sectioned LaTeX report** (compiled to PDF) and appends a row
to a shared **Excel case-log**, returning their paths (`report_pdf`, `case_log`).
Disable the report with `make_report: false`. Pass `case_name` to control the
report folder and the log row (re-running the same name updates that row, so a
design and its later CFD evaluation land on one line).

- `turbine_design_meanline` → report §1 (mean-line): input/coefficient/station
  tables, h-s diagram, annulus, constraint check.
- `turbine_design_stage` → §1 + §2 (rotor blade): also builds the airfoil and adds
  blade-parameter table, section, cascade, thickness & curvature distributions.
- `turbine_run_cascade_cfd` → §3 (CFD): geometry, mesh, SU2-config table, blade
  loading, Mach & pressure fields, pitchwise & wake, convergence history.

Reports land in `reports/<case_name>/report.pdf`; the log in `design_cases_log.xlsx`.

Rendering lives in importable modules that work on any result dict:
`turbine_design/report.py` (figures + matplotlib PDF reports), `latexreport.py`
(`write_latex_report` — sectioned `.tex` + `latexmk`/`pdflatex` compile), and
`caselog.py` (`append_case` — the Excel log). Needs `matplotlib`, `openpyxl`, and
a LaTeX toolchain. Every figure panel is guarded — an unconverged field or a
too-short probe draws a "data unavailable" note instead of failing the report.

### CFD convergence

`turbine_run_cascade_cfd` reports a `converged` flag (true only when the final
RMS-density residual reaches the `-9` target) and adds a `warning` when it does
not — a non-converged field produces unreliable loading/exit numbers, so do not
trust the post-processed results unless `converged` is true.

Earlier versions limit-cycled at ~`-5.6` (previously misattributed to GMSH
versions). The real cause was a sign error in the outlet wall extension
(`so = +tan(β3)` instead of `-tan(β3)`): the rotor exit flow leaves at `-β3`
(downward), so the outlet domain must slant down with it. With the fix the
reference case converges to `-9.0` in ~1300 iters and the CFD exit angle/Mach
match the mean-line design to ~1%. `tests/test_mesh.py` guards the sign. Note
that SU2 ignores `REYNOLDS_NUMBER` in `DIMENSIONAL` mode — the viscosity comes
from its own model.

## Install

```bash
pip install -r requirements.txt
```

`numpy`, `gmsh`, `pyvista`, `pandas`, `scipy` are only needed for the CFD tool.
The CFD tool also needs the `SU2_CFD` binary.

## Run

```bash
python3 turbine_mcp_server.py          # stdio transport
```

## Test

```bash
python3 -m pytest tests/               # fast: mean-line, blade, input validation
RUN_CFD_TESTS=1 python3 -m pytest tests/test_cfd_integration.py   # heavy GMSH+SU2 run
```

`tests/` pins the reference HP-stage numbers, the constraint pass/fail logic
(incl. off-design flips), the blade geometry and the pydantic input bounds. The
GMSH+SU2 end-to-end test is skipped unless `RUN_CFD_TESTS=1` and the `SU2_CFD`
binary and a `blade_coordinates.dat` are present.

You can also drive the live server interactively with the MCP Inspector
(`mcp dev turbine_mcp_server.py`) to browse tool schemas and fire calls by hand.

## Configuration (environment variables)

| Variable | Default |
|----------|---------|
| `TURBINE_PROJECT_DIR` | project root (working dir for generated files) |
| `SU2_BIN` | path to the `SU2_CFD` binary |
| `PARABLADE_PATH` | path to the `parablade` package |

## Use with Claude Code

A project-scoped `../.mcp.json` registers this server as `turbine`. Claude Code
auto-detects it; approve the server when prompted, then ask e.g.
*"design an HP turbine stage with ψ=1.8, φ=0.5, reaction 0.35"* and Claude will
call `turbine_design_stage`.

## Example (direct Python)

```python
from turbine_design import StageInputs, DesignCoeffs, run_meanline
r = run_meanline(StageInputs(), DesignCoeffs(psi=1.8, phi=0.5, DOR=0.35))
print(r["derived"]["U_mps"], r["station3"]["beta_deg"])
```
