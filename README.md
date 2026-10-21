# Turbine Stage Design

[![CI](https://github.com/Sanjanakamboj/turbine-stage-design/actions/workflows/ci.yml/badge.svg)](https://github.com/Sanjanakamboj/turbine-stage-design/actions/workflows/ci.yml)

Axial-turbine stage-design chain (1D mean-line → blade geometry → cascade CFD),
exposed to any MCP client as the **`turbine`** MCP server.

The repo ships a project-scoped [`.mcp.json`](.mcp.json), so once you clone it and
install the Python deps, an MCP client that reads project-scoped config can launch
the `turbine` server with **no path editing** — the server self-locates everything
relative to the repo.

## Tools

| Tool | What it does | Needs |
|------|--------------|-------|
| `turbine_design_meanline` | 1D mean-line design → stations, velocity triangles, annulus, constraints | pure Python |
| `turbine_design_stage` | mean-line **+** rotor-blade geometry | pure Python |
| `turbine_generate_airfoil` | build the airfoil with the bundled ParaBlade | pure Python |
| `turbine_run_cascade_cfd` | GMSH mesh + SU2 RANS solve + post-process | **SU2 binary** + gmsh |

The first three are pure Python and work out of the box. Only
`turbine_run_cascade_cfd` needs the heavy native toolchain.

## Design decisions

A few choices that shape the architecture — the reasoning, in 30 seconds:

- **Four tools split on a pure-Python ↔ solver boundary.** The chain is one
  pipeline (mean-line → blade → airfoil → CFD), but it's exposed as four tools cut
  at the point where the *cost and dependencies* change. The first three are
  instant, deterministic and need no native toolchain, so the design space can be
  swept cheaply and the tests run anywhere; only `turbine_run_cascade_cfd` crosses
  into GMSH + the SU2 binary. Keeping that boundary explicit means an LLM (or a CI
  job) can do all the design reasoning without ever installing a solver, and the
  expensive step is opt-in rather than load-bearing.

- **Physics lives in an importable package, not in the server.** All the
  engineering is in `turbine_design/` (`meanline.py`, `bladedesign.py`,
  `geometry.py`, `cfd.py`); `turbine_mcp_server.py` is a thin pydantic-validated
  adapter, and [`run_pipeline.py`](run_pipeline.py) is a headless CLI over the same
  functions. The MCP server, the CLI and the tests are three front-ends onto one
  core — no logic is duplicated in the tool layer.

- **Rendering is separated from computation by layer ownership.** Each physics
  module owns its own plots (`meanline.py` the h-s/triangles/annulus, `cfd.py` the
  fields/loading/wake); `report.py` is a pure orchestrator that saves PNGs and
  `latexreport.py` lays out the PDF. matplotlib is lazily imported inside the plot
  functions so the design math still runs on a headless box with no plotting stack.

- **Infeasible designs fail loudly, not silently.** The mean-line returns a
  `feasibility`/`warnings` block: if blade speed is clamped to the mechanical limit
  or the work-match doesn't converge, that surfaces instead of returning
  plausible-looking but wrong stations. The CFD tool likewise reports a `converged`
  flag tied to the SU2 residual target and refuses to vouch for an unconverged
  field. Engineering judgement is encoded, not assumed.

- **The "tribal knowledge" is a file, not folklore.** Valid input ranges, the
  design-constraint bands, where the 1D method breaks down, and how to read the CFD
  output live in [`turbine_mcp/SKILLS.md`](turbine_mcp/SKILLS.md) with pointers
  back to the code that enforces each number — so the *why* travels with the *what*.

## Setup (gets `/mcp` working)

### 1. Clone (with the ParaBlade submodule)
ParaBlade is pulled in as a git submodule, so clone recursively:
```bash
git clone --recurse-submodules https://github.com/Sanjanakamboj/turbine-stage-design.git
cd turbine-stage-design
```
Already cloned without `--recurse-submodules`? Fetch it now:
```bash
git submodule update --init parablade-master
```
(The submodule is only needed for the airfoil/stage/CFD tools; the mean-line tool
works without it.)

### 2. Install Python dependencies
Use the **same `python3`** that your MCP client will launch (the `.mcp.json` calls
`python3`). A venv is fine — just make sure it's the active `python3` when you start
the client, or point `.mcp.json`'s `"command"` at the venv's python.
```bash
python3 -m pip install -r turbine_mcp/requirements.txt
```

### 3. (CFD only) install SU2
`turbine_run_cascade_cfd` shells out to SU2. Install it from
<https://su2code.github.io/> and make sure the binary is reachable:
```bash
# either put SU2_CFD on your PATH, or set an explicit path:
export SU2_BIN=/abs/path/to/SU2_CFD
```
Skip this if you only need the design/airfoil tools.

### 4. Use it from an MCP client
Point any MCP-compatible client at [`.mcp.json`](.mcp.json) (it launches
`python3 turbine_mcp/turbine_mcp_server.py` over stdio), or run the server by hand:
```bash
python3 turbine_mcp/turbine_mcp_server.py
```
To browse the tool schemas and fire calls interactively, use the MCP Inspector:
```bash
cd turbine_mcp && mcp dev turbine_mcp_server.py
```
No MCP client? [`run_pipeline.py`](run_pipeline.py) drives the same chain from the
command line.

## Overriding paths (optional)
Everything defaults relative to the repo. Override via env vars only if needed:

| Env var | Default | Purpose |
|---------|---------|---------|
| `TURBINE_PROJECT_DIR` | repo root | where outputs are written |
| `SU2_BIN` | `SU2_CFD` (on PATH) | SU2 solver binary |
| `PARABLADE_PATH` | `parablade-master/parablade` | bundled ParaBlade library |

## Layout
- `turbine_mcp/` — the MCP server + the importable `turbine_design/` package (design logic, CFD, reports) + tests
- `parablade-master/` — **git submodule** → [ParaBlade](https://github.com/NAnand-TUD/parablade) (GPLv3), used by the airfoil/report tools
- `.mcp.json` — project MCP server registration (portable, relative paths)

## Notes
- Generated artifacts (CFD runs, meshes, `*.vtu`/`*.su2`, report PDFs, case log)
  are **git-ignored** — they're produced by the tools, not committed.
- ParaBlade is referenced as a submodule (pinned upstream commit), not redistributed
  here. It needs `numpy`/`scipy`/`matplotlib` (already in `requirements.txt`) —
  there's no need to run ParaBlade's own `setup.py`.
- `turbine_design/geometry.py` adds a small matplotlib-compat shim so upstream
  ParaBlade's `plot_*` methods work on matplotlib ≥ 3.8 (which removed `Tick.label`).
