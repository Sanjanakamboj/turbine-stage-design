# Turbine Stage Design

Axial-turbine stage-design chain (1D mean-line → blade geometry → cascade CFD),
exposed to Claude / any MCP client as the **`turbine`** MCP server.

The repo ships a project-scoped [`.mcp.json`](.mcp.json), so once you clone it and
install the Python deps, the `turbine` server shows up in Claude Code's `/mcp`
with **no path editing** — the server self-locates everything relative to the repo.

## Tools

| Tool | What it does | Needs |
|------|--------------|-------|
| `turbine_design_meanline` | 1D mean-line design → stations, velocity triangles, annulus, constraints | pure Python |
| `turbine_design_stage` | mean-line **+** rotor-blade geometry | pure Python |
| `turbine_generate_airfoil` | build the airfoil with the bundled ParaBlade | pure Python |
| `turbine_run_cascade_cfd` | GMSH mesh + SU2 RANS solve + post-process | **SU2 binary** + gmsh |

The first three are pure Python and work out of the box. Only
`turbine_run_cascade_cfd` needs the heavy native toolchain.

## Setup (gets `/mcp` working)

### 1. Clone
```bash
git clone <repo-url> turbine-stage-design
cd turbine-stage-design
```

### 2. Install Python dependencies
Use the **same `python3`** that Claude Code will launch (the `.mcp.json` calls
`python3`). A venv is fine — just make sure it's the active `python3` when you run
Claude, or point `.mcp.json`'s `"command"` at the venv's python.
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

### 4. Use it in Claude Code
Open the folder with Claude Code. Project-scoped MCP servers from `.mcp.json`
require a one-time approval — accept the prompt, then:
```
/mcp           # should list "turbine ✔ connected"
```
Or check from a terminal:
```bash
claude mcp list
# turbine: python3 turbine_mcp/turbine_mcp_server.py - ✔ Connected
```

## Overriding paths (optional)
Everything defaults relative to the repo. Override via env vars only if needed:

| Env var | Default | Purpose |
|---------|---------|---------|
| `TURBINE_PROJECT_DIR` | repo root | where outputs are written |
| `SU2_BIN` | `SU2_CFD` (on PATH) | SU2 solver binary |
| `PARABLADE_PATH` | `parablade-master/parablade` | bundled ParaBlade library |

## Layout
- `turbine_mcp/` — the MCP server + the importable `turbine_design/` package (design logic, CFD, reports) + tests
- `parablade-master/` — vendored [ParaBlade](https://github.com/RoberAgro/parablade) (GPLv3), used by the airfoil tool
- `.mcp.json` — project MCP server registration (portable, relative paths)

## Notes
- Generated artifacts (CFD runs, meshes, `*.vtu`/`*.su2`, report PDFs, case log)
  are **git-ignored** — they're produced by the tools, not committed.
- ParaBlade is bundled under GPLv3 (see `parablade-master/LICENSE`); its heavy
  `testcases/` and `docs/` are git-ignored to keep the clone small.
