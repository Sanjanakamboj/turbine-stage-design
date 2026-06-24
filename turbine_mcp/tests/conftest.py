"""Make the server module and the turbine_design package importable from tests.

Adds the turbine_mcp/ directory (parent of tests/) to sys.path so that both
`import turbine_design` and `import turbine_mcp_server` resolve when pytest is
run from anywhere.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_SERVER_DIR = os.path.dirname(_HERE)
if _SERVER_DIR not in sys.path:
    sys.path.insert(0, _SERVER_DIR)
