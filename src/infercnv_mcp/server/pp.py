"""
Preprocessing functions for infercnvpy-mcp.
"""

import os
import inspect
import scanpy as sc
from fastmcp import FastMCP , Context
from fastmcp.exceptions import ToolError
from scmcp_shared.server.pp import ScanpyPreprocessingMCP
from scmcp_shared.util import update_mcp_args

pp_mcp = ScanpyPreprocessingMCP(
    include_tools=["neighbors"]
).mcp


 
update_mcp_args(
    pp_mcp, {
        "neighbors": {
            "use_rep": {"default": "cnv_pca"},
            "key_added": {"default": "cnv_neighbors"},
            "sampleid": {"default": "adata_cnv"},
            "adtype": {"default": "cnv"}
        }
    }
)
