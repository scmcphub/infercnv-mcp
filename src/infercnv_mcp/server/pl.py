"""
Plotting functions for infercnvpy-mcp.
"""

import infercnvpy as cnv
from fastmcp import FastMCP
from ..schema.pl import *
from ..schema import CNVAdataInfo
from scmcp_shared.util import forward_request, get_ads
from scmcp_shared.logging_config import setup_logger
from scmcp_shared.util import sc_like_plot
from fastmcp.exceptions import ToolError
from scmcp_shared.server import ScanpyPlottingMCP

logger = setup_logger()


pl_mcp = ScanpyPlottingMCP(
    include_tools=["embedding"]
).mcp

@pl_mcp.tool()
async def chromosome_heatmap(
    request: ChromosomeHeatmapParam,
    adinfo: CNVAdataInfo = CNVAdataInfo()
):
    """Plot a heatmap of smoothed gene expression by chromosome.
    """
    try:
        if (res := await forward_request("pl_chromosome_heatmap", request, adinfo)) is not None:
            return res
        adata = get_ads().get_adata(adinfo=adinfo)
        fig_path = sc_like_plot(cnv.pl.chromosome_heatmap, adata, request, adinfo)
        return {"figpath": fig_path}
    except ToolError as e:
        raise ToolError(e)
    except Exception as e:
        raise ToolError(e)

