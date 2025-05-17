"""
Plotting functions for infercnvpy-mcp.
"""

import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
from fastmcp import FastMCP, Context
from pydantic import Field
from ..schema.pl import ChromosomeHeatmapModel, UMAPPlotModel
from ..util import save_figure
from scmcp_shared.util import add_op_log, forward_request
from scmcp_shared.logging_config import setup_logger

logger = setup_logger()

pl_mcp = FastMCP("InfercnvpyMCP-Plotting-Server")

@pl_mcp.tool()
async def chromosome_heatmap(
    ctx: Context,
    request: ChromosomeHeatmapModel,
    dtype: str = Field(default="exp", description="the datatype of anndata.X"),
    sampleid: str = Field(default=None, description="adata sampleid for analysis"),
):
    """
    Plot a heatmap of gene expression by chromosome.
    This visualization shows the smoothed gene expression values along chromosomes, 
    helping to identify regions with amplifications or deletions.
    It is the main visualization tool for examining CNV patterns.
    """
    try:
        result = await forward_request("pl_chromosome_heatmap", request, sampleid=sampleid, dtype=dtype)
        if result is not None:
            return result
        
        ads = ctx.request_context.lifespan_context
        adata = ads.get_adata(dtype=dtype, sampleid=sampleid)
        
        # Extract request parameters
        kwargs = request.model_dump()
        
        # Create figure
        save_param = kwargs.pop("save", None)
        show_param = kwargs.pop("show", True)
        
        # Plot chromosome heatmap
        fig = plt.figure(figsize=kwargs.pop("figsize", None))
        cnv.pl.chromosome_heatmap(adata, show=False, **kwargs)
        
        # Save figure
        fig_name = save_figure(plt.gcf(), name=save_param)
        
        # Log the operation
        add_op_log(adata, "pl_chromosome_heatmap", kwargs)
        
        # Return information
        return {
            "status": "success",
            "message": "Chromosome heatmap plotted successfully.",
            "figure": fig_name,
            "groupby": kwargs.get("groupby", None)
        }
    except KeyError as e:
        raise e
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e

@pl_mcp.tool()
async def umap(
    ctx: Context,
    request: UMAPPlotModel,
    dtype: str = Field(default="exp", description="the datatype of anndata.X"),
    sampleid: str = Field(default=None, description="adata sampleid for analysis"),
):
    """
    Plot UMAP embedding of CNV profiles.
    This visualization shows the cells in a UMAP embedding based on their CNV profiles,
    which can help identify tumor subclones and normal cells.
    It should be called after running infercnv, neighbors, and umap computation.
    """
    try:
        result = await forward_request("pl_umap", request, sampleid=sampleid, dtype=dtype)
        if result is not None:
            return result
        
        ads = ctx.request_context.lifespan_context
        adata = ads.get_adata(dtype=dtype, sampleid=sampleid)
        
        # Extract request parameters
        kwargs = request.model_dump()
        
        # Create figure
        save_param = kwargs.pop("save", None)
        show_param = kwargs.pop("show", True)
        return_fig_param = kwargs.pop("return_fig", False)
        
        # Plot UMAP
        if return_fig_param:
            fig = cnv.pl.umap(adata, show=False, return_fig=True, **kwargs)
        else:
            fig = plt.figure(figsize=kwargs.pop("figsize", None))
            cnv.pl.umap(adata, show=False, **kwargs)
        
        # Save figure
        fig_name = save_figure(plt.gcf(), name=save_param)
        
        # Log the operation
        add_op_log(adata, "pl_umap", kwargs)
        
        # Return information
        return {
            "status": "success",
            "message": "UMAP plot created successfully.",
            "figure": fig_name,
            "color": kwargs.get("color", None)
        }
    except KeyError as e:
        raise e
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e 