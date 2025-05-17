"""
Tool functions for infercnvpy-mcp.
"""

import scanpy as sc
import infercnvpy as cnv
from fastmcp import FastMCP, Context
from pydantic import Field
from ..schema.tl import InferCNVModel, PCPModel, UMAPModel, LeidenModel, CNVScoreModel
from scmcp_shared.util import add_op_log, forward_request
from scmcp_shared.logging_config import setup_logger

logger = setup_logger()

tl_mcp = FastMCP("InfercnvpyMCP-Tools-Server")

@tl_mcp.tool()
async def infercnv(
    ctx: Context,
    request: InferCNVModel,
    dtype: str = Field(default="exp", description="the datatype of anndata.X"),
    sampleid: str = Field(default=None, description="adata sampleid for analysis"),
):
    """
    Infer copy number variations (CNVs) from scRNA-seq data.
    This tool identifies regions of the genome that appear to be subject to copy number variations 
    by analyzing the gene expression profiles of single cells, comparing them to reference (normal) cells.
    """
    try:
        result = await forward_request("tl_infercnv", request, sampleid=sampleid, dtype=dtype)
        if result is not None:
            return result
        
        ads = ctx.request_context.lifespan_context
        adata = ads.get_adata(dtype=dtype, sampleid=sampleid)
        
        # Extract request parameters
        kwargs = request.model_dump()
        
        # Run infercnv
        cnv.tl.infercnv(adata, **kwargs)
        
        # Log the operation
        add_op_log(adata, "tl_infercnv", kwargs)
        
        return {
            "status": "success",
            "message": "CNV inference completed successfully.",
            "reference_key": kwargs.get("reference_key"),
            "window_size": kwargs.get("window_size")
        }
    except KeyError as e:
        raise e
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e

@tl_mcp.tool()
async def pca(
    ctx: Context,
    request: PCPModel,
    dtype: str = Field(default="exp", description="the datatype of anndata.X"),
    sampleid: str = Field(default=None, description="adata sampleid for analysis"),
):
    """
    Perform PCA on CNV profiles.
    This function is specifically designed to work on CNV profiles and should be called after running infercnv.
    """
    try:
        result = await forward_request("tl_pca", request, sampleid=sampleid, dtype=dtype)
        if result is not None:
            return result
        
        ads = ctx.request_context.lifespan_context
        adata = ads.get_adata(dtype=dtype, sampleid=sampleid)
        
        # Extract request parameters
        kwargs = request.model_dump()
        
        # Run PCA on CNV profiles
        cnv.tl.pca(adata, **kwargs)
        
        # Log the operation
        add_op_log(adata, "tl_pca", kwargs)
        
        return {
            "status": "success",
            "message": "PCA on CNV profiles computed successfully.",
            "n_comps": kwargs.get("n_comps", 50)
        }
    except KeyError as e:
        raise e
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e

@tl_mcp.tool()
async def umap(
    ctx: Context,
    request: UMAPModel,
    dtype: str = Field(default="exp", description="the datatype of anndata.X"),
    sampleid: str = Field(default=None, description="adata sampleid for analysis"),
):
    """
    Compute UMAP embedding of CNV profiles.
    This function should be called after running infercnv and neighbors computation.
    """
    try:
        result = await forward_request("tl_umap", request, sampleid=sampleid, dtype=dtype)
        if result is not None:
            return result
        
        ads = ctx.request_context.lifespan_context
        adata = ads.get_adata(dtype=dtype, sampleid=sampleid)
        
        # Extract request parameters
        kwargs = request.model_dump()
        
        # Run UMAP on CNV profiles
        cnv.tl.umap(adata, **kwargs)
        
        # Log the operation
        add_op_log(adata, "tl_umap", kwargs)
        
        return {
            "status": "success",
            "message": "UMAP embedding computed successfully.",
            "n_components": kwargs.get("n_components", 2)
        }
    except KeyError as e:
        raise e
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e

@tl_mcp.tool()
async def leiden(
    ctx: Context,
    request: LeidenModel,
    dtype: str = Field(default="exp", description="the datatype of anndata.X"),
    sampleid: str = Field(default=None, description="adata sampleid for analysis"),
):
    """
    Perform Leiden clustering on CNV profiles.
    This tool clusters cells based on their CNV profiles to identify potential tumor subclones.
    Should be called after running infercnv and neighbors computation.
    """
    try:
        result = await forward_request("tl_leiden", request, sampleid=sampleid, dtype=dtype)
        if result is not None:
            return result
        
        ads = ctx.request_context.lifespan_context
        adata = ads.get_adata(dtype=dtype, sampleid=sampleid)
        
        # Extract request parameters
        kwargs = request.model_dump()
        
        # Run Leiden clustering on CNV profiles
        cnv.tl.leiden(adata, **kwargs)
        
        # Log the operation
        add_op_log(adata, "tl_leiden", kwargs)
        
        # Get number of clusters
        key_added = kwargs.get("key_added", "cnv_leiden")
        if key_added in adata.obs:
            n_clusters = len(adata.obs[key_added].unique())
        else:
            n_clusters = "unknown"
        
        return {
            "status": "success",
            "message": "Leiden clustering on CNV profiles completed successfully.",
            "resolution": kwargs.get("resolution", 1.0),
            "n_clusters": n_clusters,
            "key_added": key_added
        }
    except KeyError as e:
        raise e
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e

@tl_mcp.tool()
async def cnv_score(
    ctx: Context,
    request: CNVScoreModel,
    dtype: str = Field(default="exp", description="the datatype of anndata.X"),
    sampleid: str = Field(default=None, description="adata sampleid for analysis"),
):
    """
    Calculate CNV scores for each cell.
    This tool quantifies the amount of CNV present in each cell, which can help identify tumor cells.
    Should be called after running infercnv.
    """
    try:
        result = await forward_request("tl_cnv_score", request, sampleid=sampleid, dtype=dtype)
        if result is not None:
            return result
        
        ads = ctx.request_context.lifespan_context
        adata = ads.get_adata(dtype=dtype, sampleid=sampleid)
        
        # Extract request parameters
        kwargs = request.model_dump()
        
        # Calculate CNV scores
        cnv.tl.cnv_score(adata, **kwargs)
        
        # Log the operation
        add_op_log(adata, "tl_cnv_score", kwargs)
        
        # Get statistics about the scores
        key_added = kwargs.get("key_added", "cnv_score")
        if key_added in adata.obs:
            min_score = float(adata.obs[key_added].min())
            max_score = float(adata.obs[key_added].max())
            mean_score = float(adata.obs[key_added].mean())
        else:
            min_score, max_score, mean_score = "unknown", "unknown", "unknown"
        
        return {
            "status": "success",
            "message": "CNV scores calculated successfully.",
            "method": kwargs.get("method", "mean_sqr"),
            "key_added": key_added,
            "min_score": min_score,
            "max_score": max_score,
            "mean_score": mean_score
        }
    except KeyError as e:
        raise e
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e 