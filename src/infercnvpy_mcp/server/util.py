"""
Utility functions for infercnvpy-mcp server.
"""

import scanpy as sc
import infercnvpy as cnv
from fastmcp import FastMCP, Context
from pydantic import Field
from ..schema.util import CNVStatusModel
from scmcp_shared.util import add_op_log, forward_request
from scmcp_shared.logging_config import setup_logger

logger = setup_logger()

ul_mcp = FastMCP("InfercnvpyMCP-Util-Server")

@ul_mcp.tool()
async def assign_cnv_status(
    ctx: Context,
    request: CNVStatusModel,
    dtype: str = Field(default="exp", description="the datatype of anndata.X"),
    sampleid: str = Field(default=None, description="adata sampleid for analysis"),
):
    """
    Assign cells to tumor or normal categories based on their CNV clusters.
    This function should be called after running infercnv, neighbors, and leiden clustering.
    It helps to identify which cells are likely tumor cells based on their CNV profiles.
    """
    try:
        result = await forward_request("ul_assign_cnv_status", request, sampleid=sampleid, dtype=dtype)
        if result is not None:
            return result
        
        ads = ctx.request_context.lifespan_context
        adata = ads.get_adata(dtype=dtype, sampleid=sampleid)
        
        # Extract request parameters
        kwargs = request.model_dump()
        tumor_clusters = kwargs.get("tumor_clusters", [])
        normal_clusters = kwargs.get("normal_clusters", None)
        cluster_key = kwargs.get("cluster_key", "cnv_leiden")
        status_key = kwargs.get("status_key", "cnv_status")
        
        # Create status column
        adata.obs[status_key] = "normal"
        
        # Mark tumor cells
        tumor_mask = adata.obs[cluster_key].isin(tumor_clusters)
        adata.obs.loc[tumor_mask, status_key] = "tumor"
        
        # Mark normal cells if specified
        if normal_clusters is not None:
            normal_mask = adata.obs[cluster_key].isin(normal_clusters)
            adata.obs.loc[normal_mask, status_key] = "normal"
        
        # Log the operation
        add_op_log(adata, "assign_cnv_status", kwargs)
        
        # Count cells in each category
        tumor_count = sum(adata.obs[status_key] == "tumor") 
        normal_count = sum(adata.obs[status_key] == "normal")
        
        return {
            "status": "success",
            "message": "CNV status assigned successfully.",
            "tumor_count": int(tumor_count),
            "normal_count": int(normal_count),
            "status_key": status_key
        }
    except KeyError as e:
        raise e
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e

@ul_mcp.tool()
async def list_available_datasets(
    ctx: Context,
):
    """
    List available example datasets from infercnvpy.
    This function can be called to see what datasets are available for testing and demonstration.
    """
    try:
        # These are the currently available datasets in infercnvpy
        available_datasets = {
            "maynard2020_3k": "Lung cancer dataset with 3,000 cells (small version)",
            "oligodendroglioma": "Oligodendroglioma tumor dataset"
        }
        
        return {
            "status": "success",
            "message": "Available datasets retrieved successfully.",
            "datasets": available_datasets
        }
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e

@ul_mcp.tool()
async def load_example_dataset(
    ctx: Context,
    dataset_name: str = Field(default="maynard2020_3k", description="Name of the dataset to load"),
    sampleid: str = Field(default=None, description="adata sampleid for the loaded dataset"),
    dtype: str = Field(default="exp", description="the datatype of anndata.X"),
):
    """
    Load an example dataset from infercnvpy.
    This function loads one of the pre-packaged datasets for testing and demonstration purposes.
    """
    try:
        ads = ctx.request_context.lifespan_context
        
        # Set sampleid if not provided
        if sampleid is None:
            sampleid = dataset_name
            
        # Load the appropriate dataset
        if dataset_name == "maynard2020_3k":
            adata = cnv.datasets.maynard2020_3k()
        elif dataset_name == "oligodendroglioma":
            adata = cnv.datasets.oligodendroglioma()
        else:
            raise ValueError(f"Dataset {dataset_name} not found. Use list_available_datasets to see available options.")
        
        # Store in AdataState
        ads.active_id = sampleid
        ads.set_adata(adata, sampleid=sampleid, sdtype=dtype)
        
        # Log basic information about the dataset
        n_obs = adata.n_obs
        n_vars = adata.n_vars
        obs_keys = list(adata.obs.columns)
        var_keys = list(adata.var.columns)
        
        return {
            "status": "success",
            "message": f"Dataset {dataset_name} loaded successfully.",
            "sampleid": sampleid,
            "n_cells": n_obs,
            "n_genes": n_vars,
            "obs_keys": obs_keys[:10] if len(obs_keys) > 10 else obs_keys,  # Limit to first 10 for readability
            "var_keys": var_keys[:10] if len(var_keys) > 10 else var_keys,  # Limit to first 10 for readability
        }
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e 