"""
Tool functions for infercnvpy-mcp.
"""
import asyncio
from fastmcp import FastMCP, Context
from fastmcp.exceptions import ToolError
from pathlib import Path
import os
from ..schema.tl import *
from ..schema import CNVAdataInfo
from scmcp_shared.util import add_op_log, savefig, filter_args, forward_request, get_ads, obsm2adata
from scmcp_shared.logging_config import setup_logger
from scmcp_shared.schema import AdataModel
import infercnvpy as cnv
from scmcp_shared.util import update_mcp_args
from scmcp_shared.server.tl import ScanpyToolsMCP

include_tools=["tsne", "umap", "leiden", "louvain", "pca"]

tl_mcp = ScanpyToolsMCP(
    include_tools=include_tools
).mcp



@tl_mcp.tool()
async def infercnv(
    request: InferCNVParam,
    adinfo: AdataModel = AdataModel()
):
    """Infer Copy Number Variation (CNV) by averaging gene expression over genomic regions.
    """
    logger = setup_logger()
    try:
        result = await forward_request("infercnv", request, adinfo)
        if result is not None:
            return result
        
        ads = get_ads()
        adata = ads.get_adata(adinfo=adinfo)
        # Validate required columns
        required_cols = ['chromosome', 'start', 'end']
        if not all(col in adata.var.columns for col in required_cols):
            msg = f"adata.var must contain columns: {required_cols}, please run load_gene_position first."
            raise ToolError(msg)

        func_kwargs = filter_args(request, cnv.tl.infercnv)
        logger.info(f"func_kwargs: {func_kwargs} 7777")
        
        cnv.tl.infercnv(adata, **func_kwargs)
        cnv_adata = obsm2adata(adata, "X_cnv")
        
        # Log the operation
        add_op_log(adata, "tl_infercnv", func_kwargs, adinfo)
        cnv_adata.uns["exp_sampleid"] = adinfo.sampleid or ads.active_id
        ads.set_adata(cnv_adata, sampleid="adata_cnv", sdtype="cnv")
        return [
            {"sampleid": adinfo.sampleid or ads.active_id, "adtype": adinfo.adtype, "adata": adata},
            {"sampleid": "adata_cnv", "adtype": "cnv", "adata": cnv_adata}
        ]
    except ToolError as e:
        raise ToolError(e)
    except Exception as e:
        raise ToolError(e)


@tl_mcp.tool()
async def cnv_score(
    request: CNVScoreParam,
    adinfo: CNVAdataInfo = CNVAdataInfo()
):
    """Calculate CNV scores for each cell."""
    try:
        result = await forward_request("tl_cnv_score", request, adinfo)
        if result is not None:
            return result
        ads = get_ads()
        adata = ads.get_adata(adinfo=adinfo)
        func_kwargs = filter_args(request, cnv.tl.cnv_score)
        cnv.tl.cnv_score(adata, **func_kwargs)
        add_op_log(adata, "tl_cnv_score", func_kwargs, adinfo)
        exp_adata = ads.get_adata(sampleid=adata.uns["exp_sampleid"], adtype="exp")
        exp_adata.obs[request.key_added] = adata.obs[request.key_added]
        return [
            {"sampleid": adinfo.sampleid or ads.active_id, "adtype": adinfo.adtype, "adata": adata}
        ]
    except ToolError as e:
        raise ToolError(e)
    except Exception as e:
        raise ToolError(e)
 

update_mcp_args(
    tl_mcp, {
        "pca": {
            "zero_center": {"default": "false"}, 
            "svd_solver": {"default": "arpack"}, 
            "key_added": {"default": "cnv_pca"},
            "sampleid": {"default": "adata_cnv"},
            "adtype": {"default": "cnv"},
        },
        "leiden": {
            "neighbors_key": {"default": "cnv_neighbors"}, 
            "key_added": {"default": "cnv_leiden"},
            "sampleid": {"default": "adata_cnv"}, 
            "adtype": {"default": "cnv"},
        },
        "umap": {
            "neighbors_key": {"default": "cnv_neighbors"},
            "key_added": {"default": "cnv_umap"},
            "sampleid": {"default": "adata_cnv"},
            "adtype": {"default": "cnv"},
        },
        "tsne": {   
            "neighbors_key": {"default": "cnv_neighbors"},
            "key_added": {"default": "cnv_tsne"},
            "sampleid": {"default": "adata_cnv"},
            "adtype": {"default": "cnv"},
        }
    }
)


# @tl_mcp.tool()
# async def cnv_score(
#     ctx: Context,
#     request: CNVScoreModel,
#     dtype: str = Field(default="exp", description="the datatype of anndata.X"),
#     sampleid: str = Field(default=None, description="adata sampleid for analysis"),
# ):
#     """
#     Calculate CNV scores for each cell.
#     This tool quantifies the amount of CNV present in each cell, which can help identify tumor cells.
#     Should be called after running infercnv.
#     """
#     try:
#         result = await forward_request("tl_cnv_score", request, sampleid=sampleid, dtype=dtype)
#         if result is not None:
#             return result
        
#         ads = ctx.request_context.lifespan_context
#         adata = ads.get_adata(dtype=dtype, sampleid=sampleid)
        
#         # Extract request parameters
#         kwargs = request.model_dump()
        
#         # Calculate CNV scores
#         cnv.tl.cnv_score(adata, **kwargs)
        
#         # Log the operation
#         add_op_log(adata, "tl_cnv_score", kwargs)
        
#         # Get statistics about the scores
#         key_added = kwargs.get("key_added", "cnv_score")
#         if key_added in adata.obs:
#             min_score = float(adata.obs[key_added].min())
#             max_score = float(adata.obs[key_added].max())
#             mean_score = float(adata.obs[key_added].mean())
#         else:
#             min_score, max_score, mean_score = "unknown", "unknown", "unknown"
        
#         return {
#             "status": "success",
#             "message": "CNV scores calculated successfully.",
#             "method": kwargs.get("method", "mean_sqr"),
#             "key_added": key_added,
#             "min_score": min_score,
#             "max_score": max_score,
#             "mean_score": mean_score
#         }
#     except KeyError as e:
#         raise e
#     except Exception as e:
#         if hasattr(e, '__context__') and e.__context__:
#             raise Exception(f"{str(e.__context__)}")
#         else:
#             raise e 