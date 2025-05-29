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
from scmcp_shared.util import add_op_log,  filter_args, get_ads, obsm2adata,forward_request
from scmcp_shared.logging_config import setup_logger
from scmcp_shared.schema import AdataInfo
import infercnvpy as cnv
from scmcp_shared.util import update_mcp_args
from scmcp_shared.server.tl import ScanpyToolsMCP

include_tools=["tsne", "umap", "leiden", "louvain", "pca"]

tl_mcp = ScanpyToolsMCP(
    include_tools=include_tools,
    AdataInfo=CNVAdataInfo
).mcp



@tl_mcp.tool()
def infercnv(
    request: InferCNVParam,
    adinfo: AdataInfo = AdataInfo()
):
    """Infer Copy Number Variation (CNV) by averaging gene expression over genomic regions.
    """
    logger = setup_logger()
    try:
        result = forward_request("infercnv", request, adinfo)
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
def cnv_score(
    request: CNVScoreParam,
    adinfo: CNVAdataInfo = CNVAdataInfo()
):
    """Calculate CNV scores for each cell."""
    try:
        result = forward_request("tl_cnv_score", request, adinfo)
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
            {"sampleid": adata.uns["exp_sampleid"], "adtype": "exp", "adata": exp_adata},
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
        },
        "leiden": {
            "neighbors_key": {"default": "cnv_neighbors"}, 
            "key_added": {"default": "cnv_leiden"},
        },
        "umap": {
            "neighbors_key": {"default": "cnv_neighbors"},
            "key_added": {"default": "cnv_umap"},
        },
        "tsne": {   
            "neighbors_key": {"default": "cnv_neighbors"},
            "key_added": {"default": "cnv_tsne"},
        }
    }
)
