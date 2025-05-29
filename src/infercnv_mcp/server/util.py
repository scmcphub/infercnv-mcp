"""
Utility functions for infercnv-mcp server.
"""
from fastmcp import FastMCP, Context
from fastmcp.exceptions import ToolError
from pathlib import Path
import os
from ..schema.util import *
from scmcp_shared.util import add_op_log, savefig, filter_args, forward_request, get_ads
from scmcp_shared.logging_config import setup_logger
from scmcp_shared.schema import AdataInfo
import pandas as pd
from scmcp_shared.server import ScanpyUtilMCP

ul_mcp = ScanpyUtilMCP(
    include_tools=["check_samples"],
).mcp


@ul_mcp.tool()
def load_gene_position(
    request: GenePosParam,
    adinfo: AdataInfo = AdataInfo()
):
    """Load gene position file and add to adata.var when finished reading adata."""
    logger = setup_logger()
    try:
        result = forward_request("load_gene_position", request, adinfo)
        if result is not None:
            return result
        adata = get_ads().get_adata(adinfo=adinfo)
        
        # Read gene position file
        gene_pos = pd.read_csv(request.file, sep=request.sep, index_col=0)
        
        # Validate required columns
        required_cols = ['chromosome', 'start', 'end']
        if not all(col in gene_pos.columns for col in required_cols):
            msg = f"Gene position file must contain columns: {required_cols}"
            raise ToolError(msg)
            
        # Check if all genes in adata.var_names exist in gene_pos
        missing_genes = set(adata.var_names) - set(gene_pos.index)
        msg_ls = []
        if missing_genes:
            msg1 = f"Some genes in adata.var_names are missing from gene position file: {missing_genes[:10]}..."
            msg_ls.append(msg1)
        # Check if all genes in gene_pos exist in adata.var_names
        extra_genes = set(gene_pos.index) - set(adata.var_names)
        if extra_genes:
            msg2 = f"Some genes in gene position file are not in adata.var_names: {extra_genes[:10]}..."
            msg_ls.append(msg2)
        # Add chromosome position information to adata.var
        for col in ['chromosome', 'start', 'end']:
            # Only add information for genes that exist in both adata.var_names and gene_pos
            common_genes = list(set(adata.var_names) & set(gene_pos.index))
            adata.var.loc[common_genes, col] = gene_pos.loc[common_genes, col]
        return {
            "msg": f"Gene position file loaded successfully.",
            "notice": msg_ls
        }
    except ToolError as e:
        raise ToolError(e)
    except Exception as e:
        raise ToolError(e)
