"""
Preprocessing functions for infercnvpy-mcp.
"""

import scanpy as sc
import infercnvpy as cnv
from fastmcp import FastMCP, Context
from pydantic import Field
from ..schema.pp import NeighborsModel
from scmcp_shared.util import add_op_log, forward_request
from scmcp_shared.logging_config import setup_logger

logger = setup_logger()

pp_mcp = FastMCP("InfercnvpyMCP-PP-Server")

@pp_mcp.tool()
async def neighbors(
    ctx: Context,
    request: NeighborsModel,
    dtype: str = Field(default="exp", description="the datatype of anndata.X"),
    sampleid: str = Field(default=None, description="adata sampleid for analysis"),
):
    """
    Compute a neighborhood graph for cells using their CNV profiles.
    This function is similar to scanpy.pp.neighbors but works on CNV profiles.
    It should be called after running infercnv and before clustering or dimensionality reduction.
    """
    try:
        result = await forward_request("pp_neighbors", request, sampleid=sampleid, dtype=dtype)
        if result is not None:
            return result
        
        ads = ctx.request_context.lifespan_context
        adata = ads.get_adata(dtype=dtype, sampleid=sampleid)
        
        # Extract request parameters
        kwargs = request.model_dump()
        
        # Run neighbors computation
        cnv.pp.neighbors(adata, **kwargs)
        
        # Log the operation
        add_op_log(adata, "pp_neighbors", kwargs)
        
        return {
            "status": "success",
            "message": "Neighborhood graph computed successfully.",
            "n_neighbors": kwargs.get("n_neighbors", 15)
        }
    except KeyError as e:
        raise e
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e 