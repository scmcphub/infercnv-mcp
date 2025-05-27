from fastmcp import FastMCP
from collections.abc import AsyncIterator
from contextlib import asynccontextmanager
from typing import Any


import scmcp_shared.server as shs

from scmcp_shared.util import filter_tools, setup_mcp
from .tl import tl_mcp
from .util import ul_mcp
from .pl import pl_mcp
from .pp import pp_mcp

ads = shs.AdataState()

@asynccontextmanager
async def adata_lifespan(server: FastMCP) -> AsyncIterator[Any]:
    yield ads


mcp = FastMCP("Infercnv-MCP-Server", lifespan=adata_lifespan)


 
module_dic = {
    "io": shs.io_mcp, 
    "tl": tl_mcp,
    "ul": ul_mcp,
    "pl": pl_mcp,
    "pp": pp_mcp
}
