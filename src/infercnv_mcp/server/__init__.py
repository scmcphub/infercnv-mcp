from scmcp_shared.server import BaseMCPManager
from scmcp_shared.server import io_mcp
from .tl import tl_mcp
from .util import ul_mcp
from .pl import pl_mcp
from .pp import pp_mcp


class InferCNVMCPManager(BaseMCPManager):
    """Manager class for Scanpy MCP modules."""
    
    def _init_modules(self):
        """Initialize available Scanpy MCP modules."""
        self.available_modules = {
            "io": io_mcp,
            "pp": pp_mcp,
            "tl": tl_mcp,
            "pl": pl_mcp,
            "ul": ul_mcp
        }
