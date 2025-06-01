"""
Command-line interface for infercnv-mcp.

This module provides a CLI entry point for the infercnv-mcp package.
"""

from scmcp_shared.cli import MCPCLI
from .server import InferCNVMCPManager

cli = MCPCLI(
    name="infercnv-mcp", 
    help_text="InferCNV MCP Server CLI",
    manager=InferCNVMCPManager
)
