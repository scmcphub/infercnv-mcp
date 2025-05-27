"""
Command-line interface for infercnv-mcp.

This module provides a CLI entry point for the infercnv-mcp package.
"""

from scmcp_shared.cli import MCPCLI

cli = MCPCLI(name="infercnv-mcp", help_text="InferCNV MCP Server CLI")

def run_cli():
    from .server import mcp, module_dic
    cli.run_cli(mcp, module_dic)
