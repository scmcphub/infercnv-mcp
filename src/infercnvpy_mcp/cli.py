"""
Command-line interface for infercnvpy-mcp.

This module provides a CLI entry point for the infercnvpy-mcp package.
"""

import asyncio
import os
import sys
import typer
from enum import Enum
from typing import Optional

app = typer.Typer(
    name="infercnvpymcp",
    help="Infercnvpy MCP Server CLI",
    add_completion=False,
    no_args_is_help=True,  # Show help if no args provided    
)

# Define enums for choices
class Module(str, Enum):
    ALL = "all"
    IO = "io"
    PP = "pp"
    PL = "pl"
    TL = "tl"
    UTIL = "util"

class Transport(str, Enum):
    STDIO = "stdio"
    SSE = "sse"
    SHTTP = "shttp"

@app.command(name="run")
def run(
    log_file: Optional[str] = typer.Option(None, "--log-file", help="Log file path, use stdout if None"),
    module: Module = typer.Option(Module.ALL, "-m", "--module", help="Specify modules to load", 
                              case_sensitive=False),
    transport: Transport = typer.Option(Transport.STDIO, "-t", "--transport", help="Specify transport type", 
                                 case_sensitive=False),
    port: int = typer.Option(8000, "-p", "--port", help="transport port"),
    host: str = typer.Option("127.0.0.1", "--host", help="transport host"),
    forward: str = typer.Option(None, "-f", "--forward", help="forward request to another server"),
):
    """Start Infercnvpy MCP Server"""
    
    # Set environment variables
    if log_file is not None:
        os.environ['SCMCP_LOG_FILE'] = log_file
    if forward is not None:
        os.environ['SCMCP_FORWARD'] = forward
    os.environ['SCMCP_TRANSPORT'] = transport.value
    os.environ['SCMCP_HOST'] = host
    os.environ['SCMCP_PORT'] = str(port)
    os.environ['SCMCP_MODULE'] = module.value
    from .server import infercnvpy_mcp, setup
    asyncio.run(setup())
    if transport == Transport.STDIO:
        infercnvpy_mcp.run()
    elif transport == Transport.SSE:
        from .util import get_figure
        from starlette.routing import Route

        infercnvpy_mcp._additional_http_routes = [Route("/figures/{figure_name}", endpoint=get_figure)]        
        infercnvpy_mcp.run(
                transport="sse",
                host=host, 
                port=port, 
                log_level="info"
            )

    elif transport == Transport.SHTTP:
        from .util import get_figure
        from starlette.routing import Route

        infercnvpy_mcp._additional_http_routes = [Route("/figures/{figure_name}", endpoint=get_figure)]
        infercnvpy_mcp.run(
                transport="streamable-http",
                host=host, 
                port=port, 
                log_level="info"
            )

@app.callback()
def main():
    """Infercnvpy MCP CLI root command."""
    pass 