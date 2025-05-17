"""
Utility functions for infercnvpy-mcp.
"""

import os
import tempfile
from pathlib import Path
from starlette.responses import FileResponse
from starlette.requests import Request
import matplotlib.pyplot as plt

# Global dictionary to store figures
FIGURES = {}

def save_figure(fig, name=None):
    """Save a matplotlib figure for later retrieval."""
    if name is None:
        name = f"figure_{len(FIGURES)}.png"
    
    # Create temp directory if doesn't exist
    temp_dir = Path(tempfile.gettempdir()) / "infercnvpy_mcp_figures"
    temp_dir.mkdir(exist_ok=True)
    
    # Save figure to file
    filename = temp_dir / name
    fig.savefig(filename)
    
    # Store in global dictionary
    FIGURES[name] = str(filename)
    
    return name

async def get_figure(request: Request):
    """Retrieve a saved figure."""
    figure_name = request.path_params["figure_name"]
    
    if figure_name in FIGURES:
        return FileResponse(FIGURES[figure_name])
    else:
        return {"error": f"Figure {figure_name} not found"} 