"""
Plotting schema definitions for infercnvpy-mcp.
"""

from pydantic import Field,BaseModel
from typing import Optional, Tuple


class ChromosomeHeatmapParam(BaseModel):
    """Input schema for chromosome heatmap plotting."""
    groupby: Optional[str] = Field(
        default="cnv_leiden",
        description="Key in adata.obs for grouping cells."
    )
    use_rep: str = Field(default="cnv", description="The representation to use for plotting.")
    cmap: str = Field(default="bwr", description="The colormap to use for plotting.")
    figsize: Optional[Tuple[float, float]] = Field(default=(16, 10), description="Figure size (width, height) in inches.")
