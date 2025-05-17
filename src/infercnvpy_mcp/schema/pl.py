"""
Plotting schema definitions for infercnvpy-mcp.
"""

from pydantic import (
    Field,
    ValidationInfo,
    computed_field,
    field_validator,
    model_validator,
    BaseModel
)
from typing import Optional, Union, List, Dict, Any, Literal, Tuple, Sequence

class ChromosomeHeatmapModel(BaseModel):
    """Input schema for chromosome heatmap plotting."""
    groupby: Optional[str] = Field(
        default=None,
        description="Key in adata.obs for grouping cells."
    )
    
    dendrogram: bool = Field(
        default=False,
        description="If True, compute a dendrogram based on the mean expression of each group."
    )
    
    
    figsize: Optional[Tuple[float, float]] = Field(
        default=None,
        description="Figure size (width, height) in inches."
    )
    
    show: bool = Field(
        default=True,
        description="If True, display the figure."
    )
    
    save: Optional[str] = Field(
        default=None,
        description="If provided, save the figure to this file path."
    )

class UMAPPlotModel(BaseModel):
    """Input schema for UMAP plotting."""
    color: Union[str, List[str]] = Field(
        default=None,
        description="Key(s) in adata.obs for coloring cells."
    )
    
    use_raw: bool = Field(
        default=False,
        description="Use raw attribute of adata if present."
    )
    
    layer: Optional[str] = Field(
        default=None,
        description="Layer in adata.layers to use."
    )
    
    sort_order: bool = Field(
        default=True,
        description="If True, sort unordered categorical in ascending order."
    )
    
    palette: Optional[Union[str, List[str], Dict[str, str]]] = Field(
        default=None,
        description="Colors to use for plotting categorical annotations."
    )
    
    size: float = Field(
        default=None,
        description="Point size for UMAP plot."
    )
    
    legend_loc: Literal['right margin', 'on data'] = Field(
        default='right margin',
        description="Location of legend, either 'on data' or 'right margin'."
    )
    
    legend_fontsize: Optional[float] = Field(
        default=None,
        description="Legend font size."
    )
    
    legend_fontoutline: float = Field(
        default=0,
        description="Width of the white line around legend text."
    )
    
    colorbar_loc: Optional[str] = Field(
        default=None,
        description="Location of colorbar, e.g., 'right', or None."
    )
    
    frameon: bool = Field(
        default=True,
        description="Draw a frame around the UMAP scatter plot."
    )
    
    ncols: int = Field(
        default=4,
        description="Number of columns for panel plots."
    )
    
    wspace: Optional[float] = Field(
        default=None,
        description="Width space between panels."
    )
    
    hspace: Optional[float] = Field(
        default=None,
        description="Height space between panels."
    )
    
    return_fig: bool = Field(
        default=False,
        description="Return the figure object."
    )
    
    show: bool = Field(
        default=True,
        description="If True, display the figure."
    )
    
    save: Optional[str] = Field(
        default=None,
        description="If provided, save the figure to this file path."
    )
    
    title: Optional[str] = Field(
        default=None,
        description="Title for the figure."
    )
    
    figsize: Optional[Tuple[float, float]] = Field(
        default=None,
        description="Figure size (width, height) in inches."
    ) 