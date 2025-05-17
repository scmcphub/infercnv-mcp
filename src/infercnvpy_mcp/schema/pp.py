"""
Preprocessing schema definitions for infercnvpy-mcp.
"""

from pydantic import (
    Field,
    ValidationInfo,
    computed_field,
    field_validator,
    model_validator,
    BaseModel
)
from typing import Optional, Union, List, Dict, Any, Literal

class NeighborsModel(BaseModel):
    """Input schema for neighbors computation."""
    n_neighbors: int = Field(
        default=15,
        description="Number of neighbors to use."
    )
    
    n_pcs: Optional[int] = Field(
        default=None,
        description="Number of PCs to use. If None, don't use PCA."
    )
    
    use_rep: Optional[str] = Field(
        default=None,
        description="Use the indicated representation. If None, use X."
    )
    
    knn: bool = Field(
        default=True,
        description="If True, use a hard threshold to restrict the number of neighbors to n_neighbors, that is, consider a knn graph. Otherwise, use a Gaussian kernel to assign low weights to neighbors more distant than the n_neighbors nearest neighbors."
    )
    
    random_state: Optional[int] = Field(
        default=0,
        description="Random seed."
    )
    
    method: Literal['umap', 'gauss', 'rapids'] = Field(
        default='umap',
        description="Use 'umap', 'gauss', or 'rapids' (if installed) for computing connectivities."
    )
    
    metric: Literal['euclidean', 'manhattan', 'chebyshev', 'minkowski', 'canberra', 'braycurtis'] = Field(
        default='euclidean',
        description="Distance metric to use for computing connectivities."
    )
    
    copy: bool = Field(
        default=False,
        description="Return a copy instead of writing to adata."
    )
    
    @field_validator('method')
    def validate_method(cls, v: str) -> str:
        if v not in ['umap', 'gauss', 'rapids']:
            raise ValueError(f"Method must be 'umap', 'gauss', or 'rapids', got {v}")
        return v
    
    @field_validator('metric')
    def validate_metric(cls, v: str) -> str:
        valid_metrics = ['euclidean', 'manhattan', 'chebyshev', 'minkowski', 'canberra', 'braycurtis']
        if v not in valid_metrics:
            raise ValueError(f"Metric must be one of {valid_metrics}, got {v}")
        return v 