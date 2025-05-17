"""
Tools schema definitions for infercnvpy-mcp.
"""

from pydantic import (
    Field,
    ValidationInfo,
    computed_field,
    field_validator,
    model_validator,
    BaseModel
)
from typing import Optional, Union, List, Dict, Any, Literal, Sequence

class InferCNVModel(BaseModel):
    """Input schema for the infercnv tool."""
    reference_key: str = Field(
        description="Key in adata.obs to find reference cells (assumed to be normal cells)."
    )
    
    reference_cat: List[str] = Field(
        description="Categories in adata.obs[reference_key] that contain reference (normal) cells."
    )
    
    window_size: int = Field(
        default=100,
        description="Size of the gene windows along the genome."
    )
    
    step: int = Field(
        default=None,
        description="Step size between windows. If None, step size is equal to window size."
    )

    
class PCPModel(BaseModel):
    """Input schema for PCA on CNV profiles."""
    n_comps: int = Field(
        default=50,
        description="Number of principal components to compute."
    )
    
    zero_center: bool = Field(
        default=True,
        description="If True, compute standard PCA from covariance matrix. If False, omit zero-centering variables, allowing to handle sparse input efficiently."
    )
    
    svd_solver: Literal['auto', 'full', 'arpack', 'randomized'] = Field(
        default='auto',
        description="SVD solver to use. Either 'auto', 'full', 'arpack', or 'randomized'."
    )
    
    random_state: Optional[int] = Field(
        default=0,
        description="Seed for random state."
    )
    
    n_iter: int = Field(
        default=5,
        description="Number of iterations for randomized PCA."
    )
    
    use_highly_variable: bool = Field(
        default=False,
        description="Whether to use highly variable regions for PCA."
    )
    
    copy: bool = Field(
        default=False,
        description="If true, return a copy instead of writing to adata."
    )

class UMAPModel(BaseModel):
    """Input schema for UMAP computation."""
    min_dist: float = Field(
        default=0.5,
        description="The effective minimum distance between embedded points."
    )
    
    spread: float = Field(
        default=1.0,
        description="The effective scale of embedded points."
    )
    
    n_components: int = Field(
        default=2,
        description="The number of dimensions of the embedding."
    )
    
    maxiter: Optional[int] = Field(
        default=None,
        description="The number of iterations (epochs) of the optimization."
    )
    
    alpha: float = Field(
        default=1.0,
        description="The initial learning rate for the embedding optimization."
    )
    
    gamma: float = Field(
        default=1.0,
        description="Weighting applied to negative samples in low dimensional embedding optimization."
    )
    
    negative_sample_rate: int = Field(
        default=5,
        description="The number of negative samples to select per positive sample in the optimization."
    )
    
    init_pos: Literal['spectral', 'random'] = Field(
        default='spectral',
        description="How to initialize the low dimensional embedding."
    )
    
    random_state: Optional[int] = Field(
        default=0,
        description="Seed for random state."
    )
    
    copy: bool = Field(
        default=False,
        description="Return a copy instead of writing to adata."
    )

class LeidenModel(BaseModel):
    """Input schema for Leiden clustering."""
    resolution: float = Field(
        default=1.0,
        description="Resolution parameter for the Leiden algorithm."
    )
    
    key_added: str = Field(
        default='cnv_leiden',
        description="Key under which to add the cluster labels."
    )
    
    directed: bool = Field(
        default=True,
        description="Whether to treat the graph as directed or undirected."
    )
    
    use_weights: bool = Field(
        default=True,
        description="Use weights from knn graph."
    )
    
    n_iterations: int = Field(
        default=-1,
        description="Number of iterations to run. If -1, run until convergence."
    )
    
    random_state: Optional[int] = Field(
        default=0,
        description="Seed for random state."
    )
    
    copy: bool = Field(
        default=False,
        description="Return a copy instead of writing to adata."
    )

class CNVScoreModel(BaseModel):
    """Input schema for CNV scoring."""
    layer: Optional[str] = Field(
        default=None,
        description="If provided, use this layer for CNV score calculation."
    )
    
    key_added: str = Field(
        default='cnv_score',
        description="Key under which to add the CNV score."
    )
    
    method: Literal['mean', 'var', 'mean_sqr'] = Field(
        default='mean_sqr',
        description="Method to use for CNV score calculation."
    )
    
    copy: bool = Field(
        default=False,
        description="Return a copy instead of writing to adata."
    ) 