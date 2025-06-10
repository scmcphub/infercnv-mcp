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


class InferCNVParam(BaseModel):
    """Input schema for the infercnv tool."""
    reference_key: Optional[str] = Field(
        default=None,
        description="Column name in adata.obs that contains tumor/normal cell cluster annotations."
    )
    reference_cat: Optional[List[str]] = Field(
        default=None,
        description="One or multiple values in adata.obs[reference_key] that annotate normal cells."
    )
    # reference: Optional[Any] = Field(
    #     default=None,
    #     description="Directly supply an array of average normal gene expression. Overrides reference_key and reference_cat."
    # )
    lfc_clip: float = Field(
        default=3.0,
        description="Clip log fold changes at this value"
    )
    window_size: int = Field(
        default=100,
        description="Size of the running window (number of genes to include in the window)"
    )
    
    step: int = Field(
        default=10,
        description="Only compute every nth running window where n = step. Set to 1 to compute all windows."
    )
    
    dynamic_threshold: Optional[float] = Field(
        default=1.5,
        description="Values < dynamic threshold * STDDEV will be set to 0, where STDDEV is the standard deviation of the smoothed gene expression. Set to None to disable this step."
    )
    
    exclude_chromosomes: Optional[Sequence[str]] = Field(
        default=('chrX', 'chrY'),
        description="List of chromosomes to exclude. The default is to exclude genosomes."
    )
    
    chunksize: int = Field(
        default=5000,
        description="Process dataset in chunks of cells. This allows to run infercnv on datasets with many cells, where the dense matrix would not fit into memory."
    )
    
    n_jobs: Optional[int] = Field(
        default=None,
        description="Number of jobs for parallel processing. Default: use all cores. Data will be submitted to workers in chunks, see chunksize."
    )
    
    layer: Optional[str] = Field(
        default=None,
        description="Layer from adata to use. If None, use X."
    )
    
    calculate_gene_values: bool = Field(
        default=False,
        description="If True per gene CNVs will be calculated and stored in adata.layers['gene_values_cnv']."
    )

class CNVScoreParam(BaseModel):

    groupby: str = Field(
        default='cnv_leiden',
        description="Key under which the clustering is stored in adata.obs. Usually the result of infercnvpy.tl.leiden(), but could also be other grouping information, e.g. sample or patient information."
    )
    key_added: str = Field(
        default='cnv_score',
        description="Key under which the score will be stored in adata.obs."
    )
