"""
Utility schema definitions for infercnvpy-mcp.
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

class CNVStatusModel(BaseModel):
    """Input schema for CNV status assignment."""
    tumor_clusters: List[str] = Field(
        description="List of cluster IDs to mark as tumor cells."
    )
    
    normal_clusters: Optional[List[str]] = Field(
        default=None,
        description="List of cluster IDs to mark as normal cells. If None, all non-tumor clusters are marked as normal."
    )
    
    cluster_key: str = Field(
        default="cnv_leiden",
        description="Key in adata.obs containing cluster labels."
    )
    
    status_key: str = Field(
        default="cnv_status",
        description="Key to add to adata.obs with tumor/normal status."
    ) 