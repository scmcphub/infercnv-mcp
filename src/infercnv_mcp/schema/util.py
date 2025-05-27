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


class GenePosParam(BaseModel):
    file: str = Field(
        ...,
        description="Path to the gene position file."
    )
    sep: str = Field(
        default=",",
        description="Separator of the gene position file. csv(,) or tsv('\\t')"
    )