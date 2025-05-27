"""Schema definitions for infercnvpy-mcp.""" 

from pydantic import Field, BaseModel,ConfigDict


class CNVAdataInfo(BaseModel):
    """Input schema for the adata tool."""
    sampleid: str = Field(default="adata_cnv", description="adata sampleid")
    adtype: str = Field(default="cnv", description="The input adata.X data type")

    model_config = ConfigDict(
        extra="ignore"
    )
