import pytest

@pytest.fixture
def mcp():
    from infercnv_mcp.server import InferCNVMCPManager

    mcp = InferCNVMCPManager("infercnv-mcp").mcp
    return mcp
