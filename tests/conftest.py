import pytest

@pytest.fixture
def mcp():
    from infercnv_mcp.server import mcp, module_dic 
    from scmcp_shared.util import setup_mcp

    mcp = setup_mcp(mcp, module_dic)
    return mcp
