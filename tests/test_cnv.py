import pytest
from fastmcp import Client
from pathlib import Path
import nest_asyncio

# Apply nest_asyncio at module level
nest_asyncio.apply()


@pytest.mark.asyncio 
async def test_activity(mcp):
    # Pass the server directly to the Client constructor
    testfile = Path(__file__).parent / "data/maynard2020_3k.h5ad"
    async with Client(mcp) as client:
        result = await client.call_tool("io_read", {"request":{"filename": testfile}})
        assert "AnnData" in result[0].text

        # result = await client.call_tool(
        #     "ul_load_gene_position", 
        #     {"request":{"file": Path(__file__).parent / "data/gene_pos.csv"}}
        # )
        # assert "success" in result[0].text

        # result = await client.call_tool(
        #     "tl_infercnv",  {'request': {'reference_key': 'cell_type',   
        #                 'reference_cat': ['B cell', 'Macrophage', 'Mast cell', 'Monocyte', 'NK cell', 'Plasma cell', 'T cell CD4', 'T cell CD8', 'T cell regulatory', 'mDC', 'pDC']}}
        # )
        # assert "adata_cnv" in result[0].text
        # result = await client.call_tool(
        #     "tl_pca",  {'request': {'n_comps': 50, 'layer': None, 'zero_center': 'false', 'svd_solver': 'arpack',        
        #     'mask_var': None, 'chunked': False, 'chunk_size': None}, 'adinfo': {'sampleid': 'adata_cnv', 'adtype': 'cnv'}}
        # )
        # assert "adata_cnv" in result[0].text

        # result = await client.call_tool(
        #     "pp_neighbors",  {'request': {'n_neighbors': 15, 'n_pcs': None, 'use_rep': 'X_pca', 'knn': True, 'method': 'umap', 'transformer': None,    
        #     'metric': 'euclidean', 'metric_kwds': {}, 'random_state': 0, 'key_added': 'cnv_neighbors'}, 'adinfo': {'sampleid': 'adata_cnv', 'adtype':
        #     'cnv'}}
        # )
        # assert "adata_cnv" in result[0].text

        # result = await client.call_tool(
        #     "tl_leiden",  {'request': {'key_added': 'cnv_leiden',    
        #     'neighbors_key': 'cnv_neighbors',  'flavor': 'igraph'}, 'adinfo': {'sampleid': 'adata_cnv',        
        #     'adtype': 'cnv'}}
        # )
        # assert "adata_cnv" in result[0].text

