# Infercnvpy-MCP

Natural language interface for Copy Number Variation (CNV) inference from scRNA-Seq data with infercnvpy through MCP.

## 🪩 What can it do?

- IO module for reading and writing scRNA-Seq data
- Preprocessing module for neighbors computation and data preparation
- Tool module for CNV inference, clustering, and dimensionality reduction
- Plotting module for chromosome heatmaps, UMAP, and t-SNE visualizations

## ❓ Who is this for?

- Researchers who want to infer CNVs from scRNA-Seq data using natural language
- Agent developers who want to integrate CNV analysis into their applications

## 🌐 Where to use it?

You can use infercnvpy-mcp in most AI clients, plugins, or agent frameworks that support the MCP:

- AI clients, like Cherry Studio
- Plugins, like Cline
- Agent frameworks, like Agno 

## 🏎️ Quickstart

### Install

Install from PyPI
```
pip install infercnvpy-mcp
```
you can test it by running
```
infercnvpy-mcp run
```

#### Run infercnvpy-server locally
Refer to the following configuration in your MCP client:

```
"mcpServers": {
  "infercnvpy-mcp": {
    "command": "infercnvpy-mcp",
    "args": [
      "run"
    ]
  }
}
```

#### Run infercnvpy-server remotely
Refer to the following configuration in your MCP client:

Run it in your server
```
infercnvpy-mcp run --transport shttp --port 8000
```

Then configure your MCP client, like this:
```
http://localhost:8000/mcp
```

## 🤝 Contributing

Contributions to the code are welcome! 