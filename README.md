# Infercnv-MCP

Natural language interface for Copy Number Variation (CNV) inference from scRNA-Seq data with infercnvpy through MCP.

## 🪩 What can it do?

- IO module for reading and writing scRNA-Seq data, load gene position 
- Preprocessing module for neighbors computation and data preparation
- Tool module for CNV inference, cnv score
- Plotting module for chromosome heatmaps, UMAP, and t-SNE visualizations

## ❓ Who is this for?

- Researchers who want to infer CNVs from scRNA-Seq data using natural language
- Agent developers who want to integrate CNV analysis into their applications

## 🌐 Where to use it?

You can use infercnv-mcp in most AI clients, plugins, or agent frameworks that support the MCP:

- AI clients, like Cherry Studio
- Plugins, like Cline
- Agent frameworks, like Agno 

## 📚 Documentation

scmcphub's complete documentation is available at https://docs.scmcphub.org

## 🏎️ Quickstart

### Install

Install from PyPI
```
pip install infercnv-mcp
```
you can test it by running
```
infercnv-mcp run
```

#### run infercnv-mcp locally
Refer to the following configuration in your MCP client:

check path
```
$ which infercnv 
/home/test/bin/infercnv-mcp
```

```
"mcpServers": {
  "infercnv-mcp": {
    "command": "/home/test/bin/infercnv-mcp",
    "args": [
      "run"
    ]
  }
}
```

#### Run infercnv-server remotely
Refer to the following configuration in your MCP client:

Run it in your server
```
infercnv-mcp run --transport shttp --port 8000
```

Then configure your MCP client, like this:
```
http://localhost:8000/mcp
```

## 🤝 Contributing

If you have any questions, welcome to submit an issue, or contact me(hsh-me@outlook.com). Contributions to the code are also welcome!

## Citing
If you use infercnv-mcp in your research, please consider citing following work: 
> https://github.com/icbi-lab/infercnvpy