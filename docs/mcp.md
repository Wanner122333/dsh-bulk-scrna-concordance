# MCP setup

`python/mcp_server.py` exposes the `bulk_scrna_concordance` tool using the
official Python MCP SDK. It supports stdio and Streamable HTTP. The analysis
parameters and interpretation are described in [usage.md](usage.md).

## Install

From the repository root, using Python 3.10 or later:

```bash
python -m venv .venv
```

Install with the virtual environment's Python:

```bash
# macOS / Linux
.venv/bin/python -m pip install -r requirements-mcp.txt

# Windows PowerShell
.venv/Scripts/python.exe -m pip install -r requirements-mcp.txt
```

Use that same interpreter in the client configuration. The optional dependency
file includes the analysis packages and the MCP SDK; DSH and direct CLI usage
do not require the MCP SDK.

## Claude Desktop and other local clients

Launch the server with stdio:

```bash
python python/mcp_server.py --transport stdio --workspace /absolute/path/to/data
```

For clients using `mcpServers` configuration, add an entry with absolute paths:

```json
{
  "mcpServers": {
    "bulk-scrna-concordance": {
      "command": "/absolute/path/to/repository/.venv/bin/python",
      "args": [
        "/absolute/path/to/repository/python/mcp_server.py",
        "--transport", "stdio",
        "--workspace", "/absolute/path/to/data"
      ]
    }
  }
}
```

On Windows, use paths such as `C:/tools/repository/.venv/Scripts/python.exe`.
Merge this entry with any existing configuration. Reload the client and check
that `bulk_scrna_concordance` appears in its tool list. See the
[MCP host configuration guide](https://py.sdk.modelcontextprotocol.io/get-started/real-host/).

For Claude Code:

```bash
claude mcp add --transport stdio bulk-scrna-concordance -- /absolute/path/to/repository/.venv/bin/python /absolute/path/to/repository/python/mcp_server.py --workspace /absolute/path/to/data
```

See [Claude Code's MCP reference](https://code.claude.com/docs/en/mcp).
Local Claude Desktop MCP configuration does not automatically provide a
connector for claude.ai or Cowork; those use a separate remote connection.
[Claude's connector guide](https://support.claude.com/en/articles/11175166-get-started-with-custom-connectors-using-remote-mcp)
describes that distinction.

## Streamable HTTP and ChatGPT

Start a private HTTP endpoint:

```bash
python python/mcp_server.py --transport streamable-http --workspace /absolute/path/to/data --host 127.0.0.1 --port 8000
```

The endpoint is `http://127.0.0.1:8000/mcp`. HTTP mode requires an existing
workspace directory. The server does not implement authentication, TLS or
per-user isolation. Keep it on loopback for local use. A remote deployment
must supply authentication and HTTPS through its hosting/proxy layer and
restrict the data available to that deployment.

ChatGPT can connect to a private server through Secure MCP Tunnel, or to a
public HTTPS endpoint that forwards to `/mcp`. In ChatGPT's custom MCP server
setup, choose the matching connection method and configure authentication.
Availability depends on the account and workspace policies. Follow the
[official OpenAI connection guide](https://developers.openai.com/plugins/deploy/connect-chatgpt)
for current setup steps. This repository supplies the server, not the tunnel,
hosting or authentication service.

## Call the tool

With the repository itself configured as the workspace, this request uses the
included synthetic data:

```json
{
  "sc_matrix_csv": "examples/synthetic/sc_expr.csv",
  "sc_meta_csv": "examples/synthetic/sc_meta.csv",
  "bulk_csv": "examples/synthetic/bulk.csv",
  "group_col": "Group",
  "out_dir": "results",
  "out_prefix": "example",
  "discriminability": true,
  "no_figure": true
}
```

For AnnData, replace the CSV inputs with `sc_h5ad` and, when appropriate,
`sc_layer: "counts"`. The default prefers `layers['counts']`, then `X`;
`sc_use_raw: true` explicitly selects `raw.X` and `raw.var`.

All paths refer to the server's filesystem. Relative paths resolve under the
configured workspace; input files and output directories must stay within it,
including after symlink resolution. `--workspace` is optional for stdio; when
omitted, absolute paths are unrestricted and relative paths use the server's
working directory. Use a workspace when the client should be limited to a
particular data directory.

Each call creates a unique subdirectory inside `out_dir`. The result contains
the analysis JSON plus `results_json`, `run_dir` and `log_path`. Failed runs
return an MCP tool error with diagnostics. The default timeout is 30 minutes;
`timeout_ms` can set a positive value up to 24 hours, subject to the client's
own timeout.

Chat attachments are not transferred to the server automatically. This server
does not provide upload, download or image-resource endpoints; returned file
paths require access to the server filesystem to open the artifacts.

## Check the connection

List the server's tools and call the synthetic example with an MCP client or
[MCP Inspector](https://modelcontextprotocol.io/docs/tools/inspector). For local
HTTP, use Streamable HTTP and `http://127.0.0.1:8000/mcp`.

Run the repository checks from the project root with the MCP environment:

```bash
python -m unittest discover -s tests -v
```

