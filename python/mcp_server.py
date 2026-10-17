#!/usr/bin/env python
"""Expose the existing analysis through official MCP stdio/Streamable HTTP.

Paths refer to files on the server machine. HTTP mode requires an operator
workspace, and every requested input/output path must stay within that root.
Each call gets a unique output directory, including calls using the same prefix.
"""

import argparse
import asyncio
from collections import deque
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import Annotated, Any, Literal

from mcp.server.fastmcp import FastMCP
from mcp.server.fastmcp.exceptions import ToolError
from pydantic import Field, StrictBool, StrictInt


SCRIPT = Path(__file__).with_name("bulk_scrna_concordance.py")
SCHEMA = "dsh.bulk-scrna-concordance/1"


def resolve_path(value: str, workspace: Path | None, *, input_file: bool) -> Path:
    if not value or not value.strip():
        raise ToolError("Input and output paths must be non-empty strings.")
    candidate = Path(value).expanduser()
    if not candidate.is_absolute():
        candidate = (workspace or Path.cwd()) / candidate
    candidate = candidate.resolve()
    if workspace is not None and not candidate.is_relative_to(workspace):
        raise ToolError("Path is outside the configured workspace: %s" % value)
    if input_file and not candidate.is_file():
        raise ToolError("Input file does not exist: %s" % value)
    return candidate


def log_tail(log_path: Path) -> str:
    with log_path.open(encoding="utf-8", errors="replace") as handle:
        return "".join(deque(handle, maxlen=30)).strip()


def run_analysis(workspace: Path | None, options: dict[str, Any]) -> dict[str, Any]:
    has_h5ad = bool(options["sc_h5ad"])
    has_csv = bool(options["sc_matrix_csv"])
    if has_h5ad == has_csv:
        raise ToolError("Provide exactly one of sc_h5ad or sc_matrix_csv.")
    if has_csv and not options["sc_meta_csv"]:
        raise ToolError("sc_matrix_csv requires sc_meta_csv.")
    if options["sc_layer"] is not None and options["sc_use_raw"]:
        raise ToolError("Choose at most one of sc_layer and sc_use_raw.")
    if has_csv and (options["sc_layer"] is not None or options["sc_use_raw"]):
        raise ToolError("sc_layer and sc_use_raw require sc_h5ad.")
    if options["sc_layer"] is not None and not options["sc_layer"].strip():
        raise ToolError("sc_layer must be a non-empty layer name.")
    prefix = options["out_prefix"]
    if not prefix or prefix in (".", "..") or any(c in prefix for c in "/\\:\x00"):
        raise ToolError("out_prefix must be a filename prefix without path separators.")

    inputs = {}
    for key in ("bulk_csv", "sc_h5ad", "sc_matrix_csv", "sc_meta_csv"):
        if options[key]:
            inputs[key] = resolve_path(options[key], workspace, input_file=True)
    out_dir = resolve_path(options["out_dir"], workspace, input_file=False)
    out_dir.mkdir(parents=True, exist_ok=True)
    run_dir = Path(tempfile.mkdtemp(prefix=prefix + "-", dir=out_dir))
    results_path = run_dir / (prefix + "_results.json")
    log_path = run_dir / (prefix + "_run.log")
    argv = [sys.executable, "-X", "utf8", str(SCRIPT), "--out-dir", str(run_dir),
            "--out-prefix", prefix, "--results-json", str(results_path)]
    for key, value in inputs.items():
        argv.extend(["--" + key.replace("_", "-"), str(value)])
    for key in ("group_col", "gene_id_key", "bulk_gene_col", "sc_layer",
                "normalize", "sc_lengths_col", "bulk_lengths_col"):
        if options[key] is not None:
            argv.extend(["--" + key.replace("_", "-"), options[key]])
    for key in ("strip_version", "sc_use_raw", "discriminability",
                "keep_missing_groups", "no_manifest", "no_figure"):
        if options[key]:
            argv.append("--" + key.replace("_", "-"))
    if options["gene_id_case"] != "keep":
        argv.append("--gene-to-" + options["gene_id_case"])

    # stdout belongs exclusively to the MCP stdio protocol; analysis output is
    # never inherited or piped into it. Keep diagnostics with the run artifacts.
    try:
        with log_path.open("w", encoding="utf-8") as log_handle:
            completed = subprocess.run(
                argv, stdin=subprocess.DEVNULL, stdout=log_handle,
                stderr=subprocess.STDOUT, cwd=workspace or Path.cwd(),
                timeout=options["timeout_ms"] / 1000, check=False,
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
            )
    except subprocess.TimeoutExpired as exc:
        raise ToolError("Analysis timed out after %d ms. Log: %s\n%s"
                        % (options["timeout_ms"], log_path, log_tail(log_path))) from exc
    except OSError as exc:
        raise ToolError("Could not start analysis: %s" % exc) from exc
    if completed.returncode != 0:
        raise ToolError("Analysis exited with code %d. Log: %s\n%s"
                        % (completed.returncode, log_path, log_tail(log_path)))
    try:
        result = json.loads(results_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ToolError("Analysis wrote no readable result JSON. Log: %s" % log_path) from exc
    if not isinstance(result, dict) or result.get("schema") != SCHEMA or result.get("status") != "ok":
        raise ToolError("Analysis returned an unexpected result schema/status.")
    return {**result, "results_json": str(results_path), "run_dir": str(run_dir),
            "log_path": str(log_path)}


def create_server(workspace: Path | None = None, host: str = "127.0.0.1", port: int = 8000) -> FastMCP:
    if workspace is not None:
        workspace = workspace.expanduser().resolve(strict=True)
        if not workspace.is_dir():
            raise ValueError("workspace must be an existing directory")
    server = FastMCP("bulk-scrna-concordance", host=host, port=port,
                     json_response=True, stateless_http=True)

    @server.tool()
    async def bulk_scrna_concordance(
        bulk_csv: str,
        out_dir: str,
        sc_h5ad: str | None = None,
        sc_matrix_csv: str | None = None,
        sc_meta_csv: str | None = None,
        out_prefix: str = "bulk_scrna_concordance",
        group_col: str = "Group",
        gene_id_key: str = "auto",
        bulk_gene_col: str | None = None,
        sc_layer: str | None = None,
        sc_use_raw: StrictBool = False,
        strip_version: StrictBool = False,
        gene_id_case: Literal["keep", "upper", "lower"] = "keep",
        discriminability: StrictBool = False,
        keep_missing_groups: StrictBool = False,
        normalize: Literal["cpm", "tpm", "none"] = "cpm",
        sc_lengths_col: str | None = None,
        bulk_lengths_col: str | None = None,
        no_manifest: StrictBool = False,
        no_figure: StrictBool = False,
        timeout_ms: Annotated[StrictInt, Field(gt=0, le=86_400_000)] = 1_800_000,
    ) -> dict[str, Any]:
        """Compare bulk and group-summed single-cell counts descriptively.

        All input/output paths are server-local. Relative paths resolve within
        the operator workspace (or working directory for unrestricted stdio).
        Select raw counts with sc_layer or explicit sc_use_raw; by default the
        counts layer is preferred, falling back to .X. Each call creates a unique
        subdirectory inside out_dir and returns the analysis JSON plus
        results_json, run_dir and log_path. Correlations
        describe pooled group profiles and do not test sample-level agreement.
        """
        options = locals().copy()
        return await asyncio.to_thread(run_analysis, workspace, options)

    return server


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--transport", choices=("stdio", "streamable-http"), default="stdio")
    parser.add_argument("--workspace", type=Path,
                        help="Operator-approved local root; required for HTTP, optional for stdio.")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args(argv)
    if args.transport == "streamable-http" and args.workspace is None:
        parser.error("--workspace is required for streamable-http")
    if not 1 <= args.port <= 65535:
        parser.error("--port must be between 1 and 65535")
    try:
        server = create_server(args.workspace, args.host, args.port)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    server.run(transport=args.transport)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

