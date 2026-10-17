"""Real MCP protocol integration tests (install requirements-mcp.txt first)."""

import asyncio
from contextlib import asynccontextmanager
from pathlib import Path
import shutil
import socket
import subprocess
import sys
import tempfile
import unittest

import httpx
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from mcp.client.streamable_http import streamable_http_client


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "python" / "mcp_server.py"


class ServerAvailabilityTests(unittest.TestCase):
    def test_server_module_exists(self):
        self.assertTrue(SCRIPT.is_file(), "the MCP server entry point must exist")

    def test_http_requires_operator_workspace(self):
        run = subprocess.run(
            [sys.executable, str(SCRIPT), "--transport", "streamable-http"],
            capture_output=True, text=True, timeout=15,
        )
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("--workspace is required", run.stderr)


class ProtocolTests(unittest.IsolatedAsyncioTestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="concordance-mcp-", dir=ROOT)
        self.addCleanup(self.temp.cleanup)
        self.workspace = Path(self.temp.name).resolve()
        examples = ROOT / "examples" / "synthetic"
        for filename in ("sc_expr.csv", "sc_meta.csv", "bulk.csv"):
            shutil.copyfile(examples / filename, self.workspace / filename)
        self.arguments = {
            "sc_matrix_csv": "sc_expr.csv", "sc_meta_csv": "sc_meta.csv",
            "bulk_csv": "bulk.csv", "out_dir": "results", "out_prefix": "same",
            "no_figure": True, "discriminability": True,
        }

    @asynccontextmanager
    async def stdio_session(self, script=SCRIPT):
        parameters = StdioServerParameters(
            command=sys.executable,
            args=[str(script), "--transport", "stdio", "--workspace", str(self.workspace)],
        )
        with (self.workspace / "stdio-server.log").open("w", encoding="utf-8") as error_log:
            async with stdio_client(parameters, errlog=error_log) as (read, write):
                async with ClientSession(read, write) as session:
                    await session.initialize()
                    yield session

    def assert_success(self, response):
        self.assertFalse(response.isError, response.content)
        data = response.structuredContent
        self.assertEqual(data["status"], "ok")
        self.assertTrue(data["pairs"])
        self.assertTrue(Path(data["results_json"]).is_file())
        self.assertTrue(Path(data["log_path"]).is_file())
        self.assertTrue(Path(data["run_dir"]).is_relative_to(self.workspace))
        # Real analysis logs are retained as artifacts, not protocol stdout.
        self.assertIn("[OK] Analysis complete.", Path(data["log_path"]).read_text(encoding="utf-8"))
        return data

    async def test_stdio_initialize_list_call_and_unique_runs(self):
        async with self.stdio_session() as session:
            listing = await session.list_tools()
            self.assertEqual([tool.name for tool in listing.tools], ["bulk_scrna_concordance"])
            schema = listing.tools[0].inputSchema
            self.assertEqual(schema["properties"]["sc_use_raw"]["type"], "boolean")
            self.assertIn("sc_layer", schema["properties"])
            responses = await asyncio.gather(*(
                session.call_tool("bulk_scrna_concordance", self.arguments) for _ in range(2)))
            first, second = map(self.assert_success, responses)
            self.assertNotEqual(first["run_dir"], second["run_dir"])

    async def test_stdio_errors_are_tool_errors_and_server_recovers(self):
        async with self.stdio_session() as session:
            cases = [
                ({"bulk_csv": "missing.csv"}, "does not exist"),
                ({"sc_layer": "counts"}, "require sc_h5ad"),
                ({"sc_h5ad": "absent.h5ad", "sc_matrix_csv": None,
                  "sc_layer": "counts", "sc_use_raw": True}, "at most one"),
                ({"out_prefix": "../escape"}, "filename prefix"),
                ({"out_dir": "../escape"}, "outside the configured workspace"),
                ({"bulk_csv": str(ROOT / "examples" / "synthetic" / "bulk.csv")}, "outside the configured workspace"),
                ({"timeout_ms": 1}, "timed out"),
                ({"sc_use_raw": "false"}, "valid boolean"),
            ]
            for overrides, expected in cases:
                with self.subTest(overrides=overrides):
                    response = await session.call_tool(
                        "bulk_scrna_concordance", {**self.arguments, **overrides})
                    self.assertTrue(response.isError)
                    self.assertIn(expected, str(response.content))
            self.assert_success(await session.call_tool("bulk_scrna_concordance", self.arguments))

    async def test_nonzero_analysis_is_tool_error(self):
        # A real malformed count matrix reaches the real core and fails there.
        source = (self.workspace / "sc_expr.csv").read_text()
        lines = source.splitlines()
        cells = lines[1].split(",")
        cells[1] = "-1"
        lines[1] = ",".join(cells)
        (self.workspace / "sc_expr.csv").write_text("\n".join(lines) + "\n")
        async with self.stdio_session() as session:
            response = await session.call_tool("bulk_scrna_concordance", self.arguments)
            self.assertTrue(response.isError)
            self.assertIn("exited with code", str(response.content))

    async def test_successful_child_without_json_is_tool_error(self):
        # Exercise the real SDK and child process with a deliberately incomplete
        # interpreter fixture, rather than mocking the MCP session or subprocess.
        fixture_dir = self.workspace / "incomplete-package"
        fixture_dir.mkdir()
        entrypoint = fixture_dir / "mcp_server.py"
        shutil.copyfile(SCRIPT, entrypoint)
        (fixture_dir / "bulk_scrna_concordance.py").write_text("raise SystemExit(0)\n")
        async with self.stdio_session(entrypoint) as session:
            response = await session.call_tool("bulk_scrna_concordance", self.arguments)
            self.assertTrue(response.isError)
            self.assertIn("no readable result JSON", str(response.content))

    async def test_http_initialize_list_call_and_workspace_boundary(self):
        with socket.socket() as reservation:
            reservation.bind(("127.0.0.1", 0))
            port = reservation.getsockname()[1]
        log_path = self.workspace / "http-server.log"
        with log_path.open("w", encoding="utf-8") as error_log:
            process = subprocess.Popen(
                [sys.executable, str(SCRIPT), "--transport", "streamable-http",
                 "--workspace", str(self.workspace), "--port", str(port)],
                stdin=subprocess.DEVNULL, stdout=error_log, stderr=subprocess.STDOUT,
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
            )
            try:
                for _ in range(100):
                    if process.poll() is not None:
                        self.fail("HTTP server exited: " + log_path.read_text())
                    try:
                        _, writer = await asyncio.wait_for(
                            asyncio.open_connection("127.0.0.1", port), timeout=0.5)
                    except (OSError, asyncio.TimeoutError):
                        await asyncio.sleep(0.1)
                    else:
                        writer.close()
                        await writer.wait_closed()
                        break
                else:
                    self.fail("HTTP server did not become ready")
                async with httpx.AsyncClient(timeout=httpx.Timeout(10, read=30)) as client:
                    async with streamable_http_client(
                        "http://127.0.0.1:%d/mcp" % port, http_client=client,
                    ) as (read, write, _):
                        async with ClientSession(read, write) as session:
                            await asyncio.wait_for(session.initialize(), timeout=15)
                            tools = await session.list_tools()
                            self.assertEqual(tools.tools[0].name, "bulk_scrna_concordance")
                            self.assert_success(await session.call_tool("bulk_scrna_concordance", self.arguments))
                            response = await session.call_tool("bulk_scrna_concordance", {
                                **self.arguments, "bulk_csv": "../outside.csv",
                            })
                            self.assertTrue(response.isError)
                            self.assertIn("outside the configured workspace", str(response.content))
            finally:
                process.terminate()
                try:
                    await asyncio.to_thread(process.wait, timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    await asyncio.to_thread(process.wait, timeout=5)


if __name__ == "__main__":
    unittest.main()

