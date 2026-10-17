"""Numerical and file-level regressions for concordance correctness."""

import importlib.util
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import linregress, pearsonr


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "python" / "bulk_scrna_concordance.py"
SPEC = importlib.util.spec_from_file_location("concordance", SCRIPT)
c = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(c)


class GeneAggregationTests(unittest.TestCase):
    def test_normalization_keeps_one_id_per_input_column(self):
        ids, duplicates = c.normalize_gene_ids(
            ["gene.1", "GENE.2", "other.1"], strip_version=True, to_upper=True
        )
        self.assertEqual(ids, ["GENE", "GENE", "OTHER"])
        self.assertEqual(duplicates, 1)

    def test_dense_and_sparse_duplicate_columns_are_summed(self):
        counts = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        labels = pd.Series(["A", "A", "B"])
        for matrix in (counts, sparse.csr_matrix(counts)):
            with self.subTest(sparse=sparse.issparse(matrix)):
                pseudo, n_cells = c.compute_pseudo_bulk(matrix, ["G", "G", "H"], labels)
                self.assertEqual(pseudo.index.tolist(), ["G", "H"])
                np.testing.assert_array_equal(pseudo.to_numpy(), [[12, 15], [9, 9]])
                self.assertEqual(n_cells, {"A": 2, "B": 1})

    def test_invalid_counts_are_rejected_before_duplicate_sums(self):
        # Negative and fractional entries could cancel into valid-looking sums.
        for bad in (-1.0, 0.5, float("nan"), float("inf")):
            for as_sparse in (False, True):
                with self.subTest(value=bad, sparse=as_sparse):
                    matrix = np.array([[bad, 2.0], [3.0, 4.0]])
                    if as_sparse:
                        matrix = sparse.csr_matrix(matrix)
                    with self.assertRaises(ValueError):
                        c.compute_pseudo_bulk(matrix, ["G", "G"], pd.Series(["A", "A"]))


class StatisticsAndFigureTests(unittest.TestCase):
    def setUp(self):
        genes = list("ABCDEFG")
        self.pseudo = pd.Series([1, 4, 16, 64, 0, 9, 0], index=genes, dtype=float)
        self.bulk = pd.Series([2, 4, 8, 32, 100, 0, 0], index=genes, dtype=float)

    def test_regression_predicts_pseudobulk_from_bulk(self):
        stats = c.correlate_pair(self.pseudo, self.bulk, normalize="none")
        x = np.log2(np.array([2, 4, 8, 32]) + 1)
        y = np.log2(np.array([1, 4, 16, 64]) + 1)
        expected = linregress(x, y)
        self.assertEqual(stats["n_genes_used"], 4)
        self.assertAlmostEqual(stats["slope"], expected.slope)
        self.assertAlmostEqual(stats["intercept"], expected.intercept)

    def test_figure_and_stats_share_normalization_population_and_fit(self):
        # Singleton zeros count toward each library, but do not become plotted
        # points. The intentionally different totals expose early filtering.
        x = np.log2(np.array([2, 4, 8, 32]) / 146.0 * 1e6 + 1)
        y = np.log2(np.array([1, 4, 16, 64]) / 94.0 * 1e6 + 1)
        expected = linregress(x, y)
        stats = c.correlate_pair(self.pseudo, self.bulk)
        self.assertAlmostEqual(stats["pearson_r"], pearsonr(x, y).statistic)
        with tempfile.TemporaryDirectory(dir=ROOT) as directory:
            paths = [str(Path(directory) / ("figure." + ext)) for ext in ("png", "pdf", "svg")]
            # Preserve the real figure until its plotted numerical data have
            # been inspected; all rendering and file writes still run.
            with patch.object(c.plt, "close") as close_spy:
                saved = c.make_figure(
                    self.pseudo.to_frame("A"), {"A": self.bulk}, ["A"], {"A": stats}, *paths
                )
            fig = close_spy.call_args.args[0]
            try:
                self.assertEqual(set(saved), set(paths))
                ax = fig.axes[0]
                np.testing.assert_allclose(ax.images[0].get_extent(), [x.min(), x.max(), y.min(), y.max()])
                fitted_line = ax.lines[1]
                np.testing.assert_allclose(
                    fitted_line.get_ydata(),
                    expected.slope * fitted_line.get_xdata() + expected.intercept,
                )
            finally:
                c.plt.close(fig)

    def test_tpm_rejects_nonpositive_or_nonfinite_lengths(self):
        for bad in (0.0, -2.0, float("nan"), float("inf")):
            with self.subTest(length=bad), self.assertRaises(ValueError):
                c.normalize_values(np.array([1.0, 2.0, 3.0]), "tpm", np.array([1.0, bad, 3.0]))

    def test_duplicate_genes_require_consistent_lengths(self):
        with self.assertRaises(ValueError):
            c.collapse_gene_lengths(pd.Series([1000, 2000, 3000]), ["G", "G", "H"])
        lengths = c.collapse_gene_lengths(pd.Series([1000, 1000, 3000]), ["G", "G", "H"])
        self.assertEqual(lengths.to_dict(), {"G": 1000, "H": 3000})


class H5adSourceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=ROOT)
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "input.h5ad"
        obs = pd.DataFrame({"Group": ["A", "A"]}, index=["c1", "c2"])
        # Raw deliberately has another gene set, so taking raw.var with X (or
        # var with raw.X) cannot silently satisfy the expected dimensions.
        raw = ad.AnnData(
            X=np.array([[10, 20, 30], [40, 50, 60]]), obs=obs.copy(),
            var=pd.DataFrame({"symbol": ["R1", "R2", "R3"], "length": [900, 1800, 2700]},
                             index=["r1", "r2", "r3"]),
        )
        data = ad.AnnData(
            X=np.array([[9, 8], [7, 6]]), obs=obs,
            var=pd.DataFrame({"symbol": ["G1", "G2"], "length": [1000, 2000]}, index=["g1", "g2"]),
        )
        data.layers["counts"] = sparse.csr_matrix([[1, 2], [3, 4]])
        data.layers["custom"] = np.array([[5, 6], [7, 8]])
        data.raw = raw
        data.write_h5ad(self.path)

    def test_default_prefers_counts_layer_over_x_and_raw(self):
        loaded = c.load_scrna_from_h5ad(str(self.path), "Group", "symbol", "length")
        self.assertEqual(len(loaded), 6)
        matrix, ids, groups, id_source, lengths, matrix_source = loaded
        np.testing.assert_array_equal(matrix.toarray(), [[1, 2], [3, 4]])
        self.assertEqual(ids, ["G1", "G2"])
        self.assertEqual(groups.tolist(), ["A", "A"])
        self.assertEqual(id_source, "symbol")
        np.testing.assert_array_equal(lengths.to_numpy(), [1000, 2000])
        self.assertIn("counts", matrix_source)

    def test_default_without_counts_uses_x_even_when_raw_exists(self):
        data = ad.read_h5ad(self.path)
        del data.layers["counts"]
        data.write_h5ad(self.path)
        loaded = c.load_scrna_from_h5ad(str(self.path), "Group", "symbol")
        self.assertEqual(len(loaded), 6)
        np.testing.assert_array_equal(loaded[0], [[9, 8], [7, 6]])
        self.assertEqual(loaded[-1], "X")

    def test_explicit_layer_and_raw_keep_matching_gene_metadata(self):
        custom = c.load_scrna_from_h5ad(str(self.path), "Group", "symbol", "length", layer="custom")
        np.testing.assert_array_equal(custom[0], [[5, 6], [7, 8]])
        self.assertEqual(custom[1], ["G1", "G2"])
        self.assertIn("custom", custom[-1])
        raw = c.load_scrna_from_h5ad(str(self.path), "Group", "symbol", "length", use_raw=True)
        np.testing.assert_array_equal(raw[0], [[10, 20, 30], [40, 50, 60]])
        self.assertEqual(raw[1], ["R1", "R2", "R3"])
        np.testing.assert_array_equal(raw[4].to_numpy(), [900, 1800, 2700])
        self.assertIn("raw", raw[-1])

    def test_missing_or_conflicting_explicit_sources_fail(self):
        with self.assertRaises((KeyError, ValueError)):
            c.load_scrna_from_h5ad(str(self.path), "Group", "symbol", layer="absent")
        with self.assertRaises(ValueError):
            c.load_scrna_from_h5ad(str(self.path), "Group", "symbol", layer="custom", use_raw=True)

    def test_loader_rejects_invalid_selected_counts(self):
        for bad in (-1.0, 0.5, float("nan"), float("inf")):
            with self.subTest(value=bad):
                data = ad.read_h5ad(self.path)
                data.layers["counts"] = sparse.csr_matrix([[bad, 2.0], [3.0, 4.0]])
                data.write_h5ad(self.path)
                with self.assertRaises(ValueError):
                    c.load_scrna_from_h5ad(str(self.path), "Group", "symbol")


class CliRegressionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=ROOT)
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)

    def run_cli(self, *extra):
        return subprocess.run(
            [sys.executable, str(SCRIPT), "--out-dir", str(self.directory / "out"),
             "--results-json", str(self.directory / "result.json"), "--no-figure", *map(str, extra)],
            capture_output=True, text=True, encoding="utf-8", timeout=60,
            env={**os.environ, "PYTHONUTF8": "1", "PYTHONIOENCODING": "utf-8"},
        )

    def write_bulk(self, ids, values):
        path = self.directory / "bulk.csv"
        pd.DataFrame({"gene": ids, "A_1": values}).to_csv(path, index=False)
        return path

    def test_cli_sums_duplicate_genes_on_both_sides(self):
        expr = self.directory / "sc.csv"
        meta = self.directory / "meta.csv"
        ids = ["g.1", "g.2", "h.1", "i.1", "j.1"]
        pd.DataFrame({"c1": [1, 2, 2, 4, 8], "c2": [2, 4, 4, 8, 16]}, index=ids).to_csv(expr)
        pd.DataFrame({"cell": ["c1", "c2"], "Group": ["A", "A"]}).to_csv(meta, index=False)
        bulk = self.write_bulk(ids, [3, 3, 4, 8, 16])
        outcome = self.run_cli("--sc-matrix-csv", expr, "--sc-meta-csv", meta,
                               "--bulk-csv", bulk, "--strip-version", "--normalize", "none")
        self.assertEqual(outcome.returncode, 0, outcome.stdout + outcome.stderr)
        result = json.loads((self.directory / "result.json").read_text(encoding="utf-8"))
        pair = result["pairs"]["A"]
        expected = linregress(np.log2(np.array([6, 4, 8, 16]) + 1),
                             np.log2(np.array([9, 6, 12, 24]) + 1))
        self.assertEqual(result["scrna"]["n_genes_unique"], 4)
        self.assertEqual(pair["n_genes_used"], 4)
        self.assertAlmostEqual(pair["slope"], expected.slope)

    def test_cli_explicit_h5ad_layer_is_used(self):
        path = self.directory / "sc.h5ad"
        data = ad.AnnData(
            X=np.array([[1, 1, 1, 1], [1, 1, 1, 1]]),
            obs=pd.DataFrame({"Group": ["A", "A"]}, index=["c1", "c2"]),
            var=pd.DataFrame(index=list("ghij")),
        )
        data.layers["custom"] = sparse.csr_matrix([[1, 2, 4, 8], [2, 4, 8, 16]])
        data.write_h5ad(path)
        bulk = self.write_bulk(list("ghij"), [2, 4, 8, 16])
        outcome = self.run_cli("--sc-h5ad", path, "--sc-layer", "custom", "--bulk-csv", bulk,
                               "--normalize", "none")
        self.assertEqual(outcome.returncode, 0, outcome.stdout + outcome.stderr)
        result = json.loads((self.directory / "result.json").read_text(encoding="utf-8"))
        expected = linregress(np.log2(np.array([2, 4, 8, 16]) + 1),
                             np.log2(np.array([3, 6, 12, 24]) + 1))
        self.assertAlmostEqual(result["pairs"]["A"]["slope"], expected.slope)

    def test_csv_and_h5ad_reject_fractional_expression(self):
        ids = list("ghij")
        bulk = self.write_bulk(ids, [2, 4, 8, 16])
        expr = self.directory / "sc.csv"
        meta = self.directory / "meta.csv"
        counts = np.array([[0.5, 2, 4, 8], [2, 4, 8, 16]])
        pd.DataFrame(counts.T, index=ids, columns=["c1", "c2"]).to_csv(expr)
        obs = pd.DataFrame({"Group": ["A", "A"]}, index=["c1", "c2"])
        obs.rename_axis("cell").to_csv(meta)
        path = self.directory / "sc.h5ad"
        ad.AnnData(X=counts, obs=obs, var=pd.DataFrame(index=ids)).write_h5ad(path)
        for source in (("--sc-matrix-csv", expr, "--sc-meta-csv", meta), ("--sc-h5ad", path)):
            with self.subTest(source=source[0]):
                outcome = self.run_cli(*source, "--bulk-csv", bulk)
                self.assertNotEqual(outcome.returncode, 0, outcome.stdout + outcome.stderr)

    def test_csv_loader_rejects_invalid_counts(self):
        expr = self.directory / "sc.csv"
        meta = self.directory / "meta.csv"
        pd.DataFrame({"cell": ["c1", "c2"], "Group": ["A", "A"]}).to_csv(meta, index=False)
        for bad in (-1.0, 0.5, float("nan"), float("inf")):
            with self.subTest(value=bad):
                pd.DataFrame({"c1": [bad, 2], "c2": [3, 4]}, index=["G", "H"]).to_csv(expr)
                with self.assertRaises(ValueError):
                    c.load_scrna_from_csv(str(expr), str(meta), "Group")

    def test_result_links_manifest_and_manifest_hashes_final_result(self):
        outcome = self.run_cli(
            "--sc-matrix-csv", ROOT / "examples/synthetic/sc_expr.csv",
            "--sc-meta-csv", ROOT / "examples/synthetic/sc_meta.csv",
            "--bulk-csv", ROOT / "examples/synthetic/bulk.csv",
        )
        self.assertEqual(outcome.returncode, 0, outcome.stdout + outcome.stderr)
        result_path = self.directory / "result.json"
        result_bytes = result_path.read_bytes()
        result = json.loads(result_bytes)
        self.assertIsNotNone(result["outputs"]["manifest_json"])
        manifest_path = Path(result["outputs"]["manifest_json"])
        self.assertTrue(manifest_path.is_file())
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(manifest["outputs"][result_path.name]["sha256"],
                         hashlib.sha256(result_bytes).hexdigest())

    def test_single_group_discriminability_has_no_mismatched_estimate(self):
        expr = self.directory / "sc.csv"
        meta = self.directory / "meta.csv"
        pd.DataFrame({"c1": [1, 2, 4, 8], "c2": [2, 4, 8, 16]}, index=list("ghij")).to_csv(expr)
        pd.DataFrame({"cell": ["c1", "c2"], "Group": ["A", "A"]}).to_csv(meta, index=False)
        bulk = self.write_bulk(list("ghij"), [2, 4, 8, 16])
        outcome = self.run_cli("--sc-matrix-csv", expr, "--sc-meta-csv", meta,
                               "--bulk-csv", bulk, "--discriminability")
        self.assertEqual(outcome.returncode, 0, outcome.stdout + outcome.stderr)
        result = json.loads((self.directory / "result.json").read_text(encoding="utf-8"))
        diagnostic = result["discriminability"]
        self.assertEqual(diagnostic["n_groups"], 1)
        self.assertEqual(diagnostic["n_correct"], 1)
        self.assertAlmostEqual(diagnostic["matched_mean_r"], 1.0)
        self.assertIsNone(diagnostic["mismatched_mean_r"])
        self.assertIsNone(diagnostic["gap"])

    def test_cli_rejects_output_prefix_traversal_before_writing(self):
        outside = self.directory / "escape_bulk_scrna_correlation_stats.csv"
        outside.write_text("keep this unrelated table", encoding="utf-8")
        outcome = self.run_cli(
            "--sc-matrix-csv", ROOT / "examples/synthetic/sc_expr.csv",
            "--sc-meta-csv", ROOT / "examples/synthetic/sc_meta.csv",
            "--bulk-csv", ROOT / "examples/synthetic/bulk.csv",
            "--out-prefix", "../escape",
        )
        self.assertEqual(outside.read_text(encoding="utf-8"), "keep this unrelated table")
        self.assertNotEqual(outcome.returncode, 0, outcome.stdout + outcome.stderr)
        self.assertFalse((self.directory / "out").exists())


if __name__ == "__main__":
    unittest.main()

