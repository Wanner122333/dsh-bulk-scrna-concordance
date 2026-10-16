#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Bulk vs scRNA pseudo-bulk concordance analysis (DSH Harness edition).

Derived from the author's TBtools plugin source `bulk_scrna_correlation_cli.py`
(which is byte-identical to the EXE inside the original .plugin package), with
the defects found during code review fixed and a discriminability assessment
added.

Fixed relative to the original
-----------------------------
1. `--strip_version` was `action="store_true", default=True`, so it could never
   be turned off. Now `--strip-version / --no-strip-version` (default: off,
   matching the original argparse author intent that the flag be optional).
2. Gene IDs collapsing to duplicates after version stripping crashed the
   original with
       ValueError: operands could not be broadcast together with shapes (n,) (m,)
   because `.loc[common_genes, g]` returned a DataFrame for a duplicated label.
   Duplicate IDs are now collapsed explicitly (values summed) and reported.
3. `infer_bulk_group()` matched sample names with loose `startswith`/`endswith`
   /substring tests, so groups that are prefixes of one another (L1 vs L10) were
   silently mis-assigned. Matching now uses exact underscore/dot/space-delimited
   tokens, and unmatched samples are reported instead of silently dropped.
4. The `nan` group (cells whose group label is missing) was treated as a real
   group and given its own pseudo-bulk column. It is now excluded and reported.

Added relative to the original
------------------------------
5. Discriminability matrix: every pseudo-bulk group against every bulk group,
   with matched-vs-mismatched summary. A single r value cannot show that a
   correlation is meaningful; this can.
6. `--results-json`: one machine-readable result object, so a caller does not
   have to scrape stdout.
7. Run provenance recorded in the JSON (parameters, input paths, versions).
8. CPM normalization is computed over the genes actually used in the
   correlation (the shared, informative set) instead of over the whole shared
   gene list as in the original. Pearson r and Spearman rho are unchanged by
   this, but the fitted regression `slope`/`intercept` are on a different scale
   than the original outputs; a linear rescaling of either axis changes the
   slope while leaving r identical.

Note on stdout/stdio
--------------------
This script prints progress to stdout. The DSH sandbox denies piped stdio to
child processes, so the caller runs it with inherited stdio and reads the JSON
result file instead of parsing stdout.
"""

import argparse
import hashlib
import importlib
import json
import os
import platform
import sys
import textwrap
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

from scipy import sparse
from scipy.stats import pearsonr, spearmanr, linregress

try:
    import anndata as ad
except ImportError:  # pragma: no cover
    ad = None

SCHEMA = "dsh.bulk-scrna-concordance/1"
# Group labels that mean "missing" rather than a real experimental group.
MISSING_GROUP_LABELS = {"nan", "na", "n/a", "none", "null", "", "unknown", "unk"}


# --------------------------------------------------------------------------
# Utilities
# --------------------------------------------------------------------------

def log(msg):
    print(msg, flush=True)


def format_p(p):
    """Format very small p-values nicely."""
    try:
        if p is None:
            return "NA"
        if p < 1e-300:
            return "<1e-300"
        return "%.1e" % p
    except (TypeError, ValueError):
        return "NA"


def normalize_gene_ids(values, strip_version=False, to_upper=False, to_lower=False):
    """Normalize gene IDs and collapse duplicates.

    Returns (ids, n_collapsed) where ids is a list of unique strings, in first
    appearance order, and n_collapsed counts how many input IDs were merged
    into an existing one.
    """
    out = []
    seen = {}
    n_collapsed = 0
    for raw in values:
        g = str(raw)
        if strip_version and "." in g:
            g = g.split(".", 1)[0]
        if to_upper:
            g = g.upper()
        elif to_lower:
            g = g.lower()
        if g in seen:
            n_collapsed += 1
            continue
        seen[g] = len(out)
        out.append(g)
    return out, n_collapsed


def collapse_matrix_by_gene(X, gene_ids):
    """Sum rows (genes) that share a normalized ID.

    X is cells x genes. Returns (X_collapsed, unique_ids).
    """
    order = []
    index_of = {}
    for g in gene_ids:
        if g not in index_of:
            index_of[g] = len(order)
            order.append(g)
    if len(order) == len(gene_ids):
        return X, gene_ids

    n_cells = X.shape[0]
    n_unique = len(order)
    mapping = np.array([index_of[g] for g in gene_ids], dtype=np.int64)

    if sparse.issparse(X):
        # group columns: X @ P where P is genes x unique, P[i, mapping[i]] = 1
        P = sparse.csr_matrix(
            (np.ones(len(gene_ids)), (np.arange(len(gene_ids)), mapping)),
            shape=(len(gene_ids), n_unique),
        )
        return (X @ P).tocsr(), order

    out = np.zeros((n_cells, n_unique), dtype=np.float64)
    for j, m in enumerate(mapping):
        out[:, m] += np.asarray(X[:, j], dtype=np.float64).ravel()
    return out, order


def tokenize_sample_name(name):
    """Split a sample name into normalized tokens.

    Splits on the separators authors actually use in sample names, so that
    matching is token-based rather than substring-based.
    """
    s = str(name)
    for ch in ("\\", "/", "|", ",", ";", ":", "(", ")", "[", "]", "+"):
        s = s.replace(ch, "_")
    s = s.replace("-", "_").replace(".", "_").replace(" ", "_")
    return [t for t in s.split("_") if t != ""]


def infer_bulk_group(sample_name, groups):
    """Return the group a bulk sample belongs to, using exact token matching.

    A sample matches group `g` when `g` appears as a whole token, or when the
    name equals `g`. This replaces the original loose substring tests, which
    mapped 'L10_a' onto group 'L1'.
    """
    tokens = tokenize_sample_name(sample_name)
    token_set = set(tokens)
    for g in groups:
        g_str = str(g)
        if g_str in token_set:
            return g_str
        # multi-token group labels, e.g. group "L1_treat" against "Rtan_L1_treat_1"
        g_tokens = tokenize_sample_name(g_str)
        if g_tokens and len(g_tokens) > 1:
            n = len(g_tokens)
            for i in range(len(tokens) - n + 1):
                if tokens[i:i + n] == g_tokens:
                    return g_str
    return None


def is_missing_label(value):
    return str(value).strip().lower() in MISSING_GROUP_LABELS


# --------------------------------------------------------------------------
# Provenance
# --------------------------------------------------------------------------

#: Nothing in this analysis is stochastic. The seed is recorded anyway so that
#: claim is explicit and checkable rather than assumed.
SEED = 42


def sha256_file(path, chunk=1 << 20):
    """Hash a file without loading it whole (inputs can be hundreds of MB)."""
    h = hashlib.sha256()
    try:
        with open(path, "rb") as fh:
            while True:
                block = fh.read(chunk)
                if not block:
                    break
                h.update(block)
    except OSError:
        return None
    return h.hexdigest()


def file_provenance(path):
    """Size, mtime and content hash for one input file."""
    try:
        st = os.stat(path)
    except OSError:
        return {"path": path, "exists": False}
    return {
        "path": path,
        "exists": True,
        "bytes": int(st.st_size),
        "mtime_utc": datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat(),
        "sha256": sha256_file(path),
    }


def environment_lock():
    """The versions that decide the numbers.

    Recording only `python`/`numpy`/`pandas` is not enough: anndata decides how
    `.X` and `.raw` are read and scipy decides every statistic computed here, so
    a silent major bump in either can move a reported r.
    """
    lock = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "seed": SEED,
    }
    for name in ("numpy", "pandas", "scipy", "anndata", "h5py", "matplotlib"):
        try:
            lock[name] = importlib.import_module(name).__version__
        except Exception:
            lock[name] = None
    return lock


# --------------------------------------------------------------------------
# Input loading
# --------------------------------------------------------------------------

def load_scrna_from_h5ad(path, group_col, gene_id_key, lengths_col=None):
    if ad is None:
        raise RuntimeError(
            "anndata is required to read .h5ad files but is not importable by "
            "this interpreter. Install anndata+h5py, or pass --sc-matrix-csv and "
            "--sc-meta-csv instead."
        )
    log("[INFO] Reading h5ad: %s" % path)
    adata = ad.read_h5ad(path)

    if adata.raw is not None:
        log("[INFO] Using adata.raw.X / adata.raw.var")
        X = adata.raw.X
        var = adata.raw.var.copy()
    else:
        log("[INFO] adata.raw is None; using adata.X / adata.var")
        X = adata.X
        var = adata.var.copy()

    obs = adata.obs.copy()
    if group_col not in obs.columns:
        raise KeyError(
            "group_col=%r not found in adata.obs. Available: %s"
            % (group_col, list(obs.columns))
        )

    var_cols = list(var.columns)
    if gene_id_key in (None, "auto"):
        chosen, used = None, None
        for ck in ("gene_symbols", "gene_symbol", "symbol", "symbols", "GeneSymbol",
                   "Gene_Symbol", "gene_name", "gene", "Gene", "GENE"):
            if ck in var_cols:
                chosen, used = var[ck].astype(str).tolist(), ck
                break
        if chosen is None:
            log("[WARN] No typical gene-ID column in var; falling back to var.index")
            chosen, used = var.index.astype(str).tolist(), "var_names"
    elif gene_id_key == "var_names":
        chosen, used = var.index.astype(str).tolist(), "var_names"
    elif gene_id_key in var_cols:
        chosen, used = var[gene_id_key].astype(str).tolist(), gene_id_key
    else:
        raise ValueError(
            "gene_id_key=%r not in var.columns (%s) and not 'var_names'/'auto'."
            % (gene_id_key, var_cols)
        )
    log("[INFO] scRNA gene-ID source: %s" % used)

    n_cells, n_genes = X.shape
    if n_genes != len(chosen):
        raise ValueError(
            "h5ad matrix has %d genes but the selected gene-ID vector has %d; "
            "var/raw.var does not match X/raw.X." % (n_genes, len(chosen))
        )
    group_series = obs[group_col]
    if len(group_series) != n_cells:
        raise ValueError(
            "obs has %d rows but X has %d cells." % (len(group_series), n_cells)
        )

    lengths = None
    if lengths_col is not None:
        if lengths_col not in var_cols:
            raise KeyError(
                "sc_lengths_col=%r not in var columns %s." % (lengths_col, var_cols)
            )
        lengths = pd.Series(
            np.asarray(var[lengths_col], dtype=np.float64),
            index=pd.Index(chosen),
        )
        log("[INFO] scRNA gene lengths from var['%s']" % lengths_col)

    return X, chosen, group_series, used, lengths


def load_scrna_from_csv(expr_path, meta_path, group_col):
    sep = "\t" if expr_path.lower().endswith((".tsv", ".txt")) else ","
    log("[INFO] Reading scRNA expression matrix: %s (sep=%r)" % (expr_path, sep))
    expr_df = pd.read_csv(expr_path, sep=sep, index_col=0)
    log("[INFO] Expression matrix: genes=%d cells=%d" % expr_df.shape)

    msep = "\t" if meta_path.lower().endswith((".tsv", ".txt")) else ","
    log("[INFO] Reading scRNA metadata: %s (sep=%r)" % (meta_path, msep))
    meta_df = pd.read_csv(meta_path, sep=msep)
    if meta_df.shape[1] < 2:
        raise ValueError(
            "Metadata file needs at least 2 columns; got %d. If the file is "
            "tab-separated, name it .tsv/.txt or make it comma-separated."
            % meta_df.shape[1]
        )
    cell_col = meta_df.columns[0]
    meta_df = meta_df.set_index(cell_col)
    if group_col not in meta_df.columns:
        raise KeyError(
            "group_col=%r not in metadata columns %s." % (group_col, list(meta_df.columns))
        )

    common = expr_df.columns.intersection(meta_df.index)
    if len(common) == 0:
        raise ValueError("No overlapping cells between matrix columns and metadata index.")
    if len(common) < expr_df.shape[1]:
        log("[WARN] %d / %d cells have metadata; using the overlapping subset."
            % (len(common), expr_df.shape[1]))

    # Guard against duplicated cell IDs, which would break positional alignment.
    if meta_df.index.duplicated().any():
        n_dup = int(meta_df.index.duplicated().sum())
        log("[WARN] Metadata has %d duplicated cell IDs; keeping the first occurrence." % n_dup)
        meta_df = meta_df[~meta_df.index.duplicated(keep="first")]
        common = expr_df.columns.intersection(meta_df.index)

    expr_df = expr_df[common]
    meta_df = meta_df.loc[common]
    return expr_df.T.values, expr_df.index.astype(str).tolist(), meta_df[group_col], "csv"


def load_bulk(path, gene_col):
    sep = "\t" if path.lower().endswith((".tsv", ".txt")) else ","
    log("[INFO] Reading bulk matrix: %s (sep=%r)" % (path, sep))
    df = pd.read_csv(path, sep=sep)
    if gene_col is None:
        df = df.set_index(df.columns[0])
    else:
        if gene_col not in df.columns:
            raise KeyError(
                "bulk_gene_col=%r not in bulk columns %s." % (gene_col, list(df.columns))
            )
        df = df.set_index(gene_col)
    return df


# --------------------------------------------------------------------------
# Core computation
# --------------------------------------------------------------------------

def compute_pseudo_bulk(X, gene_ids, group_series, exclude_missing=True):
    """Sum counts per group, after collapsing duplicate gene IDs."""
    labels = group_series.astype(str)
    groups = sorted(set(labels.tolist()))

    excluded = []
    if exclude_missing:
        keep_groups = [g for g in groups if not is_missing_label(g)]
        excluded = [g for g in groups if is_missing_label(g)]
    else:
        keep_groups = groups
    if not keep_groups:
        raise ValueError("No usable groups after excluding missing labels %s." % excluded)

    X, unique_ids = collapse_matrix_by_gene(X, gene_ids)
    n_collapsed = len(gene_ids) - len(unique_ids)
    if n_collapsed:
        log("[INFO] Collapsed %d duplicate gene IDs into %d unique genes (values summed)."
            % (n_collapsed, len(unique_ids)))

    log("[INFO] Groups used: %s" % keep_groups)
    if excluded:
        n_excl = int(labels.isin(excluded).sum())
        log("[WARN] Excluded group label(s) %s covering %d cells (treated as missing)."
            % (excluded, n_excl))

    n_genes = X.shape[1]
    pseudo = np.zeros((n_genes, len(keep_groups)), dtype=np.float64)
    cell_counts = {}
    for j, g in enumerate(keep_groups):
        mask = (labels == g).values
        n_g = int(mask.sum())
        cell_counts[g] = n_g
        if n_g == 0:
            log("[WARN] group=%s has 0 cells." % g)
            continue
        sub = X[mask, :]
        log("[INFO] group=%s cells=%d" % (g, n_g))
        if sparse.issparse(sub):
            pseudo[:, j] = np.asarray(sub.sum(axis=0)).ravel()
        else:
            pseudo[:, j] = np.asarray(sub, dtype=np.float64).sum(axis=0)

    return pd.DataFrame(pseudo, index=pd.Index(unique_ids), columns=keep_groups), cell_counts


def aggregate_bulk_by_group(bulk_df, groups, unmatched_report=True):
    """Mean bulk expression per group using token-exact sample matching."""
    assigned = {}
    unmatched = []
    for col in bulk_df.columns:
        g = infer_bulk_group(col, groups)
        if g is None:
            unmatched.append(str(col))
        else:
            assigned.setdefault(g, []).append(col)

    if unmatched_report and unmatched:
        log("[WARN] %d bulk sample(s) matched no group and were NOT used: %s"
            % (len(unmatched), unmatched))
    missing_groups = [g for g in groups if g not in assigned]
    if missing_groups:
        log("[WARN] No bulk samples matched group(s): %s" % missing_groups)

    out = {}
    for g in groups:
        cols = assigned.get(g)
        if not cols:
            continue
        log("[INFO] group=%s <- bulk samples %s" % (g, cols))
        out[g] = bulk_df[cols].mean(axis=1)
    if not out:
        raise ValueError(
            "No groups had matching bulk samples. Bulk columns: %s ; scRNA groups: %s"
            % (list(bulk_df.columns), groups)
        )
    return out, {g: assigned.get(g, []) for g in groups}, unmatched, missing_groups


def normalize_values(values, method, gene_lengths=None):
    """Scale one expression vector to a fixed library size.

    `cpm` is the historical behaviour and the default. `tpm` divides by gene
    length first, which changes what the numbers mean (and therefore the fitted
    slope) but leaves Pearson/Spearman rank agreement essentially intact.
    `none` passes values through, for inputs already on a common scale.
    """
    if method == "none":
        return values
    if method == "cpm":
        total = values.sum()
        return values / total * 1e6 if total > 0 else values
    if method == "tpm":
        if gene_lengths is None:
            raise ValueError("tpm normalization requires gene lengths")
        rate = values / gene_lengths
        total = rate.sum()
        return rate / total * 1e6 if total > 0 else rate
    raise ValueError("unknown normalization method %r" % (method,))


def correlate_pair(a, b, normalize="cpm", lengths_a=None, lengths_b=None):
    """Normalize, log2(x+1), then Pearson/Spearman/regression."""
    idx = a.index.intersection(b.index)
    x = np.asarray(a.loc[idx].values, dtype=np.float64).ravel()
    y = np.asarray(b.loc[idx].values, dtype=np.float64).ravel()
    la = None if lengths_a is None else np.asarray(lengths_a.loc[idx].values, dtype=np.float64).ravel()
    lb = None if lengths_b is None else np.asarray(lengths_b.loc[idx].values, dtype=np.float64).ravel()
    keep = (x > 0) | (y > 0)
    x, y = x[keep], y[keep]
    if la is not None:
        la = la[keep]
    if lb is not None:
        lb = lb[keep]
    if x.size == 0:
        return None
    if x.sum() <= 0 or y.sum() <= 0:
        return None
    x = normalize_values(x, normalize, la)
    y = normalize_values(y, normalize, lb)
    keep = (x > 0) & (y > 0)
    n = int(keep.sum())
    if n < 3:
        return None
    x, y = x[keep], y[keep]
    lx = np.log2(x + 1.0)
    ly = np.log2(y + 1.0)
    r, p_r = pearsonr(lx, ly)
    rho, p_s = spearmanr(lx, ly)
    slope, intercept, _, _, _ = linregress(lx, ly)
    return {
        "n_genes_used": n,
        "pearson_r": float(r),
        "pearson_p": float(p_r),
        "spearman_rho": float(rho),
        "spearman_p": float(p_s),
        "slope": float(slope),
        "intercept": float(intercept),
    }


def make_figure(pseudo_df, bulk_group_expr, groups, pairs, out_png, out_pdf, out_svg,
                normalize="cpm", lengths_a=None, lengths_b=None):
    n = len(groups)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5), dpi=150)
    if n == 1:
        axes = [axes]
    im_for_cbar = None
    for i, g in enumerate(groups):
        ax = axes[i]
        st = pairs[g]
        idx = pseudo_df.index.intersection(bulk_group_expr[g].index)
        x = np.asarray(bulk_group_expr[g].loc[idx].values, dtype=np.float64).ravel()
        y = np.asarray(pseudo_df.loc[idx, g].values, dtype=np.float64).ravel()
        la = None if lengths_b is None else np.asarray(lengths_b.loc[idx].values, dtype=np.float64).ravel()
        lb = None if lengths_a is None else np.asarray(lengths_a.loc[idx].values, dtype=np.float64).ravel()
        keep = (x > 0) & (y > 0)
        x, y = x[keep], y[keep]
        if la is not None:
            la = la[keep]
        if lb is not None:
            lb = lb[keep]
        x = normalize_values(x, normalize, la)
        y = normalize_values(y, normalize, lb)
        keep = (x > 0) & (y > 0)
        x, y = x[keep], y[keep]
        x = np.log2(x + 1.0)
        y = np.log2(y + 1.0)

        H, xe, ye = np.histogram2d(x, y, bins=200)
        H = H.T
        H[H < 1] = 0
        extent = [xe[0], xe[-1], ye[0], ye[-1]]
        vmax = max(10, H.max()) if H.max() > 0 else 10
        im = ax.imshow(H, origin="lower", extent=extent, cmap="viridis",
                       norm=LogNorm(vmin=1, vmax=vmax), aspect="auto")
        im_for_cbar = im
        lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
        ax.plot([lo, hi], [lo, hi], linestyle="--", color="gray", linewidth=1)
        xx = np.linspace(lo, hi, 100)
        ax.plot(xx, st["slope"] * xx + st["intercept"], color="blue", linewidth=1.5)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xlabel("log2(Bulk RNA-seq expression + 1)")
        ax.set_ylabel("log2(Pseudo-bulk scRNA-seq expression + 1)")
        ax.set_title("Group %s (n=%d genes)" % (g, st["n_genes_used"]))
        ax.text(0.02, 0.98,
                "Pearson r = %.3f (p=%s)\nSpearman \u03c1 = %.3f (p=%s)\ny = %.3fx + %.3f"
                % (st["pearson_r"], format_p(st["pearson_p"]),
                   st["spearman_rho"], format_p(st["spearman_p"]),
                   st["slope"], st["intercept"]),
                transform=ax.transAxes, va="top", ha="left", fontsize=9,
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))
    plt.subplots_adjust(right=0.9)
    if im_for_cbar is not None:
        cax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        cb = fig.colorbar(im_for_cbar, cax=cax)
        cb.set_label("Gene count")
    # Deliberately no tight_layout(): the manually placed colorbar axes is not
    # compatible with it, and calling it emitted a warning in the original.
    saved = []
    for path, fmt in ((out_png, "png"), (out_pdf, "pdf"), (out_svg, "svg")):
        try:
            fig.savefig(path, dpi=300)
            saved.append(path)
            log("[DONE] Figure: %s" % path)
        except Exception as e:
            log("[WARN] Failed to save %s: %s" % (fmt.upper(), e))
    plt.close(fig)
    return saved


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(
        description="Bulk vs scRNA pseudo-bulk concordance (DSH Harness edition)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
            Examples

              python bulk_scrna_concordance.py \\
                --sc-h5ad scRNA_with_celltype.h5ad \\
                --bulk-csv bulk_FPKM.csv \\
                --group-col Group \\
                --out-dir results --out-prefix my_corr

              # comma-separated input, keep version suffixes, assess discrimination
              python bulk_scrna_concordance.py \\
                --sc-matrix-csv sc_expr.csv --sc-meta-csv sc_meta.csv \\
                --bulk-csv bulk.csv --no-strip-version --discriminability
        """),
    )
    sc = p.add_argument_group("single-cell input")
    sc.add_argument("--sc-h5ad")
    sc.add_argument("--sc-matrix-csv")
    sc.add_argument("--sc-meta-csv")
    sc.add_argument("--group-col", default="Group")
    sc.add_argument("--gene-id-key", default="auto")

    b = p.add_argument_group("bulk input")
    b.add_argument("--bulk-csv", required=True)
    b.add_argument("--bulk-gene-col", default=None)

    g = p.add_argument_group("gene ID matching")
    g.add_argument("--strip-version", action=argparse.BooleanOptionalAction, default=False,
                   help="Strip the suffix after the first '.' from gene IDs on both sides.")
    g.add_argument("--gene-to-upper", action="store_true")
    g.add_argument("--gene-to-lower", action="store_true")

    o = p.add_argument_group("output")
    o.add_argument("--out-dir", required=True)
    o.add_argument("--out-prefix", default="bulk_scrna_concordance")
    o.add_argument("--results-json", default=None,
                   help="Where to write the machine-readable result object.")
    o.add_argument("--no-figure", action="store_true")

    a = p.add_argument_group("analysis")
    a.add_argument("--discriminability", action="store_true",
                   help="Correlate every pseudo-bulk group against every bulk group.")
    a.add_argument("--keep-missing-groups", action="store_true",
                   help="Do not exclude missing group labels such as 'nan'.")
    a.add_argument("--normalize", choices=("cpm", "tpm", "none"), default="cpm",
                   help="Library-size scaling before log2. 'cpm' is the historical "
                        "behaviour; 'tpm' also divides by gene length and needs "
                        "--sc-lengths-col/--bulk-lengths-col; 'none' assumes the "
                        "inputs are already on a common scale.")
    a.add_argument("--sc-lengths-col", default=None,
                   help="h5ad var column holding gene length, for --normalize tpm.")
    a.add_argument("--bulk-lengths-col", default=None,
                   help="bulk column holding gene length, for --normalize tpm. "
                        "The column is removed from the expression matrix.")
    a.add_argument("--no-manifest", action="store_true",
                   help="Skip manifest.json (provenance, hashes, environment lock).")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    if bool(args.sc_h5ad) == bool(args.sc_matrix_csv):
        raise SystemExit("Provide exactly one of --sc-h5ad or --sc-matrix-csv.")
    if args.sc_matrix_csv and not args.sc_meta_csv:
        raise SystemExit("--sc-matrix-csv requires --sc-meta-csv.")
    if args.gene_to_upper and args.gene_to_lower:
        raise SystemExit("Choose at most one of --gene-to-upper / --gene-to-lower.")

    os.makedirs(args.out_dir, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()

    # --- load single cell ---
    if args.sc_h5ad:
        X, raw_ids, group_series, id_source, sc_lengths = load_scrna_from_h5ad(
            args.sc_h5ad, args.group_col, args.gene_id_key, args.sc_lengths_col)
    else:
        X, raw_ids, group_series, id_source = load_scrna_from_csv(
            args.sc_matrix_csv, args.sc_meta_csv, args.group_col)
        sc_lengths = None
        if args.sc_lengths_col is not None:
            log("[WARN] --sc-lengths-col is only read from .h5ad var; ignored for the CSV path.")

    n_cells, n_genes = X.shape
    log("[INFO] scRNA matrix: cells=%d genes=%d" % (n_cells, n_genes))

    gene_ids, n_dup = normalize_gene_ids(
        raw_ids, strip_version=args.strip_version,
        to_upper=args.gene_to_upper, to_lower=args.gene_to_lower)
    if n_dup:
        log("[INFO] %d gene ID(s) became duplicates after normalization." % n_dup)

    pseudo_df, cell_counts = compute_pseudo_bulk(
        X, gene_ids, group_series, exclude_missing=not args.keep_missing_groups)
    sc_groups = list(pseudo_df.columns)

    # --- load bulk ---
    bulk_df = load_bulk(args.bulk_csv, args.bulk_gene_col)

    # Gene lengths for TPM. The column is kept out of the expression matrix so
    # it cannot be mistaken for a sample.
    bulk_lengths = None
    if args.bulk_lengths_col is not None:
        if args.bulk_lengths_col not in bulk_df.columns:
            raise KeyError(
                "bulk_lengths_col=%r not in bulk columns %s."
                % (args.bulk_lengths_col, list(bulk_df.columns))
            )
        bulk_lengths = pd.Series(
            np.asarray(bulk_df[args.bulk_lengths_col], dtype=np.float64),
            index=bulk_df.index,
        )
        bulk_df = bulk_df.drop(columns=[args.bulk_lengths_col])
        log("[INFO] bulk gene lengths from column '%s'" % args.bulk_lengths_col)

    log("[INFO] bulk matrix: genes=%d samples=%d" % bulk_df.shape)

    bulk_ids, bulk_dup = normalize_gene_ids(
        bulk_df.index.astype(str).tolist(), strip_version=args.strip_version,
        to_upper=args.gene_to_upper, to_lower=args.gene_to_lower)
    if bulk_dup:
        log("[WARN] %d duplicate bulk gene ID(s) after normalization; "
            "first occurrence kept." % bulk_dup)
    bulk_ids = pd.Index(bulk_ids)
    bulk_df = bulk_df.loc[~bulk_ids.duplicated()].copy()
    bulk_df.index = bulk_ids[~bulk_ids.duplicated()]
    if bulk_lengths is not None:
        bulk_lengths = bulk_lengths.loc[~bulk_ids.duplicated()]
        bulk_lengths.index = bulk_ids[~bulk_ids.duplicated()]

    if args.normalize == "tpm" and (sc_lengths is None or bulk_lengths is None):
        raise SystemExit(
            "--normalize tpm needs gene lengths on both sides: pass "
            "--sc-lengths-col and --bulk-lengths-col."
        )
    log("[INFO] normalization: %s" % args.normalize)

    bulk_group_expr, group_samples, unmatched, missing_groups = aggregate_bulk_by_group(
        bulk_df, sc_groups)

    # --- matched pairs ---
    groups = [g for g in sc_groups if g in bulk_group_expr]
    if not groups:
        raise SystemExit("No group has both pseudo-bulk and bulk data; cannot proceed.")

    def pair_stats(pseudo_series, bulk_series):
        return correlate_pair(
            pseudo_series, bulk_series, normalize=args.normalize,
            lengths_a=sc_lengths, lengths_b=bulk_lengths,
        )

    pairs = {}
    for g in groups:
        st = pair_stats(pseudo_df[g], bulk_group_expr[g])
        if st is None:
            log("[WARN] group=%s has too few informative genes; skipped." % g)
            continue
        n_overlap = len(pseudo_df.index.intersection(bulk_group_expr[g].index))
        st["n_overlap_genes"] = int(n_overlap)
        st["n_cells"] = int(cell_counts.get(g, 0))
        st["bulk_samples"] = list(group_samples.get(g, []))
        pairs[g] = st
        log("[INFO] group=%s Pearson r=%.3f (p=%s) Spearman rho=%.3f (p=%s) slope=%.3f"
            % (g, st["pearson_r"], format_p(st["pearson_p"]),
               st["spearman_rho"], format_p(st["spearman_p"]), st["slope"]))

    if not pairs:
        raise SystemExit("No group produced a usable correlation.")

    # --- discriminability ---
    discriminability = None
    if args.discriminability:
        log("\n[INFO] Computing discriminability matrix (pseudo x bulk) ...")
        matrix = {}
        for pg in sc_groups:
            row = {}
            for bg in bulk_group_expr:
                st = pair_stats(pseudo_df[pg], bulk_group_expr[bg])
                if st is not None:
                    row[bg] = {"pearson_r": st["pearson_r"], "spearman_rho": st["spearman_rho"]}
            matrix[pg] = row
        matched, mismatched, best = [], [], {}
        for pg, row in matrix.items():
            if pg not in row:
                continue
            m = row[pg]["pearson_r"]
            matched.append(m)
            offs = [v["pearson_r"] for bg, v in row.items() if bg != pg]
            mismatched.extend(offs)
            best_r = max(row.items(), key=lambda kv: kv[1]["pearson_r"])
            best[pg] = {
                "best_bulk_group": best_r[0],
                "best_r": best_r[1]["pearson_r"],
                "correct": best_r[0] == pg,
            }
        gap = (float(np.mean(matched)) - float(np.mean(mismatched))) \
            if matched and mismatched else None
        discriminability = {
            "matrix": matrix,
            "matched_mean_r": float(np.mean(matched)) if matched else None,
            "mismatched_mean_r": float(np.mean(mismatched)) if mismatched else None,
            "gap": gap,
            "best_match": best,
            "n_correct": sum(1 for v in best.values() if v["correct"]),
            "n_groups": len(best),
        }
        log("[INFO] matched mean r=%.4f  mismatched mean r=%.4f  gap=%.4f"
            % (discriminability["matched_mean_r"], discriminability["mismatched_mean_r"],
               gap if gap is not None else float("nan")))
        for pg, info in best.items():
            if not info["correct"]:
                log("[WARN] pseudo-bulk %s best-matches bulk %s (r=%.4f), not itself."
                    % (pg, info["best_bulk_group"], info["best_r"]))

    # --- outputs ---
    figure_paths = []
    if not args.no_figure:
        base = os.path.join(args.out_dir, args.out_prefix + "_bulk_scrna_correlation")
        figure_paths = make_figure(
            pseudo_df, bulk_group_expr, list(pairs.keys()), pairs,
            base + ".png", base + ".pdf", base + ".svg",
            normalize=args.normalize, lengths_a=sc_lengths, lengths_b=bulk_lengths)

    stats_df = pd.DataFrame([{
        "Group": g, "n_cells": s["n_cells"], "n_overlap_genes": s["n_overlap_genes"],
        "n_genes_used": s["n_genes_used"], "Pearson_r": s["pearson_r"],
        "Pearson_p": s["pearson_p"], "Spearman_rho": s["spearman_rho"],
        "Spearman_p": s["spearman_p"], "slope": s["slope"], "intercept": s["intercept"],
        "bulk_samples": ";".join(s["bulk_samples"]),
    } for g, s in pairs.items()])
    stats_path = os.path.join(args.out_dir, args.out_prefix + "_bulk_scrna_correlation_stats.csv")
    stats_df.to_csv(stats_path, index=False)
    log("[DONE] Stats table: %s" % stats_path)

    disc_path = None
    if discriminability is not None:
        mat = pd.DataFrame(
            {bg: {pg: discriminability["matrix"][pg][bg]["pearson_r"]
                  for pg in discriminability["matrix"] if bg in discriminability["matrix"][pg]}
             for bg in bulk_group_expr})
        disc_path = os.path.join(args.out_dir, args.out_prefix + "_discriminability.csv")
        mat.to_csv(disc_path)
        log("[DONE] Discriminability matrix: %s" % disc_path)

    result = {
        "schema": SCHEMA,
        "status": "ok",
        "generated_at_utc": started,
        "inputs": {
            "sc_h5ad": args.sc_h5ad,
            "sc_matrix_csv": args.sc_matrix_csv,
            "sc_meta_csv": args.sc_meta_csv,
            "bulk_csv": args.bulk_csv,
            "group_col": args.group_col,
            "gene_id_key": args.gene_id_key,
            "gene_id_source_used": id_source,
        },
        "parameters": {
            "strip_version": bool(args.strip_version),
            "gene_to_upper": bool(args.gene_to_upper),
            "gene_to_lower": bool(args.gene_to_lower),
            "discriminability": bool(args.discriminability),
            "keep_missing_groups": bool(args.keep_missing_groups),
            "normalize": args.normalize,
            "sc_lengths_col": args.sc_lengths_col,
            "bulk_lengths_col": args.bulk_lengths_col,
        },
        "scrna": {
            "n_cells": int(n_cells),
            "n_genes_raw": int(n_genes),
            "n_genes_unique": int(len(pseudo_df.index)),
            "n_gene_ids_collapsed": int(n_dup),
            "groups": sc_groups,
            "cells_per_group": {k: int(v) for k, v in cell_counts.items()},
        },
        "bulk": {
            "n_genes": int(bulk_df.shape[0]),
            "n_samples": int(bulk_df.shape[1]),
            "unmatched_samples": unmatched,
            "groups_without_bulk": missing_groups,
        },
        "pairs": pairs,
        "discriminability": discriminability,
        "warnings": [],
        "outputs": {
            "stats_csv": stats_path,
            "discriminability_csv": disc_path,
            "figures": figure_paths,
            "results_json": args.results_json,
            "manifest_json": None,
        },
        "environment": environment_lock(),
        "inputs_provenance": {
            key: file_provenance(value)
            for key, value in (
                ("sc_h5ad", args.sc_h5ad),
                ("sc_matrix_csv", args.sc_matrix_csv),
                ("sc_meta_csv", args.sc_meta_csv),
                ("bulk_csv", args.bulk_csv),
            )
            if value
        },
    }
    if unmatched:
        result["warnings"].append(
            "%d bulk sample(s) matched no scRNA group and were excluded: %s"
            % (len(unmatched), unmatched))
    if missing_groups:
        result["warnings"].append(
            "No bulk samples matched scRNA group(s): %s" % missing_groups)
    if discriminability and discriminability["n_correct"] < discriminability["n_groups"]:
        result["warnings"].append(
            "%d/%d pseudo-bulk groups did not best-match their own bulk group; the "
            "correlation is not reliably discriminative."
            % (discriminability["n_groups"] - discriminability["n_correct"],
               discriminability["n_groups"]))
    if n_dup:
        result["warnings"].append(
            "%d gene ID(s) collapsed onto an existing ID after normalization; their "
            "pseudo-bulk values were summed." % n_dup)

    # Write the results JSON before the manifest, so the manifest can hash it.
    if args.results_json:
        os.makedirs(os.path.dirname(os.path.abspath(args.results_json)), exist_ok=True)
        with open(args.results_json, "w", encoding="utf-8") as fh:
            json.dump(result, fh, ensure_ascii=False, indent=2)
        log("[DONE] Results JSON: %s" % args.results_json)

    # --- manifest: the one file a reader needs to reproduce this run ---
    if not args.no_manifest:
        manifest_path = os.path.join(
            args.out_dir, args.out_prefix + "_manifest.json")
        manifest = {
            "schema": "dsh.bulk-scrna-concordance.manifest/1",
            "produced_by": {
                "script": os.path.abspath(__file__),
                "script_sha256": sha256_file(os.path.abspath(__file__)),
                "seed": SEED,
            },
            "generated_at_utc": started,
            "finished_at_utc": datetime.now(timezone.utc).isoformat(),
            "command": {
                "argv": sys.argv[1:],
                "normalize": args.normalize,
                "group_col": args.group_col,
            },
            "inputs": result["inputs_provenance"],
            "environment_lock": result["environment"],
            "summary": {
                g: {
                    "pearson_r": s["pearson_r"],
                    "spearman_rho": s["spearman_rho"],
                    "slope": s["slope"],
                    "intercept": s["intercept"],
                    "n_genes_used": s["n_genes_used"],
                    "n_cells": s["n_cells"],
                }
                for g, s in pairs.items()
            },
            "discriminability": None if discriminability is None else {
                "matched_mean_r": discriminability["matched_mean_r"],
                "mismatched_mean_r": discriminability["mismatched_mean_r"],
                "gap": discriminability["gap"],
                "n_correct": discriminability["n_correct"],
                "n_groups": discriminability["n_groups"],
            },
            "outputs": {},
        }
        for path in [stats_path, disc_path, *figure_paths, args.results_json]:
            if not path:
                continue
            manifest["outputs"][os.path.basename(path)] = file_provenance(path)
        with open(manifest_path, "w", encoding="utf-8") as fh:
            json.dump(manifest, fh, ensure_ascii=False, indent=2)
        result["outputs"]["manifest_json"] = manifest_path
        log("[DONE] Manifest: %s" % manifest_path)

    log("[OK] Analysis complete.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:  # make failures machine-readable too
        log("[ERROR] %s: %s" % (type(exc).__name__, exc))
        sys.exit(1)
