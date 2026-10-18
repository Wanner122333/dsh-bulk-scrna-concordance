#!/usr/bin/env python3
"""Download and prepare the public GSE77288 real-data example.

Run from the repository root after installing requirements.txt. Source files
are checked against the SHA-256 hashes used for the published example; no
expression normalization or gene ranking is performed here.
"""

import argparse
import datetime
import gc
import gzip
import hashlib
import json
from pathlib import Path
import urllib.request


SOURCES = (
    {
        "filename": "GSE77288_reads-raw-single-per-sample.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77288/suppl/GSE77288_reads-raw-single-per-sample.txt.gz",
        "bytes": 13025232,
        "sha256": "90e6c5f93bc0e33ffd4fba84f9d5f6f2f7f234eef3a241ef629915121c4cb0c8",
    },
    {
        "filename": "GSE77288_reads-raw-bulk-per-sample.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77288/suppl/GSE77288_reads-raw-bulk-per-sample.txt.gz",
        "bytes": 330833,
        "sha256": "6e2a1a0bc59dd3159a43963bd916cbe68da921147f44c855fda1dbb2206dc787",
    },
    {
        "filename": "annotation-filter.txt",
        "url": "https://raw.githubusercontent.com/jdblischak/singleCellSeq/4c917b1e1087d292a4cddfb7a6a4548f1ca16431/data/annotation-filter.txt",
        "bytes": 23166,
        "sha256": "0a50191f9adf7be7d80f6ebbc472c774fae473b5aa471a59d09153d205b8d409",
    },
    {
        "filename": "annotation-bulk.txt",
        "url": "https://raw.githubusercontent.com/jdblischak/singleCellSeq/4c917b1e1087d292a4cddfb7a6a4548f1ca16431/data/annotation-bulk.txt",
        "bytes": 429,
        "sha256": "40e9048a1dc7b92c874f433dba69026e74216f6e5169f387ea78897835c744f4",
    },
)
DONORS = ["NA19098", "NA19101", "NA19239"]
SC_COUNTS = "GSE77288_reads-raw-single-per-sample.txt.gz"
BULK_COUNTS = "GSE77288_reads-raw-bulk-per-sample.txt.gz"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify_source(path, source):
    require(path.is_file(), f"Source file is missing: {path}")
    require(path.stat().st_size == source["bytes"], f"Unexpected byte count: {path}")
    require(sha256(path) == source["sha256"], f"SHA-256 mismatch: {path}")
    if source["filename"].endswith(".gz"):
        with gzip.open(path, "rb") as stream:
            while stream.read(1024 * 1024):
                pass  # Reaching EOF verifies gzip completeness and CRC.


def get_sources(raw_dir, reuse_only):
    if not reuse_only:
        raw_dir.mkdir(parents=True, exist_ok=True)
    for source in SOURCES:
        destination = raw_dir / source["filename"]
        if not destination.exists() and not reuse_only:
            partial = destination.with_name(destination.name + ".partial")
            request = urllib.request.Request(
                source["url"], headers={"User-Agent": "GSE77288-concordance-example/1.0"}
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                with partial.open("wb") as output:
                    for chunk in iter(lambda: response.read(1024 * 1024), b""):
                        output.write(chunk)
            verify_source(partial, source)
            partial.replace(destination)
        verify_source(destination, source)
        print(f"Verified {source['filename']}", flush=True)


def load_raw(path):
    import numpy as np
    import pandas as pd

    with gzip.open(path, "rt", encoding="utf-8") as stream:
        header = stream.readline().rstrip("\r\n").split("\t")
    require(header[:3] == ["individual", "replicate", "well"], "Unexpected count header")
    require(len(header) == len(set(header)), "Duplicate source columns")
    genes = [column for column in header[3:] if column.startswith("ENSG")]
    spikeins = [column for column in header[3:] if not column.startswith("ENSG")]
    require(len(genes) == 20327, "Expected 20,327 Ensembl genes")
    require(len(spikeins) == 92 and all(x.startswith("ERCC") for x in spikeins),
            "Unexpected spike-in columns")
    counts = pd.read_csv(
        path, sep="\t", usecols=header[:3] + genes,
        dtype={**dict.fromkeys(genes, np.int32), **dict.fromkeys(header[:3], str)},
    )
    sample_ids = counts["individual"] + "." + counts["replicate"] + "." + counts["well"]
    require(sample_ids.is_unique, "Duplicate source sample IDs")
    counts.index = pd.Index(sample_ids, name="sample_id")
    metadata = counts[header[:3]].copy()
    metadata["batch"] = metadata["individual"] + "." + metadata["replicate"]
    counts = counts[genes]
    require((counts.to_numpy() >= 0).all(), "Raw counts must be nonnegative")
    return counts, metadata, len(spikeins)


def prepare(raw_dir, out_dir):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse

    qc = pd.read_csv(raw_dir / "annotation-filter.txt", sep="\t", dtype=str).set_index("sample_id")
    bulk_metadata = pd.read_csv(raw_dir / "annotation-bulk.txt", sep="\t", dtype=str).set_index("sample_id")
    require(qc.index.is_unique and len(qc) == 564, "Expected 564 author-QC single cells")
    require(bulk_metadata.index.is_unique and len(bulk_metadata) == 9, "Expected 9 raw bulk libraries")
    require(sorted(qc.individual.unique()) == DONORS, "Unexpected donor labels")
    require(qc.batch.nunique() == 8 and "NA19098.r2" not in set(qc.batch),
            "Unexpected author-QC batches")

    sc, sc_metadata, sc_spikeins = load_raw(raw_dir / SC_COUNTS)
    require(len(sc) == 864, "Expected 864 raw single cells")
    require(set(qc.index) <= set(sc.index), "QC metadata has unknown cells")
    pd.testing.assert_frame_equal(sc_metadata.loc[qc.index, qc.columns], qc, check_dtype=False)
    sc = sc.loc[qc.index].copy()
    del sc_metadata
    gc.collect()

    bulk, raw_bulk_metadata, bulk_spikeins = load_raw(raw_dir / BULK_COUNTS)
    pd.testing.assert_frame_equal(
        raw_bulk_metadata.loc[bulk_metadata.index, bulk_metadata.columns], bulk_metadata,
        check_dtype=False,
    )
    require(list(sc.columns) == list(bulk.columns), "Single-cell and bulk genes differ")
    keep = bulk_metadata.batch.isin(set(qc.batch))
    excluded_bulk = bulk_metadata.index[~keep].tolist()
    require(excluded_bulk == ["NA19098.r2.bulk"], "Unexpected excluded bulk libraries")
    bulk_metadata = bulk_metadata.loc[keep].copy()
    require(len(bulk_metadata) == 8 and set(bulk_metadata.batch) == set(qc.batch),
            "The bulk and author-QC single-cell batches do not match")
    bulk = bulk.loc[bulk_metadata.index].T
    bulk.index.name = "gene_id"
    bulk.to_csv(out_dir / "bulk_read_counts.csv", lineterminator="\n")
    bulk_metadata.to_csv(out_dir / "bulk_metadata.csv", lineterminator="\n")
    qc.to_csv(out_dir / "sc_metadata.csv", lineterminator="\n")

    obs = qc.copy()
    obs["individual"] = pd.Categorical(obs.individual, categories=DONORS)
    adata = ad.AnnData(
        X=sparse.csr_matrix(sc.to_numpy(dtype=np.int32)), obs=obs,
        var=pd.DataFrame(index=pd.Index(sc.columns, name="gene_id")),
    )
    adata.uns["source"] = "GSE77288 raw read counts; author annotation-filter.txt QC"
    adata.uns["expression_units"] = "untransformed integer read counts, ERCC excluded"
    adata.write_h5ad(out_dir / "sc_read_counts.h5ad", compression="gzip")
    summary = {
        "study": "Tung et al. 2017 / GSE77288; human iPSC",
        "raw_single_cells": 864,
        "author_qc_single_cells": len(qc),
        "excluded_single_cells": 864 - len(qc),
        "raw_bulk_libraries": 9,
        "matched_bulk_libraries": len(bulk_metadata),
        "excluded_bulk_libraries": excluded_bulk,
        "genes_before_expression_filter": adata.n_vars,
        "excluded_ERCC_columns": {"single_cell": sc_spikeins, "bulk": bulk_spikeins},
        "single_cells_by_donor": qc.groupby("individual").size().to_dict(),
        "bulk_libraries_by_donor": bulk_metadata.groupby("individual").size().to_dict(),
        "single_cells_by_batch": qc.groupby("batch").size().to_dict(),
        "sc_matrix_source": "AnnData.X raw integer read counts (cells x genes)",
        "bulk_matrix_source": "independently measured bulk raw read counts (genes x samples)",
        "preparation": "Author fixed QC list; exclude ERCC; align 8 measured batches; no CPM, log, imputation, or gene selection before tool execution.",
        "analysis": "Donor-level sc sum and bulk raw-count mean, then default CPM over pairwise shared genes positive on either side; retain both-positive and log2(CPM+1).",
        "interpretation": "Three donors of the same iPSC cell type. Gene-wise correlations are descriptive, not donor-level replicate tests or cell-type discrimination.",
        "inputs": {},
    }
    for name in ("bulk_metadata.csv", "bulk_read_counts.csv", "sc_metadata.csv", "sc_read_counts.h5ad"):
        path = out_dir / name
        summary["inputs"][name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    (out_dir / "preparation_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("results/gse77288/input-data"),
                        help="Prepared inputs and provenance directory (default: %(default)s)")
    parser.add_argument("--source-dir", type=Path,
                        help="Reuse these four existing source files without downloading; all hashes are still checked")
    args = parser.parse_args()
    out_dir = args.out_dir.resolve()
    raw_dir = args.source_dir.resolve() if args.source_dir else out_dir / "raw"
    get_sources(raw_dir, reuse_only=args.source_dir is not None)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = prepare(raw_dir, out_dir)
    provenance = {
        "checked_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "study": "Tung et al. 2017 / GSE77288; human iPSC",
        "sources": list(SOURCES),
    }
    (out_dir / "download_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output_directory": str(out_dir), **summary}, indent=2), flush=True)


if __name__ == "__main__":
    main()
