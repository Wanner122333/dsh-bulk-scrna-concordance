# GSE77288 real-data example

Tung et al. measured single-cell and independent bulk RNA-seq in human iPSCs
from three donors. The example uses the authors' fixed list of 564 QC-passing
cells and the 8 bulk libraries whose batches remain after that QC. The ninth
bulk library, `NA19098.r2.bulk`, is omitted because that single-cell batch has
no retained cells.

Sources: [GEO GSE77288](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77288)
and the [authors' count preparation workflow](https://jdblischak.github.io/singleCellSeq/analysis/prepare-counts.html).

## Run from the repository root

```bash
python -m pip install -r requirements.txt
python examples/gse77288/prepare.py --out-dir results/gse77288/input-data
python python/bulk_scrna_concordance.py --sc-h5ad results/gse77288/input-data/sc_read_counts.h5ad --bulk-csv results/gse77288/input-data/bulk_read_counts.csv --group-col individual --gene-id-key var_names --discriminability --out-dir results/gse77288 --out-prefix GSE77288_real --results-json results/gse77288/GSE77288_real_results.json
```

These commands also work in PowerShell. The preparation step downloads four
public source files (about 13 MB total), verifies their fixed SHA-256 hashes
and gzip integrity, excludes the 92 ERCC spike-ins, and writes 20,327 Ensembl
genes as untransformed integer read counts. It does not normalize, impute or
rank genes. Preparation can use several hundred MB of RAM.

To reuse a downloaded source directory without network access:

```bash
python examples/gse77288/prepare.py --source-dir /path/to/raw --out-dir results/gse77288/input-data
```

The source directory must contain the two `GSE77288_reads-raw-*-per-sample.txt.gz`
files and `annotation-filter.txt` / `annotation-bulk.txt`. Existing files are
still checked against the pinned hashes. Missing or changed sources cause an
error rather than silently changing the example.

The two annotation URLs are pinned to author-repository commit
`4c917b1e1087d292a4cddfb7a6a4548f1ca16431`. Count files come from the GEO
supplement directory. The fixed hashes identify the exact downloaded content;
they are not an independently supplied GEO checksum certification.

## Prepared inputs

| File | Contents |
|---|---|
| `sc_read_counts.h5ad` | 564 cells × 20,327 genes; raw counts in `X`, donor in `obs['individual']` |
| `bulk_read_counts.csv` | 20,327 genes × 8 independent bulk libraries |
| `sc_metadata.csv`, `bulk_metadata.csv` | Retained sample IDs, donors and batches |
| `download_provenance.json` | Public URLs, sizes and fixed source hashes |
| `preparation_summary.json` | Selection rules, sample counts and prepared-input hashes |

The analysis sums single-cell counts by donor and averages bulk counts by
donor, using raw-count means for bulk before normalization. It applies the
tool's default CPM normalization and `log2(x+1)` transform to each comparison.
The CPM denominator includes shared genes positive on either side, before
filtering to genes positive on both sides. Only those both-positive genes enter the plotted
correlations and regression; the number of plotted genes therefore differs
between donors.

These are gene-wise descriptive comparisons among three donors of the same
cell type. High correlations do not establish donor-level statistical
significance or cell-type discrimination. Numerical output can be compared
with the [published result table](results/GSE77288_real_bulk_scrna_correlation_stats.csv).
Rendering details may vary with package versions.

## Reference run

The [reference summary](reference_run.json) records analysis commit
`7478bdef898cdca66ba05b101f6a5cbfe239035e`, full-precision results, all source
hashes, parameters, environment versions and figure/CSV hashes. Its paths are
relative to the example directory; it is a portable summary rather than an
unchanged runtime manifest. The five files in `results/` are unchanged copies
of that run's PNG, PDF, SVG, statistics and cross-donor correlation matrix.

The reference environment was Python 3.13.7, NumPy 2.2.6, pandas 2.3.2,
SciPy 1.16.1, Matplotlib 3.10.6, AnnData 0.12.2 and h5py 3.14.0. The repository
requirements specify supported minimum versions; they are not an exact lock
of that environment. New runs write their own JSON and manifest with hashes
and installed versions.

Raw counts were independently compared entry by entry with the prepared
inputs, and correlations/regression coefficients were independently recomputed
within a tolerance of 1e-12. The very small p-values were not independently
recomputed: their numeric zero is floating-point underflow, displayed as
`p < 1e-300` rather than interpreted as an exact zero.
