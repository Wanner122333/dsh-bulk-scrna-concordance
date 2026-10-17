# Usage

The DSH and MCP tool is `bulk_scrna_concordance`. See [MCP setup](mcp.md) for
ChatGPT and Claude connections. The same analysis is available as a Python CLI:

```bash
python python/bulk_scrna_concordance.py --help
```

## Analysis

For each group, the tool:

1. Sums single-cell raw counts across cells and averages bulk expression across
   the samples assigned to that group.
2. Intersects gene IDs and excludes genes that are zero on both sides.
3. Applies the selected normalization to each vector, then retains genes
   positive on both sides and applies `log2(x+1)`.
4. Computes Pearson `r`, Spearman `rho` and an ordinary least-squares fit
   `pseudo = slope * bulk + intercept`.

The figure uses the same genes and transformed values as the statistics:
bulk on the x-axis and pseudo-bulk on the y-axis. A pair needs at least three
retained genes.

## Single-cell inputs

### AnnData (`.h5ad`)

```bash
python python/bulk_scrna_concordance.py \
  --sc-h5ad scRNA.h5ad --sc-layer counts \
  --bulk-csv bulk.csv --group-col Group \
  --out-dir results --out-prefix my_run
```

The default matrix is `adata.layers['counts']` when present, otherwise
`adata.X`. Use `--sc-layer <name>` for a named layer or `--sc-use-raw`
to request `adata.raw.X` explicitly. These options are mutually exclusive,
apply only to `.h5ad` input, and fail if the requested source is absent.
The default never selects `.raw` implicitly.

Gene IDs are read from the selected matrix's metadata: `adata.var` for
`X`/layers, or `adata.raw.var` for explicit `.raw` selection. Automatic
selection tries recognized columns such as `gene_symbols`, `symbol` and
`gene_name`, then falls back to the gene index. Override it with
`--gene-id-key`. Gene lengths use that same metadata.

### Expression CSV and cell metadata

```bash
python python/bulk_scrna_concordance.py \
  --sc-matrix-csv sc_counts.csv --sc-meta-csv sc_meta.csv \
  --bulk-csv bulk.csv --group-col Group \
  --out-dir results --out-prefix my_run
```

- Expression: raw counts, genes as rows, cells as columns; the first column is
  the gene ID.
- Metadata: the first column is the cell ID; `--group-col` identifies the
  group column.
- `.tsv` and `.txt` files are tab-separated; other extensions are read as CSV.

Both input routes require finite, nonnegative, integer-valued raw counts
before aggregation. Validation catches common normalized/log-transformed
matrices but cannot prove that integer-valued data are untransformed. Confirm
the source provenance; a `counts` layer name or `.raw` label is not sufficient.

Cells with missing group labels such as `nan`, `NA` and `none` are
excluded and reported by default.

## Bulk input

Genes are rows and samples are columns. The first column is the gene ID unless
`--bulk-gene-col` selects another column. Values may be nonnegative continuous
expression measures such as counts, CPM, TPM or FPKM. Choose normalization with
the units and prior processing in mind; applying CPM to FPKM does not turn it
into raw-count CPM.

Samples are assigned by exact delimited group tokens. With groups `A` and
`B`, both `A_1` and `sample_A_rep1` match `A`. Unmatched sample names
are reported in the log and `bulk.unmatched_samples`. Use unambiguous names
containing only the intended group token.

## Parameters

| Tool argument (DSH and MCP) | CLI flag | Default | Meaning |
|---|---|---|---|
| `sc_h5ad` | `--sc-h5ad` | — | AnnData input; exactly one of this or `sc_matrix_csv` |
| `sc_layer` | `--sc-layer` | — | Named raw-count layer; excludes `sc_use_raw` |
| `sc_use_raw` | `--sc-use-raw` | off | Explicitly select `raw.X` and `raw.var`; excludes `sc_layer` |
| `sc_matrix_csv` | `--sc-matrix-csv` | — | Raw-count genes × cells matrix |
| `sc_meta_csv` | `--sc-meta-csv` | — | Required with `sc_matrix_csv` |
| `bulk_csv` | `--bulk-csv` | — | Required bulk expression matrix |
| `bulk_gene_col` (MCP only) | `--bulk-gene-col` | first column | Bulk gene-ID column |
| `out_dir` | `--out-dir` | — | Required output directory; created if absent |
| `out_prefix` | `--out-prefix` | `bulk_scrna_concordance` | Output filename prefix |
| — | `--results-json` | off | Results-JSON path; DSH and MCP supply this automatically |
| `group_col` | `--group-col` | `Group` | Cell metadata group column |
| `gene_id_key` | `--gene-id-key` | `auto` | AnnData gene-ID column or `var_names` |
| `strip_version` | `--strip-version` / `--no-strip-version` | off | Drop the suffix after the first `.` |
| `gene_id_case` | `--gene-to-upper` / `--gene-to-lower` | keep | Normalize gene-ID case |
| `normalize` | `--normalize` | `cpm` | `cpm`, `tpm` or `none` |
| `sc_lengths_col` | `--sc-lengths-col` | — | Selected AnnData gene metadata column for TPM |
| `bulk_lengths_col` | `--bulk-lengths-col` | — | Bulk gene-length column for TPM; removed from expression |
| `discriminability` | `--discriminability` | off | Compute all pseudo-bulk × bulk group comparisons |
| `keep_missing_groups` (MCP only) | `--keep-missing-groups` | off | Retain missing group labels |
| `no_figure` | `--no-figure` | off | Write tables without figures |
| `no_manifest` | `--no-manifest` | off | Skip the provenance manifest |
| `timeout_ms` | — | 1800000 | Child-process timeout in milliseconds; configurable in DSH |

`strip_version` truncates at the first period, so use it only when that suffix
is a version identifier. IDs that become identical are summed on both the
single-cell and bulk sides and reported.

## Normalization

| Mode | Operation before log transform |
|---|---|
| `cpm` | Divide each vector by its sum and multiply by 1e6 |
| `tpm` | Divide by gene length, then scale the rates to sum to 1e6 |
| `none` | Leave values unchanged |

Scaling uses the shared genes positive on at least one side, before the final
both-positive filter. Different pairs can therefore use different gene sets
and denominators.

CPM scaling followed by `log2(x+1)` can change Pearson `r`, slope and
intercept when the denominator changes. For an unchanged gene set, positive
scalar scaling followed by a strictly increasing log transform preserves ranks
and Spearman `rho`, apart from numerical precision. TPM divides by
gene-specific lengths and can change ranks and both correlations. Results
from different modes are not interchangeable.

TPM requires `--sc-lengths-col` and `--bulk-lengths-col` with finite,
positive lengths. When normalized IDs merge, all lengths for the same ID must
agree; otherwise the run fails. The single-cell length column comes from
`raw.var` when `--sc-use-raw` is selected. Lengths must describe comparable
features and use consistent units. Do not apply length normalization again to
inputs already corrected for gene length.

## Group-matching diagnostics

With `--discriminability`, every pseudo-bulk group is compared with every
bulk group:

- `matched_mean_r`: mean Pearson r for groups with the same label.
- `mismatched_mean_r`: mean Pearson r for different labels.
- `gap`: matched mean minus mismatched mean.
- `n_correct / n_groups`: groups whose best-matching bulk group has the
  expected label.

The gap and best match flag unexpected pairings. A high correlation or positive
gap does not prove biological agreement. This analysis pools cells and
averages bulk samples by group; it does not estimate sample-to-sample
variation. Correlation p-values across genes do not provide inference across
independent biological replicates.

## Outputs and provenance

| File | Contents |
|---|---|
| `<prefix>_results.json` | Schema `dsh.bulk-scrna-concordance/1`; parameters, statistics, diagnostics, warnings and environment; CLI requires `--results-json` |
| `<prefix>_bulk_scrna_correlation_stats.csv` | One row per group |
| `<prefix>_discriminability.csv` | Pseudo-bulk × bulk Pearson matrix, when requested |
| `<prefix>_bulk_scrna_correlation.{png,pdf,svg}` | Scatter plots with regression and y=x lines |
| `<prefix>_manifest.json` | Script/input/output hashes, command, timestamps and versions |
| `<prefix>_run.log` | Child stdout/stderr when called through DSH or MCP |

The `environment_lock` field records versions, including Python, NumPy,
pandas, SciPy, AnnData, h5py and Matplotlib. It is not an installable dependency
lockfile. Reproducing a run requires the actual inputs and code plus a
compatible runtime; recorded versions alone do not guarantee byte-identical
figures across fonts, platforms or rendering environments.

## Checks and failures

Run from the repository root:

```bash
python -m unittest discover -s tests -v
node --check index.js
```

The tests use synthetic inputs. They check numerical and input-handling
behavior rather than biological validity.

The CLI exits with code 0 on success and 1 on analysis failure; argument-parser
errors can exit with code 2. In DSH, failed or timed-out runs return a structured
`failed` or `killed` result with the log tail. MCP reports a tool error with
diagnostics; each call writes to a new subdirectory within `out_dir`.

