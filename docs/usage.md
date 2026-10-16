# Usage

The plugin registers one tool, `bulk_scrna_concordance`. Everything below is
also reachable from the shell by calling the analysis script directly — useful
for scripting and for reproducing a figure outside the Harness.

```bash
python python/bulk_scrna_concordance.py --help
```

## What it computes

For each group:

1. sum single-cell counts across that group's cells → one pseudo-bulk profile;
2. scale both sides (default CPM), take `log2(x + 1)`;
3. over the genes present and positive on both sides, compute Pearson `r`,
   Spearman `rho`, and a least-squares fit `y = slope·x + intercept`.

`r` answers "do these two look like the same transcriptome". It does **not**
answer "is this the right match" — see [Discriminability](#discriminability).

## Inputs

Two ways to supply the single-cell side.

### `.h5ad` (preferred for real data)

```bash
python python/bulk_scrna_concordance.py \
  --sc-h5ad scRNA.h5ad \
  --bulk-csv bulk_FPKM.csv \
  --group-col Group \
  --out-dir results --out-prefix my_run
```

- Gene IDs are read from a `var` column (`gene_symbols`, `symbol`, `gene_name`,
  … in that order) or fall back to `var.index`; override with `--gene-id-key`.
- `adata.raw` is used when present, otherwise `adata.X`. **Pseudo-bulk sums the
  matrix as-is**, so if `.X` holds normalized or log-transformed values rather
  than counts, the pseudo-bulk is wrong. Check this first on a new dataset.
- Cells whose group label is missing (`nan`, `NA`, `none`, …) are excluded and
  reported, rather than being analysed as their own group.

### Expression CSV + metadata CSV

```bash
python python/bulk_scrna_concordance.py \
  --sc-matrix-csv sc_expr.csv --sc-meta-csv sc_meta.csv \
  --bulk-csv bulk.csv \
  --out-dir results --out-prefix my_run
```

- `sc_expr.csv` — genes as rows, cells as columns, first column is the gene ID.
- `sc_meta.csv` — first column is the cell ID, and one column holds the group.
- A path ending in `.tsv`/`.txt` is read as tab-separated; **anything else is
  read as comma-separated**, whatever its extension says.

### Bulk matrix

Genes as rows, samples as columns, first column is the gene ID. Values may be
TPM/FPKM/CPM/counts — any continuous expression measure.

Samples are assigned to groups by **exact delimited token** match on the sample
name, not by substring. With groups `A` and `B`, both `A_1` and `Rtan_A_rep1`
map to `A`. Samples that match no group are listed in the log and in
`bulk.unmatched_samples`; they are never silently dropped.

## Parameters

| Tool argument | CLI flag | Default | Meaning |
|---|---|---|---|
| `sc_h5ad` | `--sc-h5ad` | — | AnnData input; exactly one of this or `sc_matrix_csv` |
| `sc_matrix_csv` | `--sc-matrix-csv` | — | genes × cells matrix |
| `sc_meta_csv` | `--sc-meta-csv` | — | required with `sc_matrix_csv` |
| `bulk_csv` | `--bulk-csv` | — | **required** |
| `out_dir` | `--out-dir` | — | **required**; created if absent |
| `out_prefix` | `--out-prefix` | `bulk_scrna_concordance` | output filename prefix |
| `group_col` | `--group-col` | `Group` | metadata column holding the group |
| `gene_id_key` | `--gene-id-key` | `auto` | `.h5ad` var column, or `auto` / `var_names` |
| `strip_version` | `--strip-version` / `--no-strip-version` | **off** | drop the suffix after the first `.` |
| `gene_id_case` | `--gene-to-upper` / `--gene-to-lower` | keep | case-normalize gene IDs |
| `normalize` | `--normalize` | `cpm` | `cpm` \| `tpm` \| `none` |
| `sc_lengths_col` | `--sc-lengths-col` | — | `.h5ad` var column of gene length (tpm) |
| `bulk_lengths_col` | `--bulk-lengths-col` | — | bulk column of gene length (tpm) |
| `discriminability` | `--discriminability` | off | also compute the full matrix |
| `no_figure` | `--no-figure` | off | tables only |
| `no_manifest` | `--no-manifest` | off | skip `manifest.json` |
| `timeout_ms` | — | plugin setting | kill the child after N ms |

`strip_version` defaults to off. Turn it on only when IDs genuinely carry a
version suffix; it truncates at the first `.`, so `ENSG00000141510.17` becomes
`ENSG00000141510` — but so would a gene literally named `AT1G01010.1`. IDs that
collapse onto each other are summed and reported, never silently dropped.

## Normalization

| Value | What it does | When |
|---|---|---|
| `cpm` | `x / sum(x) · 1e6` | default; comparable to the historical numbers |
| `tpm` | divides by gene length first | when length bias matters; needs both length columns |
| `none` | passes values through | inputs already on a common scale |

Scale affects the fitted `slope`, barely affects `r`/`rho`. Whatever you pick is
recorded in the results JSON and the manifest.

## Discriminability

Pass `--discriminability` and the run also correlates **every** pseudo-bulk group
against **every** bulk group, then reports:

- `matched_mean_r` — mean r of each group against its own bulk samples;
- `mismatched_mean_r` — mean r against the others;
- `gap` — the difference;
- `n_correct / n_groups` — how often a pseudo-bulk group's best-matching bulk
  group is its own.

**Read this before quoting an r.** Any two transcriptomes from the same
experiment correlate strongly, so a high `r` on its own is not evidence that
bulk and single-cell agree — it is largely a statement that both are plant (or
human, or whatever) transcriptomes. What carries information is the gap, and
whether the intended pairing actually wins.

Two very different outcomes on two real inputs:

| Input | matched | mismatched | gap | correct |
|---|---|---|---|---|
| Synthetic markers, `examples/synthetic/` | 0.9925 | −0.7784 | **1.7709** | 2/2 |
| The authors' Rheum dataset | 0.8781 | 0.8213 | **0.0568** | 2/3 |

In the second, `pseudo-L1` best-matches `bulk-L2`, not its own samples. A run
like that should not be reported as "bulk and scRNA agree"; the tool says so in
its warnings rather than leaving it to the reader.

## Outputs

| File | Content |
|---|---|
| `<prefix>_manifest.json` | script hash, input/output SHA-256, params, seed, environment lock |
| `<prefix>_results.json` | full result object, schema `dsh.bulk-scrna-concordance/1` |
| `<prefix>_bulk_scrna_correlation_stats.csv` | one row per group |
| `<prefix>_discriminability.csv` | pseudo × bulk Pearson matrix (with the flag) |
| `<prefix>_bulk_scrna_correlation.{png,pdf,svg}` | density scatter, regression and y=x lines |
| `<prefix>_run.log` | full stdout/stderr of the child process |

The figure is a file on disk. A tool result cannot carry a picture — Harness
production adapters declare text-only assistant output, so only user messages
hold images — which is why the tool prints the path and a Markdown image
reference and says that showing it as a card needs a separate `present` call.

## Reproducing a published figure

`<prefix>_manifest.json` is written to be the only file a reader needs:

```jsonc
{
  "produced_by": { "script_sha256": "...", "seed": 42 },
  "command":     { "argv": [ /* the exact invocation */ ] },
  "inputs": {
    "sc_h5ad":  { "bytes": 720779019, "sha256": "9e2f70e5..." },
    "bulk_csv": { "bytes": 4641941,   "sha256": "d3979353..." }
  },
  "environment_lock": {
    "python": "3.13.7", "scipy": "1.16.1", "anndata": "0.12.2",
    "numpy": "2.2.6", "pandas": "2.3.2", "h5py": "3.14.0", "matplotlib": "3.10.6"
  },
  "summary": { "L1": { "pearson_r": 0.8726, ... } },
  "outputs": { "....png": { "sha256": "c8271e23..." } }
}
```

`scipy` and `anndata` are in the lock because they decide the numbers: anndata
controls how `.X`/`.raw` are read, scipy computes every statistic. Recording
only numpy/pandas would miss both.

With the same inputs and the same environment the PNG is byte-identical across
runs, so `outputs.*.sha256` is a real check, not decoration. (PDF and SVG embed
timestamps and will differ.)

## Exit codes

| Code | Meaning |
|---|---|
| 0 | analysis completed |
| 1 | failed; the reason is the `[ERROR]` line in the log |

Inside the Harness the tool never throws for these: it returns a result object
with `status: "failed"` or `"killed"` plus the log tail, so a bad path or a
missing interpreter is a readable message rather than a stack trace.
