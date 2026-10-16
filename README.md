# Bulk–scRNA Concordance — DSH Host plugin

Turns the TBtools `Transcriptome Integration — Bulk & scRNA Correlation` plugin
into an agent-callable analysis tool in the DeepSeek Harness.

The original `.plugin` package contained a 118 MB PyInstaller-encrypted EXE plus
a TBtools Java GUI shell. This bundle replaces both: the numerics live in a
readable Python script inside the package, and the Harness registers a tool the
agent can call directly.

- Full parameter reference, input formats and interpretation: **[docs/usage.md](docs/usage.md)**
- Runnable input fixtures: **[examples/synthetic/](examples/synthetic/)** (600 genes × 6 cells, disjoint group markers)
- Survey of the neighbouring DSH plugin ecosystem: [docs/github-survey.md](docs/github-survey.md)
- CI runs the analysis on those fixtures on Linux and Windows: [.github/workflows/checks.yml](.github/workflows/checks.yml)

---

## What it adds

One tool: **`bulk_scrna_concordance`**.

For each group (L1/L2/L3 …) it sums single-cell counts into a pseudo-bulk
profile, CPM-normalizes both sides over the shared genes, applies `log2(x+1)`,
and reports Pearson r, Spearman rho, and a linear fit.

With `discriminability: true` it also correlates **every** pseudo-bulk group
against **every** bulk group and reports:

- matched vs mismatched mean r and their gap
- whether each pseudo-bulk group best-matches its own bulk group

Do this whenever the correlation is meant as evidence. A single high r cannot
show that a correlation is meaningful, because any two transcriptomes from the
same experiment correlate strongly; the discriminability block can.

## Files

| File | Role |
|---|---|
| `package.json` | bundle manifest (`dsh.bundle.patch`), display metadata, dependency declarations |
| `cordis.patch.yml` | inserts the Host plugin row and its config |
| `index.js` | registers the tool, validates arguments, launches the analysis, reads the result JSON |
| `python/bulk_scrna_concordance.py` | the analysis itself |
| `locale/{en,zh}.json` | Plugin Manager card text |

**This bundle declares no npm dependencies and imports nothing outside `node:`.**
A Host-only bundle needs none, and a dependency would make installation require
the npm registry. The tool definition is therefore a plain registry-ready object
with JSON Schema parameters, plus the local argument check that `defineTool`
would otherwise supply. `Config` is declared as JSON Schema for the same reason.

## Prerequisite: a Python interpreter with the scientific stack

The analysis needs **`anndata`, `h5py`, `scipy`, `pandas`, `numpy`,
`matplotlib`**. The interpreter that ships inside DSH has only `numpy` and
`pandas`, so it is *not* sufficient for `.h5ad` input.

Check a candidate interpreter:

```bash
python -c "import anndata, h5py, scipy, pandas, numpy, matplotlib; print('ok')"
```

Then point the plugin at it. The shipped `cordis.patch.yml` assumes `python` on
`PATH`; override it in your profile's own `cordis.patch.yml`:

```yaml
- id: bulk-scrna-concordance
  name: "@local/dsh-bulk-scrna-concordance"
  config:
    pythonPath: /usr/bin/python3        # or C:/Python313/python.exe
    defaultTimeoutMs: 1800000
    maxLogLines: 40
```

An absolute `pythonPath` is recommended: the Host's `PATH` is not guaranteed to
match your shell's.

If the interpreter is wrong the tool returns a clean failure naming the path and
listing the packages it needs — it does not throw.

## Install

```bash
git clone https://github.com/<you>/dsh-bulk-scrna-concordance
```

Then install the checkout into a profile. In the Harness Web UI:

1. Open **Settings → Plugins**.
2. Install a bundle, with the **absolute path of the directory containing
   `package.json`** as the target (e.g. `/home/me/dsh-bulk-scrna-concordance`
   or `C:\src\dsh-bulk-scrna-concordance`).
3. Read the returned `application` field. `applied` means it is live;
   `restart-required` means it is not yet.

From a shell, where a `dsh` CLI is available:

```bash
dsh plugin --profile <profile> add /abs/path/to/dsh-bulk-scrna-concordance
```

That installs the package, links it into `node_modules` and adds it to
`dsh.profile.bundles`. It writes **no row** into the profile's
`cordis.patch.yml`, so add the config row yourself if you need a non-default
`pythonPath` (see above).

Notes:

- `@local/...` is a deliberate scope: this package is meant to be installed from
  a path or a git checkout, not published to npm.
- The installer performs package installation and bundle selection. Do not
  hand-edit the profile's `package.json` or `cordis.patch.yml` to select a
  bundle, and do not run `pnpm` in the profile directory.
- A fresh bundle can activate through HMR. **Replacing an already-installed
  package requires a restart** to load the new JavaScript generation.
- The Desktop app owns the `desktop` profile and the `dsh` CLI refuses to write
  it; install from inside the app there.

## Verify

1. The tool row exists and is enabled (Plugin Manager, or `cordis_inspect_query`).
2. Run it on real data and confirm the tool result, not just the log:

   ```
   bulk_scrna_concordance with
     sc_h5ad        /path/to/scRNA_with_celltype.h5ad
     bulk_csv       /path/to/bulk_FPKM.csv
     group_col      Group
     discriminability true
     out_dir        <somewhere writable>
   ```

   Expected on the dataset this was developed against (verified — see below):

   | group | Pearson r | Spearman rho | n_genes |
   |---|---|---|---|
   | L1 | 0.8726 | 0.8606 | 21843 |
   | L2 | 0.8844 | 0.8782 | 22206 |
   | L3 | 0.8772 | 0.8748 | 22104 |

   and `discriminability`: matched 0.8781, mismatched 0.8213, gap 0.0568,
   `correct=2/3`, `misassigned: ["L1->L2"]`.

## What was fixed relative to the original

1. **`--strip_version` was a dead switch.** It was `action="store_true"` with
   `default=True`, so it could never be turned off. It is now
   `--strip-version / --no-strip-version`, default off.
2. **Duplicate gene IDs after version stripping crashed the analysis.** The
   original truncated at the first `.` with no dedup, then
   `pseudo_df.loc[common_genes, g]` returned a DataFrame for a duplicated label
   and died with
   `ValueError: operands could not be broadcast together with shapes (n,) (m,)`.
   Duplicate IDs are now collapsed explicitly (values summed) and the count is
   reported.
3. **Bulk sample-to-group matching was ambiguous.** `startswith`/`endswith`
   /substring tests silently mis-assigned groups that are prefixes of one
   another (L1 vs L10). Matching now uses exact delimited tokens, and unmatched
   samples are reported instead of dropped quietly.
4. **The `nan` group was treated as a real group.** In the original dataset 3440
   cells have no group label, and the original built a pseudo-bulk column for
   them. They are now excluded and reported.
5. **Added the discriminability assessment** described above.
6. **Added `--results-json`** so a caller reads one structured object instead of
   scraping stdout (stdout is a log file here, for the sandbox reason below).
7. **Run parameters are recorded** in the result JSON for reproducibility.
8. **The figure no longer emits a `tight_layout` warning.**

### One deliberate numeric change

CPM normalization is computed over the genes actually used in the correlation
(the shared, informative set), not over the whole shared gene list as in the
original. **Pearson r and Spearman rho are identical** to the original for this
reason — they are invariant to a linear rescaling of either axis. The fitted
`slope`/`intercept` are on a different scale: e.g. L1 slope is 0.901 here vs
0.845 in the original. Scale-free reporting should quote r, not the slope.

### Normalization is now explicit

`--normalize {cpm,tpm,none}` (tool argument `normalize`), default `cpm`:

- **`cpm`** — the historical behaviour, and what every published number above used.
- **`tpm`** — also divides by gene length. Needs `--sc-lengths-col` (an `.h5ad`
  `var` column) and `--bulk-lengths-col` (a bulk column, excluded from the
  expression matrix). Lengths must be present on both sides or the run stops,
  rather than silently falling back.
- **`none`** — inputs already on a common scale.

The choice is recorded in the results JSON and the manifest, because it changes
the fitted slope even though it barely moves r.

## Why stdout is a file rather than a pipe

Under the Harness file sandbox, child processes may not use piped stdio
(`spawnSync ... EPERM`); `stdio: 'inherit'` and raw file-descriptor redirection
both work. The plugin opens a descriptor on `<out_dir>/<prefix>_run.log`, gives
the child a duplicate of it, then reads the tail back for diagnostics. The
machine-readable result does not travel over stdout at all — it is read from the
JSON file the script writes.

## Outputs

Everything lands in `out_dir`:

| File | Content |
|---|---|
| `<prefix>_manifest.json` | **provenance**: script hash, input/output SHA-256, params, seed, environment lock |
| `<prefix>_results.json` | schema `dsh.bulk-scrna-concordance/1`: parameters, per-group stats, discriminability matrix, warnings, environment |
| `<prefix>_bulk_scrna_correlation_stats.csv` | one row per group |
| `<prefix>_discriminability.csv` | the pseudo × bulk Pearson matrix |
| `<prefix>_bulk_scrna_correlation.{png,pdf,svg}` | density scatter with regression and y=x lines |
| `<prefix>_run.log` | full analysis log |

### Provenance

Every run writes a manifest answering *"this figure = which code + which inputs +
which parameters + which environment"*:

- **script hash** — the analysis script's own SHA-256, so a figure produced by a
  modified script is detectable;
- **input hashes** — SHA-256, byte size and mtime for the `.h5ad` and bulk matrix.
  An input that changed silently between runs is otherwise invisible;
- **environment lock** — `python`, `numpy`, `pandas`, **`scipy`**, **`anndata`**,
  `h5py`, `matplotlib` and the seed. `anndata` decides how `.X`/`.raw` are read and
  `scipy` computes every statistic here, so those two are the versions that can
  actually move a reported `r`; recording only numpy/pandas would miss them;
- **output hashes** — every produced file, so a figure can be checked against the
  run that claims to have made it;
- **command line and timestamps**.

Nothing in this analysis is stochastic. The manifest records `seed: 42` anyway, so
that is a checkable claim rather than an assumption.

This follows the provenance model of
[`poplarity/dsh-science-workbench`](https://github.com/poplarity/dsh-science-workbench)
(see [docs/github-survey.md](docs/github-survey.md)), minus the project/cell/feedback
machinery that belongs to a workbench rather than to one analysis.

Pass `--no-manifest` to skip it.

### The figure does not appear in the conversation, by platform design

A tool result **cannot** carry a picture. The Harness attachment documentation is
explicit that the role-neutral `ImageBlock` leaves assistant-side image output as
forward compatibility, and that **current production adapters declare text-only
output, so only user content carries images**. An `ImageBlock` also needs an
`attachmentId` minted by the `attachments` service, and that service is not
mounted in every profile.

The render therefore emits the path plus a Markdown image reference, and states
that surfacing the figure as a card in the Web UI needs a separate `present`
call with the file path. Without that second call a user sees only a path and may
reasonably conclude no figure was produced at all.

## Verification status

Performed, and reproducible:

- **Python layer** — full run on the real 720 MB `.h5ad` and the real bulk
  matrix; per-group r reproduced to 1e-9 against an independent recomputation;
  discriminability matrix matches an independent 3×3 computation exactly.
- **Node layer** — the registered tool was exercised through a stub `ctx.tools`
  in five cases: invalid input combinations rejected; a bad `pythonPath`
  returned a clean `failed` object; a real analysis returned `status: ok` with
  the expected statistics; `render()` produced the expected model-facing text
  for both success and failure; and disposing the registration removes the tool.
- **Manifest** — `package.json` parses, all `exports` targets exist, every file
  in the bundle validates as UTF-8, `node --check` passes on `index.js`,
  `py_compile` passes on the script.
- **Composition against the real Loader** — the bundle's patch was applied to a
  scratch profile through `dsh --profile <name> --patch … --dump-config`; the row
  appears in the composed tree with its `id`, resolved module URL and full
  `config`, and the tree composes without error.
- **Installed into the `desktop` profile** with
  `dsh plugin --profile desktop add <bundle dir>` (pnpm `link:`, no registry
  access needed). Confirmed afterwards: the bundle name is in
  `dsh.profile.bundles`, the dependency is in `package.json`, and
  `node_modules/@local/dsh-bulk-scrna-concordance` resolves to this directory.
  The profile patch carries a matching `bulk-scrna-concordance` row whose row
  `id` overrides the bundle's default `config`, so `pythonPath` points at
  `D:/Python/python.exe`.
- **API contract read from the installed Harness**, not assumed:
  `ctx.tools.register()` was read at `dsh-tools/lib/index.js:2878` and requires
  `output.{schema,render}` plus a schema accepted by
  `assertSupportedJsonSchema` — it does **not** require a `defineTool()` result,
  which is what makes the dependency-free form valid. The supported keyword set
  is `type/oneOf/properties/required/additionalProperties/items/enum/const` plus
  `description/title/default/examples`; this bundle uses only those.
- **Live activation and a live tool call**, through a real Harness agent. A
  scratch headless profile (`base` + `dsh-headless` + this bundle) was booted
  with `dsh --profile <name> "<task>"`. The agent called
  `bulk_scrna_concordance` once on the real `.h5ad` and bulk matrix and reported
  `status: ok` with r = 0.873 / 0.884 / 0.877 and gap 0.0568 (`L1→L2`
  misassigned) — matching the independently computed values exactly. All six
  output files were written to the requested directory.

## Two traps worth knowing, both hit while installing this

1. **Do not export `Config` unless its schema library resolvable from the
   plugin.** Cordis resolves an exported `Config` through schemastery
   (`resolveConfig` calls `Config.validate`). `@deepseek-ai/schemastery` is not
   resolvable from a bundle installed by path, so an exported `Config` makes the
   plugin fail to activate with:

   ```
   bulk-scrna-concordance (@local/dsh-bulk-scrna-concordance):
   TypeError: Cannot read properties of undefined (reading 'validate')
   ```

   `Config` is optional, so it is omitted here and the row's config is read
   defensively in `apply()`. The symptom is a **silent** missing tool: the row
   composes, but only a boot-time warning says `1 entry did not activate`.

2. **A newly installed bundle needs a restart.** `dsh plugin --profile desktop
   add <path>` installs and selects the bundle, but writes no row into the
   profile's `cordis.patch.yml` (unlike `install_bundle`), and the running app
   does not pick up the new JavaScript generation. Both were needed here: the
   install, then a restart.

## Scratch profile used for verification

`~/.dsh/profiles/concordlive` — `dsh-base` + the shipped `dsh-headless` bundle +
this bundle, with a patch row supplying `pythonPath`. It is what produced the
live verification above, and is safe to delete:

```powershell
dsh plugin --profile concordlive remove @local/dsh-bulk-scrna-concordance
```

Note `dsh plugin … add @deepseek-ai/dsh-headless` installs the **registry**
copy, which is `0.0.1-rc.1` and is rejected as incompatible with dsh
`0.2.0-rc.2`; the shipped bundle inside the installation is `0.2.0-rc.2` and is
selected by name from `dsh.profile.bundles`, which is what this profile does.


`exec.agent.session.cwd` is read defensively for relative-path resolution; if
that field is absent in the installed version, paths fall back to the Host's
working directory, so passing absolute paths always works.
