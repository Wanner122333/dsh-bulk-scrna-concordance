# GitHub survey — DSH plugins adjacent to bulk–scRNA concordance

Survey method: GitHub search API (`web_fetch`, since the shell has no network),
queries `dsh-plugin {transcriptome,rnaseq,single-cell,bioinformatics}`,
`topic:dsh-plugin bioinformatics`, `topic:dsh-plugin rna`, plus
`topic:dsh-plugin` sorted by stars as the ecosystem baseline.

## Findings

**No plugin does bulk RNA-seq vs scRNA pseudo-bulk correlation.** The closest
domain neighbours:

| Repo | ★ | What it is | Overlap |
|---|---|---|---|
| [sagirimo/BioDSH](https://github.com/sagirimo/BioDSH) | 59 | Bioinformatics agent *desktop* (Tauri), skill store, offline mode | Packaging/distribution, not analysis |
| [poplarity/dsh-science-workbench](https://github.com/poplarity/dsh-science-workbench) | 11 | Reproducible science workbench: cells, inline figures, provenance, git | **Provenance/reproducibility model** |
| [moonbowterfly/dsh-bio-genie](https://github.com/moonbowterfly/dsh-bio-genie) | 3 | 57 semantic bio tools + self-bootstrapping private Python env | **Env bootstrapping, CN mirrors** |
| [gsh150801/dsh-bioinf-verify](https://github.com/gsh150801/dsh-bioinf-verify) | 2 | Bioinformatics verification plugin | Verification approach |
| [aiyacharley/dsh-pubmed](https://github.com/aiyacharley/dsh-pubmed) | 2 | 25 tools over PubMed/Europe PMC/PubTator3/Semantic Scholar | Tool-count/UX scale |
| [Sculptor815/scrna-seq-workbench](https://github.com/Sculptor815/scrna-seq-workbench) | 1 | Single-cell plugin with 6 skills + bundled knowledge base + benchmark | Domain, pseudobulk DE, KB pattern |

Ecosystem baseline (`topic:dsh-plugin`, stars descending) is desktop apps and
generic tooling — e.g. [deepseek-ai/deepseek-harness](https://github.com/deepseek-ai/deepseek-harness)
(250,610), which is the platform itself and the authoritative architecture
reference. The first-party figure/tool conventions in the installed runtime came
from there.

So the **highest-starred genuinely comparable plugin is `dsh-science-workbench`
(11★)**. Combined with `scrna-seq-workbench` (1★) it still gives the useful ideas,
because star count in this niche tracks age and marketing, not quality.

## What to adopt

Ranked by value to this plugin:

1. **Provenance per artifact** (science-workbench). Every output should record the
   producing script, input SHA-256, output SHA-256, params, seed and time, so
   "this figure = this code + these inputs + these params" is answerable. Our
   results JSON records parameters and input paths but **no hashes**, so an input
   that changed silently is undetectable.
2. **Environment lock** (science-workbench `environment.lock`). We record
   `python`, `platform`, `numpy`, `pandas` — but not `anndata`, `h5py` or
   `scipy`, which are exactly the versions that change numerical results.
3. **`manifest.json` as single source of truth** (science-workbench). Ours is one
   JSON per run; the idea worth taking is that it is the *only* thing a reader
   needs.
4. **Explicit normalization choice.** `--normalize {cpm,tpm,none}` — TPM needs
   gene lengths, but the biological default should be stated and recorded rather
   than implicit.
5. **Declared seed.** Nothing here is stochastic today, but recording a seed
   makes that an explicit, checkable claim.

## What NOT to adopt

- **A whole `bio_*` cell/notebook tool family.** Nine tools managing projects,
  cells and feedback are a workbench product; this plugin is one analysis.
  Adopting that surface would trade a working tool for a platform.
- **The Client/UI half.** science-workbench needs `lib/client.js` plus
  `webServer` routes. The figure is already written to disk and surfaced through
  `present`; a bespoke UI is disproportionate here.
- **Knowledge-base bundling.** scrna-seq-workbench ships a KB with per-run policy
  hashes — right for a 6-stage pipeline with annotation choices, overkill for one
  deterministic statistic.

## Where this plugin is ahead

Neither neighbour checks whether a correlation is *meaningful*. The
`discriminability` block (matched vs mismatched mean r, gap, best-match
correctness) is the part of this tool that a reviewer will actually ask about,
and it is not something to trade away for parity with either repo.
