# Bundled GECKO 3.0 example datasets

This folder ships ready-to-use GECKO ecModel inputs for two organisms so
CNApy users can build an enzyme-constrained model with zero external
downloads.

## Datasets

| Species | Directory | Reactions | kcat source | Size |
|---------|-----------|:--------:|-------------|-----:|
| Human (*Homo sapiens*) | `human/` | 11,983 | DLKcat predictions (NaN cleaned) | ~22 MB |
| Yeast (*S. cerevisiae*) | `yeast/` | 4,063  | BRENDA + DLKcat + manual curation | ~7 MB |

Each subdirectory contains:

* `*-GEM.yml` — the genome-scale metabolic model in GECKO YAML format
* `DLKcat.tsv` or `customKcats.tsv` — kcat table (`proteins/genes/gene_name/kcat/rxns/notes/stoicho`)
* `uniprot.tsv` — enzyme MW + sequence table (`Entry/Gene Names/.../Mass/Sequence`)
* `manifest.json` — metadata consumed by CNApy (`File → New project from GECKO example` and the GECKO dialog's "Example…" buttons)
* `LICENSE` — per-file upstream licenses

## Usage in CNApy

1. `File → New project from GECKO example → Human / Yeast` loads the GEM YAML.
2. `Configure → GECKO ecModel` opens the Build page.
3. The `kcat` and `UniProt` rows each have an **Example…** button. Click it,
   pick the species from the dropdown, and the bundled TSV loads automatically.
4. Choose **Full** or **Light** mode and press **Convert to ecModel**.

## Sources and licenses

| Asset | Upstream | License |
|-------|---------|---------|
| Human-GEM | <https://github.com/SysBioChalmers/Human-GEM> (v1.15.0) | MIT |
| yeast-GEM | <https://github.com/SysBioChalmers/yeast-GEM> (v8.6.2) | MIT |
| DLKcat predictions | <https://github.com/SysBioChalmers/DLKcat> | Apache 2.0 |
| UniProt extracts | <https://www.uniprot.org/> | CC-BY 4.0 |

Per-species `LICENSE` files carry the original text. This bundle is a
redistribution in compliance with those licenses; no modifications to the
scientific content were made beyond removing NaN kcat rows (DLKcat) that
would otherwise break ecModel builds on Gurobi/GLPK.

## Updating the bundle

Replace the files in place (keeping file names) and refresh the
`version` / `coverage` fields of each `manifest.json`. The GUI reads the
manifest at load-time, so no code change is needed for straightforward
updates.
