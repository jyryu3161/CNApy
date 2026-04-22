# CNApy Changelog

This document records all additions and changes made in this fork on top of the original CNApy [[cnapy-org/CNApy v1.2.7]](https://github.com/cnapy-org/CNApy), listed in reverse chronological order.
Change types follow the [Keep a Changelog](https://keepachangelog.com/) convention: **Added / Changed / Fixed / Removed**.

---

## [Unreleased]

> This fork does not yet have a formal release (no git tag or PyPI distribution).
> All entries below are changes accumulated on `master`, grouped into feature milestones
> in reverse chronological order for readability.

### LLM-based Strain Analysis removed — 2026-04-22

#### Removed
- **LLM / agent layer removed**: The ChatGPT / Gemini-based strain analysis feature originally introduced under "AI integration" in "Initial extensions — 2025-12" has been removed. It was dead code — no menu in the GUI ever instantiated it, so end users could not reach it.
  - Files removed: `cnapy/gui_elements/llm_analysis_dialog.py`, `cnapy/gui_elements/agent_dialog.py`, entire `cnapy/agents/`, entire `cnapy/tests/test_agents/`
  - The `[llm]` and `[ai-agents]` optional dependencies in `pyproject.toml` and the `requires_crewai` pytest marker were dropped as well.

### GECKO dead dialog cleanup — 2026-04-21

#### Removed
- **GECKO dead dialog files removed**: Six dialog files that had been superseded by the unified dialog were deleted (−1,774 lines).
  - `gecko_dialog.py`, `enzyme_usage_dialog.py`, `proteomics_dialog.py`, `kcat_whatif_dialog.py`, `unikp_dialog.py`, `unikp_predictor.py`

#### Changed
- Corrected stale docstrings in `ecmodel_data.py` to reference `gecko_unified_dialog`.

### GECKO enzyme-constrained modeling — 2026-04-20

#### Added
- **GECKO 3.0 ecModel build**: *Analysis > "Enzyme-Constrained Model (GECKO)…"*
  - Both Full and Light formulations are supported.
  - GECKO 3.0 sign convention (v2 sign convention) is implemented correctly.
  - Handles isozymes and pseudoreactions.
  - Automatically generates protein-pool exchange and usage reactions.
  - Files: `cnapy/ecmodel/ecmodel_builder.py`, `expansion.py`, `ec_structure.py`, `ecmodel_data.py`, `exceptions.py`
- **Unified workflow dialog**: A single sidebar-navigation dialog integrates four pages.
  - Build ecModel / Enzyme Usage Report / Proteomics Integration / kcat What-if Analysis
  - File: `cnapy/gui_elements/gecko_unified_dialog.py`
- **GECKO YAML save/load**: `!!omap`-tagged ordered mapping, supports both plain GEM and ecModel YAML, persisted inside the `.cna` project ZIP as `ecmodel_data.json`.
  - Includes automatic v1→v2 sign-convention migration.
  - File: `cnapy/ecmodel/yaml_io.py`
- **File > New project from YAML**: start a new project directly from a GECKO YAML file.
- **File > New project from GECKO example**: launch instantly with bundled Human-GEM / Yeast-GEM datasets.
  - Human-GEM: `human-GEM.yml`, DLKcat.tsv, uniprot.tsv
  - Yeast-GEM: `yeast-GEM.yml`, customKcats.tsv, uniprot.tsv
  - Files: `cnapy/data/examples/gecko/`
- **GECKO test suite** (124 tests): Stage 1 signs/isozymes/pseudoreactions, YAML round-trip, CNA migration, paper parity, and bundled-example build verification.
  - Files: `cnapy/tests/test_ecmodel/`

#### Fixed
- Reject NaN kcat / MW values before they enter the stoichiometric computation.
- Tolerate externally produced GECKO YAMLs and prevent duplicate additions of usage reactions.
- Fix isozyme / pseudoreaction bugs and adopt the GECKO 3.0 sign convention.
- Align GECKO proteomics integration and kcat constraints with the paper formulation.

### UI Refinement — 2026-03-20

#### Changed
- Apply **FilterableComboBox** to the Flux Response Analysis (FRA) reaction selector to support substring search.

#### Fixed
- Prevent accidental value changes caused by mouse-wheel scrolling over comboboxes and spin boxes.
  - File: `cnapy/utils.py` (new `no_scroll` helper)

### Dynamic FBA improvements — 2026-03-19

#### Changed
- Set the default number of dFBA substrates to 3.
- Add a **reaction picker** to the dFBA dialog.
- Improve the dFBA table UI and add **bulk checkbox toggling**.
- File: `cnapy/gui_elements/dynamic_fba_dialog.py`

### Omics Gene Knockout — 2026-02-25

#### Added
- **MOMA/ROOM reference flux selection dialog**: a shared UI for selecting the reference flux during batch analyses.
  - File: `cnapy/gui_elements/moma_room_reference_dialog.py`

### Plot Customization — 2026-02-19

#### Added
- **Plot Customization dialog**: a common customization UI for `FigureCanvasQTAgg`-based plots.
  - Adjust title, axis labels, axis scale (log/linear), and axis ranges.
  - "Customize Plot" button integrated in 9 dialogs:
    - Flux Response Analysis, FSEOF, FVSEOF, Gene Essentiality, Robustness, Flux Sampling, Flux Optimization, Yield Optimization, and Phase plane / Yield space.
  - File: `cnapy/gui_elements/plot_customization_dialog.py`

### E-Flux2 and dependency cleanup — 2026-02-13

#### Added
- **E-Flux2 algorithm** (alongside LAD in Omics Integration): true L2 norm (QP)-based flux prediction with pFBA / FBA fallback.
  - Menu label: "Transcriptome-based Flux Prediction (LAD/E-Flux2)..."

#### Fixed
- **Recursive GPR-rule evaluation**: replace the previous flat gene aggregation with a proper recursive tree walk.
  - OR → max, AND → min traversal.
- **Restore missing dependencies**: add back dependencies from upstream CNApy to `pyproject.toml` (e.g., `gurobipy`).

### Multi-condition Omics and UI — 2026-02-05

#### Added
- **Multi-condition omics integration**: load multiple condition files simultaneously for comparative analysis.
- **Sorting and fold-change columns** in the omics results table.
- **Real-time text filtering** on reaction-selector comboboxes (substring search).
- **Map → Reactions tab synchronization**: clicking a reaction on the map automatically selects it in the Reactions tab.

#### Changed
- Make the omics-integration dialog **non-modal** so that other operations can run in parallel.

#### Fixed
- Fix the optlang variable-access error in LAD analysis.

### Strain Design extensions — 2026-01-22

#### Added
- **FSEOF Analysis** (Flux Scanning based on Enforced Objective Flux): identify reactions correlated with a target production flux to suggest over-expression / knockout targets for metabolic engineering.
  - File: `cnapy/gui_elements/fseof_dialog.py`
- **FVSEOF Analysis** (FVA-based FSEOF): perform FVA at each scan point for stricter target identification.
  - File: `cnapy/gui_elements/fvseof_dialog.py`
- **Batch MOMA/ROOM Analysis**: run MOMA / ROOM in batch across many knockout scenarios, with constraint reset support.
  - File: `cnapy/gui_elements/batch_moma_room_dialog.py`
- **Gene Essentiality Analysis**: systematic screening of gene essentiality.
  - File: `cnapy/gui_elements/gene_essentiality_dialog.py`
- **Robustness Analysis**: evaluate model robustness against parameter variations.
  - File: `cnapy/gui_elements/robustness_analysis_dialog.py`
- **MOMA/ROOM template flux selector**: UI for selecting the reference flux.
- **Reaction list enhancements**: added equation column and direction controls.
- **Analysis dialog UX improvements**.

#### Fixed
- Multiprocessing bug in Batch MOMA/ROOM.
- Remove quote text from background SVG files for maps.
- Replace residual console references with `print()`.

#### Removed
- Python console and related UI removed (refactor preparing for future extensions).

### Initial extensions — 2025-12

#### Added

##### Analysis features
- **ROOM (Regulatory On/Off Minimization)**: minimize flux changes after gene knockouts using a MILP solver.
  - Requirements: a MILP solver such as CPLEX, Gurobi, or GLPK.
  - File: `cnapy/moma.py`
- **Linear MOMA**: linear MOMA analysis with optional handling of external dependencies.
- **Flux Sampling**: GUI dialog for flux sampling.
  - Files: `cnapy/flux_sampling.py`, `cnapy/gui_elements/flux_sampling_dialog.py`
- **Flux Response Analysis**: scan the flux of a target reaction and plot the maximum production rate of a product.
  - File: `cnapy/gui_elements/flux_response_dialog.py`
- **Dynamic FBA (dFBA)**: FBA + ODE coupling for time-course simulation.
  - References: Mahadevan et al. 2002, Varma & Palsson 1994.
  - File: `cnapy/gui_elements/dynamic_fba_dialog.py`
- **Omics Integration (LAD)**: transcriptome-driven Least Absolute Deviation flux prediction.
  - Gene expression data loading (CSV / TSV / Excel).
  - Gene-to-reaction mapping (GPR rules).
  - Various aggregation methods (min / max / mean / sum).
  - File: `cnapy/gui_elements/omics_integration_dialog.py`
- **Configurable Auto Analysis method** (FBA / MOMA).

##### Model and scenario management
- **Model Management** tools:
  - GPR cleanup (automatically detect and consolidate duplicate genes).
  - Dead-end metabolite detection.
  - Blocked-reaction detection (FVA-based).
  - Orphan-reaction detection.
  - Model Validation (mass / charge balance, bound errors, etc.).
  - File: `cnapy/gui_elements/model_management_dialog.py`
- **External Flux Data Loading**: load reaction-flux data from CSV / TSV, compare multiple conditions, and display a Log2 fold-change heatmap (green = up, red = down).
  - File: `cnapy/gui_elements/flux_data_dialog.py`
- **Scenario Templates & Bookmarks** (Ctrl+T): predefined culture condition templates, quick knockout creation, and scenario bookmarks.
  - File: `cnapy/gui_elements/scenario_templates_dialog.py`
- **Media Management** (Ctrl+M): manage culture media configurations.
  - File: `cnapy/gui_elements/media_management_dialog.py`

##### AI integration
- **LLM-based Strain Analysis**: analyze reaction and gene plausibility using ChatGPT / Google Gemini.
  - Supports OpenAI GPT-4o and Google Gemini Flash.
  - Real-time information via web search.
  - Local storage of API keys and export of results to JSON / CSV.
  - File: `cnapy/gui_elements/llm_analysis_dialog.py`

  > **Note**: This feature was later removed in "LLM-based Strain Analysis removed — 2026-04-22".

##### Map features
- **Create maps from PNG / SVG images alone**: build CNApy maps from image files without a JSON file.
- **Custom reaction boxes**: add boxes for reaction IDs not in the model to display flux values on the map.
- File: `cnapy/gui_elements/central_widget.py`

##### UI/UX
- **Alt + Left-click or context menu to toggle reaction knockouts**.
- **Improved OptKnock explanation** (Strain Design dialog):
  - Outer Objective example: `EX_succ_e` (succinate production).
  - Inner Objective example: `BIOMASS` (growth).

#### Fixed
- Improved exception handling and code quality.
- Various minor bug fixes.

---

## [Base] Upstream CNApy 1.2.7 — 2025-11-04

The base of this fork: [`cnapy-org/CNApy` v1.2.7](https://github.com/cnapy-org/CNApy/releases/tag/v1.2.7).
Features that were part of the upstream release are not recorded in this changelog; please refer to the upstream release notes.

---

## License

These changes are distributed under the **Apache License 2.0** and are part of the original CNApy project. All modifications are compatible with the license of the upstream CNApy project.
