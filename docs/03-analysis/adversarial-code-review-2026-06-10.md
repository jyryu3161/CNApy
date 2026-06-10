# CNApy — Adversarial Code Review

**Date:** 2026-06-10  
**Scope:** `cnapy/` package (~37k hand-written LoC; generated `resources.py` excluded)  
**Method:** Static adversarial review — 52 module units, each reviewed by an independent skeptical agent, then every medium-or-higher finding re-checked by a second independent verifier agent to eliminate false positives. 5/5 sampled HIGH findings were additionally confirmed by hand (signatures, call sites, undefined attributes).  
**Tooling caveat:** Dynamic analysis (full `pytest` / `ruff`) was unreliable on this machine — the project `.venv/bin/python` is a symlink to the Anaconda base interpreter (not an isolated env), Gurobi's local license is version-mismatched (license v12.0 vs runtime v13.0), and example models are absent. Findings are therefore code-evidenced, not runtime-reproduced, except where noted.

## Result summary

- **236** raw findings → **195** confirmed after independent verification; **41** rejected as false positives.
- Confirmed by severity: **35 High**, **70 Medium**, **90 Low** (0 Critical — no silent-corruption-on-default-path or RCE).
- Recurring themes:
  1. **QThread lifecycle** — several analysis dialogs never stop/await their worker thread on close (`QThread: Destroyed while thread is still running` → process abort).
  2. **Unguarded `float()`/`int()`** on user-entered text (`ValueError` crashes).
  3. **Silently wrong scientific results** (objective / bounds / units) — the most dangerous class for a scientific tool, because no error is raised.
  4. **Method/attribute typos** that crash on real user actions (`dict` called as a function, undefined `pop_menu`, wrong arg counts).
  5. **YAML round-trip data loss** in the new `ecmodel` subsystem (objective coefficients dropped).

---

## A. Project-wide / environment issues

_Found via a repository-wide corpus scan, outside the per-module pass._

| # | Issue | Location | Severity | Notes |
|---|-------|----------|----------|-------|
| A1 | `numpy.bool` (removed in numpy ≥1.24) used as a dtype | `core.py:35`, `gui_elements/mode_navigator.py:269`, `flux_vector_container.py:23` | Medium (latent) | Works only because `pyproject.toml` pins `numpy==1.23`; any numpy bump → `AttributeError` on EFM computation, mode 'Select all', and flux-vector load. Use `bool` / `numpy.bool_`. |
| A2 | `PyYAML` imported but undeclared as a dependency | `ecmodel/yaml_io.py:36` vs `pyproject.toml` | Medium | Present only transitively in `uv.lock`; a clean install under a different resolver leaves the whole `ecmodel` feature un-importable. Add `pyyaml` to `[project.dependencies]`. |
| A3 | `subprocess.check_call(command, shell=True)` with interpolated user path | `configuration_gurobi.py:113-115`, `configuration_cplex.py:136-137` | Low | Local user-initiated solver setup; a path containing `"`/`&` breaks or mis-runs the command. Pass an argument list and drop `shell=True`. |
| A4 | `yaml.load(fp, Loader=_EcLoader)` | `ecmodel/yaml_io.py:270` | None (verified safe) | `_EcLoader` subclasses `yaml.SafeLoader`; **not** unsafe deserialization. Listed only to pre-empt a common false positive. |

---

## B. High-severity findings (35)

_Wrong results on a real path, or a crash on a real user action._

#### H1. compute_color_heat passes 0-255 integers to QColor.fromRgbF which expects 0.0-1.0, collapsing the flux heatmap

- **Location:** `cnapy/appdata.py:224-231` · **Category:** logic · **Verifier confidence:** high
- **Problem:** h is computed as round(mean*255/high) (or .../low), yielding a value in the 0-255 range. Lines 225 and 231 then call QColor.fromRgbF(255 - h, 255, 255 - h) / QColor.fromRgbF(255, 255 - h, 255 - h). fromRgbF is the FLOAT variant of the API and expects each channel in [0.0, 1.0]; values >1.0 are clamped to 1.0. Since 255 - h is >= 2 for essentially all h < 254, the red/blue channel is always clamped to full (1.0 = white). The intended green-to-white / red-to-white heat gradient therefore collapses: almost every reaction is rendered near-white regardless of its flux magnitude, so the heatmap conveys no information. git blame confirms this is a regression: commit 46e7a36b changed the original QColor.fromRgb(255-h, ...) (the integer 0-255 API, correct) to QColor.fromRgbF(...) without rescaling h, and the positive branch (line 225) was changed to fromRgbF in the most recent commit cc98275f.
- **Trigger:** User enables heat-map flux coloring on a map/reaction list (compute_color_heat is called from central_widget.py lines 1125-1150). Any model with computed flux values renders all reactions near-white instead of a graded color.
- **Fix:** Use the integer API: `return QColor.fromRgb(255 - h, 255, 255 - h)` and `return QColor.fromRgb(255, 255 - h, 255 - h)` (matching the original pre-regression code and the fromRgb calls in compute_color_onoff). Alternatively keep fromRgbF but rescale: divide each channel by 255.0, e.g. `QColor.fromRgbF((255 - h) / 255.0, 1.0, (255 - h) / 255.0)`.

#### H2. Isozyme splitting mis-handles complexes that contain a nested OR group, under-splitting and emitting non-gene tokens

- **Location:** `cnapy/ecmodel/expansion.py:186-195` · **Category:** logic · **Verifier confidence:** high
- **Problem:** parse_gpr_to_isozyme_sets only descends one level: it splits the top-level on OR, strips one outer paren pair from each part, then splits on AND. Any OR group nested inside an AND complex is never expanded. For a GPR like '(g1 or g2) and g3' (extremely common in genome-scale models such as yeast-GEM), the top level has no OR, so _split_top_level on OR returns the whole string as one part, _strip_outer_parens does nothing (the outer chars are not a matched pair), and the AND split yields ['(g1 or g2)', 'g3']. The result is a SINGLE isozyme set [['(g1 or g2)', 'g3']] in which '(g1 or g2)' is treated as one gene/protein token. GECKO 3.0 expandModel requires this to be distributed into two isozyme reactions catalysed by complexes (g1 and g3) and (g2 and g3). Consequences: (1) expand_isozymes sees len(isozyme_sets)==1 and does NOT split (line 289), so the reaction keeps a multi-isozyme GPR that later kcat application treats as a single enzyme complex — wrong enzyme-cost constraint; (2) the literal token '(g1 or g2)' is not a valid UniProt/gene ID, so downstream protein/kcat mapping finds no enzyme for it, silently dropping a real catalytic alternative. This is silently-incorrect scientific output, not a crash.
- **Trigger:** expand_isozymes (or parse_gpr_to_isozyme_sets) is called on a reaction whose gene_reaction_rule embeds an OR inside an AND complex, e.g. '(g1 or g2) and g3' — a standard pattern in yeast-GEM / iML1515-class models.
- **Fix:** Implement a real recursive GPR distribution to disjunctive normal form: parse the rule into an AST (or recurse), expand each AND of (OR-groups) into the Cartesian product of its alternatives, so '(g1 or g2) and g3' yields [['g1','g3'],['g2','g3']]. cobra already provides this via cobra.core.gene.GPR / ast2str; reusing GPR.from_string(...).body and walking the AST avoids the hand-rolled one-level parser entirely.

#### H3. Reaction objective coefficients are silently lost on YAML round-trip

- **Location:** `cnapy/ecmodel/yaml_io.py:117-132, 325-355` · **Category:** data-loss · **Verifier confidence:** high
- **Problem:** save_ecmodel never serialises each reaction's objective_coefficient (the objective function), and load_ecmodel never restores it. cobra's own save_yaml_model/load_yaml_model persists objective_coefficient per reaction, so a model loaded back through this module has an empty/zero objective. After a save->load round-trip, FBA on the reloaded model optimises a null objective (all-zero), producing wrong/degenerate flux solutions with no error raised. This is the central scientific output of the tool.
- **Trigger:** Save any ecModel/GEM whose objective is a normal reaction (e.g. biomass) via save_ecmodel, then reload with load_ecmodel and run FBA — the objective is gone and optimization returns a meaningless solution.
- **Fix:** In _reaction_to_omap append ('objective_coefficient', _as_number(rxn.objective_coefficient)) when nonzero, and in load_ecmodel set rxn.objective_coefficient = _as_float(r.get('objective_coefficient', 0.0)) (or set model.objective from the collected coefficients) before/after add_reactions.

#### H4. Bound-swap fabricates an invalid flux interval when the reference range does not intersect the model's existing bounds, silently relaxing hard constraints

- **Location:** `cnapy/flux_sampling.py:134-139` · **Category:** math · **Verifier confidence:** high
- **Problem:** When the requested constraint interval (from the reference flux / near-zero rule) has an empty intersection with the reaction's existing bounds, new_lb > new_ub. Instead of recognizing this as infeasible, lines 135-136 simply swap the two values and assign them as the reaction's bounds. The swapped interval can lie partly or wholly OUTSIDE the model's original feasible bounds, so a hard model constraint is silently relaxed. The applied_bounds dict (returned to the GUI and reported to the user) then advertises bounds the original model never permitted. Sampling proceeds over a feasible space that is scientifically wrong.
- **Trigger:** constraint_mode=='bounds' and a reference flux whose [min_fraction*flux, max_fraction*flux] (or [-0.1,0.1] near-zero) window does not overlap the reaction's existing lower/upper bounds — e.g. near-zero reference on a reaction with a strictly-positive lower bound, or a reference flux of opposite sign to the reaction's feasible direction.
- **Fix:** Treat new_lb > new_ub as an empty/infeasible intersection rather than swapping. Either skip constraining that reaction (leave its original bounds) and record a warning, or clamp the window into the existing bounds (new_lb = clamp(ref_lb, rxn.lower_bound, rxn.upper_bound); new_ub = clamp(ref_ub, rxn.lower_bound, rxn.upper_bound)) so the result always stays within [rxn.lower_bound, rxn.upper_bound]. Never produce a range outside the model's own bounds.

#### H5. Failed .npz load leaves a half-constructed object that crashes callers with AttributeError

- **Location:** `cnapy/flux_vector_container.py:15-22` · **Category:** exception · **Verifier confidence:** high
- **Problem:** In FluxVectorContainer.__init__, when numpy.load (or the subsequent dict access) raises, the except block shows a QMessageBox and `return`s from __init__ BEFORE any of self.fv_mat, self.reac_id, self.irreversible, self.unbounded are assigned. Python still returns the half-constructed instance. Callers (main_window.load_modes/load_mcs) immediately assign it to self.appdata.project.modes and call set_to_efm()/update_mode(), which access modes via __len__ (self.fv_mat.shape[0]) and modes.reac_id. Because fv_mat was never set, this raises AttributeError: 'FluxVectorContainer' object has no attribute 'fv_mat', and the previously valid project.modes has already been overwritten -> corrupted application state, not a graceful error.
- **Trigger:** User selects File > Load EFMs/MCS (load_modes/load_mcs) and picks a .npz that is corrupted or not a valid EFM file (e.g. an unrelated .npz). The except branch fires, the broken object is stored in project.modes, and the immediately-following update_mode/set_to_efm call raises AttributeError.
- **Fix:** Do not leave the object half-built. Either re-raise after notifying the user, or initialize the object to a valid empty state in the except branch (self.fv_mat = numpy.zeros((0,0)); self.reac_id = []; self.irreversible = numpy.array(0); self.unbounded = numpy.array(0)) before returning, and have callers check for an empty/invalid result before overwriting project.modes.

#### H6. Map template silently uses FVA/bounds lower bound as reference flux (scientifically wrong)

- **Location:** `cnapy/gui_elements/batch_moma_room_dialog.py:131` · **Category:** math · **Verifier confidence:** high
- **Problem:** When template_type=='map', the worker builds reference_fluxes via {rid: vals[0] for rid, vals in self.comp_values.items()}, treating each comp_values entry as a (flux, precision) pair and taking element [0] as the flux. But comp_values can also hold FVA / 'show model bounds' results where each entry is a (lower_bound, upper_bound) tuple (appdata sets project.comp_values_type = 1 in that case; see main_window.py show_model_bounds line 2688 and the fva path). In that mode vals[0] is the LOWER BOUND, not an actual flux. MOMA/ROOM then minimize distance to the lower-bound vector instead of a real flux distribution, producing silently incorrect results. The dialog never inspects comp_values_type: _update_template_availability (line 570) only checks bool(self.appdata.project.comp_values), so the 'Current Map' option is enabled and selectable even when the stored data is FVA bounds.
- **Trigger:** User runs FVA or 'Show model bounds' (comp_values_type becomes 1), then opens this dialog and selects Template Flux = 'Current Map'. Analysis runs without error but uses lower bounds as the reference flux distribution.
- **Fix:** In _update_template_availability, disable/relabel the 'Current Map' option when self.appdata.project.comp_values_type != 0; and/or in the worker guard `if self.appdata.project.comp_values_type != 0: emit error` before line 131. Only use Map fluxes for simple flux-vector results.

#### H7. Worker QThread not stopped/awaited on dialog close — destroyed-while-running crash

- **Location:** `cnapy/gui_elements/batch_moma_room_dialog.py:558-559` · **Category:** concurrency · **Verifier confidence:** high
- **Problem:** The Close button connects to self.accept and there is no closeEvent/reject override, no request_cancel(), and no worker_thread.wait() anywhere. If the user clicks Close (or the window 'X') while a batch run is in progress, the QDialog and its slots can be torn down while BatchMomaRoomWorkerThread is still executing model.optimize()/linear_moma/room. The running thread continues to emit progress_update/result_ready/error_occurred into slots of a dialog being destroyed, and Qt raises 'QThread: Destroyed while thread is still running', which typically aborts/crashes the process. The thread also keeps mutating its model copy and references self.appdata after the dialog is gone.
- **Trigger:** Start a long batch analysis (large model, MOMA/ROOM over all genes/reactions) and close the dialog before it finishes.
- **Fix:** Override closeEvent: if self.worker_thread and self.worker_thread.isRunning(): call request_cancel(), then self.worker_thread.wait() (optionally with a timeout) before accepting; or ignore the close until the thread finishes. Ensure the thread is joined before the dialog is destroyed.

#### H8. Box position dialog does item-assignment on a tuple, crashing after a box has been dragged

- **Location:** `cnapy/gui_elements/box_position_dialog.py:72-73` · **Category:** crash · **Verifier confidence:** high
- **Problem:** set_position() mutates the stored box position via index assignment: maps[name]['boxes'][id][0] = x_float and [1] = y_float. The box position is stored in two different forms depending on how it was last set: map_view.py stores it as a tuple when the box is dragged (lines 107, 117, 121: '= (point_item.x(), point_item.y())'), and as a list when freshly deserialized from a saved project's JSON (JSON arrays load as lists). Tuples do not support item assignment, so once the user drags a box (or any session where a box position became a tuple), opening the 'Set reaction box position' dialog and clicking 'Set position' raises TypeError: 'tuple' object does not support item assignment, aborting the operation and leaving the box unmoved.
- **Trigger:** User drags a reaction box on a map (storing its position as a tuple), then opens the box's 'Set reaction box position' dialog and clicks 'Set position'.
- **Fix:** Replace the two element assignments with a single whole-value reassignment to a new tuple/list, e.g.: self.map.appdata.project.maps[self.map.name]['boxes'][self.reaction_box.id] = (x_float, y_float). Also standardize the stored type across map_view.py and the dialog.

#### H9. Gurobi setup success path calls undefined method get_and_set_environmental_variable -> AttributeError

- **Location:** `cnapy/gui_elements/configuration_gurobi.py:135` · **Category:** crash · **Verifier confidence:** high
- **Problem:** After the Gurobi setup.py connection script runs successfully, the success branch calls self.get_and_set_environmental_variable(). This method is NOT defined anywhere in GurobiConfigurationDialog (it has only __init__, folder_error, choose_gurobi_directory, run_python_connection_script). The method exists only in configuration_cplex.py. Calling it raises AttributeError, so the success path crashes precisely when the user has done everything correctly and Gurobi was actually configured.
- **Trigger:** User sets a valid Gurobi folder, presses 'Run Python connection script', and subprocess.check_call returns 0 (setup succeeded). The success message shows, then the call to the missing method raises AttributeError, aborting and leaving the dialog/operation in an error state.
- **Fix:** Either implement get_and_set_environmental_variable() in GurobiConfigurationDialog (mirroring the CPLEX one for Gurobi env vars) or remove/guard the call. If env var handling is intended, factor it into a shared base class or module-level helper.

#### H10. Growth rate and biomass flux taken from objective_value, not the selected biomass reaction

- **Location:** `cnapy/gui_elements/dynamic_fba_dialog.py:243-251` · **Category:** math · **Verifier confidence:** high
- **Problem:** run_dfba never sets the model objective to params.biomass_reaction. It uses `mu = solution.objective_value`, the value of whatever objective the (scenario-loaded) model already has, then uses mu for dX/dt = mu*X (line 259) and records it as the biomass reaction's flux (line 251). The biomass_combo is editable and lists multiple candidates plus the current objective; if the user picks a biomass reaction that is not the model's objective (or the loaded scenario set a different/zero objective), the growth rate used to integrate biomass and the recorded biomass flux are simply wrong (could be a non-growth objective or 0). The simulation then reports physically incorrect biomass trajectories without any warning.
- **Trigger:** User selects (or types) a biomass reaction that differs from the model's current objective, or load_scenario_into_model changes the objective; growth dynamics are computed from the wrong reaction.
- **Fix:** Inside the `with model as m:` block, set `m.objective = m.reactions.get_by_id(params.biomass_reaction)` before optimize, and use `mu = m.reactions.get_by_id(params.biomass_reaction).flux` (or solution.fluxes[biomass_reaction]) rather than objective_value.

#### H11. modes_coloring left True after log2FC apply corrupts subsequent default flux coloring

- **Location:** `cnapy/gui_elements/flux_data_dialog.py:548` · **Category:** logic · **Verifier confidence:** high
- **Problem:** In the log2FC branch of _apply, self.appdata.modes_coloring is set to True and is never reset within this flow. Because the custom flux_value_display override is (incorrectly) reverted synchronously (see other finding), the actual rendering uses appdata.flux_value_display with modes_coloring=True. In that path (appdata.py:155-159) any reaction whose stored value is exactly 0 is colored red and every nonzero value is colored green — a binary mode-coloring scheme, not a log2FC gradient. Worse, modes_coloring stays True globally after the dialog acts, so all later normal flux displays elsewhere in the app are wrongly colored until something else resets it. _update_visualization (single mode) sets it back to False, but the log2FC branch never does.
- **Trigger:** Apply log2 fold change to the map; then view any flux values anywhere in the app — they are colored with the binary modes scheme instead of normal flux coloring.
- **Fix:** Reset appdata.modes_coloring = False after the log2FC visualization completes (mirroring _update_visualization), and do not rely on modes_coloring for log2FC gradient coloring at all.

#### H12. UnboundLocalError: fixed_growth_rate used outside the block that assigns it (biomass-only adjustment crashes)

- **Location:** `cnapy/gui_elements/flux_feasibility_dialog.py:338-373` · **Category:** crash · **Verifier confidence:** high
- **Problem:** fixed_growth_rate is assigned only at line 338, which lives inside two nested guards: `if len(reactions_in_objective) > 0:` (line 323) AND `if bm_is_modified and self.bm_mod_scenario.isChecked():` (line 337). It is later READ at lines 372-373 under a different, weaker guard: `if bm_is_modified:` (line 345) -> `if self.bm_mod_scenario.isChecked():` (line 368), which does NOT require reactions_in_objective to be non-empty. When the biomass equation is modified but NO flux corrections were applied, reactions_in_objective is empty, so line 338 never executes and line 372 raises UnboundLocalError.
- **Trigger:** Uncheck the 'Allow corrections to given fluxes' group (so flux_weight_scale=0 and reactions_in_objective comes back empty from make_scenario_feasible), enable biomass adjustment, check 'Add modified biomass reaction to scenario', and run Compute on a scenario where the solver returns an optimal solution that modifies the biomass equation. This is a fully supported usage (adjust only the biomass equation, not fluxes).
- **Fix:** Compute fixed_growth_rate where bm_reac_id is established (e.g. right after `bm_reac_id = self.bm_reac_id` near line 250, as `fixed_growth_rate = self.appdata.project.scen_values[bm_reac_id][0]`), or assign it once inside the `if bm_is_modified:` block before it is used at lines 372-373, independent of reactions_in_objective.

#### H13. UnboundLocalError: gam_base returned without being assigned when GAM metabolites are invalid

- **Location:** `cnapy/gui_elements/flux_feasibility_dialog.py:536-544` · **Category:** crash · **Verifier confidence:** high
- **Problem:** get_gam_removal_parameters assigns gam_base only inside `if valid:` (line 539) after validate_gam_mets. If validate_gam_mets returns valid=False (a metabolite in the gam_mets_edit text is not part of the biomass reaction), the `if valid:` block is skipped entirely and gam_base is never bound, yet line 544 unconditionally returns `valid, gam_mets, gam_base`, raising UnboundLocalError before the caller can inspect `valid`.
- **Trigger:** The gam_mets_edit field contains a metabolite id that is not in the selected biomass reaction (e.g. set programmatically, pre-populated, or entered without re-triggering the editingFinished isModified validation), then Compute (line 278) or any path that calls update_bm_constituents_table -> get_gam_removal_parameters (line 402, e.g. toggling the GAM group) executes while adjust_gam is checked.
- **Fix:** Initialize `gam_base = 0` at the top of get_gam_removal_parameters (and `gam_mets`/return early) so the function returns a well-defined tuple when validation fails, letting callers handle `valid is False` as intended.

#### H14. Worker QThread is never stopped/awaited when the dialog closes, leaving a thread emitting signals into a destroyed dialog

- **Location:** `cnapy/gui_elements/flux_response_dialog.py:464-470, 626-641` · **Category:** concurrency · **Verifier confidence:** high
- **Problem:** The dialog wires the worker thread's signals (progress_update, result_ready, error_occurred, finished) to its own slots, but there is no closeEvent/reject override and no call to request_cancel/quit/wait/deleteLater anywhere (grep for closeEvent/reject/wait/terminate/quit/deleteLater returns nothing). The Close button is connected to self.accept(). If the user closes the dialog while an analysis is running, the QDialog (and its Qt C++ object) is destroyed while FluxResponseWorkerThread.run() is still executing and still holds connected signals. When the worker later emits result_ready/error_occurred/finished, Qt delivers them to slots on a deleted object, which causes 'wrapped C/C++ object has been deleted' RuntimeError or a hard crash. The running optimization (model.optimize()) also continues consuming the solver and CPU after the dialog is gone.
- **Trigger:** Start an analysis (Run Analysis) on a non-trivial model, then click Close (or press Esc/window X) before it finishes. The worker keeps running and emits into the destroyed dialog.
- **Fix:** Override closeEvent(): if self.worker_thread and self.worker_thread.isRunning(): self.worker_thread.request_cancel(); self.worker_thread.wait(timeout) (or block close until finished). Disconnect signals / call deleteLater on the worker in _on_finished, and reject/accept only after the thread has stopped.

#### H15. model.objective_direction is left at 'min' for the entire flux scan, maximizing nothing / minimizing biomass

- **Location:** `cnapy/gui_elements/fseof_dialog.py:135-169, 184` · **Category:** math · **Verifier confidence:** high
- **Problem:** The target-range computation runs a 'min' block (lines 135-142) and then a 'max' block (lines 145-152), each inside 'with model:'. cobra reverts objective_direction set inside a context on exit, so after both blocks the direction returns to the model's pre-block value. At line 169 only model.objective is reassigned to obj_rxn (biomass); objective_direction is NOT explicitly set. The scan loop (line 184) then calls model.optimize() relying on an implicit 'max'. This is only correct if the model's original objective_direction was 'max'. cobra models created from some formats / after prior operations can carry direction 'min'; in that case every scan FBA MINIMIZES biomass, producing systematically wrong fluxes and meaningless correlations, with no error raised.
- **Trigger:** Run FSEOF on any model whose solver objective direction was not 'max' at dialog start (e.g. a model previously configured to minimize, or imported with direction min). All scan optimizations then minimize biomass.
- **Fix:** Explicitly set model.objective_direction = "max" at line 169 right after model.objective = obj_rxn, so the scan deterministically maximizes the biomass objective regardless of prior state.

#### H16. Dialog can be closed while worker QThread is running, destroying a live QThread (crash)

- **Location:** `cnapy/gui_elements/fseof_dialog.py:255-268, 441-443` · **Category:** concurrency · **Verifier confidence:** high
- **Problem:** FSEOFDialog spawns an FSEOFWorkerThread (self.worker_thread.start() at line 540) but never overrides closeEvent/reject and never calls wait()/quit() on the thread before the dialog is destroyed. The Close button is wired to self.accept() (line 442), and the window-manager close box / Escape key trigger QDialog.reject(). Either of these closes and schedules destruction of the QDialog while the child QThread (parented implicitly to the dialog as its Python owner) is still executing model.optimize() in run(). Qt then emits 'QThread: Destroyed while thread is still running' and aborts the process.
- **Trigger:** Start an FSEOF analysis on a non-trivial model, then click the Close button (or press Escape / click the window X) before the scan finishes.
- **Fix:** Override closeEvent (and/or reject): if self.worker_thread and self.worker_thread.isRunning(): self.worker_thread.request_cancel(); self.worker_thread.wait(); then accept the close. Also keep a reference and call deleteLater on finished.

#### H17. XLSX export crashes with KeyError: 'up_candidates' / 'down_candidates'

- **Location:** `cnapy/gui_elements/fvseof_dialog.py:1059,1076,1102-1103` · **Category:** crash · **Verifier confidence:** high
- **Problem:** _export_xlsx has the same root cause as the CSV defect: it reads self.last_results["up_candidates"] (line 1059), self.last_results["down_candidates"] (line 1076), and len() of both in the Summary sheet (lines 1102-1103). Since these keys are never stored in last_results, every XLSX export raises KeyError. The whole export fails after the 'All Reactions' sheet has been built in memory; nothing is saved to disk and the user only sees 'Failed to export'.
- **Trigger:** User runs analysis successfully and clicks 'Export XLSX' (with openpyxl installed). Building the Up-regulation sheet hits self.last_results['up_candidates'] -> KeyError; wb.save is never reached, so no file is produced.
- **Fix:** Store up_candidates/down_candidates on self.last_results in _apply_cutoffs (see CSV fix), or recompute candidate lists inside _export_xlsx from regression_results using the cutoff spin values.

#### H18. Worker QThread not stopped/awaited on dialog close — abandoned thread keeps mutating a copied model and emitting to a destroyed dialog

- **Location:** `cnapy/gui_elements/fvseof_dialog.py:740-746,914-918` · **Category:** concurrency · **Verifier confidence:** high
- **Problem:** There is no closeEvent/reject override that cancels and waits for the worker thread. _cancel_analysis only sets the cancel flag and returns immediately (it does not call wait()). If the user clicks 'Close' (close_btn -> self.accept) or the window 'X' while the analysis is running, the QDialog is destroyed but the FVSEOFWorkerThread keeps running. The running thread continues to emit signals (progress_update, result_ready, etc.) bound to slots on the now-deleted dialog, which can cross-thread-touch deleted Qt C++ objects and crash the process, and the thread may outlive expectations holding a full copied COBRA model. The cancel flag is only polled at coarse points (lines 192, 233), so even cancellation takes a long time.
- **Trigger:** Start an FVSEOF run on a non-trivial model (long-running), then close the dialog (Close button or window close) before completion. The worker thread continues optimizing and emits result_ready/progress_update into deleted slots, risking a crash or silent resource retention.
- **Fix:** Override closeEvent/reject to call self.worker_thread.request_cancel() and self.worker_thread.wait() (with a timeout) before accepting the close, and disconnect signals; ensure the thread is fully stopped before the dialog is destroyed.

#### H19. CSV export crashes with KeyError: 'up_candidates' / 'down_candidates'

- **Location:** `cnapy/gui_elements/fvseof_dialog.py:979,987` · **Category:** crash · **Verifier confidence:** high
- **Problem:** _export_csv reads self.last_results["up_candidates"] and self.last_results["down_candidates"], but these keys are never written into last_results. The result dict emitted by the worker (lines 357-372) contains only 'regression_results' (plus metadata), and _apply_cutoffs computes up_candidates/down_candidates purely as local variables (lines 793, 806) without ever assigning them back to self.last_results. Accessing a missing key on a plain dict raises KeyError, which is not the generic 'export failed' situation the code anticipates — it aborts the entire export after the main regression CSV may already be partially written.
- **Trigger:** User runs an analysis successfully, then clicks 'Export CSV'. The first CSV (regression_results) is written, then 'up_filename' is opened and the KeyError fires while iterating self.last_results['up_candidates'], leaving partial/empty candidate files and showing a misleading 'Failed to export' message.
- **Fix:** In _apply_cutoffs, after computing the lists, persist them: self.last_results["up_candidates"] = up_candidates; self.last_results["down_candidates"] = down_candidates. Alternatively, recompute the candidate lists from regression_results inside the export functions using the current cutoff spin values.

#### H20. Right-click on gene list crashes: self.pop_menu is never defined

- **Location:** `cnapy/gui_elements/gene_list.py:75-77` · **Category:** crash · **Verifier confidence:** high
- **Problem:** on_context_menu calls self.pop_menu.exec_(...) but GeneList never creates a pop_menu attribute. The comment at line 46 ('create context menu') is a leftover; unlike metabolite_list.py (which builds self.pop_menu = QMenu(...) at line 86), GeneList has no such construction. The context-menu policy IS wired up (lines 43-44 connect customContextMenuRequested to on_context_menu), so the handler runs on every right-click.
- **Trigger:** User right-clicks (context menu) anywhere on the gene tree widget while the model has at least one gene. Raises AttributeError: 'GeneList' object has no attribute 'pop_menu'.
- **Fix:** Either remove the context-menu policy/connect (lines 43-44) and on_context_menu, or construct self.pop_menu (a QMenu with the intended actions) in __init__ before connecting, mirroring metabolite_list.py.

#### H21. InOutFluxDialog 'Plot' calls a non-existent method with wrong arg count -> always crashes

- **Location:** `cnapy/gui_elements/in_out_flux_dialog.py:42-44` · **Category:** crash · **Verifier confidence:** high
- **Problem:** compute() invokes self.appdata.window.centralWidget().in_out_fluxes(metabolite). centralWidget() returns a CentralWidget instance, but in_out_fluxes is a method of MainWindow, not CentralWidget (grep confirms there is no in_out_fluxes attribute anywhere in central_widget.py). This raises AttributeError ('CentralWidget' object has no attribute 'in_out_fluxes'). Even if it resolved to MainWindow's method, that method is defined as 'def in_out_fluxes(self, metabolite_id, soldict)' and requires TWO arguments, while compute() passes only one, so it would instead raise TypeError. The intended entry point is MainWindow.print_in_out_fluxes(metabolite) (line 2682), which builds soldict from comp_values and then calls in_out_fluxes(metabolite, soldict); that helper is defined but currently wired to nothing else, confirming the dialog was supposed to call it. As written, clicking the dialog's 'Plot' button always fails.
- **Trigger:** Open Compute in/out fluxes dialog, pick a metabolite, click 'Plot'. Raises AttributeError (and would be TypeError even if the receiver were correct).
- **Fix:** Call the helper on the MainWindow instead of CentralWidget and let it build soldict: replace line 44 with `self.appdata.window.print_in_out_fluxes(metabolite)`.

#### H22. Concentration loaders call _set_concentrations() without required replace_all arg -> TypeError

- **Location:** `cnapy/gui_elements/main_window.py:3071` · **Category:** crash · **Verifier confidence:** high
- **Problem:** _set_concentrations is defined as `def _set_concentrations(self, concentrations, replace_all)` (line 3028) and requires two arguments. Both call sites in the concentration loaders pass only one argument: `_load_concentrations_json` calls `self._set_concentrations(concentrations)` (line 3071) and `_load_concentrations_xlsx` calls `self._set_concentrations(concentrations)` (line 3135). Every one of the four 'Load concentration ranges [in M]' menu actions (JSON/XLSX x amend/replace-all, wired at lines 665-683) therefore raises `TypeError: _set_concentrations() missing 1 required positional argument: 'replace_all'`. The dG0 variants correctly pass both args (line 3145/3202), which confirms the concentration calls are the defect. As a secondary consequence, even if the call were fixed by defaulting, the replace_all flag selected by the user (amend vs replace-all) is silently dropped, so the two menu items would behave identically.
- **Trigger:** User invokes any of the Thermodynamics > 'Load concentration ranges [in M] (replacing/amending)' JSON or XLSX menu items. The TypeError fires regardless of file content (it even fires when the user cancels the file dialog, because _load_json returns {} and _set_concentrations({}) is still called with one arg).
- **Fix:** Pass replace_all through: in `_load_concentrations_json` call `self._set_concentrations(concentrations, replace_all)`, and in `_load_concentrations_xlsx` call `self._set_concentrations(concentrations, replace_all)`.

#### H23. _save_fluxes crashes with TypeError when user cancels the save dialog

- **Location:** `cnapy/gui_elements/main_window.py:3219-3241` · **Category:** crash · **Verifier confidence:** high
- **Problem:** `_get_filename` returns None (bare `return` at line 3214) when the user cancels the save dialog. `_save_fluxes` assigns `filename = self._get_filename(filetype)` (line 3220) and then unconditionally uses it: `open(filename, 'w', ...)` for CSV (line 3225) or `wb.save(filename)` for XLSX (line 3241). With filename=None this raises `TypeError: expected str, bytes or os.PathLike object, not NoneType` (open) or an openpyxl error (save). There is no guard for the cancel case. Note that the sibling method all_in_out_fluxes (line 2757) DOES guard against None, confirming the omission here is a bug.
- **Trigger:** Clipboard > 'Save flux solution as CSV...' or 'Save flux solution as Excel...' followed by pressing Cancel in the file-save dialog.
- **Fix:** After `filename = self._get_filename(filetype)` add `if filename is None: return`. Also annotate `_get_filename` return type as `str | None` since it can return None.

#### H24. Secretion template entries (positive value) become forced-flux lb=ub=1000, making model infeasible

- **Location:** `cnapy/gui_elements/media_management_dialog.py:856-861` · **Category:** math · **Verifier confidence:** high
- **Problem:** In _apply_template the template value is always treated as the lower bound: `for pattern, lb in all_components.items()`. The upper bound is computed as `ub = rxn.upper_bound if lb < 0 else 1000`. Several templates intentionally encode a POSITIVE value to represent allowed secretion/production (Heterotrophic Plant Culture `EX_co2_e: 1000.0` at line 260 commented "CO2 release", and Autotrophic `EX_o2_e: 1000.0` at line 278 commented "O2 production"). For these, lb=1000.0 and since lb>=0, ub=1000. The reaction is therefore constrained to (1000, 1000), i.e. forced to carry exactly +1000 flux, instead of being allowed to secrete up to 1000. This over-constrains the model and will typically render FBA infeasible.
- **Trigger:** User selects the 'Heterotrophic Plant Culture' or 'Autotrophic (Light + CO2)' template (or any custom media with a positive component value) and clicks 'Apply Media' on a model containing that exchange reaction; the resulting scenario forces lb=ub=1000 and FBA becomes infeasible.
- **Fix:** Interpret the value sign correctly: for a positive value, set it as the upper bound and keep lb at its existing model lower bound (or 0), e.g. `if lb < 0: ub = rxn.upper_bound; bounds=(lb, ub)` else `bounds=(rxn.lower_bound if rxn.lower_bound<0 else 0, lb)`. Alternatively document that positive values mean secretion upper bound and apply accordingly.

#### H25. Operator-precedence bug causes KeyError when filtering strain designs by must_occur reactions

- **Location:** `cnapy/gui_elements/mode_navigator.py:326-329` · **Category:** logic · **Verifier confidence:** high
- **Problem:** In the strain-design branch (mode_type==2), the condition `if selected and r not in s or numpy.any(numpy.isnan(s[r])) or numpy.all(s[r] == 0):` is parsed by Python as `(selected and (r not in s)) or numpy.any(numpy.isnan(s[r])) or numpy.all(s[r] == 0)` because `and` binds tighter than `or`. When a mode `s` (a dict returned by modes[i]) does NOT contain reaction `r`, the first conjunct `selected and (r not in s)` is True only if `selected` is True; but if `selected` is False (mode already deselected by a prior reaction in must_occur), the expression falls through to `numpy.any(numpy.isnan(s[r]))`, which indexes `s[r]` for a key that does not exist, raising KeyError. Even when no crash occurs, the intended logic (deselect modes where r is absent OR nan OR all-zero) is wrong: a mode is deselected merely because some OTHER reaction was nan/zero, since the nan/all-zero tests are not gated by `selected`.
- **Trigger:** Strain Design Navigation active; user types a must-occur reaction ID in the selector and presses Enter, when at least one navigated strain design does not contain that reaction AND was already deselected by an earlier must-occur reaction (so `selected` is False for it). The KeyError propagates; it is NOT caught (the except in apply_selection only catches ValueError/IndexError).
- **Fix:** Parenthesize and gate every test by `selected`, and guard the dict lookup: `if selected and (r not in s or numpy.any(numpy.isnan(s[r])) or numpy.all(s[r] == 0)):` — and only index `s[r]` after confirming `r in s`, e.g. `if selected and (r not in s or numpy.any(numpy.isnan(s[r])) or numpy.all(s[r] == 0)):` with a short-circuit `r not in s or ...` so `s[r]` is only evaluated when present.

#### H26. Operator-precedence bug causes KeyError and inverted logic when filtering strain designs by must_not_occur

- **Location:** `cnapy/gui_elements/mode_navigator.py:332-335` · **Category:** logic · **Verifier confidence:** high
- **Problem:** The must_not_occur strain-design condition `if selected and r in s and not numpy.any(numpy.isnan(s[r])) or numpy.all(s[r] == 0):` parses as `(selected and r in s and not numpy.any(numpy.isnan(s[r]))) or numpy.all(s[r] == 0)`. The trailing `or numpy.all(s[r] == 0)` is evaluated whenever the first group is False, including when `r not in s`, so it indexes `s[r]` for an absent key → KeyError. Logically it is also wrong: a mode where `r` is absent should be KEPT (it does not contain the forbidden reaction), but the dangling `or numpy.all(s[r]==0)` would deselect it if it didn't crash. The all-zero clause is meant to be part of the AND chain, not a standalone OR.
- **Trigger:** Strain Design Navigation active; user enters a must-not-occur token (prefixed with '!') for a reaction that is absent from at least one navigated strain design, then presses Enter. KeyError propagates uncaught (apply_selection only catches ValueError/IndexError).
- **Fix:** Group correctly and guard lookups: `if selected and r in s and not numpy.any(numpy.isnan(s[r])) and not numpy.all(s[r] == 0):` so the all-zero test is part of the AND chain and `s[r]` is only read when `r in s`.

#### H27. Spinbox range hard-clamped to +/-1000 silently corrupts the analyzed flux sweep

- **Location:** `cnapy/gui_elements/robustness_analysis_dialog.py:281,296,496-497` · **Category:** logic · **Verifier confidence:** high
- **Problem:** min_spin and max_spin are created with setRange(-1000, 1000). When the user selects 'Use scenario bounds as range', _update_bounds_from_scenario does self.min_spin.setValue(lb) / self.max_spin.setValue(ub). QDoubleSpinBox.setValue silently clamps any value outside [-1000, 1000] to the nearest bound. COBRApy / SBML models very commonly use bounds of +/-1e6, +/-1e3 is only the cobra default; many imported models use 1000000 or float('inf') for effectively-unbounded exchange/internal reactions. A reaction whose true bound is 1e6 (or inf) is clamped to 1000, so the sweep silently runs over the wrong range and the user is never told. The same clamp applies if the user types e.g. 5000 manually. The resulting xs = linspace(min, max, steps) and the entire robustness curve are computed over a wrong, truncated interval, producing scientifically incorrect results with no warning.
- **Trigger:** User analyzes a reaction whose model/scenario bound magnitude exceeds 1000 (e.g. an SBML model using 1e6 or inf bounds), either by checking 'Use scenario bounds' or by typing a larger value; the sweep range is silently truncated to +/-1000.
- **Fix:** Widen setRange to cover realistic COBRA bounds (e.g. -1e7..1e7 or based on appdata default bounds), and detect/replace non-finite bounds. When populating from scenario bounds, if the value exceeds the spinbox range, expand the range first or warn the user instead of silently clamping.

#### H28. Annotation reaction_id edits are never persisted (no signal wired for cellWidget)

- **Location:** `cnapy/gui_elements/scenario_tab.py:182, 289-299, 449-455` · **Category:** data-loss · **Verifier confidence:** high
- **Problem:** The Reaction_id column of the annotations table is a QComplReceivLineEdit cellWidget (set in new_annotation_row at line 451). The only handler that writes annotation data into scen_values.annotations is cell_content_changed_annotations, which is connected solely to self.annotations.cellChanged (line 182). QTableWidget.cellChanged fires only when a QTableWidgetItem's data changes (the Key/Value columns), NOT when the text of a cellWidget changes. Therefore typing/changing a reaction ID in the annotation row never triggers cell_content_changed_annotations and is never stored. The reaction_id is only captured opportunistically when the user later edits the Key or Value item, reading whatever is in the reaction_id field at that moment. If the user sets the reaction_id after editing Key/Value (or never touches Key/Value again), the annotation is saved with an empty/stale reaction_id, so the annotation is silently dropped at apply time (appdata.py line 549-552 skips annotations whose reaction_id is missing/unknown).
- **Trigger:** User adds a scenario reaction annotation, types a reaction ID into the first column, but does not (re)edit the Key or Value cell afterward. The reaction_id is lost; on save/apply the annotation is skipped.
- **Fix:** Connect the reaction_id QComplReceivLineEdit's editingFinished/textChanged signal (in new_annotation_row) to a handler that updates scen_values.annotations for that row, or store the row index on the widget and update on its own signal rather than relying on cellChanged.

#### H29. Direct clear_flux_values() bypasses undo history, desyncing undo/redo across all apply paths

- **Location:** `cnapy/gui_elements/scenario_templates_dialog.py:764, 837, 928, 1020` · **Category:** logic · **Verifier confidence:** high
- **Problem:** All four apply methods (_apply_template, _apply_knockout, _apply_bookmark, _apply_custom) clear the scenario via self.appdata.project.scen_values.clear_flux_values(), which only calls dict.clear() and does NOT append a 'clear' record to appdata.scenario_past. The subsequent value writes go through scen_values_set / scen_values_set_multiple, which DO append 'set' records. The undo/redo machinery in main_window (undo_scenario_edit -> recreate_scenario_from_history) rebuilds scen_values by replaying scenario_past from an empty dict and applying every recorded 'set' in order. Because the clear was never recorded, all previously recorded 'set' entries that the clear was meant to discard get replayed again, so the reconstructed scenario contains stale reaction bounds that should have been removed. AppData already exposes scen_values_clear() precisely to keep history consistent; this code bypasses it.
- **Trigger:** Apply any template/knockout/bookmark/custom scenario without merge so the dict is cleared, then press Undo and Redo (or any action that calls recreate_scenario_from_history). Reactions that were cleared reappear, producing an incorrect scenario / wrong FBA bounds.
- **Fix:** Route the clear through self.appdata.scen_values_clear() (which records the 'clear' in scenario_past) instead of calling scen_values.clear_flux_values() directly.

#### H30. load() crashes with KeyError on a setup saved for a different model

- **Location:** `cnapy/gui_elements/strain_design_dialog.py:1532-1561` · **Category:** crash · **Verifier confidence:** high
- **Problem:** After warning (but not aborting) when MODEL_ID does not match (lines 1467-1473), load() iterates the saved KOCOST/KICOST/GKOCOST/GKICOST and indexes self.reaction_itv[r] / self.gene_itv[g] directly. If a saved reaction or gene id is absent in the current model, the dict lookup raises KeyError and aborts loading mid-way, leaving the dialog in a partially-populated, inconsistent state.
- **Trigger:** Load a .sdc setup file that references reaction/gene ids not present in the currently open model (e.g., a setup built for another model, exactly the scenario the line-1467 warning anticipates).
- **Fix:** Guard each lookup, e.g. `if r in self.reaction_itv:` (and analogously for genes), skipping or collecting unknown ids and reporting them, instead of indexing unconditionally.

#### H31. self.module_edit() calls a dict as if it were a method -> TypeError

- **Location:** `cnapy/gui_elements/strain_design_dialog.py:1584, 1645` · **Category:** crash · **Verifier confidence:** high
- **Problem:** self.module_edit is a dict (initialized at line 230 as `self.module_edit = {}`), not a callable. The intent was clearly to call the method self.edit_module() (defined at line 811). Calling self.module_edit() raises `TypeError: 'dict' object is not callable`, crashing the Compute action.
- **Trigger:** Click 'Compute' when at least one module row was added via '+' but never set up (m is None) — line 1575 branch is entered and line 1584 runs. Or when more than one bilevel module (OptKnock/RobustKnock/OptCouple) is defined with a non-OPTLANG solver — line 1636 branch runs line 1645.
- **Fix:** Replace `self.module_edit()` with `self.edit_module()` at lines 1584 and 1645.

#### H32. Failed strain design emits a 3-tuple instead of SDSolutions, crashing conclude_computation

- **Location:** `cnapy/gui_elements/strain_design_dialog.py:1910-1911` · **Category:** crash · **Verifier confidence:** high
- **Problem:** On the exception path of SDComputationThread.run, line 1910 builds a proper SDSolutions object into `sd_solutions`, but line 1911 ignores it and emits `pickle.dumps(([], [], ERROR))` — a plain tuple. The connected slot SDComputationViewer.conclude_computation (line 1829, wired in main_window.py:1006) unpickles this and calls `self.solutions.get_num_sols()` (line 1832). A tuple has no get_num_sols method, so the GUI slot raises AttributeError every time a computation fails.
- **Trigger:** Any strain design computation that raises (e.g., infeasible problem, missing/invalid solver license — the project's known broken Gurobi license guarantees this) drives run() into the except block, emitting the bad tuple and crashing the result slot.
- **Fix:** Emit the already-constructed object: `self.finished_computation.emit(pickle.dumps(sd_solutions))` at line 1911.

#### H33. Empty/blank scenario constraint causes AttributeError on constraint[0].items()

- **Location:** `cnapy/gui_elements/thermodynamics_dialog.py:220-221` · **Category:** crash · **Verifier confidence:** high
- **Problem:** The constraint loop assumes constraint[0] is always a dict: `stoichiometries={key: value for key, value in constraint[0].items()}`. However Scenario defines `empty_constraint = (None, "", "")` (cnapy/appdata.py:258) and appends it into scen_values.constraints in the incompatible/placeholder path (cnapy/appdata.py:336). When such an empty constraint is present, constraint[0] is None and `None.items()` raises AttributeError. The try/except is commented out (lines 322-324), so the exception propagates, the dialog hangs, and the BusyCursor is never reset.
- **Trigger:** Load/run a scenario that contains an empty constraint placeholder (None, '', '') — e.g. a constraint row that was incompatible with the current model on scenario load — then press Compute.
- **Fix:** Skip empty constraints before processing: `if constraint[0] is None or constraint[1] == '': continue` (or check direction in {'>=','<=','='}) at the top of the loop body.

#### H34. OPTIMAL result box always reports "Maximum yield" even when minimizing

- **Location:** `cnapy/gui_elements/yield_optimization_dialog.py:143-160` · **Category:** logic · **Verifier confidence:** high
- **Problem:** In compute(), the optimization sense is read at lines 107-110 into `sense` ("Maximum" or "Minimum"). The unbounded and undefined branches correctly use `sense` (lines 127, 137). But the OPTIMAL branch hardcodes the literal string "Maximum yield (" (line 152) and never uses `sense`. When the user selects "minimize" in the sense_combo, yopt correctly computes a minimized yield, yet the result dialog mislabels it as the maximum. This is a silently incorrect scientific result presented to the user.
- **Trigger:** User selects "minimize" in the sense combo and the problem is OPTIMAL; the resulting message box claims "Maximum yield ..." while the displayed value is actually the minimum.
- **Fix:** Use the computed `sense` variable: replace the hardcoded "Maximum yield (" with `sense + " yield ("` so the label matches the chosen optimization direction.

#### H35. ROOM Big-M becomes infinite for reactions with unbounded (inf) bounds, crashing the solver

- **Location:** `cnapy/moma.py:186-187` · **Category:** math · **Verifier confidence:** high
- **Problem:** M_upper and M_lower are computed as max(abs(ub - w_upper), 1000) / max(abs(w_lower - lb), 1000). When a reaction's upper_bound is +inf or lower_bound is -inf (very common for exchange/demand/sink/transport reactions, where users routinely set unbounded bounds), abs(inf - w_upper) evaluates to inf, so M_upper/M_lower becomes inf. The resulting constraint expression rxn.flux_expression - inf*y (or + inf*y) contains an infinite coefficient, which the LP/MILP backend rejects or mis-scales. The reference implementation in test_data/Simulator.py:300-301 explicitly caps this with `ub - w_upper if ub < 1000 else 2000`, but moma.py never caps — it floors at 1000 and lets inf propagate upward.
- **Trigger:** Call room(model, reference_fluxes) where any reaction listed in reference_fluxes has upper_bound == +inf or lower_bound == -inf (e.g. an exchange reaction opened to unbounded uptake/secretion). The Big-M coefficient becomes inf and the solver errors out or returns garbage.
- **Fix:** Cap the Big-M at a finite value when bounds are large/infinite, mirroring the reference: e.g. `M_upper = (ub - w_upper) if (ub is not None and ub < 1000) else 2000` and `M_lower = (w_lower - lb) if (lb is not None and lb > -1000) else 2000`, ensuring M is finite and >= the needed relaxation (ub - w_upper) / (w_lower - lb).

---

## C. Medium-severity findings (70)

_Edge-case crashes, wrong results on non-default inputs, resource leaks, or data loss on uncommon paths._

#### M1. format_flux_value corrupts large values in scientific notation by stripping the exponent's trailing zero

- **Location:** `cnapy/appdata.py:149-150` · **Category:** data-loss · **Verifier confidence:** high
- **Problem:** format_flux_value does `str(round(float(flux_value), self.rounding)).rstrip("0").rstrip(".")`. For values whose str() uses scientific notation with an exponent ending in 0 (e.g. 1.5e20 -> '1.5e+20'), rstrip('0') removes the trailing '0' of the EXPONENT, producing '1.5e+2'. The displayed number is then 1.5e2 instead of 1.5e20 - an 18-orders-of-magnitude error silently shown to the user. The rstrip logic only safely strips trailing fractional zeros for normal decimal notation; it does not account for the exponent field.
- **Trigger:** A flux/bound value displayed via format_flux_value is >= ~1e16 and its string form ends in a zero (e.g. 1e16, 1e20, 1.5e20, 5e30). Such magnitudes can arise from very large or effectively-unbounded reaction bounds, enzyme-constraint capacities, or pathological FVA results.
- **Fix:** Avoid blind rstrip; only trim the fractional part when there is no exponent. E.g. detect 'e' in the string and skip stripping, or format explicitly: `s = format(round(float(flux_value), self.rounding), 'g')` (g-format already drops insignificant trailing zeros and handles exponents correctly), or strip only on the mantissa portion.

#### M2. Malformed numeric/color config values crash startup (ValueError/SyntaxError not caught)

- **Location:** `cnapy/application.py:220-291` · **Category:** exception · **Verifier confidence:** high
- **Problem:** In read_config, each config value is parsed via int()/float()/ast.literal_eval() inside a try/except that only catches (KeyError, NoOptionError). The conversions raise ValueError (int/float on bad text, ast.literal_eval on an invalid name) or SyntaxError (ast.literal_eval on malformed brackets). These are NOT among the caught exceptions, and the only outer handler is `except NoSectionError`. A single corrupted value in cnapy-config.txt (e.g. scen_color, font_size, box_width, rounding, abs_tol, or recent_cna_files) therefore propagates out of read_config and out of __init__, crashing the entire application at startup instead of falling back to defaults as the surrounding code clearly intends.
- **Trigger:** User's cnapy-config.txt has any numeric/color/recent_cna_files field with a non-parseable value (manual edit, partial write, version skew, or corruption). E.g. font_size=abc or scen_color=#FF0000 or recent_cna_files=[1,2 . CNApy then fails to launch with an uncaught traceback.
- **Fix:** Broaden each inner handler to also catch ValueError/SyntaxError (e.g. `except (KeyError, NoOptionError, ValueError, SyntaxError):`), or wrap the conversion separately and fall back to the existing appdata default while logging the bad value.

#### M3. Non-boolean use_results_cache value raises uncaught ValueError at startup

- **Location:** `cnapy/application.py:282-284` · **Category:** exception · **Verifier confidence:** high
- **Problem:** getboolean(..., fallback=...) only uses the fallback when the option/section is missing; when the option exists but holds a value configparser cannot interpret as boolean it raises ValueError. This call sits in the outer try whose only handler is `except NoSectionError`, so a malformed use_results_cache value propagates and crashes startup rather than using the intended fallback.
- **Trigger:** cnapy-config.txt contains use_results_cache set to any string that is not a recognized boolean literal (e.g. use_results_cache = maybe). Application crashes on launch.
- **Fix:** Wrap the getboolean call in its own try/except ValueError that falls back to self.appdata.use_results_cache, mirroring the per-option handling used for the other fields.

#### M4. model_optimization_with_exceptions returns None on swallowed exception, crashing net_conversion caller

- **Location:** `cnapy/core_gui.py:47-55` · **Category:** crash · **Verifier confidence:** high
- **Problem:** model_optimization_with_exceptions returns model.optimize() on success (a cobra Solution), but on any exception it falls into the except block which does NOT re-raise and has no return statement, so it implicitly returns None. One of its callers, MainWindow.net_conversion() (cnapy/gui_elements/main_window.py:2594-2595), uses the result directly as `solution = model_optimization_with_exceptions(model); if solution.status == "optimal":`. When optimization raises any exception, the function returns None and `None.status` raises AttributeError: 'NoneType' object has no attribute 'status', which propagates uncaught out of net_conversion (no surrounding try/except).
- **Trigger:** User invokes 'Net conversion of external metabolites' (MainWindow.net_conversion) while the active solver raises during model.optimize() — e.g. infeasible/unbounded model raising OptimizationError, a misconfigured/unavailable solver, or any non-community solver error. The function returns None and `solution.status` raises AttributeError.
- **Fix:** Make the function's failure contract explicit and consistent: either re-raise after showing the dialog (so callers' own try/except handle it) or return a sentinel that all callers check. Minimal fix: add `raise` at the end of the except block (after the optional dialog) so the caller's error handling triggers, and update net_conversion to guard against None / catch the exception. Also correct the `-> None` annotation to `-> Optional[cobra.Solution]`.

#### M5. Non-community-edition solver/model errors are silently swallowed with no user feedback

- **Location:** `cnapy/core_gui.py:50-54` · **Category:** logic · **Verifier confidence:** high
- **Problem:** The bare `except Exception` only surfaces an error to the user when has_community_error_substring(exstr) is True. For every other failure (infeasible/unbounded raising an error, solver crash, license problems not matching the hardcoded substrings, numerical errors), the exception is caught, NOT logged, NOT printed, NOT re-raised, and no message box is shown — the function just returns None. For the fba()/fba_optimize_reaction() callers (which route through process_fba_solution, where hasattr(None,'status') is False), this means the FBA action silently does nothing: the previous/stale solution and comp_values are left untouched and the user gets no indication the optimization failed. This violates fail-fast/never-suppress-silently and can present stale flux results as if current.
- **Trigger:** Run FBA (or any caller of this function) on a model/solver combination that raises a non-community error during optimize() — e.g. GLPK numerical failure, an unexpected optlang/solver exception, or a license error whose text doesn't contain the hardcoded substrings ('10012','10013','license 12', etc.). No error is shown; the GUI silently keeps the old solution.
- **Fix:** Add an else branch that surfaces the error: print/log exstr and call utils.show_unknown_error_box(exstr) (as pfba() already does at main_window.py:2560-2562), and/or re-raise so callers can react. Do not swallow unrecognized exceptions silently.

#### M6. updateReactionStoichiometry JS handler reuses loop variable 'i' for nested loops, corrupting the outer iteration over search records

- **Location:** `cnapy/data/escher_cnapy.html:251-282` · **Category:** logic · **Verifier confidence:** high
- **Problem:** The handler updateReactionStoichiometry declares the outer loop with 'for (i=0; i<records.length; i++)' (implicit global i, line 251) and then inside the matching branch declares an inner loop ALSO using 'for (i=0; i<reaction.metabolites.length; i++)' (line 255) over the same global i. After the inner loop completes, i equals reaction.metabolites.length, which the outer loop's i++ then continues from, skipping records and/or terminating the outer loop early. When the same reaction id appears on the map multiple times (Escher explicitly allows this, see comment at escher_map_view.py:168), only the first matching record gets correctly processed; subsequent duplicate occurrences are skipped, so their on-map stoichiometry is not updated.
- **Trigger:** A reaction that is drawn more than once on the same Escher map has its stoichiometry/reversibility changed; the second and later instances on the map are not updated because the shared loop counter 'i' was advanced by the inner metabolite loop.
- **Fix:** Declare the inner loop counter with a distinct local variable (e.g. 'for (var j = 0; j < reaction.metabolites.length; j++)') and declare the outer 'i' with let/var to avoid the shared global.

#### M7. Duplicate enzyme entries in one ec-rxn are silently summed, corrupting subunit copy counts

- **Location:** `cnapy/ecmodel/ec_structure.py:262-264` · **Category:** math · **Verifier confidence:** high
- **Problem:** from_yaml_omap iterates the enzyme pairs of each reaction and appends one COO triple per pair into rxn_enz_rows, then builds the matrix via _triples_to_csr -> coo_matrix(...).tocsr(). scipy's COO->CSR conversion SUMS values at duplicate (row, col) coordinates. If a single reaction's enzymes omap lists the same UniProt id more than once (a realistic malformed/relaxed input, which from_yaml_omap explicitly promises to 'tolerate' per its docstring and the missing-enzyme handling at lines 254-261), the two subunit-copy counts are silently added together instead of being deduplicated or treated as an error. The entry[i, j] value (documented as 'the subunit-copy count of enzyme j in reaction i') is then wrong, silently corrupting the EC stoichiometry used downstream in enzyme-constrained FBA.
- **Trigger:** Loading a GECKO YAML where an ec-rxns entry's 'enzymes' omap contains the same UniProt accession twice (duplicate key in a relaxed/hand-edited or partially-exported model). The two subunit counts are added rather than reported or overwritten.
- **Fix:** Detect duplicate (rxn_index, enz_index) pairs before building the matrix and either raise a ValueError, log a warning, or keep last/first deterministically. E.g. accumulate into a dict keyed by (i, j) with explicit overwrite-or-raise semantics instead of relying on COO's implicit summation.

#### M8. from_dict crashes (TypeError/ValueError) on None or empty numeric fields, unlike the tolerant string fields

- **Location:** `cnapy/ecmodel/ec_structure.py:308-316` · **Category:** crash · **Verifier confidence:** high
- **Problem:** from_dict is documented as the inverse of to_dict and is used to reload EcStructure from the .cna JSON archive. The string fields (source, notes, eccodes, genes, enzymes, sequence) are converted with the None-tolerant `_as_str` helper, but the numeric fields kcat, mw, and concs are converted with a raw `float(v)` inside a list comprehension. If any element is None (JSON null) or an empty string, `float(v)` raises TypeError (for None) or ValueError (for ''), aborting the entire model load. This is an inconsistency: the author clearly intended tolerant parsing for the string fields but left the numeric fields brittle. A null/empty kcat is plausible in a hand-edited .cna file or one produced by a partial/older export, and there is no try/except around the load.
- **Trigger:** Loading a .cna archive whose ec-structure has a null or empty-string entry in the kcat, mw, or concs arrays (e.g. hand-edited JSON, or a model where a concentration was left unmeasured and serialized as null).
- **Fix:** Use the existing `_as_float` helper instead of raw `float`: `kcat=[_as_float(v) for v in data.get('kcat', [])]` (and likewise for mw and concs), matching the tolerant behavior already applied to the string fields.

#### M9. revert_to_gem corrupts/removes native reactions whose ID legitimately ends in _REV

- **Location:** `cnapy/ecmodel/ecmodel_builder.py:844-853` · **Category:** data-loss · **Verifier confidence:** medium
- **Problem:** The _REV cleanup selects reactions purely by id suffix: rev_rxns = [r for r in gem.reactions if r.id.endswith('_REV')]. _split_reversible_reactions records the real reverse-split reactions in ec_data.split_rxn_map, but revert ignores that map and matches on the string suffix instead. Any reaction in the model whose identifier natively ends in '_REV' (reverse-direction reactions exist in real BiGG/Human-GEM naming) is misclassified as a split artifact: base_id = r.id[:-4]; if a reaction with that base id exists, its lower_bound is overwritten (base_rxn.lower_bound = -rev_rxn.upper_bound) and then the '_REV' reaction is permanently removed by gem.remove_reactions(rev_rxns, remove_orphans=True). This silently deletes a genuine model reaction and corrupts another's bounds on revert.
- **Trigger:** Build an ecModel from a GEM that contains a reaction whose ID ends in '_REV' (or any user reaction named *_REV), then call revert_to_gem. That reaction is removed and a sibling's lower bound is clobbered.
- **Fix:** Restrict the _REV cleanup to reactions actually recorded in ec_data.split_rxn_map (the rev ids stored as the second element of each split list) instead of matching on the '_REV' string suffix.

#### M10. apply_proteomics over-debits prot_pool when a usage reaction is constrained twice or already connected, corrupting pool budget

- **Location:** `cnapy/ecmodel/ecmodel_builder.py:972-992` · **Category:** math · **Verifier confidence:** high
- **Problem:** apply_proteomics subtracts every applied level from the pool capacity unconditionally (total_applied_mg += level for each uid in prot_data). It does not check whether the usage reaction was ALREADY disconnected from the pool by a prior apply_proteomics call. If apply_proteomics is invoked twice (e.g. user re-applies proteomics, or applies a second proteomics set without calling remove_proteomics first), the first run already reduced pool_rxn.lower_bound (via current_cap = usage_capacity(pool_rxn)). The second run reads the ALREADY-reduced current_cap and subtracts the levels again, double-debiting the pool. The pool capacity collapses toward 0 and the model becomes artificially over-constrained or infeasible, silently producing wrong FBA results.
- **Trigger:** Call apply_proteomics twice without an intervening remove_proteomics (e.g. user loads proteomics, then loads a corrected/second proteomics file). The pool capacity is reduced by the applied total on each call, compounding the reduction.
- **Fix:** Before debiting, reset the pool to its full bound (pool_rxn.lower_bound = -ec_data.pool_bound()) and only add an enzyme's level to total_applied_mg if its usage reaction was still connected to the pool (i.e. its prot_pool coefficient was non-zero) prior to this call, so re-application is idempotent.

#### M11. Duplicate metabolite/reaction ids silently collapse, losing entries instead of erroring

- **Location:** `cnapy/ecmodel/yaml_io.py:293-317, 70-71` · **Category:** data-loss · **Verifier confidence:** high
- **Problem:** Metabolites are collected into mets_by_id keyed by id (line 315) then bulk-added (line 317). If two metabolite entries share an id (corrupt/foreign YAML), the second silently overwrites the first in the dict, so the model ends up with fewer metabolites than the file declared, and reaction coefficients referencing the lost variant bind to the survivor. Likewise the !!omap constructor flattens any duplicated keys via out.update (ec_structure.py _to_dict / construct_mapping). No validation detects the count mismatch.
- **Trigger:** Load a YAML (e.g. a partially merged/foreign GECKO export) containing two metabolites or reactions with the same id; one is silently dropped and the model is built with wrong stoichiometry.
- **Fix:** Detect duplicate ids while collecting (if met.id in mets_by_id: raise EcYamlError) and verify len(mets_by_id) == len(document['metabolites']).

#### M12. Malformed numeric fields silently coerced to 0.0, masking corrupt bounds/kcat/MW

- **Location:** `cnapy/ecmodel/yaml_io.py:329-330, 340, 369-371` · **Category:** data-loss · **Verifier confidence:** high
- **Problem:** load_ecmodel routes lower_bound/upper_bound/metabolite coefficients/sigma/f/ptot through _as_float (ec_structure.py:380-386), which returns 0.0 for any value that is None, '' or not float-coercible, swallowing the error via except (TypeError, ValueError): return 0.0. A bound or stoichiometric coefficient that is malformed in the YAML (e.g. an accidental string, a list, 'NA') becomes 0.0 silently. A lower_bound silently becoming 0.0 turns a reversible reaction irreversible; a stoichiometric coefficient becoming 0.0 drops a metabolite from the reaction — both change the science with no error.
- **Trigger:** Load a YAML produced by a non-CNApy tool (or a hand-edited file) where a coefficient/bound is non-numeric; the reaction is silently built with a 0 coefficient/closed bound instead of failing.
- **Fix:** Distinguish 'absent' (use default) from 'present but unparseable' (raise EcYamlError with the reaction/metabolite id and offending value) instead of universally returning 0.0.

#### M13. Sampling fallback silently returns far fewer samples than requested without notifying the caller

- **Location:** `cnapy/flux_sampling.py:152-159` · **Category:** data-loss · **Verifier confidence:** high
- **Problem:** If the first sample() call raises any Exception, the code retries with min(n, 100) samples and processes=1. When the retry succeeds, the function returns that reduced result with no flag, warning, or error. A user who requested n=5000 (the dialog default) can receive only 100 samples; the GUI then reports 'Total samples: 100' as if it were a complete success. Statistics/correlations computed downstream are based on a silently truncated sample set, which is scientifically misleading for a sampling tool.
- **Trigger:** n > 100 and the first sample() call fails (e.g. transient solver error, multiprocessing failure with processes>1, OptGP convergence issue) while the reduced retry succeeds.
- **Fix:** Either propagate the original failure, or surface the degradation explicitly (return a status/metadata flag or log a warning) so the caller and user know fewer than n samples were produced. Do not silently substitute min(n,100).

#### M14. is_integer_vector_rounded crashes when fv_mat is a sparse/list-of-lists matrix

- **Location:** `cnapy/flux_vector_container.py:45-48` · **Category:** exception · **Verifier confidence:** high
- **Problem:** is_integer_vector_rounded iterates `for val in self.fv_mat[idx, :]`. When fv_mat is a scipy.sparse matrix (the documented storage mode after loading a sparse-saved .npz at lines 23-24, where tolist() yields a scipy.sparse lil_matrix), `self.fv_mat[idx, :]` is a 1xN sparse matrix and iterating it yields the whole 1-row matrix object rather than scalar values, so round(val, decimals) / .is_integer() fails. It is called from central_widget.update_mode for mode_type==0 (EFM display). A user who saves EFMs in sparse form and reloads them hits this code path and crashes. The TODO comment acknowledges the limitation, but the method still executes and raises in production.
- **Trigger:** EFMs whose fv_mat is stored as a scipy.sparse matrix (e.g. loaded from a .npz saved with a sparse fv_mat, which the object-dtype branch at lines 23-24 restores as scipy.sparse) are displayed in EFM mode (mode_type==0); update_mode calls is_integer_vector_rounded -> AttributeError/TypeError.
- **Fix:** Handle the sparse case explicitly: densify the single row first, e.g. `row = self.fv_mat[idx, :]; row = row.toarray().ravel() if scipy.sparse.issparse(self.fv_mat) else numpy.asarray(self.fv_mat[idx, :]).ravel()` and then test integrality on the dense 1D array.

#### M15. Header read crashes with IndexError on empty/truncated efmtool binary file

- **Location:** `cnapy/flux_vector_container.py:80-82` · **Category:** crash · **Verifier confidence:** high
- **Problem:** The memmap header is read with numpy.fromfile(...)[0]. If the file is empty or shorter than the expected header (8-byte >i8 efm count + 4-byte >i4 reac count), numpy.fromfile returns a zero-length array and indexing [0] raises IndexError. There is no length/validity check, so a truncated or zero-byte efms.bin (e.g. efmtool aborted mid-write, or disk full) produces an uncaught IndexError instead of a meaningful error.
- **Trigger:** FluxVectorMemmap is constructed on an efms.bin that is empty or truncated (efmtool process killed/aborted, partial write, 0 EFMs written without a full header). The [0] index on the empty fromfile result raises IndexError.
- **Fix:** Read the header with explicit size validation: check that each fromfile result has the expected length (e.g. assert len(arr) == 1) and raise a clear, caught exception (or return an empty container) when the file is too short.

#### M16. FluxVectorMemmap.__del__ raises AttributeError when __init__ fails before fv_mat is set

- **Location:** `cnapy/flux_vector_container.py:95-96` · **Category:** exception · **Verifier confidence:** high
- **Problem:** FluxVectorMemmap.__del__ unconditionally does `del self.fv_mat`. If __init__ raises before super().__init__() assigns self.fv_mat (e.g. the open() at line 80 fails, the fromfile header reads at lines 81-82 raise IndexError on a truncated file, or numpy.memmap construction at line 84 fails), the partially-constructed instance still gets finalized, and __del__ raises AttributeError: fv_mat. This is reported as 'Exception ignored in __del__' and masks/obscures the real initialization error during garbage collection.
- **Trigger:** Constructing FluxVectorMemmap with a missing/locked/truncated efms.bin (open fails, or fromfile([0]) raises IndexError on an empty/short file) so __init__ raises before super().__init__ runs; the subsequent GC of the partial object triggers __del__ -> AttributeError.
- **Fix:** Guard the cleanup: `if hasattr(self, 'fv_mat'): del self.fv_mat` in both __del__ and clear(); or initialize self.fv_mat = None at the very start of __init__ before the file operations.

#### M17. open_in_browser uses unbound variable after AttributeError is swallowed

- **Location:** `cnapy/gui_elements/annotation_widget.py:101-109` · **Category:** crash · **Verifier confidence:** high
- **Problem:** In open_in_browser the assignments to identifier_type/identifier_value are inside a try/except AttributeError: pass. If either table cell item is None (empty key or empty value cell), .text() raises AttributeError which is swallowed, but execution then falls through to `if identifier_value.startswith("[")` which references the never-assigned name -> NameError, crashing the handler. The except clause should abort, not pass.
- **Trigger:** Select an annotation row whose value (column 1) or key (column 0) cell has no QTableWidgetItem (e.g. a row created with '+' where only the key was typed, or a key with empty value) and click 'Open chosen in browser'.
- **Fix:** Replace `except AttributeError: pass` with `except AttributeError: return` (and optionally show an information message), so the function does not proceed with unbound variables.

#### M18. ast.literal_eval on annotation value can raise uncaught ValueError/SyntaxError/IndexError

- **Location:** `cnapy/gui_elements/annotation_widget.py:106-107` · **Category:** crash · **Verifier confidence:** high
- **Problem:** open_in_browser does `if identifier_value.startswith("["): identifier_value = ast.literal_eval(identifier_value)[0]` with no try/except. If the value text starts with '[' but is not a valid Python list literal (e.g. '[foo' or '[a,b'), ast.literal_eval raises ValueError/SyntaxError. If the value is the literal '[]' (empty list), `ast.literal_eval('[]')[0]` raises IndexError. Any of these crashes the handler. Note apply_annotation (line 76-80) guards the same parse with try/except, but this path does not.
- **Trigger:** Select an annotation row whose value cell contains text starting with '[' that is not a valid non-empty list literal (e.g. '[]' or '[abc'), then click 'Open chosen in browser'.
- **Fix:** Wrap the parse in try/except (ValueError, SyntaxError, IndexError) and fall back to the raw string or abort with a message, mirroring the guard already used in apply_annotation.

#### M19. Per-KO exception isolation swallows all errors and reports KO biomass = 0 (false essentiality)

- **Location:** `cnapy/gui_elements/batch_moma_room_dialog.py:289-290` · **Category:** exception · **Verifier confidence:** medium
- **Problem:** _analyze_gene_knockout and _analyze_reaction_knockout wrap the whole solve in a bare `except Exception:` that returns (0.0, 0.0/None). This conflates a genuine infeasible/zero-growth knockout with any unexpected failure (solver error, KeyError because biomass_reaction not in solution.fluxes for that sub-problem, ROOM big-M failure, transient solver exception). Every such failure is silently recorded as ko_biomass = 0.0, which then yields growth_ratio 0.0 and is_essential = True. The user gets a results table that reports those targets as essential with zero growth, with no indication that the computation actually failed — silently incorrect scientific output.
- **Trigger:** Any knockout where the MOMA/ROOM solve raises (e.g., solver hiccup, or solution.fluxes lookup error). That target is reported as essential (KO growth 0) rather than flagged as failed.
- **Fix:** Catch narrowly and record a failure marker (e.g., return None/NaN or add result['status']='failed') so failed solves are visually distinguished from true zero-growth knockouts, instead of mapping every error to essential.

#### M20. convert_current_escher_to_cnapy busy-waits forever if the Escher callback never fires

- **Location:** `cnapy/gui_elements/central_widget.py:1013-1019` · **Category:** concurrency · **Verifier confidence:** high
- **Problem:** After requesting map data from the Escher (QWebEngine) view, the code spins `while sem[0] < 1: QApplication.processEvents()` with no timeout. retrieve_map_data relies on an asynchronous JS->Python bridge callback to set sem[0]=1. If the web page errored, the bridge is not connected, or the JS callback never resolves (load failure, JS exception), this loop never terminates and the GUI hangs in a 100% CPU busy-spin, requiring force-kill. There is no timeout, no iteration cap, and no error path.
- **Trigger:** Converting an Escher map whose JS bridge callback does not fire (page load failure, JS error, or bridge not yet ready) — the while loop spins indefinitely, hanging the application.
- **Fix:** Bound the wait with a deadline (e.g. time.monotonic()+timeout) and break with a QMessageBox warning if exceeded; or make retrieval fully asynchronous with a completion callback instead of a synchronous processEvents spin.

#### M21. Clipboard calculator raises KeyError when clipboard/current value sets have mismatched keys

- **Location:** `cnapy/gui_elements/clipboard_calculator.py:115-123` · **Category:** crash · **Verifier confidence:** high
- **Problem:** The compute loop iterates `for key in self.appdata.project.comp_values` and then indexes `l_comp[key]` (line 119) and `r_comp[key]` (line 123) with raw subscription. When the left or right operand is 'Clipboard values', l_comp/r_comp is self.appdata.clipboard_comp_values, which is a snapshot taken at an earlier time and is NOT guaranteed to contain the same reaction ids as the current comp_values (model edited, reactions added, or different scenario). Any key present in comp_values but absent from the clipboard dict raises KeyError and crashes compute().
- **Trigger:** Store values to clipboard, then add/remove reactions or switch to a model/scenario whose comp_values keys differ from the clipboard's, then open the calculator with one side = 'Clipboard values' and compute. KeyError on the first missing reaction id.
- **Fix:** Guard each access, e.g. skip keys missing from the operand dict or default to (0,0): `lv = l_comp.get(key)`; if None, `continue`. Optionally warn the user that key sets differ.

#### M22. Division operator can raise ZeroDivisionError on user/clipboard values

- **Location:** `cnapy/gui_elements/clipboard_calculator.py:143-144` · **Category:** math · **Verifier confidence:** high
- **Problem:** combine() performs `llb / rlb` and `lub / rub` for the '/' operator with no zero check. Reaction flux bounds frequently include 0 (e.g., a blocked reaction lb=ub=0, or an irreversible reaction lb=0). Dividing by such values raises ZeroDivisionError (for int/int) or yields inf/nan, propagating an exception out of compute() (which is connected directly to the Compute button) and crashing the operation.
- **Trigger:** Select '/' operator with the right operand being a reaction whose lb or ub is 0 (very common), or type 0 into the right value field, then Compute. ZeroDivisionError propagates and aborts; for the manual-entry int path it crashes the dialog.
- **Fix:** Catch ZeroDivisionError in compute()/combine() and report via QMessageBox, or define a policy (skip key, produce inf/nan deliberately) and apply it consistently.

#### M23. Cancelling the results-cache directory chooser silently sets cache dir to the current working directory

- **Location:** `cnapy/gui_elements/config_dialog.py:196-199` · **Category:** logic · **Verifier confidence:** high
- **Problem:** QFileDialog.getExistingDirectory() returns an empty string when the user cancels. `Path("")` evaluates to `Path('.')`, whose `.exists()` is always True (the process CWD). So when the user cancels the dialog, the guard at line 197 passes and line 199 sets the results cache directory button text to the current working directory instead of leaving it unchanged. This contrasts with choose_work_directory (line 190) which correctly guards `if not directory or len(directory) == 0`. Result: FVA/result caches get written to an unintended location (CWD).
- **Trigger:** User clicks the results-cache directory button then presses Cancel in the file dialog.
- **Fix:** Capture the raw string first and bail on empty: `selected = dialog.getExistingDirectory(); if not selected: return; directory = Path(selected)`.

#### M24. Font size field has no validator; empty/non-numeric input crashes apply() with ValueError

- **Location:** `cnapy/gui_elements/config_dialog.py:252` · **Category:** exception · **Verifier confidence:** high
- **Problem:** The font_size QLineEdit (line 39) is created without any validator, unlike rounding (QIntValidator, line 124) and abs_tol (QDoubleValidator, line 135). In apply(), line 252 does `new_fontsize = float(self.font_size.text())` with no try/except. If the user clears the field or types non-numeric text, float() raises ValueError, which propagates out of the clicked-signal slot and is logged as an unhandled exception; the Apply button silently fails and no settings are saved.
- **Trigger:** User opens Config dialog, deletes the contents of the Font size field (or types e.g. 'big'), then clicks 'Apply Changes'.
- **Fix:** Attach a QDoubleValidator/QIntValidator to self.font_size, and/or wrap the parse: `try: new_fontsize = float(self.font_size.text())` and on ValueError show a QMessageBox and abort apply().

#### M25. Box width field has no validator; empty/non-integer input crashes apply() with ValueError

- **Location:** `cnapy/gui_elements/config_dialog.py:262` · **Category:** exception · **Verifier confidence:** high
- **Problem:** The box_width QLineEdit (line 49) has no validator. apply() calls `int(self.box_width.text())` at line 262 (and again at line 282) with no try/except. Empty or non-integer text raises ValueError, which aborts apply() before any settings (including colors, work_directory, etc.) are saved, so a user editing colors loses those edits too because the exception fires mid-apply.
- **Trigger:** User clears the Box width field or enters non-integer text and clicks 'Apply Changes'.
- **Fix:** Add a QIntValidator to self.box_width and/or guard the int() conversion with try/except, validating before mutating appdata.

#### M26. Connection script runs the bare interpreter name, resolving via PATH instead of sys.executable

- **Location:** `cnapy/gui_elements/configuration_cplex.py:130-138` · **Category:** logic · **Verifier confidence:** high
- **Problem:** The command is built as 'cd "{python_dir}" && {python_exe_name} "{...}setup.py" install' where python_exe_name is just the basename of sys.executable (e.g. 'python3' or 'python'). It is then run with shell=True. The cd into python_dir does NOT put that directory on the executable search path on POSIX shells (cwd is not on PATH), so the shell resolves the bare name 'python3'/'python' via the existing PATH — which may be a DIFFERENT interpreter than CNApy's sys.executable, or may not exist at all (CalledProcessError -> 'Run Error'). On a venv whose bin dir is not on PATH, the install fails entirely; on a system with multiple Pythons it can install the CPLEX bindings into the wrong interpreter, silently leaving CNApy without CPLEX while reporting nothing useful.
- **Trigger:** On Linux/macOS where CNApy's interpreter directory is not on the spawned shell's PATH (typical for an un-activated venv), or where PATH's 'python'/'python3' points to a different interpreter than sys.executable.
- **Fix:** Invoke the full interpreter path directly and avoid the cd hack: subprocess.check_call([sys.executable, os.path.join(self.cplex_directory.text(), 'python', 'setup.py'), 'install']) with shell=False, or at minimum use the quoted python_exe_path instead of python_exe_name.

#### M27. Hardcoded Python '3.10' path causes uncaught FileNotFoundError after successful install

- **Location:** `cnapy/gui_elements/configuration_cplex.py:160-162` · **Category:** crash · **Verifier confidence:** high
- **Problem:** get_and_set_environmental_variable() builds base_path = self.cplex_directory.text() + "cplex/python/3.10/" with the Python minor version '3.10' hardcoded, then immediately calls os.listdir(base_path). CNApy frequently runs on Python 3.9/3.11/3.12, in which case the CPLEX layout is <cplex>/cplex/python/3.11/ (or 3.9, 3.12), so the '3.10' directory does not exist. os.listdir then raises FileNotFoundError. This call happens at line 157 in the else-branch, which is OUTSIDE the try/except (the except at line 139 only covers subprocess.check_call), so the exception propagates uncaught and crashes the handler right after the user was shown a 'Success' message box. Even on Python 3.10 the same crash occurs if the installed CPLEX uses a different sub-folder layout.
- **Trigger:** User runs CNApy on any Python version other than 3.10 (e.g. 3.11/3.12/3.9), completes step 5 successfully, so get_and_set_environmental_variable() runs and os.listdir('<cplex>/cplex/python/3.10/') hits a non-existent directory.
- **Fix:** Derive the version folder from the running interpreter, e.g. ver = f"{sys.version_info.major}.{sys.version_info.minor}"; base_path = self.cplex_directory.text() + f"cplex/python/{ver}/". Also guard with if not os.path.isdir(base_path): show an informative message instead of calling os.listdir, and wrap the call in try/except FileNotFoundError.

#### M28. Network/zip download runs on the GUI thread with unguarded I/O, hanging or crashing the UI on failure

- **Location:** `cnapy/gui_elements/download_dialog.py:54-79` · **Category:** exception · **Verifier confidence:** high
- **Problem:** download() performs urllib.request.urlretrieve(), ZipFile.extractall(), os.mkdir/os.remove synchronously on the GUI thread with no try/except. The dialog text explicitly promises graceful handling ('If a working directory error occurs, you can solve by setting a working directory'), but none of the failure modes are caught: no internet / DNS failure / HTTP 404 raises URLError/HTTPError; a corrupt or partial download raises BadZipFile from ZipFile; an unwritable work_directory raises OSError from os.mkdir. Any of these propagates as an unhandled exception out of the slot, and because the call is on the GUI thread, the whole window is frozen for the duration of the (potentially large) download.
- **Trigger:** User clicks 'Yes, download...' with no internet connection, a 404 on the release asset, a corrupted/partial zip, or a non-writable working directory.
- **Fix:** Wrap the download/extract in try/except (URLError, HTTPError, OSError, BadZipFile) and show a QMessageBox.warning on failure; move the blocking work to a QThread/worker so the GUI stays responsive, and only show the success message after the worker completes.

#### M29. Partial/corrupt zip leaves a stale file that blocks all future download attempts

- **Location:** `cnapy/gui_elements/download_dialog.py:64-77` · **Category:** data-loss · **Verifier confidence:** high
- **Problem:** The download is gated by 'if not os.path.exists(target)'. urlretrieve writes directly to 'target'. If the download is interrupted (network drop, app killed) after the file is partially written but before os.remove(zip_path) at line 77 runs, the partial zip remains on disk. On the next attempt, 'if not os.path.exists(target)' is False, so the code SKIPS the download AND the extraction entirely, silently producing no projects and then showing the 'downloaded successfully' message box. The user can never recover via the dialog because the stale file always short-circuits the path; they must manually delete the file.
- **Trigger:** A download interrupted after the file begins/finishes writing but before extraction+removal complete (network drop or app exit), followed by re-opening the dialog and clicking download again.
- **Fix:** Download to a temporary path and atomically rename only after successful extraction, or always remove a pre-existing target before downloading; do not show the success message unless extraction actually occurred.

#### M30. Flux interpolation uses non-monotonic time points, silently producing wrong flux output

- **Location:** `cnapy/gui_elements/dynamic_fba_dialog.py:321-329` · **Category:** math · **Verifier confidence:** high
- **Problem:** flux_time_points is built by appending t every time the ODE RHS dfba_rhs is called (line 250). Adaptive solvers (RK45/Radau/BDF/LSODA, all selectable) evaluate the RHS at intermediate stages, rejected steps, and generally out of chronological order, so flux_time_arr is NOT monotonically increasing. np.interp(time, flux_time_arr, flux_arr) requires xp to be monotonically increasing; when it is not, numpy does NOT raise but returns silently incorrect values. Therefore every reaction flux reported in DFBAResult.fluxes (used for analysis/CSV) is wrong, not just imprecise.
- **Trigger:** Any successful dFBA run (the default RK45 is adaptive). The flux arrays in the Data tab / exported CSV / DFBAResult.fluxes are silently wrong.
- **Fix:** Sort the (time, value) pairs by time before interpolating, and de-duplicate repeated t values, e.g. `order = np.argsort(flux_time_arr); ftp = flux_time_arr[order]; fa = flux_arr[order]; uniq = np.concatenate(([True], np.diff(ftp) > 0)); fluxes[rid] = np.interp(time, ftp[uniq], fa[uniq])`. Better: capture fluxes only at accepted t_eval points via dense output or a post-solve pass.

#### M31. QThread.terminate() with no wait(); stale finished slot mutates GUI after reset

- **Location:** `cnapy/gui_elements/dynamic_fba_dialog.py:798-803` · **Category:** concurrency · **Verifier confidence:** high
- **Problem:** _stop_simulation calls self.simulation_thread.terminate() and immediately _reset_ui(), without disconnecting the finished signal or calling wait(). terminate() forcibly kills the thread at an arbitrary point (possibly mid cobra optimize / mid model.copy), which can corrupt the solver/model state and leak resources. If the simulation actually completed in the small race window, the custom finished signal is still connected and _on_simulation_finished will fire afterwards, re-enabling buttons and popping a 'Simulation Complete' QMessageBox for a run the user asked to stop. There is also no guard preventing starting a new thread while an old one is being torn down, and the old thread object is overwritten in _run_simulation, potentially dropping a still-running QThread (Qt warns/crashes if a QThread is destroyed while running).
- **Trigger:** User clicks Stop during a simulation (especially a long one); or clicks Run again quickly. Risk of spurious completion dialog, GUI state corruption, or Qt 'QThread destroyed while still running' crash.
- **Fix:** Use a cooperative cancellation flag checked inside run_dfba instead of terminate(); on stop, disconnect the finished signal, call requestInterruption()/wait(), and guard against launching a new thread while one is active.

#### M32. O(R*M) membership test rebuilds reaction-id list for every reaction on apply

- **Location:** `cnapy/gui_elements/flux_data_dialog.py:515` · **Category:** logic · **Verifier confidence:** high
- **Problem:** In _apply, the membership check `if rxn_id in [r.id for r in model.reactions]` rebuilds a full Python list of all model reaction ids on every iteration of the flux_data loop, then does a linear scan. For a genome-scale model (thousands of reactions) and a file with thousands of entries this is O(R*M) and can freeze the GUI thread for many seconds. The same pattern repeats at line 543 in the log2FC branch. _preview already does the right thing using a set comprehension once (lines 450, 478).
- **Trigger:** Apply flux data for a genome-scale model (e.g. ~2500+ reactions) with a large flux file; the Apply button blocks the GUI thread for a noticeable hang.
- **Fix:** Build the reaction-id set once before the loop: `model_rxn_ids = {r.id for r in model.reactions}` and test `if rxn_id in model_rxn_ids`.

#### M33. IndexError on Ctrl+C when no table cell is selected

- **Location:** `cnapy/gui_elements/flux_feasibility_dialog.py:515` · **Category:** crash · **Verifier confidence:** high
- **Problem:** copy_table_selection indexes selectedRanges()[0] without checking for an empty result. QTableWidget.selectedRanges() returns an empty list when nothing is selected, so `[0]` raises IndexError. The slot is wired to a Ctrl+C QAction on the table (lines 138-141), so pressing Ctrl+C while the table has focus but no active selection crashes.
- **Trigger:** Give the bm_constituents table focus with no cell selected (e.g. immediately after the table is (re)populated by update_bm_constituents_table which resets rows) and press Ctrl+C.
- **Fix:** Guard the slot: `ranges = self.bm_constituents.selectedRanges(); if not ranges: return; selection_range = ranges[0]`.

#### M34. Scenario-loading failures are silently swallowed, producing wrong results from an unconstrained model

- **Location:** `cnapy/gui_elements/flux_response_dialog.py:137-139, 559-562` · **Category:** exception · **Verifier confidence:** medium
- **Problem:** Both the worker (run) and _autodetect_range wrap load_scenario_into_model in 'except Exception: pass'. If applying scenario constraints raises (e.g. a scenario reaction/constraint references an id not present after model.copy(), or constraint construction fails), the exception is discarded and analysis proceeds on the model WITHOUT the scenario bounds the user explicitly requested via 'Use scenario constraints'. The FVA range and the whole flux-response scan are then computed against the wrong feasible space and silently reported as valid results, with no indication the scenario was not applied. This violates fail-explicitly and yields silently incorrect scientific output.
- **Trigger:** Enable 'Use scenario constraints' with a scenario whose application raises (malformed/added-reaction scenario, constraint referencing missing reaction). Results are computed and shown as if scenario applied, but it was not.
- **Fix:** Do not swallow; emit error_occurred (worker) / show QMessageBox (autodetect) with the failure, or at least set a visible warning flag in the results so the user knows constraints were not applied. Catch only the specific expected exception.

#### M35. Empty biomass reaction passes validation but crashes the worker via get_by_id('')

- **Location:** `cnapy/gui_elements/flux_response_dialog.py:594, 601-604, 142` · **Category:** crash · **Verifier confidence:** high
- **Problem:** _run_analysis only guards 'if not target_id or not product_id' (line 594); biomass_id may be empty. The validation loop (lines 601-604) uses 'if rxn_id and rxn_id not in model.reactions', so an empty biomass_id is skipped (the 'rxn_id and ...' short-circuits) and passes validation. The worker then unconditionally calls model.reactions.get_by_id(self.biomass_reaction) at line 142, which raises KeyError('') for an empty string (verified: cobra get_by_id('') -> KeyError). The whole analysis aborts with a confusing 'Flux response analysis failed' traceback even though the user only failed to pick a biomass reaction.
- **Trigger:** Open the dialog on a model where the biomass selector ends up empty (e.g. no objective/biomass-pattern reaction found so the combo has no auto-selected biomass, and the user clears it), select valid target+product, click Run. Worker crashes at line 142.
- **Fix:** Validate biomass_id explicitly before starting (require it, or skip the WT-biomass step when empty). At minimum check 'if not biomass_id' and warn, or guard the get_by_id call in the worker and emit a clear error.

#### M36. Predicted sampling treats FVA/bounds intervals as point reference fluxes, yielding wrong (often zero) reference and silently corrupting the sampling space

- **Location:** `cnapy/gui_elements/flux_sampling_dialog.py:80-85, 218-224` · **Category:** math · **Verifier confidence:** high
- **Problem:** The dialog enables 'Predicted Flux-Based Sampling' whenever appdata.project.comp_values is non-empty (line 80), without checking appdata.project.comp_values_type. comp_values_type==1 means the stored tuples are genuine flux INTERVALS from an FVA / show-model-bounds result (main_window.py:2688-2689 stores (reaction.lower_bound, reaction.upper_bound); fva() at main_window.py:2746 does the same), not a single predicted flux per reaction. get_reference_fluxes() (lines 218-224) then computes reference[rid] = (lb + ub) / 2 for every entry, treating that interval midpoint as a 'reference flux'. For a symmetric/unbounded reaction the interval is e.g. (-1000, 1000), so the midpoint is 0 even though the reaction may carry substantial flux. That bogus 0 reference is passed to perform_predicted_flux_sampling (main_window.py:2443-2457), which in 'bounds' constraint mode treats abs(flux)<1e-6 as near-zero and clamps the reaction to [-0.1, 0.1] (flux_sampling.py:121-124). The result is a silently incorrect, massively over-constrained sampling space with no error raised. The mode tooltip explicitly claims the reference comes 'from FBA, MOMA, etc.' (point solutions), confirming the midpoint-of-interval handling is unintended for FVA-type computed values.
- **Trigger:** User runs FVA or 'Show model bounds' (sets comp_values_type=1 with flux intervals), then opens Flux Sampling, selects 'Predicted Flux-Based Sampling' with constraint mode 'Bounds', and starts sampling. Reactions whose interval midpoint is ~0 (e.g. symmetric reversible reactions) get clamped to [-0.1,0.1], producing wrong sampling results with no warning.
- **Fix:** Only enable predicted mode for point-flux results: gate on `len(comp_values) > 0 and self.appdata.project.comp_values_type == 0` at line 80 (and disable with an explanatory tooltip otherwise). Alternatively, document/handle FVA intervals explicitly rather than silently using the midpoint as a reference flux.

#### M37. Biomass constrained as exact equality at fraction of max can over-constrain and bias FVA

- **Location:** `cnapy/gui_elements/fvseof_dialog.py:223,229` · **Category:** math · **Verifier confidence:** high
- **Problem:** The algorithm docstring says 'Constrain biomass to 95% (configurable) of max' (typically a lower bound: biomass >= fraction*max). The code instead fixes biomass to an exact equality at the fraction: biomass_rxn.bounds = (constrained_biomass, constrained_biomass) where constrained_biomass = max_biomass * biomass_fraction. Fixing biomass to exactly 95% of max (rather than >=95%) is a stronger, different constraint than the standard FSEOF/FVSEOF convention and changes which flux distributions are feasible, yielding different (and arguably incorrect per the stated method) FVA ranges and correlations.
- **Trigger:** Every feasible scan step: FVA is run with biomass pinned exactly at fraction*max instead of allowed to range from fraction*max up to max, narrowing feasible flux space and altering reported min/max ranges and the resulting regression-based candidate lists.
- **Fix:** Use a lower bound on biomass: biomass_rxn.lower_bound = constrained_biomass (leaving the upper bound at max_biomass or its original value), matching the documented and conventional FSEOF/FVSEOF constraint.

#### M38. FVA min/max failure silently substituted with 0.0, corrupting flux ranges and regressions

- **Location:** `cnapy/gui_elements/fvseof_dialog.py:249,251,259,261,266-267` · **Category:** math · **Verifier confidence:** medium
- **Problem:** When an individual reaction's min or max optimization is non-optimal (or raises), min_flux/max_flux is set to 0.0 instead of being marked missing. 0.0 is a valid, meaningful flux, so a solver failure is indistinguishable from a true zero. This poisons avg_flux=(min+max)/2 and sol_range=abs(max-min), and these fabricated zeros then enter the linear regressions (lines 301-304 read these values, with 'or 0.0' further masking). The result is silently incorrect correlation coefficients and candidate classification — the core scientific output of the tool — with no warning to the user.
- **Trigger:** Any scan step where biomass is fixed at constrained level and a particular reaction's directional optimization returns non-optimal (e.g., numerical infeasibility for that single-reaction objective under the equality biomass constraint). That reaction gets a spurious 0 flux for that step, skewing its regression slope/R and possibly mis-listing it as (or excluding it from) a candidate.
- **Fix:** Treat non-optimal sub-optimizations as missing (None) and exclude those (reaction, step) points from the per-reaction regression arrays, rather than substituting 0.0; or abort/record the step as unreliable and surface a warning.

#### M39. _run_whatif text-fallback uses substring match and silently selects the wrong enzyme

- **Location:** `cnapy/gui_elements/gecko_unified_dialog.py:1176-1187` · **Category:** logic · **Verifier confidence:** high
- **Problem:** When currentData() returns no uid, the fallback loop matches the typed text against combo items with 'if typed == item_uid or typed in item_text'. The 'typed in item_text' branch is a substring containment test, so typing a short or partial string (e.g. 'P12') matches the FIRST item whose label contains that substring and applies the kcat multiplier to that enzyme. The label format is f"{uid}  ({gene})", so a typed gene fragment or partial UID can match an unintended enzyme. The user receives no warning that a fuzzy/partial match was used; apply_kcat_multiplier then scales the wrong enzyme's kcat, producing a silently incorrect what-if FBA result.
- **Trigger:** User types a partial UniProt ID or gene fragment (rather than selecting from the dropdown) that is a substring of multiple enzyme labels; the first matching enzyme — not necessarily the intended one — is silently used for the kcat multiplier.
- **Fix:** Require an exact match (typed == item_uid or typed == item_text), or collect all candidate matches and refuse to proceed / prompt the user when the match is ambiguous, instead of using substring containment with first-match-wins.

#### M40. Cancel button and per-gene progress are non-functional during the actual deletion computation

- **Location:** `cnapy/gui_elements/gene_essentiality_dialog.py:94-146` · **Category:** logic · **Verifier confidence:** high
- **Problem:** All of the heavy work happens inside the single blocking call single_gene_deletion(model) on line 98. The progress bar is only advanced to 0/total before that call (line 96) and then per-gene during the fast result-processing loop (line 146), which merely iterates an already-computed DataFrame. Likewise, the _cancel_requested checks on lines 100 and 108 are only reached AFTER single_gene_deletion has fully returned. Consequently, while the analysis is actually running (the slow part), the progress bar sits frozen at 0% and pressing Cancel has no effect until the whole computation is already finished. For a model with many genes this presents as a hung, uncancellable UI. The status label even says 'Analyzing gene N of M' (line 359) during what is really just table/result post-processing, misrepresenting actual progress.
- **Trigger:** Run the analysis on any non-trivial model and press Cancel while the computation is in progress; the cancel is ignored until single_gene_deletion() finishes, and the progress bar stays at 0% the entire time the solver is actually working.
- **Fix:** Drive single_gene_deletion in a way that allows progress/cancellation, e.g. iterate genes manually (knock out each gene with model context, optimize, restore) emitting progress and checking _cancel_requested per gene, or pass a cobra progress/processes hook. At minimum, do not advertise per-gene progress that does not reflect the real workload.

#### M41. paint() crashes with KeyError for a box whose reaction is not in the model but has a scenario value

- **Location:** `cnapy/gui_elements/map_view.py:1014-1017` · **Category:** crash · **Verifier confidence:** high
- **Problem:** ReactionBox explicitly supports reactions that are not in cobra_py_model (rebuild_scene lines 357-359 and ReactionBox.__init__ line 769 handle the off-model case). However paint() unconditionally calls cobra_py_model.reactions.get_by_id(self.id) whenever self.id is in scen_values, without guarding that the reaction is actually in the model. get_by_id raises KeyError for an unknown id. Because paint() is invoked by Qt during every repaint, the exception is raised inside the rendering path and is not caught here.
- **Trigger:** A box exists on the map for a reaction id that is not in the current cobra model (e.g. after a model was changed/reloaded while the map kept the box, or an imported map referencing an unknown id), and that id has an entry in scen_values. On the next repaint, get_by_id raises KeyError.
- **Fix:** Guard with `if self.id in self.map.appdata.project.cobra_py_model.reactions:` before calling get_by_id, and fall back to the non-scenario drawing branch (or skip bound comparison) when the reaction is absent from the model.

#### M42. Unvalidated float()/int() on Max. Solutions / Max. Size / Time Limit crashes compute

- **Location:** `cnapy/gui_elements/mcs_dialog.py:293-295` · **Category:** exception · **Verifier confidence:** high
- **Problem:** compute_optlang() parses three free-text QLineEdits with float(self.max_solu.text()), int(self.max_size.text()), and float(self.time_limit.text()) with no try/except. These fields are never validated by check_for_mcs_equation_errors() (which only inspects column-3 cells of the target/desired tables). Any non-numeric, empty, or out-of-domain text raises ValueError, which propagates out of the clicked slot as an unhandled exception. Notably int() cannot parse a float-like or 'inf' string: the default for max_size is '7' but the adjacent max_solu/time_limit fields default to 'inf', so a user who types 'inf' (a natural 'no limit' value) into Max. Size triggers int('inf') -> ValueError. Empty string in any of the three fields also raises.
- **Trigger:** User enters any non-integer (e.g. 'inf', '1e3', '', '7.5', 'abc') in the Max. Size field, or empty/non-numeric text in Max. Solutions or Time Limit, then clicks 'Compute MCS'.
- **Fix:** Wrap the three conversions in try/except ValueError and surface a QMessageBox.warning (as done for region parse errors), or validate these fields inside check_for_mcs_equation_errors before calling compute_optlang. For max_size, accept 'inf' by parsing as float and converting, or document that only finite integers are allowed and validate accordingly.

#### M43. rxn.name may be None → TypeError on name[:30] / name[:30] slicing crashes current-media display

- **Location:** `cnapy/gui_elements/media_management_dialog.py:772` · **Category:** crash · **Verifier confidence:** high
- **Problem:** In _update_current_media_display the reaction name is sliced with `name[:30]` (line 772) and also stored raw at line 764. In COBRApy a reaction's `.name` attribute can be an empty string but, for models loaded from certain formats or constructed programmatically, `.name` can be None. Slicing None raises `TypeError: 'NoneType' object is not subscriptable`, crashing the dialog's refresh.
- **Trigger:** Model contains an EX_ reaction whose .name is None (e.g. imported model without reaction names). Opening the dialog or clicking 'Refresh Current Status' raises TypeError.
- **Fix:** Coerce to string: `name = rxn.name or rxn.id` before slicing, or `QTableWidgetItem((name or '')[:30])`.

#### M44. Non-numeric bound from imported/custom JSON crashes table render via f"{bound:.2f}"

- **Location:** `cnapy/gui_elements/media_management_dialog.py:888` · **Category:** crash · **Verifier confidence:** high
- **Problem:** _on_custom_media_selected formats each component bound with `f"{bound:.2f}"` (line 888) and has NO surrounding try/except. Custom media are loaded from custom-media.json and can also originate from arbitrary user-imported JSON (via _import_media_from_json, which stores whatever 'components'/'patterns' values it finds without numeric validation). If a component value is a string or null, `format(bound, '.2f')` raises ValueError/TypeError when the user merely selects the item in the list, crashing the GUI.
- **Trigger:** User imports a JSON file whose components contain a non-float value (e.g. `"EX_glc__D_e": "-10"` or null), then clicks that media in the Saved Custom Media list. The currentItemChanged slot runs `f"{bound:.2f}"` on a str/None and raises.
- **Fix:** Validate/coerce values to float on import (reject or float() with try) and/or guard formatting: `try: txt = f"{float(bound):.2f}" except (TypeError, ValueError): txt = str(bound)`.

#### M45. Context-menu 'compute in/out fluxes' crashes when no metabolite row is selected

- **Location:** `cnapy/gui_elements/metabolite_list.py:212` · **Category:** crash · **Verifier confidence:** high
- **Problem:** emit_in_out_fluxes_action() calls self.metabolite_list.currentItem().text(0) with no None guard. The context menu that triggers this action is shown by on_context_menu() which only checks that the model has >0 metabolites (line 123), NOT that a row is currently selected. QTreeWidget.currentItem() returns None when there is no current item, so .text(0) raises AttributeError: 'NoneType' object has no attribute 'text', crashing the action.
- **Trigger:** Open a model with metabolites, ensure no row is selected (e.g. click empty area of the tree or after a programmatic clear of selection), right-click to open the context menu, and choose 'compute in/out fluxes for this metabolite'.
- **Fix:** Guard against None: `item = self.metabolite_list.currentItem(); if item is not None: self.computeInOutFlux.emit(item.text(0))`, and/or only enable the action when currentItem() is not None.

#### M46. prev()/next() infinite-loop hang (and IndexError on empty modes) when no mode is selected

- **Location:** `cnapy/gui_elements/mode_navigator.py:234-252` · **Category:** logic · **Verifier confidence:** high
- **Problem:** Both prev() and next() use `while True:` and only break when `self.selection[self.current]` is True. If the current selection array is all-False, neither loop ever breaks → the GUI thread hangs (application freeze). Additionally, if `len(self.appdata.project.modes) == 0`, the first iteration sets `self.current = len(modes) - 1 = -1` (prev) and then evaluates `self.selection[-1]` on a zero-length array, raising IndexError. An all-False selection is reachable: apply_selection() can leave `self.selection` partially/fully deselected if select() raises mid-iteration (num_selected is then stale and the `num_selected == 0` guard at line 296 uses the old value), or via the reset path interacting with a stale selection.
- **Trigger:** Press the '<' or '>' navigation button when `self.selection` contains no True entries (e.g., after a selection filter left no surviving modes due to the strain-design KeyError path partially deselecting, or an empty modes set). Result: infinite loop / frozen UI, or IndexError on empty modes.
- **Fix:** Guard before looping: `if self.num_selected == 0 or len(self.appdata.project.modes) == 0: return`. Bound the loop by the number of modes and break if no True is found to avoid an unbounded spin.

#### M47. Partial deselection and stale num_selected when select() raises mid-iteration

- **Location:** `cnapy/gui_elements/mode_navigator.py:287-307` · **Category:** logic · **Verifier confidence:** high
- **Problem:** apply_selection() calls select() inside a try/except. select() mutates self.selection in place as it iterates reactions and only recomputes self.num_selected at its final line (349). If a reaction ID is not found, `self.appdata.project.modes.reac_id.index(r)` raises ValueError partway through, leaving self.selection partially deselected while self.num_selected still holds the value from a PREVIOUS selection. The except shows a dialog but does not restore selection or recompute num_selected. The subsequent `if self.num_selected == 0` check (line 296) then operates on the stale count, so the empty-selection guard can be bypassed (num_selected > 0 while selection has fewer or zero True entries), and display_mode()/next() may navigate into an inconsistent selection (feeding the prev/next hang).
- **Trigger:** User enters a comma-separated selection where a valid reaction precedes an invalid/unknown reaction ID (e.g., `R1,bogus`). select() deselects modes for R1, then raises ValueError on `index('bogus')`; selection is left partially filtered with a stale num_selected, and the empty-selection guard uses the wrong count.
- **Fix:** On exception in apply_selection, call self.reset_selection() to restore a consistent state (all True, num_selected = len(modes)). Alternatively, validate all reaction IDs before mutating self.selection, and recompute num_selected in a finally block.

#### M48. find_blocked_reactions swallows all exceptions and falls back to a near-useless lb==ub==0 check, hiding solver/infeasibility failures and reporting wrong results

- **Location:** `cnapy/gui_elements/model_management_dialog.py:178-182` · **Category:** exception · **Verifier confidence:** high
- **Problem:** The FVA call is wrapped in `except Exception:` (line 178) which silently swallows ANY error — solver-not-installed, model infeasibility, optlang errors, MemoryError on large models — and falls back to scanning for reactions with lower_bound == 0 and upper_bound == 0. This fallback is not equivalent to FVA: it only catches reactions hard-bounded to zero, missing the vast majority of truly blocked reactions (which have non-zero bounds but cannot carry flux due to network structure). The user is shown an incomplete 'blocked reactions' list with no indication that FVA failed, producing silently incorrect scientific results. _scan_blocked's outer try/except (line 652) never sees the error because it was already swallowed here.
- **Trigger:** Run Find Blocked Reactions on a model where FVA fails (no LP solver configured, infeasible model, or solver error). Instead of an error, the user silently gets only reactions with bounds exactly [0,0], a scientifically wrong 'blocked' set.
- **Fix:** Narrow the except to expected solver exceptions and re-raise/report others; at minimum surface that FVA failed (so _scan_blocked can warn the user) rather than silently substituting a non-equivalent bounds check. Log the swallowed exception.

#### M49. FC reference defaults to first condition without honoring a user-chosen '(None)' selection

- **Location:** `cnapy/gui_elements/omics_integration_dialog.py:1140-1161` · **Category:** logic · **Verifier confidence:** high
- **Problem:** On every results refresh the FC reference combo is cleared and rebuilt. The previous selection is restored only if current_ref is in conditions (line 1146); '(None)' is NOT in conditions, so a user who explicitly selected '(None)' to disable fold-change has it overwritten by setCurrentIndex(1) (first real condition) on the next refresh (e.g. after _on_fc_reference_changed re-enters _update_results_display, or after re-running). This silently re-enables FC columns the user turned off and changes the displayed table semantics.
- **Trigger:** User selects '(None)' in FC Reference, then any action that calls _update_results_display (re-run analysis, or the signal round-trip) resets it to the first condition, recomputing FC columns unexpectedly.
- **Fix:** Treat '(None)' as a valid restorable selection: if current_ref == '(None)' or current_ref in conditions, restore it; only default to index 1 when current_ref is empty.

#### M50. Objective reaction id taken from first objective variable can be the reverse-variable, silently skipping the objective fix

- **Location:** `cnapy/gui_elements/omics_integration_dialog.py:442-449` · **Category:** api-misuse · **Verifier confidence:** high
- **Problem:** objective_rxn_id is derived from the first element of model.objective.variables via .name and then break. In cobra/optlang each reaction contributes TWO optlang variables (forward and reverse, named '<rid>' and '<rid>_reverse_xxxxx'). model.objective.variables is a set, so iteration order is non-deterministic; the first variable may be the reverse variable whose name is not a reaction id. objective_rxn_id then fails the 'in m.reactions' test at line 510, so the intended objective-fix path is skipped and the fragile objective_coefficient fallback (line 515-518) runs instead. The same wrong id is also added to excluded_reactions (line 448), so the real objective reaction may NOT be excluded from E-Flux scaling, letting expression data shrink the growth/objective bounds.
- **Trigger:** Any model whose objective variable set yields the reverse variable first (set ordering / hash dependent); the objective reaction is then neither correctly fixed in stage 2 nor excluded from E-Flux scaling.
- **Fix:** Derive the objective reaction(s) from model.reactions where reaction.objective_coefficient != 0 (or use cobra.util.solver.linear_reaction_coefficients(model)), rather than parsing optlang variable names.

#### M51. E-Flux2 objective-fix uses lower_bound, infeasible for negative or minimized objective values

- **Location:** `cnapy/gui_elements/omics_integration_dialog.py:509-518` · **Category:** math · **Verifier confidence:** high
- **Problem:** Step 2 fixes the objective at a fraction of optimum by setting the objective reaction's lower_bound to optimal_value * objective_fraction. This is only correct when the objective is a MAXIMIZE problem with a positive optimum. If optimal_value is negative (a maximize objective whose optimum is negative, e.g. a net-consumption objective, or a model whose objective direction is 'min' where optlang reports a signed value), then optimal_value * 0.99 is GREATER (less negative) than optimal_value. Forcing lower_bound above the actual optimum makes the QP/pFBA step infeasible or silently changes the constraint meaning, yielding wrong fluxes or a status that is not 'optimal'. For minimization objectives the correct cap is an upper_bound, not a lower_bound.
- **Trigger:** Run E-Flux2 on a model whose objective optimum is negative or whose objective direction is minimization (e.g. minimize a maintenance flux). Stage-2 QP/pFBA becomes infeasible or constrains incorrectly.
- **Fix:** Branch on model.objective.direction: for 'max' set lower_bound = optimal_value*frac when optimal_value>=0 (and handle negative optima explicitly); for 'min' set upper_bound = optimal_value/frac (or *something>=1). Generally constrain the objective expression with a Constraint rather than mutating a single reaction bound, and guard the sign of optimal_value.

#### M52. UnboundLocalError on old_id when reaction not found in list

- **Location:** `cnapy/gui_elements/reactions_list.py:411-424` · **Category:** crash · **Verifier confidence:** high
- **Problem:** In handle_changed_reaction, old_id is assigned ONLY inside the `if item.reaction == reaction:` branch of the loop. If the loop completes without finding a matching item, old_id is never bound, and line 424 `self.reactionChanged.emit(old_id, reaction)` raises UnboundLocalError. cobra.Reaction.__eq__ is identity-based (object.__eq__, confirmed), so the match relies on the exact same Reaction object instance being present as item.reaction in the tree. Any desync (item rebuilt/removed, reaction replaced, or the list cleared/repopulated between mask edit and signal delivery) leaves old_id unbound and crashes the slot.
- **Trigger:** reactionChanged is emitted (e.g. after applying an edit) for a reaction whose item is not present in the tree widget (list out of sync with the model), so the `if item.reaction == reaction` branch never executes and old_id stays undefined.
- **Fix:** Initialize `old_id = reaction.id` (or None) before the loop, and/or only emit when a matching item was found. e.g. set `old_id = None` before the loop and guard the emit: `if old_id is not None: self.reactionChanged.emit(old_id, reaction)`.

#### M53. RenameMapDialog silently corrupts map name on empty input / silently no-ops on collision

- **Location:** `cnapy/gui_elements/rename_map_dialog.py:40-48` · **Category:** data-loss · **Verifier confidence:** high
- **Problem:** apply() does no validation of new_name. (1) Empty name: if the user clears the field, new_name='' which is not a key in maps, so the guard passes and the code renames the current map to an empty string ('' becomes the dict key and the tab text becomes ''). This produces a map with an empty/invalid name and a blank tab, and may collide/serialize badly later. (2) Collision: if new_name already exists, the if-branch is skipped entirely but self.accept() still runs, so the dialog closes as if the rename succeeded while nothing changed and no error/warning is shown to the user; the rename is silently lost. Note: the dict key swap (maps[new_name] = maps.pop(old_name)) loses the original mapping with no recovery if anything in the if-body subsequently raised.
- **Trigger:** Open 'Change map name'; either clear the text and click Rename (creates a map keyed by '' with a blank tab) or type the name of an existing map (dialog closes, rename silently does nothing).
- **Fix:** Validate new_name: strip it, reject empty (`if not new_name.strip(): warn and return`), and on collision show a QMessageBox warning and keep the dialog open instead of silently accept().

#### M54. 'Use scenario bounds' on a fixed reaction (lb==ub) makes the analysis impossible to run

- **Location:** `cnapy/gui_elements/robustness_analysis_dialog.py:489-497,513-515` · **Category:** logic · **Verifier confidence:** high
- **Problem:** When 'Use scenario bounds as range' is enabled, _update_bounds_from_scenario reads scen_lb, scen_ub from scen_values[rxn_id] (or the reaction's own bounds) and sets min_spin=lb, max_spin=ub. For a reaction that is fixed in the scenario (a very common case, e.g. a defined glucose uptake fixed at -10, or any reaction the user pinned to a single value), lb == ub. _run_analysis then rejects the run with 'Min flux must be less than max flux' because of the x_min >= x_max guard. The feature that is meant to auto-populate a useful range instead produces an un-runnable state for exactly the fixed reactions users most often want to scan. The user must manually edit the values, defeating the checkbox, and may be confused since the dialog just refuses to run.
- **Trigger:** Enable 'Use scenario bounds as range' for any reaction that is fixed (lower_bound == upper_bound) in the current scenario, then press Run Analysis.
- **Fix:** When auto-populating from a fixed reaction, expand the range around the fixed value (e.g. [val - delta, val + delta] or use the model default bounds) rather than setting equal min/max, or special-case lb==ub to widen the sweep.

#### M55. Editing one bound silently discards a valid edit when the partner bound cell is non-numeric

- **Location:** `cnapy/gui_elements/scenario_tab.py:271-286` · **Category:** logic · **Verifier confidence:** high
- **Problem:** When the LB or UB cell of a scenario reaction is edited, both LB and UB cells are re-parsed via verify_bound (lines 272-273). The model update at lines 276-277 is gated on `not (isnan(lb) or isnan(ub))`. If the partner cell currently holds any non-parseable text, verify_bound returns NaN for it, so the entire update branch is skipped and the just-entered valid value is never written into scen_values.reactions, even though the user committed a correct number. The cell is repainted but the value is silently dropped; the model and scenario stay stale with no error shown to the user.
- **Trigger:** A scenario reaction has a non-numeric value in one bound cell (e.g. empty or a typo); the user enters a valid number in the other bound cell. The valid edit is discarded silently.
- **Fix:** Validate and persist each bound independently; only block the comparison/cross-validation (lb <= ub) when both are numeric, but still store the individually-valid edited cell and flag only the invalid one.

#### M56. Reaction bound sign-preservation skips bounds equal to zero, leaving directionality inconsistent with equation

- **Location:** `cnapy/gui_elements/scenario_tab.py:332-339` · **Category:** math · **Verifier confidence:** high
- **Problem:** After parsing a new reaction equation, the code attempts to copy the parsed reaction's lower/upper bound into the stored scenario bounds only when a sign conflict is detected. The conditions are `(stored < 0 and parsed >= 0) or (stored > 0 and parsed <= 0)`. When the stored bound is exactly 0.0, neither sub-condition is true, so the parsed bound is never propagated. Consequently, changing an equation to a reversible form (`<=>`, which makes cobra set lower_bound to a negative default) while the stored lower bound is 0 leaves the scenario reaction's lower bound at 0, contradicting the reversible equation the user typed. The reaction is silently treated as irreversible despite the entered equation, producing scientifically wrong flux bounds.
- **Trigger:** A scenario reaction whose stored lower bound is 0 has its equation edited to a reversible form (using `<=>`). The negative lower bound implied by reversibility is not propagated; the bound stays at 0.
- **Fix:** Adopt the parsed reaction's reversibility/bounds directly when the directionality implied by the equation differs from the stored bounds, including the boundary case of 0; or compare reversibility flags rather than strict sign inequalities.

#### M57. constraint_edited writes to wrong constraint via currentRow() instead of sender widget

- **Location:** `cnapy/gui_elements/scenario_tab.py:513-524` · **Category:** logic · **Verifier confidence:** medium
- **Problem:** constraint_edited is connected to every constraint row's textCorrect signal (line 488). textCorrect is emitted from check_text, which is invoked on every textChanged (utils.py line 220) and on focusOut (utils.py line 239), i.e. it can fire for a widget that is not the table's 'current' cell. The slot, however, identifies the target row with self.constraints.currentRow() (line 515) and indexes self.appdata.project.scen_values.constraints[row] (lines 519/523) using that. The signal carries no row identity and self.sender() is ignored. When the emitting widget's row differs from currentRow() — e.g. focusOut of one row while another is current, or programmatic text changes — the edited text is parsed and written into the wrong constraint entry, silently corrupting a different constraint.
- **Trigger:** Two or more scenario constraints exist; editing one constraint's field while the table's current cell points at a different row (e.g. focus-out events, or completer-driven text changes) causes the parsed constraint to overwrite the wrong constraints[] entry.
- **Fix:** Resolve the emitting widget via self.sender() and map it to its row (e.g. self.constraints.indexAt(widget.pos()) or store the row on the widget), then index constraints by that row instead of currentRow().

#### M58. Corrupt/non-list bookmarks JSON is accepted silently, then crashes the refresh loop

- **Location:** `cnapy/gui_elements/scenario_templates_dialog.py:250-257` · **Category:** crash · **Verifier confidence:** high
- **Problem:** load_bookmarks reads the bookmarks file with json.load and only resets to [] inside a bare 'except Exception'. If the file is valid JSON but not a list (e.g. it was corrupted to a JSON object {...}), json.load succeeds and self.bookmarks becomes a dict. No exception is raised, so the guard does not trigger. _refresh_bookmarks_list then does 'for bm in self.bookmarks' (iterating dict keys -> str) and 'bm['name']' (line 858) -> TypeError: string indices must be integers, crashing the dialog. The same applies if the JSON is a list of non-dict items.
- **Trigger:** scenario-bookmarks.json becomes a JSON object or a list of non-dicts (partial write, manual edit, or an older/incompatible format). Opening the Bookmarks tab raises TypeError.
- **Fix:** After json.load, validate isinstance(data, list) and that each element is a dict with a 'name' key; otherwise reset to [].

#### M59. Custom-template spin boxes clamp bounds to +/-1000, silently truncating larger model bounds

- **Location:** `cnapy/gui_elements/scenario_templates_dialog.py:648, 654` · **Category:** math · **Verifier confidence:** high
- **Problem:** The lower/upper bound QDoubleSpinBoxes used to add a custom reaction row are configured with setRange(-1000, 1000). Many COBRA / genome-scale models use bound magnitudes larger than 1000 (e.g. 100000 or 999999 for effectively-unbounded reactions). When the user dials in such a value via the spin box, Qt silently clamps it to +/-1000, so _add_custom_row writes a wrong, more-restrictive bound into the table (str(self.custom_ub_spin.value())). The subsequently applied flux bound is incorrect, changing FBA results without any warning.
- **Trigger:** User adds a custom reaction row intending a bound such as 100000 (a common unbounded value); the spin box caps it at 1000 and the wrong bound is applied to the scenario.
- **Fix:** Set the spin box range to match the model's bound magnitude (e.g. +/-1e6 or the cobra Configuration bounds) or allow free-form numeric entry instead of a clamped spin box.

#### M60. Non-merge apply with zero matched reactions silently wipes the scenario and cannot be undone

- **Location:** `cnapy/gui_elements/scenario_templates_dialog.py:763-782` · **Category:** data-loss · **Verifier confidence:** medium
- **Problem:** In _apply_template, when the 'Merge' checkbox is unchecked, the existing scenario is cleared via self.appdata.project.scen_values.clear_flux_values() BEFORE any template reaction is matched. If none of the template's pattern/reaction IDs exist in the loaded model (applied_count == 0), the method shows a 'No Matches' warning and returns at line 780-782 without restoring anything. The user's entire prior scenario is destroyed even though nothing was applied. Worse, clear_flux_values() mutates the dict directly and does NOT append a 'clear' entry to appdata.scenario_past (compare with the proper AppData.scen_values_clear() which records ('clear','all',0)). Because the destructive clear is invisible to the undo history, recreate_scenario_from_history()/undo cannot restore the lost values.
- **Trigger:** User has a populated scenario, opens the dialog, leaves 'Merge' unchecked, and applies a template whose oxygen/carbon/etc. exchange IDs do not match the current model's reaction naming (e.g. a model not using BiGG IDs). The scenario is cleared and unrecoverable via undo.
- **Fix:** Compute matched reactions first; only clear (and prefer routing the clear through appdata.scen_values_clear() so it is recorded in history) when applied_count>0, or restore the previous scen_values if no matches are found.

#### M61. Unguarded float()/int() on user-entered QLineEdit text raises ValueError

- **Location:** `cnapy/gui_elements/strain_design_dialog.py:1045, 1396, 1657-1659` · **Category:** exception · **Verifier confidence:** high
- **Problem:** Several numeric fields read free-form QLineEdit/QTableItem text and convert without try/except or validation. Non-numeric or locale-formatted input (e.g. comma decimal, stray text) raises an uncaught ValueError that propagates out of the GUI handler. Examples: min_gcp at line 1045 (`float(min_gcp)`), regulatory cost at line 1396 (`float(self.regulatory_itv_list.item(i,1).text())`), and the optlang path max-solutions/max-cost/time-limit at lines 1657-1659.
- **Trigger:** User types a non-numeric value into the Min. growth-coupling potential, a regulatory cost cell, Max. Solutions, Max. cost, or Time Limit field (e.g. '0,2' with a comma, or 'abc'), then triggers Check module / Compute. The conversion raises ValueError and the operation aborts with an unhandled traceback (and, in compute, a stuck override cursor).
- **Fix:** Validate/parse these inputs inside try/except and show a QMessageBox on failure, or use a QDoubleValidator/QIntValidator on the line edits.

#### M62. gene_names empty-check uses set("") (empty set) so the fallback to gene_ids never triggers

- **Location:** `cnapy/gui_elements/strain_design_dialog.py:156-159` · **Category:** logic · **Verifier confidence:** high
- **Problem:** `set(self.appdata.project.cobra_py_model.genes.list_attr('name')) != set('')` is intended to detect the case where all gene names are empty. But `set('')` evaluates to the empty set `set()`, not `{''}`. For any model with >=1 gene, the names-set is non-empty and thus `!= set()` is always True, so self.gene_names is always set to the raw names list (which may contain empty strings) and the fallback `self.gene_names = self.gene_ids` is effectively dead. Genes with empty-string names then propagate `''` into gene_names.
- **Trigger:** A COBRA model whose genes have empty 'name' attributes. gene_names then holds '' entries; line 1415 `gkoCost.update({self.gene_names[i]: ...})` keys cost dicts by '', collapsing distinct genes to one key (data loss), and line 1547 `self.gene_names.index(k)` returns a wrong/first index on reload.
- **Fix:** Compare against `{''}` (or test `if any(self.appdata...list_attr('name'))`): e.g. `if set(...list_attr('name')) != {''}: ... else: self.gene_names = self.gene_ids`.

#### M63. Default Cmin/Cmax line edits parsed with float() and no error handling

- **Location:** `cnapy/gui_elements/thermodynamics_dialog.py:278-279` · **Category:** exception · **Verifier confidence:** high
- **Problem:** min_default_conc and max_default_conc are user-editable QLineEdits (defaults '1e-6' and '0.2'). They are passed directly through `float(self.min_default_conc.text())` and `float(self.max_default_conc.text())` with no try/except, in contrast to the min_mdf field which IS guarded with a ValueError handler (lines 241-250). If the user clears the field or types a non-numeric value, float() raises ValueError. With the surrounding try/except commented out, the exception propagates uncaught, leaving the dialog open with a stuck BusyCursor and no user feedback.
- **Trigger:** Edit the 'Default Cmin'/'Default Cmax' field to an empty string or any non-numeric text (e.g. '0,2' under a comma locale, or 'abc'), then press Compute.
- **Fix:** Wrap both float() conversions in try/except ValueError that shows a QMessageBox.warning, resets the cursor to ArrowCursor, and returns — mirroring the existing min_mdf handling.

#### M64. Bottleneck solution flagged ALL_OK=True regardless of actual solver termination status

- **Location:** `cnapy/gui_elements/thermodynamics_dialog.py:319-320` · **Category:** logic · **Verifier confidence:** high
- **Problem:** For BOTTLENECK_ANALYSIS, any non-empty solution dict is unconditionally marked successful: `if solution != {}: solution[ALL_OK_KEY] = True`. This overrides whatever ALL_OK_KEY / TERMINATION_CONDITION_KEY the solver actually produced. If perform_lp_thermodynamic_bottleneck_analysis returns a populated dict that nonetheless represents a failed/suboptimal solve (time limit, infeasible-but-with-partial-data, etc.), the code forces ALL_OK_KEY=True and proceeds to set_boxes, silently writing potentially invalid bottleneck/flux/concentration results into the project instead of showing the appropriate warning that the OPTMDFPATHWAY/FBA paths would show.
- **Trigger:** A bottleneck analysis that returns a non-empty result dict but with a non-optimal termination condition (e.g. solver time/iteration limit) — results are accepted as valid and written into the model.
- **Fix:** Respect the solver-reported status: only set ALL_OK_KEY if not already present, and do not overwrite a False ALL_OK_KEY / non-success TERMINATION_CONDITION_KEY; let get_solution_from_thread decide based on the real status.

#### M65. No None-guard on yopt() return value; sol.status raises AttributeError

- **Location:** `cnapy/gui_elements/yield_optimization_dialog.py:116-123` · **Category:** crash · **Verifier confidence:** high
- **Problem:** `sol = yopt(...)` is immediately dereferenced via `sol.status` (line 123) with no check for None. straindesign.yopt (lptools.py) has a code path where the final `else: status = INFEASIBLE` branch (lptools.py line 720-721) sets a local variable but never executes a `return`, so yopt falls through and returns None. This occurs when the LFP solve in the den_sign loop yields a status that is neither OPTIMAL nor UNBOUNDED (e.g. solver-reported INFEASIBLE/ERROR/TIME_LIMIT on the scaled linear-fractional problem). The dialog then crashes with AttributeError: 'NoneType' object has no attribute 'status', and because BusyCursor is never reset, the application is left with a stuck busy cursor.
- **Trigger:** Run yield optimization where the LFP problem solver returns a status other than OPTIMAL/UNBOUNDED (e.g. numerically infeasible scaled LFP, solver time limit, or solver error). yopt returns None and the dialog dereferences it.
- **Fix:** After calling yopt, guard for None (and unexpected statuses): `if sol is None or sol.status not in (UNBOUNDED, OPTIMAL): show infeasible warning and return` before accessing sol.status / sol.objective_value. Also reset the cursor in a finally block.

#### M66. IndexError on Ctrl+C with no cell selected in QTableCopyable

- **Location:** `cnapy/utils.py:289-291` · **Category:** crash · **Verifier confidence:** high
- **Problem:** keyPressEvent handles Ctrl+C by sorting selectedIndexes() and immediately indexing the last element to compute max_column. When no cell is selected, selectedIndexes() returns an empty list, so copied_cells is [] and copied_cells[-1] raises IndexError, crashing the key handler (and propagating up through Qt's event dispatch).
- **Trigger:** User focuses a QTableCopyable table and presses Ctrl+C while no cell/row is selected (e.g. after a click that cleared the selection, or on an empty/freshly-populated table). selectedIndexes() returns [] and copied_cells[-1] raises IndexError.
- **Fix:** Guard against an empty selection before indexing, e.g. `if not copied_cells: return` immediately after computing copied_cells (still calling super().keyPressEvent first), so Ctrl+C with no selection is a no-op instead of a crash.

#### M67. Busy cursor never reset when an exception occurs mid-check (no try/finally)

- **Location:** `cnapy/utils_for_cnapy_api.py:19-73` · **Category:** resource-leak · **Verifier confidence:** high
- **Problem:** check_in_identifiers_org sets the widget to a busy cursor at line 19 and only resets it to the arrow cursor at line 73, the very last statement, with no try/finally. Any exception raised in between leaves the widget permanently stuck displaying the busy cursor. Because check_identifiers_org_entry can raise uncaught exceptions on plausible inputs (JSONDecodeError at line 91/106, KeyError at line 93/108 per the other findings), this is reachable in normal operation whenever the remote API misbehaves, degrading the UI until the widget is recreated.
- **Trigger:** Any exception during the loop body (e.g., the remote API returns non-JSON or JSON missing 'errorMessage'), which prevents reaching line 73.
- **Fix:** Wrap the loop body in try/finally and move widget.setCursor(Qt.ArrowCursor) into the finally block so the cursor is always restored.

#### M68. Connection error during the colon-split re-check is silently ignored

- **Location:** `cnapy/utils_for_cnapy_api.py:61-72` · **Category:** logic · **Verifier confidence:** high
- **Problem:** The connection_error flag is only inspected once per value, at line 50, immediately after the first check_identifiers_org_entry call. If that first check reports the pair as invalid AND the value contains ':', a second check_identifiers_org_entry is invoked at line 63. This second call can itself set connection_error = True (lines 89/104) when the network/server fails on the retry, but the result of that second call is never tested for connection_error. The code then proceeds to lines 65-71 and interprets is_key_valid / is_key_value_pair_valid (both default False on a connection error) as a genuine validation failure, marking the annotation cells red. A transient network failure on the retry is thus silently converted into a false 'invalid' result that corrupts the displayed validation state.
- **Trigger:** The first check marks the key:value pair invalid, the value contains a colon, and the network/server fails (RequestException) only on the second check_identifiers_org_entry call.
- **Fix:** After the re-check at line 63, re-test identifiers_org_result.connection_error (show the connection-error dialog and break) before applying the is_key_valid / is_key_value_pair_valid styling at lines 65-71.

#### M69. Unguarded .json() call crashes on non-JSON HTTP responses

- **Location:** `cnapy/utils_for_cnapy_api.py:91` · **Category:** crash · **Verifier confidence:** high
- **Problem:** result_object.json() is called immediately after a successful requests.get() without any guard. The try/except at lines 85-90 only catches requests.exceptions.RequestException, which covers connection/timeout errors but NOT a successful HTTP response whose body is not valid JSON. If identifiers.org returns an HTML error page, a 5xx/4xx with a non-JSON body, a redirect/maintenance page, or any malformed body, .json() raises requests.exceptions.JSONDecodeError (a subclass of ValueError, not RequestException), which propagates uncaught out of check_identifiers_org_entry and up through check_in_identifiers_org, crashing the validation operation. The same defect exists at line 106 for the key:value pair check.
- **Trigger:** identifiers.org returns a non-JSON response body (server error page, gateway 502/503 HTML, rate-limit page, or any malformed payload) while the TCP/HTTP request itself succeeds.
- **Fix:** Wrap result_object.json() in try/except (catching ValueError/requests.exceptions.JSONDecodeError) and treat a decode failure as a connection_error (or otherwise as an invalid result), e.g. set identifiers_org_result.connection_error = True and return. Optionally also check result_object.ok / status_code before parsing.

#### M70. KeyError when response JSON lacks the 'errorMessage' key

- **Location:** `cnapy/utils_for_cnapy_api.py:93` · **Category:** crash · **Verifier confidence:** high
- **Problem:** After parsing the response, the code indexes result_json["errorMessage"] directly. If the identifiers.org resolver returns a valid JSON object that does not contain the 'errorMessage' field (different schema for some endpoints, API version changes, or an unexpected success/error envelope), this raises KeyError, which is not caught and propagates out of the function, crashing the annotation-check workflow. The same unguarded access occurs at line 108 for the key:value pair response.
- **Trigger:** identifiers.org returns JSON without an 'errorMessage' key (e.g., an API schema change, a non-error payload, or an alternate error envelope).
- **Fix:** Use result_json.get("errorMessage") instead of result_json["errorMessage"], or guard with 'if "errorMessage" in result_json'. Treat a missing key as an invalid/unknown result rather than crashing.

---

## D. Low-severity findings (90)

_Minor robustness / correctness nits; unlikely triggers or cosmetic-but-wrong behavior._

| # | Location | Category | Title |
|---|----------|----------|-------|
| L1 | `cnapy/__main__.py:53-61` | logic | JAVA_HOME left pointing at a bogus path when no JVM is found in any site-packages |
| L2 | `cnapy/__main__.py:68-79` | logic | Installed `cnapy` console-script ignores command-line file arguments (argv never parsed) |
| L3 | `cnapy/appdata.py:177-202` | resource-leak | save_cnapy_config leaks the file handle and can truncate the config to empty if writing fails |
| L4 | `cnapy/appdata.py:452-475` | exception | except gurobipy.GurobiError raises AttributeError when gurobipy is not installed, masking the real model-construction error |
| L5 | `cnapy/application.py:111-130` | compat | make_light_palette uses deprecated/enum-mismatched QPalette role access vs make_dark_palette |
| L6 | `cnapy/application.py:222-242` | logic | QColor.fromRgb(int(color)) loses the alpha/format and assumes a plain decimal RGB int |
| L7 | `cnapy/core.py:184` | exception | KeyError on metabolites not in biomass reaction when variable_constituents passed directly |
| L8 | `cnapy/core.py:35` | compat | numpy.bool used — removed in NumPy >= 1.24, AttributeError crash |
| L9 | `cnapy/core.py:364` | exception | Unguarded get_by_id over scen_values.items() raises KeyError on stale scenario keys |
| L10 | `cnapy/ecmodel/ecmodel_builder.py:223-226` | logic | _build_gene_to_uniprot maps gene→protein from kcat entries by positional zip, silently dropping or mismatching when proteins/genes lengths differ |
| L11 | `cnapy/ecmodel/ecmodel_builder.py:726-730` | exception | _populate_ec_structure int() conversion of stoicho can raise on non-integer or NaN subunit counts |
| L12 | `cnapy/ecmodel/ecmodel_data.py:158-176` | logic | schema_version values below 1 trigger a fake migration that never completes, re-firing on every load |
| L13 | `cnapy/ecmodel/ecmodel_data.py:192` | logic | reset() discards self.ec entirely, dropping the ec-substructure's gecko_light flag while the top-level flag is preserved |
| L14 | `cnapy/ecmodel/ecmodel_data.py:96-107` | exception | to_dict requires every kcat entry to be a dict carrying all 7 keys; KcatEntry dataclass instances or sparse dicts break project save |
| L15 | `cnapy/ecmodel/expansion.py:61-71` | logic | clear_pseudoreaction_gpr silently ignores a missing extra_tsv path, dropping user-specified pseudoreaction IDs |
| L16 | `cnapy/ecmodel/yaml_io.py:167-168` | type | _as_number returns booleans unchanged, emitting YAML true/false for numeric fields |
| L17 | `cnapy/ecmodel/yaml_io.py:306-311` | data-loss | Metabolite charge truncated to int on load, dropping fractional/decimal charges and silently discarding non-numeric values |
| L18 | `cnapy/ecmodel/yaml_io.py:366-367` | logic | Legacy v1 schema_version in YAML is read but never migrated, so v1 sign convention is not flipped |
| L19 | `cnapy/flux_sampling.py:216-220` | math | Gaussian-noise post-processing can push sampled fluxes outside their feasible bounds, producing invalid 'samples' |
| L20 | `cnapy/flux_sampling.py:218` | math | std_fraction reinterpreted as an absolute standard deviation for near-zero reference fluxes (unit inconsistency) |
| L21 | `cnapy/gui_elements/batch_moma_room_dialog.py:670-671` | logic | Progress bar maximum miscomputed for reactions (excludes EX_ targets) so it never reaches 100% |
| L22 | `cnapy/gui_elements/batch_moma_room_dialog.py:712` | logic | Target/Genes radio buttons not disabled during run; status text uses live (mutable) radio state |
| L23 | `cnapy/gui_elements/central_widget.py:545-557` | exception | reaction_participation reads modes.fv_mat for mode_type<=1 without verifying modes is a FluxVectorContainer |
| L24 | `cnapy/gui_elements/central_widget.py:568-571` | dead-code | Dead/duplicated relative_participation flatten block can raise NameError for mode_type==2 path |
| L25 | `cnapy/gui_elements/config_cobrapy_dialog.py:127-133` | compat | Config file written with default platform encoding (no encoding specified) |
| L26 | `cnapy/gui_elements/config_cobrapy_dialog.py:98-135` | logic | current solver/tolerance partially applied and persisted even when current-model solver set fails |
| L27 | `cnapy/gui_elements/config_dialog.py:166` | api-misuse | self.close attribute shadows inherited QDialog.close() method, breaking programmatic close() |
| L28 | `cnapy/gui_elements/configuration_cplex.py:136-138` | concurrency | Synchronous subprocess.check_call (setup.py install) blocks the GUI event loop |
| L29 | `cnapy/gui_elements/configuration_cplex.py:88-89` | api-misuse | self.close push button shadows the inherited QDialog.close() method |
| L30 | `cnapy/gui_elements/download_dialog.py:42` | api-misuse | QPushButton assigned to self.close shadows QDialog.close(), making the inherited method uncallable |
| L31 | `cnapy/gui_elements/dynamic_fba_dialog.py:1073-1077` | exception | Header 'toggle all' assumes every Use cell is a checkbox; reads .isChecked() unguarded in all() |
| L32 | `cnapy/gui_elements/efmtool_dialog.py:61-70` | logic | efmtool error path leaves Compute button hidden with no way to retry, and abort wiring not reset |
| L33 | `cnapy/gui_elements/escher_map_view.py:149-165` | concurrency | retrieve_pos_and_zoom increments the save-completion semaphore synchronously before the async pos/zoom callbacks run, allowing the project to be saved with stale geometry |
| L34 | `cnapy/gui_elements/escher_map_view.py:290-301` | logic | last_accepted_value stored as an HTML attribute string in set_reaction_box_scenario_value but compared against a raw numeric string in value_changed, defeating the redundant-call guard |
| L35 | `cnapy/gui_elements/escher_map_view.py:73-80` | type | set_geometry concatenates non-string zoom/pos values directly into a JavaScript string, raising TypeError for default-initialized maps |
| L36 | `cnapy/gui_elements/flux_data_dialog.py:160-161` | logic | QSpinBox column range capped at 100 silently skips valid wider files |
| L37 | `cnapy/gui_elements/flux_data_dialog.py:459-462` | logic | Preview/stats compute min/max/mean over reactions absent from the model, misrepresenting matched data |
| L38 | `cnapy/gui_elements/flux_feasibility_dialog.py:617-625` | api-misuse | Aliasing of project scen_values.reactions into a transient Scenario in show_elemental_balance |
| L39 | `cnapy/gui_elements/flux_response_dialog.py:557-579` | resource-leak | FVA auto-detect runs solver synchronously on the GUI thread, freezing the UI |
| L40 | `cnapy/gui_elements/flux_response_dialog.py:588, 626-641` | concurrency | Re-running analysis does not check for an already-running worker, allowing thread reference to be overwritten |
| L41 | `cnapy/gui_elements/flux_response_dialog.py:724-731` | math | Fold-change division guard only excludes zero, so negative wild-type product flux yields a misleading negative fold |
| L42 | `cnapy/gui_elements/fseof_dialog.py:191-192, 214` | exception | Iterating model.reactions while building row dict relies on sol.fluxes containing every reaction id; KeyError on partial solution objects |
| L43 | `cnapy/gui_elements/fseof_dialog.py:602` | logic | 'Top positive' summary uses correlations[0] (sorted by \|r\|), so it is None when the strongest reaction is negatively correlated |
| L44 | `cnapy/gui_elements/fvseof_dialog.py:118-119` | math | Time-estimate label cleared on every stage but stale estimate persists across analyses isn't reset; minor — main numeric guard fine |
| L45 | `cnapy/gui_elements/gecko_unified_dialog.py:909-922` | api-misuse | kcat-file column header guidance vs parser mismatch / proteomics float() header heuristic can misclassify |
| L46 | `cnapy/gui_elements/gecko_unified_dialog.py:974-983` | crash | Unguarded FBA solver call on GUI thread in _flexibilize can raise an unhandled exception |
| L47 | `cnapy/gui_elements/gene_essentiality_dialog.py:449-459` | data-loss | Highlight-on-map silently wipes existing computed flux values and corrupts map coloring scale |
| L48 | `cnapy/gui_elements/gene_essentiality_dialog.py:458-459` | logic | central_widget.update() called for map refresh updates the currently active list tab instead of the map only |
| L49 | `cnapy/gui_elements/gene_essentiality_dialog.py:90-121` | math | Exact wild-type growth of 1e-9 makes every gene appear essential |
| L50 | `cnapy/gui_elements/gene_list.py:248-253` | exception | delete_selected_annotation catches IndexError but dict deletion raises KeyError |
| L51 | `cnapy/gui_elements/gene_list.py:85-98` | exception | handle_changed_gene raises NameError when the gene is not found in the tree |
| L52 | `cnapy/gui_elements/main_window.py:1872-1875` | resource-leak | open_project leaves window stuck on BusyCursor when model SBML inside the zip is invalid |
| L53 | `cnapy/gui_elements/map_view.py:156-166` | logic | Pinch-zoom applies scale but never updates _zoom tracking, desyncing zoom math |
| L54 | `cnapy/gui_elements/map_view.py:382-398` | exception | update_reaction looks up new reaction name without guarding for an off-model reaction |
| L55 | `cnapy/gui_elements/map_view.py:73-79` | logic | Initial zoom restoration is off-by-one, scale does not match saved zoom level |
| L56 | `cnapy/gui_elements/map_view.py:912-913` | logic | set_default_style mutates the shared appdata.default_color QColor in place |
| L57 | `cnapy/gui_elements/mcs_dialog.py:299-306` | logic | enum_method may be left unbound if no MCS-search radio button is checked |
| L58 | `cnapy/gui_elements/mcs_dialog.py:399-406` | exception | Bare 'except Exception' prints traceback but discards real solver/setup failures as a generic error box |
| L59 | `cnapy/gui_elements/mcs_dialog.py:454-462` | logic | Target 't' field validated unconditionally while empty/blank rows are silently skipped during computation, rejecting valid input |
| L60 | `cnapy/gui_elements/media_management_dialog.py:1025-1030` | data-loss | Import preserves arbitrary JSON values into custom media without numeric validation, propagating bad data |
| L61 | `cnapy/gui_elements/media_management_dialog.py:883` | crash | components/patterns explicitly null in JSON → {**None} TypeError in unguarded selection handler |
| L62 | `cnapy/gui_elements/media_management_dialog.py:897-899` | logic | Description QInputDialog cancel ignored: 'ok' overwritten and never checked, save proceeds with empty desc |
| L63 | `cnapy/gui_elements/metabolite_list.py:135` | logic | Metabolite identity matched with '==' instead of 'is', risking false matches on cobra equality semantics |
| L64 | `cnapy/gui_elements/metabolite_list.py:138-141` | data-loss | Renaming a metabolite that has a concentration leaves a stale concentration in the list and orphans the value |
| L65 | `cnapy/gui_elements/metabolite_list.py:434` | logic | validate_compartment validates the compartment string as a Metabolite name, not as a compartment |
| L66 | `cnapy/gui_elements/model_management_dialog.py:611-623` | logic | find_dead_end_metabolites and _scan_dead_ends use divergent reversibility branch ordering, risking inconsistent classification for reactions with lb<0 and ub<=0 vs lb>=0 and ub>0 |
| L67 | `cnapy/gui_elements/omics_integration_dialog.py:1002-1091` | resource-leak | run_btn left disabled and progress bar left visible if results-display/update path raises |
| L68 | `cnapy/gui_elements/omics_integration_dialog.py:253-261` | logic | Gene case-insensitive matching can map two model genes to the same expression entry / mismatched value |
| L69 | `cnapy/gui_elements/plot_customization_dialog.py:249-259` | exception | set_xscale/set_yscale exceptions are silently swallowed, leaving scale unapplied with no user feedback |
| L70 | `cnapy/gui_elements/plot_customization_dialog.py:261-275` | logic | Inverted or equal axis limits (min >= max) accepted silently, producing a reversed or degenerate axis |
| L71 | `cnapy/gui_elements/plot_customization_dialog.py:277-293` | exception | Rollback canvas.draw() is unguarded — a second draw failure crashes the slot |
| L72 | `cnapy/gui_elements/plot_space_dialog.py:166-195` | resource-leak | Busy cursor permanently stuck if load_scenario_into_model raises |
| L73 | `cnapy/gui_elements/reactions_list.py:288-290` | logic | df_val left stale (-inf) for reactions without a df value, corrupting DF-column sort |
| L74 | `cnapy/gui_elements/reactions_list.py:300-305` | logic | Flux column foreground forced to black overrides dark-mode white |
| L75 | `cnapy/gui_elements/reactions_list.py:477-498` | exception | update() leaves itemChanged disconnected and sorting disabled if an update raises |
| L76 | `cnapy/gui_elements/rename_map_dialog.py:21-22, 43` | crash | RenameMapDialog.apply pops a possibly-missing key (KeyError) when no map tab is active |
| L77 | `cnapy/gui_elements/robustness_analysis_dialog.py:197-203` | math | Bottleneck 'insensitive' decision uses an absolute 1e-6 gradient threshold, misclassifying small-magnitude objectives |
| L78 | `cnapy/gui_elements/robustness_analysis_dialog.py:89-100` | logic | original_bounds are read after the scenario has already overwritten them, so reported 'original bounds' are wrong |
| L79 | `cnapy/gui_elements/scenario_tab.py:204-214` | exception | update() dereferences possibly-None table items for flux column |
| L80 | `cnapy/gui_elements/scenario_templates_dialog.py:1006-1011` | data-loss | Duplicate bookmark name discards the entire custom template with no save/merge |
| L81 | `cnapy/gui_elements/strain_design_dialog.py:1425-1428` | resource-leak | save() leaves widget busy cursor set when the current module is invalid |
| L82 | `cnapy/gui_elements/strain_design_dialog.py:1565-1650` | resource-leak | compute() leaves the application override cursor stuck (busy) on every early-return path |
| L83 | `cnapy/gui_elements/strain_design_dialog.py:358` | api-misuse | clicked.connect(self.rem_module, True) passes True as the Qt connection type, not a slot argument |
| L84 | `cnapy/gui_elements/thermodynamics_dialog.py:212-321` | resource-leak | BusyCursor not reset when no reaction has dG0 / when match falls through to no objective |
| L85 | `cnapy/gui_elements/yield_optimization_dialog.py:105-168` | resource-leak | BusyCursor never reset on early return or exception, leaving GUI stuck |
| L86 | `cnapy/gui_elements/yield_optimization_dialog.py:173-176` | dead-code | Dead idx counter in set_boxes signals incomplete/abandoned logic |
| L87 | `cnapy/moma.py:106-108, 156` | api-misuse | has_milp_solver() inspects a fresh default cobra.Model, not the model actually being solved |
| L88 | `cnapy/utils.py:246` | logic | is_valid set to None instead of a boolean on empty-string rejection |
| L89 | `cnapy/utils.py:62-69` | exception | TypeError in update_selected when a model element has a None name/id |
| L90 | `cnapy/utils_for_cnapy_api.py:62-63` | logic | Colon-containing values are mis-split, dropping everything after the second colon |

---

## E. Rejected (false positives, 41)

_Claims the verifier overturned after reading the surrounding code — recorded for transparency. A few (e.g. MOMA/ROOM solver-status checking, all-zero-mode normalization) are defensible hardening ideas even though they were judged not-a-bug as written._

| Location | Rejected claim |
|----------|----------------|
| `cnapy/appdata.py:405-407` | Scenario.clear() overrides dict.clear() to reset ALL scenario state, surprising callers that expect only fluxes cleared |
| `cnapy/moma.py:91-92, 209-210` | Solver status never checked: infeasible/unbounded MOMA and ROOM silently return garbage (NaN/stale) fluxes |
| `cnapy/application.py:42-57` | excepthook calls a Qt GUI message box from any thread / before app init, risking secondary crash |
| `cnapy/ecmodel/ecmodel_builder.py:1129-1144` | apply_kcat_multiplier raises KeyError/ValueError on GECKO-light models and on reactions not connected to the enzyme metabolite |
| `cnapy/ecmodel/ecmodel_data.py:163-176` | upgrade() stamps schema_version=2 without flipping signs when called with cobra_model=None, mislabeling v1 data as v2 |
| `cnapy/gui_elements/main_window.py:3069-3077` | load_concentrations_json/xlsx replace_all flag silently ignored (amend == replace-all) |
| `cnapy/gui_elements/strain_design_dialog.py:1999, 2092-2093` | modes are built from set(self.assoc) ordering but selected by equivalence-class value, risking IndexError/wrong mode |
| `cnapy/gui_elements/central_widget.py:378-382, 1401-1406` | ConfirmMapDeleteDialog deletes wrong map after tab reorder due to captured stale index |
| `cnapy/gui_elements/central_widget.py:433-439` | ZeroDivisionError when normalizing an all-zero flux mode vector |
| `cnapy/gui_elements/central_widget.py:404-410` | Genes-tab search crashes when a found gene id is not a model gene (get_by_id KeyError) |
| `cnapy/appdata.py:128-132 (called from central_widget.py set_scen_value line 295-298)` | set_comp_value_as_scen_value silently drops a zero computed flux value |
| `cnapy/gui_elements/gecko_unified_dialog.py:957-989` | Busy cursor never reset on success path in _flexibilize; dialog left with permanent BusyCursor |
| `cnapy/gui_elements/gecko_unified_dialog.py:502` | Unguarded self.dialog.parent().centralWidget() chain crashes after model already mutated |
| `cnapy/gui_elements/omics_integration_dialog.py:1011-1014` | Scenario flux_constraints unpacked as (lb, ub) but scen_values stores single values, causing crash or wrong bounds |
| `cnapy/gui_elements/reactions_list.py:833-838` | apply() parses coefficient and bounds with float() without try, can raise after partial mutation |
| `cnapy/gui_elements/fvseof_dialog.py:157,167,169,177` | linspace over infinite reaction bounds when target FVA endpoints fall back to bounds produces NaN scan points |
| `cnapy/gui_elements/map_view.py:461-464` | remove_box uses unguarded del on reaction_boxes, raising KeyError when box was never built |
| `cnapy/gui_elements/dynamic_fba_dialog.py:356-369` | Custom Signal named 'finished' shadows QThread's built-in finished signal |
| `cnapy/gui_elements/dynamic_fba_dialog.py:253-269` | fluxes[rid] indexing raises KeyError for tracked reactions absent from the solution |
| `cnapy/gui_elements/robustness_analysis_dialog.py:89-90,116-117` | load_scenario_into_model called from worker QThread reads GUI-owned shared scen_values without synchronization |
| `cnapy/gui_elements/model_management_dialog.py:54-78` | simplify_gpr silently corrupts GPR rules with mixed AND/OR operators, producing logically wrong gene rules |
| `cnapy/gui_elements/model_management_dialog.py:260-273` | find_duplicate_gpr_genes false-positives on gene names containing duplicated word tokens, flagging reactions that have no actual duplicate genes |
| `cnapy/gui_elements/model_management_dialog.py:510, 554, 584, 602, 642, 663, 679` | _scan_gpr_duplicates / scan slots access self.appdata.project.cobra_py_model with no None/empty-model guard, raising AttributeError |
| `cnapy/gui_elements/model_management_dialog.py:213-214` | find_unbalanced_reactions skips multi-metabolite exchange/sink reactions and mislabels real multi-species reactions as exchanges |
| `cnapy/gui_elements/model_management_dialog.py:501` | _fix_all_gpr iterates and writes GPR but never disables Fix-All afterward correctly relies on rescan; class-level mutable dict _gpr_duplicates_data shared across all dialog instances |
| `cnapy/gui_elements/model_management_dialog.py:747` | _run_validation uses rxn.objective_coefficient which is deprecated/removed in modern cobrapy, raising AttributeError |
| `cnapy/gui_elements/flux_data_dialog.py:567-615` | log2FC custom coloring is restored before async repaint runs, so colors never apply |
| `cnapy/gui_elements/flux_data_dialog.py:393` | Remove buttons use stale captured row index, deleting the wrong condition after any removal |
| `cnapy/gui_elements/mode_navigator.py:269` | numpy.bool used as dtype raises AttributeError on NumPy >= 1.24 |
| `cnapy/gui_elements/mode_navigator.py:356-363` | Operator-precedence bug inflates strain-design mode size in size_histogram (zeros counted as present) |
| `cnapy/gui_elements/thermodynamics_dialog.py:319-321` | BOTTLENECK_ANALYSIS with empty solver result raises KeyError on ALL_OK_KEY |
| `cnapy/data/escher_cnapy.html:257-259` | updateReactionStoichiometry JS bridge drops stoichiometric coefficients equal to 0 due to truthy check, silently producing wrong stoichiometry on the map |
| `cnapy/gui_elements/config_dialog.py:284` | abs_tol uses unbounded QDoubleValidator; empty or locale-formatted text crashes apply() with ValueError |
| `cnapy/gui_elements/gene_list.py:91-95` | GPR rename silently fails for genes adjacent to ')' or at end of rule (no trailing space) |
| `cnapy/gui_elements/gene_list.py:242-245` | delete_gene removes wrong row when gene list is sorted (view index vs. take index) |
| `cnapy/gui_elements/plot_space_dialog.py:23-32, 187-194` | Plot button callable when straindesign missing; real cause hidden as a calculation error |
| `cnapy/gui_elements/yield_optimization_dialog.py:171-178` | set_boxes() can KeyError / write NaN comp_values from unbounded/undefined solutions |
| `cnapy/gui_elements/flux_optimization_dialog.py:139-140` | set_boxes assumes solution contains every reaction id -> KeyError on partial flux vectors |
| `cnapy/gui_elements/clipboard_calculator.py:91-95` | Left 'Current values' branch mutates comp_values in place, polluting it with scenario values |
| `cnapy/gui_elements/moma_room_reference_dialog.py:114-120` | get_reference indexes _conditions with checkedId-1, IndexError when no button is checked |
| `cnapy/gui_elements/solver_buttons.py:95-99` | select_solver result used as dict key can raise KeyError when solver sets diverge |
