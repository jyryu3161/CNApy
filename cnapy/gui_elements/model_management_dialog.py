"""Model Management Dialog for CNApy

This module provides model management utilities including:
- GPR simplification (removing duplicate genes in GPR rules)
- Finding dead-end metabolites
- Finding blocked reactions
- Finding orphan reactions
- Model validation
"""

from qtpy.QtCore import Qt, Slot
from qtpy.QtWidgets import (QDialog, QHBoxLayout, QLabel, QVBoxLayout,
                            QPushButton, QGroupBox, QTextEdit, QTabWidget,
                            QWidget, QMessageBox, QCheckBox, QListWidget,
                            QListWidgetItem, QAbstractItemView, QProgressBar,
                            QSplitter)
from qtpy.QtGui import QFont
import cobra
import re
from typing import List, Set, Dict, Tuple, Optional
from cnapy.appdata import AppData


def simplify_gpr(gpr_rule: str) -> str:
    """
    Simplify a GPR rule by removing duplicate genes while preserving logic.

    Examples:
        "(gene1 or gene1)" -> "gene1"
        "(gene1 and gene1)" -> "gene1"
        "(gene1 or gene2 or gene1)" -> "(gene1 or gene2)"
        "((gene1 or gene2) and (gene1 or gene3))" -> "((gene1 or gene2) and (gene1 or gene3))"
    """
    if not gpr_rule or gpr_rule.strip() == "":
        return ""

    def simplify_group(match):
        """Simplify a single group (within parentheses or the whole expression)."""
        content = match.group(1) if match.lastindex else match.group(0)

        # Check if this is an 'or' group or 'and' group
        if ' or ' in content and ' and ' not in content:
            # Split by 'or', remove duplicates while preserving order
            parts = [p.strip() for p in content.split(' or ')]
            seen = set()
            unique_parts = []
            for part in parts:
                if part not in seen:
                    seen.add(part)
                    unique_parts.append(part)
            if len(unique_parts) == 1:
                return unique_parts[0]
            return '(' + ' or '.join(unique_parts) + ')'
        elif ' and ' in content and ' or ' not in content:
            # Split by 'and', remove duplicates while preserving order
            parts = [p.strip() for p in content.split(' and ')]
            seen = set()
            unique_parts = []
            for part in parts:
                if part not in seen:
                    seen.add(part)
                    unique_parts.append(part)
            if len(unique_parts) == 1:
                return unique_parts[0]
            return '(' + ' and '.join(unique_parts) + ')'
        return '(' + content + ')'

    # Process from innermost parentheses outward
    result = gpr_rule
    prev_result = None

    while prev_result != result:
        prev_result = result
        # Match innermost parentheses containing only genes and operators (no nested parens)
        result = re.sub(r'\(([^()]+)\)', simplify_group, result)

    # Handle the case where there are no parentheses
    if '(' not in result:
        result = simplify_group(type('Match', (), {'group': lambda self, x=0: result, 'lastindex': 0})())
        # Remove outer parentheses if they were added
        if result.startswith('(') and result.endswith(')'):
            inner = result[1:-1]
            if inner.count('(') == inner.count(')'):
                result = inner

    # Clean up: remove redundant outer parentheses
    while result.startswith('(') and result.endswith(')'):
        inner = result[1:-1]
        # Check if the parentheses are matching and redundant
        depth = 0
        is_redundant = True
        for i, c in enumerate(inner):
            if c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
            if depth < 0:
                is_redundant = False
                break
        if is_redundant and depth == 0:
            result = inner
        else:
            break

    return result.strip()


def find_dead_end_metabolites(model: cobra.Model) -> List[cobra.Metabolite]:
    """
    Find dead-end metabolites (metabolites that can only be produced or only consumed).

    A dead-end metabolite is one where all reactions involving it have the metabolite
    on the same side (all produce or all consume).
    """
    dead_ends = []

    for met in model.metabolites:
        producing = 0
        consuming = 0

        for rxn in met.reactions:
            coef = rxn.metabolites[met]
            # Consider reversibility
            if rxn.lower_bound < 0 and rxn.upper_bound > 0:
                # Reversible - can both produce and consume
                producing += 1
                consuming += 1
            elif rxn.lower_bound >= 0:
                # Forward only
                if coef > 0:
                    producing += 1
                else:
                    consuming += 1
            else:
                # Reverse only
                if coef > 0:
                    consuming += 1
                else:
                    producing += 1

        # Dead-end if only produced or only consumed
        if (producing == 0 and consuming > 0) or (producing > 0 and consuming == 0):
            dead_ends.append(met)

    return dead_ends


def find_blocked_reactions(model: cobra.Model, tolerance: float = 1e-9) -> List[cobra.Reaction]:
    """
    Find blocked reactions (reactions that cannot carry flux).
    Uses Flux Variability Analysis approach.
    """
    blocked = []

    try:
        from cobra.flux_analysis import flux_variability_analysis
        fva_result = flux_variability_analysis(model, fraction_of_optimum=0.0)

        for rxn_id in fva_result.index:
            min_flux = fva_result.loc[rxn_id, 'minimum']
            max_flux = fva_result.loc[rxn_id, 'maximum']

            if abs(min_flux) < tolerance and abs(max_flux) < tolerance:
                blocked.append(model.reactions.get_by_id(rxn_id))
    except Exception as e:
        # Fallback: check if bounds are both zero
        for rxn in model.reactions:
            if rxn.lower_bound == 0 and rxn.upper_bound == 0:
                blocked.append(rxn)

    return blocked


def find_orphan_reactions(model: cobra.Model) -> List[cobra.Reaction]:
    """
    Find orphan reactions (reactions with metabolites that don't appear in any other reaction).
    """
    orphans = []

    for rxn in model.reactions:
        has_orphan_metabolite = False
        for met in rxn.metabolites:
            if len(met.reactions) == 1:
                has_orphan_metabolite = True
                break
        if has_orphan_metabolite:
            orphans.append(rxn)

    return orphans


def find_unbalanced_reactions(model: cobra.Model) -> List[Tuple[cobra.Reaction, Dict[str, float]]]:
    """
    Find reactions with unbalanced mass/charge.
    Returns list of (reaction, imbalance_dict) tuples.
    """
    unbalanced = []

    for rxn in model.reactions:
        # Skip exchange reactions
        if len(rxn.metabolites) == 1:
            continue

        try:
            imbalance = {}

            # Check mass balance for each element
            for met, coef in rxn.metabolites.items():
                if met.elements:
                    for element, count in met.elements.items():
                        if element not in imbalance:
                            imbalance[element] = 0
                        imbalance[element] += coef * count

            # Check charge balance
            charge_balance = 0
            for met, coef in rxn.metabolites.items():
                if met.charge is not None:
                    charge_balance += coef * met.charge

            if abs(charge_balance) > 1e-6:
                imbalance['charge'] = charge_balance

            # Remove balanced elements
            imbalance = {k: v for k, v in imbalance.items() if abs(v) > 1e-6}

            if imbalance:
                unbalanced.append((rxn, imbalance))
        except:
            pass  # Skip reactions that can't be checked

    return unbalanced


def find_duplicate_gpr_genes(model: cobra.Model) -> Dict[str, List[str]]:
    """
    Find reactions with duplicate genes in their GPR rules.
    Returns dict mapping reaction_id to list of duplicate genes.
    """
    duplicates = {}

    for rxn in model.reactions:
        if not rxn.gene_reaction_rule:
            continue

        # Extract all gene mentions
        genes = re.findall(r'\b(\w+)\b', rxn.gene_reaction_rule)
        # Filter out operators
        genes = [g for g in genes if g not in ('and', 'or')]

        # Find duplicates
        seen = set()
        dups = set()
        for g in genes:
            if g in seen:
                dups.add(g)
            seen.add(g)

        if dups:
            duplicates[rxn.id] = list(dups)

    return duplicates


class ModelManagementDialog(QDialog):
    """Dialog for model management utilities."""

    def __init__(self, appdata: AppData):
        QDialog.__init__(self)
        self.setWindowTitle("Model Management")
        self.appdata = appdata
        self.setMinimumSize(800, 600)

        self.layout = QVBoxLayout()

        # Create tab widget
        self.tabs = QTabWidget()

        # Tab 1: GPR Cleanup
        self.gpr_tab = self._create_gpr_tab()
        self.tabs.addTab(self.gpr_tab, "GPR Cleanup")

        # Tab 2: Dead-end Metabolites
        self.dead_end_tab = self._create_dead_end_tab()
        self.tabs.addTab(self.dead_end_tab, "Dead-end Metabolites")

        # Tab 3: Blocked Reactions
        self.blocked_tab = self._create_blocked_tab()
        self.tabs.addTab(self.blocked_tab, "Blocked Reactions")

        # Tab 4: Orphan Reactions
        self.orphan_tab = self._create_orphan_tab()
        self.tabs.addTab(self.orphan_tab, "Orphan Reactions")

        # Tab 5: Model Validation
        self.validation_tab = self._create_validation_tab()
        self.tabs.addTab(self.validation_tab, "Model Validation")

        self.layout.addWidget(self.tabs)

        # Close button
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()
        self.close_btn = QPushButton("Close")
        self.close_btn.clicked.connect(self.accept)
        btn_layout.addWidget(self.close_btn)
        self.layout.addLayout(btn_layout)

        self.setLayout(self.layout)

    def _create_gpr_tab(self) -> QWidget:
        """Create the GPR cleanup tab."""
        widget = QWidget()
        layout = QVBoxLayout()

        # Description
        desc = QLabel(
            "This tool finds and simplifies GPR (Gene-Protein-Reaction) rules that contain "
            "duplicate genes. For example, '(gene1 or gene1 or gene2)' will be simplified to "
            "'(gene1 or gene2)'."
        )
        desc.setWordWrap(True)
        layout.addWidget(desc)

        # Scan button
        btn_layout = QHBoxLayout()
        self.gpr_scan_btn = QPushButton("Scan for Duplicate Genes in GPR")
        self.gpr_scan_btn.clicked.connect(self._scan_gpr_duplicates)
        btn_layout.addWidget(self.gpr_scan_btn)
        btn_layout.addStretch()
        layout.addLayout(btn_layout)

        # Results list
        layout.addWidget(QLabel("Reactions with duplicate genes in GPR:"))

        splitter = QSplitter(Qt.Vertical)

        self.gpr_list = QListWidget()
        self.gpr_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.gpr_list.itemSelectionChanged.connect(self._show_gpr_details)
        splitter.addWidget(self.gpr_list)

        # Details panel
        details_widget = QWidget()
        details_layout = QVBoxLayout()
        details_layout.addWidget(QLabel("GPR Details:"))
        self.gpr_details = QTextEdit()
        self.gpr_details.setReadOnly(True)
        self.gpr_details.setFont(QFont("Courier"))
        details_layout.addWidget(self.gpr_details)
        details_widget.setLayout(details_layout)
        splitter.addWidget(details_widget)

        layout.addWidget(splitter)

        # Action buttons
        action_layout = QHBoxLayout()
        self.gpr_fix_selected_btn = QPushButton("Fix Selected")
        self.gpr_fix_selected_btn.clicked.connect(self._fix_selected_gpr)
        self.gpr_fix_selected_btn.setEnabled(False)
        action_layout.addWidget(self.gpr_fix_selected_btn)

        self.gpr_fix_all_btn = QPushButton("Fix All")
        self.gpr_fix_all_btn.clicked.connect(self._fix_all_gpr)
        self.gpr_fix_all_btn.setEnabled(False)
        action_layout.addWidget(self.gpr_fix_all_btn)
        action_layout.addStretch()
        layout.addLayout(action_layout)

        widget.setLayout(layout)
        return widget

    def _create_dead_end_tab(self) -> QWidget:
        """Create the dead-end metabolites tab."""
        widget = QWidget()
        layout = QVBoxLayout()

        desc = QLabel(
            "Dead-end metabolites are metabolites that can only be produced (no consumption) "
            "or only consumed (no production). These often indicate incomplete pathways."
        )
        desc.setWordWrap(True)
        layout.addWidget(desc)

        btn_layout = QHBoxLayout()
        self.dead_end_scan_btn = QPushButton("Find Dead-end Metabolites")
        self.dead_end_scan_btn.clicked.connect(self._scan_dead_ends)
        btn_layout.addWidget(self.dead_end_scan_btn)
        btn_layout.addStretch()
        layout.addLayout(btn_layout)

        layout.addWidget(QLabel("Dead-end metabolites:"))
        self.dead_end_list = QListWidget()
        self.dead_end_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        layout.addWidget(self.dead_end_list)

        widget.setLayout(layout)
        return widget

    def _create_blocked_tab(self) -> QWidget:
        """Create the blocked reactions tab."""
        widget = QWidget()
        layout = QVBoxLayout()

        desc = QLabel(
            "Blocked reactions are reactions that cannot carry any flux under the current "
            "model constraints. This analysis uses Flux Variability Analysis (FVA)."
        )
        desc.setWordWrap(True)
        layout.addWidget(desc)

        btn_layout = QHBoxLayout()
        self.blocked_scan_btn = QPushButton("Find Blocked Reactions (FVA)")
        self.blocked_scan_btn.clicked.connect(self._scan_blocked)
        btn_layout.addWidget(self.blocked_scan_btn)
        btn_layout.addStretch()
        layout.addLayout(btn_layout)

        layout.addWidget(QLabel("Blocked reactions:"))
        self.blocked_list = QListWidget()
        self.blocked_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        layout.addWidget(self.blocked_list)

        # Progress bar for FVA
        self.blocked_progress = QProgressBar()
        self.blocked_progress.setVisible(False)
        layout.addWidget(self.blocked_progress)

        widget.setLayout(layout)
        return widget

    def _create_orphan_tab(self) -> QWidget:
        """Create the orphan reactions tab."""
        widget = QWidget()
        layout = QVBoxLayout()

        desc = QLabel(
            "Orphan reactions are reactions that have at least one metabolite that doesn't "
            "appear in any other reaction. These may indicate disconnected pathways."
        )
        desc.setWordWrap(True)
        layout.addWidget(desc)

        btn_layout = QHBoxLayout()
        self.orphan_scan_btn = QPushButton("Find Orphan Reactions")
        self.orphan_scan_btn.clicked.connect(self._scan_orphans)
        btn_layout.addWidget(self.orphan_scan_btn)
        btn_layout.addStretch()
        layout.addLayout(btn_layout)

        layout.addWidget(QLabel("Orphan reactions:"))
        self.orphan_list = QListWidget()
        self.orphan_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        layout.addWidget(self.orphan_list)

        widget.setLayout(layout)
        return widget

    def _create_validation_tab(self) -> QWidget:
        """Create the model validation tab."""
        widget = QWidget()
        layout = QVBoxLayout()

        desc = QLabel(
            "Model validation checks for common issues such as unbalanced reactions "
            "(mass/charge imbalance), reactions with invalid bounds, and other problems."
        )
        desc.setWordWrap(True)
        layout.addWidget(desc)

        btn_layout = QHBoxLayout()
        self.validate_btn = QPushButton("Run Validation")
        self.validate_btn.clicked.connect(self._run_validation)
        btn_layout.addWidget(self.validate_btn)
        btn_layout.addStretch()
        layout.addLayout(btn_layout)

        layout.addWidget(QLabel("Validation results:"))
        self.validation_results = QTextEdit()
        self.validation_results.setReadOnly(True)
        self.validation_results.setFont(QFont("Courier"))
        layout.addWidget(self.validation_results)

        widget.setLayout(layout)
        return widget

    # Store duplicates data for fixing
    _gpr_duplicates_data: Dict[str, Tuple[str, str]] = {}  # rxn_id -> (original, simplified)

    @Slot()
    def _scan_gpr_duplicates(self):
        """Scan for duplicate genes in GPR rules."""
        self.gpr_list.clear()
        self.gpr_details.clear()
        self._gpr_duplicates_data = {}

        model = self.appdata.project.cobra_py_model
        duplicates = find_duplicate_gpr_genes(model)

        for rxn_id, dup_genes in duplicates.items():
            rxn = model.reactions.get_by_id(rxn_id)
            original = rxn.gene_reaction_rule
            simplified = simplify_gpr(original)

            if original != simplified:
                self._gpr_duplicates_data[rxn_id] = (original, simplified)
                item = QListWidgetItem(f"{rxn_id}: {', '.join(dup_genes)}")
                item.setData(Qt.UserRole, rxn_id)
                self.gpr_list.addItem(item)

        count = len(self._gpr_duplicates_data)
        if count == 0:
            self.gpr_details.setText("No duplicate genes found in GPR rules.")
        else:
            self.gpr_details.setText(f"Found {count} reactions with duplicate genes.")

        self.gpr_fix_all_btn.setEnabled(count > 0)
        self.gpr_fix_selected_btn.setEnabled(False)

    @Slot()
    def _show_gpr_details(self):
        """Show details for selected GPR."""
        selected = self.gpr_list.selectedItems()
        self.gpr_fix_selected_btn.setEnabled(len(selected) > 0)

        if len(selected) == 1:
            rxn_id = selected[0].data(Qt.UserRole)
            if rxn_id in self._gpr_duplicates_data:
                original, simplified = self._gpr_duplicates_data[rxn_id]
                self.gpr_details.setText(
                    f"Reaction: {rxn_id}\n\n"
                    f"Original GPR:\n{original}\n\n"
                    f"Simplified GPR:\n{simplified}"
                )

    @Slot()
    def _fix_selected_gpr(self):
        """Fix GPR for selected reactions."""
        selected = self.gpr_list.selectedItems()
        if not selected:
            return

        model = self.appdata.project.cobra_py_model
        fixed_count = 0

        for item in selected:
            rxn_id = item.data(Qt.UserRole)
            if rxn_id in self._gpr_duplicates_data:
                _, simplified = self._gpr_duplicates_data[rxn_id]
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.gene_reaction_rule = simplified
                fixed_count += 1

        QMessageBox.information(
            self, "GPR Fixed",
            f"Fixed GPR rules for {fixed_count} reaction(s)."
        )

        self.appdata.unsaved = True
        self._scan_gpr_duplicates()  # Refresh list

    @Slot()
    def _fix_all_gpr(self):
        """Fix all GPR rules with duplicates."""
        if not self._gpr_duplicates_data:
            return

        reply = QMessageBox.question(
            self, "Confirm Fix All",
            f"This will fix GPR rules for {len(self._gpr_duplicates_data)} reaction(s). Continue?",
            QMessageBox.Yes | QMessageBox.No
        )

        if reply == QMessageBox.Yes:
            model = self.appdata.project.cobra_py_model

            for rxn_id, (_, simplified) in self._gpr_duplicates_data.items():
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.gene_reaction_rule = simplified

            QMessageBox.information(
                self, "GPR Fixed",
                f"Fixed GPR rules for {len(self._gpr_duplicates_data)} reaction(s)."
            )

            self.appdata.unsaved = True
            self._scan_gpr_duplicates()  # Refresh list

    @Slot()
    def _scan_dead_ends(self):
        """Scan for dead-end metabolites."""
        self.dead_end_list.clear()

        model = self.appdata.project.cobra_py_model
        dead_ends = find_dead_end_metabolites(model)

        for met in dead_ends:
            # Determine if it's only produced or only consumed
            producing = 0
            consuming = 0
            for rxn in met.reactions:
                coef = rxn.metabolites[met]
                if rxn.lower_bound >= 0:  # Forward only
                    if coef > 0:
                        producing += 1
                    else:
                        consuming += 1
                elif rxn.upper_bound <= 0:  # Reverse only
                    if coef > 0:
                        consuming += 1
                    else:
                        producing += 1
                else:  # Reversible
                    producing += 1
                    consuming += 1

            status = "only produced" if consuming == 0 else "only consumed"
            item = QListWidgetItem(f"{met.id} ({met.name}) - {status}")
            item.setData(Qt.UserRole, met.id)
            self.dead_end_list.addItem(item)

        if len(dead_ends) == 0:
            self.dead_end_list.addItem("No dead-end metabolites found.")

    @Slot()
    def _scan_blocked(self):
        """Scan for blocked reactions."""
        self.blocked_list.clear()
        self.blocked_progress.setVisible(True)
        self.blocked_progress.setRange(0, 0)  # Indeterminate
        self.blocked_scan_btn.setEnabled(False)

        try:
            model = self.appdata.project.cobra_py_model
            blocked = find_blocked_reactions(model)

            for rxn in blocked:
                item = QListWidgetItem(f"{rxn.id}: {rxn.name}")
                item.setData(Qt.UserRole, rxn.id)
                self.blocked_list.addItem(item)

            if len(blocked) == 0:
                self.blocked_list.addItem("No blocked reactions found.")
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Error during FVA: {str(e)}")
        finally:
            self.blocked_progress.setVisible(False)
            self.blocked_scan_btn.setEnabled(True)

    @Slot()
    def _scan_orphans(self):
        """Scan for orphan reactions."""
        self.orphan_list.clear()

        model = self.appdata.project.cobra_py_model
        orphans = find_orphan_reactions(model)

        for rxn in orphans:
            # Find which metabolites are orphans
            orphan_mets = [met.id for met in rxn.metabolites if len(met.reactions) == 1]
            item = QListWidgetItem(f"{rxn.id}: {rxn.name} (orphan metabolites: {', '.join(orphan_mets)})")
            item.setData(Qt.UserRole, rxn.id)
            self.orphan_list.addItem(item)

        if len(orphans) == 0:
            self.orphan_list.addItem("No orphan reactions found.")

    @Slot()
    def _run_validation(self):
        """Run model validation."""
        model = self.appdata.project.cobra_py_model
        results = []

        results.append("=" * 60)
        results.append("MODEL VALIDATION REPORT")
        results.append("=" * 60)
        results.append("")

        # 1. Basic statistics
        results.append("BASIC STATISTICS:")
        results.append(f"  Reactions: {len(model.reactions)}")
        results.append(f"  Metabolites: {len(model.metabolites)}")
        results.append(f"  Genes: {len(model.genes)}")
        results.append("")

        # 2. Check for invalid bounds
        results.append("INVALID BOUNDS CHECK:")
        invalid_bounds = []
        for rxn in model.reactions:
            if rxn.lower_bound > rxn.upper_bound:
                invalid_bounds.append(f"  {rxn.id}: lower_bound ({rxn.lower_bound}) > upper_bound ({rxn.upper_bound})")

        if invalid_bounds:
            results.extend(invalid_bounds)
        else:
            results.append("  No invalid bounds found.")
        results.append("")

        # 3. Check for unbalanced reactions
        results.append("UNBALANCED REACTIONS:")
        unbalanced = find_unbalanced_reactions(model)

        if unbalanced:
            for rxn, imbalance in unbalanced[:20]:  # Show first 20
                imbalance_str = ", ".join([f"{k}: {v:+.2f}" for k, v in imbalance.items()])
                results.append(f"  {rxn.id}: {imbalance_str}")
            if len(unbalanced) > 20:
                results.append(f"  ... and {len(unbalanced) - 20} more")
        else:
            results.append("  No unbalanced reactions found.")
        results.append("")

        # 4. Check for reactions without genes
        results.append("REACTIONS WITHOUT GPR:")
        no_gpr = [rxn for rxn in model.reactions if not rxn.gene_reaction_rule]
        results.append(f"  Count: {len(no_gpr)} reactions")
        if no_gpr and len(no_gpr) <= 10:
            for rxn in no_gpr:
                results.append(f"    - {rxn.id}")
        results.append("")

        # 5. Check for duplicate reaction IDs (shouldn't happen, but check)
        results.append("DUPLICATE CHECK:")
        rxn_ids = [rxn.id for rxn in model.reactions]
        met_ids = [met.id for met in model.metabolites]
        gene_ids = [g.id for g in model.genes]

        dup_rxns = len(rxn_ids) - len(set(rxn_ids))
        dup_mets = len(met_ids) - len(set(met_ids))
        dup_genes = len(gene_ids) - len(set(gene_ids))

        results.append(f"  Duplicate reaction IDs: {dup_rxns}")
        results.append(f"  Duplicate metabolite IDs: {dup_mets}")
        results.append(f"  Duplicate gene IDs: {dup_genes}")
        results.append("")

        # 6. Check objective
        results.append("OBJECTIVE FUNCTION:")
        obj_rxns = [rxn.id for rxn in model.reactions if rxn.objective_coefficient != 0]
        if obj_rxns:
            results.append(f"  Objective reactions: {', '.join(obj_rxns)}")
        else:
            results.append("  No objective function defined.")
        results.append("")

        results.append("=" * 60)
        results.append("VALIDATION COMPLETE")
        results.append("=" * 60)

        self.validation_results.setText("\n".join(results))
