#!/usr/bin/env python3
#
# Copyright 2022 CNApy organization
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# -*- coding: utf-8 -*-

"""Omics integration dialog for transcriptome-based flux prediction.

This module provides two methods for integrating gene expression data:

1. LAD (Least Absolute Deviation) method:
   - Predicts flux distributions by minimizing the absolute deviation
     between predicted fluxes and expression-based targets.

2. E-Flux2 method:
   - Constrains reaction upper bounds based on gene expression levels.
   - Performs FBA followed by L2-norm minimization (QP, min Σv²) for a
     unique flux solution.  Falls back to pFBA (L1) if the solver does
     not support quadratic objectives.
   - Reference: Kim & Lun (2016), PLOS Computational Biology.

This module was added/modified as part of CNApy enhancements.
"""

import ast

import cobra
import numpy as np
import pandas as pd
from qtpy.QtCore import Qt, Slot
from qtpy.QtGui import QColor
from qtpy.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from cnapy.appdata import AppData
from cnapy.utils import no_scroll


class NumericTableWidgetItem(QTableWidgetItem):
    """QTableWidgetItem that sorts by absolute value (flux magnitude matters more than sign)."""

    def __init__(self, text: str, value: float):
        super().__init__(text)
        self._value = value

    def __lt__(self, other):
        if isinstance(other, NumericTableWidgetItem):
            # Handle NaN values - put them at the end
            if np.isnan(self._value):
                return False
            if np.isnan(other._value):
                return True
            # Sort by absolute value (magnitude)
            return abs(self._value) < abs(other._value)
        return super().__lt__(other)


def read_expression_data(filename: str) -> dict[str, float]:
    """
    Read gene expression data from a file.

    Supported formats:
    - CSV/TSV: First column = gene ID, second column = expression value
    - Excel: First column = gene ID, second column = expression value

    Parameters:
    -----------
    filename : str
        Path to the expression data file

    Returns:
    --------
    Dict[str, float]
        Dictionary mapping gene IDs to expression values
    """
    expression_info = {}

    if filename.endswith(".xlsx") or filename.endswith(".xls"):
        df = pd.read_excel(filename, header=0)
    elif filename.endswith(".csv"):
        df = pd.read_csv(filename, header=0)
    else:
        # Try tab-separated
        df = pd.read_csv(filename, sep="\t", header=0)

    # Assume first column is gene ID, second column is expression value
    if len(df.columns) >= 2:
        for idx, row in df.iterrows():
            gene_id = str(row.iloc[0]).strip()
            try:
                value = float(row.iloc[1])
                expression_info[gene_id] = value
            except (ValueError, TypeError):
                continue

    return expression_info


def read_multi_condition_expression_data(filename: str) -> tuple[list[str], dict[str, dict[str, float]]]:
    """
    Read gene expression data with multiple conditions from a file.

    Supported formats:
    - CSV/TSV: First column = gene ID, remaining columns = expression values per condition
    - Excel: First column = gene ID, remaining columns = expression values per condition

    Parameters:
    -----------
    filename : str
        Path to the expression data file

    Returns:
    --------
    tuple[list[str], dict[str, dict[str, float]]]
        - conditions: List of condition names (column headers)
        - data: Dictionary mapping gene IDs to dict of {condition: expression_value}

    Example file format:
        gene_id,WT,Drug1,Drug2
        b0001,5.23,3.12,4.56
        b0002,8.45,7.89,6.01
    """
    if filename.endswith(".xlsx") or filename.endswith(".xls"):
        df = pd.read_excel(filename, header=0, index_col=0)
    elif filename.endswith(".csv"):
        df = pd.read_csv(filename, header=0, index_col=0)
    else:
        # Try tab-separated
        df = pd.read_csv(filename, sep="\t", header=0, index_col=0)

    conditions = list(df.columns)
    data = {}

    for gene_id in df.index:
        gene_id_str = str(gene_id).strip()
        gene_values = {}
        for col in conditions:
            try:
                value = float(df.loc[gene_id, col])
                if not np.isnan(value):
                    gene_values[col] = value
            except (ValueError, TypeError):
                continue
        if gene_values:
            data[gene_id_str] = gene_values

    return conditions, data


def evaluate_gpr_expression(
    gpr_node: ast.AST, gene_expression: dict[str, float]
) -> float | None:
    """Recursively evaluate a GPR AST node with continuous expression values.

    AND (enzyme complex) → min of subunits
    OR  (isozymes)       → max of alternatives
    Gene leaf            → expression value (None if missing)

    Missing gene handling:
    - AND with any None → None (complex cannot form)
    - OR with some None → max of available values; all None → None

    Parameters
    ----------
    gpr_node : ast.AST
        A node from cobra's GPR AST (reaction.gpr.body).
    gene_expression : dict[str, float]
        Case-insensitive lookup dict mapping gene IDs to expression values.

    Returns
    -------
    float | None
        Evaluated expression level, or None if not evaluable.
    """
    if isinstance(gpr_node, ast.BoolOp):
        child_values = [
            evaluate_gpr_expression(child, gene_expression)
            for child in gpr_node.values
        ]

        if isinstance(gpr_node.op, ast.And):
            # Enzyme complex: all subunits required → min
            if any(v is None for v in child_values):
                return None
            return min(child_values)
        else:
            # Isozymes: any alternative suffices → max
            valid = [v for v in child_values if v is not None]
            if not valid:
                return None
            return max(valid)

    elif isinstance(gpr_node, ast.Name):
        return gene_expression.get(gpr_node.id)

    return None


def gene_expression_to_reaction_weights(
    model: cobra.Model, gene_expression: dict[str, float]
) -> dict[str, float]:
    """Convert gene expression values to reaction weights using GPR rules.

    Evaluates GPR boolean rules recursively:
    AND (enzyme complex) → min, OR (isozymes) → max.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model.
    gene_expression : dict[str, float]
        Dictionary mapping gene IDs to expression values.

    Returns
    -------
    dict[str, float]
        Dictionary mapping reaction IDs to weights.
    """
    # Build case-insensitive lookup for gene expression
    expr_lookup: dict[str, float] = {}
    for gene in model.genes:
        gene_id = gene.id
        if gene_id in gene_expression:
            expr_lookup[gene_id] = gene_expression[gene_id]
        elif gene_id.upper() in gene_expression:
            expr_lookup[gene_id] = gene_expression[gene_id.upper()]
        elif gene_id.lower() in gene_expression:
            expr_lookup[gene_id] = gene_expression[gene_id.lower()]

    reaction_weights = {}

    for reaction in model.reactions:
        if not reaction.genes:
            continue

        gpr = reaction.gpr
        if gpr is None or gpr.body is None:
            continue

        weight = evaluate_gpr_expression(gpr.body, expr_lookup)
        if weight is not None:
            reaction_weights[reaction.id] = weight

    return reaction_weights


def run_lad_fitting(
    model: cobra.Model,
    reaction_weights: dict[str, float],
    flux_constraints: dict[str, tuple[float, float]] | None = None,
    weight_threshold: float = 0.01,
    scaling_factor: float = 1.0,
) -> tuple[str, float | None, dict[str, float] | None]:
    """
    Run LAD (Least Absolute Deviation) fitting to predict fluxes from expression data.

    Parameters:
    -----------
    model : cobra.Model
        The metabolic model
    reaction_weights : Dict[str, float]
        Dictionary mapping reaction IDs to target weights (from expression data)
    flux_constraints : Dict[str, Tuple[float, float]], optional
        Additional flux constraints as {reaction_id: (lb, ub)}
    weight_threshold : float
        Minimum weight to include a reaction in the fitting
    scaling_factor : float
        Scaling factor for target fluxes

    Returns:
    --------
    Tuple[str, Optional[float], Optional[Dict[str, float]]]
        (status, objective_value, flux_distribution)
    """
    from optlang.symbolics import Add

    # Filter reactions by weight threshold
    target_reactions = {
        rid: w for rid, w in reaction_weights.items() if rid in model.reactions and abs(w) > weight_threshold
    }

    if not target_reactions:
        return "no_targets", None, None

    with model as m:
        # Apply flux constraints if provided
        if flux_constraints:
            for rid, (lb, ub) in flux_constraints.items():
                if rid in m.reactions:
                    rxn = m.reactions.get_by_id(rid)
                    rxn.bounds = (lb, ub)

        # Create deviation variables for target reactions
        vars_to_add = []
        cons_to_add = []

        # Store deviation variables for later access
        delta_pos_vars = {}
        delta_neg_vars = {}

        for rid, weight in target_reactions.items():
            rxn = m.reactions.get_by_id(rid)
            target_flux = abs(weight) * scaling_factor

            # Create positive and negative deviation variables
            # flux = target + delta_pos - delta_neg
            # Total deviation = delta_pos + delta_neg
            delta_pos = m.problem.Variable(f"delta_pos_{rid}", lb=0)
            delta_neg = m.problem.Variable(f"delta_neg_{rid}", lb=0)
            delta_pos_vars[rid] = delta_pos
            delta_neg_vars[rid] = delta_neg
            vars_to_add.extend([delta_pos, delta_neg])

            # flux = target + delta_pos - delta_neg
            cons_decomp = m.problem.Constraint(
                rxn.flux_expression - delta_pos + delta_neg, lb=target_flux, ub=target_flux, name=f"lad_decomp_{rid}"
            )
            cons_to_add.append(cons_decomp)

        m.add_cons_vars(vars_to_add)
        m.add_cons_vars(cons_to_add)

        # Objective: minimize sum of deviations
        objective_terms = []
        for rid in target_reactions:
            objective_terms.extend([delta_pos_vars[rid], delta_neg_vars[rid]])

        m.objective = m.problem.Objective(Add(*objective_terms), direction="min")

        # Optimize
        solution = m.optimize()

        if solution.status == "optimal":
            flux_distribution = {rxn.id: solution.fluxes[rxn.id] for rxn in m.reactions}
            return "optimal", solution.objective_value, flux_distribution
        else:
            return solution.status, None, None


def run_eflux2(
    model: cobra.Model,
    reaction_weights: dict[str, float],
    flux_constraints: dict[str, tuple[float, float]] | None = None,
    weight_threshold: float = 0.0,
    min_scale: float = 0.001,
    objective_fraction: float = 0.99,
    exclude_exchange: bool = True,
    additional_excluded_reactions: list[str] | None = None,
    normalization_percentile: float = 100.0,
) -> tuple[str, float | None, dict[str, float] | None, str]:
    """Run E-Flux2 analysis to predict fluxes from expression data.

    E-Flux2 scales each constrained reaction's bounds by its normalized
    expression ratio (expr / max_expr), then performs FBA to find the optimal
    objective value, fixes the objective as a linear constraint at
    ``objective_fraction * optimum``, and finally minimizes the L2-norm of the
    flux vector (QP, ``min Σ vᵢ²``) for a unique solution.  If the solver does
    not support quadratic objectives, the function falls back to pFBA (L1) and
    finally to plain FBA — restoring the original biological objective in
    between so the fallback solves the right problem.

    Reference: Kim MK, Lane A, Kelley JJ, Lun DS (2016). E-Flux2 and SPOT:
    Validated Methods for Inferring Intracellular Metabolic Flux Distributions
    from Transcriptomic Data. PLoS ONE 11(6): e0157101.

    Parameters
    ----------
    model : cobra.Model
        Metabolic model. Not mutated — all changes happen inside a ``with model``
        context manager.
    reaction_weights : dict[str, float]
        Reaction-level expression weights from
        :func:`gene_expression_to_reaction_weights`.
    flux_constraints : dict[str, tuple[float, float]], optional
        Per-reaction ``(lb, ub)`` overrides (e.g. scenario constraints).
        Reactions listed here keep these bounds — expression-based scaling is
        NOT applied to them.
    weight_threshold : float, default 0.0
        Raw-expression noise filter (in the same units as the input weights,
        e.g. TPM).  Reactions whose raw expression is at or below this
        threshold are treated as minimally expressed — their bounds are
        scaled by ``min_scale`` (NOT left unconstrained, which would defeat
        E-Flux's continuous scaling).  Above-threshold reactions are scaled by
        the usual normalized ratio.  The threshold also gates which weights
        participate in the max-normalization, so a handful of noise spikes
        cannot inflate ``max_weight`` and crush every other reaction.
    min_scale : float, default 0.001
        Numerical-stability floor on the normalized expression ratio (also the
        scaling applied to below-threshold "noise" reactions).
    objective_fraction : float, default 0.99
        Fraction of the FBA optimum that the L2 step must still achieve.
    exclude_exchange : bool, default True
        Exclude EX_-prefixed exchange reactions and reactions whose id
        contains "biomass"/"growth" from expression-based scaling.
    additional_excluded_reactions : list[str], optional
        Extra reaction IDs to exclude from scaling.
    normalization_percentile : float, default 100.0
        Percentile (1–100) of the above-threshold weights used as the
        normalization denominator.  ``100`` is the paper-faithful max
        normalization.  Lower values make the denominator robust against
        extreme-expression outliers — recommended ``99`` for real RNA-seq
        data where a single highly-expressed gene/complex would otherwise
        crush every other reaction's bound.  Normalized values above 1.0
        (the outliers themselves) are clipped to 1.0.

    Returns
    -------
    tuple[str, float | None, dict[str, float] | None, str]
        ``(status, objective_value, flux_distribution, l2_method)``.
        ``l2_method`` is one of ``"qp"``, ``"pfba"``, ``"fba"``, or ``""``.
    """
    from cobra.flux_analysis import pfba
    from cobra.util import linear_reaction_coefficients
    from optlang.symbolics import Add

    if not reaction_weights:
        return "no_targets", None, None, ""

    # Keep only reactions that actually exist in the model.  Negative weights
    # are interpreted by magnitude.
    valid_weights = {
        rid: abs(w)
        for rid, w in reaction_weights.items()
        if rid in model.reactions
    }
    if not valid_weights:
        return "no_targets", None, None, ""

    # Separate above-threshold (real signal) from below-threshold (noise).
    # IMPORTANT: ``weight_threshold`` is interpreted in the *raw* units of the
    # incoming expression data, not as a normalized [0,1] floor.  Below-threshold
    # reactions are NOT left unconstrained — they get scaled by ``min_scale``,
    # matching E-Flux's "low expression → tight bound" semantics.
    above_threshold = {rid: w for rid, w in valid_weights.items() if w > weight_threshold}
    if not above_threshold:
        return "no_targets", None, None, ""

    # Normalize against the percentile of above-threshold weights.  At
    # percentile=100 this is the paper-faithful max; lower percentiles make
    # the denominator robust to a handful of extreme outliers (e.g. a single
    # very-highly-expressed gene/complex) that would otherwise crush every
    # other reaction's bound.
    if not 1.0 <= normalization_percentile <= 100.0:
        raise ValueError(
            f"normalization_percentile must be in [1, 100], got {normalization_percentile}"
        )
    if normalization_percentile >= 100.0:
        max_weight = max(above_threshold.values())
    else:
        max_weight = float(
            np.percentile(list(above_threshold.values()), normalization_percentile)
        )
    if max_weight <= 0:
        max_weight = 1.0

    # Identify ALL objective reactions (linear-combination objectives are
    # supported) and exclude them from bound scaling.
    objective_coeffs = linear_reaction_coefficients(model)  # {Reaction: coef}
    excluded_reactions: set[str] = {rxn.id for rxn in objective_coeffs}
    if exclude_exchange:
        for rxn in model.reactions:
            rid_lower = rxn.id.lower()
            if rxn.id.startswith("EX_") or "biomass" in rid_lower or "growth" in rid_lower:
                excluded_reactions.add(rxn.id)
    if additional_excluded_reactions:
        excluded_reactions.update(additional_excluded_reactions)

    with model as m:
        # Scenario constraints take precedence and skip expression scaling.
        if flux_constraints:
            for rid, (lb, ub) in flux_constraints.items():
                if rid in m.reactions:
                    m.reactions.get_by_id(rid).bounds = (lb, ub)

        constrained_count = 0
        for rxn in m.reactions:
            if rxn.id in excluded_reactions:
                continue
            if flux_constraints and rxn.id in flux_constraints:
                continue
            if rxn.id not in valid_weights:
                # No expression data → keep original bounds (documented choice).
                continue
            # Skip forced-flow reactions (e.g. ATP maintenance with lb > 0).
            # Their bounds are explicit modeling constraints, not catalytic
            # capacity, so expression scaling should not override them.
            if rxn.lower_bound > 0 or rxn.upper_bound < 0:
                continue

            raw = valid_weights[rxn.id]
            if raw <= weight_threshold:
                # Noise-floor reactions: tighten bounds to ``min_scale × orig``.
                normalized = min_scale
            else:
                normalized = raw / max_weight
                # Numerical-stability floor / ceiling on the normalized ratio.
                normalized = max(min_scale, min(1.0, normalized))

            if rxn.upper_bound > 0:
                rxn.upper_bound = normalized * rxn.upper_bound
            if rxn.lower_bound < 0:
                rxn.lower_bound = normalized * rxn.lower_bound
            constrained_count += 1

        # Snapshot the original biological objective so we can restore it
        # before any non-QP fallback.  We capture the optlang expression and
        # the direction; cobra's `with` context will revert anything we add.
        original_obj_expression = m.objective.expression
        original_obj_direction = m.objective.direction

        # Step 1 — FBA to find the optimal biological objective under the
        # expression-scaled bounds.
        try:
            solution1 = m.optimize()
        except Exception as e:
            return f"optimization_error: {e!s}", None, None, ""
        if solution1.status != "optimal":
            return solution1.status, None, None, ""

        optimal_value = solution1.objective_value
        stage1_fluxes = {rxn.id: solution1.fluxes[rxn.id] for rxn in m.reactions}

        # Step 2 — fix the biological objective via a linear constraint, then
        # try QP (true L2) → pFBA (L1) → plain FBA.
        if objective_coeffs:
            obj_linear_expr = Add(*[
                coef * m.reactions.get_by_id(rxn.id).flux_expression
                for rxn, coef in objective_coeffs.items()
            ])
            tolerance = abs(optimal_value) * (1.0 - objective_fraction)
            if original_obj_direction == "max":
                obj_floor = m.problem.Constraint(
                    obj_linear_expr,
                    lb=optimal_value - tolerance,
                    name="_eflux2_obj_floor",
                )
            else:
                obj_floor = m.problem.Constraint(
                    obj_linear_expr,
                    ub=optimal_value + tolerance,
                    name="_eflux2_obj_floor",
                )
            m.add_cons_vars([obj_floor])

        # --- QP (true L2: min Σ vᵢ²) ---
        flux_distribution: dict[str, float] | None = None
        l2_method = ""
        try:
            quadratic_terms = Add(*[rxn.flux_expression ** 2 for rxn in m.reactions])
            m.objective = m.problem.Objective(quadratic_terms, direction="min")
            qp_solution = m.optimize()
            if qp_solution.status == "optimal":
                flux_distribution = {
                    rxn.id: qp_solution.fluxes[rxn.id] for rxn in m.reactions
                }
                l2_method = "qp"
        except Exception:
            pass  # QP not supported — fall through to pFBA.

        # --- Fallback: pFBA (L1) ---
        # CRITICAL: restore the biological objective first so that pfba()'s
        # internal `fix_objective_as_constraint` step fixes the *biological*
        # objective at its optimum, not the leftover QP objective from above.
        if flux_distribution is None:
            m.objective = m.problem.Objective(
                original_obj_expression, direction=original_obj_direction
            )
            try:
                pfba_solution = pfba(m)
                if pfba_solution.status == "optimal":
                    flux_distribution = {
                        rxn.id: pfba_solution.fluxes[rxn.id] for rxn in m.reactions
                    }
                    l2_method = "pfba"
            except Exception:
                pass

        if flux_distribution is None:
            # Final fallback: Stage-1 FBA fluxes.
            return "optimal", optimal_value, stage1_fluxes, "fba"

        return "optimal", optimal_value, flux_distribution, l2_method


class OmicsIntegrationDialog(QDialog):
    """Dialog for Omics integration using LAD or E-Flux2 methods.

    Supports both single-condition and multi-condition expression data files.
    Multi-condition files allow batch analysis across multiple experimental conditions.
    """

    def __init__(self, appdata: AppData, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Omics Integration - Multi-Condition Flux Prediction")
        self.setMinimumWidth(900)
        self.setMinimumHeight(900)
        self.resize(1200, 950)
        self.appdata = appdata

        # Single-condition data (legacy support)
        self.expression_data: dict[str, float] = {}
        self.reaction_weights: dict[str, float] = {}

        # Multi-condition data
        self.conditions: list[str] = []
        self.multi_expression_data: dict[str, dict[str, float]] = {}
        self.multi_results: dict[str, dict[str, float]] = {}
        self.is_multi_condition = False

        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout()

        # File selection group
        file_group = QGroupBox("Expression Data")
        file_layout = QVBoxLayout()

        file_row = QHBoxLayout()
        file_row.addWidget(QLabel("Expression file:"))
        self.file_edit = QLineEdit()
        self.file_edit.setPlaceholderText("Select CSV, TSV, or Excel file with gene expression data")
        self.file_edit.setReadOnly(True)
        file_row.addWidget(self.file_edit)
        self.browse_btn = QPushButton("Browse...")
        self.browse_btn.setAutoDefault(False)
        self.browse_btn.clicked.connect(self._browse_file)
        file_row.addWidget(self.browse_btn)
        file_layout.addLayout(file_row)

        info_row = QHBoxLayout()
        self.gene_count_label = QLabel("Genes loaded: 0")
        info_row.addWidget(self.gene_count_label)
        self.matched_count_label = QLabel("Matched to model: 0")
        info_row.addWidget(self.matched_count_label)
        self.conditions_label = QLabel("Conditions: 0")
        info_row.addWidget(self.conditions_label)
        info_row.addStretch()
        file_layout.addLayout(info_row)

        # Example file format display
        example_label = QLabel(
            "<b>Expected File Format (Multi-Condition):</b><br>"
            "<pre style='background-color: #f0f0f0; padding: 8px; font-family: monospace;'>"
            "gene_id,WT,Drug1,Drug2\n"
            "b0001,5.23,3.12,4.56\n"
            "b0002,8.45,7.89,6.01\n"
            "...</pre>"
            "<b>Requirements:</b><br>"
            "• <b>Column 1</b>: Gene ID (must match model gene IDs)<br>"
            "• <b>Column 2+</b>: Expression values for each condition<br>"
            "• <b>Header row</b>: Condition names (WT, Drug1, Drug2, etc.)"
        )
        example_label.setWordWrap(True)
        example_label.setStyleSheet(
            "QLabel { background-color: #fafafa; padding: 10px; border: 1px solid #ddd; border-radius: 4px; }"
        )
        file_layout.addWidget(example_label)

        # Data preview label (shown after loading)
        self.data_preview_label = QLabel("")
        self.data_preview_label.setWordWrap(True)
        self.data_preview_label.setVisible(False)
        self.data_preview_label.setStyleSheet(
            "QLabel { background-color: #e8f5e9; padding: 10px; border: 1px solid #c8e6c9; border-radius: 4px; }"
        )
        file_layout.addWidget(self.data_preview_label)

        file_group.setLayout(file_layout)
        layout.addWidget(file_group)

        # Main content area with splitter
        content_splitter = QSplitter(Qt.Horizontal)

        # Left panel: Conditions and Parameters
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        left_layout.setContentsMargins(0, 0, 0, 0)

        # Condition selection group
        condition_group = QGroupBox("Conditions to Analyze")
        condition_layout = QVBoxLayout()

        self.condition_list = QListWidget()
        self.condition_list.setSelectionMode(QAbstractItemView.MultiSelection)
        self.condition_list.setMaximumHeight(120)
        condition_layout.addWidget(self.condition_list)

        cond_btn_row = QHBoxLayout()
        self.select_all_btn = QPushButton("Select All")
        self.select_all_btn.setAutoDefault(False)
        self.select_all_btn.clicked.connect(self._select_all_conditions)
        self.deselect_all_btn = QPushButton("Deselect All")
        self.deselect_all_btn.setAutoDefault(False)
        self.deselect_all_btn.clicked.connect(self._deselect_all_conditions)
        cond_btn_row.addWidget(self.select_all_btn)
        cond_btn_row.addWidget(self.deselect_all_btn)
        condition_layout.addLayout(cond_btn_row)

        condition_group.setLayout(condition_layout)
        left_layout.addWidget(condition_group)

        # Method selection group
        method_group = QGroupBox("Integration Method")
        method_layout = QVBoxLayout()

        method_row = QHBoxLayout()
        method_row.addWidget(QLabel("Method:"))
        self.method_combo = no_scroll(QComboBox())
        self.method_combo.addItems(["LAD (Least Absolute Deviation)", "E-Flux2 (Expression-based FBA)"])
        self.method_combo.setToolTip(
            "LAD: Minimizes absolute deviation between fluxes and expression targets.\n"
            "E-Flux2: Constrains upper bounds based on expression, then optimizes growth."
        )
        self.method_combo.currentTextChanged.connect(self._on_method_changed)
        method_row.addWidget(self.method_combo)
        method_layout.addLayout(method_row)

        self.method_desc_label = QLabel(
            "LAD minimizes the deviation between predicted fluxes and expression-derived targets."
        )
        self.method_desc_label.setWordWrap(True)
        self.method_desc_label.setStyleSheet("color: gray; font-style: italic;")
        method_layout.addWidget(self.method_desc_label)

        method_group.setLayout(method_layout)
        left_layout.addWidget(method_group)

        # Parameters group
        params_group = QGroupBox("Analysis Parameters")
        params_layout = QVBoxLayout()

        agg_row = QHBoxLayout()
        agg_row.addWidget(QLabel("Gene aggregation:"))
        agg_label = QLabel("GPR-based (AND→min, OR→max)")
        agg_label.setToolTip(
            "Gene-Protein-Reaction rules are evaluated recursively:\n"
            "AND (enzyme complex) → min of subunit expression\n"
            "OR (isozymes) → max of alternative expression"
        )
        agg_label.setStyleSheet("font-weight: bold;")
        agg_row.addWidget(agg_label)
        params_layout.addLayout(agg_row)

        threshold_row = QHBoxLayout()
        threshold_row.addWidget(QLabel("Weight threshold:"))
        self.threshold_spin = no_scroll(QDoubleSpinBox())
        self.threshold_spin.setDecimals(4)
        self.threshold_spin.setRange(0, 1000)
        self.threshold_spin.setValue(0.01)
        threshold_row.addWidget(self.threshold_spin)
        params_layout.addLayout(threshold_row)

        scale_row = QHBoxLayout()
        scale_row.addWidget(QLabel("Scaling factor:"))
        self.scale_spin = no_scroll(QDoubleSpinBox())
        self.scale_spin.setDecimals(4)
        self.scale_spin.setRange(0.0001, 10000)
        self.scale_spin.setValue(1.0)
        scale_row.addWidget(self.scale_spin)
        params_layout.addLayout(scale_row)

        self.use_scenario_check = QCheckBox("Apply scenario constraints")
        self.use_scenario_check.setChecked(True)
        params_layout.addWidget(self.use_scenario_check)

        # E-Flux2 specific parameters
        self.eflux2_params_widget = QGroupBox("E-Flux2 Options")
        eflux2_layout = QVBoxLayout()

        min_scale_row = QHBoxLayout()
        min_scale_row.addWidget(QLabel("Minimum scale:"))
        self.min_scale_spin = no_scroll(QDoubleSpinBox())
        self.min_scale_spin.setDecimals(4)
        self.min_scale_spin.setRange(0.0001, 1.0)
        self.min_scale_spin.setValue(0.001)
        self.min_scale_spin.setToolTip(
            "Minimum scaling factor for numerical stability (0.001 = 0.1%).\n"
            "Prevents reaction bounds from becoming too small."
        )
        min_scale_row.addWidget(self.min_scale_spin)
        eflux2_layout.addLayout(min_scale_row)

        obj_frac_row = QHBoxLayout()
        obj_frac_row.addWidget(QLabel("Objective fraction:"))
        self.obj_fraction_spin = no_scroll(QDoubleSpinBox())
        self.obj_fraction_spin.setDecimals(2)
        self.obj_fraction_spin.setRange(0.0, 1.0)
        self.obj_fraction_spin.setValue(0.99)
        self.obj_fraction_spin.setToolTip(
            "Fraction of optimal objective to maintain during L2 minimization.\n"
            "Default 0.99 means 99% of maximum growth rate is preserved."
        )
        obj_frac_row.addWidget(self.obj_fraction_spin)
        eflux2_layout.addLayout(obj_frac_row)

        norm_pctl_row = QHBoxLayout()
        norm_pctl_row.addWidget(QLabel("Normalization percentile:"))
        self.norm_percentile_spin = no_scroll(QDoubleSpinBox())
        self.norm_percentile_spin.setDecimals(1)
        self.norm_percentile_spin.setRange(50.0, 100.0)
        self.norm_percentile_spin.setSingleStep(1.0)
        self.norm_percentile_spin.setValue(100.0)
        self.norm_percentile_spin.setToolTip(
            "Percentile of above-threshold weights used as the normalization denominator.\n"
            "100 = max normalization (paper-faithful).\n"
            "99 = recommended robust mode for real RNA-seq data — reduces sensitivity\n"
            "      to extreme expression outliers that would otherwise crush every\n"
            "      other reaction's bound.\n"
            "95 = aggressive robustness, useful for sensitivity analysis."
        )
        norm_pctl_row.addWidget(self.norm_percentile_spin)
        eflux2_layout.addLayout(norm_pctl_row)

        self.eflux2_params_widget.setLayout(eflux2_layout)
        self.eflux2_params_widget.setVisible(False)
        params_layout.addWidget(self.eflux2_params_widget)

        params_group.setLayout(params_layout)
        left_layout.addWidget(params_group)

        left_layout.addStretch()
        content_splitter.addWidget(left_panel)

        # Right panel: Results
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        right_layout.setContentsMargins(0, 0, 0, 0)

        results_group = QGroupBox("Analysis Results")
        results_layout = QVBoxLayout()

        # FC (Fold Change) reference condition selector
        fc_row = QHBoxLayout()
        fc_row.addWidget(QLabel("FC Reference:"))
        self.fc_reference_combo = no_scroll(QComboBox())
        self.fc_reference_combo.setToolTip(
            "Select reference condition for Fold Change calculation.\n"
            "FC = log2(condition / reference)"
        )
        self.fc_reference_combo.currentTextChanged.connect(self._on_fc_reference_changed)
        fc_row.addWidget(self.fc_reference_combo)
        fc_row.addStretch()
        results_layout.addLayout(fc_row)

        self.results_tabs = QTabWidget()

        # Comparison table tab
        self.comparison_table = QTableWidget()
        self.comparison_table.setAlternatingRowColors(True)
        self.comparison_table.setSortingEnabled(True)  # Enable column sorting
        self.results_tabs.addTab(self.comparison_table, "Comparison")

        results_layout.addWidget(self.results_tabs)

        # Result actions
        result_actions = QHBoxLayout()
        self.apply_btn = QPushButton("Apply Selected to Main View")
        self.apply_btn.setAutoDefault(False)
        self.apply_btn.setToolTip("Apply the selected condition's flux values to the main view")
        self.apply_btn.clicked.connect(self._apply_selected_to_main)
        self.apply_btn.setEnabled(False)
        result_actions.addWidget(self.apply_btn)

        self.export_btn = QPushButton("Export Results...")
        self.export_btn.setAutoDefault(False)
        self.export_btn.setToolTip("Export all results to CSV file")
        self.export_btn.clicked.connect(self._export_results)
        self.export_btn.setEnabled(False)
        result_actions.addWidget(self.export_btn)

        results_layout.addLayout(result_actions)

        results_group.setLayout(results_layout)
        right_layout.addWidget(results_group)

        content_splitter.addWidget(right_panel)
        content_splitter.setSizes([350, 550])

        layout.addWidget(content_splitter)

        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)

        # Bottom buttons
        btn_layout = QHBoxLayout()

        self.run_btn = QPushButton("Run Analysis")
        self.run_btn.setAutoDefault(False)
        self.run_btn.clicked.connect(self._run_multi_analysis)
        self.run_btn.setEnabled(False)
        btn_layout.addWidget(self.run_btn)

        self.close_btn = QPushButton("Close")
        self.close_btn.setAutoDefault(False)
        self.close_btn.clicked.connect(self.reject)
        btn_layout.addWidget(self.close_btn)

        layout.addLayout(btn_layout)

        self.setLayout(layout)
        self._on_method_changed(self.method_combo.currentText())

    def _on_method_changed(self, method_text: str):
        """Handle method selection change."""
        is_eflux2 = "E-Flux2" in method_text
        self.eflux2_params_widget.setVisible(is_eflux2)

        if is_eflux2:
            self.method_desc_label.setText(
                "E-Flux2 constrains reaction bounds based on gene expression levels\n"
                "(excluding exchange/biomass reactions), then performs FBA followed by\n"
                "L2-norm minimization (QP) for a unique flux solution.\n"
                "Falls back to pFBA (L1) if QP is unavailable."
            )
            self.run_btn.setText("Run E-Flux2 Analysis")
        else:
            self.method_desc_label.setText(
                "LAD minimizes the deviation between predicted fluxes and expression targets."
            )
            self.run_btn.setText("Run LAD Analysis")

    def _browse_file(self):
        """Open file dialog to select expression data file."""
        filename, _ = QFileDialog.getOpenFileName(
            self,
            "Select Expression Data File",
            self.appdata.work_directory,
            "Data files (*.csv *.tsv *.txt *.xlsx *.xls);;All files (*.*)",
        )

        if filename:
            self.file_edit.setText(filename)
            self._load_expression_data(filename)

    def _load_expression_data(self, filename: str):
        """Load expression data from the selected file (auto-detect single or multi-condition)."""
        try:
            # Try to load as multi-condition first
            self.conditions, self.multi_expression_data = read_multi_condition_expression_data(filename)

            # Check if it's really multi-condition or single-condition
            if len(self.conditions) == 1:
                # Single condition - convert to legacy format for backward compatibility
                self.is_multi_condition = False
                self.expression_data = {
                    gene: values[self.conditions[0]]
                    for gene, values in self.multi_expression_data.items()
                }
            else:
                self.is_multi_condition = True
                # Also create legacy single-condition data using first condition
                if self.conditions:
                    self.expression_data = {
                        gene: values.get(self.conditions[0], 0.0)
                        for gene, values in self.multi_expression_data.items()
                    }

            # Update condition list widget
            self.condition_list.clear()
            for cond in self.conditions:
                item = QListWidgetItem(cond)
                item.setSelected(True)
                self.condition_list.addItem(item)

            # Update info labels
            self.gene_count_label.setText(f"Genes loaded: {len(self.multi_expression_data)}")
            self.conditions_label.setText(f"Conditions: {len(self.conditions)}")

            # Count genes that match model genes
            model_genes = {g.id for g in self.appdata.project.cobra_py_model.genes}
            matched = sum(1 for g in self.multi_expression_data if g in model_genes)
            self.matched_count_label.setText(f"Matched to model: {matched}")

            # Show data preview
            if len(self.multi_expression_data) > 0:
                preview_html = (
                    f"<b>Data Preview:</b><br>"
                    f"Total genes: {len(self.multi_expression_data)}<br>"
                    f"Matched to model: {matched}<br>"
                    f"Conditions: {', '.join(self.conditions[:5])}"
                    f"{'...' if len(self.conditions) > 5 else ''}<br><br>"
                    f"<b>Sample data (first 3 genes):</b><br>"
                )
                for i, (gene_id, values) in enumerate(list(self.multi_expression_data.items())[:3]):
                    match_status = "+" if gene_id in model_genes else "-"
                    val_str = ", ".join(f"{c}:{values.get(c, 0):.2f}" for c in self.conditions[:3])
                    preview_html += f"&nbsp;&nbsp;{match_status} {gene_id}: {val_str}<br>"

                self.data_preview_label.setText(preview_html)
                self.data_preview_label.setVisible(True)
            else:
                self.data_preview_label.setVisible(False)

            # Enable run button if we have data
            self.run_btn.setEnabled(len(self.multi_expression_data) > 0 and matched > 0)

            if len(self.multi_expression_data) == 0:
                QMessageBox.warning(self, "No data", "No expression data could be loaded from the file.")
            elif matched == 0:
                QMessageBox.warning(
                    self,
                    "No matches",
                    "No genes from the expression data match genes in the model.\n"
                    "Check that gene IDs use the same format.",
                )

        except Exception as e:
            QMessageBox.critical(self, "Load error", f"Failed to load expression data:\n{str(e)}")
            self.multi_expression_data = {}
            self.conditions = []
            self.gene_count_label.setText("Genes loaded: 0")
            self.matched_count_label.setText("Matched to model: 0")
            self.conditions_label.setText("Conditions: 0")
            self.data_preview_label.setVisible(False)
            self.run_btn.setEnabled(False)

    def _select_all_conditions(self):
        """Select all conditions in the list."""
        for i in range(self.condition_list.count()):
            self.condition_list.item(i).setSelected(True)

    def _deselect_all_conditions(self):
        """Deselect all conditions in the list."""
        for i in range(self.condition_list.count()):
            self.condition_list.item(i).setSelected(False)

    def _is_eflux2(self) -> bool:
        """Check if E-Flux2 method is selected."""
        return "E-Flux2" in self.method_combo.currentText()

    def _run_multi_analysis(self):
        """Run analysis for all selected conditions."""
        selected_conditions = [
            self.condition_list.item(i).text()
            for i in range(self.condition_list.count())
            if self.condition_list.item(i).isSelected()
        ]

        if not selected_conditions:
            QMessageBox.warning(self, "No conditions", "Please select at least one condition to analyze.")
            return

        if not self.multi_expression_data:
            QMessageBox.warning(self, "No data", "Please load expression data first.")
            return

        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, len(selected_conditions))
        self.run_btn.setEnabled(False)

        is_eflux2 = self._is_eflux2()
        method_name = "E-Flux2" if is_eflux2 else "LAD"
        threshold = self.threshold_spin.value()

        # Get scenario constraints if enabled
        flux_constraints = {}
        if self.use_scenario_check.isChecked():
            for rid, (lb, ub) in self.appdata.project.scen_values.items():
                flux_constraints[rid] = (lb, ub)

        self.multi_results.clear()
        failed_conditions = []
        successful_conditions = []
        pfba_fallback_conditions: list[str] = []  # conditions where pFBA was used instead of QP
        fba_fallback_conditions: list[str] = []   # conditions where plain FBA was used

        for i, condition in enumerate(selected_conditions):
            self.progress_bar.setValue(i)
            self.progress_bar.setFormat(f"Analyzing {condition}... (%p%)")

            try:
                # Extract expression data for this condition
                condition_expression = {
                    gene: values[condition]
                    for gene, values in self.multi_expression_data.items()
                    if condition in values
                }

                # Compute reaction weights
                reaction_weights = gene_expression_to_reaction_weights(
                    self.appdata.project.cobra_py_model, condition_expression
                )

                if not reaction_weights:
                    failed_conditions.append((condition, "No reaction weights computed"))
                    continue

                # Run analysis
                model = self.appdata.project.cobra_py_model.copy()

                if is_eflux2:
                    min_scale = self.min_scale_spin.value()
                    obj_fraction = self.obj_fraction_spin.value()
                    norm_percentile = self.norm_percentile_spin.value()
                    status, obj_value, flux_dist, l2_method = run_eflux2(
                        model,
                        reaction_weights,
                        flux_constraints,
                        weight_threshold=threshold,
                        min_scale=min_scale,
                        objective_fraction=obj_fraction,
                        normalization_percentile=norm_percentile,
                    )
                    if l2_method == "pfba":
                        pfba_fallback_conditions.append(condition)
                    elif l2_method == "fba":
                        fba_fallback_conditions.append(condition)
                else:
                    scaling = self.scale_spin.value()
                    status, obj_value, flux_dist = run_lad_fitting(
                        model, reaction_weights, flux_constraints, threshold, scaling
                    )

                if status == "optimal" and flux_dist:
                    self.multi_results[condition] = flux_dist
                    successful_conditions.append(condition)
                else:
                    failed_conditions.append((condition, f"Optimization status: {status}"))

            except Exception as e:
                failed_conditions.append((condition, str(e)))

        self.progress_bar.setValue(len(selected_conditions))
        self.progress_bar.setVisible(False)
        self.run_btn.setEnabled(True)

        # Update results display
        self._update_results_display()

        # Store results in appdata for main view access
        self.appdata.project.omics_results = self.multi_results.copy()
        self.appdata.project.omics_conditions = list(self.multi_results.keys())

        # Update central widget if available
        if hasattr(self.parent(), "centralWidget"):
            central = self.parent().centralWidget()
            if hasattr(central, "update_omics_selector"):
                central.update_omics_selector()

        # Show summary
        msg = f"{method_name} analysis completed.\n\n"
        if is_eflux2:
            p = self.norm_percentile_spin.value()
            note = " (paper-faithful max)" if p >= 100.0 else " (robust — outlier-resistant)"
            msg += f"Normalization percentile: {p:.1f}{note}\n\n"
        msg += f"Successful: {len(successful_conditions)} conditions\n"
        if successful_conditions:
            msg += f"  ({', '.join(successful_conditions[:5])}"
            if len(successful_conditions) > 5:
                msg += f"... +{len(successful_conditions)-5} more"
            msg += ")\n"

        if failed_conditions:
            msg += f"\nFailed: {len(failed_conditions)} conditions\n"
            for cond, reason in failed_conditions[:3]:
                msg += f"  - {cond}: {reason}\n"
            if len(failed_conditions) > 3:
                msg += f"  ... and {len(failed_conditions)-3} more\n"

        # Warn about L2 fallback methods
        if pfba_fallback_conditions:
            msg += (
                f"\nWarning: {len(pfba_fallback_conditions)} condition(s) used pFBA (L1) "
                f"instead of QP (L2) — solver does not support quadratic objectives.\n"
            )
        if fba_fallback_conditions:
            msg += (
                f"\nWarning: {len(fba_fallback_conditions)} condition(s) fell back to "
                f"plain FBA (no L2/L1 minimization).\n"
            )

        if successful_conditions:
            msg += "\nResults are available in the comparison table."
            QMessageBox.information(self, "Analysis Complete", msg)
        else:
            QMessageBox.warning(self, "Analysis Failed", msg)

    def _update_results_display(self):
        """Update the comparison table with multi-condition results."""
        if not self.multi_results:
            self.comparison_table.setRowCount(0)
            self.comparison_table.setColumnCount(0)
            self.fc_reference_combo.clear()
            self.apply_btn.setEnabled(False)
            self.export_btn.setEnabled(False)
            return

        conditions = list(self.multi_results.keys())

        # Update FC reference combo
        current_ref = self.fc_reference_combo.currentText()
        self.fc_reference_combo.blockSignals(True)
        self.fc_reference_combo.clear()
        self.fc_reference_combo.addItem("(None)")
        self.fc_reference_combo.addItems(conditions)
        # Restore previous selection if still valid
        if current_ref and current_ref in conditions:
            self.fc_reference_combo.setCurrentText(current_ref)
        elif conditions:
            self.fc_reference_combo.setCurrentIndex(1)  # Select first condition as default
        self.fc_reference_combo.blockSignals(False)

        # Collect all reactions
        all_reactions = set()
        for flux_dist in self.multi_results.values():
            all_reactions.update(flux_dist.keys())
        reactions = sorted(all_reactions)

        # Determine FC reference condition
        fc_ref = self.fc_reference_combo.currentText()
        show_fc = fc_ref and fc_ref != "(None)" and fc_ref in self.multi_results
        other_conditions = [c for c in conditions if c != fc_ref] if show_fc else []

        # Calculate column count: Reaction + conditions + FC columns
        num_fc_cols = len(other_conditions) if show_fc else 0
        total_cols = 1 + len(conditions) + num_fc_cols

        # Build header labels
        headers = ["Reaction"] + conditions
        if show_fc:
            headers += [f"FC({c}/{fc_ref})" for c in other_conditions]

        # Disable sorting while updating to prevent issues
        self.comparison_table.setSortingEnabled(False)

        # Setup table
        self.comparison_table.setRowCount(len(reactions))
        self.comparison_table.setColumnCount(total_cols)
        self.comparison_table.setHorizontalHeaderLabels(headers)

        # Fill table
        for row, rxn_id in enumerate(reactions):
            # Reaction ID
            item = QTableWidgetItem(rxn_id)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            self.comparison_table.setItem(row, 0, item)

            # Flux values for each condition
            for col, cond in enumerate(conditions, 1):
                flux = self.multi_results[cond].get(rxn_id, 0.0)
                item = NumericTableWidgetItem(f"{flux:.4f}", flux)
                item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.comparison_table.setItem(row, col, item)

            # FC columns
            if show_fc:
                ref_flux = self.multi_results[fc_ref].get(rxn_id, 0.0)
                for fc_col, cond in enumerate(other_conditions):
                    cond_flux = self.multi_results[cond].get(rxn_id, 0.0)
                    fc_value = self._calculate_fc(cond_flux, ref_flux)
                    col_idx = 1 + len(conditions) + fc_col

                    if fc_value is not None and not np.isinf(fc_value):
                        item = NumericTableWidgetItem(f"{fc_value:.3f}", fc_value)
                        # Color code: red for down, green for up
                        if fc_value < -1:
                            item.setBackground(QColor(255, 200, 200))  # Light red
                        elif fc_value > 1:
                            item.setBackground(QColor(200, 255, 200))  # Light green
                    else:
                        item = NumericTableWidgetItem("N/A", float('nan'))

                    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                    self.comparison_table.setItem(row, col_idx, item)

        # Re-enable sorting
        self.comparison_table.setSortingEnabled(True)

        # Resize columns
        self.comparison_table.resizeColumnsToContents()

        # Enable action buttons
        self.apply_btn.setEnabled(True)
        self.export_btn.setEnabled(True)

    def _calculate_fc(self, value: float, reference: float) -> float | None:
        """Calculate log2 fold change.

        Returns None if calculation is not possible (division by zero or negative values).
        """
        # Handle edge cases
        if reference == 0 or value < 0 or reference < 0:
            return None
        if value == 0:
            return float('-inf') if reference > 0 else None

        try:
            fc = np.log2(value / reference)
            return fc
        except (ValueError, ZeroDivisionError):
            return None

    @Slot(str)
    def _on_fc_reference_changed(self, _reference: str):
        """Handle FC reference condition change."""
        if self.multi_results:
            self._update_results_display()

    def _apply_selected_to_main(self):
        """Apply the first selected condition's results to the main view."""
        if not self.multi_results:
            return

        conditions = list(self.multi_results.keys())
        if not conditions:
            return

        # Show condition selection dialog
        from qtpy.QtWidgets import QInputDialog

        condition, ok = QInputDialog.getItem(
            self,
            "Select Condition",
            "Choose condition to apply to main view:",
            conditions,
            0,
            False,
        )

        if not ok or condition not in self.multi_results:
            return

        # Apply to comp_values
        self.appdata.project.comp_values.clear()
        for rxn_id, flux in self.multi_results[condition].items():
            self.appdata.project.comp_values[rxn_id] = (flux, flux)

        # Update central widget
        if hasattr(self.parent(), "centralWidget"):
            central = self.parent().centralWidget()
            central.update()
            if hasattr(central, "set_current_omics_condition"):
                central.set_current_omics_condition(condition)

        QMessageBox.information(
            self,
            "Applied",
            f"Flux values from condition '{condition}' have been applied to the main view.",
        )

    def _export_results(self):
        """Export all results to a CSV file."""
        if not self.multi_results:
            return

        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Export Results",
            self.appdata.work_directory,
            "CSV files (*.csv);;All files (*.*)",
        )

        if not filename:
            return

        if not filename.endswith(".csv"):
            filename += ".csv"

        try:
            conditions = list(self.multi_results.keys())
            all_reactions = set()
            for flux_dist in self.multi_results.values():
                all_reactions.update(flux_dist.keys())
            reactions = sorted(all_reactions)

            # Create DataFrame
            data = {"Reaction": reactions}
            for cond in conditions:
                data[cond] = [self.multi_results[cond].get(rxn, 0.0) for rxn in reactions]

            df = pd.DataFrame(data)
            df.to_csv(filename, index=False)

            QMessageBox.information(
                self,
                "Export Complete",
                f"Results exported to:\n{filename}\n\n"
                f"Exported {len(reactions)} reactions across {len(conditions)} conditions.",
            )

        except Exception as e:
            QMessageBox.critical(self, "Export Error", f"Failed to export results:\n{str(e)}")
