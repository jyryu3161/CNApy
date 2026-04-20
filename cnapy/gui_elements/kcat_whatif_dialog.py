"""
kcat_whatif_dialog.py
=====================
kcat What-if / Sensitivity Analysis dialog.

Lets the user:
  1. Select an enzyme (by UniProt ID or gene name)
  2. Enter a kcat multiplier (e.g. 2× = kcat doubled = 50% less protein needed)
  3. Run FBA on the modified model (non-destructive: uses a copy)
  4. Compare original vs. modified objective value and pool usage

Equivalent to GECKO's sensitivityTuning / kcat sensitivity analysis.
"""

from __future__ import annotations

from qtpy.QtCore import Qt, Slot
from qtpy.QtWidgets import (
    QDialog,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QHeaderView,
    QAbstractItemView,
)

from cnapy.appdata import AppData
from cnapy.ecmodel.ecmodel_builder import (
    apply_kcat_multiplier,
    get_enzyme_usage,
    run_fba_on_ecmodel,
)
from cnapy.utils import no_scroll


class KcatWhatifDialog(QDialog):
    """kcat sensitivity / what-if analysis dialog."""

    def __init__(self, appdata: AppData, parent=None):
        super().__init__(parent)
        self.appdata = appdata
        self.setWindowTitle("kcat What-if Analysis")
        self.setMinimumWidth(750)
        self.setMinimumHeight(560)
        self._baseline: dict | None = None   # baseline FBA result
        self._build_ui()
        self._run_baseline()

    # ── UI ────────────────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)

        info = QLabel(
            "<b>kcat What-if Analysis</b><br>"
            "Simulates the effect of engineering an enzyme's kcat on the objective "
            "(e.g. growth rate). A higher kcat means less protein needed per unit flux, "
            "so the same protein budget can support more growth."
        )
        info.setWordWrap(True)
        info.setTextFormat(Qt.RichText)
        info.setStyleSheet(
            "QLabel { background:#f0fff4; padding:8px; "
            "border:1px solid #cdc; border-radius:4px; }"
        )
        root.addWidget(info)

        # Baseline group
        baseline_group = QGroupBox("Baseline (current ecModel)")
        baseline_layout = QHBoxLayout()
        self.baseline_obj_label  = QLabel("Objective: —")
        self.baseline_pool_label = QLabel("Pool used: —  /  UB: —")
        baseline_layout.addWidget(self.baseline_obj_label)
        baseline_layout.addWidget(self.baseline_pool_label)
        baseline_layout.addStretch()
        self.refresh_baseline_btn = QPushButton("Re-run Baseline")
        self.refresh_baseline_btn.setAutoDefault(False)
        self.refresh_baseline_btn.clicked.connect(self._run_baseline)
        baseline_layout.addWidget(self.refresh_baseline_btn)
        baseline_group.setLayout(baseline_layout)
        root.addWidget(baseline_group)

        # What-if settings
        settings_group = QGroupBox("What-if Settings")
        form = QFormLayout()

        self.enzyme_edit = QLineEdit()
        self.enzyme_edit.setPlaceholderText("UniProt ID (e.g. P0A9P0) or gene name (e.g. thrA)")
        form.addRow("Enzyme:", self.enzyme_edit)

        self.multiplier_spin = QDoubleSpinBox()
        self.multiplier_spin.setRange(0.01, 1000.0)
        self.multiplier_spin.setSingleStep(0.5)
        self.multiplier_spin.setDecimals(2)
        self.multiplier_spin.setValue(2.0)
        self.multiplier_spin.setToolTip(
            "2× = kcat doubled (enzyme twice as efficient, uses half the protein)\n"
            "0.5× = kcat halved (enzyme less efficient, needs twice the protein)"
        )
        no_scroll(self.multiplier_spin)
        form.addRow("kcat multiplier (×):", self.multiplier_spin)

        settings_group.setLayout(form)
        root.addWidget(settings_group)

        # Run button
        self.run_btn = QPushButton("Run What-if FBA")
        self.run_btn.setAutoDefault(False)
        self.run_btn.clicked.connect(self._run_whatif)
        root.addWidget(self.run_btn)

        # Result group
        result_group = QGroupBox("Result Comparison")
        result_layout = QVBoxLayout()

        comp_row = QHBoxLayout()
        self.result_obj_label  = QLabel("Modified objective: —")
        self.result_diff_label = QLabel("Change: —")
        self.result_pool_label = QLabel("Pool used: —")
        comp_row.addWidget(self.result_obj_label)
        comp_row.addWidget(self.result_diff_label)
        comp_row.addWidget(self.result_pool_label)
        comp_row.addStretch()
        result_layout.addLayout(comp_row)

        result_layout.addWidget(QLabel("Top 10 enzyme usage changes:"))
        self.diff_table = QTableWidget(0, 4)
        self.diff_table.setHorizontalHeaderLabels([
            "UniProt ID", "Baseline cap. (%)", "Modified cap. (%)", "Δ cap. (%)"
        ])
        self.diff_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.diff_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.diff_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.diff_table.setMaximumHeight(200)
        result_layout.addWidget(self.diff_table)

        result_group.setLayout(result_layout)
        root.addWidget(result_group)

        # Close
        btn_row = QHBoxLayout()
        btn_row.addStretch()
        close_btn = QPushButton("Close")
        close_btn.setAutoDefault(False)
        close_btn.clicked.connect(self.close)
        btn_row.addWidget(close_btn)
        root.addLayout(btn_row)

    # ── helpers ───────────────────────────────────────────────────────────────

    def _resolve_uniprot(self, text: str) -> str | None:
        """Resolve a UniProt ID or gene name to a UniProt ID."""
        text = text.strip()
        if not text:
            return None
        ec = self.appdata.project.ec_model_data
        # Direct UniProt ID
        if text in ec.uniprot_data:
            return text
        # Gene name lookup
        for uid, info in ec.uniprot_data.items():
            gene = info.get("gene", "")
            if gene and text.lower() in gene.lower().split():
                return uid
        # Check kcat_entries gene_name field
        for entry in ec.kcat_entries:
            if text.lower() == entry.get("gene_name", "").lower():
                return entry["proteins"][0] if entry["proteins"] else None
        return None

    # ── slots ─────────────────────────────────────────────────────────────────

    @Slot()
    def _run_baseline(self):
        ec = self.appdata.project.ec_model_data
        if not ec.is_ecmodel:
            self.baseline_obj_label.setText("Objective: (not an ecModel)")
            return

        ecmodel = self.appdata.project.cobra_py_model
        self.setCursor(Qt.BusyCursor)
        try:
            sol = run_fba_on_ecmodel(ecmodel)
        except Exception as exc:
            self.setCursor(Qt.ArrowCursor)
            self.baseline_obj_label.setText(f"FBA failed: {exc}")
            return
        self.setCursor(Qt.ArrowCursor)

        if sol.status != "optimal":
            self.baseline_obj_label.setText(f"Status: {sol.status}")
            return

        self._baseline = {
            "objective": sol.objective_value,
            "usage": get_enzyme_usage(ecmodel, sol),
            "solution": sol,
        }

        pool_rxn_id = ec.protein_pool_rxn_id
        try:
            pool_ub = ecmodel.reactions.get_by_id(pool_rxn_id).upper_bound
            pool_used = abs(sol.fluxes.get(pool_rxn_id, 0.0))
        except KeyError:
            pool_ub = pool_used = 0.0

        self.baseline_obj_label.setText(f"Objective: {sol.objective_value:.6f}")
        self.baseline_pool_label.setText(
            f"Pool used: {pool_used:.4f}  /  UB: {pool_ub:.4f} mg/gDCW"
        )

    @Slot()
    def _run_whatif(self):
        ec = self.appdata.project.ec_model_data
        if not ec.is_ecmodel:
            QMessageBox.warning(self, "Not an ecModel", "Convert to ecModel first.")
            return
        if self._baseline is None:
            QMessageBox.warning(self, "No baseline", "Run baseline FBA first.")
            return

        uid = self._resolve_uniprot(self.enzyme_edit.text())
        if uid is None:
            QMessageBox.warning(self, "Enzyme not found",
                                f"Could not find '{self.enzyme_edit.text()}' "
                                f"in UniProt data or kcat entries.")
            return

        multiplier = self.multiplier_spin.value()
        ecmodel = self.appdata.project.cobra_py_model

        # Work on a copy — never modify the live model
        self.setCursor(Qt.BusyCursor)
        try:
            modified = ecmodel.copy()
            n_rxns = apply_kcat_multiplier(modified, ec, uid, multiplier)
            if n_rxns == 0:
                self.setCursor(Qt.ArrowCursor)
                QMessageBox.warning(self, "No reactions updated",
                                    f"Enzyme {uid} has no constrained reactions in the ecModel.")
                return
            sol = run_fba_on_ecmodel(modified)
        except Exception as exc:
            self.setCursor(Qt.ArrowCursor)
            QMessageBox.critical(self, "What-if FBA failed", str(exc))
            return
        self.setCursor(Qt.ArrowCursor)

        if sol.status != "optimal":
            self.result_obj_label.setText(f"Status: {sol.status}")
            return

        new_obj = sol.objective_value
        old_obj = self._baseline["objective"]
        diff = new_obj - old_obj
        pct  = (diff / old_obj * 100) if old_obj != 0 else 0.0

        self.result_obj_label.setText(f"Modified objective: {new_obj:.6f}")
        sign = "+" if diff >= 0 else ""
        color = "green" if diff >= 0 else "red"
        self.result_diff_label.setText(
            f"<span style='color:{color}'>"
            f"Change: {sign}{diff:.6f} ({sign}{pct:.2f}%)</span>"
        )
        self.result_diff_label.setTextFormat(Qt.RichText)

        pool_rxn_id = ec.protein_pool_rxn_id
        try:
            pool_used = abs(sol.fluxes.get(pool_rxn_id, 0.0))
            self.result_pool_label.setText(f"Pool used: {pool_used:.4f}")
        except Exception:
            pass

        # Build diff table
        new_usage = get_enzyme_usage(modified, sol)
        base_usage = self._baseline["usage"]

        diffs = []
        all_ids = set(base_usage) | set(new_usage)
        for u in all_ids:
            base_cap = base_usage.get(u, {}).get("cap_usage", 0.0) * 100
            new_cap  = new_usage.get(u, {}).get("cap_usage", 0.0) * 100
            delta    = new_cap - base_cap
            diffs.append((u, base_cap, new_cap, delta))

        diffs.sort(key=lambda x: abs(x[3]), reverse=True)
        top = diffs[:10]

        self.diff_table.setRowCount(len(top))
        for i, (u, base_c, new_c, delta) in enumerate(top):
            self.diff_table.setItem(i, 0, QTableWidgetItem(u))
            self.diff_table.setItem(i, 1, QTableWidgetItem(f"{base_c:.1f}%"))
            self.diff_table.setItem(i, 2, QTableWidgetItem(f"{new_c:.1f}%"))
            sign = "+" if delta >= 0 else ""
            delta_item = QTableWidgetItem(f"{sign}{delta:.1f}%")
            delta_item.setForeground(Qt.darkGreen if delta < 0 else Qt.darkRed)
            self.diff_table.setItem(i, 3, delta_item)
