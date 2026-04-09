"""
enzyme_usage_dialog.py
======================
Enzyme Usage Viewer — shows per-enzyme protein usage after FBA.

Equivalent to GECKO's reportEnzymeUsage / enzymeUsage MATLAB functions.

Displays:
  - Absolute usage (mg/gDCW) per enzyme
  - Capacity usage (%) = abs_usage / UB × 100
  - Bottleneck enzymes highlighted (capUsage ≥ 90%)
  - Catalysed reactions for each enzyme
"""

from __future__ import annotations

import cobra
from qtpy.QtCore import Qt, Slot
from qtpy.QtGui import QColor
from qtpy.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from cnapy.appdata import AppData
from cnapy.ecmodel.ecmodel_builder import get_enzyme_usage, run_fba_on_ecmodel


class _NumItem(QTableWidgetItem):
    """Sortable numeric table item."""
    def __init__(self, value: float, fmt: str = "{:.4f}"):
        super().__init__(fmt.format(value))
        self._value = value
        self.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)

    def __lt__(self, other):
        if isinstance(other, _NumItem):
            return self._value < other._value
        return super().__lt__(other)


class EnzymeUsageDialog(QDialog):
    """Dialog showing per-enzyme protein usage from an FBA solution."""

    def __init__(self, appdata: AppData, parent=None):
        super().__init__(parent)
        self.appdata = appdata
        self.setWindowTitle("Enzyme Usage Report")
        self.setMinimumWidth(900)
        self.setMinimumHeight(600)
        self._usage: dict = {}
        self._build_ui()

    # ── UI ────────────────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)

        # Summary bar
        summary_group = QGroupBox("Protein Pool Summary")
        summary_layout = QHBoxLayout()
        self.pool_ub_label    = QLabel("Pool UB: —")
        self.pool_used_label  = QLabel("Used: —")
        self.pool_pct_label   = QLabel("Pool usage: —")
        self.bottleneck_label = QLabel("Bottleneck enzymes (≥90%): —")
        for lbl in [self.pool_ub_label, self.pool_used_label,
                    self.pool_pct_label, self.bottleneck_label]:
            summary_layout.addWidget(lbl)
        summary_layout.addStretch()
        summary_group.setLayout(summary_layout)
        root.addWidget(summary_group)

        # Main table
        table_group = QGroupBox("Enzyme Usage (sorted by capacity usage)")
        table_layout = QVBoxLayout()

        self.table = QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels([
            "UniProt ID", "Gene", "Abs. Usage\n(mg/gDCW)",
            "UB (mg/gDCW)", "Cap. Usage (%)", "Catalysed Reactions"
        ])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(5, QHeaderView.Stretch)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSortingEnabled(True)
        table_layout.addWidget(self.table)

        legend_row = QHBoxLayout()
        red_lbl = QLabel("  ")
        red_lbl.setStyleSheet("background-color: #ffcccc;")
        legend_row.addWidget(red_lbl)
        legend_row.addWidget(QLabel("Bottleneck (≥90%)"))
        legend_row.addSpacing(12)
        yellow_lbl = QLabel("  ")
        yellow_lbl.setStyleSheet("background-color: #fff3cd;")
        legend_row.addWidget(yellow_lbl)
        legend_row.addWidget(QLabel("High usage (70–90%)"))
        legend_row.addStretch()
        table_layout.addLayout(legend_row)

        table_group.setLayout(table_layout)
        root.addWidget(table_group)

        # Buttons
        btn_row = QHBoxLayout()
        self.run_btn = QPushButton("Run FBA & Compute Usage")
        self.run_btn.setAutoDefault(False)
        self.run_btn.clicked.connect(self._run_fba)

        self.export_btn = QPushButton("Export CSV…")
        self.export_btn.setAutoDefault(False)
        self.export_btn.setEnabled(False)
        self.export_btn.clicked.connect(self._export_csv)

        close_btn = QPushButton("Close")
        close_btn.setAutoDefault(False)
        close_btn.clicked.connect(self.close)

        btn_row.addWidget(self.run_btn)
        btn_row.addWidget(self.export_btn)
        btn_row.addStretch()
        btn_row.addWidget(close_btn)
        root.addLayout(btn_row)

    # ── slots ─────────────────────────────────────────────────────────────────

    @Slot()
    def _run_fba(self):
        ec = self.appdata.project.ec_model_data
        if not ec.is_ecmodel:
            from qtpy.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Not an ecModel",
                                "Please convert the model to an ecModel first.")
            return

        ecmodel = self.appdata.project.cobra_py_model
        self.setCursor(Qt.BusyCursor)
        try:
            solution = run_fba_on_ecmodel(ecmodel)
        except Exception as exc:
            self.setCursor(Qt.ArrowCursor)
            from qtpy.QtWidgets import QMessageBox
            QMessageBox.critical(self, "FBA failed", str(exc))
            return
        self.setCursor(Qt.ArrowCursor)

        if solution.status != "optimal":
            from qtpy.QtWidgets import QMessageBox
            QMessageBox.warning(self, "No optimal solution",
                                f"FBA status: {solution.status}")
            return

        self._usage = get_enzyme_usage(ecmodel, solution)
        self._populate_table(ecmodel, solution)
        self.export_btn.setEnabled(True)

    def _populate_table(self, ecmodel: cobra.Model, solution: cobra.Solution):
        ec = self.appdata.project.ec_model_data

        # Pool summary
        pool_rxn_id = ec.protein_pool_rxn_id
        try:
            pool_rxn = ecmodel.reactions.get_by_id(pool_rxn_id)
            pool_ub = pool_rxn.upper_bound
            pool_used = abs(solution.fluxes.get(pool_rxn_id, 0.0))
            pool_pct = (pool_used / pool_ub * 100) if pool_ub > 0 else 0.0
            self.pool_ub_label.setText(f"Pool UB: {pool_ub:.2f} mg/gDCW")
            self.pool_used_label.setText(f"Used: {pool_used:.4f} mg/gDCW")
            self.pool_pct_label.setText(f"Pool usage: {pool_pct:.1f}%")
        except KeyError:
            pool_ub = 1.0

        # Build rows: uid → gene, reactions
        rows = []
        for uid, data in self._usage.items():
            gene = ec.uniprot_data.get(uid, {}).get("gene", "")
            # Find catalysed reactions (not usage itself, not REV unless relevant)
            met_id = f"prot_{uid}"
            try:
                met = ecmodel.metabolites.get_by_id(met_id)
                rxn_ids = [r.id for r in met.reactions
                           if not r.id.startswith("usage_prot_")
                           and abs(solution.fluxes.get(r.id, 0.0)) > 1e-9]
            except KeyError:
                rxn_ids = []
            rows.append((uid, gene, data["abs_usage"], data["ub"],
                         data["cap_usage"] * 100, ", ".join(rxn_ids)))

        # Sort by cap_usage descending
        rows.sort(key=lambda x: x[4], reverse=True)

        bottleneck_n = sum(1 for r in rows if r[4] >= 90)
        self.bottleneck_label.setText(f"Bottleneck enzymes (≥90%): {bottleneck_n}")

        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(rows))
        for i, (uid, gene, abs_u, ub, cap_pct, rxns) in enumerate(rows):
            self.table.setItem(i, 0, QTableWidgetItem(uid))
            self.table.setItem(i, 1, QTableWidgetItem(gene))
            self.table.setItem(i, 2, _NumItem(abs_u))
            self.table.setItem(i, 3, _NumItem(ub))
            self.table.setItem(i, 4, _NumItem(cap_pct, "{:.1f}%"))
            self.table.setItem(i, 5, QTableWidgetItem(rxns))

            # Colour coding
            if cap_pct >= 90:
                color = QColor("#ffcccc")
            elif cap_pct >= 70:
                color = QColor("#fff3cd")
            else:
                color = QColor("#ffffff")
            for col in range(6):
                item = self.table.item(i, col)
                if item:
                    item.setBackground(color)

        self.table.setSortingEnabled(True)

    @Slot()
    def _export_csv(self):
        if not self._usage:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Enzyme Usage",
            self.appdata.work_directory,
            "CSV files (*.csv)",
            options=QFileDialog.DontUseNativeDialog,
        )
        if not path:
            return
        if not path.endswith(".csv"):
            path += ".csv"
        try:
            with open(path, "w", encoding="utf-8") as f:
                f.write("UniProt ID,Gene,Abs Usage (mg/gDCW),UB (mg/gDCW),Cap Usage (%),Reactions\n")
                for row in range(self.table.rowCount()):
                    vals = [self.table.item(row, c).text() for c in range(6)]
                    f.write(",".join(f'"{v}"' for v in vals) + "\n")
        except Exception as exc:
            from qtpy.QtWidgets import QMessageBox
            QMessageBox.critical(self, "Export failed", str(exc))
