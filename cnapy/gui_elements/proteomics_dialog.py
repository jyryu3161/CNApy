"""
proteomics_dialog.py
====================
Proteomics Integration Dialog for GECKO ecModel.

Allows loading measured protein abundances (mg/gDCW) to individually
constrain enzyme usage reactions, replacing the shared protein pool
constraint for measured enzymes.

Equivalent to GECKO's fillEnzConcs + constrainEnzConcs.
"""

from __future__ import annotations

from qtpy.QtCore import Qt, Slot
from qtpy.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)

from cnapy.appdata import AppData
from cnapy.ecmodel.ecmodel_builder import (
    apply_proteomics,
    parse_proteomics_file,
    remove_proteomics,
)


class ProteomicsDialog(QDialog):
    """Dialog for loading proteomics data to constrain individual enzymes."""

    def __init__(self, appdata: AppData, parent=None):
        super().__init__(parent)
        self.appdata = appdata
        self.setWindowTitle("Proteomics Integration")
        self.setMinimumWidth(750)
        self.setMinimumHeight(520)
        self._prot_data: dict = {}
        self._applied = False
        self._build_ui()
        self._check_state()

    # ── UI ────────────────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)

        # Info box
        info = QLabel(
            "<b>Proteomics Integration</b><br>"
            "Measured protein abundances (mg/gDCW) are used to set individual "
            "upper bounds on each enzyme's usage reaction, replacing the shared "
            "protein pool constraint for that enzyme.<br>"
            "Enzymes without measured data continue to draw from the protein pool.<br>"
            "File format: two columns — <tt>UniProt_ID</tt> | <tt>level (mg/gDCW)</tt>"
        )
        info.setWordWrap(True)
        info.setTextFormat(Qt.RichText)
        info.setStyleSheet(
            "QLabel { background:#f0f4ff; padding:8px; "
            "border:1px solid #ccd; border-radius:4px; }"
        )
        root.addWidget(info)

        # File load group
        file_group = QGroupBox("Proteomics Data File")
        file_layout = QVBoxLayout()

        file_row = QHBoxLayout()
        self.file_edit = QLineEdit()
        self.file_edit.setPlaceholderText("proteomics.tsv / .csv / .xlsx — UniProt ID | level (mg/gDCW)")
        self.file_edit.setReadOnly(True)
        file_row.addWidget(self.file_edit)
        browse_btn = QPushButton("Browse…")
        browse_btn.setAutoDefault(False)
        browse_btn.clicked.connect(self._browse)
        file_row.addWidget(browse_btn)
        file_layout.addLayout(file_row)

        self.file_status = QLabel("Entries loaded: —")
        file_layout.addWidget(self.file_status)
        file_group.setLayout(file_layout)
        root.addWidget(file_group)

        # Preview table
        preview_group = QGroupBox("Preview (matched to ecModel)")
        preview_layout = QVBoxLayout()
        self.match_label = QLabel("Matched: — / —")
        preview_layout.addWidget(self.match_label)

        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(["UniProt ID", "Level (mg/gDCW)", "In ecModel?"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setMaximumHeight(220)
        preview_layout.addWidget(self.table)
        preview_group.setLayout(preview_layout)
        root.addWidget(preview_group)

        # Status
        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        root.addWidget(self.status_label)

        # Buttons
        btn_row = QHBoxLayout()
        self.apply_btn = QPushButton("Apply Proteomics Constraints")
        self.apply_btn.setAutoDefault(False)
        self.apply_btn.setEnabled(False)
        self.apply_btn.clicked.connect(self._apply)

        self.remove_btn = QPushButton("Remove Proteomics Constraints")
        self.remove_btn.setAutoDefault(False)
        self.remove_btn.setEnabled(False)
        self.remove_btn.clicked.connect(self._remove)

        close_btn = QPushButton("Close")
        close_btn.setAutoDefault(False)
        close_btn.clicked.connect(self.close)

        btn_row.addWidget(self.apply_btn)
        btn_row.addWidget(self.remove_btn)
        btn_row.addStretch()
        btn_row.addWidget(close_btn)
        root.addLayout(btn_row)

    # ── helpers ───────────────────────────────────────────────────────────────

    def _check_state(self):
        ec = self.appdata.project.ec_model_data
        if not ec.is_ecmodel:
            self.status_label.setText(
                "ecModel이 아닙니다. 먼저 GECKO 다이얼로그에서 Convert to ecModel을 실행하세요."
            )
            self.status_label.setStyleSheet("color: red;")
            self.apply_btn.setEnabled(False)

    def _populate_table(self):
        ec = self.appdata.project.ec_model_data
        ecmodel = self.appdata.project.cobra_py_model

        self.table.setRowCount(len(self._prot_data))
        matched = 0
        for i, (uid, level) in enumerate(sorted(self._prot_data.items())):
            in_model = f"usage_prot_{uid}" in {r.id for r in ecmodel.reactions}
            if in_model:
                matched += 1
            self.table.setItem(i, 0, QTableWidgetItem(uid))
            self.table.setItem(i, 1, QTableWidgetItem(f"{level:.6f}"))
            status_item = QTableWidgetItem("Yes" if in_model else "No")
            status_item.setForeground(Qt.darkGreen if in_model else Qt.red)
            self.table.setItem(i, 2, status_item)

        self.match_label.setText(
            f"Matched: {matched} / {len(self._prot_data)} enzymes found in ecModel"
        )
        self.apply_btn.setEnabled(matched > 0 and ec.is_ecmodel)

    # ── slots ─────────────────────────────────────────────────────────────────

    @Slot()
    def _browse(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open Proteomics File",
            self.appdata.work_directory,
            "Data files (*.tsv *.csv *.xlsx *.xls);;All files (*)",
            options=QFileDialog.DontUseNativeDialog,
        )
        if not path:
            return
        try:
            self._prot_data = parse_proteomics_file(path)
            self.file_edit.setText(path)
            self.file_status.setText(f"Entries loaded: {len(self._prot_data)}")
            self._populate_table()
        except Exception as exc:
            QMessageBox.critical(self, "Error reading file", str(exc))

    @Slot()
    def _apply(self):
        if not self._prot_data:
            return
        ecmodel = self.appdata.project.cobra_py_model
        ec = self.appdata.project.ec_model_data
        self.setCursor(Qt.BusyCursor)
        try:
            n_applied, missing = apply_proteomics(ecmodel, ec, self._prot_data)
        except Exception as exc:
            self.setCursor(Qt.ArrowCursor)
            QMessageBox.critical(self, "Error applying proteomics", str(exc))
            return
        self.setCursor(Qt.ArrowCursor)

        self._applied = True
        self.remove_btn.setEnabled(True)
        self.apply_btn.setEnabled(False)

        msg = f"Applied proteomics constraints to {n_applied} enzyme(s)."
        if missing:
            msg += f"\n{len(missing)} UniProt ID(s) not found in ecModel: {', '.join(missing[:5])}"
            if len(missing) > 5:
                msg += f" …(+{len(missing)-5})"
        self.status_label.setText(msg)
        self.status_label.setStyleSheet("color: green;")
        QMessageBox.information(self, "Proteomics applied", msg)

        # Refresh main window
        self.parent().centralWidget().update(rebuild_all_tabs=True)

    @Slot()
    def _remove(self):
        ret = QMessageBox.question(
            self, "Remove Proteomics Constraints",
            "This will reconnect all enzymes to the protein pool "
            "and restore default upper bounds.\nContinue?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No
        )
        if ret != QMessageBox.Yes:
            return
        ecmodel = self.appdata.project.cobra_py_model
        ec = self.appdata.project.ec_model_data
        n = remove_proteomics(ecmodel, ec)
        self._applied = False
        self.remove_btn.setEnabled(False)
        self.apply_btn.setEnabled(bool(self._prot_data))
        self.status_label.setText(f"Removed proteomics constraints from {n} enzyme(s).")
        self.status_label.setStyleSheet("color: blue;")
