"""
gecko_dialog.py
===============
Main GECKO enzyme-constrained model (ecModel) dialog for CNApy.

Allows the user to:
1. Load customKcats.tsv and uniprot.tsv files
2. Set sigma / f / Ptot parameters
3. Choose Full or Light formulation
4. Preview coverage (which reactions are constrained, which are missing)
5. Convert the current GEM to an ecModel (or revert)
6. Optionally launch the UniKP kcat prediction sub-dialog
"""

from __future__ import annotations

import os

from qtpy.QtCore import Qt, Slot
from qtpy.QtWidgets import (
    QAbstractItemView,
    QButtonGroup,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QRadioButton,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from cnapy.appdata import AppData
from cnapy.ecmodel.ecmodel_builder import (
    build_ecmodel,
    compute_coverage,
    parse_customkcats_file,
    parse_uniprot_file,
    revert_to_gem,
)
from cnapy.ecmodel.ecmodel_data import ECModelData
from cnapy.utils import no_scroll


class GeckoDialog(QDialog):
    """Dialog for building GECKO enzyme-constrained models in CNApy."""

    def __init__(self, appdata: AppData, parent=None):
        super().__init__(parent)
        self.appdata = appdata
        self.setWindowTitle("GECKO – Enzyme-Constrained Model Builder")
        self.setMinimumWidth(950)
        self.setMinimumHeight(650)
        self._build_ui()
        self._refresh_state()

    # ── UI construction ────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)

        splitter = QSplitter(Qt.Horizontal)

        # ── Left panel ────────────────────────────────────────────────────────
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 4, 0)

        # --- Input files group ---
        files_group = QGroupBox("Input Files")
        files_layout = QVBoxLayout()

        # customKcats file
        kcat_row = QHBoxLayout()
        kcat_row.addWidget(QLabel("customKcats file:"))
        self.kcat_edit = QLineEdit()
        self.kcat_edit.setPlaceholderText("TSV / CSV / Excel with kcat data")
        self.kcat_edit.setReadOnly(True)
        kcat_row.addWidget(self.kcat_edit)
        self.kcat_browse = QPushButton("Browse…")
        self.kcat_browse.setAutoDefault(False)
        self.kcat_browse.clicked.connect(self._browse_kcat)
        kcat_row.addWidget(self.kcat_browse)
        files_layout.addLayout(kcat_row)

        # uniprot file
        uni_row = QHBoxLayout()
        uni_row.addWidget(QLabel("UniProt data file:"))
        self.uni_edit = QLineEdit()
        self.uni_edit.setPlaceholderText("TSV / CSV / Excel with Entry, Mass, Sequence…")
        self.uni_edit.setReadOnly(True)
        uni_row.addWidget(self.uni_edit)
        self.uni_browse = QPushButton("Browse…")
        self.uni_browse.setAutoDefault(False)
        self.uni_browse.clicked.connect(self._browse_uniprot)
        uni_row.addWidget(self.uni_browse)
        files_layout.addLayout(uni_row)

        # status labels
        self.kcat_status = QLabel("kcat entries: —")
        self.uni_status = QLabel("UniProt entries: —")
        files_layout.addWidget(self.kcat_status)
        files_layout.addWidget(self.uni_status)

        files_group.setLayout(files_layout)
        left_layout.addWidget(files_group)

        # --- Protein pool parameters ---
        params_group = QGroupBox("Protein Pool Parameters")
        params_layout = QVBoxLayout()

        def _spin(label_text, default, lo, hi, tip):
            row = QHBoxLayout()
            lbl = QLabel(label_text)
            lbl.setToolTip(tip)
            spin = QDoubleSpinBox()
            spin.setRange(lo, hi)
            spin.setSingleStep(0.05)
            spin.setDecimals(4)
            spin.setValue(default)
            no_scroll(spin)
            row.addWidget(lbl)
            row.addWidget(spin)
            params_layout.addLayout(row)
            return spin

        self.sigma_spin = _spin(
            "σ (avg. enzyme saturation):", 0.5, 0.0001, 1.0,
            "Average in-vivo enzyme saturation (0–1). Default: 0.5"
        )
        self.f_spin = _spin(
            "f (metabolic enzyme fraction):", 0.5, 0.0001, 1.0,
            "Fraction of proteome mass that is metabolic enzymes (0–1). Default: 0.5"
        )
        self.ptot_spin = _spin(
            "Ptot (total protein, g/gDCW):", 0.5, 0.0001, 10.0,
            "Total cellular protein content (g protein / g dry cell weight). Default: 0.5"
        )
        self.pool_bound_label = QLabel("Pool UB: —")
        params_layout.addWidget(self.pool_bound_label)

        self.sigma_spin.valueChanged.connect(self._update_pool_label)
        self.f_spin.valueChanged.connect(self._update_pool_label)
        self.ptot_spin.valueChanged.connect(self._update_pool_label)
        self._update_pool_label()

        params_group.setLayout(params_layout)
        left_layout.addWidget(params_group)

        # --- Formulation mode ---
        mode_group = QGroupBox("Formulation Mode")
        mode_layout = QVBoxLayout()

        self.full_radio = QRadioButton("Full  (per-enzyme prot_{id} metabolites + usage reactions)")
        self.light_radio = QRadioButton("Light (direct prot_pool stoichiometry, fewer variables)")
        self.full_radio.setChecked(True)

        self._mode_group = QButtonGroup(self)
        self._mode_group.addButton(self.full_radio)
        self._mode_group.addButton(self.light_radio)

        mode_layout.addWidget(self.full_radio)
        mode_layout.addWidget(self.light_radio)
        mode_group.setLayout(mode_layout)
        left_layout.addWidget(mode_group)

        # --- UniKP button ---
        unikp_group = QGroupBox("Optional: kcat Prediction")
        unikp_layout = QVBoxLayout()
        unikp_info = QLabel(
            "Use UniKP (ProtT5-XL + Extra Trees) to predict kcat values\n"
            "for reactions that lack experimental data."
        )
        unikp_info.setWordWrap(True)
        unikp_layout.addWidget(unikp_info)
        self.unikp_btn = QPushButton("Predict kcat with UniKP…")
        self.unikp_btn.setAutoDefault(False)
        self.unikp_btn.clicked.connect(self._open_unikp_dialog)
        unikp_layout.addWidget(self.unikp_btn)
        unikp_group.setLayout(unikp_layout)
        left_layout.addWidget(unikp_group)

        left_layout.addStretch()

        # ── Right panel ───────────────────────────────────────────────────────
        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(4, 0, 0, 0)

        coverage_group = QGroupBox("Coverage Preview")
        cov_layout = QVBoxLayout()

        self.cov_total_label = QLabel("Enzyme-associated reactions: —")
        self.cov_covered_label = QLabel("Constrained (have kcat): —")
        self.cov_missing_label = QLabel("Unconstrained (no kcat): —")
        self.cov_mw_label = QLabel("Proteins with missing MW: —")
        for lbl in [self.cov_total_label, self.cov_covered_label,
                    self.cov_missing_label, self.cov_mw_label]:
            cov_layout.addWidget(lbl)

        cov_layout.addWidget(QLabel("Unconstrained reaction IDs:"))
        self.missing_rxns_table = QTableWidget(0, 1)
        self.missing_rxns_table.setHorizontalHeaderLabels(["Reaction ID"])
        self.missing_rxns_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.missing_rxns_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.missing_rxns_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.missing_rxns_table.setMaximumHeight(180)
        cov_layout.addWidget(self.missing_rxns_table)

        cov_layout.addWidget(QLabel("Proteins with missing MW or Sequence:"))
        self.missing_mw_table = QTableWidget(0, 1)
        self.missing_mw_table.setHorizontalHeaderLabels(["UniProt ID"])
        self.missing_mw_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.missing_mw_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.missing_mw_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.missing_mw_table.setMaximumHeight(120)
        cov_layout.addWidget(self.missing_mw_table)

        coverage_group.setLayout(cov_layout)
        right_layout.addWidget(coverage_group)
        right_layout.addStretch()

        # ── Assemble splitter ─────────────────────────────────────────────────
        splitter.addWidget(left)
        splitter.addWidget(right)
        splitter.setSizes([380, 570])
        root.addWidget(splitter)

        # ── Bottom button bar ─────────────────────────────────────────────────
        btn_bar = QHBoxLayout()

        self.preview_btn = QPushButton("Refresh Coverage Preview")
        self.preview_btn.setAutoDefault(False)
        self.preview_btn.clicked.connect(self._refresh_coverage)

        self.convert_btn = QPushButton("Convert to ecModel")
        self.convert_btn.setAutoDefault(False)
        self.convert_btn.clicked.connect(self._convert)

        self.revert_btn = QPushButton("Revert to GEM")
        self.revert_btn.setAutoDefault(False)
        self.revert_btn.clicked.connect(self._revert)

        self.close_btn = QPushButton("Close")
        self.close_btn.setAutoDefault(False)
        self.close_btn.clicked.connect(self.close)

        btn_bar.addWidget(self.preview_btn)
        btn_bar.addWidget(self.convert_btn)
        btn_bar.addWidget(self.revert_btn)
        btn_bar.addStretch()
        btn_bar.addWidget(self.close_btn)
        root.addLayout(btn_bar)

    # ── helpers ───────────────────────────────────────────────────────────────

    def _update_pool_label(self):
        ptot = self.ptot_spin.value()
        f = self.f_spin.value()
        sigma = self.sigma_spin.value()
        bound = ptot * f * sigma * 1000
        self.pool_bound_label.setText(f"Pool UB: {bound:.4f} mg/gDCW")

    def _refresh_state(self):
        """Sync UI from existing ECModelData (e.g. after project load)."""
        ec = self.appdata.project.ec_model_data
        self.sigma_spin.setValue(ec.sigma)
        self.f_spin.setValue(ec.f)
        self.ptot_spin.setValue(ec.ptot)
        self.light_radio.setChecked(ec.gecko_light)
        self.full_radio.setChecked(not ec.gecko_light)

        is_ec = ec.is_ecmodel
        self.convert_btn.setEnabled(not is_ec)
        self.revert_btn.setEnabled(is_ec)
        self.kcat_browse.setEnabled(not is_ec)
        self.uni_browse.setEnabled(not is_ec)
        self.full_radio.setEnabled(not is_ec)
        self.light_radio.setEnabled(not is_ec)
        self.sigma_spin.setEnabled(not is_ec)
        self.f_spin.setEnabled(not is_ec)
        self.ptot_spin.setEnabled(not is_ec)

        if ec.kcat_entries:
            self.kcat_status.setText(f"kcat entries: {len(ec.kcat_entries)}")
        if ec.uniprot_data:
            self.uni_status.setText(f"UniProt entries: {len(ec.uniprot_data)}")

    # ── slots ─────────────────────────────────────────────────────────────────

    @Slot()
    def _browse_kcat(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open customKcats file",
            self.appdata.work_directory,
            "Data files (*.tsv *.csv *.xlsx *.xls);;All files (*)",
            options=QFileDialog.DontUseNativeDialog,
        )
        if not path:
            return
        try:
            entries = parse_customkcats_file(path)
            self.appdata.project.ec_model_data.kcat_entries = entries
            self.kcat_edit.setText(path)
            self.kcat_status.setText(f"kcat entries: {len(entries)}")
            self._auto_refresh_coverage()
        except Exception as exc:
            QMessageBox.critical(self, "Error reading kcat file", str(exc))

    @Slot()
    def _browse_uniprot(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open UniProt data file",
            self.appdata.work_directory,
            "Data files (*.tsv *.csv *.xlsx *.xls);;All files (*)",
            options=QFileDialog.DontUseNativeDialog,
        )
        if not path:
            return
        try:
            data = parse_uniprot_file(path)
            self.appdata.project.ec_model_data.uniprot_data = data
            self.uni_edit.setText(path)
            self.uni_status.setText(f"UniProt entries: {len(data)}")
            self._auto_refresh_coverage()
        except Exception as exc:
            QMessageBox.critical(self, "Error reading UniProt file", str(exc))

    def _auto_refresh_coverage(self):
        """kcat 데이터가 있을 때만 자동으로 coverage를 갱신 (팝업 없음)."""
        if not self.appdata.project.ec_model_data.kcat_entries:
            return
        try:
            ec = self.appdata.project.ec_model_data
            model = self.appdata.project.cobra_py_model
            cov = compute_coverage(model, ec)
            self.cov_total_label.setText(f"Enzyme-associated reactions: {cov['total_enzyme_rxns']}")
            self.cov_covered_label.setText(f"Constrained (have kcat): {cov['covered_rxns']}")
            missing_n = cov['total_enzyme_rxns'] - cov['covered_rxns']
            self.cov_missing_label.setText(f"Unconstrained (no kcat): {missing_n}")
            self.cov_mw_label.setText(f"Proteins with missing MW: {len(cov['missing_mw_proteins'])}")
            missing_ids = cov.get("missing_rxn_ids", [])
            self.missing_rxns_table.setRowCount(len(missing_ids))
            for i, rxn_id in enumerate(sorted(missing_ids)):
                self.missing_rxns_table.setItem(i, 0, QTableWidgetItem(rxn_id))
            missing_mw = cov.get("missing_mw_proteins", [])
            self.missing_mw_table.setRowCount(len(missing_mw))
            for i, uid in enumerate(sorted(missing_mw)):
                self.missing_mw_table.setItem(i, 0, QTableWidgetItem(uid))
        except Exception:
            pass  # 자동 갱신 실패는 무시 (수동 버튼으로 에러 확인 가능)

    @Slot()
    def _refresh_coverage(self):
        ec = self.appdata.project.ec_model_data
        model = self.appdata.project.cobra_py_model
        if not ec.kcat_entries:
            QMessageBox.information(self, "No kcat data", "Please load a customKcats file first.")
            return
        try:
            cov = compute_coverage(model, ec)
        except Exception as exc:
            QMessageBox.critical(self, "Error computing coverage", str(exc))
            return

        self.cov_total_label.setText(f"Enzyme-associated reactions: {cov['total_enzyme_rxns']}")
        self.cov_covered_label.setText(f"Constrained (have kcat): {cov['covered_rxns']}")
        missing_n = cov['total_enzyme_rxns'] - cov['covered_rxns']
        self.cov_missing_label.setText(f"Unconstrained (no kcat): {missing_n}")
        self.cov_mw_label.setText(f"Proteins with missing MW: {len(cov['missing_mw_proteins'])}")

        # Populate missing reactions table
        missing_ids = cov.get("missing_rxn_ids", [])
        self.missing_rxns_table.setRowCount(len(missing_ids))
        for i, rxn_id in enumerate(sorted(missing_ids)):
            self.missing_rxns_table.setItem(i, 0, QTableWidgetItem(rxn_id))

        # Populate missing MW table
        missing_mw = cov.get("missing_mw_proteins", [])
        self.missing_mw_table.setRowCount(len(missing_mw))
        for i, uid in enumerate(sorted(missing_mw)):
            self.missing_mw_table.setItem(i, 0, QTableWidgetItem(uid))

    @Slot()
    def _convert(self):
        ec = self.appdata.project.ec_model_data
        if not ec.kcat_entries:
            QMessageBox.warning(self, "No kcat data", "Please load a customKcats file first.")
            return
        if not ec.uniprot_data:
            ret = QMessageBox.question(
                self, "No UniProt data",
                "No UniProt file loaded.  Molecular weights will be unavailable.\n"
                "Continue anyway?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No
            )
            if ret != QMessageBox.Yes:
                return

        # Save parameters from UI into ECModelData
        ec.sigma = self.sigma_spin.value()
        ec.f = self.f_spin.value()
        ec.ptot = self.ptot_spin.value()
        ec.gecko_light = self.light_radio.isChecked()

        model = self.appdata.project.cobra_py_model
        self.setCursor(Qt.BusyCursor)
        try:
            ecmodel, unconstrained, missing_mw = build_ecmodel(model, ec)
        except Exception as exc:
            self.setCursor(Qt.ArrowCursor)
            QMessageBox.critical(self, "ecModel build failed", str(exc))
            return
        self.setCursor(Qt.ArrowCursor)

        self.appdata.project.cobra_py_model = ecmodel
        self.appdata.project.store_original_bounds()
        self._refresh_state()
        self._refresh_coverage()
        self.parent().centralWidget().update(rebuild_all_tabs=True)

        msg_parts = [f"ecModel built successfully.\n"]
        if unconstrained:
            msg_parts.append(
                f"{len(unconstrained)} reaction(s) could not be constrained "
                f"(no kcat match).\n"
            )
        if missing_mw:
            msg_parts.append(
                f"{len(missing_mw)} protein(s) had no molecular weight data.\n"
            )
        msg_parts.append(
            "\nNote: Save the project (.cna) to persist the ecModel state."
        )
        QMessageBox.information(self, "ecModel built", "".join(msg_parts))

    @Slot()
    def _revert(self):
        ret = QMessageBox.question(
            self, "Revert to GEM",
            "This will remove all GECKO enzyme additions and restore the original GEM.\n"
            "Are you sure?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No
        )
        if ret != QMessageBox.Yes:
            return

        ec = self.appdata.project.ec_model_data
        model = self.appdata.project.cobra_py_model
        self.setCursor(Qt.BusyCursor)
        try:
            gem = revert_to_gem(model, ec)
        except Exception as exc:
            self.setCursor(Qt.ArrowCursor)
            QMessageBox.critical(self, "Revert failed", str(exc))
            return
        self.setCursor(Qt.ArrowCursor)

        self.appdata.project.cobra_py_model = gem
        self.appdata.project.store_original_bounds()
        self._refresh_state()
        # Clear tables
        self.missing_rxns_table.setRowCount(0)
        self.missing_mw_table.setRowCount(0)
        self.cov_total_label.setText("Enzyme-associated reactions: —")
        self.cov_covered_label.setText("Constrained (have kcat): —")
        self.cov_missing_label.setText("Unconstrained (no kcat): —")
        self.cov_mw_label.setText("Proteins with missing MW: —")
        self.parent().centralWidget().update(rebuild_all_tabs=True)
        QMessageBox.information(self, "Reverted", "Model has been reverted to the original GEM.")

    @Slot()
    def _open_unikp_dialog(self):
        from cnapy.gui_elements.unikp_dialog import UniKPDialog
        dlg = UniKPDialog(self.appdata, self)
        dlg.exec()
        # After UniKP dialog closes, refresh coverage if new entries were added
        ec = self.appdata.project.ec_model_data
        if ec.kcat_entries:
            self._refresh_coverage()
