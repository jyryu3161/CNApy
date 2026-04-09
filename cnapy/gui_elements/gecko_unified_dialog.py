"""
gecko_unified_dialog.py
=======================
Unified GECKO workflow dialog.

Left sidebar:
  [전제 조건]
    Build ecModel
  [분석 도구]
    Enzyme Usage Report
    Proteomics Integration
    kcat What-if Analysis

Right panel: QStackedWidget — content switches on sidebar click.
Analysis tools are disabled until ecModel is built.
"""

from __future__ import annotations

import cobra
from qtpy.QtCore import Qt, Slot
from qtpy.QtGui import QColor, QFont
from qtpy.QtWidgets import (
    QAbstractItemView,
    QButtonGroup,
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSizePolicy,
    QSplitter,
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from cnapy.appdata import AppData
from cnapy.gui_elements.filterable_combobox import FilterableComboBox
from cnapy.ecmodel.ecmodel_builder import (
    apply_kcat_multiplier,
    apply_proteomics,
    build_ecmodel,
    compute_coverage,
    flexibilize_enz_concs,
    get_enzyme_usage,
    parse_customkcats_file,
    parse_proteomics_file,
    parse_uniprot_file,
    remove_proteomics,
    revert_to_gem,
    run_fba_on_ecmodel,
)
from cnapy.utils import no_scroll


# ── helpers ───────────────────────────────────────────────────────────────────

class _NumItem(QTableWidgetItem):
    def __init__(self, value: float, fmt: str = "{:.4f}"):
        super().__init__(fmt.format(value))
        self._value = value
        self.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)

    def __lt__(self, other):
        if isinstance(other, _NumItem):
            return self._value < other._value
        return super().__lt__(other)


# ── sidebar nav button ────────────────────────────────────────────────────────

class _NavButton(QPushButton):
    """Sidebar navigation button with status indicator."""

    STYLE_NORMAL = """
        QPushButton {
            text-align: left; padding: 8px 12px;
            border: none; border-radius: 4px;
            background: transparent; color: #333;
        }
        QPushButton:hover { background: #e8edf5; }
        QPushButton:checked { background: #d0dff7; color: #1a4fa0; font-weight: bold; }
        QPushButton:disabled { color: #aaa; }
    """

    def __init__(self, text: str, parent=None):
        super().__init__(text, parent)
        self.setCheckable(True)
        self.setStyleSheet(self.STYLE_NORMAL)
        self.setAutoDefault(False)

    def set_status(self, status: str):
        """status: 'done' | 'active' | 'pending' | 'optional'"""
        icons = {"done": "✓  ", "active": "→  ", "pending": "○  ", "optional": "○  "}
        base = self.text().lstrip("✓→○ ")
        self.setText(icons.get(status, "○  ") + base)


# ── Page 0: Build ecModel ─────────────────────────────────────────────────────

class _BuildPage(QWidget):
    def __init__(self, appdata: AppData, dialog: "GeckoUnifiedDialog"):
        super().__init__()
        self.appdata = appdata
        self.dialog = dialog
        self._build_ui()
        self._refresh()

    def _build_ui(self):
        layout = QVBoxLayout(self)

        # ── kcat file ────────────────────────────────────────────────────────
        kcat_group = QGroupBox("kcat File  (customKcats)")
        kfl = QVBoxLayout()

        kcat_row = QHBoxLayout()
        kcat_row.addWidget(QLabel("File:"))
        self.kcat_edit = QLineEdit()
        self.kcat_edit.setPlaceholderText("TSV / CSV / Excel")
        self.kcat_edit.setReadOnly(True)
        kcat_row.addWidget(self.kcat_edit)
        self.kcat_btn = QPushButton("Browse…")
        self.kcat_btn.setAutoDefault(False)
        self.kcat_btn.clicked.connect(self._browse_kcat)
        kcat_row.addWidget(self.kcat_btn)
        kfl.addLayout(kcat_row)
        self.kcat_status = QLabel("kcat entries: —")
        kfl.addWidget(self.kcat_status)

        kcat_guide = QLabel(
            "<b>Required columns</b> (column names are case-insensitive):<br>"
            "&nbsp;&nbsp;• <tt>kcat</tt> — turnover number in 1/s, must be &gt; 0<br>"
            "&nbsp;&nbsp;• <tt>proteins</tt> — UniProt ID(s); use <tt>+</tt> for enzyme complexes "
            "<i>(e.g. P0A9P0+P12345)</i><br>"
            "&nbsp;&nbsp;&nbsp;&nbsp;<b>or</b> <tt>rxns</tt> — model reaction ID(s), comma-separated "
            "<i>(when proteins is absent)</i><br>"
            "<b>Optional columns</b>:<br>"
            "&nbsp;&nbsp;• <tt>stoicho</tt> — subunit copies per protein, <tt>+</tt>-separated "
            "<i>(default: 1)</i><br>"
            "&nbsp;&nbsp;• <tt>genes</tt> — gene IDs for GPR-based reaction matching<br>"
            "<br>"
            "<b>Example</b> (TSV/CSV):<br>"
            "<tt>proteins&nbsp;&nbsp;&nbsp;&nbsp;kcat&nbsp;&nbsp;&nbsp;rxns&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;stoicho<br>"
            "P0A9P0&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;45.2&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1<br>"
            "P0A9P0+P1B2&nbsp;12.0&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2+1<br>"
            "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;8.5&nbsp;&nbsp;&nbsp;PFK,PGI</tt>"
        )
        kcat_guide.setTextFormat(Qt.RichText)
        kcat_guide.setWordWrap(True)
        kcat_guide.setStyleSheet(
            "QLabel { background:#f8f9fa; padding:8px; "
            "border:1px solid #dde; border-radius:4px; font-size:11px; }"
        )
        kcat_guide_scroll = QScrollArea()
        kcat_guide_scroll.setWidgetResizable(True)
        kcat_guide_scroll.setWidget(kcat_guide)
        kcat_guide_scroll.setMaximumHeight(150)
        kcat_guide_scroll.setFrameShape(QScrollArea.NoFrame)
        kfl.addWidget(kcat_guide_scroll)
        kcat_group.setLayout(kfl)
        layout.addWidget(kcat_group)

        # ── UniProt file ──────────────────────────────────────────────────────
        uni_group = QGroupBox("UniProt Data File")
        ufl = QVBoxLayout()

        uni_row = QHBoxLayout()
        uni_row.addWidget(QLabel("File:"))
        self.uni_edit = QLineEdit()
        self.uni_edit.setPlaceholderText("TSV / CSV / Excel")
        self.uni_edit.setReadOnly(True)
        uni_row.addWidget(self.uni_edit)
        self.uni_btn = QPushButton("Browse…")
        self.uni_btn.setAutoDefault(False)
        self.uni_btn.clicked.connect(self._browse_uniprot)
        uni_row.addWidget(self.uni_btn)
        ufl.addLayout(uni_row)
        self.uni_status = QLabel("UniProt entries: —")
        ufl.addWidget(self.uni_status)

        uni_guide = QLabel(
            "<b>Required columns</b>:<br>"
            "&nbsp;&nbsp;• <tt>Entry</tt> — UniProt accession ID <i>(e.g. P0A9P0)</i><br>"
            "&nbsp;&nbsp;• <tt>Mass</tt> — molecular weight in Da <i>(e.g. 86376)</i>; "
            "also accepted: <tt>Molecular weight</tt>, <tt>MW</tt><br>"
            "<b>Optional columns</b>:<br>"
            "&nbsp;&nbsp;• <tt>Gene Names</tt> — gene name for display &amp; GPR matching "
            "<i>(any column containing \"gene\" and \"name\")</i><br>"
            "&nbsp;&nbsp;• <tt>Sequence</tt> — amino acid sequence <i>(for UniKP kcat prediction only)</i><br>"
            "<br>"
            "<b>Example</b> (TSV/CSV):<br>"
            "<tt>Entry&nbsp;&nbsp;&nbsp;Mass&nbsp;&nbsp;&nbsp;Gene Names<br>"
            "P0A9P0&nbsp;&nbsp;86376&nbsp;&nbsp;thrA<br>"
            "P12345&nbsp;&nbsp;34200&nbsp;&nbsp;aceE</tt>"
        )
        uni_guide.setTextFormat(Qt.RichText)
        uni_guide.setWordWrap(True)
        uni_guide.setStyleSheet(
            "QLabel { background:#f8f9fa; padding:8px; "
            "border:1px solid #dde; border-radius:4px; font-size:11px; }"
        )
        uni_guide_scroll = QScrollArea()
        uni_guide_scroll.setWidgetResizable(True)
        uni_guide_scroll.setWidget(uni_guide)
        uni_guide_scroll.setMaximumHeight(150)
        uni_guide_scroll.setFrameShape(QScrollArea.NoFrame)
        ufl.addWidget(uni_guide_scroll)
        uni_group.setLayout(ufl)
        layout.addWidget(uni_group)

        # Parameters
        params_group = QGroupBox("Protein Pool Parameters")
        pl = QVBoxLayout()

        def _spin(label, default, lo, hi, tip):
            row = QHBoxLayout()
            lbl = QLabel(label)
            lbl.setToolTip(tip)
            spin = QDoubleSpinBox()
            spin.setRange(lo, hi)
            spin.setSingleStep(0.05)
            spin.setDecimals(4)
            spin.setValue(default)
            no_scroll(spin)
            row.addWidget(lbl)
            row.addWidget(spin)
            pl.addLayout(row)
            return spin

        self.ptot_spin  = _spin("Ptot (g/gDCW):", 0.5, 0.0001, 10.0,
            "Total protein content per gram of cell dry weight.\n"
            "Typical range: 0.3–0.6 g/gDCW  |  Default: 0.5")
        self.f_spin     = _spin("f (metabolic enzyme fraction):", 0.5, 0.0001, 1.0,
            "Fraction of total protein allocated to metabolic enzymes.\n"
            "Excludes structural, ribosomal, and other non-metabolic proteins.\n"
            "Typical range: 0.3–0.7  |  Default: 0.5")
        self.sigma_spin = _spin("σ (avg. enzyme saturation):", 0.5, 0.0001, 1.0,
            "Average in vivo saturation of enzymes (fraction of Vmax achieved).\n"
            "σ = 1 means all enzymes work at full capacity; σ = 0.5 means half capacity.\n"
            "Typical range: 0.3–0.7  |  Default: 0.5")

        param_info = QLabel(
            "<small><b>Pool UB</b> = Ptot × f × σ × 1000  [mg/gDCW]<br>"
            "Hover over each parameter label for details.</small>"
        )
        param_info.setTextFormat(Qt.RichText)
        param_info.setStyleSheet("color: #666; padding-top: 2px;")
        pl.addWidget(param_info)

        self.pool_label = QLabel("Pool UB: —")
        pl.addWidget(self.pool_label)
        for sp in [self.sigma_spin, self.f_spin, self.ptot_spin]:
            sp.valueChanged.connect(self._update_pool_label)
        self._update_pool_label()
        params_group.setLayout(pl)
        layout.addWidget(params_group)

        # Mode
        mode_group = QGroupBox("Formulation Mode")
        ml = QVBoxLayout()
        self.full_radio  = QRadioButton("Full  (prot_{id} metabolites + usage reactions)")
        self.light_radio = QRadioButton("Light (direct prot_pool stoichiometry)")
        self.full_radio.setChecked(True)
        bg = QButtonGroup(self)
        bg.addButton(self.full_radio)
        bg.addButton(self.light_radio)
        ml.addWidget(self.full_radio)
        ml.addWidget(self.light_radio)
        mode_group.setLayout(ml)
        layout.addWidget(mode_group)

        # Coverage summary
        cov_group = QGroupBox("Coverage Preview")
        cl = QVBoxLayout()
        self.cov_total   = QLabel("Enzyme reactions: —")
        self.cov_covered = QLabel("Constrained: —")
        self.cov_missing = QLabel("Unconstrained: —")
        for lbl in [self.cov_total, self.cov_covered, self.cov_missing]:
            cl.addWidget(lbl)
        cov_group.setLayout(cl)
        layout.addWidget(cov_group)

        layout.addStretch()

        btn_row = QHBoxLayout()
        self.convert_btn = QPushButton("Convert to ecModel")
        self.convert_btn.setAutoDefault(False)
        self.convert_btn.clicked.connect(self._convert)
        self.revert_btn = QPushButton("Revert to GEM")
        self.revert_btn.setAutoDefault(False)
        self.revert_btn.clicked.connect(self._revert)
        btn_row.addWidget(self.convert_btn)
        btn_row.addWidget(self.revert_btn)
        btn_row.addStretch()
        layout.addLayout(btn_row)

    def _update_pool_label(self):
        b = self.ptot_spin.value() * self.f_spin.value() * self.sigma_spin.value() * 1000
        self.pool_label.setText(f"Pool UB: {b:.4f} mg/gDCW")

    def _refresh(self):
        ec = self.appdata.project.ec_model_data
        self.sigma_spin.setValue(ec.sigma)
        self.f_spin.setValue(ec.f)
        self.ptot_spin.setValue(ec.ptot)
        self.light_radio.setChecked(ec.gecko_light)
        self.full_radio.setChecked(not ec.gecko_light)
        is_ec = ec.is_ecmodel
        self.convert_btn.setEnabled(not is_ec)
        self.revert_btn.setEnabled(is_ec)
        for w in [self.kcat_btn, self.uni_btn, self.sigma_spin,
                  self.f_spin, self.ptot_spin, self.full_radio, self.light_radio]:
            w.setEnabled(not is_ec)
        if ec.kcat_entries:
            self.kcat_status.setText(f"kcat entries: {len(ec.kcat_entries)}")
        if ec.uniprot_data:
            self.uni_status.setText(f"UniProt entries: {len(ec.uniprot_data)}")

    def _auto_coverage(self):
        ec = self.appdata.project.ec_model_data
        if not ec.kcat_entries:
            return
        try:
            cov = compute_coverage(self.appdata.project.cobra_py_model, ec)
            self.cov_total.setText(f"Enzyme reactions: {cov['total_enzyme_rxns']}")
            self.cov_covered.setText(f"Constrained: {cov['covered_rxns']}")
            self.cov_missing.setText(
                f"Unconstrained: {cov['total_enzyme_rxns'] - cov['covered_rxns']}")
        except Exception:
            pass

    @Slot()
    def _browse_kcat(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open customKcats file", self.appdata.work_directory,
            "Data files (*.tsv *.csv *.xlsx *.xls);;All files (*)",
            options=QFileDialog.DontUseNativeDialog)
        if not path:
            return
        try:
            entries = parse_customkcats_file(path)
            self.appdata.project.ec_model_data.kcat_entries = entries
            self.kcat_edit.setText(path)
            self.kcat_status.setText(f"kcat entries: {len(entries)}")
            self._auto_coverage()
        except Exception as exc:
            QMessageBox.critical(self, "Error", str(exc))

    @Slot()
    def _browse_uniprot(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open UniProt file", self.appdata.work_directory,
            "Data files (*.tsv *.csv *.xlsx *.xls);;All files (*)",
            options=QFileDialog.DontUseNativeDialog)
        if not path:
            return
        try:
            data = parse_uniprot_file(path)
            self.appdata.project.ec_model_data.uniprot_data = data
            self.uni_edit.setText(path)
            self.uni_status.setText(f"UniProt entries: {len(data)}")
            self._auto_coverage()
        except Exception as exc:
            QMessageBox.critical(self, "Error", str(exc))

    @Slot()
    def _convert(self):
        ec = self.appdata.project.ec_model_data
        if not ec.kcat_entries:
            QMessageBox.warning(self, "No kcat data", "Please load a customKcats file first.")
            return
        ec.sigma = self.sigma_spin.value()
        ec.f = self.f_spin.value()
        ec.ptot = self.ptot_spin.value()
        ec.gecko_light = self.light_radio.isChecked()
        self.dialog.setCursor(Qt.BusyCursor)
        try:
            ecmodel, unconstrained, missing_mw = build_ecmodel(
                self.appdata.project.cobra_py_model, ec)
        except Exception as exc:
            self.dialog.setCursor(Qt.ArrowCursor)
            QMessageBox.critical(self, "Build failed", str(exc))
            return
        self.dialog.setCursor(Qt.ArrowCursor)
        self.appdata.project.cobra_py_model = ecmodel
        self.appdata.project.store_original_bounds()
        self._refresh()
        self._auto_coverage()
        self.dialog._refresh_nav()
        self.dialog._update_status_bar()
        self.dialog.parent().centralWidget().update(rebuild_all_tabs=True)
        msg = "ecModel built successfully."
        if unconstrained:
            msg += f"\nUnconstrained reactions (GPR but no kcat): {len(unconstrained)}"
        if missing_mw:
            msg += f"\nProteins with missing MW: {len(missing_mw)}"
        QMessageBox.information(self, "Done", msg)

    @Slot()
    def _revert(self):
        ret = QMessageBox.question(self, "Revert to GEM",
            "All GECKO additions will be removed and the original GEM will be restored.\nContinue?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if ret != QMessageBox.Yes:
            return
        ec = self.appdata.project.ec_model_data
        self.dialog.setCursor(Qt.BusyCursor)
        try:
            gem = revert_to_gem(self.appdata.project.cobra_py_model, ec)
        except Exception as exc:
            self.dialog.setCursor(Qt.ArrowCursor)
            QMessageBox.critical(self, "Revert failed", str(exc))
            return
        self.dialog.setCursor(Qt.ArrowCursor)
        self.appdata.project.cobra_py_model = gem
        self.appdata.project.store_original_bounds()
        self._refresh()
        self.dialog._refresh_nav()
        self.dialog._update_status_bar()
        self.dialog.parent().centralWidget().update(rebuild_all_tabs=True)
        QMessageBox.information(self, "Done", "Model reverted to the original GEM.")


# ── Page 1: Enzyme Usage ──────────────────────────────────────────────────────

class _UsagePage(QWidget):
    def __init__(self, appdata: AppData, dialog: "GeckoUnifiedDialog"):
        super().__init__()
        self.appdata = appdata
        self.dialog = dialog
        self._usage: dict = {}
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)

        # Summary
        sum_group = QGroupBox("Protein Pool Summary")
        sl = QHBoxLayout()
        self.pool_ub_lbl   = QLabel("Pool UB: —")
        self.pool_used_lbl = QLabel("Used: —")
        self.pool_pct_lbl  = QLabel("Pool usage: —")
        self.bottleneck_lbl = QLabel("Bottleneck (≥90%): —")
        for lbl in [self.pool_ub_lbl, self.pool_used_lbl,
                    self.pool_pct_lbl, self.bottleneck_lbl]:
            sl.addWidget(lbl)
        sl.addStretch()
        sum_group.setLayout(sl)
        layout.addWidget(sum_group)

        # Table
        self.table = QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels([
            "UniProt ID", "Gene", "Abs. Usage\n(mg/gDCW)",
            "UB (mg/gDCW)", "Cap. Usage (%)", "Catalysed Reactions"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(5, QHeaderView.Stretch)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSortingEnabled(True)
        layout.addWidget(self.table)

        legend = QHBoxLayout()
        for color, text in [("#ffcccc", "Bottleneck (≥90%)"),
                             ("#fff3cd", "High usage (70–90%)")]:
            lbl = QLabel("  ")
            lbl.setStyleSheet(f"background:{color};")
            legend.addWidget(lbl)
            legend.addWidget(QLabel(text))
            legend.addSpacing(8)
        legend.addStretch()
        layout.addLayout(legend)

        btn_row = QHBoxLayout()
        self.run_btn = QPushButton("Run FBA & Compute Usage")
        self.run_btn.setAutoDefault(False)
        self.run_btn.clicked.connect(self._run)
        self.export_btn = QPushButton("Export CSV…")
        self.export_btn.setAutoDefault(False)
        self.export_btn.setEnabled(False)
        self.export_btn.clicked.connect(self._export)
        btn_row.addWidget(self.run_btn)
        btn_row.addWidget(self.export_btn)
        btn_row.addStretch()
        layout.addLayout(btn_row)

    @Slot()
    def _run(self):
        ec = self.appdata.project.ec_model_data
        if not ec.is_ecmodel:
            QMessageBox.warning(self, "Not an ecModel", "Please build an ecModel first.")
            return
        ecmodel = self.appdata.project.cobra_py_model
        self.dialog.setCursor(Qt.BusyCursor)
        try:
            sol = run_fba_on_ecmodel(ecmodel)
        except Exception as exc:
            self.dialog.setCursor(Qt.ArrowCursor)
            QMessageBox.critical(self, "FBA failed", str(exc))
            return
        self.dialog.setCursor(Qt.ArrowCursor)
        if sol.status != "optimal":
            QMessageBox.warning(self, "No solution", f"FBA status: {sol.status}")
            return
        self._usage = get_enzyme_usage(ecmodel, sol)
        self._populate(ecmodel, sol, ec)
        self.export_btn.setEnabled(True)
        self.dialog._update_status_bar(sol)

    def _populate(self, ecmodel, sol, ec):
        try:
            pool_rxn = ecmodel.reactions.get_by_id(ec.protein_pool_rxn_id)
            pool_ub   = pool_rxn.upper_bound
            pool_used = abs(sol.fluxes.get(ec.protein_pool_rxn_id, 0.0))
            pool_pct  = (pool_used / pool_ub * 100) if pool_ub > 0 else 0.0
            self.pool_ub_lbl.setText(f"Pool UB: {pool_ub:.2f}")
            self.pool_used_lbl.setText(f"Used: {pool_used:.4f}")
            self.pool_pct_lbl.setText(f"Pool usage: {pool_pct:.1f}%")
        except KeyError:
            pool_ub = 1.0

        rows = []
        for uid, data in self._usage.items():
            gene = ec.uniprot_data.get(uid, {}).get("gene", "")
            met_id = f"prot_{uid}"
            try:
                met = ecmodel.metabolites.get_by_id(met_id)
                rxn_ids = [r.id for r in met.reactions
                           if not r.id.startswith("usage_prot_")
                           and abs(sol.fluxes.get(r.id, 0.0)) > 1e-9]
            except KeyError:
                rxn_ids = []
            rows.append((uid, gene, data["abs_usage"], data["ub"],
                         data["cap_usage"] * 100, ", ".join(rxn_ids)))

        rows.sort(key=lambda x: x[4], reverse=True)
        self.bottleneck_lbl.setText(
            f"Bottleneck (≥90%): {sum(1 for r in rows if r[4] >= 90)}")

        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(rows))
        for i, (uid, gene, abs_u, ub, cap_pct, rxns) in enumerate(rows):
            self.table.setItem(i, 0, QTableWidgetItem(uid))
            self.table.setItem(i, 1, QTableWidgetItem(gene))
            self.table.setItem(i, 2, _NumItem(abs_u))
            self.table.setItem(i, 3, _NumItem(ub))
            self.table.setItem(i, 4, _NumItem(cap_pct, "{:.1f}%"))
            self.table.setItem(i, 5, QTableWidgetItem(rxns))
            color = ("#ffcccc" if cap_pct >= 90
                     else "#fff3cd" if cap_pct >= 70 else "#ffffff")
            for col in range(6):
                item = self.table.item(i, col)
                if item:
                    item.setBackground(QColor(color))
        self.table.setSortingEnabled(True)

    @Slot()
    def _export(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Export CSV", self.appdata.work_directory,
            "CSV (*.csv)", options=QFileDialog.DontUseNativeDialog)
        if not path:
            return
        if not path.endswith(".csv"):
            path += ".csv"
        try:
            with open(path, "w", encoding="utf-8") as f:
                f.write("UniProt ID,Gene,Abs Usage,UB,Cap Usage (%),Reactions\n")
                for row in range(self.table.rowCount()):
                    vals = [self.table.item(row, c).text() for c in range(6)]
                    f.write(",".join(f'"{v}"' for v in vals) + "\n")
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc))


# ── Page 2: Proteomics ────────────────────────────────────────────────────────

class _ProteomicsPage(QWidget):
    def __init__(self, appdata: AppData, dialog: "GeckoUnifiedDialog"):
        super().__init__()
        self.appdata = appdata
        self.dialog = dialog
        self._prot_data: dict = {}
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)

        info = QLabel(
            "Load experimentally measured protein abundances (mg/gDCW) to set individual "
            "upper bounds on each enzyme's usage reaction. "
            "Enzymes without measured data continue to draw from the shared protein pool."
        )
        info.setWordWrap(True)
        info.setStyleSheet(
            "QLabel { background:#f0f4ff; padding:8px; "
            "border:1px solid #ccd; border-radius:4px; }")
        layout.addWidget(info)

        file_group = QGroupBox("Proteomics File")
        fl = QVBoxLayout()
        file_row = QHBoxLayout()
        self.file_edit = QLineEdit()
        self.file_edit.setPlaceholderText("TSV / CSV / Excel")
        self.file_edit.setReadOnly(True)
        file_row.addWidget(self.file_edit)
        browse_btn = QPushButton("Browse…")
        browse_btn.setAutoDefault(False)
        browse_btn.clicked.connect(self._browse)
        file_row.addWidget(browse_btn)
        fl.addLayout(file_row)
        self.file_status = QLabel("Entries: —")
        fl.addWidget(self.file_status)

        prot_guide = QLabel(
            "<b>Required columns</b> (two columns, header is optional):<br>"
            "&nbsp;&nbsp;• Column 1 — UniProt accession ID <i>(e.g. P0A9P0)</i><br>"
            "&nbsp;&nbsp;• Column 2 — protein abundance in mg/gDCW, must be &gt; 0<br>"
            "<br>"
            "<b>Example</b> (TSV/CSV):<br>"
            "<tt>UniProt&nbsp;&nbsp;level<br>"
            "P0A9P0&nbsp;&nbsp;0.00312<br>"
            "P12345&nbsp;&nbsp;0.00087</tt>"
        )
        prot_guide.setTextFormat(Qt.RichText)
        prot_guide.setWordWrap(True)
        prot_guide.setStyleSheet(
            "QLabel { background:#f8f9fa; padding:8px; "
            "border:1px solid #dde; border-radius:4px; font-size:11px; }"
        )
        prot_guide_scroll = QScrollArea()
        prot_guide_scroll.setWidgetResizable(True)
        prot_guide_scroll.setWidget(prot_guide)
        prot_guide_scroll.setMaximumHeight(150)
        prot_guide_scroll.setFrameShape(QScrollArea.NoFrame)
        fl.addWidget(prot_guide_scroll)
        file_group.setLayout(fl)
        layout.addWidget(file_group)

        preview_group = QGroupBox("Preview (ecModel Matching)")
        pl = QVBoxLayout()
        self.match_label = QLabel("Matched: — / —")
        pl.addWidget(self.match_label)
        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(["UniProt ID", "Level (mg/gDCW)", "In ecModel?"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setMaximumHeight(200)
        pl.addWidget(self.table)
        preview_group.setLayout(pl)
        layout.addWidget(preview_group)

        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)
        layout.addStretch()

        btn_row = QHBoxLayout()
        self.apply_btn = QPushButton("Apply Proteomics Constraints")
        self.apply_btn.setAutoDefault(False)
        self.apply_btn.setEnabled(False)
        self.apply_btn.clicked.connect(self._apply)
        self.flex_btn = QPushButton("Flexibilize & Re-run FBA")
        self.flex_btn.setAutoDefault(False)
        self.flex_btn.setEnabled(False)
        self.flex_btn.setToolTip(
            "Use when FBA is infeasible: relaxes individual enzyme constraints one by one\n"
            "until FBA becomes feasible, keeping only the minimum necessary relaxations.")
        self.flex_btn.clicked.connect(self._flexibilize)
        self.remove_btn = QPushButton("Remove Constraints")
        self.remove_btn.setAutoDefault(False)
        self.remove_btn.setEnabled(False)
        self.remove_btn.clicked.connect(self._remove)
        btn_row.addWidget(self.apply_btn)
        btn_row.addWidget(self.flex_btn)
        btn_row.addWidget(self.remove_btn)
        btn_row.addStretch()
        layout.addLayout(btn_row)

    def _populate_table(self):
        ec = self.appdata.project.ec_model_data
        ecmodel = self.appdata.project.cobra_py_model
        rxn_ids = {r.id for r in ecmodel.reactions}
        self.table.setRowCount(len(self._prot_data))
        matched = 0
        for i, (uid, level) in enumerate(sorted(self._prot_data.items())):
            in_model = f"usage_prot_{uid}" in rxn_ids
            if in_model:
                matched += 1
            self.table.setItem(i, 0, QTableWidgetItem(uid))
            self.table.setItem(i, 1, QTableWidgetItem(f"{level:.6f}"))
            item = QTableWidgetItem("Yes" if in_model else "No")
            item.setForeground(Qt.darkGreen if in_model else Qt.red)
            self.table.setItem(i, 2, item)
        self.match_label.setText(
            f"Matched: {matched} / {len(self._prot_data)}")
        self.apply_btn.setEnabled(matched > 0 and ec.is_ecmodel)

    @Slot()
    def _browse(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open Proteomics File", self.appdata.work_directory,
            "Data files (*.tsv *.csv *.xlsx *.xls);;All files (*)",
            options=QFileDialog.DontUseNativeDialog)
        if not path:
            return
        try:
            self._prot_data = parse_proteomics_file(path)
            self.file_edit.setText(path)
            self.file_status.setText(f"Entries: {len(self._prot_data)}")
            self._populate_table()
        except Exception as exc:
            QMessageBox.critical(self, "Error", str(exc))

    @Slot()
    def _apply(self):
        ec = self.appdata.project.ec_model_data
        ecmodel = self.appdata.project.cobra_py_model
        self.dialog.setCursor(Qt.BusyCursor)
        try:
            n, missing = apply_proteomics(ecmodel, ec, self._prot_data)
        except Exception as exc:
            self.dialog.setCursor(Qt.ArrowCursor)
            QMessageBox.critical(self, "Error", str(exc))
            return
        self.dialog.setCursor(Qt.ArrowCursor)
        self.remove_btn.setEnabled(True)
        self.flex_btn.setEnabled(True)
        self.apply_btn.setEnabled(False)
        msg = f"Proteomics constraints applied to {n} enzyme(s)."
        if missing:
            msg += f"\nNot found in ecModel: {', '.join(missing[:5])}"
            if len(missing) > 5:
                msg += f" …(+{len(missing)-5})"
        msg += "\nIf FBA becomes infeasible, click 'Flexibilize & Re-run FBA' to relax constraints."
        self.status_label.setText(msg)
        self.status_label.setStyleSheet("color: green;")
        # Refresh main window so UB changes are visible immediately
        try:
            self.dialog.parent().centralWidget().update(rebuild_all_tabs=True)
        except Exception:
            pass

    @Slot()
    def _flexibilize(self):
        ec = self.appdata.project.ec_model_data
        ecmodel = self.appdata.project.cobra_py_model
        self.dialog.setCursor(Qt.BusyCursor)
        try:
            relaxed = flexibilize_enz_concs(ecmodel, self._prot_data)
        except Exception as exc:
            self.dialog.setCursor(Qt.ArrowCursor)
            QMessageBox.critical(self, "Error", str(exc))
            return
        self.dialog.setCursor(Qt.ArrowCursor)

        if relaxed:
            msg = (f"Relaxed constraints on {len(relaxed)} enzyme(s) — FBA is now feasible.\n"
                   f"Relaxed: {', '.join(relaxed[:10])}")
            if len(relaxed) > 10:
                msg += f" …(+{len(relaxed)-10})"
            self.status_label.setText(msg)
            self.status_label.setStyleSheet("color: darkorange;")
        else:
            # Check if FBA is still infeasible
            sol = run_fba_on_ecmodel(ecmodel)
            if sol.status != "optimal":
                self.status_label.setText(
                    "FBA remains infeasible after flexibilization. "
                    "Consider increasing the protein pool UB or checking model constraints.")
                self.status_label.setStyleSheet("color: red;")
                return
            self.status_label.setText("FBA is feasible with all constraints retained.")
            self.status_label.setStyleSheet("color: green;")

        self.flex_btn.setEnabled(False)
        try:
            self.dialog.parent().centralWidget().update(rebuild_all_tabs=True)
        except Exception:
            pass

    @Slot()
    def _remove(self):
        ret = QMessageBox.question(self, "Remove Constraints",
            "All enzymes will be reconnected to prot_pool and individual bounds restored.\nContinue?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if ret != QMessageBox.Yes:
            return
        ec = self.appdata.project.ec_model_data
        n = remove_proteomics(self.appdata.project.cobra_py_model, ec)
        self.remove_btn.setEnabled(False)
        self.flex_btn.setEnabled(False)
        self.apply_btn.setEnabled(bool(self._prot_data))
        self.status_label.setText(f"Removed constraints from {n} enzyme(s). All reconnected to prot_pool.")
        self.status_label.setStyleSheet("color: blue;")
        try:
            self.dialog.parent().centralWidget().update(rebuild_all_tabs=True)
        except Exception:
            pass


# ── Page 3: kcat What-if ──────────────────────────────────────────────────────

class _WhatifPage(QWidget):
    def __init__(self, appdata: AppData, dialog: "GeckoUnifiedDialog"):
        super().__init__()
        self.appdata = appdata
        self.dialog = dialog
        self._baseline: dict | None = None
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)

        info = QLabel(
            "Simulates the effect of scaling an enzyme's kcat on the FBA objective (e.g. growth rate). "
            "The original model is never modified — analysis runs on an internal copy."
        )
        info.setWordWrap(True)
        info.setStyleSheet(
            "QLabel { background:#f0fff4; padding:8px; "
            "border:1px solid #cdc; border-radius:4px; }")
        layout.addWidget(info)

        # Baseline
        bl_group = QGroupBox("Baseline (current ecModel FBA)")
        bl_layout = QHBoxLayout()
        self.bl_obj_lbl  = QLabel("Objective: —")
        self.bl_pool_lbl = QLabel("Pool used: —")
        bl_layout.addWidget(self.bl_obj_lbl)
        bl_layout.addWidget(self.bl_pool_lbl)
        bl_layout.addStretch()
        self.bl_btn = QPushButton("Run Baseline FBA")
        self.bl_btn.setAutoDefault(False)
        self.bl_btn.clicked.connect(self._run_baseline)
        bl_layout.addWidget(self.bl_btn)
        bl_group.setLayout(bl_layout)
        layout.addWidget(bl_group)

        # Settings
        settings_group = QGroupBox("What-if Settings")
        form = QFormLayout()

        enzyme_row = QHBoxLayout()
        self.enzyme_combo = FilterableComboBox()
        self.enzyme_combo.setMinimumWidth(280)
        self.enzyme_combo.setMaxVisibleItems(15)
        self.enzyme_combo.lineEdit().setPlaceholderText("Type UniProt ID or gene name to search…")
        enzyme_row.addWidget(self.enzyme_combo)
        self.refresh_enzyme_btn = QPushButton("↺")
        self.refresh_enzyme_btn.setFixedWidth(28)
        self.refresh_enzyme_btn.setToolTip("Reload enzyme list from ecModel")
        self.refresh_enzyme_btn.setAutoDefault(False)
        self.refresh_enzyme_btn.clicked.connect(self._populate_enzyme_combo)
        enzyme_row.addWidget(self.refresh_enzyme_btn)
        enzyme_widget = QWidget()
        enzyme_widget.setLayout(enzyme_row)
        form.addRow("Enzyme:", enzyme_widget)

        self.mult_spin = QDoubleSpinBox()
        self.mult_spin.setRange(0.01, 1000.0)
        self.mult_spin.setSingleStep(0.5)
        self.mult_spin.setDecimals(2)
        self.mult_spin.setValue(2.0)
        self.mult_spin.setToolTip("2× = kcat doubled (half the protein cost per unit flux)")
        no_scroll(self.mult_spin)
        form.addRow("kcat multiplier (×):", self.mult_spin)
        settings_group.setLayout(form)
        layout.addWidget(settings_group)

        self.run_btn = QPushButton("Run What-if FBA")
        self.run_btn.setAutoDefault(False)
        self.run_btn.clicked.connect(self._run_whatif)
        layout.addWidget(self.run_btn)

        # Result
        result_group = QGroupBox("Result Comparison")
        rl = QVBoxLayout()
        rl.setSpacing(4)
        rl.setContentsMargins(8, 6, 8, 6)
        comp_row = QHBoxLayout()
        comp_row.setSpacing(16)
        self.res_obj_lbl  = QLabel("Modified objective: —")
        self.res_diff_lbl = QLabel("Change: —")
        self.res_pool_lbl = QLabel("Pool used: —")
        comp_row.addWidget(self.res_obj_lbl)
        comp_row.addWidget(self.res_diff_lbl)
        comp_row.addWidget(self.res_pool_lbl)
        comp_row.addStretch()
        rl.addLayout(comp_row)

        rl.addWidget(QLabel("Top 10 enzyme usage changes:"))
        self.diff_table = QTableWidget(0, 4)
        self.diff_table.setHorizontalHeaderLabels([
            "UniProt ID", "Baseline cap. (%)", "Modified cap. (%)", "Δ (%)"])
        self.diff_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.diff_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.diff_table.setMaximumHeight(200)
        rl.addWidget(self.diff_table)
        result_group.setLayout(rl)

        result_group.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)

        layout.addWidget(result_group)
        layout.addStretch()

    def _populate_enzyme_combo(self):
        """Populate enzyme combo box from usage_prot_* reactions in the ecModel."""
        ec = self.appdata.project.ec_model_data
        self.enzyme_combo.clear()
        if not ec.is_ecmodel:
            return
        ecmodel = self.appdata.project.cobra_py_model
        enzymes = []
        for rxn in ecmodel.reactions:
            if rxn.id.startswith("usage_prot_"):
                uid = rxn.id[len("usage_prot_"):]
                gene = ec.uniprot_data.get(uid, {}).get("gene", "")
                label = f"{uid}  ({gene})" if gene else uid
                enzymes.append((label, uid))
        enzymes.sort(key=lambda x: x[0])
        for label, uid in enzymes:
            self.enzyme_combo.addItem(label, userData=uid)

    @Slot()
    def _run_baseline(self):
        ec = self.appdata.project.ec_model_data
        if not ec.is_ecmodel:
            QMessageBox.warning(self, "Not an ecModel", "Please build an ecModel first.")
            return
        ecmodel = self.appdata.project.cobra_py_model
        self.dialog.setCursor(Qt.BusyCursor)
        try:
            sol = run_fba_on_ecmodel(ecmodel)
        except Exception as exc:
            self.dialog.setCursor(Qt.ArrowCursor)
            self.bl_obj_lbl.setText(f"FBA failed: {exc}")
            return
        self.dialog.setCursor(Qt.ArrowCursor)
        if sol.status != "optimal":
            self.bl_obj_lbl.setText(f"Status: {sol.status}")
            return
        self._baseline = {
            "objective": sol.objective_value,
            "usage": get_enzyme_usage(ecmodel, sol),
        }
        pool_used = abs(sol.fluxes.get(ec.protein_pool_rxn_id, 0.0))
        self.bl_obj_lbl.setText(f"Objective: {sol.objective_value:.6f}")
        self.bl_pool_lbl.setText(f"Pool used: {pool_used:.4f}")
        self.dialog._update_status_bar(sol)

    @Slot()
    def _run_whatif(self):
        ec = self.appdata.project.ec_model_data
        if not ec.is_ecmodel:
            QMessageBox.warning(self, "Not an ecModel", "Please build an ecModel first.")
            return
        if self._baseline is None:
            QMessageBox.warning(self, "No baseline", "Please run the baseline FBA first.")
            return
        uid = self.enzyme_combo.currentData()
        if not uid:
            # Fallback: match typed text against combo labels or uid directly
            typed = self.enzyme_combo.currentText().strip()
            for i in range(self.enzyme_combo.count()):
                item_uid = self.enzyme_combo.itemData(i)
                item_text = self.enzyme_combo.itemText(i)
                if typed == item_uid or typed in item_text:
                    uid = item_uid
                    break
        if not uid:
            QMessageBox.warning(self, "No enzyme selected",
                                "Please select an enzyme from the list or type a UniProt ID / gene name.")
            return
        ecmodel = self.appdata.project.cobra_py_model
        mult = self.mult_spin.value()
        self.dialog.setCursor(Qt.BusyCursor)
        try:
            modified = ecmodel.copy()
            n = apply_kcat_multiplier(modified, ec, uid, mult)
            if n == 0:
                self.dialog.setCursor(Qt.ArrowCursor)
                QMessageBox.warning(self, "No reactions updated",
                    f"Enzyme {uid} has no constrained reactions in the ecModel.")
                return
            sol = run_fba_on_ecmodel(modified)
        except Exception as exc:
            self.dialog.setCursor(Qt.ArrowCursor)
            QMessageBox.critical(self, "FBA failed", str(exc))
            return
        self.dialog.setCursor(Qt.ArrowCursor)
        if sol.status != "optimal":
            self.res_obj_lbl.setText(f"Status: {sol.status}")
            return

        new_obj = sol.objective_value
        old_obj = self._baseline["objective"]
        diff = new_obj - old_obj
        pct  = (diff / old_obj * 100) if old_obj != 0 else 0.0
        sign = "+" if diff >= 0 else ""
        color = "green" if diff >= 0 else "red"

        self.res_obj_lbl.setText(f"Modified: {new_obj:.6f}")
        self.res_diff_lbl.setText(
            f"<span style='color:{color}'>{sign}{diff:.6f} ({sign}{pct:.2f}%)</span>")
        self.res_diff_lbl.setTextFormat(Qt.RichText)
        pool_used = abs(sol.fluxes.get(ec.protein_pool_rxn_id, 0.0))
        self.res_pool_lbl.setText(f"Pool used: {pool_used:.4f}")

        new_usage = get_enzyme_usage(modified, sol)
        base_usage = self._baseline["usage"]
        diffs = sorted(
            [(u,
              base_usage.get(u, {}).get("cap_usage", 0.0) * 100,
              new_usage.get(u, {}).get("cap_usage", 0.0) * 100,
              (new_usage.get(u, {}).get("cap_usage", 0.0) -
               base_usage.get(u, {}).get("cap_usage", 0.0)) * 100)
             for u in set(base_usage) | set(new_usage)],
            key=lambda x: abs(x[3]), reverse=True)[:10]

        self.diff_table.setRowCount(len(diffs))
        for i, (u, bc, nc, delta) in enumerate(diffs):
            self.diff_table.setItem(i, 0, QTableWidgetItem(u))
            self.diff_table.setItem(i, 1, QTableWidgetItem(f"{bc:.1f}%"))
            self.diff_table.setItem(i, 2, QTableWidgetItem(f"{nc:.1f}%"))
            s = "+" if delta >= 0 else ""
            d_item = QTableWidgetItem(f"{s}{delta:.1f}%")
            d_item.setForeground(Qt.darkGreen if delta < 0 else Qt.darkRed)
            self.diff_table.setItem(i, 3, d_item)


# ── Unified dialog ────────────────────────────────────────────────────────────

class GeckoUnifiedDialog(QDialog):
    """Unified GECKO workflow dialog with sidebar navigation."""

    _NAV_ITEMS = [
        # (page_index, key, label_base)
        (0, "build",      "Build ecModel"),
        (1, "usage",      "Enzyme Usage Report"),
        (2, "proteomics", "Proteomics Integration"),
        (3, "whatif",     "kcat What-if Analysis"),
    ]

    def __init__(self, appdata: AppData, parent=None):
        super().__init__(parent)
        self.appdata = appdata
        self.setWindowTitle("GECKO Enzyme-Constrained Modeling")
        self.setMinimumWidth(1000)
        self.setMinimumHeight(720)
        self.resize(1100, 820)
        self._nav_btns: dict[str, _NavButton] = {}
        self._build_ui()
        self._refresh_nav()
        self._update_status_bar()

    # ── UI ────────────────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        splitter = QSplitter(Qt.Horizontal)

        # ── Left sidebar ──────────────────────────────────────────────────────
        sidebar = QWidget()
        sidebar.setMinimumWidth(190)
        sidebar.setMaximumWidth(220)
        sidebar.setStyleSheet("QWidget { background: #f4f6fa; }")
        sb_layout = QVBoxLayout(sidebar)
        sb_layout.setContentsMargins(8, 12, 8, 8)
        sb_layout.setSpacing(2)

        # Section: Prerequisites
        prereq_lbl = QLabel("Prerequisites")
        prereq_lbl.setStyleSheet(
            "color:#888; font-size:11px; font-weight:bold; "
            "padding:4px 4px 2px 4px;")
        sb_layout.addWidget(prereq_lbl)

        btn_build = _NavButton("○  Build ecModel")
        self._nav_btns["build"] = btn_build
        btn_build.clicked.connect(lambda: self._switch(0, "build"))
        sb_layout.addWidget(btn_build)

        sb_layout.addSpacing(10)

        # Section: Analysis Tools
        tools_lbl = QLabel("Analysis Tools")
        tools_lbl.setStyleSheet(
            "color:#888; font-size:11px; font-weight:bold; "
            "padding:4px 4px 2px 4px;")
        sb_layout.addWidget(tools_lbl)

        for key, label in [("usage",      "○  Enzyme Usage Report"),
                            ("proteomics", "○  Proteomics Integration"),
                            ("whatif",     "○  kcat What-if Analysis")]:
            btn = _NavButton(label)
            self._nav_btns[key] = btn
            idx = {"usage": 1, "proteomics": 2, "whatif": 3}[key]
            btn.clicked.connect(lambda checked, i=idx, k=key: self._switch(i, k))
            sb_layout.addWidget(btn)

        sb_layout.addStretch()

        # Status bar
        sep = QWidget()
        sep.setFixedHeight(1)
        sep.setStyleSheet("background:#ddd;")
        sb_layout.addWidget(sep)

        self.status_model_lbl = QLabel("Model: GEM")
        self.status_pool_lbl  = QLabel("")
        for lbl in [self.status_model_lbl, self.status_pool_lbl]:
            lbl.setStyleSheet("color:#555; font-size:11px; padding:2px 4px;")
            lbl.setWordWrap(True)
            sb_layout.addWidget(lbl)

        # Close button
        close_btn = QPushButton("Close")
        close_btn.setAutoDefault(False)
        close_btn.clicked.connect(self.close)
        sb_layout.addWidget(close_btn)

        # ── Right: stacked pages ──────────────────────────────────────────────
        self.stack = QStackedWidget()

        self._page_build = _BuildPage(self.appdata, self)
        self._page_usage = _UsagePage(self.appdata, self)
        self._page_prot  = _ProteomicsPage(self.appdata, self)
        self._page_wif   = _WhatifPage(self.appdata, self)

        for page in [self._page_build, self._page_usage,
                     self._page_prot, self._page_wif]:
            container = QWidget()
            cl = QVBoxLayout(container)
            cl.setContentsMargins(16, 16, 16, 16)
            cl.addWidget(page)
            self.stack.addWidget(container)

        splitter.addWidget(sidebar)
        splitter.addWidget(self.stack)
        splitter.setSizes([200, 800])
        splitter.setHandleWidth(1)

        root.addWidget(splitter)

        # Start on build page
        self._switch(0, "build")

    # ── helpers ───────────────────────────────────────────────────────────────

    def _switch(self, index: int, key: str):
        self.stack.setCurrentIndex(index)
        for k, btn in self._nav_btns.items():
            btn.setChecked(k == key)

    def _refresh_nav(self):
        ec = self.appdata.project.ec_model_data
        is_ec = ec.is_ecmodel

        # Build button
        self._nav_btns["build"].set_status("done" if is_ec else "active")

        # Analysis tool buttons
        for key in ["usage", "proteomics", "whatif"]:
            btn = self._nav_btns[key]
            btn.setEnabled(is_ec)
            btn.set_status("pending" if not is_ec else "optional")

        # Populate enzyme combo whenever ecModel state changes
        self._page_wif._populate_enzyme_combo()

    def _update_status_bar(self, solution=None):
        ec = self.appdata.project.ec_model_data
        if ec.is_ecmodel:
            self.status_model_lbl.setText("Model: ecModel ✓")
            self.status_model_lbl.setStyleSheet(
                "color:#1a7a3a; font-size:11px; padding:2px 4px;")
        else:
            self.status_model_lbl.setText("Model: GEM")
            self.status_model_lbl.setStyleSheet(
                "color:#555; font-size:11px; padding:2px 4px;")

        if solution is not None and ec.is_ecmodel:
            try:
                model = self.appdata.project.cobra_py_model
                pool_rxn = model.reactions.get_by_id(ec.protein_pool_rxn_id)
                pool_ub   = pool_rxn.upper_bound
                pool_used = abs(solution.fluxes.get(ec.protein_pool_rxn_id, 0.0))
                self.status_pool_lbl.setText(
                    f"Pool: {pool_used:.2f} / {pool_ub:.2f} mg/gDCW")
            except Exception:
                pass
        elif not ec.is_ecmodel:
            self.status_pool_lbl.setText("")
