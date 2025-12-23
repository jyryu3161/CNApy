"""External Flux Data Loading Dialog for CNApy

This module provides functionality to:
- Load external flux data (CSV/TSV) for visualization
- Support multiple conditions/files
- Calculate and visualize log2 fold change between conditions
"""

import os
import csv
import numpy as np
from qtpy.QtCore import Qt, Slot, Signal
from qtpy.QtWidgets import (QDialog, QHBoxLayout, QLabel, QVBoxLayout,
                            QPushButton, QGroupBox, QTextEdit, QTabWidget,
                            QWidget, QMessageBox, QCheckBox, QListWidget,
                            QListWidgetItem, QAbstractItemView, QFileDialog,
                            QComboBox, QDoubleSpinBox, QTableWidget,
                            QTableWidgetItem, QHeaderView, QRadioButton,
                            QButtonGroup, QSpinBox, QLineEdit, QGridLayout)
from qtpy.QtGui import QColor
from typing import Dict, List, Tuple, Optional
import math
from cnapy.appdata import AppData


class FluxDataCondition:
    """Represents a single flux data condition loaded from a file."""

    def __init__(self, name: str, file_path: str = ""):
        self.name = name
        self.file_path = file_path
        self.flux_data: Dict[str, float] = {}  # reaction_id -> flux value

    def load_from_file(self, file_path: str, reaction_col: int = 0,
                       flux_col: int = 1, has_header: bool = True,
                       delimiter: str = ',') -> Tuple[bool, str]:
        """Load flux data from CSV/TSV file.

        Returns:
            Tuple of (success, message)
        """
        self.file_path = file_path
        self.flux_data = {}

        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                reader = csv.reader(f, delimiter=delimiter)

                if has_header:
                    next(reader)  # Skip header

                for row in reader:
                    if len(row) > max(reaction_col, flux_col):
                        rxn_id = row[reaction_col].strip()
                        try:
                            flux_val = float(row[flux_col].strip())
                            self.flux_data[rxn_id] = flux_val
                        except ValueError:
                            pass  # Skip rows with non-numeric flux

            return True, f"Loaded {len(self.flux_data)} reactions from {os.path.basename(file_path)}"
        except Exception as e:
            return False, f"Error loading file: {str(e)}"


def calculate_log2_fold_change(condition1: FluxDataCondition,
                               condition2: FluxDataCondition,
                               pseudocount: float = 1e-6) -> Dict[str, float]:
    """
    Calculate log2 fold change between two conditions.

    log2FC = log2(condition2 / condition1)

    Args:
        condition1: Reference condition (denominator)
        condition2: Target condition (numerator)
        pseudocount: Small value added to avoid division by zero

    Returns:
        Dict mapping reaction_id to log2 fold change value
    """
    log2fc = {}

    all_reactions = set(condition1.flux_data.keys()) | set(condition2.flux_data.keys())

    for rxn_id in all_reactions:
        val1 = abs(condition1.flux_data.get(rxn_id, 0)) + pseudocount
        val2 = abs(condition2.flux_data.get(rxn_id, 0)) + pseudocount

        try:
            log2fc[rxn_id] = math.log2(val2 / val1)
        except (ValueError, ZeroDivisionError):
            log2fc[rxn_id] = 0.0

    return log2fc


class FluxDataDialog(QDialog):
    """Dialog for loading and visualizing external flux data."""

    flux_applied = Signal()  # Emitted when flux data is applied

    def __init__(self, appdata: AppData, central_widget=None):
        QDialog.__init__(self)
        self.setWindowTitle("Load External Flux Data")
        self.appdata = appdata
        self.central_widget = central_widget
        self.setMinimumSize(900, 700)

        self.conditions: List[FluxDataCondition] = []
        self._setup_ui()

    def _setup_ui(self):
        """Setup the dialog UI."""
        main_layout = QVBoxLayout()

        # Instructions
        instructions = QLabel(
            "Load flux data from CSV or TSV files to visualize on the map. "
            "You can load multiple conditions and calculate log2 fold change between them."
        )
        instructions.setWordWrap(True)
        main_layout.addWidget(instructions)

        # File loading section
        file_group = QGroupBox("Load Flux Data Files")
        file_layout = QVBoxLayout()

        # File format options
        format_layout = QGridLayout()

        format_layout.addWidget(QLabel("Delimiter:"), 0, 0)
        self.delimiter_combo = QComboBox()
        self.delimiter_combo.addItem("Comma (,)", ",")
        self.delimiter_combo.addItem("Tab (\\t)", "\t")
        self.delimiter_combo.addItem("Semicolon (;)", ";")
        format_layout.addWidget(self.delimiter_combo, 0, 1)

        format_layout.addWidget(QLabel("Reaction ID column:"), 0, 2)
        self.rxn_col_spin = QSpinBox()
        self.rxn_col_spin.setRange(0, 100)
        self.rxn_col_spin.setValue(0)
        format_layout.addWidget(self.rxn_col_spin, 0, 3)

        format_layout.addWidget(QLabel("Flux value column:"), 1, 0)
        self.flux_col_spin = QSpinBox()
        self.flux_col_spin.setRange(0, 100)
        self.flux_col_spin.setValue(1)
        format_layout.addWidget(self.flux_col_spin, 1, 1)

        self.has_header_check = QCheckBox("File has header row")
        self.has_header_check.setChecked(True)
        format_layout.addWidget(self.has_header_check, 1, 2, 1, 2)

        file_layout.addLayout(format_layout)

        # Condition name input
        name_layout = QHBoxLayout()
        name_layout.addWidget(QLabel("Condition name:"))
        self.condition_name_edit = QLineEdit()
        self.condition_name_edit.setPlaceholderText("e.g., Control, Treatment, WT, KO...")
        name_layout.addWidget(self.condition_name_edit)
        file_layout.addLayout(name_layout)

        # Load button
        btn_layout = QHBoxLayout()
        self.load_file_btn = QPushButton("Load File...")
        self.load_file_btn.clicked.connect(self._load_file)
        btn_layout.addWidget(self.load_file_btn)
        btn_layout.addStretch()
        file_layout.addLayout(btn_layout)

        file_group.setLayout(file_layout)
        main_layout.addWidget(file_group)

        # Loaded conditions section
        conditions_group = QGroupBox("Loaded Conditions")
        conditions_layout = QVBoxLayout()

        self.conditions_table = QTableWidget()
        self.conditions_table.setColumnCount(4)
        self.conditions_table.setHorizontalHeaderLabels(["Condition Name", "File", "# Reactions", "Actions"])
        self.conditions_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.conditions_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.conditions_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        conditions_layout.addWidget(self.conditions_table)

        conditions_group.setLayout(conditions_layout)
        main_layout.addWidget(conditions_group)

        # Visualization options
        viz_group = QGroupBox("Visualization Options")
        viz_layout = QVBoxLayout()

        # Mode selection
        mode_layout = QHBoxLayout()
        self.mode_group = QButtonGroup(self)

        self.single_mode_radio = QRadioButton("Single condition")
        self.single_mode_radio.setChecked(True)
        self.mode_group.addButton(self.single_mode_radio)
        mode_layout.addWidget(self.single_mode_radio)

        self.log2fc_mode_radio = QRadioButton("Log2 Fold Change")
        self.mode_group.addButton(self.log2fc_mode_radio)
        mode_layout.addWidget(self.log2fc_mode_radio)

        mode_layout.addStretch()
        viz_layout.addLayout(mode_layout)

        # Single condition selection
        self.single_widget = QWidget()
        single_layout = QHBoxLayout()
        single_layout.setContentsMargins(0, 0, 0, 0)
        single_layout.addWidget(QLabel("Select condition:"))
        self.single_condition_combo = QComboBox()
        single_layout.addWidget(self.single_condition_combo)
        single_layout.addStretch()
        self.single_widget.setLayout(single_layout)
        viz_layout.addWidget(self.single_widget)

        # Log2 fold change options
        self.log2fc_widget = QWidget()
        log2fc_layout = QGridLayout()
        log2fc_layout.setContentsMargins(0, 0, 0, 0)

        log2fc_layout.addWidget(QLabel("Reference (denominator):"), 0, 0)
        self.ref_condition_combo = QComboBox()
        log2fc_layout.addWidget(self.ref_condition_combo, 0, 1)

        log2fc_layout.addWidget(QLabel("Target (numerator):"), 0, 2)
        self.target_condition_combo = QComboBox()
        log2fc_layout.addWidget(self.target_condition_combo, 0, 3)

        log2fc_layout.addWidget(QLabel("Pseudocount:"), 1, 0)
        self.pseudocount_spin = QDoubleSpinBox()
        self.pseudocount_spin.setDecimals(10)
        self.pseudocount_spin.setRange(0, 1)
        self.pseudocount_spin.setValue(1e-6)
        self.pseudocount_spin.setSingleStep(1e-6)
        log2fc_layout.addWidget(self.pseudocount_spin, 1, 1)

        log2fc_layout.addWidget(QLabel("log2FC = log2(Target / Reference)"), 1, 2, 1, 2)

        self.log2fc_widget.setLayout(log2fc_layout)
        self.log2fc_widget.setVisible(False)
        viz_layout.addWidget(self.log2fc_widget)

        # Coloring options for log2FC
        self.color_options_widget = QWidget()
        color_layout = QHBoxLayout()
        color_layout.setContentsMargins(0, 0, 0, 0)

        color_layout.addWidget(QLabel("Color scale range (symmetric):"))
        self.color_range_spin = QDoubleSpinBox()
        self.color_range_spin.setRange(0.1, 10)
        self.color_range_spin.setValue(2.0)
        self.color_range_spin.setSingleStep(0.5)
        color_layout.addWidget(self.color_range_spin)
        color_layout.addWidget(QLabel("(±log2FC for full color saturation)"))
        color_layout.addStretch()

        self.color_options_widget.setLayout(color_layout)
        self.color_options_widget.setVisible(False)
        viz_layout.addWidget(self.color_options_widget)

        viz_group.setLayout(viz_layout)
        main_layout.addWidget(viz_group)

        # Connect mode radio buttons
        self.single_mode_radio.toggled.connect(self._update_mode_visibility)
        self.log2fc_mode_radio.toggled.connect(self._update_mode_visibility)

        # Preview section
        preview_group = QGroupBox("Preview / Statistics")
        preview_layout = QVBoxLayout()
        self.preview_text = QTextEdit()
        self.preview_text.setReadOnly(True)
        self.preview_text.setMaximumHeight(150)
        preview_layout.addWidget(self.preview_text)
        preview_group.setLayout(preview_layout)
        main_layout.addWidget(preview_group)

        # Action buttons
        btn_layout = QHBoxLayout()

        self.preview_btn = QPushButton("Preview")
        self.preview_btn.clicked.connect(self._preview)
        btn_layout.addWidget(self.preview_btn)

        self.apply_btn = QPushButton("Apply to Map")
        self.apply_btn.clicked.connect(self._apply)
        btn_layout.addWidget(self.apply_btn)

        self.clear_btn = QPushButton("Clear Visualization")
        self.clear_btn.clicked.connect(self._clear_visualization)
        btn_layout.addWidget(self.clear_btn)

        btn_layout.addStretch()

        self.close_btn = QPushButton("Close")
        self.close_btn.clicked.connect(self.accept)
        btn_layout.addWidget(self.close_btn)

        main_layout.addLayout(btn_layout)

        self.setLayout(main_layout)

    @Slot()
    def _update_mode_visibility(self):
        """Update visibility of mode-specific widgets."""
        is_log2fc = self.log2fc_mode_radio.isChecked()
        self.single_widget.setVisible(not is_log2fc)
        self.log2fc_widget.setVisible(is_log2fc)
        self.color_options_widget.setVisible(is_log2fc)

    @Slot()
    def _load_file(self):
        """Load a flux data file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Open Flux Data File",
            self.appdata.work_directory,
            "Data files (*.csv *.tsv *.txt);;CSV files (*.csv);;TSV files (*.tsv);;All files (*.*)"
        )

        if not file_path:
            return

        # Get condition name
        condition_name = self.condition_name_edit.text().strip()
        if not condition_name:
            condition_name = os.path.splitext(os.path.basename(file_path))[0]

        # Check for duplicate names
        existing_names = [c.name for c in self.conditions]
        if condition_name in existing_names:
            base_name = condition_name
            counter = 1
            while condition_name in existing_names:
                condition_name = f"{base_name}_{counter}"
                counter += 1

        # Create and load condition
        condition = FluxDataCondition(condition_name, file_path)

        delimiter = self.delimiter_combo.currentData()
        rxn_col = self.rxn_col_spin.value()
        flux_col = self.flux_col_spin.value()
        has_header = self.has_header_check.isChecked()

        success, message = condition.load_from_file(
            file_path, rxn_col, flux_col, has_header, delimiter
        )

        if success:
            self.conditions.append(condition)
            self._update_conditions_table()
            self._update_condition_combos()
            self.preview_text.setText(message)
            self.condition_name_edit.clear()
        else:
            QMessageBox.warning(self, "Load Error", message)

    def _update_conditions_table(self):
        """Update the conditions table."""
        self.conditions_table.setRowCount(len(self.conditions))

        for i, cond in enumerate(self.conditions):
            self.conditions_table.setItem(i, 0, QTableWidgetItem(cond.name))
            self.conditions_table.setItem(i, 1, QTableWidgetItem(os.path.basename(cond.file_path)))
            self.conditions_table.setItem(i, 2, QTableWidgetItem(str(len(cond.flux_data))))

            # Remove button
            remove_btn = QPushButton("Remove")
            remove_btn.clicked.connect(lambda checked, idx=i: self._remove_condition(idx))
            self.conditions_table.setCellWidget(i, 3, remove_btn)

    def _remove_condition(self, index: int):
        """Remove a condition."""
        if 0 <= index < len(self.conditions):
            del self.conditions[index]
            self._update_conditions_table()
            self._update_condition_combos()

    def _update_condition_combos(self):
        """Update all condition combo boxes."""
        current_single = self.single_condition_combo.currentText()
        current_ref = self.ref_condition_combo.currentText()
        current_target = self.target_condition_combo.currentText()

        self.single_condition_combo.clear()
        self.ref_condition_combo.clear()
        self.target_condition_combo.clear()

        for cond in self.conditions:
            self.single_condition_combo.addItem(cond.name)
            self.ref_condition_combo.addItem(cond.name)
            self.target_condition_combo.addItem(cond.name)

        # Restore selections if possible
        idx = self.single_condition_combo.findText(current_single)
        if idx >= 0:
            self.single_condition_combo.setCurrentIndex(idx)

        idx = self.ref_condition_combo.findText(current_ref)
        if idx >= 0:
            self.ref_condition_combo.setCurrentIndex(idx)

        idx = self.target_condition_combo.findText(current_target)
        if idx >= 0:
            self.target_condition_combo.setCurrentIndex(idx)

    def _get_selected_condition(self, name: str) -> Optional[FluxDataCondition]:
        """Get condition by name."""
        for cond in self.conditions:
            if cond.name == name:
                return cond
        return None

    @Slot()
    def _preview(self):
        """Preview the selected visualization."""
        if not self.conditions:
            self.preview_text.setText("No conditions loaded.")
            return

        if self.single_mode_radio.isChecked():
            cond = self._get_selected_condition(self.single_condition_combo.currentText())
            if cond:
                # Match with model reactions
                model = self.appdata.project.cobra_py_model
                model_rxn_ids = {rxn.id for rxn in model.reactions}
                matched = set(cond.flux_data.keys()) & model_rxn_ids
                unmatched = set(cond.flux_data.keys()) - model_rxn_ids

                stats = []
                stats.append(f"Condition: {cond.name}")
                stats.append(f"Total reactions in file: {len(cond.flux_data)}")
                stats.append(f"Matched to model: {len(matched)}")
                stats.append(f"Unmatched: {len(unmatched)}")
                if cond.flux_data:
                    vals = list(cond.flux_data.values())
                    stats.append(f"Flux range: [{min(vals):.4f}, {max(vals):.4f}]")
                    stats.append(f"Mean flux: {np.mean(vals):.4f}")

                self.preview_text.setText("\n".join(stats))
            else:
                self.preview_text.setText("No condition selected.")

        elif self.log2fc_mode_radio.isChecked():
            ref_cond = self._get_selected_condition(self.ref_condition_combo.currentText())
            target_cond = self._get_selected_condition(self.target_condition_combo.currentText())

            if ref_cond and target_cond:
                pseudocount = self.pseudocount_spin.value()
                log2fc = calculate_log2_fold_change(ref_cond, target_cond, pseudocount)

                # Match with model reactions
                model = self.appdata.project.cobra_py_model
                model_rxn_ids = {rxn.id for rxn in model.reactions}
                matched = set(log2fc.keys()) & model_rxn_ids

                stats = []
                stats.append(f"Log2 Fold Change: {target_cond.name} / {ref_cond.name}")
                stats.append(f"Total reactions: {len(log2fc)}")
                stats.append(f"Matched to model: {len(matched)}")

                if log2fc:
                    vals = list(log2fc.values())
                    stats.append(f"log2FC range: [{min(vals):.4f}, {max(vals):.4f}]")
                    stats.append(f"Mean log2FC: {np.mean(vals):.4f}")
                    stats.append(f"Up-regulated (log2FC > 0): {sum(1 for v in vals if v > 0)}")
                    stats.append(f"Down-regulated (log2FC < 0): {sum(1 for v in vals if v < 0)}")

                self.preview_text.setText("\n".join(stats))
            else:
                self.preview_text.setText("Please select both reference and target conditions.")

    @Slot()
    def _apply(self):
        """Apply flux data to map visualization."""
        if not self.conditions:
            QMessageBox.warning(self, "No Data", "Please load at least one condition first.")
            return

        model = self.appdata.project.cobra_py_model

        if self.single_mode_radio.isChecked():
            cond = self._get_selected_condition(self.single_condition_combo.currentText())
            if not cond:
                QMessageBox.warning(self, "No Condition", "Please select a condition.")
                return

            # Apply single condition flux values
            self.appdata.project.comp_values.clear()
            for rxn_id, flux_val in cond.flux_data.items():
                if rxn_id in [r.id for r in model.reactions]:
                    self.appdata.project.comp_values[rxn_id] = (flux_val, flux_val)

            self._update_visualization()
            self.preview_text.setText(
                f"Applied {len(self.appdata.project.comp_values)} flux values from '{cond.name}' to map."
            )

        elif self.log2fc_mode_radio.isChecked():
            ref_cond = self._get_selected_condition(self.ref_condition_combo.currentText())
            target_cond = self._get_selected_condition(self.target_condition_combo.currentText())

            if not ref_cond or not target_cond:
                QMessageBox.warning(self, "Missing Conditions",
                                    "Please select both reference and target conditions.")
                return

            if ref_cond.name == target_cond.name:
                QMessageBox.warning(self, "Same Condition",
                                    "Please select different conditions for comparison.")
                return

            pseudocount = self.pseudocount_spin.value()
            log2fc = calculate_log2_fold_change(ref_cond, target_cond, pseudocount)

            # Apply log2FC values as comp_values for visualization
            self.appdata.project.comp_values.clear()
            color_range = self.color_range_spin.value()

            for rxn_id, fc_val in log2fc.items():
                if rxn_id in [r.id for r in model.reactions]:
                    # Store as tuple for compatibility
                    self.appdata.project.comp_values[rxn_id] = (fc_val, fc_val)

            # Enable modes coloring for log2FC visualization
            self.appdata.modes_coloring = True

            self._update_visualization_log2fc(color_range)
            self.preview_text.setText(
                f"Applied log2 fold change ({target_cond.name} / {ref_cond.name}) to {len(self.appdata.project.comp_values)} reactions."
            )

    def _update_visualization(self):
        """Update map visualization with standard flux coloring."""
        self.appdata.modes_coloring = False

        if self.central_widget:
            self.central_widget.update()
            # Update all map tabs
            for i in range(self.central_widget.map_tabs.count()):
                widget = self.central_widget.map_tabs.widget(i)
                if hasattr(widget, 'update'):
                    widget.update()

    def _update_visualization_log2fc(self, color_range: float):
        """Update map visualization with log2FC coloring."""
        # For log2FC, we use a custom coloring scheme
        # This will be handled by overriding the color computation temporarily

        if self.central_widget:
            # Store original color computation
            original_flux_value_display = self.appdata.flux_value_display

            def log2fc_display(vl, vu):
                """Custom display for log2FC values."""
                fc_val = vl  # Both bounds are the same for log2FC

                # Scale to [-1, 1] based on color_range
                scaled = fc_val / color_range
                scaled = max(-1, min(1, scaled))

                # Color: red for negative (down), green for positive (up)
                if scaled > 0:
                    # Green for up-regulation
                    r = int(255 * (1 - scaled))
                    g = 255
                    b = int(255 * (1 - scaled))
                    background_color = QColor(r, g, b)
                elif scaled < 0:
                    # Red for down-regulation
                    r = 255
                    g = int(255 * (1 + scaled))
                    b = int(255 * (1 + scaled))
                    background_color = QColor(r, g, b)
                else:
                    # White for no change
                    background_color = QColor(255, 255, 255)

                flux_text = f"{fc_val:.3f}"
                return flux_text, background_color, True

            # Temporarily override
            self.appdata.flux_value_display = log2fc_display

            # Update visualization
            self.central_widget.update()
            for i in range(self.central_widget.map_tabs.count()):
                widget = self.central_widget.map_tabs.widget(i)
                if hasattr(widget, 'update'):
                    widget.update()

            # Restore original function
            self.appdata.flux_value_display = original_flux_value_display

    @Slot()
    def _clear_visualization(self):
        """Clear the current visualization."""
        self.appdata.project.comp_values.clear()
        self.appdata.modes_coloring = False

        if self.central_widget:
            self.central_widget.update()
            for i in range(self.central_widget.map_tabs.count()):
                widget = self.central_widget.map_tabs.widget(i)
                if hasattr(widget, 'update'):
                    widget.update()

        self.preview_text.setText("Visualization cleared.")


class FluxDataVisualizationManager:
    """Manager class to handle log2FC color state for persistent visualization."""

    def __init__(self, appdata: AppData):
        self.appdata = appdata
        self.log2fc_mode = False
        self.color_range = 2.0
        self.log2fc_values: Dict[str, float] = {}

    def set_log2fc_data(self, log2fc_values: Dict[str, float], color_range: float = 2.0):
        """Set log2FC data for visualization."""
        self.log2fc_mode = True
        self.log2fc_values = log2fc_values
        self.color_range = color_range

        # Apply to comp_values
        self.appdata.project.comp_values.clear()
        for rxn_id, fc_val in log2fc_values.items():
            self.appdata.project.comp_values[rxn_id] = (fc_val, fc_val)

    def clear(self):
        """Clear log2FC visualization."""
        self.log2fc_mode = False
        self.log2fc_values = {}
        self.appdata.project.comp_values.clear()

    def get_color_for_value(self, value: float) -> QColor:
        """Get color for a log2FC value."""
        scaled = value / self.color_range
        scaled = max(-1, min(1, scaled))

        if scaled > 0:
            r = int(255 * (1 - scaled))
            g = 255
            b = int(255 * (1 - scaled))
        elif scaled < 0:
            r = 255
            g = int(255 * (1 + scaled))
            b = int(255 * (1 + scaled))
        else:
            r, g, b = 255, 255, 255

        return QColor(r, g, b)
