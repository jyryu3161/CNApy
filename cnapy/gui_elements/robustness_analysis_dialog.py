"""Robustness Analysis Dialog for CNApy

This dialog performs robustness analysis by varying the flux of a target reaction
and observing the effect on the objective function.

This helps identify:
- Bottleneck reactions (where objective function drops sharply)
- Flux ranges where the system remains viable
- Sensitivity of the objective to specific reactions
"""

import csv

import numpy as np
from qtpy.QtCore import Qt, QThread, Signal, Slot
from qtpy.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QSlider,
    QSpinBox,
    QSplitter,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from cnapy.appdata import AppData
from cnapy.gui_elements.plot_customization_dialog import PlotCustomizationDialog
from cnapy.utils import no_scroll

# Check for matplotlib availability
try:
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure

    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    FigureCanvas = None
    Figure = None


class RobustnessWorkerThread(QThread):
    """Worker thread for running robustness analysis in background."""

    progress_update = Signal(int, int)  # (current, total)
    result_ready = Signal(object)  # dict with results
    error_occurred = Signal(str)

    def __init__(self, model, appdata: AppData, rxn_id: str, x_min: float, x_max: float, steps: int):
        """Initialize the worker thread.

        Args:
            model: COBRApy model (copy)
            appdata: Application data for scenario loading
            rxn_id: Target reaction ID
            x_min: Minimum flux value
            x_max: Maximum flux value
            steps: Number of steps
        """
        super().__init__()
        self.model = model
        self.appdata = appdata
        self.rxn_id = rxn_id
        self.x_min = x_min
        self.x_max = x_max
        self.steps = steps
        self._cancel_requested = False

    def request_cancel(self):
        """Request cancellation of the analysis."""
        self._cancel_requested = True

    def run(self):
        """Execute the robustness analysis."""
        try:
            model = self.model

            # Load scenario into model
            if self.appdata and self.appdata.project:
                self.appdata.project.load_scenario_into_model(model)

            # Check if reaction exists
            if self.rxn_id not in model.reactions:
                self.error_occurred.emit(f"Reaction '{self.rxn_id}' not found in model.")
                return

            # Get the target reaction
            target_rxn = model.reactions.get_by_id(self.rxn_id)
            original_lb = target_rxn.lower_bound
            original_ub = target_rxn.upper_bound

            # Generate flux values to test
            xs = np.linspace(self.x_min, self.x_max, self.steps)
            ys = np.full_like(xs, np.nan)
            statuses = []
            flux_distributions = []

            self.progress_update.emit(0, self.steps)

            for i, x in enumerate(xs):
                if self._cancel_requested:
                    return

                with model as m:
                    # Load scenario
                    if self.appdata and self.appdata.project:
                        self.appdata.project.load_scenario_into_model(m)

                    # Fix the target reaction flux
                    rxn = m.reactions.get_by_id(self.rxn_id)
                    rxn.lower_bound = x
                    rxn.upper_bound = x

                    # Optimize
                    try:
                        sol = m.optimize()
                        if sol.status == "optimal":
                            ys[i] = sol.objective_value
                            statuses.append("optimal")
                            # Store flux distribution
                            flux_distributions.append({r.id: sol.fluxes[r.id] for r in m.reactions})
                        else:
                            statuses.append(sol.status)
                            flux_distributions.append(None)
                    except Exception:
                        statuses.append("error")
                        flux_distributions.append(None)

                self.progress_update.emit(i + 1, self.steps)

            # Analyze bottleneck
            bottleneck_info = self._analyze_bottleneck(xs, ys)

            self.result_ready.emit(
                {
                    "rxn_id": self.rxn_id,
                    "xs": xs,
                    "ys": ys,
                    "statuses": statuses,
                    "flux_distributions": flux_distributions,
                    "original_bounds": (original_lb, original_ub),
                    "bottleneck": bottleneck_info,
                }
            )

        except Exception as e:
            self.error_occurred.emit(f"Analysis failed: {str(e)}")

    def _analyze_bottleneck(self, xs: np.ndarray, ys: np.ndarray) -> dict:
        """Analyze bottleneck in the robustness curve.

        Args:
            xs: Flux values
            ys: Objective values

        Returns:
            dict with bottleneck analysis results
        """
        # Remove NaN values for analysis
        valid_mask = ~np.isnan(ys)
        xs_valid = xs[valid_mask]
        ys_valid = ys[valid_mask]

        if len(xs_valid) < 2:
            return {
                "found": False,
                "message": "Not enough valid data points for bottleneck analysis.",
            }

        # Calculate gradient
        dy_dx = np.diff(ys_valid) / np.diff(xs_valid)

        if len(dy_dx) == 0:
            return {
                "found": False,
                "message": "Cannot calculate gradient.",
            }

        # Find steepest decline (most negative gradient)
        min_gradient_idx = np.argmin(dy_dx)
        min_gradient = dy_dx[min_gradient_idx]

        # Find steepest increase (most positive gradient)
        max_gradient_idx = np.argmax(dy_dx)
        max_gradient = dy_dx[max_gradient_idx]

        # Determine if there's a significant bottleneck
        gradient_range = max_gradient - min_gradient
        if gradient_range < 1e-6:
            return {
                "found": False,
                "message": "Objective is relatively insensitive to this reaction flux.",
            }

        # Get bottleneck location
        bottleneck_x = (xs_valid[min_gradient_idx] + xs_valid[min_gradient_idx + 1]) / 2
        bottleneck_y = (ys_valid[min_gradient_idx] + ys_valid[min_gradient_idx + 1]) / 2

        return {
            "found": True,
            "bottleneck_flux": bottleneck_x,
            "bottleneck_objective": bottleneck_y,
            "steepest_decline_gradient": min_gradient,
            "steepest_decline_range": (xs_valid[min_gradient_idx], xs_valid[min_gradient_idx + 1]),
            "max_objective": np.nanmax(ys),
            "max_objective_flux": xs_valid[np.argmax(ys_valid)],
            "sensitivity": abs(min_gradient),
        }


class RobustnessAnalysisDialog(QDialog):
    """Dialog for robustness analysis."""

    def __init__(self, appdata: AppData, central_widget=None):
        super().__init__()
        self.setWindowTitle("Robustness Analysis")
        self.setMinimumSize(1000, 700)

        self.appdata = appdata
        self.central_widget = central_widget
        self.worker_thread: RobustnessWorkerThread | None = None
        self.last_results: dict | None = None

        self._setup_ui()
        self._populate_reactions()

    def _setup_ui(self):
        """Setup the dialog UI."""
        main_layout = QHBoxLayout()

        # Left panel: Parameters
        left_widget = QWidget()
        left_layout = QVBoxLayout()

        # Description - simple and concise
        desc_label = QLabel("Analyze how the objective changes when varying a reaction's flux.")
        desc_label.setWordWrap(True)
        left_layout.addWidget(desc_label)

        # Reaction selection
        rxn_group = QGroupBox("Target Reaction")
        rxn_layout = QVBoxLayout()

        self.rxn_combo = no_scroll(QComboBox())
        self.rxn_combo.setEditable(True)
        self.rxn_combo.setInsertPolicy(QComboBox.NoInsert)
        self.rxn_combo.setToolTip(
            "Select a reaction to analyze.\n"
            "Typically exchange reactions (starting with EX_) are used.\n"
            "Examples: EX_glc__D_e (glucose uptake), EX_o2_e (oxygen uptake)"
        )
        self.rxn_combo.currentTextChanged.connect(self._on_reaction_changed)
        rxn_layout.addWidget(self.rxn_combo)

        # Reaction info
        self.rxn_info_label = QLabel("")
        self.rxn_info_label.setWordWrap(True)
        self.rxn_info_label.setStyleSheet("color: gray; font-size: 10px;")
        rxn_layout.addWidget(self.rxn_info_label)

        rxn_group.setLayout(rxn_layout)
        left_layout.addWidget(rxn_group)

        # Flux range
        range_group = QGroupBox("Flux Range")
        range_layout = QVBoxLayout()

        min_row = QHBoxLayout()
        min_row.addWidget(QLabel("Min flux:"))
        self.min_spin = no_scroll(QDoubleSpinBox())
        self.min_spin.setRange(-1e7, 1e7)
        self.min_spin.setValue(-10)
        self.min_spin.setDecimals(3)
        self.min_spin.setToolTip(
            "Minimum value of the flux range to analyze.\n"
            "For exchange reactions:\n"
            "• Negative (-): Uptake into the model\n"
            "• Positive (+): Secretion from the model"
        )
        min_row.addWidget(self.min_spin)
        range_layout.addLayout(min_row)

        max_row = QHBoxLayout()
        max_row.addWidget(QLabel("Max flux:"))
        self.max_spin = no_scroll(QDoubleSpinBox())
        self.max_spin.setRange(-1e7, 1e7)
        self.max_spin.setValue(0)
        self.max_spin.setDecimals(3)
        self.max_spin.setToolTip("Maximum value of the flux range to analyze.\n" "Must be greater than Min.")
        max_row.addWidget(self.max_spin)
        range_layout.addLayout(max_row)

        steps_row = QHBoxLayout()
        steps_row.addWidget(QLabel("Steps:"))
        self.steps_spin = no_scroll(QSpinBox())
        self.steps_spin.setRange(5, 500)
        self.steps_spin.setValue(31)
        self.steps_spin.setToolTip(
            "Number of points to sample in the flux range.\n"
            "• Higher values give finer resolution\n"
            "• Higher values take longer to compute\n"
            "• Recommended: 20-50 (quick), 100+ (detailed)"
        )
        steps_row.addWidget(self.steps_spin)
        range_layout.addLayout(steps_row)

        # Use scenario bounds checkbox
        self.use_scenario_bounds_cb = QCheckBox("Use scenario bounds as range")
        self.use_scenario_bounds_cb.setChecked(False)
        self.use_scenario_bounds_cb.setToolTip(
            "When checked, automatically fills Min/Max\n" "with the bounds from the current scenario."
        )
        self.use_scenario_bounds_cb.toggled.connect(self._on_use_scenario_bounds)
        range_layout.addWidget(self.use_scenario_bounds_cb)

        range_group.setLayout(range_layout)
        left_layout.addWidget(range_group)

        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        left_layout.addWidget(self.progress_bar)

        # Run button
        btn_layout = QHBoxLayout()
        self.run_btn = QPushButton("Run Analysis")
        self.run_btn.setStyleSheet("font-size: 14px; padding: 10px; background-color: #4CAF50; color: white;")
        self.run_btn.clicked.connect(self._run_analysis)
        btn_layout.addWidget(self.run_btn)

        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setEnabled(False)
        self.cancel_btn.clicked.connect(self._cancel_analysis)
        btn_layout.addWidget(self.cancel_btn)

        left_layout.addLayout(btn_layout)

        # Bottleneck info
        bottleneck_group = QGroupBox("Bottleneck Analysis")
        bottleneck_layout = QVBoxLayout()

        self.bottleneck_text = QTextEdit()
        self.bottleneck_text.setReadOnly(True)
        self.bottleneck_text.setMaximumHeight(150)
        self.bottleneck_text.setPlaceholderText("Run analysis to see bottleneck information...")
        bottleneck_layout.addWidget(self.bottleneck_text)

        bottleneck_group.setLayout(bottleneck_layout)
        left_layout.addWidget(bottleneck_group)

        # Map slider
        slider_group = QGroupBox("Apply to Map")
        slider_layout = QVBoxLayout()

        self.flux_slider = QSlider(Qt.Horizontal)
        self.flux_slider.setEnabled(False)
        self.flux_slider.valueChanged.connect(self._on_slider_changed)
        slider_layout.addWidget(self.flux_slider)

        self.slider_value_label = QLabel("Flux: --")
        slider_layout.addWidget(self.slider_value_label)

        self.apply_btn = QPushButton("Apply to Map")
        self.apply_btn.setEnabled(False)
        self.apply_btn.clicked.connect(self._apply_to_map)
        slider_layout.addWidget(self.apply_btn)

        slider_group.setLayout(slider_layout)
        left_layout.addWidget(slider_group)

        left_layout.addStretch()

        # Close and export buttons
        bottom_btn_layout = QHBoxLayout()

        export_btn = QPushButton("Export...")
        export_btn.clicked.connect(self._export_results)
        bottom_btn_layout.addWidget(export_btn)

        bottom_btn_layout.addStretch()

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.accept)
        bottom_btn_layout.addWidget(close_btn)

        left_layout.addLayout(bottom_btn_layout)

        left_widget.setLayout(left_layout)
        left_widget.setMaximumWidth(350)

        # Right panel: Plot
        right_widget = QWidget()
        right_layout = QVBoxLayout()

        if MATPLOTLIB_AVAILABLE:
            self.figure = Figure(figsize=(8, 6))
            self.canvas = FigureCanvas(self.figure)
            right_layout.addWidget(self.canvas)

            # Plot controls
            plot_controls = QHBoxLayout()
            self.customize_btn = QPushButton("Customize Plot")
            self.customize_btn.setEnabled(False)
            self.customize_btn.clicked.connect(self._customize_plot)
            plot_controls.addWidget(self.customize_btn)
            save_plot_btn = QPushButton("Save Plot...")
            save_plot_btn.clicked.connect(self._save_plot)
            plot_controls.addWidget(save_plot_btn)
            plot_controls.addStretch()
            right_layout.addLayout(plot_controls)
        else:
            no_plot_label = QLabel(
                "Matplotlib not available.\n" "Install with: pip install matplotlib\n" "to see the robustness plot."
            )
            no_plot_label.setAlignment(Qt.AlignCenter)
            right_layout.addWidget(no_plot_label)

        right_widget.setLayout(right_layout)

        # Splitter
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(left_widget)
        splitter.addWidget(right_widget)
        splitter.setSizes([350, 650])

        main_layout.addWidget(splitter)
        self.setLayout(main_layout)

    def _populate_reactions(self):
        """Populate the reaction combo box."""
        model = self.appdata.project.cobra_py_model

        # Add all reactions
        for rxn in sorted(model.reactions, key=lambda r: r.id):
            self.rxn_combo.addItem(rxn.id)

        # Try to select a common exchange reaction
        common_exchanges = ["EX_glc__D_e", "EX_glc_D_e", "EX_o2_e"]
        for ex in common_exchanges:
            if ex in model.reactions:
                self.rxn_combo.setCurrentText(ex)
                break

    @Slot(str)
    def _on_reaction_changed(self, rxn_id: str):
        """Handle reaction selection change."""
        model = self.appdata.project.cobra_py_model

        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            info = f"Name: {rxn.name or 'N/A'}\n" f"Bounds: [{rxn.lower_bound}, {rxn.upper_bound}]\n" f"{rxn.reaction}"
            self.rxn_info_label.setText(info)

            # Update min/max spin boxes if using scenario bounds
            if self.use_scenario_bounds_cb.isChecked():
                self._update_bounds_from_scenario()
        else:
            self.rxn_info_label.setText("")

    @Slot(bool)
    def _on_use_scenario_bounds(self, checked: bool):
        """Handle use scenario bounds checkbox toggle."""
        if checked:
            self._update_bounds_from_scenario()

    def _update_bounds_from_scenario(self):
        """Update flux range from scenario bounds."""
        rxn_id = self.rxn_combo.currentText()
        model = self.appdata.project.cobra_py_model

        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)

            # Check if there are scenario constraints
            lb = rxn.lower_bound
            ub = rxn.upper_bound

            # Look for scenario-specific bounds
            if rxn_id in self.appdata.project.scen_values:
                scen_lb, scen_ub = self.appdata.project.scen_values[rxn_id]
                if scen_lb is not None:
                    lb = scen_lb
                if scen_ub is not None:
                    ub = scen_ub

            # Replace non-finite bounds with the spinbox limits so an
            # effectively-unbounded reaction (e.g. +/-inf) is not silently
            # clamped to a misleading value.
            if not np.isfinite(lb):
                lb = self.min_spin.minimum() if lb < 0 else self.min_spin.maximum()
            if not np.isfinite(ub):
                ub = self.max_spin.minimum() if ub < 0 else self.max_spin.maximum()

            self.min_spin.setValue(lb)
            self.max_spin.setValue(ub)

    @Slot()
    def _run_analysis(self):
        """Start the robustness analysis."""
        rxn_id = self.rxn_combo.currentText()
        model = self.appdata.project.cobra_py_model

        if rxn_id not in model.reactions:
            QMessageBox.warning(self, "Invalid Reaction", f"Reaction '{rxn_id}' not found in the model.")
            return

        x_min = self.min_spin.value()
        x_max = self.max_spin.value()
        steps = self.steps_spin.value()

        if x_min >= x_max:
            QMessageBox.warning(self, "Invalid Range", "Min flux must be less than max flux.")
            return

        # Disable run, enable cancel
        self.run_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, steps)
        self.progress_bar.setValue(0)

        # Create model copy
        model_copy = model.copy()

        # Start worker thread
        self.worker_thread = RobustnessWorkerThread(model_copy, self.appdata, rxn_id, x_min, x_max, steps)
        self.worker_thread.progress_update.connect(self._on_progress)
        self.worker_thread.result_ready.connect(self._on_results)
        self.worker_thread.error_occurred.connect(self._on_error)
        self.worker_thread.finished.connect(self._on_finished)
        self.worker_thread.start()

    @Slot()
    def _cancel_analysis(self):
        """Cancel the running analysis."""
        if self.worker_thread and self.worker_thread.isRunning():
            self.worker_thread.request_cancel()
            self.cancel_btn.setEnabled(False)

    @Slot(int, int)
    def _on_progress(self, current: int, total: int):
        """Update progress bar."""
        self.progress_bar.setValue(current)

    @Slot(object)
    def _on_results(self, data: dict):
        """Handle analysis results."""
        self.last_results = data

        # Update plot
        if MATPLOTLIB_AVAILABLE:
            self._update_plot(data)

        # Update bottleneck info
        self._update_bottleneck_info(data["bottleneck"])

        # Enable slider
        self._setup_slider(data)

    def _update_plot(self, data: dict):
        """Update the matplotlib plot."""
        self.figure.clear()
        ax = self.figure.add_subplot(111)

        xs = data["xs"]
        ys = data["ys"]

        # Plot the main curve
        ax.plot(xs, ys, "b-", linewidth=2, label="Objective Value")

        # Mark valid and invalid points
        valid_mask = ~np.isnan(ys)
        ax.scatter(xs[valid_mask], ys[valid_mask], c="blue", s=20, zorder=5)

        # Mark infeasible points
        invalid_mask = np.isnan(ys)
        if np.any(invalid_mask):
            ax.scatter(
                xs[invalid_mask],
                np.zeros(np.sum(invalid_mask)),
                c="red",
                marker="x",
                s=50,
                label="Infeasible",
                zorder=5,
            )

        # Mark bottleneck if found
        bottleneck = data["bottleneck"]
        if bottleneck.get("found", False):
            ax.axvline(bottleneck["bottleneck_flux"], color="orange", linestyle="--", linewidth=1.5, label="Bottleneck")
            ax.scatter(
                [bottleneck["bottleneck_flux"]], [bottleneck["bottleneck_objective"]], c="orange", s=100, zorder=6
            )

        # Labels and title
        ax.set_xlabel(f"Flux of {data['rxn_id']}")
        ax.set_ylabel("Objective Value")
        ax.set_title(f"Robustness Analysis: {data['rxn_id']}")
        ax.legend(loc="best")
        ax.grid(True, alpha=0.3)

        self.figure.tight_layout()
        self.canvas.draw()
        self.customize_btn.setEnabled(True)

    def _update_bottleneck_info(self, bottleneck: dict):
        """Update the bottleneck information text."""
        if not bottleneck.get("found", False):
            info_html = (
                "<b>Bottleneck Analysis Results</b><br><br>"
                "<b>Bottleneck Found:</b> No<br><br>"
                f"{bottleneck.get('message', 'No sharp change point was found in the analyzed range.')}<br><br>"
                "The objective function may change linearly with flux,<br>"
                "or the changes may be minimal."
            )
            self.bottleneck_text.setHtml(info_html)
            return

        info_html = (
            f"<b>Bottleneck Analysis Results</b><br><br>"
            f"<b>Bottleneck Found:</b> Yes<br>"
            f"<b>Location:</b> flux = {bottleneck['bottleneck_flux']:.4f}<br>"
            f"&nbsp;&nbsp;→ The objective function changes sharply at this point.<br><br>"
            f"<b>Objective at Bottleneck:</b> {bottleneck['bottleneck_objective']:.4f}<br>"
            f"<b>Maximum Objective:</b> {bottleneck['max_objective']:.4f} "
            f"(at flux = {bottleneck['max_objective_flux']:.4f})<br><br>"
            f"<b>Sensitivity:</b> {bottleneck['sensitivity']:.4f}<br>"
            f"&nbsp;&nbsp;→ Higher values indicate greater sensitivity to flux changes."
        )
        self.bottleneck_text.setHtml(info_html)

    def _setup_slider(self, data: dict):
        """Setup the flux slider based on results."""
        xs = data["xs"]
        self.flux_slider.setRange(0, len(xs) - 1)
        self.flux_slider.setValue(len(xs) // 2)
        self.flux_slider.setEnabled(True)
        self.apply_btn.setEnabled(True)
        self._on_slider_changed(self.flux_slider.value())

    @Slot(int)
    def _on_slider_changed(self, value: int):
        """Handle slider value change."""
        if not self.last_results:
            return

        xs = self.last_results["xs"]
        ys = self.last_results["ys"]

        flux = xs[value]
        obj = ys[value] if not np.isnan(ys[value]) else "Infeasible"

        if isinstance(obj, float):
            self.slider_value_label.setText(f"Flux: {flux:.4f}, Objective: {obj:.4f}")
        else:
            self.slider_value_label.setText(f"Flux: {flux:.4f}, Objective: {obj}")

    @Slot()
    def _apply_to_map(self):
        """Apply the selected flux distribution to the map."""
        if not self.last_results:
            return

        idx = self.flux_slider.value()
        flux_dist = self.last_results["flux_distributions"][idx]

        if flux_dist is None:
            QMessageBox.warning(self, "No Solution", "No flux distribution available at this point (infeasible).")
            return

        # Apply to comp_values
        self.appdata.project.comp_values.clear()
        for rxn_id, flux in flux_dist.items():
            self.appdata.project.comp_values[rxn_id] = (flux, flux)

        self.appdata.project.comp_values_type = 0

        # Update the map view
        if self.central_widget:
            self.central_widget.update()

        flux = self.last_results["xs"][idx]
        QMessageBox.information(
            self, "Applied", f"Flux distribution at {self.last_results['rxn_id']} = {flux:.4f} applied to map."
        )

    @Slot(str)
    def _on_error(self, error: str):
        """Handle analysis error."""
        QMessageBox.warning(self, "Analysis Error", error)

    @Slot()
    def _on_finished(self):
        """Handle worker thread completion."""
        self.run_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)
        self.progress_bar.setVisible(False)

    @Slot()
    def _customize_plot(self):
        """Open the plot customization dialog."""
        dialog = PlotCustomizationDialog(self, self.figure, self.canvas)
        dialog.exec()

    @Slot()
    def _save_plot(self):
        """Save the plot to file."""
        if not MATPLOTLIB_AVAILABLE:
            return

        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Save Plot",
            self.appdata.work_directory,
            "PNG files (*.png);;PDF files (*.pdf);;SVG files (*.svg)",
        )

        if filename:
            self.figure.savefig(filename, dpi=150, bbox_inches="tight")
            QMessageBox.information(self, "Saved", f"Plot saved to {filename}")

    @Slot()
    def _export_results(self):
        """Export results to CSV."""
        if not self.last_results:
            QMessageBox.warning(self, "No Results", "No results to export. Run analysis first.")
            return

        filename, _ = QFileDialog.getSaveFileName(
            self, "Export Results", self.appdata.work_directory, "CSV files (*.csv)"
        )

        if not filename:
            return

        try:
            with open(filename, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(["Flux", "Objective Value", "Status"])

                for i, (x, y, status) in enumerate(
                    zip(
                        self.last_results["xs"],
                        self.last_results["ys"],
                        self.last_results["statuses"],
                        strict=True,
                    )
                ):
                    writer.writerow([x, y if not np.isnan(y) else "", status])

            QMessageBox.information(self, "Exported", f"Results exported to {filename}")
        except Exception as e:
            QMessageBox.warning(self, "Export Error", f"Failed to export: {str(e)}")

    def closeEvent(self, event):
        """Handle dialog close."""
        if self.worker_thread and self.worker_thread.isRunning():
            self.worker_thread.request_cancel()
            self.worker_thread.wait(2000)
        event.accept()
