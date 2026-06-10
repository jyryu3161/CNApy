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

"""Dynamic FBA (dFBA) Dialog for CNApy

Dynamic Flux Balance Analysis simulates time-course behavior of metabolic
systems by coupling FBA with ordinary differential equations (ODEs).

The simulation tracks:
- Biomass concentration over time
- External metabolite concentrations (substrates, products)
- Flux distributions at each time point

References:
- Mahadevan et al. (2002). Dynamic Flux Balance Analysis of Diauxic Growth
  in Escherichia coli. Biophysical Journal.
- Varma & Palsson (1994). Stoichiometric flux balance models quantitatively
  predict growth and metabolic by-product secretion in wild-type E. coli W3110.
"""

from collections.abc import Callable
from dataclasses import dataclass, field

import cobra
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
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from cnapy.appdata import AppData
from cnapy.gui_elements.plot_customization_dialog import PlotCustomizationDialog
from cnapy.utils import no_scroll

# Check for scipy availability
try:
    from scipy.integrate import solve_ivp

    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    solve_ivp = None

# Check for matplotlib availability
try:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.figure import Figure

    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    plt = None
    FigureCanvas = None
    Figure = None


# ============================================================================
# dFBA Core Implementation
# ============================================================================


@dataclass
class DFBAResult:
    """Results from a dFBA simulation."""

    time: np.ndarray
    biomass: np.ndarray
    metabolites: dict[str, np.ndarray]  # metabolite_id -> concentration array
    fluxes: dict[str, np.ndarray]  # reaction_id -> flux array (sampled at time points)
    success: bool
    message: str


@dataclass
class DFBAParameters:
    """Parameters for dFBA simulation."""

    initial_biomass: float = 0.1  # g/L
    time_start: float = 0.0  # hours
    time_end: float = 10.0  # hours
    time_steps: int = 100  # number of output time points

    # Substrate parameters: {exchange_rxn_id: (initial_conc, Km, Vmax)}
    substrates: dict[str, tuple[float, float, float]] = field(default_factory=dict)

    # Products to track: {exchange_rxn_id: initial_conc}
    products: dict[str, float] = field(default_factory=dict)

    # Biomass reaction ID
    biomass_reaction: str = ""

    # ODE solver settings
    method: str = "RK45"  # ODE solver method
    max_step: float = 0.1  # Maximum step size


def run_dfba(
    model: cobra.Model, params: DFBAParameters, progress_callback: Callable[[float], None] | None = None
) -> DFBAResult:
    """
    Run Dynamic Flux Balance Analysis simulation.

    The dFBA couples FBA optimization with ODE integration:

    dX/dt = μ * X  (biomass growth)
    dS/dt = -v_uptake * X  (substrate consumption)
    dP/dt = v_production * X  (product formation)

    where:
    - X is biomass concentration
    - S is substrate concentration
    - P is product concentration
    - μ is specific growth rate (from FBA objective)
    - v_uptake follows Michaelis-Menten kinetics: Vmax * S / (Km + S)

    Parameters:
    -----------
    model : cobra.Model
        The metabolic model
    params : DFBAParameters
        Simulation parameters
    progress_callback : Callable, optional
        Function to report progress (0-100)

    Returns:
    --------
    DFBAResult
        Simulation results
    """
    if not SCIPY_AVAILABLE:
        return DFBAResult(
            time=np.array([]),
            biomass=np.array([]),
            metabolites={},
            fluxes={},
            success=False,
            message="SciPy is required for dFBA simulation. Install it with: pip install scipy",
        )

    if not params.biomass_reaction:
        return DFBAResult(
            time=np.array([]),
            biomass=np.array([]),
            metabolites={},
            fluxes={},
            success=False,
            message="No biomass reaction specified.",
        )

    # Get substrate and product reaction IDs
    substrate_ids = list(params.substrates.keys())
    product_ids = list(params.products.keys())

    # Initial conditions: [biomass, substrate1, substrate2, ..., product1, product2, ...]
    y0 = [params.initial_biomass]
    for sid in substrate_ids:
        y0.append(params.substrates[sid][0])  # initial concentration
    for pid in product_ids:
        y0.append(params.products[pid])  # initial concentration

    y0 = np.array(y0)

    # Time span
    t_span = (params.time_start, params.time_end)
    t_eval = np.linspace(params.time_start, params.time_end, params.time_steps)

    # Store flux history
    flux_history = {params.biomass_reaction: []}
    for rid in substrate_ids + product_ids:
        flux_history[rid] = []

    flux_time_points = []

    def dfba_rhs(t, y):
        """Right-hand side of the dFBA ODE system."""
        # Unpack state variables
        biomass = max(y[0], 1e-10)  # Prevent negative biomass

        substrates_conc = {}
        for i, sid in enumerate(substrate_ids):
            substrates_conc[sid] = max(y[1 + i], 0)  # Prevent negative concentrations

        products_conc = {}
        for i, pid in enumerate(product_ids):
            products_conc[pid] = max(y[1 + len(substrate_ids) + i], 0)

        # Update model bounds based on substrate availability (Michaelis-Menten)
        with model as m:
            for sid, (_, km, vmax) in params.substrates.items():
                if sid in m.reactions:
                    S = substrates_conc.get(sid, 0)
                    # Michaelis-Menten uptake rate
                    uptake_rate = vmax * S / (km + S) if (km + S) > 0 else 0
                    rxn = m.reactions.get_by_id(sid)
                    # Exchange reactions: negative flux = uptake
                    rxn.lower_bound = -uptake_rate

            # Set the objective to the selected biomass reaction so that the
            # growth rate reflects that reaction, regardless of the model's
            # (or loaded scenario's) current objective.
            if params.biomass_reaction in m.reactions:
                m.objective = m.reactions.get_by_id(params.biomass_reaction)

            # Perform FBA
            try:
                solution = m.optimize()
                if solution.status != "optimal":
                    # If infeasible, set zero growth
                    mu = 0
                    fluxes = {rid: 0 for rid in [params.biomass_reaction] + substrate_ids + product_ids}
                else:
                    # Growth rate is the flux through the selected biomass
                    # reaction, not the raw objective value.
                    mu = solution.fluxes[params.biomass_reaction]
                    fluxes = solution.fluxes
            except Exception:
                mu = 0
                fluxes = {rid: 0 for rid in [params.biomass_reaction] + substrate_ids + product_ids}

        # Record fluxes
        flux_time_points.append(t)
        flux_history[params.biomass_reaction].append(mu)
        for rid in substrate_ids + product_ids:
            flux_history[rid].append(fluxes.get(rid, 0) if isinstance(fluxes, dict) else fluxes[rid])

        # Compute derivatives
        dydt = np.zeros_like(y)

        # dX/dt = μ * X (biomass growth)
        dydt[0] = mu * biomass

        # dS/dt = v_exchange * X (substrate consumption)
        for i, sid in enumerate(substrate_ids):
            v_exchange = fluxes.get(sid, 0) if isinstance(fluxes, dict) else fluxes[sid]
            # Exchange flux is negative for uptake
            dydt[1 + i] = v_exchange * biomass

        # dP/dt = v_exchange * X (product formation)
        for i, pid in enumerate(product_ids):
            v_exchange = fluxes.get(pid, 0) if isinstance(fluxes, dict) else fluxes[pid]
            # Exchange flux is positive for secretion
            dydt[1 + len(substrate_ids) + i] = v_exchange * biomass

        return dydt

    # Define termination event (substrate depletion)
    def substrate_depleted(t, y):
        # Check if all substrates are below threshold
        for i, sid in enumerate(substrate_ids):
            if y[1 + i] > 0.01:  # At least one substrate remains
                return 1
        return -1  # All substrates depleted

    substrate_depleted.terminal = True
    substrate_depleted.direction = -1

    # Solve ODE system
    try:
        sol = solve_ivp(
            dfba_rhs,
            t_span,
            y0,
            method=params.method,
            t_eval=t_eval,
            events=[substrate_depleted] if substrate_ids else None,
            max_step=params.max_step,
            dense_output=True,
        )

        if not sol.success:
            return DFBAResult(
                time=np.array([]),
                biomass=np.array([]),
                metabolites={},
                fluxes={},
                success=False,
                message=f"ODE solver failed: {sol.message}",
            )

        # Extract results
        time = sol.t
        biomass = sol.y[0, :]

        metabolites = {}
        for i, sid in enumerate(substrate_ids):
            metabolites[sid] = sol.y[1 + i, :]
        for i, pid in enumerate(product_ids):
            metabolites[pid] = sol.y[1 + len(substrate_ids) + i, :]

        # Interpolate fluxes to match output time points
        fluxes = {}
        if flux_time_points:
            flux_time_arr = np.array(flux_time_points)
            for rid, flux_vals in flux_history.items():
                flux_arr = np.array(flux_vals)
                # Simple interpolation to output time points
                if len(flux_arr) > 1:
                    fluxes[rid] = np.interp(time, flux_time_arr, flux_arr)
                else:
                    fluxes[rid] = np.zeros_like(time)

        return DFBAResult(
            time=time,
            biomass=biomass,
            metabolites=metabolites,
            fluxes=fluxes,
            success=True,
            message="Simulation completed successfully.",
        )

    except Exception as e:
        return DFBAResult(
            time=np.array([]),
            biomass=np.array([]),
            metabolites={},
            fluxes={},
            success=False,
            message=f"Simulation error: {str(e)}",
        )


# ============================================================================
# dFBA Dialog
# ============================================================================


class DFBASimulationThread(QThread):
    """Thread for running dFBA simulation."""

    progress = Signal(float)
    finished = Signal(object)  # DFBAResult

    def __init__(self, model: cobra.Model, params: DFBAParameters):
        super().__init__()
        self.model = model.copy()
        self.params = params

    def run(self):
        result = run_dfba(self.model, self.params, self.progress.emit)
        self.finished.emit(result)


class DynamicFBADialog(QDialog):
    """Dialog for Dynamic FBA simulation and visualization."""

    def __init__(self, appdata: AppData, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Dynamic FBA (dFBA) Simulation")
        self.setMinimumSize(1200, 900)
        self.appdata = appdata

        self.simulation_thread: DFBASimulationThread | None = None
        self.last_result: DFBAResult | None = None

        self._setup_ui()
        self._populate_reactions()

    def _setup_ui(self):
        """Setup the dialog UI."""
        main_layout = QHBoxLayout()

        # Left panel: Parameters
        left_widget = QWidget()
        left_layout = QVBoxLayout()

        # Biomass reaction selection
        biomass_group = QGroupBox("Biomass Reaction")
        biomass_layout = QVBoxLayout()

        self.biomass_combo = no_scroll(QComboBox())
        self.biomass_combo.setEditable(True)
        biomass_layout.addWidget(QLabel("Select biomass/growth reaction:"))
        biomass_layout.addWidget(self.biomass_combo)

        biomass_group.setLayout(biomass_layout)
        left_layout.addWidget(biomass_group)

        # Initial conditions
        init_group = QGroupBox("Initial Conditions")
        init_layout = QVBoxLayout()

        biomass_row = QHBoxLayout()
        biomass_row.addWidget(QLabel("Initial biomass (g/L):"))
        self.initial_biomass_spin = no_scroll(QDoubleSpinBox())
        self.initial_biomass_spin.setRange(0.001, 100)
        self.initial_biomass_spin.setValue(0.1)
        self.initial_biomass_spin.setDecimals(4)
        biomass_row.addWidget(self.initial_biomass_spin)
        init_layout.addLayout(biomass_row)

        init_group.setLayout(init_layout)
        left_layout.addWidget(init_group)

        # Time settings
        time_group = QGroupBox("Time Settings")
        time_layout = QVBoxLayout()

        start_row = QHBoxLayout()
        start_row.addWidget(QLabel("Start time (h):"))
        self.time_start_spin = no_scroll(QDoubleSpinBox())
        self.time_start_spin.setRange(0, 1000)
        self.time_start_spin.setValue(0)
        start_row.addWidget(self.time_start_spin)
        time_layout.addLayout(start_row)

        end_row = QHBoxLayout()
        end_row.addWidget(QLabel("End time (h):"))
        self.time_end_spin = no_scroll(QDoubleSpinBox())
        self.time_end_spin.setRange(0.1, 1000)
        self.time_end_spin.setValue(10)
        end_row.addWidget(self.time_end_spin)
        time_layout.addLayout(end_row)

        steps_row = QHBoxLayout()
        steps_row.addWidget(QLabel("Time steps:"))
        self.time_steps_spin = no_scroll(QSpinBox())
        self.time_steps_spin.setRange(10, 10000)
        self.time_steps_spin.setValue(100)
        steps_row.addWidget(self.time_steps_spin)
        time_layout.addLayout(steps_row)

        time_group.setLayout(time_layout)
        left_layout.addWidget(time_group)

        # Substrates table
        substrates_group = QGroupBox("Substrates (Carbon/Nitrogen Sources)")
        substrates_layout = QVBoxLayout()

        self.substrates_table = QTableWidget()
        self.substrates_table.setColumnCount(5)
        self.substrates_table.setHorizontalHeaderLabels(
            ["Use", "EX Reaction", "Initial (mM)", "Km (mM)", "Vmax (mmol/gDW/h)"]
        )
        hdr_s = self.substrates_table.horizontalHeader()
        hdr_s.setSectionResizeMode(0, QHeaderView.Fixed)
        self.substrates_table.setColumnWidth(0, 36)
        hdr_s.setSectionResizeMode(1, QHeaderView.Stretch)
        hdr_s.setSectionResizeMode(2, QHeaderView.Fixed)
        self.substrates_table.setColumnWidth(2, 80)
        hdr_s.setSectionResizeMode(3, QHeaderView.Fixed)
        self.substrates_table.setColumnWidth(3, 65)
        hdr_s.setSectionResizeMode(4, QHeaderView.Fixed)
        self.substrates_table.setColumnWidth(4, 130)
        self.substrates_table.setMaximumHeight(220)
        self.substrates_table.horizontalHeader().sectionClicked.connect(self._on_substrates_header_clicked)
        substrates_layout.addWidget(self.substrates_table)

        substrates_btn_layout = QHBoxLayout()
        add_substrate_btn = QPushButton("Add")
        add_substrate_btn.clicked.connect(self._add_substrate_row)
        substrates_btn_layout.addWidget(add_substrate_btn)
        remove_substrate_btn = QPushButton("Remove")
        remove_substrate_btn.clicked.connect(self._remove_substrate_row)
        substrates_btn_layout.addWidget(remove_substrate_btn)
        substrates_btn_layout.addStretch()
        substrates_layout.addLayout(substrates_btn_layout)

        substrates_group.setLayout(substrates_layout)
        left_layout.addWidget(substrates_group)

        # Products table
        products_group = QGroupBox("Products to Track")
        products_layout = QVBoxLayout()

        self.products_table = QTableWidget()
        self.products_table.setColumnCount(3)
        self.products_table.setHorizontalHeaderLabels(["Use", "EX Reaction", "Initial (mM)"])
        hdr_p = self.products_table.horizontalHeader()
        hdr_p.setSectionResizeMode(0, QHeaderView.Fixed)
        self.products_table.setColumnWidth(0, 36)
        hdr_p.setSectionResizeMode(1, QHeaderView.Stretch)
        hdr_p.setSectionResizeMode(2, QHeaderView.Fixed)
        self.products_table.setColumnWidth(2, 90)
        self.products_table.setMaximumHeight(200)
        self.products_table.horizontalHeader().sectionClicked.connect(self._on_products_header_clicked)
        products_layout.addWidget(self.products_table)

        products_btn_layout = QHBoxLayout()
        add_product_btn = QPushButton("Add")
        add_product_btn.clicked.connect(self._add_product_row)
        products_btn_layout.addWidget(add_product_btn)
        remove_product_btn = QPushButton("Remove")
        remove_product_btn.clicked.connect(self._remove_product_row)
        products_btn_layout.addWidget(remove_product_btn)
        products_btn_layout.addStretch()
        products_layout.addLayout(products_btn_layout)

        products_group.setLayout(products_layout)
        left_layout.addWidget(products_group)

        # Solver settings
        solver_group = QGroupBox("Solver Settings")
        solver_layout = QVBoxLayout()

        method_row = QHBoxLayout()
        method_row.addWidget(QLabel("ODE Method:"))
        self.method_combo = no_scroll(QComboBox())
        self.method_combo.addItems(["RK45", "RK23", "DOP853", "Radau", "BDF", "LSODA"])
        self.method_combo.setCurrentText("RK45")
        method_row.addWidget(self.method_combo)
        solver_layout.addLayout(method_row)

        maxstep_row = QHBoxLayout()
        maxstep_row.addWidget(QLabel("Max step size (h):"))
        self.max_step_spin = no_scroll(QDoubleSpinBox())
        self.max_step_spin.setRange(0.001, 10)
        self.max_step_spin.setValue(0.1)
        self.max_step_spin.setDecimals(4)
        maxstep_row.addWidget(self.max_step_spin)
        solver_layout.addLayout(maxstep_row)

        solver_group.setLayout(solver_layout)
        left_layout.addWidget(solver_group)

        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        left_layout.addWidget(self.progress_bar)

        # Buttons
        btn_layout = QHBoxLayout()

        self.run_btn = QPushButton("Run Simulation")
        self.run_btn.setStyleSheet("font-size: 14px; padding: 10px; background-color: #4CAF50; color: white;")
        self.run_btn.clicked.connect(self._run_simulation)
        btn_layout.addWidget(self.run_btn)

        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setEnabled(False)
        self.stop_btn.clicked.connect(self._stop_simulation)
        btn_layout.addWidget(self.stop_btn)

        left_layout.addLayout(btn_layout)

        left_layout.addStretch()

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.accept)
        left_layout.addWidget(close_btn)

        left_widget.setLayout(left_layout)

        # Right panel: Results
        right_widget = QWidget()
        right_layout = QVBoxLayout()

        # Results tabs
        self.results_tabs = QTabWidget()

        # Plot tab
        if MATPLOTLIB_AVAILABLE:
            self.plot_widget = QWidget()
            plot_layout = QVBoxLayout()

            self.figure = Figure(figsize=(8, 6))
            self.canvas = FigureCanvas(self.figure)
            plot_layout.addWidget(self.canvas)

            plot_btn_layout = QHBoxLayout()
            self.customize_btn = QPushButton("Customize Plot")
            self.customize_btn.setEnabled(False)
            self.customize_btn.clicked.connect(self._customize_plot)
            plot_btn_layout.addWidget(self.customize_btn)
            save_plot_btn = QPushButton("Save Plot...")
            save_plot_btn.clicked.connect(self._save_plot)
            plot_btn_layout.addWidget(save_plot_btn)
            plot_btn_layout.addStretch()
            plot_layout.addLayout(plot_btn_layout)

            self.plot_widget.setLayout(plot_layout)
            self.results_tabs.addTab(self.plot_widget, "📈 Plot")
        else:
            no_plot_label = QLabel("Matplotlib not available. Install with: pip install matplotlib")
            self.results_tabs.addTab(no_plot_label, "📈 Plot")

        # Data tab
        self.data_widget = QWidget()
        data_layout = QVBoxLayout()

        self.data_table = QTableWidget()
        data_layout.addWidget(self.data_table)

        export_btn = QPushButton("Export Data to CSV...")
        export_btn.clicked.connect(self._export_data)
        data_layout.addWidget(export_btn)

        self.data_widget.setLayout(data_layout)
        self.results_tabs.addTab(self.data_widget, "📊 Data")

        # Log tab
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.results_tabs.addTab(self.log_text, "📝 Log")

        right_layout.addWidget(self.results_tabs)
        right_widget.setLayout(right_layout)

        # Splitter
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(left_widget)
        splitter.addWidget(right_widget)
        splitter.setSizes([500, 700])

        main_layout.addWidget(splitter)
        self.setLayout(main_layout)

        # Check dependencies
        if not SCIPY_AVAILABLE:
            self.run_btn.setEnabled(False)
            self._log("Warning: SciPy is not installed. dFBA requires SciPy for ODE integration.")
            self._log("Install it with: pip install scipy")

    def _populate_reactions(self):
        """Populate reaction lists."""
        model = self.appdata.project.cobra_py_model

        # Find biomass reactions
        biomass_candidates = []
        for rxn in model.reactions:
            name_lower = rxn.name.lower() if rxn.name else ""
            id_lower = rxn.id.lower()
            if "biomass" in name_lower or "biomass" in id_lower or "growth" in name_lower:
                biomass_candidates.append(rxn.id)

        # Also add current objective
        obj_rxn = None
        for rxn in model.reactions:
            if rxn.objective_coefficient != 0:
                obj_rxn = rxn.id
                break

        if obj_rxn and obj_rxn not in biomass_candidates:
            biomass_candidates.insert(0, obj_rxn)

        for rid in biomass_candidates:
            self.biomass_combo.addItem(rid)

        # Find exchange reactions for substrates and products
        # Default substrates: only these 3
        default_substrates = ["EX_glc__D_e", "EX_o2_e", "EX_nh4_e"]
        # Common products
        common_products = ["EX_ac_e", "EX_etoh_e", "EX_lac__D_e", "EX_co2_e", "EX_for_e"]

        # Populate substrates table with only the 3 default substrates (if present in model)
        substrate_rxns = [rid for rid in default_substrates if rid in model.reactions]

        self.substrates_table.setRowCount(len(substrate_rxns))
        for row, rid in enumerate(substrate_rxns):
            # Checkbox
            cb = QCheckBox()
            cb.setChecked(True)  # All 3 defaults checked
            self.substrates_table.setCellWidget(row, 0, cb)

            # Reaction ID
            self.substrates_table.setItem(row, 1, QTableWidgetItem(rid))

            # Initial concentration (mM)
            init_spin = no_scroll(QDoubleSpinBox())
            init_spin.setRange(0, 10000)
            init_spin.setValue(10.0 if "glc" in rid else 1000.0)
            self.substrates_table.setCellWidget(row, 2, init_spin)

            # Km (mM)
            km_spin = no_scroll(QDoubleSpinBox())
            km_spin.setRange(0.001, 1000)
            km_spin.setValue(0.5 if "glc" in rid else 0.1)
            self.substrates_table.setCellWidget(row, 3, km_spin)

            # Vmax (mmol/gDW/h)
            vmax_spin = no_scroll(QDoubleSpinBox())
            vmax_spin.setRange(0.1, 100)
            vmax_spin.setValue(10.0 if "glc" in rid else 20.0)
            self.substrates_table.setCellWidget(row, 4, vmax_spin)

        # Populate products table
        product_rxns = []
        for rid in common_products:
            if rid in model.reactions:
                product_rxns.append(rid)

        self.products_table.setRowCount(len(product_rxns))
        for row, rid in enumerate(product_rxns):
            # Checkbox
            cb = QCheckBox()
            cb.setChecked(True)
            self.products_table.setCellWidget(row, 0, cb)

            # Reaction ID
            self.products_table.setItem(row, 1, QTableWidgetItem(rid))

            # Initial concentration (mM)
            init_spin = no_scroll(QDoubleSpinBox())
            init_spin.setRange(0, 10000)
            init_spin.setValue(0.0)
            self.products_table.setCellWidget(row, 2, init_spin)

    def _log(self, message: str):
        """Add message to log."""
        self.log_text.append(message)

    def _get_parameters(self) -> DFBAParameters:
        """Collect parameters from UI."""
        params = DFBAParameters()

        params.biomass_reaction = self.biomass_combo.currentText()
        params.initial_biomass = self.initial_biomass_spin.value()
        params.time_start = self.time_start_spin.value()
        params.time_end = self.time_end_spin.value()
        params.time_steps = self.time_steps_spin.value()
        params.method = self.method_combo.currentText()
        params.max_step = self.max_step_spin.value()

        # Collect substrates
        for row in range(self.substrates_table.rowCount()):
            cb = self.substrates_table.cellWidget(row, 0)
            if cb and cb.isChecked():
                rid = self.substrates_table.item(row, 1).text()
                init = self.substrates_table.cellWidget(row, 2).value()
                km = self.substrates_table.cellWidget(row, 3).value()
                vmax = self.substrates_table.cellWidget(row, 4).value()
                params.substrates[rid] = (init, km, vmax)

        # Collect products
        for row in range(self.products_table.rowCount()):
            cb = self.products_table.cellWidget(row, 0)
            if cb and cb.isChecked():
                rid = self.products_table.item(row, 1).text()
                init = self.products_table.cellWidget(row, 2).value()
                params.products[rid] = init

        return params

    @Slot()
    def _run_simulation(self):
        """Start dFBA simulation."""
        params = self._get_parameters()

        if not params.biomass_reaction:
            QMessageBox.warning(self, "No Biomass", "Please select a biomass reaction.")
            return

        if not params.substrates:
            QMessageBox.warning(self, "No Substrates", "Please select at least one substrate.")
            return

        self._log("Starting dFBA simulation...")
        self._log(f"  Biomass reaction: {params.biomass_reaction}")
        self._log(f"  Initial biomass: {params.initial_biomass} g/L")
        self._log(f"  Time: {params.time_start} - {params.time_end} h")
        self._log(f"  Substrates: {list(params.substrates.keys())}")
        self._log(f"  Products: {list(params.products.keys())}")

        self.run_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # Indeterminate

        # Run in thread
        model = self.appdata.project.cobra_py_model.copy()

        # Apply scenario
        self.appdata.project.load_scenario_into_model(model)

        self.simulation_thread = DFBASimulationThread(model, params)
        self.simulation_thread.finished.connect(self._on_simulation_finished)
        self.simulation_thread.start()

    @Slot()
    def _stop_simulation(self):
        """Stop running simulation."""
        if self.simulation_thread and self.simulation_thread.isRunning():
            self.simulation_thread.terminate()
            self._log("Simulation stopped by user.")
            self._reset_ui()

    def _reset_ui(self):
        """Reset UI after simulation."""
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        self.progress_bar.setVisible(False)

    @Slot(object)
    def _on_simulation_finished(self, result: DFBAResult):
        """Handle simulation completion."""
        self._reset_ui()
        self.last_result = result

        if result.success:
            self._log(f"✓ {result.message}")
            self._log(f"  Simulation time points: {len(result.time)}")
            self._log(f"  Final biomass: {result.biomass[-1]:.4f} g/L")

            for mid, conc in result.metabolites.items():
                self._log(f"  {mid}: {conc[0]:.2f} -> {conc[-1]:.2f} mM")

            self._update_plot(result)
            self._update_data_table(result)

            QMessageBox.information(
                self,
                "Simulation Complete",
                f"dFBA simulation completed.\n\nFinal biomass: {result.biomass[-1]:.4f} g/L",
            )
        else:
            self._log(f"✗ Simulation failed: {result.message}")
            QMessageBox.warning(self, "Simulation Failed", result.message)

    def _update_plot(self, result: DFBAResult):
        """Update the plot with simulation results."""
        if not MATPLOTLIB_AVAILABLE:
            return

        self.figure.clear()

        # Create subplots
        ax1 = self.figure.add_subplot(211)
        ax2 = self.figure.add_subplot(212)

        # Plot biomass
        ax1.plot(result.time, result.biomass, "b-", linewidth=2, label="Biomass")
        ax1.set_xlabel("Time (h)")
        ax1.set_ylabel("Biomass (g/L)", color="b")
        ax1.tick_params(axis="y", labelcolor="b")
        ax1.set_title("Dynamic FBA Simulation")
        ax1.legend(loc="upper left")
        ax1.grid(True, alpha=0.3)

        # Plot metabolites
        colors = plt.cm.tab10.colors
        for i, (mid, conc) in enumerate(result.metabolites.items()):
            ax2.plot(result.time, conc, color=colors[i % len(colors)], linewidth=1.5, label=mid)

        ax2.set_xlabel("Time (h)")
        ax2.set_ylabel("Concentration (mM)")
        ax2.legend(loc="upper right", fontsize=8)
        ax2.grid(True, alpha=0.3)

        self.figure.tight_layout()
        self.canvas.draw()
        self.customize_btn.setEnabled(True)

    def _update_data_table(self, result: DFBAResult):
        """Update the data table with simulation results."""
        # Columns: Time, Biomass, Metabolite1, Metabolite2, ...
        n_cols = 2 + len(result.metabolites)
        n_rows = len(result.time)

        self.data_table.setColumnCount(n_cols)
        self.data_table.setRowCount(n_rows)

        headers = ["Time (h)", "Biomass (g/L)"] + list(result.metabolites.keys())
        self.data_table.setHorizontalHeaderLabels(headers)

        for row in range(n_rows):
            self.data_table.setItem(row, 0, QTableWidgetItem(f"{result.time[row]:.4f}"))
            self.data_table.setItem(row, 1, QTableWidgetItem(f"{result.biomass[row]:.6f}"))

            for col, mid in enumerate(result.metabolites.keys()):
                self.data_table.setItem(row, 2 + col, QTableWidgetItem(f"{result.metabolites[mid][row]:.6f}"))

        self.data_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

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
            self, "Save Plot", self.appdata.work_directory, "PNG files (*.png);;PDF files (*.pdf);;SVG files (*.svg)"
        )

        if filename:
            self.figure.savefig(filename, dpi=150, bbox_inches="tight")
            self._log(f"Plot saved to {filename}")

    @Slot()
    def _export_data(self):
        """Export data to CSV."""
        if not self.last_result or not self.last_result.success:
            QMessageBox.warning(self, "No Data", "No simulation results to export.")
            return

        filename, _ = QFileDialog.getSaveFileName(self, "Export Data", self.appdata.work_directory, "CSV files (*.csv)")

        if filename:
            import csv

            result = self.last_result
            headers = ["Time (h)", "Biomass (g/L)"] + list(result.metabolites.keys())

            with open(filename, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(headers)

                for i in range(len(result.time)):
                    row = [result.time[i], result.biomass[i]]
                    for mid in result.metabolites:
                        row.append(result.metabolites[mid][i])
                    writer.writerow(row)

            self._log(f"Data exported to {filename}")
            QMessageBox.information(self, "Exported", f"Data exported to {filename}")

    def _pick_reaction(self, title: str) -> str | None:
        """Open a searchable dialog to pick a reaction from the model. Returns reaction ID or None."""
        model = self.appdata.project.cobra_py_model
        all_reaction_ids = [rxn.id for rxn in model.reactions]

        dialog = QDialog(self)
        dialog.setWindowTitle(title)
        dialog.setMinimumSize(400, 500)
        layout = QVBoxLayout()

        search_box = QLineEdit()
        search_box.setPlaceholderText("Search reactions...")
        layout.addWidget(search_box)

        list_widget = QListWidget()
        for rid in all_reaction_ids:
            rxn = model.reactions.get_by_id(rid)
            label = f"{rid}  —  {rxn.name}" if rxn.name else rid
            item = QListWidgetItem(label)
            item.setData(Qt.UserRole, rid)
            list_widget.addItem(item)
        layout.addWidget(list_widget)

        def _filter(text):
            text_lower = text.lower()
            for i in range(list_widget.count()):
                item = list_widget.item(i)
                item.setHidden(text_lower not in item.text().lower())

        search_box.textChanged.connect(_filter)

        btn_layout = QHBoxLayout()
        ok_btn = QPushButton("OK")
        ok_btn.clicked.connect(dialog.accept)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(dialog.reject)
        btn_layout.addWidget(ok_btn)
        btn_layout.addWidget(cancel_btn)
        layout.addLayout(btn_layout)

        list_widget.itemDoubleClicked.connect(lambda _: dialog.accept())

        dialog.setLayout(layout)

        if dialog.exec() == QDialog.Accepted:
            selected = list_widget.currentItem()
            if selected:
                return selected.data(Qt.UserRole)
        return None

    @Slot()
    def _add_substrate_row(self):
        """Add a reaction to the substrates table."""
        rid = self._pick_reaction("Add Substrate Reaction")
        if not rid:
            return

        # Check for duplicates
        for row in range(self.substrates_table.rowCount()):
            item = self.substrates_table.item(row, 1)
            if item and item.text() == rid:
                QMessageBox.information(self, "Already Added", f"{rid} is already in the substrates list.")
                return

        row = self.substrates_table.rowCount()
        self.substrates_table.setRowCount(row + 1)

        cb = QCheckBox()
        cb.setChecked(True)
        self.substrates_table.setCellWidget(row, 0, cb)
        self.substrates_table.setItem(row, 1, QTableWidgetItem(rid))

        init_spin = no_scroll(QDoubleSpinBox())
        init_spin.setRange(0, 10000)
        init_spin.setValue(10.0 if "glc" in rid else 1.0)
        self.substrates_table.setCellWidget(row, 2, init_spin)

        km_spin = no_scroll(QDoubleSpinBox())
        km_spin.setRange(0.001, 1000)
        km_spin.setValue(0.5 if "glc" in rid else 0.1)
        self.substrates_table.setCellWidget(row, 3, km_spin)

        vmax_spin = no_scroll(QDoubleSpinBox())
        vmax_spin.setRange(0.1, 100)
        vmax_spin.setValue(10.0 if "glc" in rid else 20.0)
        self.substrates_table.setCellWidget(row, 4, vmax_spin)

    @Slot()
    def _remove_substrate_row(self):
        """Remove the selected row from the substrates table."""
        selected_rows = sorted(set(idx.row() for idx in self.substrates_table.selectedIndexes()), reverse=True)
        for row in selected_rows:
            self.substrates_table.removeRow(row)

    @Slot()
    def _add_product_row(self):
        """Add a reaction to the products table."""
        rid = self._pick_reaction("Add Product Reaction")
        if not rid:
            return

        # Check for duplicates
        for row in range(self.products_table.rowCount()):
            item = self.products_table.item(row, 1)
            if item and item.text() == rid:
                QMessageBox.information(self, "Already Added", f"{rid} is already in the products list.")
                return

        row = self.products_table.rowCount()
        self.products_table.setRowCount(row + 1)

        cb = QCheckBox()
        cb.setChecked(True)
        self.products_table.setCellWidget(row, 0, cb)
        self.products_table.setItem(row, 1, QTableWidgetItem(rid))

        init_spin = no_scroll(QDoubleSpinBox())
        init_spin.setRange(0, 10000)
        init_spin.setValue(0.0)
        self.products_table.setCellWidget(row, 2, init_spin)

    @Slot()
    def _remove_product_row(self):
        """Remove the selected row from the products table."""
        selected_rows = sorted(set(idx.row() for idx in self.products_table.selectedIndexes()), reverse=True)
        for row in selected_rows:
            self.products_table.removeRow(row)

    @Slot(int)
    def _on_substrates_header_clicked(self, col: int):
        """Toggle all checkboxes in the substrates table when the Use header is clicked."""
        if col != 0:
            return
        all_checked = all(
            self.substrates_table.cellWidget(row, 0).isChecked()
            for row in range(self.substrates_table.rowCount())
            if self.substrates_table.cellWidget(row, 0)
        )
        new_state = not all_checked
        for row in range(self.substrates_table.rowCount()):
            cb = self.substrates_table.cellWidget(row, 0)
            if cb:
                cb.setChecked(new_state)

    @Slot(int)
    def _on_products_header_clicked(self, col: int):
        """Toggle all checkboxes in the products table when the Use header is clicked."""
        if col != 0:
            return
        all_checked = all(
            self.products_table.cellWidget(row, 0).isChecked()
            for row in range(self.products_table.rowCount())
            if self.products_table.cellWidget(row, 0)
        )
        new_state = not all_checked
        for row in range(self.products_table.rowCount()):
            cb = self.products_table.cellWidget(row, 0)
            if cb:
                cb.setChecked(new_state)
