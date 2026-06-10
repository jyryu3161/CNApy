"""Scenario Templates Dialog for CNApy

This module provides:
- Predefined scenario templates for common culture conditions
- Quick knockout generation
- Scenario bookmarks/favorites
"""

import json
import os
from dataclasses import asdict, dataclass

import appdirs
from qtpy.QtCore import Qt, Signal, Slot
from qtpy.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QInputDialog,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
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

# ============================================================================
# Scenario Template Definitions
# ============================================================================


@dataclass
class ScenarioTemplate:
    """A scenario template definition."""

    name: str
    description: str
    category: str
    reactions: dict[str, tuple[float, float]]  # reaction_id -> (lb, ub)
    # Common reaction patterns (will be matched against model)
    patterns: dict[str, tuple[float, float]]  # pattern -> (lb, ub)

    def to_dict(self) -> dict:
        return asdict(self)

    @staticmethod
    def from_dict(data: dict) -> "ScenarioTemplate":
        return ScenarioTemplate(**data)


# Default templates for common culture conditions
DEFAULT_TEMPLATES = [
    # ==================== Culture Conditions ====================
    ScenarioTemplate(
        name="Aerobic",
        description="Standard aerobic growth conditions with unlimited oxygen uptake.",
        category="Culture Conditions",
        reactions={},
        patterns={
            "EX_o2_e": (-1000, 1000),  # Oxygen exchange (unlimited)
            "EX_o2(e)": (-1000, 1000),
            "O2t": (-1000, 1000),
        },
    ),
    ScenarioTemplate(
        name="Anaerobic",
        description="Anaerobic conditions with no oxygen uptake.",
        category="Culture Conditions",
        reactions={},
        patterns={
            "EX_o2_e": (0, 0),  # No oxygen
            "EX_o2(e)": (0, 0),
            "O2t": (0, 0),
        },
    ),
    ScenarioTemplate(
        name="Microaerobic",
        description="Limited oxygen uptake (microaerobic conditions).",
        category="Culture Conditions",
        reactions={},
        patterns={
            "EX_o2_e": (-2, 0),  # Limited oxygen uptake
            "EX_o2(e)": (-2, 0),
            "O2t": (-2, 1000),
        },
    ),
    # ==================== Carbon Sources ====================
    ScenarioTemplate(
        name="Glucose",
        description="Glucose as primary carbon source (10 mmol/gDW/h uptake).",
        category="Carbon Sources",
        reactions={},
        patterns={
            "EX_glc__D_e": (-10, 0),
            "EX_glc_D_e": (-10, 0),
            "EX_glc(e)": (-10, 0),
            "GLCt1": (0, 10),
            "GLCpts": (0, 10),
        },
    ),
    ScenarioTemplate(
        name="Glucose (Limited)",
        description="Limited glucose uptake (5 mmol/gDW/h).",
        category="Carbon Sources",
        reactions={},
        patterns={
            "EX_glc__D_e": (-5, 0),
            "EX_glc_D_e": (-5, 0),
            "EX_glc(e)": (-5, 0),
        },
    ),
    ScenarioTemplate(
        name="Xylose",
        description="Xylose as primary carbon source.",
        category="Carbon Sources",
        reactions={},
        patterns={
            "EX_xyl__D_e": (-10, 0),
            "EX_xyl_D_e": (-10, 0),
            "EX_xyl(e)": (-10, 0),
            "EX_glc__D_e": (0, 0),  # No glucose
            "EX_glc_D_e": (0, 0),
            "EX_glc(e)": (0, 0),
        },
    ),
    ScenarioTemplate(
        name="Glycerol",
        description="Glycerol as primary carbon source.",
        category="Carbon Sources",
        reactions={},
        patterns={
            "EX_glyc_e": (-10, 0),
            "EX_glyc(e)": (-10, 0),
            "GLYCt": (0, 10),
            "EX_glc__D_e": (0, 0),  # No glucose
            "EX_glc_D_e": (0, 0),
        },
    ),
    ScenarioTemplate(
        name="Acetate",
        description="Acetate as primary carbon source.",
        category="Carbon Sources",
        reactions={},
        patterns={
            "EX_ac_e": (-10, 0),
            "EX_ac(e)": (-10, 0),
            "ACt2r": (-10, 10),
            "EX_glc__D_e": (0, 0),  # No glucose
            "EX_glc_D_e": (0, 0),
        },
    ),
    ScenarioTemplate(
        name="Succinate",
        description="Succinate as primary carbon source.",
        category="Carbon Sources",
        reactions={},
        patterns={
            "EX_succ_e": (-10, 0),
            "EX_succ(e)": (-10, 0),
            "EX_glc__D_e": (0, 0),
            "EX_glc_D_e": (0, 0),
        },
    ),
    # ==================== Production Scenarios ====================
    ScenarioTemplate(
        name="Ethanol Production",
        description="Optimize for ethanol production.",
        category="Production",
        reactions={},
        patterns={
            "EX_etoh_e": (0, 1000),  # Allow ethanol export
            "EX_etoh(e)": (0, 1000),
        },
    ),
    ScenarioTemplate(
        name="Lactate Production",
        description="Optimize for lactate production (fermentation).",
        category="Production",
        reactions={},
        patterns={
            "EX_lac__D_e": (0, 1000),
            "EX_lac_D_e": (0, 1000),
            "EX_lac(e)": (0, 1000),
        },
    ),
    ScenarioTemplate(
        name="Succinate Production",
        description="Optimize for succinate production.",
        category="Production",
        reactions={},
        patterns={
            "EX_succ_e": (0, 1000),
            "EX_succ(e)": (0, 1000),
        },
    ),
    # ==================== Stress Conditions ====================
    ScenarioTemplate(
        name="Nitrogen Limitation",
        description="Limited nitrogen source (ammonia).",
        category="Stress Conditions",
        reactions={},
        patterns={
            "EX_nh4_e": (-1, 0),  # Limited ammonia uptake
            "EX_nh4(e)": (-1, 0),
        },
    ),
    ScenarioTemplate(
        name="Phosphate Limitation",
        description="Limited phosphate availability.",
        category="Stress Conditions",
        reactions={},
        patterns={
            "EX_pi_e": (-0.5, 0),  # Limited phosphate
            "EX_pi(e)": (-0.5, 0),
        },
    ),
]


# ============================================================================
# Scenario Bookmarks Manager
# ============================================================================


class ScenarioBookmarksManager:
    """Manages scenario bookmarks/favorites."""

    def __init__(self):
        self.config_path = os.path.join(
            appdirs.user_config_dir("cnapy", roaming=True, appauthor=False), "scenario-bookmarks.json"
        )
        self.bookmarks: list[dict] = []
        self.load_bookmarks()

    def load_bookmarks(self):
        """Load bookmarks from file."""
        if os.path.exists(self.config_path):
            try:
                with open(self.config_path, encoding="utf-8") as f:
                    self.bookmarks = json.load(f)
            except Exception:
                self.bookmarks = []

    def save_bookmarks(self):
        """Save bookmarks to file."""
        try:
            os.makedirs(os.path.dirname(self.config_path), exist_ok=True)
            with open(self.config_path, "w", encoding="utf-8") as f:
                json.dump(self.bookmarks, f, indent=2, ensure_ascii=False)
        except Exception as e:
            print(f"Error saving bookmarks: {e}")

    def add_bookmark(self, name: str, description: str, scenario_data: dict[str, tuple[float, float]]) -> bool:
        """Add a new bookmark."""
        # Check for duplicate name
        for bm in self.bookmarks:
            if bm["name"] == name:
                return False

        self.bookmarks.append(
            {"name": name, "description": description, "scenario": {k: list(v) for k, v in scenario_data.items()}}
        )
        self.save_bookmarks()
        return True

    def remove_bookmark(self, name: str):
        """Remove a bookmark by name."""
        self.bookmarks = [bm for bm in self.bookmarks if bm["name"] != name]
        self.save_bookmarks()

    def get_bookmark(self, name: str) -> dict | None:
        """Get a bookmark by name."""
        for bm in self.bookmarks:
            if bm["name"] == name:
                return bm
        return None

    def update_bookmark(self, name: str, new_data: dict):
        """Update an existing bookmark."""
        for i, bm in enumerate(self.bookmarks):
            if bm["name"] == name:
                self.bookmarks[i] = new_data
                self.save_bookmarks()
                return True
        return False


# ============================================================================
# Scenario Templates Dialog
# ============================================================================


class ScenarioTemplatesDialog(QDialog):
    """Dialog for managing scenario templates and bookmarks."""

    scenarioApplied = Signal()

    def __init__(self, appdata: AppData):
        super().__init__()
        self.setWindowTitle("Scenario Templates & Bookmarks")
        self.appdata = appdata
        self.bookmarks_manager = ScenarioBookmarksManager()

        self.setMinimumSize(900, 700)
        self._setup_ui()
        self._load_templates()

    def _setup_ui(self):
        """Setup the dialog UI."""
        main_layout = QVBoxLayout()

        # Create tabs
        self.tabs = QTabWidget()

        # Tab 1: Templates
        self.templates_tab = self._create_templates_tab()
        self.tabs.addTab(self.templates_tab, "📋 Templates")

        # Tab 2: Quick Knockout
        self.knockout_tab = self._create_knockout_tab()
        self.tabs.addTab(self.knockout_tab, "🔧 Quick Knockout")

        # Tab 3: Bookmarks
        self.bookmarks_tab = self._create_bookmarks_tab()
        self.tabs.addTab(self.bookmarks_tab, "⭐ Bookmarks")

        # Tab 4: Custom Template
        self.custom_tab = self._create_custom_tab()
        self.tabs.addTab(self.custom_tab, "➕ Create Custom")

        main_layout.addWidget(self.tabs)

        # Bottom buttons
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()

        self.close_btn = QPushButton("Close")
        self.close_btn.clicked.connect(self.accept)
        btn_layout.addWidget(self.close_btn)

        main_layout.addLayout(btn_layout)
        self.setLayout(main_layout)

    def _create_templates_tab(self) -> QWidget:
        """Create the templates tab."""
        widget = QWidget()
        layout = QHBoxLayout()

        # Left: Category and template list
        left_widget = QWidget()
        left_layout = QVBoxLayout()

        # Category filter
        cat_layout = QHBoxLayout()
        cat_layout.addWidget(QLabel("Category:"))
        self.category_combo = no_scroll(QComboBox())
        self.category_combo.addItem("All Categories")
        self.category_combo.currentTextChanged.connect(self._filter_templates)
        cat_layout.addWidget(self.category_combo)
        left_layout.addLayout(cat_layout)

        # Template list
        self.template_list = QListWidget()
        self.template_list.currentItemChanged.connect(self._on_template_selected)
        left_layout.addWidget(self.template_list)

        left_widget.setLayout(left_layout)

        # Right: Template details and apply
        right_widget = QWidget()
        right_layout = QVBoxLayout()

        # Template info
        info_group = QGroupBox("Template Details")
        info_layout = QVBoxLayout()

        self.template_name_label = QLabel("")
        self.template_name_label.setStyleSheet("font-weight: bold; font-size: 14px;")
        info_layout.addWidget(self.template_name_label)

        self.template_desc_label = QLabel("")
        self.template_desc_label.setWordWrap(True)
        info_layout.addWidget(self.template_desc_label)

        # Reactions preview
        info_layout.addWidget(QLabel("Reactions to set:"))
        self.template_reactions_table = QTableWidget()
        self.template_reactions_table.setColumnCount(3)
        self.template_reactions_table.setHorizontalHeaderLabels(["Reaction", "Lower Bound", "Upper Bound"])
        self.template_reactions_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.template_reactions_table.setMaximumHeight(200)
        info_layout.addWidget(self.template_reactions_table)

        info_group.setLayout(info_layout)
        right_layout.addWidget(info_group)

        # Apply options
        options_group = QGroupBox("Apply Options")
        options_layout = QVBoxLayout()

        self.merge_checkbox = QCheckBox("Merge with current scenario (don't clear existing values)")
        options_layout.addWidget(self.merge_checkbox)

        self.run_fba_checkbox = QCheckBox("Run FBA after applying")
        self.run_fba_checkbox.setChecked(True)
        options_layout.addWidget(self.run_fba_checkbox)

        options_group.setLayout(options_layout)
        right_layout.addWidget(options_group)

        # Apply button
        self.apply_template_btn = QPushButton("Apply Template")
        self.apply_template_btn.setStyleSheet("font-size: 14px; padding: 10px;")
        self.apply_template_btn.clicked.connect(self._apply_template)
        right_layout.addWidget(self.apply_template_btn)

        right_layout.addStretch()
        right_widget.setLayout(right_layout)

        # Add to splitter
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(left_widget)
        splitter.addWidget(right_widget)
        splitter.setSizes([300, 500])

        layout.addWidget(splitter)
        widget.setLayout(layout)
        return widget

    def _create_knockout_tab(self) -> QWidget:
        """Create the quick knockout tab."""
        widget = QWidget()
        layout = QVBoxLayout()

        # Instructions
        desc = QLabel(
            "Quickly create knockout scenarios by selecting reactions to disable. "
            "Selected reactions will have their flux bounds set to (0, 0)."
        )
        desc.setWordWrap(True)
        layout.addWidget(desc)

        # Search
        search_layout = QHBoxLayout()
        search_layout.addWidget(QLabel("Search:"))
        self.knockout_search = QLineEdit()
        self.knockout_search.setPlaceholderText("Filter reactions by ID or name...")
        self.knockout_search.textChanged.connect(self._filter_knockout_list)
        search_layout.addWidget(self.knockout_search)
        layout.addLayout(search_layout)

        # Reactions list with checkboxes
        list_layout = QHBoxLayout()

        # Available reactions
        available_widget = QWidget()
        available_layout = QVBoxLayout()
        available_layout.addWidget(QLabel("Available Reactions:"))
        self.knockout_available_list = QListWidget()
        self.knockout_available_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        available_layout.addWidget(self.knockout_available_list)
        available_widget.setLayout(available_layout)

        # Buttons
        btn_widget = QWidget()
        btn_layout = QVBoxLayout()
        btn_layout.addStretch()

        add_btn = QPushButton(">>")
        add_btn.clicked.connect(self._add_to_knockout)
        btn_layout.addWidget(add_btn)

        remove_btn = QPushButton("<<")
        remove_btn.clicked.connect(self._remove_from_knockout)
        btn_layout.addWidget(remove_btn)

        btn_layout.addStretch()
        btn_widget.setLayout(btn_layout)

        # Selected for knockout
        knockout_widget = QWidget()
        knockout_layout = QVBoxLayout()
        knockout_layout.addWidget(QLabel("Selected for Knockout:"))
        self.knockout_selected_list = QListWidget()
        self.knockout_selected_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        knockout_layout.addWidget(self.knockout_selected_list)
        knockout_widget.setLayout(knockout_layout)

        list_layout.addWidget(available_widget)
        list_layout.addWidget(btn_widget)
        list_layout.addWidget(knockout_widget)
        layout.addLayout(list_layout)

        # Apply knockout
        apply_layout = QHBoxLayout()
        apply_layout.addStretch()

        self.knockout_merge_cb = QCheckBox("Merge with current scenario")
        apply_layout.addWidget(self.knockout_merge_cb)

        apply_knockout_btn = QPushButton("Apply Knockout")
        apply_knockout_btn.setStyleSheet("font-size: 14px; padding: 8px 20px;")
        apply_knockout_btn.clicked.connect(self._apply_knockout)
        apply_layout.addWidget(apply_knockout_btn)

        layout.addLayout(apply_layout)
        widget.setLayout(layout)
        return widget

    def _create_bookmarks_tab(self) -> QWidget:
        """Create the bookmarks tab."""
        widget = QWidget()
        layout = QVBoxLayout()

        # Instructions
        desc = QLabel(
            "Save frequently used scenarios as bookmarks for quick access. "
            "Bookmarks are stored locally and persist across sessions."
        )
        desc.setWordWrap(True)
        layout.addWidget(desc)

        # Bookmark list
        list_layout = QHBoxLayout()

        # Bookmarks list
        left_widget = QWidget()
        left_layout = QVBoxLayout()
        left_layout.addWidget(QLabel("Saved Bookmarks:"))
        self.bookmarks_list = QListWidget()
        self.bookmarks_list.currentItemChanged.connect(self._on_bookmark_selected)
        left_layout.addWidget(self.bookmarks_list)

        # Bookmark actions
        bm_btn_layout = QHBoxLayout()
        save_current_btn = QPushButton("Save Current Scenario")
        save_current_btn.clicked.connect(self._save_current_as_bookmark)
        bm_btn_layout.addWidget(save_current_btn)

        delete_bm_btn = QPushButton("Delete Selected")
        delete_bm_btn.clicked.connect(self._delete_bookmark)
        bm_btn_layout.addWidget(delete_bm_btn)

        left_layout.addLayout(bm_btn_layout)
        left_widget.setLayout(left_layout)

        # Bookmark details
        right_widget = QWidget()
        right_layout = QVBoxLayout()

        details_group = QGroupBox("Bookmark Details")
        details_layout = QVBoxLayout()

        self.bookmark_name_label = QLabel("")
        self.bookmark_name_label.setStyleSheet("font-weight: bold;")
        details_layout.addWidget(self.bookmark_name_label)

        self.bookmark_desc_label = QLabel("")
        self.bookmark_desc_label.setWordWrap(True)
        details_layout.addWidget(self.bookmark_desc_label)

        details_layout.addWidget(QLabel("Scenario Values:"))
        self.bookmark_values_table = QTableWidget()
        self.bookmark_values_table.setColumnCount(3)
        self.bookmark_values_table.setHorizontalHeaderLabels(["Reaction", "LB", "UB"])
        self.bookmark_values_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        details_layout.addWidget(self.bookmark_values_table)

        details_group.setLayout(details_layout)
        right_layout.addWidget(details_group)

        # Apply button
        apply_bm_btn = QPushButton("Apply Bookmark")
        apply_bm_btn.setStyleSheet("font-size: 14px; padding: 10px;")
        apply_bm_btn.clicked.connect(self._apply_bookmark)
        right_layout.addWidget(apply_bm_btn)

        right_widget.setLayout(right_layout)

        list_layout.addWidget(left_widget)
        list_layout.addWidget(right_widget)
        layout.addLayout(list_layout)

        widget.setLayout(layout)
        return widget

    def _create_custom_tab(self) -> QWidget:
        """Create custom template creation tab."""
        widget = QWidget()
        layout = QVBoxLayout()

        desc = QLabel(
            "Create a custom scenario template by specifying reaction bounds. "
            "Custom templates can be saved as bookmarks for reuse."
        )
        desc.setWordWrap(True)
        layout.addWidget(desc)

        # Template name
        name_layout = QHBoxLayout()
        name_layout.addWidget(QLabel("Template Name:"))
        self.custom_name_edit = QLineEdit()
        self.custom_name_edit.setPlaceholderText("My Custom Scenario")
        name_layout.addWidget(self.custom_name_edit)
        layout.addLayout(name_layout)

        # Description
        desc_layout = QHBoxLayout()
        desc_layout.addWidget(QLabel("Description:"))
        self.custom_desc_edit = QLineEdit()
        self.custom_desc_edit.setPlaceholderText("Description of the scenario...")
        desc_layout.addWidget(self.custom_desc_edit)
        layout.addLayout(desc_layout)

        # Reactions table
        layout.addWidget(QLabel("Reaction Bounds:"))

        self.custom_reactions_table = QTableWidget()
        self.custom_reactions_table.setColumnCount(4)
        self.custom_reactions_table.setHorizontalHeaderLabels(["Reaction ID", "Name", "Lower Bound", "Upper Bound"])
        self.custom_reactions_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.custom_reactions_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        layout.addWidget(self.custom_reactions_table)

        # Add reaction row
        add_row_layout = QHBoxLayout()
        self.custom_reaction_combo = no_scroll(QComboBox())
        self.custom_reaction_combo.setEditable(True)
        self.custom_reaction_combo.setPlaceholderText("Select or type reaction ID...")
        add_row_layout.addWidget(self.custom_reaction_combo)

        self.custom_lb_spin = no_scroll(QDoubleSpinBox())
        self.custom_lb_spin.setRange(-1000, 1000)
        self.custom_lb_spin.setValue(0)
        self.custom_lb_spin.setPrefix("LB: ")
        add_row_layout.addWidget(self.custom_lb_spin)

        self.custom_ub_spin = no_scroll(QDoubleSpinBox())
        self.custom_ub_spin.setRange(-1000, 1000)
        self.custom_ub_spin.setValue(0)
        self.custom_ub_spin.setPrefix("UB: ")
        add_row_layout.addWidget(self.custom_ub_spin)

        add_row_btn = QPushButton("Add")
        add_row_btn.clicked.connect(self._add_custom_row)
        add_row_layout.addWidget(add_row_btn)

        remove_row_btn = QPushButton("Remove Selected")
        remove_row_btn.clicked.connect(self._remove_custom_row)
        add_row_layout.addWidget(remove_row_btn)

        layout.addLayout(add_row_layout)

        # Save and Apply buttons
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()

        save_as_bookmark_btn = QPushButton("Save as Bookmark")
        save_as_bookmark_btn.clicked.connect(self._save_custom_as_bookmark)
        btn_layout.addWidget(save_as_bookmark_btn)

        apply_custom_btn = QPushButton("Apply")
        apply_custom_btn.setStyleSheet("font-size: 14px; padding: 8px 20px;")
        apply_custom_btn.clicked.connect(self._apply_custom)
        btn_layout.addWidget(apply_custom_btn)

        layout.addLayout(btn_layout)
        widget.setLayout(layout)
        return widget

    def _load_templates(self):
        """Load templates and populate lists."""
        # Get categories
        categories = set()
        for template in DEFAULT_TEMPLATES:
            categories.add(template.category)

        for cat in sorted(categories):
            self.category_combo.addItem(cat)

        # Populate template list
        self._filter_templates("All Categories")

        # Load reactions for knockout tab
        self._load_knockout_reactions()

        # Load bookmarks
        self._refresh_bookmarks_list()

        # Load reactions for custom tab
        self._load_custom_reactions()

    def _filter_templates(self, category: str):
        """Filter templates by category."""
        self.template_list.clear()

        for template in DEFAULT_TEMPLATES:
            if category == "All Categories" or template.category == category:
                item = QListWidgetItem(f"{template.name}")
                item.setData(Qt.UserRole, template)
                item.setToolTip(template.description)
                self.template_list.addItem(item)

    def _on_template_selected(self, current: QListWidgetItem, previous: QListWidgetItem):
        """Handle template selection."""
        if current is None:
            return

        template: ScenarioTemplate = current.data(Qt.UserRole)
        self.template_name_label.setText(template.name)
        self.template_desc_label.setText(template.description)

        # Show matching reactions
        self.template_reactions_table.setRowCount(0)
        model = self.appdata.project.cobra_py_model

        matched_reactions = []
        for pattern, bounds in template.patterns.items():
            if pattern in model.reactions:
                matched_reactions.append((pattern, bounds))

        for reac_id, bounds in template.reactions.items():
            if reac_id in model.reactions:
                matched_reactions.append((reac_id, bounds))

        self.template_reactions_table.setRowCount(len(matched_reactions))
        for row, (reac_id, bounds) in enumerate(matched_reactions):
            self.template_reactions_table.setItem(row, 0, QTableWidgetItem(reac_id))
            self.template_reactions_table.setItem(row, 1, QTableWidgetItem(str(bounds[0])))
            self.template_reactions_table.setItem(row, 2, QTableWidgetItem(str(bounds[1])))

        if len(matched_reactions) == 0:
            self.template_desc_label.setText(
                template.description + "\n\n⚠️ No matching reactions found in current model."
            )

    @Slot()
    def _apply_template(self):
        """Apply selected template."""
        current = self.template_list.currentItem()
        if current is None:
            QMessageBox.warning(self, "No Selection", "Please select a template to apply.")
            return

        template: ScenarioTemplate = current.data(Qt.UserRole)
        model = self.appdata.project.cobra_py_model

        if not self.merge_checkbox.isChecked():
            self.appdata.scen_values_clear()

        applied_count = 0

        # Apply pattern-matched reactions
        for pattern, bounds in template.patterns.items():
            if pattern in model.reactions:
                self.appdata.scen_values_set(pattern, (bounds[0], bounds[1]))
                applied_count += 1

        # Apply explicit reactions
        for reac_id, bounds in template.reactions.items():
            if reac_id in model.reactions:
                self.appdata.scen_values_set(reac_id, (bounds[0], bounds[1]))
                applied_count += 1

        if applied_count == 0:
            QMessageBox.warning(self, "No Matches", "No reactions from the template were found in the current model.")
            return

        QMessageBox.information(
            self, "Template Applied", f"Applied {applied_count} reaction bound(s) from template '{template.name}'."
        )

        self.scenarioApplied.emit()

    def _load_knockout_reactions(self):
        """Load reactions for knockout tab."""
        self.knockout_available_list.clear()
        model = self.appdata.project.cobra_py_model

        for rxn in model.reactions:
            item = QListWidgetItem(f"{rxn.id}: {rxn.name}")
            item.setData(Qt.UserRole, rxn.id)
            self.knockout_available_list.addItem(item)

    def _filter_knockout_list(self, text: str):
        """Filter knockout reaction list."""
        text = text.lower()
        for i in range(self.knockout_available_list.count()):
            item = self.knockout_available_list.item(i)
            visible = text in item.text().lower()
            item.setHidden(not visible)

    @Slot()
    def _add_to_knockout(self):
        """Add selected reactions to knockout list."""
        selected = self.knockout_available_list.selectedItems()
        existing = set()
        for i in range(self.knockout_selected_list.count()):
            existing.add(self.knockout_selected_list.item(i).data(Qt.UserRole))

        for item in selected:
            reac_id = item.data(Qt.UserRole)
            if reac_id not in existing:
                new_item = QListWidgetItem(item.text())
                new_item.setData(Qt.UserRole, reac_id)
                self.knockout_selected_list.addItem(new_item)

    @Slot()
    def _remove_from_knockout(self):
        """Remove selected from knockout list."""
        for item in self.knockout_selected_list.selectedItems():
            self.knockout_selected_list.takeItem(self.knockout_selected_list.row(item))

    @Slot()
    def _apply_knockout(self):
        """Apply knockout scenario."""
        if self.knockout_selected_list.count() == 0:
            QMessageBox.warning(self, "No Selection", "Please select reactions to knockout.")
            return

        if not self.knockout_merge_cb.isChecked():
            self.appdata.scen_values_clear()

        reactions = []
        values = []
        for i in range(self.knockout_selected_list.count()):
            reac_id = self.knockout_selected_list.item(i).data(Qt.UserRole)
            reactions.append(reac_id)
            values.append((0.0, 0.0))

        self.appdata.scen_values_set_multiple(reactions, values)

        QMessageBox.information(self, "Knockout Applied", f"Applied knockout to {len(reactions)} reaction(s).")

        self.scenarioApplied.emit()

    def _refresh_bookmarks_list(self):
        """Refresh the bookmarks list."""
        self.bookmarks_list.clear()
        self.bookmarks_manager.load_bookmarks()

        for bm in self.bookmarks_manager.bookmarks:
            item = QListWidgetItem(f"⭐ {bm['name']}")
            item.setData(Qt.UserRole, bm)
            item.setToolTip(bm.get("description", ""))
            self.bookmarks_list.addItem(item)

    def _on_bookmark_selected(self, current: QListWidgetItem, previous: QListWidgetItem):
        """Handle bookmark selection."""
        if current is None:
            return

        bm = current.data(Qt.UserRole)
        self.bookmark_name_label.setText(bm["name"])
        self.bookmark_desc_label.setText(bm.get("description", ""))

        scenario = bm.get("scenario", {})
        self.bookmark_values_table.setRowCount(len(scenario))

        for row, (reac_id, bounds) in enumerate(scenario.items()):
            self.bookmark_values_table.setItem(row, 0, QTableWidgetItem(reac_id))
            self.bookmark_values_table.setItem(row, 1, QTableWidgetItem(str(bounds[0])))
            self.bookmark_values_table.setItem(row, 2, QTableWidgetItem(str(bounds[1])))

    @Slot()
    def _save_current_as_bookmark(self):
        """Save current scenario as bookmark."""
        name, ok = QInputDialog.getText(self, "Save Bookmark", "Enter bookmark name:", QLineEdit.Normal, "")

        if not ok or not name.strip():
            return

        desc, ok = QInputDialog.getText(self, "Save Bookmark", "Enter description (optional):", QLineEdit.Normal, "")

        scenario_data = dict(self.appdata.project.scen_values)

        if self.bookmarks_manager.add_bookmark(name.strip(), desc, scenario_data):
            QMessageBox.information(self, "Saved", f"Bookmark '{name}' saved successfully.")
            self._refresh_bookmarks_list()
        else:
            QMessageBox.warning(self, "Error", f"Bookmark with name '{name}' already exists.")

    @Slot()
    def _delete_bookmark(self):
        """Delete selected bookmark."""
        current = self.bookmarks_list.currentItem()
        if current is None:
            return

        bm = current.data(Qt.UserRole)
        reply = QMessageBox.question(
            self,
            "Delete Bookmark",
            f"Are you sure you want to delete bookmark '{bm['name']}'?",
            QMessageBox.Yes | QMessageBox.No,
        )

        if reply == QMessageBox.Yes:
            self.bookmarks_manager.remove_bookmark(bm["name"])
            self._refresh_bookmarks_list()

    @Slot()
    def _apply_bookmark(self):
        """Apply selected bookmark."""
        current = self.bookmarks_list.currentItem()
        if current is None:
            QMessageBox.warning(self, "No Selection", "Please select a bookmark to apply.")
            return

        bm = current.data(Qt.UserRole)
        scenario = bm.get("scenario", {})

        self.appdata.scen_values_clear()

        reactions = []
        values = []
        model = self.appdata.project.cobra_py_model

        for reac_id, bounds in scenario.items():
            if reac_id in model.reactions:
                reactions.append(reac_id)
                values.append((bounds[0], bounds[1]))

        if reactions:
            self.appdata.scen_values_set_multiple(reactions, values)

        QMessageBox.information(
            self, "Bookmark Applied", f"Applied bookmark '{bm['name']}' with {len(reactions)} reaction(s)."
        )

        self.scenarioApplied.emit()

    def _load_custom_reactions(self):
        """Load reactions for custom template."""
        self.custom_reaction_combo.clear()
        model = self.appdata.project.cobra_py_model

        for rxn in model.reactions:
            self.custom_reaction_combo.addItem(f"{rxn.id}", rxn.id)

    @Slot()
    def _add_custom_row(self):
        """Add row to custom template table."""
        reac_id = self.custom_reaction_combo.currentData()
        if not reac_id:
            reac_id = self.custom_reaction_combo.currentText()

        if not reac_id:
            return

        model = self.appdata.project.cobra_py_model
        rxn_name = ""
        if reac_id in model.reactions:
            rxn_name = model.reactions.get_by_id(reac_id).name

        row = self.custom_reactions_table.rowCount()
        self.custom_reactions_table.setRowCount(row + 1)

        self.custom_reactions_table.setItem(row, 0, QTableWidgetItem(reac_id))
        self.custom_reactions_table.setItem(row, 1, QTableWidgetItem(rxn_name))
        self.custom_reactions_table.setItem(row, 2, QTableWidgetItem(str(self.custom_lb_spin.value())))
        self.custom_reactions_table.setItem(row, 3, QTableWidgetItem(str(self.custom_ub_spin.value())))

    @Slot()
    def _remove_custom_row(self):
        """Remove selected row from custom template table."""
        current_row = self.custom_reactions_table.currentRow()
        if current_row >= 0:
            self.custom_reactions_table.removeRow(current_row)

    @Slot()
    def _save_custom_as_bookmark(self):
        """Save custom template as bookmark."""
        name = self.custom_name_edit.text().strip()
        if not name:
            QMessageBox.warning(self, "Error", "Please enter a template name.")
            return

        scenario_data = {}
        for row in range(self.custom_reactions_table.rowCount()):
            reac_id = self.custom_reactions_table.item(row, 0).text()
            try:
                lb = float(self.custom_reactions_table.item(row, 2).text())
                ub = float(self.custom_reactions_table.item(row, 3).text())
                scenario_data[reac_id] = (lb, ub)
            except ValueError:
                continue

        desc = self.custom_desc_edit.text()

        if self.bookmarks_manager.add_bookmark(name, desc, scenario_data):
            QMessageBox.information(self, "Saved", f"Custom template '{name}' saved as bookmark.")
            self._refresh_bookmarks_list()
            self.tabs.setCurrentWidget(self.bookmarks_tab)
        else:
            QMessageBox.warning(self, "Error", f"Bookmark with name '{name}' already exists.")

    @Slot()
    def _apply_custom(self):
        """Apply custom template."""
        if self.custom_reactions_table.rowCount() == 0:
            QMessageBox.warning(self, "No Reactions", "Please add reactions to the custom template.")
            return

        self.appdata.scen_values_clear()

        reactions = []
        values = []
        model = self.appdata.project.cobra_py_model

        for row in range(self.custom_reactions_table.rowCount()):
            reac_id = self.custom_reactions_table.item(row, 0).text()
            try:
                lb = float(self.custom_reactions_table.item(row, 2).text())
                ub = float(self.custom_reactions_table.item(row, 3).text())
                if reac_id in model.reactions:
                    reactions.append(reac_id)
                    values.append((lb, ub))
            except ValueError:
                continue

        if reactions:
            self.appdata.scen_values_set_multiple(reactions, values)

        QMessageBox.information(self, "Applied", f"Applied custom template with {len(reactions)} reaction(s).")

        self.scenarioApplied.emit()
