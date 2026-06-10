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

"""Media Management Dialog for CNApy

This module provides functionality for:
- Predefined media templates for bacteria, plant cells, and animal cells
- Save/load custom media compositions
- Apply media to current model (set exchange reaction bounds)
- Display current media information
"""

import json
import os
from dataclasses import asdict, dataclass

import appdirs
from qtpy.QtCore import Qt, Signal
from qtpy.QtGui import QColor
from qtpy.QtWidgets import (
    QCheckBox,
    QDialog,
    QFileDialog,
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

# ============================================================================
# Media Template Definitions
# ============================================================================


@dataclass
class MediaComponent:
    """A single component in a media definition."""

    metabolite_id: str  # Exchange reaction ID pattern (e.g., EX_glc__D_e)
    name: str
    default_uptake: float  # Default uptake rate (negative for uptake)
    is_essential: bool = False  # Whether this is an essential component


@dataclass
class MediaTemplate:
    """A complete media template definition."""

    name: str
    description: str
    category: str  # "Bacteria", "Plant", "Animal", "Custom"
    organism_type: str  # More specific organism type
    components: dict[str, float]  # exchange_reaction_id -> uptake_rate (lb)
    patterns: dict[str, float]  # pattern -> uptake_rate for flexible matching

    def to_dict(self) -> dict:
        return asdict(self)

    @staticmethod
    def from_dict(data: dict) -> "MediaTemplate":
        return MediaTemplate(**data)


# ============================================================================
# Default Media Templates
# ============================================================================

# Common exchange reaction ID patterns for different naming conventions
GLUCOSE_PATTERNS = ["EX_glc__D_e", "EX_glc_D_e", "EX_glc(e)", "EX_glc_e", "GLCt1", "GLCpts"]
OXYGEN_PATTERNS = ["EX_o2_e", "EX_o2(e)", "O2t"]
AMMONIA_PATTERNS = ["EX_nh4_e", "EX_nh4(e)", "NH4t"]
PHOSPHATE_PATTERNS = ["EX_pi_e", "EX_pi(e)", "PIt"]
SULFATE_PATTERNS = ["EX_so4_e", "EX_so4(e)", "SO4t"]
CO2_PATTERNS = ["EX_co2_e", "EX_co2(e)", "CO2t"]
H2O_PATTERNS = ["EX_h2o_e", "EX_h2o(e)", "H2Ot"]
PROTON_PATTERNS = ["EX_h_e", "EX_h(e)", "Ht"]

# Bacteria Media Templates
BACTERIA_TEMPLATES = [
    MediaTemplate(
        name="M9 Minimal (Glucose)",
        description="M9 minimal medium with glucose as carbon source. Standard for E. coli and other bacteria.",
        category="Bacteria",
        organism_type="E. coli / General Bacteria",
        components={},
        patterns={
            "EX_glc__D_e": -10.0,  # Glucose uptake
            "EX_glc_D_e": -10.0,
            "EX_nh4_e": -1000.0,  # Ammonia (unlimited)
            "EX_nh4(e)": -1000.0,
            "EX_pi_e": -1000.0,  # Phosphate
            "EX_pi(e)": -1000.0,
            "EX_so4_e": -1000.0,  # Sulfate
            "EX_so4(e)": -1000.0,
            "EX_o2_e": -1000.0,  # Oxygen (aerobic)
            "EX_o2(e)": -1000.0,
            "EX_h2o_e": -1000.0,  # Water
            "EX_h_e": -1000.0,  # Protons
            "EX_fe2_e": -1000.0,  # Iron
            "EX_mg2_e": -1000.0,  # Magnesium
            "EX_ca2_e": -1000.0,  # Calcium
            "EX_k_e": -1000.0,  # Potassium
            "EX_cl_e": -1000.0,  # Chloride
            "EX_na1_e": -1000.0,  # Sodium
            "EX_co2_e": -1000.0,  # CO2
        },
    ),
    MediaTemplate(
        name="M9 Minimal (Glycerol)",
        description="M9 minimal medium with glycerol as carbon source.",
        category="Bacteria",
        organism_type="E. coli / General Bacteria",
        components={},
        patterns={
            "EX_glyc_e": -10.0,  # Glycerol uptake
            "EX_glyc(e)": -10.0,
            "EX_glc__D_e": 0.0,  # No glucose
            "EX_glc_D_e": 0.0,
            "EX_nh4_e": -1000.0,
            "EX_nh4(e)": -1000.0,
            "EX_pi_e": -1000.0,
            "EX_o2_e": -1000.0,
            "EX_h2o_e": -1000.0,
            "EX_fe2_e": -1000.0,
            "EX_mg2_e": -1000.0,
            "EX_so4_e": -1000.0,
        },
    ),
    MediaTemplate(
        name="LB (Lysogeny Broth)",
        description="Rich medium for rapid bacterial growth. Contains amino acids and vitamins.",
        category="Bacteria",
        organism_type="E. coli / General Bacteria",
        components={},
        patterns={
            "EX_glc__D_e": -10.0,  # Glucose
            "EX_nh4_e": -1000.0,  # Ammonia
            "EX_o2_e": -1000.0,
            "EX_pi_e": -1000.0,
            # Amino acids (typical in LB)
            "EX_ala__L_e": -1.0,
            "EX_arg__L_e": -1.0,
            "EX_asn__L_e": -1.0,
            "EX_asp__L_e": -1.0,
            "EX_cys__L_e": -1.0,
            "EX_gln__L_e": -1.0,
            "EX_glu__L_e": -1.0,
            "EX_gly_e": -1.0,
            "EX_his__L_e": -1.0,
            "EX_ile__L_e": -1.0,
            "EX_leu__L_e": -1.0,
            "EX_lys__L_e": -1.0,
            "EX_met__L_e": -1.0,
            "EX_phe__L_e": -1.0,
            "EX_pro__L_e": -1.0,
            "EX_ser__L_e": -1.0,
            "EX_thr__L_e": -1.0,
            "EX_trp__L_e": -1.0,
            "EX_tyr__L_e": -1.0,
            "EX_val__L_e": -1.0,
        },
    ),
    MediaTemplate(
        name="Anaerobic M9 (Glucose)",
        description="M9 minimal medium without oxygen for anaerobic growth.",
        category="Bacteria",
        organism_type="E. coli / General Bacteria",
        components={},
        patterns={
            "EX_glc__D_e": -10.0,
            "EX_glc_D_e": -10.0,
            "EX_nh4_e": -1000.0,
            "EX_pi_e": -1000.0,
            "EX_so4_e": -1000.0,
            "EX_o2_e": 0.0,  # No oxygen!
            "EX_o2(e)": 0.0,
            "EX_h2o_e": -1000.0,
            "EX_co2_e": -1000.0,
        },
    ),
]

# Plant Cell Media Templates
PLANT_TEMPLATES = [
    MediaTemplate(
        name="MS Medium (Murashige & Skoog)",
        description="Standard plant cell culture medium with sucrose. Used for most plant tissue cultures.",
        category="Plant",
        organism_type="Plant Cells / Arabidopsis",
        components={},
        patterns={
            # Sucrose as carbon source
            "EX_sucr_e": -10.0,
            "EX_sucr(e)": -10.0,
            # Or glucose/fructose
            "EX_glc__D_e": -5.0,
            "EX_fru_e": -5.0,
            # Nitrogen sources
            "EX_no3_e": -20.0,  # Nitrate (primary N source for plants)
            "EX_nh4_e": -5.0,  # Ammonium
            # Phosphate
            "EX_pi_e": -1000.0,
            # Sulfate
            "EX_so4_e": -1000.0,
            # CO2 (for photosynthesis if applicable)
            "EX_co2_e": -1000.0,
            # Light energy (photons) if modeled
            "EX_photon_e": -100.0,
            # Water and ions
            "EX_h2o_e": -1000.0,
            "EX_k_e": -1000.0,
            "EX_ca2_e": -1000.0,
            "EX_mg2_e": -1000.0,
            "EX_fe2_e": -1000.0,
            "EX_cl_e": -1000.0,
        },
    ),
    MediaTemplate(
        name="Heterotrophic Plant Culture",
        description="Plant cell culture in dark (heterotrophic) conditions with glucose.",
        category="Plant",
        organism_type="Plant Cells",
        components={},
        patterns={
            "EX_glc__D_e": -10.0,
            "EX_sucr_e": -5.0,
            "EX_no3_e": -20.0,
            "EX_nh4_e": -5.0,
            "EX_pi_e": -1000.0,
            "EX_so4_e": -1000.0,
            "EX_o2_e": -1000.0,  # Oxygen for respiration
            "EX_co2_e": 1000.0,  # CO2 release (no photosynthesis)
            "EX_photon_e": 0.0,  # No light
        },
    ),
    MediaTemplate(
        name="Autotrophic (Light + CO2)",
        description="Photosynthetic conditions with CO2 and light as energy source.",
        category="Plant",
        organism_type="Plant Cells / Algae",
        components={},
        patterns={
            "EX_co2_e": -10.0,  # CO2 uptake
            "EX_photon_e": -100.0,  # Light energy
            "EX_no3_e": -10.0,
            "EX_nh4_e": -5.0,
            "EX_pi_e": -1000.0,
            "EX_so4_e": -1000.0,
            "EX_h2o_e": -1000.0,
            "EX_o2_e": 1000.0,  # O2 production
            "EX_glc__D_e": 0.0,  # No glucose needed
        },
    ),
]

# Animal Cell Media Templates
ANIMAL_TEMPLATES = [
    MediaTemplate(
        name="DMEM (Dulbecco's Modified Eagle Medium)",
        description="Standard medium for mammalian cell culture. Contains glucose, amino acids, and vitamins.",
        category="Animal",
        organism_type="Mammalian Cells / CHO / HEK293",
        components={},
        patterns={
            # Glucose (high glucose DMEM)
            "EX_glc__D_e": -25.0,
            "EX_glc_D_e": -25.0,
            # Glutamine (essential for mammalian cells)
            "EX_gln__L_e": -4.0,
            "EX_gln_L_e": -4.0,
            # Essential amino acids
            "EX_arg__L_e": -1.0,
            "EX_cys__L_e": -0.5,
            "EX_his__L_e": -0.5,
            "EX_ile__L_e": -1.0,
            "EX_leu__L_e": -1.0,
            "EX_lys__L_e": -1.0,
            "EX_met__L_e": -0.5,
            "EX_phe__L_e": -0.5,
            "EX_ser__L_e": -1.0,
            "EX_thr__L_e": -1.0,
            "EX_trp__L_e": -0.2,
            "EX_tyr__L_e": -0.5,
            "EX_val__L_e": -1.0,
            # Vitamins
            "EX_pnto__R_e": -0.01,  # Pantothenate
            "EX_chol_e": -0.01,  # Choline
            "EX_ncam_e": -0.01,  # Nicotinamide
            "EX_pydx_e": -0.01,  # Pyridoxine
            "EX_ribflv_e": -0.01,  # Riboflavin
            "EX_thm_e": -0.01,  # Thiamine
            # Oxygen
            "EX_o2_e": -1000.0,
            # Inorganic
            "EX_pi_e": -1000.0,
            "EX_so4_e": -1000.0,
            "EX_fe2_e": -1000.0,
            "EX_ca2_e": -1000.0,
            "EX_mg2_e": -1000.0,
            "EX_k_e": -1000.0,
            "EX_na1_e": -1000.0,
            "EX_cl_e": -1000.0,
        },
    ),
    MediaTemplate(
        name="RPMI 1640",
        description="Medium for suspension cells (lymphocytes, hybridomas). Lower glucose than DMEM.",
        category="Animal",
        organism_type="Mammalian Cells / Immune Cells",
        components={},
        patterns={
            "EX_glc__D_e": -11.0,  # Lower glucose
            "EX_gln__L_e": -2.0,
            "EX_arg__L_e": -1.0,
            "EX_asn__L_e": -0.5,
            "EX_asp__L_e": -0.5,
            "EX_cys__L_e": -0.5,
            "EX_glu__L_e": -0.5,
            "EX_gly_e": -0.5,
            "EX_his__L_e": -0.5,
            "EX_ile__L_e": -0.5,
            "EX_leu__L_e": -0.5,
            "EX_lys__L_e": -0.5,
            "EX_met__L_e": -0.3,
            "EX_phe__L_e": -0.3,
            "EX_pro__L_e": -0.5,
            "EX_ser__L_e": -0.5,
            "EX_thr__L_e": -0.5,
            "EX_trp__L_e": -0.1,
            "EX_tyr__L_e": -0.3,
            "EX_val__L_e": -0.5,
            "EX_o2_e": -1000.0,
            "EX_pi_e": -1000.0,
        },
    ),
    MediaTemplate(
        name="Serum-Free Medium",
        description="Chemically defined medium without serum. Used for biopharmaceutical production.",
        category="Animal",
        organism_type="CHO / Industrial Cells",
        components={},
        patterns={
            "EX_glc__D_e": -30.0,
            "EX_gln__L_e": -6.0,
            # All amino acids
            "EX_ala__L_e": -1.0,
            "EX_arg__L_e": -1.5,
            "EX_asn__L_e": -2.0,
            "EX_asp__L_e": -1.0,
            "EX_cys__L_e": -1.0,
            "EX_glu__L_e": -1.0,
            "EX_gly_e": -1.0,
            "EX_his__L_e": -0.5,
            "EX_ile__L_e": -1.0,
            "EX_leu__L_e": -1.5,
            "EX_lys__L_e": -1.5,
            "EX_met__L_e": -0.5,
            "EX_phe__L_e": -0.5,
            "EX_pro__L_e": -1.0,
            "EX_ser__L_e": -2.0,
            "EX_thr__L_e": -1.0,
            "EX_trp__L_e": -0.3,
            "EX_tyr__L_e": -0.5,
            "EX_val__L_e": -1.0,
            "EX_o2_e": -1000.0,
            "EX_pi_e": -1000.0,
            "EX_so4_e": -1000.0,
        },
    ),
]

# Combine all default templates
DEFAULT_MEDIA_TEMPLATES = BACTERIA_TEMPLATES + PLANT_TEMPLATES + ANIMAL_TEMPLATES


# ============================================================================
# Custom Media Manager
# ============================================================================


class CustomMediaManager:
    """Manages custom (user-defined) media."""

    def __init__(self):
        self.config_path = os.path.join(
            appdirs.user_config_dir("cnapy", roaming=True, appauthor=False), "custom-media.json"
        )
        self.custom_media: list[dict] = []
        self.load_media()

    def load_media(self):
        """Load custom media from file."""
        if os.path.exists(self.config_path):
            try:
                with open(self.config_path, encoding="utf-8") as f:
                    self.custom_media = json.load(f)
            except Exception:
                self.custom_media = []

    def save_media(self):
        """Save custom media to file."""
        try:
            os.makedirs(os.path.dirname(self.config_path), exist_ok=True)
            with open(self.config_path, "w", encoding="utf-8") as f:
                json.dump(self.custom_media, f, indent=2, ensure_ascii=False)
        except Exception as e:
            print(f"Error saving custom media: {e}")

    def add_media(self, name: str, description: str, category: str, components: dict[str, float]) -> bool:
        """Add a new custom media."""
        for m in self.custom_media:
            if m["name"] == name:
                return False

        self.custom_media.append(
            {
                "name": name,
                "description": description,
                "category": category,
                "organism_type": "Custom",
                "components": components,
                "patterns": {},
            }
        )
        self.save_media()
        return True

    def remove_media(self, name: str):
        """Remove a custom media by name."""
        self.custom_media = [m for m in self.custom_media if m["name"] != name]
        self.save_media()

    def get_media(self, name: str) -> dict | None:
        """Get a custom media by name."""
        for m in self.custom_media:
            if m["name"] == name:
                return m
        return None

    def update_media(self, name: str, new_data: dict) -> bool:
        """Update an existing custom media."""
        for i, m in enumerate(self.custom_media):
            if m["name"] == name:
                self.custom_media[i] = new_data
                self.save_media()
                return True
        return False


# ============================================================================
# Media Management Dialog
# ============================================================================


class MediaManagementDialog(QDialog):
    """Dialog for managing and applying media to metabolic models."""

    mediaApplied = Signal()

    def __init__(self, appdata: AppData, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Media Management")
        self.setMinimumSize(1000, 700)
        self.appdata = appdata
        self.custom_media_manager = CustomMediaManager()

        self._setup_ui()
        self._load_templates()
        self._update_current_media_display()

    def _setup_ui(self):
        """Setup the dialog UI."""
        main_layout = QVBoxLayout()

        # Current media display at top
        current_group = QGroupBox("Current Media Status")
        current_layout = QVBoxLayout()

        self.current_media_label = QLabel("No media information available.")
        self.current_media_label.setWordWrap(True)
        current_layout.addWidget(self.current_media_label)

        # Exchange reactions table (compact)
        self.current_exchanges_table = QTableWidget()
        self.current_exchanges_table.setColumnCount(4)
        self.current_exchanges_table.setHorizontalHeaderLabels(["Reaction", "Metabolite", "Lower Bound", "Upper Bound"])
        self.current_exchanges_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.current_exchanges_table.setMaximumHeight(150)
        current_layout.addWidget(self.current_exchanges_table)

        refresh_btn = QPushButton("Refresh Current Status")
        refresh_btn.clicked.connect(self._update_current_media_display)
        current_layout.addWidget(refresh_btn)

        current_group.setLayout(current_layout)
        main_layout.addWidget(current_group)

        # Tabs for template selection
        self.tabs = QTabWidget()

        # Tab 1: Bacteria
        self.bacteria_tab = self._create_category_tab("Bacteria", BACTERIA_TEMPLATES)
        self.tabs.addTab(self.bacteria_tab, "🦠 Bacteria")

        # Tab 2: Plant
        self.plant_tab = self._create_category_tab("Plant", PLANT_TEMPLATES)
        self.tabs.addTab(self.plant_tab, "🌱 Plant")

        # Tab 3: Animal
        self.animal_tab = self._create_category_tab("Animal", ANIMAL_TEMPLATES)
        self.tabs.addTab(self.animal_tab, "🐁 Animal")

        # Tab 4: Custom
        self.custom_tab = self._create_custom_tab()
        self.tabs.addTab(self.custom_tab, "⚙️ Custom")

        main_layout.addWidget(self.tabs)

        # Bottom buttons
        btn_layout = QHBoxLayout()

        export_btn = QPushButton("Export Current as JSON...")
        export_btn.clicked.connect(self._export_current_media)
        btn_layout.addWidget(export_btn)

        import_btn = QPushButton("Import Media from JSON...")
        import_btn.clicked.connect(self._import_media_from_json)
        btn_layout.addWidget(import_btn)

        btn_layout.addStretch()

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.accept)
        btn_layout.addWidget(close_btn)

        main_layout.addLayout(btn_layout)
        self.setLayout(main_layout)

    def _create_category_tab(self, category: str, templates: list[MediaTemplate]) -> QWidget:
        """Create a tab for a specific category of media."""
        widget = QWidget()
        layout = QHBoxLayout()

        # Left: Template list
        left_widget = QWidget()
        left_layout = QVBoxLayout()
        left_layout.addWidget(QLabel(f"{category} Media Templates:"))

        template_list = QListWidget()
        template_list.setObjectName(f"{category.lower()}_list")

        for template in templates:
            item = QListWidgetItem(f"📋 {template.name}")
            item.setData(Qt.UserRole, template)
            item.setToolTip(template.description)
            template_list.addItem(item)

        template_list.currentItemChanged.connect(lambda current, prev: self._on_template_selected(current, category))
        left_layout.addWidget(template_list)
        left_widget.setLayout(left_layout)

        # Right: Details and apply
        right_widget = QWidget()
        right_layout = QVBoxLayout()

        # Details group
        details_group = QGroupBox("Template Details")
        details_layout = QVBoxLayout()

        name_label = QLabel("")
        name_label.setObjectName(f"{category.lower()}_name")
        name_label.setStyleSheet("font-weight: bold; font-size: 14px;")
        details_layout.addWidget(name_label)

        desc_label = QLabel("")
        desc_label.setObjectName(f"{category.lower()}_desc")
        desc_label.setWordWrap(True)
        details_layout.addWidget(desc_label)

        org_label = QLabel("")
        org_label.setObjectName(f"{category.lower()}_org")
        org_label.setStyleSheet("color: gray;")
        details_layout.addWidget(org_label)

        # Components table
        details_layout.addWidget(QLabel("Media Components:"))
        comp_table = QTableWidget()
        comp_table.setObjectName(f"{category.lower()}_components")
        comp_table.setColumnCount(3)
        comp_table.setHorizontalHeaderLabels(["Exchange Reaction", "Bound (LB)", "Status"])
        comp_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        comp_table.setMaximumHeight(200)
        details_layout.addWidget(comp_table)

        details_group.setLayout(details_layout)
        right_layout.addWidget(details_group)

        # Apply options
        options_group = QGroupBox("Apply Options")
        options_layout = QVBoxLayout()

        clear_cb = QCheckBox("Clear scenario before applying")
        clear_cb.setObjectName(f"{category.lower()}_clear")
        clear_cb.setChecked(True)
        options_layout.addWidget(clear_cb)

        only_exchange_cb = QCheckBox("Only set exchange reactions (leave internal reactions unchanged)")
        only_exchange_cb.setObjectName(f"{category.lower()}_only_exchange")
        only_exchange_cb.setChecked(True)
        options_layout.addWidget(only_exchange_cb)

        options_group.setLayout(options_layout)
        right_layout.addWidget(options_group)

        # Apply button
        apply_btn = QPushButton("Apply Media")
        apply_btn.setObjectName(f"{category.lower()}_apply")
        apply_btn.setStyleSheet("font-size: 14px; padding: 10px; background-color: #4CAF50; color: white;")
        apply_btn.clicked.connect(lambda: self._apply_template(category))
        right_layout.addWidget(apply_btn)

        right_layout.addStretch()
        right_widget.setLayout(right_layout)

        # Splitter
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(left_widget)
        splitter.addWidget(right_widget)
        splitter.setSizes([300, 500])

        layout.addWidget(splitter)
        widget.setLayout(layout)
        return widget

    def _create_custom_tab(self) -> QWidget:
        """Create the custom media tab."""
        widget = QWidget()
        layout = QVBoxLayout()

        desc = QLabel("Create custom media compositions or save current exchange reaction settings.")
        desc.setWordWrap(True)
        layout.addWidget(desc)

        # Two columns
        content_layout = QHBoxLayout()

        # Left: Saved custom media
        left_widget = QWidget()
        left_layout = QVBoxLayout()
        left_layout.addWidget(QLabel("Saved Custom Media:"))

        self.custom_media_list = QListWidget()
        self.custom_media_list.currentItemChanged.connect(self._on_custom_media_selected)
        left_layout.addWidget(self.custom_media_list)

        btn_layout = QHBoxLayout()
        save_current_btn = QPushButton("Save Current as Custom")
        save_current_btn.clicked.connect(self._save_current_as_custom)
        btn_layout.addWidget(save_current_btn)

        delete_btn = QPushButton("Delete Selected")
        delete_btn.clicked.connect(self._delete_custom_media)
        btn_layout.addWidget(delete_btn)

        left_layout.addLayout(btn_layout)
        left_widget.setLayout(left_layout)

        # Right: Custom media details
        right_widget = QWidget()
        right_layout = QVBoxLayout()

        details_group = QGroupBox("Custom Media Details")
        details_layout = QVBoxLayout()

        self.custom_name_label = QLabel("")
        self.custom_name_label.setStyleSheet("font-weight: bold;")
        details_layout.addWidget(self.custom_name_label)

        self.custom_desc_label = QLabel("")
        self.custom_desc_label.setWordWrap(True)
        details_layout.addWidget(self.custom_desc_label)

        details_layout.addWidget(QLabel("Components:"))
        self.custom_comp_table = QTableWidget()
        self.custom_comp_table.setColumnCount(2)
        self.custom_comp_table.setHorizontalHeaderLabels(["Exchange Reaction", "Bound (LB)"])
        self.custom_comp_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        details_layout.addWidget(self.custom_comp_table)

        details_group.setLayout(details_layout)
        right_layout.addWidget(details_group)

        apply_custom_btn = QPushButton("Apply Selected Custom Media")
        apply_custom_btn.setStyleSheet("font-size: 14px; padding: 10px;")
        apply_custom_btn.clicked.connect(self._apply_custom_media)
        right_layout.addWidget(apply_custom_btn)

        right_layout.addStretch()
        right_widget.setLayout(right_layout)

        content_layout.addWidget(left_widget)
        content_layout.addWidget(right_widget)

        layout.addLayout(content_layout)
        widget.setLayout(layout)
        return widget

    def _load_templates(self):
        """Load templates and custom media."""
        self._refresh_custom_media_list()

    def _refresh_custom_media_list(self):
        """Refresh the custom media list."""
        self.custom_media_list.clear()
        self.custom_media_manager.load_media()

        for media in self.custom_media_manager.custom_media:
            item = QListWidgetItem(f"⭐ {media['name']}")
            item.setData(Qt.UserRole, media)
            item.setToolTip(media.get("description", ""))
            self.custom_media_list.addItem(item)

    def _update_current_media_display(self):
        """Update the display of current media/exchange reactions."""
        model = self.appdata.project.cobra_py_model

        # Find all exchange reactions
        exchange_rxns = [r for r in model.reactions if r.id.startswith("EX_")]

        # Filter to those with non-default bounds or in scenario
        active_exchanges = []
        for rxn in exchange_rxns:
            # Check if in scenario
            if rxn.id in self.appdata.project.scen_values:
                lb, ub = self.appdata.project.scen_values[rxn.id]
                active_exchanges.append((rxn.id, rxn.name, lb, ub, "Scenario"))
            elif rxn.lower_bound != 0 or rxn.upper_bound != 0:
                active_exchanges.append((rxn.id, rxn.name, rxn.lower_bound, rxn.upper_bound, "Model"))

        # Update table
        self.current_exchanges_table.setRowCount(len(active_exchanges[:20]))  # Show top 20
        for row, (rid, name, lb, ub, source) in enumerate(active_exchanges[:20]):
            self.current_exchanges_table.setItem(row, 0, QTableWidgetItem(rid))
            self.current_exchanges_table.setItem(row, 1, QTableWidgetItem(name[:30]))
            self.current_exchanges_table.setItem(row, 2, QTableWidgetItem(f"{lb:.2f}"))
            self.current_exchanges_table.setItem(row, 3, QTableWidgetItem(f"{ub:.2f}"))

        # Update label
        uptake_count = sum(1 for _, _, lb, _, _ in active_exchanges if lb < 0)
        secrete_count = sum(1 for _, _, _, ub, _ in active_exchanges if ub > 0)

        self.current_media_label.setText(
            f"Exchange reactions: {len(exchange_rxns)} total, "
            f"{uptake_count} with uptake, {secrete_count} with secretion. "
            f"(Showing top 20 active exchanges)"
        )

    def _on_template_selected(self, current: QListWidgetItem, category: str):
        """Handle template selection."""
        if current is None:
            return

        template: MediaTemplate = current.data(Qt.UserRole)
        cat_lower = category.lower()

        # Update labels
        name_label = self.findChild(QLabel, f"{cat_lower}_name")
        if name_label:
            name_label.setText(template.name)

        desc_label = self.findChild(QLabel, f"{cat_lower}_desc")
        if desc_label:
            desc_label.setText(template.description)

        org_label = self.findChild(QLabel, f"{cat_lower}_org")
        if org_label:
            org_label.setText(f"Organism: {template.organism_type}")

        # Update components table
        comp_table = self.findChild(QTableWidget, f"{cat_lower}_components")
        if comp_table:
            model = self.appdata.project.cobra_py_model
            all_components = {**template.components, **template.patterns}

            # Find matching reactions
            matched = []
            for pattern, bound in all_components.items():
                if pattern in model.reactions:
                    matched.append((pattern, bound, "✓ Found"))
                else:
                    matched.append((pattern, bound, "✗ Not in model"))

            comp_table.setRowCount(len(matched))
            for row, (rxn_id, bound, status) in enumerate(matched):
                comp_table.setItem(row, 0, QTableWidgetItem(rxn_id))
                comp_table.setItem(row, 1, QTableWidgetItem(f"{bound:.2f}"))
                status_item = QTableWidgetItem(status)
                if "Found" in status:
                    status_item.setForeground(QColor(0, 150, 0))
                else:
                    status_item.setForeground(QColor(150, 0, 0))
                comp_table.setItem(row, 2, status_item)

    def _apply_template(self, category: str):
        """Apply selected template."""
        cat_lower = category.lower()
        template_list = self.findChild(QListWidget, f"{cat_lower}_list")

        if not template_list or not template_list.currentItem():
            QMessageBox.warning(self, "No Selection", "Please select a media template.")
            return

        template: MediaTemplate = template_list.currentItem().data(Qt.UserRole)
        model = self.appdata.project.cobra_py_model

        # Check clear option
        clear_cb = self.findChild(QCheckBox, f"{cat_lower}_clear")
        if clear_cb and clear_cb.isChecked():
            self.appdata.project.scen_values.clear_flux_values()

        # Apply components
        applied_count = 0
        all_components = {**template.components, **template.patterns}

        reactions = []
        values = []

        for pattern, value in all_components.items():
            if pattern in model.reactions:
                rxn = model.reactions.get_by_id(pattern)
                if value < 0:
                    # Negative value = uptake: use it as the lower bound, keep upper bound
                    new_lb = value
                    new_ub = rxn.upper_bound
                else:
                    # Positive value = allowed secretion/production: use it as the upper
                    # bound and keep the existing (negative) lower bound, or 0 otherwise
                    new_lb = rxn.lower_bound if rxn.lower_bound < 0 else 0
                    new_ub = value
                reactions.append(pattern)
                values.append((new_lb, new_ub))
                applied_count += 1

        if reactions:
            self.appdata.scen_values_set_multiple(reactions, values)

        QMessageBox.information(
            self, "Media Applied", f"Applied media '{template.name}' with {applied_count} exchange reaction(s)."
        )

        self._update_current_media_display()
        self.mediaApplied.emit()

    def _on_custom_media_selected(self, current: QListWidgetItem, previous: QListWidgetItem):
        """Handle custom media selection."""
        if current is None:
            return

        media = current.data(Qt.UserRole)
        self.custom_name_label.setText(media["name"])
        self.custom_desc_label.setText(media.get("description", ""))

        components = {**media.get("components", {}), **media.get("patterns", {})}
        self.custom_comp_table.setRowCount(len(components))

        for row, (rxn_id, bound) in enumerate(components.items()):
            self.custom_comp_table.setItem(row, 0, QTableWidgetItem(rxn_id))
            self.custom_comp_table.setItem(row, 1, QTableWidgetItem(f"{bound:.2f}"))

    def _save_current_as_custom(self):
        """Save current exchange settings as custom media."""
        name, ok = QInputDialog.getText(self, "Save Custom Media", "Enter media name:", QLineEdit.Normal, "")

        if not ok or not name.strip():
            return

        desc, ok = QInputDialog.getText(
            self, "Save Custom Media", "Enter description (optional):", QLineEdit.Normal, ""
        )

        # Collect current exchange reaction bounds from scenario
        components = {}
        for rxn_id, (lb, ub) in self.appdata.project.scen_values.items():
            if rxn_id.startswith("EX_"):
                components[rxn_id] = lb

        if not components:
            # Also include exchange reactions with non-zero bounds from model
            for rxn in self.appdata.project.cobra_py_model.reactions:
                if rxn.id.startswith("EX_") and rxn.lower_bound != 0:
                    components[rxn.id] = rxn.lower_bound

        if not components:
            QMessageBox.warning(self, "No Exchange Reactions", "No exchange reactions found in scenario or model.")
            return

        if self.custom_media_manager.add_media(name.strip(), desc, "Custom", components):
            QMessageBox.information(self, "Saved", f"Custom media '{name}' saved successfully.")
            self._refresh_custom_media_list()
        else:
            QMessageBox.warning(self, "Error", f"Media with name '{name}' already exists.")

    def _delete_custom_media(self):
        """Delete selected custom media."""
        current = self.custom_media_list.currentItem()
        if not current:
            return

        media = current.data(Qt.UserRole)
        reply = QMessageBox.question(
            self,
            "Delete Media",
            f"Are you sure you want to delete '{media['name']}'?",
            QMessageBox.Yes | QMessageBox.No,
        )

        if reply == QMessageBox.Yes:
            self.custom_media_manager.remove_media(media["name"])
            self._refresh_custom_media_list()

    def _apply_custom_media(self):
        """Apply selected custom media."""
        current = self.custom_media_list.currentItem()
        if not current:
            QMessageBox.warning(self, "No Selection", "Please select a custom media.")
            return

        media = current.data(Qt.UserRole)
        model = self.appdata.project.cobra_py_model

        # Clear scenario
        self.appdata.project.scen_values.clear_flux_values()

        # Apply components
        components = {**media.get("components", {}), **media.get("patterns", {})}
        applied_count = 0

        reactions = []
        values = []

        for rxn_id, lb in components.items():
            if rxn_id in model.reactions:
                rxn = model.reactions.get_by_id(rxn_id)
                ub = rxn.upper_bound
                reactions.append(rxn_id)
                values.append((lb, ub))
                applied_count += 1

        if reactions:
            self.appdata.scen_values_set_multiple(reactions, values)

        QMessageBox.information(
            self, "Media Applied", f"Applied custom media '{media['name']}' with {applied_count} reaction(s)."
        )

        self._update_current_media_display()
        self.mediaApplied.emit()

    def _export_current_media(self):
        """Export current media settings to JSON."""
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export Media", self.appdata.work_directory, "JSON files (*.json)"
        )

        if not filename:
            return

        # Collect current settings
        components = {}
        for rxn_id, (lb, ub) in self.appdata.project.scen_values.items():
            if rxn_id.startswith("EX_"):
                components[rxn_id] = lb

        media_data = {
            "name": "Exported Media",
            "description": "Exported from CNApy",
            "category": "Custom",
            "organism_type": "Custom",
            "components": components,
            "patterns": {},
        }

        try:
            with open(filename, "w", encoding="utf-8") as f:
                json.dump(media_data, f, indent=2, ensure_ascii=False)
            QMessageBox.information(self, "Exported", f"Media exported to {filename}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to export: {str(e)}")

    def _import_media_from_json(self):
        """Import media from JSON file."""
        filename, _ = QFileDialog.getOpenFileName(
            self, "Import Media", self.appdata.work_directory, "JSON files (*.json)"
        )

        if not filename:
            return

        try:
            with open(filename, encoding="utf-8") as f:
                media_data = json.load(f)

            name = media_data.get("name", "Imported Media")
            desc = media_data.get("description", "")
            components = media_data.get("components", {})
            patterns = media_data.get("patterns", {})

            all_components = {**components, **patterns}

            if self.custom_media_manager.add_media(name, desc, "Custom", all_components):
                QMessageBox.information(self, "Imported", f"Media '{name}' imported successfully.")
                self._refresh_custom_media_list()
            else:
                QMessageBox.warning(self, "Error", f"Media with name '{name}' already exists.")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import: {str(e)}")
