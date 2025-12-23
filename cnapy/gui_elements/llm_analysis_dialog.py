"""LLM-based Strain Analysis Dialog for CNApy

This module provides functionality to:
- Analyze reactions and genes existence in specific strains using LLM APIs
- Support OpenAI (ChatGPT) and Google Gemini Flash APIs
- Manage API keys per project or globally
- Utilize web search capabilities for real-time information
"""

import os
import json
import re
from typing import Dict, List, Optional, Tuple, Any
from configparser import ConfigParser

from qtpy.QtCore import Qt, Slot, Signal, QThread
from qtpy.QtWidgets import (QDialog, QHBoxLayout, QLabel, QVBoxLayout,
                            QPushButton, QGroupBox, QTextEdit, QTabWidget,
                            QWidget, QMessageBox, QCheckBox, QListWidget,
                            QListWidgetItem, QAbstractItemView, QFileDialog,
                            QComboBox, QLineEdit, QGridLayout, QProgressBar,
                            QRadioButton, QButtonGroup, QSpinBox, QSplitter,
                            QTableWidget, QTableWidgetItem, QHeaderView,
                            QApplication)
from qtpy.QtGui import QFont, QColor

import appdirs

from cnapy.appdata import AppData


class LLMConfig:
    """Manages LLM API configurations."""

    def __init__(self):
        self.config_path = os.path.join(
            appdirs.user_config_dir("cnapy", roaming=True, appauthor=False),
            "llm-config.json"
        )
        self.openai_api_key = ""
        self.gemini_api_key = ""
        self.default_provider = "openai"  # "openai" or "gemini"
        self.default_model_openai = "gpt-4o"
        self.default_model_gemini = "gemini-2.0-flash-exp"
        self.load_config()

    def load_config(self):
        """Load configuration from file."""
        if os.path.exists(self.config_path):
            try:
                with open(self.config_path, 'r') as f:
                    data = json.load(f)
                    self.openai_api_key = data.get('openai_api_key', '')
                    self.gemini_api_key = data.get('gemini_api_key', '')
                    self.default_provider = data.get('default_provider', 'openai')
                    self.default_model_openai = data.get('default_model_openai', 'gpt-4o')
                    self.default_model_gemini = data.get('default_model_gemini', 'gemini-2.0-flash-exp')
            except Exception:
                pass

    def save_config(self):
        """Save configuration to file."""
        try:
            os.makedirs(os.path.dirname(self.config_path), exist_ok=True)
            with open(self.config_path, 'w') as f:
                json.dump({
                    'openai_api_key': self.openai_api_key,
                    'gemini_api_key': self.gemini_api_key,
                    'default_provider': self.default_provider,
                    'default_model_openai': self.default_model_openai,
                    'default_model_gemini': self.default_model_gemini
                }, f, indent=2)
        except Exception as e:
            print(f"Error saving LLM config: {e}")


class LLMAnalysisThread(QThread):
    """Background thread for LLM API calls."""

    result_ready = Signal(str, str, dict)  # item_id, item_type, result
    error_occurred = Signal(str, str)  # item_id, error message
    progress_update = Signal(int, int)  # current, total

    def __init__(self, config: LLMConfig, items: List[Tuple[str, str, str]],
                 strain_name: str, provider: str, model: str,
                 use_web_search: bool = True):
        """
        Args:
            config: LLM configuration
            items: List of (item_id, item_name, item_type) tuples
            strain_name: Name of the strain to analyze
            provider: "openai" or "gemini"
            model: Model name to use
            use_web_search: Whether to use web search (for OpenAI)
        """
        super().__init__()
        self.config = config
        self.items = items
        self.strain_name = strain_name
        self.provider = provider
        self.model = model
        self.use_web_search = use_web_search
        self._stop_requested = False

    def stop(self):
        self._stop_requested = True

    def run(self):
        """Execute LLM analysis for all items."""
        total = len(self.items)

        for i, (item_id, item_name, item_type) in enumerate(self.items):
            if self._stop_requested:
                break

            self.progress_update.emit(i + 1, total)

            try:
                if self.provider == "openai":
                    result = self._analyze_with_openai(item_id, item_name, item_type)
                else:
                    result = self._analyze_with_gemini(item_id, item_name, item_type)

                self.result_ready.emit(item_id, item_type, result)
            except Exception as e:
                self.error_occurred.emit(item_id, str(e))

    def _build_prompt(self, item_id: str, item_name: str, item_type: str) -> str:
        """Build the analysis prompt."""
        if item_type == "reaction":
            return f"""Analyze whether the following metabolic reaction exists or is likely to exist in {self.strain_name}.

Reaction ID: {item_id}
Reaction Name: {item_name}

Please provide:
1. Existence status: Does this reaction exist in {self.strain_name}? (Yes/No/Unknown/Likely)
2. Confidence level: How confident is this assessment? (High/Medium/Low)
3. Evidence: What evidence supports this conclusion?
4. Alternative: If the reaction doesn't exist, is there an alternative pathway?
5. References: Any relevant database entries or literature references.

Respond in JSON format:
{{
    "exists": "Yes/No/Unknown/Likely",
    "confidence": "High/Medium/Low",
    "evidence": "explanation",
    "alternative": "alternative pathway or null",
    "references": ["ref1", "ref2"]
}}"""
        else:  # gene
            return f"""Analyze whether the following gene exists or has an ortholog in {self.strain_name}.

Gene ID: {item_id}
Gene Name: {item_name}

Please provide:
1. Existence status: Does this gene or its ortholog exist in {self.strain_name}? (Yes/No/Unknown/Likely)
2. Confidence level: How confident is this assessment? (High/Medium/Low)
3. Evidence: What evidence supports this conclusion?
4. Ortholog: If there's an ortholog, what is its ID/name in {self.strain_name}?
5. Function: Brief description of the gene's function.
6. References: Any relevant database entries or literature references.

Respond in JSON format:
{{
    "exists": "Yes/No/Unknown/Likely",
    "confidence": "High/Medium/Low",
    "evidence": "explanation",
    "ortholog": "ortholog gene ID or null",
    "function": "gene function description",
    "references": ["ref1", "ref2"]
}}"""

    def _analyze_with_openai(self, item_id: str, item_name: str, item_type: str) -> dict:
        """Analyze using OpenAI API."""
        try:
            import openai
        except ImportError:
            raise ImportError("OpenAI package not installed. Install with: pip install openai")

        client = openai.OpenAI(api_key=self.config.openai_api_key)

        prompt = self._build_prompt(item_id, item_name, item_type)

        messages = [
            {"role": "system", "content": "You are a bioinformatics expert specializing in metabolic network analysis and comparative genomics. Provide accurate, evidence-based analysis."},
            {"role": "user", "content": prompt}
        ]

        # Use web search if available and requested
        kwargs = {
            "model": self.model,
            "messages": messages,
            "temperature": 0.3,
        }

        # For models that support web search
        if self.use_web_search and "gpt-4" in self.model:
            try:
                # Try with web_search_options for newer API versions
                kwargs["web_search_options"] = {"search_context_size": "medium"}
            except:
                pass

        try:
            response = client.chat.completions.create(**kwargs)
        except Exception as e:
            # Retry without web search if it fails
            if "web_search_options" in kwargs:
                del kwargs["web_search_options"]
                response = client.chat.completions.create(**kwargs)
            else:
                raise e

        content = response.choices[0].message.content

        # Parse JSON from response
        return self._parse_json_response(content)

    def _analyze_with_gemini(self, item_id: str, item_name: str, item_type: str) -> dict:
        """Analyze using Google Gemini API."""
        try:
            import google.generativeai as genai
        except ImportError:
            raise ImportError("Google Generative AI package not installed. Install with: pip install google-generativeai")

        genai.configure(api_key=self.config.gemini_api_key)

        model = genai.GenerativeModel(
            model_name=self.model,
            generation_config={
                "temperature": 0.3,
                "max_output_tokens": 2048,
            }
        )

        prompt = self._build_prompt(item_id, item_name, item_type)
        system_prompt = "You are a bioinformatics expert specializing in metabolic network analysis and comparative genomics. Provide accurate, evidence-based analysis."

        full_prompt = f"{system_prompt}\n\n{prompt}"

        # Gemini 2.0 Flash supports grounding with Google Search
        try:
            from google.generativeai.types import Tool
            tools = [Tool(google_search={})]
            response = model.generate_content(full_prompt, tools=tools)
        except Exception:
            # Fall back to regular generation
            response = model.generate_content(full_prompt)

        content = response.text

        return self._parse_json_response(content)

    def _parse_json_response(self, content: str) -> dict:
        """Parse JSON from LLM response."""
        # Try to extract JSON from the response
        try:
            # Find JSON block in markdown code block
            json_match = re.search(r'```(?:json)?\s*(\{.*?\})\s*```', content, re.DOTALL)
            if json_match:
                return json.loads(json_match.group(1))

            # Try to find raw JSON
            json_match = re.search(r'\{[^{}]*\}', content, re.DOTALL)
            if json_match:
                return json.loads(json_match.group(0))

            # Return raw content if JSON parsing fails
            return {
                "exists": "Unknown",
                "confidence": "Low",
                "evidence": content,
                "raw_response": True
            }
        except json.JSONDecodeError:
            return {
                "exists": "Unknown",
                "confidence": "Low",
                "evidence": content,
                "raw_response": True
            }


class LLMAnalysisDialog(QDialog):
    """Dialog for LLM-based strain analysis."""

    def __init__(self, appdata: AppData):
        super().__init__()
        self.setWindowTitle("LLM-based Strain Analysis")
        self.appdata = appdata
        self.config = LLMConfig()
        self.analysis_thread: Optional[LLMAnalysisThread] = None
        self.results: Dict[str, dict] = {}  # item_id -> result

        self.setMinimumSize(1000, 800)
        self._setup_ui()
        self._load_model_items()

    def _setup_ui(self):
        """Setup the dialog UI."""
        main_layout = QVBoxLayout()

        # Create tabs
        self.tabs = QTabWidget()

        # Tab 1: API Configuration
        self.config_tab = self._create_config_tab()
        self.tabs.addTab(self.config_tab, "API Configuration")

        # Tab 2: Analysis
        self.analysis_tab = self._create_analysis_tab()
        self.tabs.addTab(self.analysis_tab, "Strain Analysis")

        # Tab 3: Results
        self.results_tab = self._create_results_tab()
        self.tabs.addTab(self.results_tab, "Results")

        main_layout.addWidget(self.tabs)

        # Close button
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()
        self.close_btn = QPushButton("Close")
        self.close_btn.clicked.connect(self.accept)
        btn_layout.addWidget(self.close_btn)
        main_layout.addLayout(btn_layout)

        self.setLayout(main_layout)

    def _create_config_tab(self) -> QWidget:
        """Create API configuration tab."""
        widget = QWidget()
        layout = QVBoxLayout()

        # Instructions
        desc = QLabel(
            "Configure your LLM API keys. Keys are stored securely in your user config directory. "
            "You need at least one API key to use the analysis feature."
        )
        desc.setWordWrap(True)
        layout.addWidget(desc)

        # OpenAI Configuration
        openai_group = QGroupBox("OpenAI (ChatGPT) Configuration")
        openai_layout = QGridLayout()

        openai_layout.addWidget(QLabel("API Key:"), 0, 0)
        self.openai_key_edit = QLineEdit()
        self.openai_key_edit.setPlaceholderText("sk-...")
        self.openai_key_edit.setEchoMode(QLineEdit.Password)
        self.openai_key_edit.setText(self.config.openai_api_key)
        openai_layout.addWidget(self.openai_key_edit, 0, 1)

        self.openai_show_key = QCheckBox("Show")
        self.openai_show_key.toggled.connect(
            lambda checked: self.openai_key_edit.setEchoMode(
                QLineEdit.Normal if checked else QLineEdit.Password
            )
        )
        openai_layout.addWidget(self.openai_show_key, 0, 2)

        openai_layout.addWidget(QLabel("Model:"), 1, 0)
        self.openai_model_combo = QComboBox()
        self.openai_model_combo.addItems([
            "gpt-4o",
            "gpt-4o-mini",
            "gpt-4-turbo",
            "gpt-4",
            "gpt-3.5-turbo"
        ])
        idx = self.openai_model_combo.findText(self.config.default_model_openai)
        if idx >= 0:
            self.openai_model_combo.setCurrentIndex(idx)
        openai_layout.addWidget(self.openai_model_combo, 1, 1, 1, 2)

        self.openai_web_search = QCheckBox("Enable Web Search (if supported by model)")
        self.openai_web_search.setChecked(True)
        openai_layout.addWidget(self.openai_web_search, 2, 0, 1, 3)

        openai_group.setLayout(openai_layout)
        layout.addWidget(openai_group)

        # Gemini Configuration
        gemini_group = QGroupBox("Google Gemini Configuration")
        gemini_layout = QGridLayout()

        gemini_layout.addWidget(QLabel("API Key:"), 0, 0)
        self.gemini_key_edit = QLineEdit()
        self.gemini_key_edit.setPlaceholderText("AI...")
        self.gemini_key_edit.setEchoMode(QLineEdit.Password)
        self.gemini_key_edit.setText(self.config.gemini_api_key)
        gemini_layout.addWidget(self.gemini_key_edit, 0, 1)

        self.gemini_show_key = QCheckBox("Show")
        self.gemini_show_key.toggled.connect(
            lambda checked: self.gemini_key_edit.setEchoMode(
                QLineEdit.Normal if checked else QLineEdit.Password
            )
        )
        gemini_layout.addWidget(self.gemini_show_key, 0, 2)

        gemini_layout.addWidget(QLabel("Model:"), 1, 0)
        self.gemini_model_combo = QComboBox()
        self.gemini_model_combo.addItems([
            "gemini-2.0-flash-exp",
            "gemini-1.5-flash",
            "gemini-1.5-flash-8b",
            "gemini-1.5-pro"
        ])
        idx = self.gemini_model_combo.findText(self.config.default_model_gemini)
        if idx >= 0:
            self.gemini_model_combo.setCurrentIndex(idx)
        gemini_layout.addWidget(self.gemini_model_combo, 1, 1, 1, 2)

        gemini_group.setLayout(gemini_layout)
        layout.addWidget(gemini_group)

        # Default provider
        provider_group = QGroupBox("Default Provider")
        provider_layout = QHBoxLayout()

        self.provider_group = QButtonGroup(self)
        self.openai_radio = QRadioButton("OpenAI (ChatGPT)")
        self.gemini_radio = QRadioButton("Google Gemini")
        self.provider_group.addButton(self.openai_radio)
        self.provider_group.addButton(self.gemini_radio)

        if self.config.default_provider == "gemini":
            self.gemini_radio.setChecked(True)
        else:
            self.openai_radio.setChecked(True)

        provider_layout.addWidget(self.openai_radio)
        provider_layout.addWidget(self.gemini_radio)
        provider_layout.addStretch()
        provider_group.setLayout(provider_layout)
        layout.addWidget(provider_group)

        # Save button
        btn_layout = QHBoxLayout()
        self.save_config_btn = QPushButton("Save Configuration")
        self.save_config_btn.clicked.connect(self._save_config)
        btn_layout.addWidget(self.save_config_btn)

        self.test_api_btn = QPushButton("Test API Connection")
        self.test_api_btn.clicked.connect(self._test_api)
        btn_layout.addWidget(self.test_api_btn)

        btn_layout.addStretch()
        layout.addLayout(btn_layout)

        layout.addStretch()
        widget.setLayout(layout)
        return widget

    def _create_analysis_tab(self) -> QWidget:
        """Create the analysis tab."""
        widget = QWidget()
        layout = QVBoxLayout()

        # Strain name input
        strain_layout = QHBoxLayout()
        strain_layout.addWidget(QLabel("Target Strain:"))
        self.strain_edit = QLineEdit()
        self.strain_edit.setPlaceholderText("e.g., Escherichia coli K-12, Saccharomyces cerevisiae S288C...")
        strain_layout.addWidget(self.strain_edit)
        layout.addLayout(strain_layout)

        # Item selection
        selection_group = QGroupBox("Select Items to Analyze")
        selection_layout = QVBoxLayout()

        # Analysis type
        type_layout = QHBoxLayout()
        self.analysis_type_group = QButtonGroup(self)
        self.reactions_radio = QRadioButton("Reactions")
        self.reactions_radio.setChecked(True)
        self.genes_radio = QRadioButton("Genes")
        self.both_radio = QRadioButton("Both")
        self.analysis_type_group.addButton(self.reactions_radio)
        self.analysis_type_group.addButton(self.genes_radio)
        self.analysis_type_group.addButton(self.both_radio)

        type_layout.addWidget(self.reactions_radio)
        type_layout.addWidget(self.genes_radio)
        type_layout.addWidget(self.both_radio)
        type_layout.addStretch()

        self.reactions_radio.toggled.connect(self._update_item_list)
        self.genes_radio.toggled.connect(self._update_item_list)
        self.both_radio.toggled.connect(self._update_item_list)

        selection_layout.addLayout(type_layout)

        # Search filter
        filter_layout = QHBoxLayout()
        filter_layout.addWidget(QLabel("Filter:"))
        self.filter_edit = QLineEdit()
        self.filter_edit.setPlaceholderText("Search by ID or name...")
        self.filter_edit.textChanged.connect(self._filter_items)
        filter_layout.addWidget(self.filter_edit)
        selection_layout.addLayout(filter_layout)

        # Item list
        list_layout = QHBoxLayout()

        # Available items
        available_widget = QWidget()
        available_layout = QVBoxLayout()
        available_layout.setContentsMargins(0, 0, 0, 0)
        available_layout.addWidget(QLabel("Available Items:"))
        self.available_list = QListWidget()
        self.available_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        available_layout.addWidget(self.available_list)
        available_widget.setLayout(available_layout)
        list_layout.addWidget(available_widget)

        # Add/Remove buttons
        btn_widget = QWidget()
        btn_layout_v = QVBoxLayout()
        btn_layout_v.addStretch()
        self.add_btn = QPushButton(">>")
        self.add_btn.clicked.connect(self._add_items)
        btn_layout_v.addWidget(self.add_btn)
        self.add_all_btn = QPushButton("Add All")
        self.add_all_btn.clicked.connect(self._add_all_items)
        btn_layout_v.addWidget(self.add_all_btn)
        self.remove_btn = QPushButton("<<")
        self.remove_btn.clicked.connect(self._remove_items)
        btn_layout_v.addWidget(self.remove_btn)
        self.clear_selection_btn = QPushButton("Clear")
        self.clear_selection_btn.clicked.connect(self._clear_selection)
        btn_layout_v.addWidget(self.clear_selection_btn)
        btn_layout_v.addStretch()
        btn_widget.setLayout(btn_layout_v)
        list_layout.addWidget(btn_widget)

        # Selected items
        selected_widget = QWidget()
        selected_layout = QVBoxLayout()
        selected_layout.setContentsMargins(0, 0, 0, 0)
        selected_layout.addWidget(QLabel("Selected for Analysis:"))
        self.selected_list = QListWidget()
        self.selected_list.setSelectionMode(QAbstractItemView.ExtendedSelection)
        selected_layout.addWidget(self.selected_list)
        selected_widget.setLayout(selected_layout)
        list_layout.addWidget(selected_widget)

        selection_layout.addLayout(list_layout)
        selection_group.setLayout(selection_layout)
        layout.addWidget(selection_group)

        # Progress
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)

        self.status_label = QLabel("")
        layout.addWidget(self.status_label)

        # Analysis buttons
        analysis_btn_layout = QHBoxLayout()
        self.start_analysis_btn = QPushButton("Start Analysis")
        self.start_analysis_btn.clicked.connect(self._start_analysis)
        analysis_btn_layout.addWidget(self.start_analysis_btn)

        self.stop_analysis_btn = QPushButton("Stop")
        self.stop_analysis_btn.clicked.connect(self._stop_analysis)
        self.stop_analysis_btn.setEnabled(False)
        analysis_btn_layout.addWidget(self.stop_analysis_btn)

        analysis_btn_layout.addStretch()
        layout.addLayout(analysis_btn_layout)

        widget.setLayout(layout)
        return widget

    def _create_results_tab(self) -> QWidget:
        """Create results tab."""
        widget = QWidget()
        layout = QVBoxLayout()

        # Results table
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(6)
        self.results_table.setHorizontalHeaderLabels([
            "ID", "Name", "Type", "Exists", "Confidence", "Evidence"
        ])
        self.results_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.results_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.results_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self.results_table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeToContents)
        self.results_table.horizontalHeader().setSectionResizeMode(4, QHeaderView.ResizeToContents)
        self.results_table.horizontalHeader().setSectionResizeMode(5, QHeaderView.Stretch)
        self.results_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.results_table.itemSelectionChanged.connect(self._show_result_details)
        layout.addWidget(self.results_table)

        # Details panel
        details_group = QGroupBox("Details")
        details_layout = QVBoxLayout()
        self.details_text = QTextEdit()
        self.details_text.setReadOnly(True)
        self.details_text.setMaximumHeight(200)
        details_layout.addWidget(self.details_text)
        details_group.setLayout(details_layout)
        layout.addWidget(details_group)

        # Export buttons
        export_layout = QHBoxLayout()
        self.export_json_btn = QPushButton("Export as JSON")
        self.export_json_btn.clicked.connect(self._export_json)
        export_layout.addWidget(self.export_json_btn)

        self.export_csv_btn = QPushButton("Export as CSV")
        self.export_csv_btn.clicked.connect(self._export_csv)
        export_layout.addWidget(self.export_csv_btn)

        self.clear_results_btn = QPushButton("Clear Results")
        self.clear_results_btn.clicked.connect(self._clear_results)
        export_layout.addWidget(self.clear_results_btn)

        export_layout.addStretch()
        layout.addLayout(export_layout)

        widget.setLayout(layout)
        return widget

    def _load_model_items(self):
        """Load reactions and genes from the model."""
        self._all_reactions = []
        self._all_genes = []

        model = self.appdata.project.cobra_py_model

        for rxn in model.reactions:
            self._all_reactions.append((rxn.id, rxn.name, "reaction"))

        for gene in model.genes:
            self._all_genes.append((gene.id, gene.name, "gene"))

        self._update_item_list()

    @Slot()
    def _update_item_list(self):
        """Update the available items list based on selection type."""
        self.available_list.clear()

        items = []
        if self.reactions_radio.isChecked():
            items = self._all_reactions
        elif self.genes_radio.isChecked():
            items = self._all_genes
        else:  # both
            items = self._all_reactions + self._all_genes

        for item_id, item_name, item_type in items:
            display_text = f"[{item_type[0].upper()}] {item_id}: {item_name}"
            item = QListWidgetItem(display_text)
            item.setData(Qt.UserRole, (item_id, item_name, item_type))
            self.available_list.addItem(item)

    @Slot()
    def _filter_items(self):
        """Filter available items based on search text."""
        filter_text = self.filter_edit.text().lower()

        for i in range(self.available_list.count()):
            item = self.available_list.item(i)
            data = item.data(Qt.UserRole)
            if data:
                item_id, item_name, _ = data
                visible = (filter_text in item_id.lower() or
                          filter_text in item_name.lower())
                item.setHidden(not visible)

    @Slot()
    def _add_items(self):
        """Add selected items to the analysis list."""
        selected = self.available_list.selectedItems()
        existing_ids = set()
        for i in range(self.selected_list.count()):
            data = self.selected_list.item(i).data(Qt.UserRole)
            if data:
                existing_ids.add(data[0])

        for item in selected:
            data = item.data(Qt.UserRole)
            if data and data[0] not in existing_ids:
                new_item = QListWidgetItem(item.text())
                new_item.setData(Qt.UserRole, data)
                self.selected_list.addItem(new_item)
                existing_ids.add(data[0])

    @Slot()
    def _add_all_items(self):
        """Add all visible items to the analysis list."""
        existing_ids = set()
        for i in range(self.selected_list.count()):
            data = self.selected_list.item(i).data(Qt.UserRole)
            if data:
                existing_ids.add(data[0])

        for i in range(self.available_list.count()):
            item = self.available_list.item(i)
            if not item.isHidden():
                data = item.data(Qt.UserRole)
                if data and data[0] not in existing_ids:
                    new_item = QListWidgetItem(item.text())
                    new_item.setData(Qt.UserRole, data)
                    self.selected_list.addItem(new_item)
                    existing_ids.add(data[0])

    @Slot()
    def _remove_items(self):
        """Remove selected items from the analysis list."""
        for item in self.selected_list.selectedItems():
            self.selected_list.takeItem(self.selected_list.row(item))

    @Slot()
    def _clear_selection(self):
        """Clear all selected items."""
        self.selected_list.clear()

    @Slot()
    def _save_config(self):
        """Save API configuration."""
        self.config.openai_api_key = self.openai_key_edit.text()
        self.config.gemini_api_key = self.gemini_key_edit.text()
        self.config.default_provider = "gemini" if self.gemini_radio.isChecked() else "openai"
        self.config.default_model_openai = self.openai_model_combo.currentText()
        self.config.default_model_gemini = self.gemini_model_combo.currentText()
        self.config.save_config()

        QMessageBox.information(self, "Saved", "API configuration saved successfully.")

    @Slot()
    def _test_api(self):
        """Test API connection."""
        provider = "gemini" if self.gemini_radio.isChecked() else "openai"

        self.test_api_btn.setEnabled(False)
        self.test_api_btn.setText("Testing...")
        QApplication.processEvents()

        try:
            if provider == "openai":
                import openai
                client = openai.OpenAI(api_key=self.openai_key_edit.text())
                response = client.chat.completions.create(
                    model=self.openai_model_combo.currentText(),
                    messages=[{"role": "user", "content": "Say 'API connection successful' in exactly those words."}],
                    max_tokens=20
                )
                QMessageBox.information(self, "Success", f"OpenAI API connection successful!\nModel: {self.openai_model_combo.currentText()}")
            else:
                import google.generativeai as genai
                genai.configure(api_key=self.gemini_key_edit.text())
                model = genai.GenerativeModel(self.gemini_model_combo.currentText())
                response = model.generate_content("Say 'API connection successful'")
                QMessageBox.information(self, "Success", f"Gemini API connection successful!\nModel: {self.gemini_model_combo.currentText()}")
        except ImportError as e:
            QMessageBox.critical(self, "Error", f"Required package not installed:\n{str(e)}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"API test failed:\n{str(e)}")
        finally:
            self.test_api_btn.setEnabled(True)
            self.test_api_btn.setText("Test API Connection")

    @Slot()
    def _start_analysis(self):
        """Start the analysis."""
        # Validate inputs
        strain_name = self.strain_edit.text().strip()
        if not strain_name:
            QMessageBox.warning(self, "Error", "Please enter a target strain name.")
            return

        if self.selected_list.count() == 0:
            QMessageBox.warning(self, "Error", "Please select at least one item to analyze.")
            return

        # Get provider and check API key
        provider = "gemini" if self.gemini_radio.isChecked() else "openai"

        if provider == "openai" and not self.openai_key_edit.text():
            QMessageBox.warning(self, "Error", "Please enter OpenAI API key.")
            return

        if provider == "gemini" and not self.gemini_key_edit.text():
            QMessageBox.warning(self, "Error", "Please enter Gemini API key.")
            return

        # Update config with current values
        self.config.openai_api_key = self.openai_key_edit.text()
        self.config.gemini_api_key = self.gemini_key_edit.text()

        # Collect items
        items = []
        for i in range(self.selected_list.count()):
            data = self.selected_list.item(i).data(Qt.UserRole)
            if data:
                items.append(data)

        # Get model
        if provider == "openai":
            model = self.openai_model_combo.currentText()
            use_web_search = self.openai_web_search.isChecked()
        else:
            model = self.gemini_model_combo.currentText()
            use_web_search = True  # Gemini always uses grounding when available

        # Start analysis thread
        self.analysis_thread = LLMAnalysisThread(
            self.config, items, strain_name, provider, model, use_web_search
        )
        self.analysis_thread.result_ready.connect(self._on_result)
        self.analysis_thread.error_occurred.connect(self._on_error)
        self.analysis_thread.progress_update.connect(self._on_progress)
        self.analysis_thread.finished.connect(self._on_finished)

        # Update UI
        self.start_analysis_btn.setEnabled(False)
        self.stop_analysis_btn.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, len(items))
        self.progress_bar.setValue(0)
        self.status_label.setText(f"Analyzing {len(items)} items...")

        self.analysis_thread.start()

    @Slot()
    def _stop_analysis(self):
        """Stop the analysis."""
        if self.analysis_thread:
            self.analysis_thread.stop()
            self.status_label.setText("Stopping...")

    @Slot(str, str, dict)
    def _on_result(self, item_id: str, item_type: str, result: dict):
        """Handle analysis result."""
        self.results[item_id] = {
            "item_type": item_type,
            **result
        }
        self._add_result_to_table(item_id, result)

    @Slot(str, str)
    def _on_error(self, item_id: str, error: str):
        """Handle analysis error."""
        self.results[item_id] = {
            "exists": "Error",
            "confidence": "N/A",
            "evidence": error,
            "error": True
        }
        self._add_result_to_table(item_id, self.results[item_id])

    @Slot(int, int)
    def _on_progress(self, current: int, total: int):
        """Handle progress update."""
        self.progress_bar.setValue(current)
        self.status_label.setText(f"Analyzing item {current}/{total}...")

    @Slot()
    def _on_finished(self):
        """Handle analysis completion."""
        self.start_analysis_btn.setEnabled(True)
        self.stop_analysis_btn.setEnabled(False)
        self.progress_bar.setVisible(False)
        self.status_label.setText(f"Analysis complete. {len(self.results)} items analyzed.")
        self.tabs.setCurrentWidget(self.results_tab)

    def _add_result_to_table(self, item_id: str, result: dict):
        """Add a result to the results table."""
        row = self.results_table.rowCount()
        self.results_table.insertRow(row)

        # Find item name
        item_name = ""
        item_type = result.get("item_type", "unknown")

        for rxn_id, rxn_name, _ in self._all_reactions:
            if rxn_id == item_id:
                item_name = rxn_name
                break

        if not item_name:
            for gene_id, gene_name, _ in self._all_genes:
                if gene_id == item_id:
                    item_name = gene_name
                    break

        self.results_table.setItem(row, 0, QTableWidgetItem(item_id))
        self.results_table.setItem(row, 1, QTableWidgetItem(item_name))
        self.results_table.setItem(row, 2, QTableWidgetItem(item_type))

        exists_item = QTableWidgetItem(result.get("exists", "Unknown"))
        exists = result.get("exists", "").lower()
        if exists == "yes":
            exists_item.setBackground(QColor(144, 238, 144))  # Light green
        elif exists == "no":
            exists_item.setBackground(QColor(255, 182, 193))  # Light pink
        elif exists == "likely":
            exists_item.setBackground(QColor(255, 255, 150))  # Light yellow
        self.results_table.setItem(row, 3, exists_item)

        self.results_table.setItem(row, 4, QTableWidgetItem(result.get("confidence", "Unknown")))

        evidence = result.get("evidence", "")
        if len(evidence) > 100:
            evidence = evidence[:100] + "..."
        self.results_table.setItem(row, 5, QTableWidgetItem(evidence))

        # Store full result in first column item
        self.results_table.item(row, 0).setData(Qt.UserRole, result)

    @Slot()
    def _show_result_details(self):
        """Show details for selected result."""
        selected = self.results_table.selectedItems()
        if not selected:
            return

        row = selected[0].row()
        result = self.results_table.item(row, 0).data(Qt.UserRole)

        if result:
            details = json.dumps(result, indent=2, ensure_ascii=False)
            self.details_text.setText(details)

    @Slot()
    def _export_json(self):
        """Export results as JSON."""
        if not self.results:
            QMessageBox.warning(self, "No Data", "No results to export.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export Results as JSON",
            self.appdata.work_directory,
            "JSON files (*.json)"
        )

        if file_path:
            if not file_path.endswith('.json'):
                file_path += '.json'

            with open(file_path, 'w', encoding='utf-8') as f:
                json.dump(self.results, f, indent=2, ensure_ascii=False)

            QMessageBox.information(self, "Exported", f"Results exported to:\n{file_path}")

    @Slot()
    def _export_csv(self):
        """Export results as CSV."""
        if not self.results:
            QMessageBox.warning(self, "No Data", "No results to export.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export Results as CSV",
            self.appdata.work_directory,
            "CSV files (*.csv)"
        )

        if file_path:
            if not file_path.endswith('.csv'):
                file_path += '.csv'

            import csv
            with open(file_path, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(["ID", "Type", "Exists", "Confidence", "Evidence", "References"])

                for item_id, result in self.results.items():
                    writer.writerow([
                        item_id,
                        result.get("item_type", ""),
                        result.get("exists", ""),
                        result.get("confidence", ""),
                        result.get("evidence", ""),
                        "; ".join(result.get("references", []))
                    ])

            QMessageBox.information(self, "Exported", f"Results exported to:\n{file_path}")

    @Slot()
    def _clear_results(self):
        """Clear all results."""
        self.results.clear()
        self.results_table.setRowCount(0)
        self.details_text.clear()
