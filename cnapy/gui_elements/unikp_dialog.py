"""
unikp_dialog.py
===============
Sub-dialog for kcat prediction using UniKP.

Lets the user supply:
  - A FASTA file (protein sequences) or paste sequences manually
  - A SMILES CSV/TSV file (substrate SMILES) or paste manually

After prediction, injects the predicted kcat values into the project's
ECModelData.kcat_entries so they are picked up by the GECKO builder.
"""

from __future__ import annotations

import io

from qtpy.QtCore import Qt, QThread, Signal, Slot
from qtpy.QtWidgets import (
    QDialog,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from cnapy.appdata import AppData
from cnapy.ecmodel.unikp_predictor import check_unikp_available, predict_kcat_unikp


# ── worker thread ──────────────────────────────────────────────────────────────

class _UniKPWorker(QThread):
    """Runs UniKP prediction in a background thread."""

    progress = Signal(int, int)     # (current, total)
    finished = Signal(dict, list)   # (predictions, failed)
    error = Signal(str)

    def __init__(self, sequences: dict[str, str], smiles: dict[str, str]):
        super().__init__()
        self.sequences = sequences
        self.smiles = smiles

    def run(self):
        try:
            preds, failed = predict_kcat_unikp(
                self.sequences,
                self.smiles,
                progress_callback=lambda c, t: self.progress.emit(c, t),
            )
            self.finished.emit(preds, failed)
        except Exception as exc:  # noqa: BLE001
            self.error.emit(str(exc))


# ── dialog ─────────────────────────────────────────────────────────────────────

class UniKPDialog(QDialog):
    """Dialog for UniKP-based kcat prediction."""

    def __init__(self, appdata: AppData, parent=None):
        super().__init__(parent)
        self.appdata = appdata
        self.setWindowTitle("UniKP – kcat Prediction")
        self.setMinimumWidth(700)
        self.setMinimumHeight(560)
        self._worker: _UniKPWorker | None = None
        self._sequences: dict[str, str] = {}
        self._smiles: dict[str, str] = {}
        self._build_ui()
        self._check_deps()

    # ── UI ────────────────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)

        # --- Dependency status ---
        dep_group = QGroupBox("Dependency Status")
        dep_layout = QVBoxLayout()
        self.dep_label = QLabel("Checking…")
        self.dep_label.setWordWrap(True)
        dep_layout.addWidget(self.dep_label)
        dep_group.setLayout(dep_layout)
        root.addWidget(dep_group)

        # --- Sequences (FASTA) ---
        seq_group = QGroupBox("Protein Sequences (FASTA)")
        seq_layout = QVBoxLayout()

        seq_file_row = QHBoxLayout()
        self.seq_edit = QLineEdit()
        self.seq_edit.setPlaceholderText("Select FASTA file…")
        self.seq_edit.setReadOnly(True)
        seq_file_row.addWidget(self.seq_edit)
        seq_browse = QPushButton("Browse…")
        seq_browse.setAutoDefault(False)
        seq_browse.clicked.connect(self._browse_fasta)
        seq_file_row.addWidget(seq_browse)
        seq_layout.addLayout(seq_file_row)

        seq_layout.addWidget(QLabel(
            "Format: standard FASTA – header line starting with '>' followed by\n"
            "the identifier (used to match SMILES), then the amino-acid sequence."
        ))
        self.seq_text = QPlainTextEdit()
        self.seq_text.setPlaceholderText(
            ">enzyme_1\nMKILLTFLVVVTIAACSAKIEEGK…\n>enzyme_2\nMAGSSHHHHHH…"
        )
        self.seq_text.setMaximumHeight(130)
        seq_layout.addWidget(self.seq_text)
        self.seq_status = QLabel("Sequences loaded: 0")
        seq_layout.addWidget(self.seq_status)
        seq_group.setLayout(seq_layout)
        root.addWidget(seq_group)

        # --- SMILES ---
        smi_group = QGroupBox("Substrate SMILES")
        smi_layout = QVBoxLayout()

        smi_file_row = QHBoxLayout()
        self.smi_edit = QLineEdit()
        self.smi_edit.setPlaceholderText("Select CSV/TSV file with identifier, SMILES columns…")
        self.smi_edit.setReadOnly(True)
        smi_file_row.addWidget(self.smi_edit)
        smi_browse = QPushButton("Browse…")
        smi_browse.setAutoDefault(False)
        smi_browse.clicked.connect(self._browse_smiles)
        smi_file_row.addWidget(smi_browse)
        smi_layout.addLayout(smi_file_row)

        smi_layout.addWidget(QLabel(
            "Format: two-column file (no header required): identifier  SMILES\n"
            "The identifier must match the FASTA header identifier above."
        ))
        self.smi_text = QPlainTextEdit()
        self.smi_text.setPlaceholderText("enzyme_1\tCC(=O)Nc1ccc(O)cc1\nenzyme_2\tO=C(O)CCCC…")
        self.smi_text.setMaximumHeight(100)
        smi_layout.addWidget(self.smi_text)
        self.smi_status = QLabel("SMILES loaded: 0")
        smi_layout.addWidget(self.smi_status)
        smi_group.setLayout(smi_layout)
        root.addWidget(smi_group)

        # --- Progress bar ---
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        root.addWidget(self.progress_bar)

        # --- Result label ---
        self.result_label = QLabel("")
        self.result_label.setWordWrap(True)
        root.addWidget(self.result_label)

        # --- Buttons ---
        btn_row = QHBoxLayout()
        self.parse_btn = QPushButton("Parse Input")
        self.parse_btn.setAutoDefault(False)
        self.parse_btn.clicked.connect(self._parse_input)

        self.predict_btn = QPushButton("Run UniKP Prediction")
        self.predict_btn.setAutoDefault(False)
        self.predict_btn.setEnabled(False)
        self.predict_btn.clicked.connect(self._run_prediction)

        self.close_btn = QPushButton("Close")
        self.close_btn.setAutoDefault(False)
        self.close_btn.clicked.connect(self.close)

        btn_row.addWidget(self.parse_btn)
        btn_row.addWidget(self.predict_btn)
        btn_row.addStretch()
        btn_row.addWidget(self.close_btn)
        root.addLayout(btn_row)

    # ── helpers ───────────────────────────────────────────────────────────────

    def _check_deps(self):
        available, msg = check_unikp_available()
        if available:
            self.dep_label.setText("✓ All UniKP dependencies are available.")
            self.dep_label.setStyleSheet("color: green;")
        else:
            self.dep_label.setText(msg)
            self.dep_label.setStyleSheet("color: red;")
            self.predict_btn.setEnabled(False)
        self._deps_ok = available

    def _parse_fasta(self, text: str) -> dict[str, str]:
        """Parse FASTA-format text → {identifier: sequence}."""
        seqs = {}
        current_id = None
        current_seq: list[str] = []
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                if current_id:
                    current_seq.append(line)
        if current_id:
            seqs[current_id] = "".join(current_seq)
        return seqs

    def _parse_smiles_text(self, text: str) -> dict[str, str]:
        """Parse two-column (id, SMILES) text → {identifier: smiles}."""
        result = {}
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t", 1) if "\t" in line else line.split(None, 1)
            if len(parts) == 2:
                result[parts[0].strip()] = parts[1].strip()
        return result

    # ── slots ─────────────────────────────────────────────────────────────────

    @Slot()
    def _browse_fasta(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open FASTA file",
            self.appdata.work_directory,
            "FASTA files (*.fasta *.fa *.faa);;Text files (*.txt);;All files (*)",
            options=QFileDialog.DontUseNativeDialog,
        )
        if not path:
            return
        try:
            with open(path, encoding="utf-8", errors="replace") as fh:
                content = fh.read()
            self.seq_edit.setText(path)
            self.seq_text.setPlainText(content)
        except Exception as exc:
            QMessageBox.critical(self, "Error reading FASTA", str(exc))

    @Slot()
    def _browse_smiles(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open SMILES file",
            self.appdata.work_directory,
            "Data files (*.tsv *.csv *.txt);;All files (*)",
            options=QFileDialog.DontUseNativeDialog,
        )
        if not path:
            return
        try:
            with open(path, encoding="utf-8", errors="replace") as fh:
                content = fh.read()
            self.smi_edit.setText(path)
            self.smi_text.setPlainText(content)
        except Exception as exc:
            QMessageBox.critical(self, "Error reading SMILES file", str(exc))

    @Slot()
    def _parse_input(self):
        seq_text = self.seq_text.toPlainText().strip()
        smi_text = self.smi_text.toPlainText().strip()

        if seq_text:
            self._sequences = self._parse_fasta(seq_text)
        else:
            self._sequences = {}
        if smi_text:
            self._smiles = self._parse_smiles_text(smi_text)
        else:
            self._smiles = {}

        self.seq_status.setText(f"Sequences loaded: {len(self._sequences)}")
        self.smi_status.setText(f"SMILES loaded: {len(self._smiles)}")

        common = set(self._sequences) & set(self._smiles)
        self.result_label.setText(
            f"Matched pairs (sequence + SMILES): {len(common)}"
        )
        can_predict = bool(common) and self._deps_ok
        self.predict_btn.setEnabled(can_predict)

    @Slot()
    def _run_prediction(self):
        if not self._sequences or not self._smiles:
            QMessageBox.warning(self, "No input", "Please parse input first.")
            return

        common = set(self._sequences) & set(self._smiles)
        if not common:
            QMessageBox.warning(self, "No matched pairs",
                                "No identifiers match between sequences and SMILES.")
            return

        self.predict_btn.setEnabled(False)
        self.parse_btn.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, len(common))
        self.progress_bar.setValue(0)
        self.result_label.setText("Running UniKP prediction…")

        self._worker = _UniKPWorker(self._sequences, self._smiles)
        self._worker.progress.connect(self._on_progress)
        self._worker.finished.connect(self._on_finished)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    @Slot(int, int)
    def _on_progress(self, current: int, total: int):
        self.progress_bar.setMaximum(max(total, 1))
        self.progress_bar.setValue(current)

    @Slot(dict, list)
    def _on_finished(self, predictions: dict, failed: list):
        self.progress_bar.setVisible(False)
        self.predict_btn.setEnabled(True)
        self.parse_btn.setEnabled(True)

        if not predictions:
            self.result_label.setText("No predictions returned.")
            return

        # Inject predicted kcat values into ECModelData.kcat_entries
        ec = self.appdata.project.ec_model_data
        added = 0
        for uid, kcat in predictions.items():
            seq = self._sequences.get(uid, "")
            # Create a minimal kcat entry for each predicted value
            entry = {
                "proteins": [uid],
                "genes": [],
                "gene_name": uid,
                "kcat": kcat,
                "rxns": [],
                "notes": f"UniKP prediction (1/s)",
                "stoicho": [1],
            }
            # Only add if not already present
            existing_ids = {
                p
                for e in ec.kcat_entries
                for p in e.get("proteins", [])
            }
            if uid not in existing_ids:
                ec.kcat_entries.append(entry)
                added += 1

        msg = (
            f"Prediction complete.\n"
            f"  Predicted: {len(predictions)}\n"
            f"  Failed: {len(failed)}\n"
            f"  Added to kcat entries: {added}\n"
        )
        if failed:
            msg += f"\nFailed IDs: {', '.join(failed[:10])}"
            if len(failed) > 10:
                msg += f" … (+{len(failed)-10} more)"
        self.result_label.setText(msg)
        QMessageBox.information(self, "UniKP Prediction Done", msg)

    @Slot(str)
    def _on_error(self, message: str):
        self.progress_bar.setVisible(False)
        self.predict_btn.setEnabled(True)
        self.parse_btn.setEnabled(True)
        self.result_label.setText(f"Error: {message}")
        QMessageBox.critical(self, "UniKP Error", message)
