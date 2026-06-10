"""The in out flux dialog"""

from qtpy.QtWidgets import QComboBox, QDialog, QHBoxLayout, QLabel, QPushButton, QVBoxLayout

from cnapy.appdata import AppData
from cnapy.utils import no_scroll


class InOutFluxDialog(QDialog):
    """A dialog to plot in out fluxes on a metabolite"""

    def __init__(self, appdata: AppData):
        QDialog.__init__(self)
        self.setWindowTitle("Compute in/out fluxes")

        self.appdata = appdata

        self.layout = QVBoxLayout()

        t1 = QLabel("Choose metabolite")
        self.layout.addWidget(t1)
        self.metabolite_chooser = no_scroll(QComboBox())
        self.layout.addWidget(self.metabolite_chooser)
        l = QHBoxLayout()
        self.button = QPushButton("Plot")
        self.cancel = QPushButton("Close")
        l.addWidget(self.button)
        l.addWidget(self.cancel)
        self.layout.addItem(l)
        self.setLayout(self.layout)

        # Connecting the signal
        self.cancel.clicked.connect(self.reject)
        self.button.clicked.connect(self.compute)
        self.update()

    def update(self):
        sorted_metabolite_ids = sorted([x.id for x in self.appdata.project.cobra_py_model.metabolites], reverse=True)
        for metabolite_id in sorted_metabolite_ids:
            self.metabolite_chooser.insertItem(0, metabolite_id)

    def compute(self):
        metabolite = self.metabolite_chooser.currentText()
        self.appdata.window.print_in_out_fluxes(metabolite)
