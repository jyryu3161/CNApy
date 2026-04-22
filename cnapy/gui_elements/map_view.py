"""The PyNetAnalyzer map view"""

import importlib.resources as resources
import math
from ast import literal_eval as make_tuple
from math import isclose

from qtpy.QtCore import QEvent, QMimeData, QRectF, Qt, Signal
from qtpy.QtGui import QColor, QDrag, QFont, QKeyEvent, QMouseEvent, QPainter, QPalette, QPen, QPixmap
from qtpy.QtSvg import QGraphicsSvgItem
from qtpy.QtWidgets import (
    QAction,
    QApplication,
    QCheckBox,
    QDoubleSpinBox,
    QGraphicsItem,
    QGraphicsPixmapItem,
    QGraphicsProxyWidget,
    QGraphicsScene,
    QGraphicsSceneDragDropEvent,
    QGraphicsSceneMouseEvent,
    QGraphicsView,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMenu,
    QMessageBox,
    QPinchGesture,
    QPushButton,
    QSlider,
    QTreeWidget,
    QVBoxLayout,
    QWidget,
)

from cnapy.appdata import AppData
from cnapy.gui_elements.box_position_dialog import BoxPositionDialog
from cnapy.utils import no_scroll

INCREASE_FACTOR = 1.1
DECREASE_FACTOR = 1 / INCREASE_FACTOR


class MapView(QGraphicsView):
    """A map of reaction boxes"""

    def __init__(self, appdata: AppData, central_widget, name: str):
        self.scene: QGraphicsScene = QGraphicsScene()
        QGraphicsView.__init__(self, self.scene)
        self.background: QGraphicsItem = None  # Can be QGraphicsSvgItem or QGraphicsPixmapItem
        palette = self.palette()
        if appdata.is_in_dark_mode:
            palette.setColor(QPalette.Base, QColor(90, 90, 90))  # Map etc. backgrounds
        else:
            palette.setColor(QPalette.Base, QColor(250, 250, 250))  # Map etc. backgrounds
        self.setPalette(palette)
        self.setInteractive(True)
        self.setDragMode(QGraphicsView.ScrollHandDrag)
        self.appdata = appdata
        self.central_widget = central_widget
        self.name: str = name
        self.setAcceptDrops(True)
        self.drag_map = False
        self.reaction_boxes: dict[str, ReactionBox] = {}
        self.constraint_slider = None
        self._zoom = 0
        self.previous_point = None
        self.select = False
        self.select_start = None
        self._last_mouse_pos = None

        # initial scale
        self._zoom = self.appdata.project.maps[self.name]["zoom"]
        if self._zoom > 0:
            for _ in range(1, self._zoom):
                self.scale(INCREASE_FACTOR, INCREASE_FACTOR)
        if self._zoom < 0:
            for _ in range(self._zoom, -1):
                self.scale(DECREASE_FACTOR, DECREASE_FACTOR)

        # connect events to methods
        self.horizontalScrollBar().valueChanged.connect(self.on_hbar_change)
        self.verticalScrollBar().valueChanged.connect(self.on_vbar_change)
        # enable pinch-zoom on trackpad
        self.grabGesture(Qt.PinchGesture)

        self.rebuild_scene()
        self.update()

    def on_hbar_change(self, x):
        self.appdata.project.maps[self.name]["pos"] = (x, self.verticalScrollBar().value())

    def on_vbar_change(self, y):
        self.appdata.project.maps[self.name]["pos"] = (self.horizontalScrollBar().value(), y)

    def dragEnterEvent(self, event: QGraphicsSceneDragDropEvent):
        self.previous_point = self.mapToScene(event.pos())
        event.acceptProposedAction()

    def dragMoveEvent(self, event: QGraphicsSceneDragDropEvent):
        event.setAccepted(True)
        point_item = self.mapToScene(event.pos())
        r_id = event.mimeData().text()

        if r_id in self.appdata.project.maps[self.name]["boxes"].keys():
            if isinstance(event.source(), QTreeWidget):  # existing/continued drag from reaction list
                self.appdata.project.maps[self.name]["boxes"][r_id] = (point_item.x(), point_item.y())
                self.mapChanged.emit(r_id)
            else:
                move_x = point_item.x() - self.previous_point.x()
                move_y = point_item.y() - self.previous_point.y()
                self.previous_point = point_item
                selected = self.scene.selectedItems()
                for item in selected:
                    pos = self.appdata.project.maps[self.name]["boxes"][item.id]

                    self.appdata.project.maps[self.name]["boxes"][item.id] = (pos[0] + move_x, pos[1] + move_y)
                    self.mapChanged.emit(item.id)

        else:  # drag reaction from list that has not yet a box on this map
            self.appdata.project.maps[self.name]["boxes"][r_id] = (point_item.x(), point_item.y())
            self.reactionAdded.emit(r_id)
            self.rebuild_scene()  # TODO don't rebuild the whole scene only add one item

        self.update()

    def dropEvent(self, event: QGraphicsSceneDragDropEvent):
        self.drag_map = False
        identifier = event.mimeData().text()
        self.mapChanged.emit(identifier)
        self.scene.setSceneRect(self.scene.itemsBoundingRect())
        self.viewport().setCursor(Qt.OpenHandCursor)
        self.update()

    def wheelEvent(self, event):
        modifiers = QApplication.queryKeyboardModifiers()
        if modifiers == Qt.ControlModifier:
            if event.angleDelta().y() > 0:
                self.appdata.project.maps[self.name]["box-size"] *= INCREASE_FACTOR
            else:
                self.appdata.project.maps[self.name]["box-size"] *= DECREASE_FACTOR

            self.mapChanged.emit("dummy")
            self.update()
        else:
            if event.angleDelta().y() > 0:
                self.zoom_in()
            else:
                self.zoom_out()

    def event(self, event):
        if event.type() == QEvent.Gesture:
            return self.gestureEvent(event)
        return super().event(event)

    def gestureEvent(self, event):
        pinch = event.gesture(Qt.PinchGesture)
        if pinch and isinstance(pinch, QPinchGesture):
            factor = pinch.scaleFactor()
            if factor > 0:
                self.scale(factor, factor)
                # keep scrollbar positions coherent
                self.appdata.project.maps[self.name]["zoom"] = self._zoom
            event.accept()
            return True
        return False

    def fit(self):
        self.fitInView(self.scene.sceneRect(), Qt.KeepAspectRatio)

    def zoom_in(self):
        self._zoom += 1

        self.appdata.project.maps[self.name]["zoom"] = self._zoom
        self.scale(INCREASE_FACTOR, INCREASE_FACTOR)

    def zoom_out(self):
        self._zoom -= 1

        self.appdata.project.maps[self.name]["zoom"] = self._zoom
        self.scale(DECREASE_FACTOR, DECREASE_FACTOR)

    def keyPressEvent(self, event: QKeyEvent):
        if not self.drag_map and event.key() in (Qt.Key_Control, Qt.Key_Shift):
            self.viewport().setCursor(Qt.ArrowCursor)
            self.select = True
        elif event.key() in (Qt.Key_Delete, Qt.Key_Backspace):
            # Remove selected reaction boxes
            selected_items = self.scene.selectedItems()
            for item in selected_items:
                if isinstance(item, ReactionBox):
                    self.remove_box(item.id)
            event.accept()
        else:
            super().keyPressEvent(event)

    def keyReleaseEvent(self, event: QKeyEvent):
        if (
            self.select
            and QApplication.mouseButtons() != Qt.LeftButton
            and event.key() in (Qt.Key_Control, Qt.Key_Shift)
        ):
            self.viewport().setCursor(Qt.OpenHandCursor)
            self.select = False
        else:
            super().keyReleaseEvent(event)

    def mousePressEvent(self, event: QMouseEvent):
        self._last_mouse_pos = event.pos()
        if self.select:  # select multiple boxes
            self.setDragMode(QGraphicsView.RubberBandDrag)  # switches to ArrowCursor
            self.select_start = self.mapToScene(event.pos())
            super(MapView, self).mousePressEvent(event)
            return

        # If clicking empty space, clear selection and start panning.
        target_item = self.itemAt(event.pos())
        if target_item is None:
            self.scene.clearSelection()
            self.scene.clearFocus()
            self.viewport().setCursor(Qt.ClosedHandCursor)
            self.setDragMode(QGraphicsView.ScrollHandDrag)
            self.drag_map = True
            super(MapView, self).mousePressEvent(event)
            return

        # Otherwise, let items handle the click; do not start map drag here.
        self.drag_map = False
        super(MapView, self).mousePressEvent(event)  # generates events for the graphics scene items

    def mouseMoveEvent(self, event: QMouseEvent):
        if self.drag_map and self._last_mouse_pos is not None:
            delta = event.pos() - self._last_mouse_pos
            self._last_mouse_pos = event.pos()
            self.horizontalScrollBar().setValue(self.horizontalScrollBar().value() - delta.x())
            self.verticalScrollBar().setValue(self.verticalScrollBar().value() - delta.y())
            event.accept()
            return
        super(MapView, self).mouseMoveEvent(event)

    def mouseReleaseEvent(self, event: QMouseEvent):
        if self.drag_map:
            self.viewport().setCursor(Qt.OpenHandCursor)
            self.drag_map = False
            self._last_mouse_pos = None
        if self.select:
            modifiers = QApplication.keyboardModifiers()
            if modifiers != Qt.ControlModifier and modifiers != Qt.ShiftModifier:
                self.viewport().setCursor(Qt.OpenHandCursor)
                self.select = False
            point = self.mapToScene(event.pos())

            width = point.x() - self.select_start.x()
            height = point.y() - self.select_start.y()
            selected = self.scene.items(QRectF(self.select_start.x(), self.select_start.y(), width, height))

            for item in selected:
                if isinstance(item, QGraphicsProxyWidget):
                    item.widget().parent.setSelected(True)

        painter = QPainter()
        self.render(painter)

        super(MapView, self).mouseReleaseEvent(event)
        event.accept()

    def focusOutEvent(self, event):
        super(MapView, self).focusOutEvent(event)
        self.viewport().setCursor(Qt.OpenHandCursor)
        self.select = False

    def enterEvent(self, event) -> None:
        super().enterEvent(event)
        if not isinstance(QApplication.focusWidget(), QLineEdit):
            # only take focus if no QlineEdit is active to prevent
            # editingFinished signals there
            if len(self.scene.selectedItems()) == 1:
                self.scene.selectedItems()[0].item.setFocus()
            else:
                self.scene.setFocus()  # to capture Shift/Ctrl keys

    def leaveEvent(self, event) -> None:
        super().leaveEvent(event)
        self.scene.clearFocus()  # finishes editing of potentially active ReactionBox

    def update_selected(self, found_ids):
        for r_id, box in self.reaction_boxes.items():
            box.item.setHidden(True)
            for found_id in found_ids:
                if found_id.lower() in r_id.lower():
                    box.item.setHidden(False)
                elif found_id.lower() in box.name.lower():
                    box.item.setHidden(False)

    def focus_reaction(self, reaction: str):
        x = self.appdata.project.maps[self.name]["boxes"][reaction][0]
        y = self.appdata.project.maps[self.name]["boxes"][reaction][1]
        self.centerOn(x, y)
        self.zoom_in_reaction()

    def zoom_in_reaction(self):
        bg_size = self.appdata.project.maps[self.name]["bg-size"]
        x = (INCREASE_FACTOR**self._zoom) / bg_size
        while x < 1:
            x = (INCREASE_FACTOR**self._zoom) / bg_size
            self._zoom += 1
            self.appdata.project.maps[self.name]["zoom"] = self._zoom
            self.scale(INCREASE_FACTOR, INCREASE_FACTOR)

    def highlight_reaction(self, string):
        treffer = self.reaction_boxes[string]
        treffer.item.setHidden(False)
        treffer.item.setFocus()

    def select_single_reaction(self, reac_id: str):
        box: ReactionBox = self.reaction_boxes.get(reac_id, None)
        if box is not None:
            self.scene.clearSelection()
            self.scene.clearFocus()
            box.setSelected(True)

    def set_background(self):
        if self.background is not None:
            self.scene.removeItem(self.background)
        bg_path = self.appdata.project.maps[self.name]["background"]
        bg_size = self.appdata.project.maps[self.name]["bg-size"]

        # Support both SVG and PNG files
        if bg_path.lower().endswith(".png"):
            pixmap = QPixmap(bg_path)
            self.background = QGraphicsPixmapItem(pixmap)
            self.background.setFlags(QGraphicsItem.ItemClipsToShape)
            self.background.setScale(bg_size)
        else:
            # Default to SVG
            self.background = QGraphicsSvgItem(bg_path)
            self.background.setFlags(QGraphicsItem.ItemClipsToShape)
            self.background.setScale(bg_size)
        self.scene.addItem(self.background)

    def rebuild_scene(self):
        self.scene.clear()
        self.background = None

        if (len(self.appdata.project.maps[self.name]["boxes"]) > 0) and self.appdata.project.maps[self.name][
            "background"
        ].replace("\\", "/").endswith("/data/default-bg.svg"):
            with resources.as_file(resources.files("cnapy") / "data" / "blank.svg") as path:
                self.appdata.project.maps[self.name]["background"] = str(path)

        self.set_background()

        for r_id in self.appdata.project.maps[self.name]["boxes"]:
            try:
                if r_id in self.appdata.project.cobra_py_model.reactions:
                    name = self.appdata.project.cobra_py_model.reactions.get_by_id(r_id).name
                else:
                    # Use reaction ID as name if not in model
                    name = r_id
                box = ReactionBox(self, r_id, name)

                self.scene.addItem(box)
                box.add_line_widget()
                self.reaction_boxes[r_id] = box
            except (KeyError, AttributeError) as e:
                print(f"failed to add reaction box for {r_id}: {e}")

    def delete_box(self, reaction_id: str) -> bool:
        box = self.reaction_boxes.get(reaction_id, None)
        if box is not None:
            lineedit = box.proxy
            self.scene.removeItem(lineedit)
            self.scene.removeItem(box)
            return True
        else:
            # print(f"Reaction {reaction_id} does not occur on map {self.name}")
            return False

    def update_reaction(self, old_reaction_id: str, new_reaction_id: str):
        if not self.delete_box(old_reaction_id):  # reaction is not on map
            return
        try:
            name = self.appdata.project.cobra_py_model.reactions.get_by_id(new_reaction_id).name
            box = ReactionBox(self, new_reaction_id, name)

            self.scene.addItem(box)
            box.add_line_widget()
            self.reaction_boxes[new_reaction_id] = box

            box.setScale(self.appdata.project.maps[self.name]["box-size"])
            box.proxy.setScale(self.appdata.project.maps[self.name]["box-size"])
            box.setPos(
                self.appdata.project.maps[self.name]["boxes"][box.id][0],
                self.appdata.project.maps[self.name]["boxes"][box.id][1],
            )

        except KeyError:
            print(f"Failed to add reaction box for {new_reaction_id} on map {self.name}")

    def update(self):
        for item in self.scene.items():
            if isinstance(item, QGraphicsSvgItem):
                item.setScale(self.appdata.project.maps[self.name]["bg-size"])
            elif isinstance(item, ReactionBox):
                item.setScale(self.appdata.project.maps[self.name]["box-size"])
                item.proxy.setScale(self.appdata.project.maps[self.name]["box-size"])
                try:
                    item.setPos(
                        self.appdata.project.maps[self.name]["boxes"][item.id][0],
                        self.appdata.project.maps[self.name]["boxes"][item.id][1],
                    )
                except KeyError:
                    print(f"{item.id} not found as box")
            else:
                pass

        self.set_values()
        self.recolor_all()

        # set scrollbars
        self.horizontalScrollBar().setValue(self.appdata.project.maps[self.name]["pos"][0])
        self.verticalScrollBar().setValue(self.appdata.project.maps[self.name]["pos"][1])

    def recolor_all(self):
        for r_id in self.appdata.project.maps[self.name]["boxes"]:
            box = self.reaction_boxes.get(r_id)
            if box is None:
                continue
            box.recolor()

    def show_flux_slider(self, reaction_box: "ReactionBox"):
        """Show or reuse the floating flux constraint slider for the given reaction box."""
        if self.constraint_slider is None:
            self.constraint_slider = FluxConstraintSlider(self)

        self.constraint_slider.configure(reaction_box)

        # place the slider near the clicked box in global coordinates
        scene_point = reaction_box.mapToScene(reaction_box.boundingRect().topRight())
        view_point = self.mapFromScene(scene_point)
        global_point = self.mapToGlobal(view_point)
        global_point.setX(global_point.x() + 10)
        global_point.setY(global_point.y() + 10)
        self.constraint_slider.move(global_point)
        self.constraint_slider.show()
        self.constraint_slider.raise_()

    def set_values(self):
        for r_id in self.appdata.project.maps[self.name]["boxes"]:
            box = self.reaction_boxes.get(r_id)
            if box is None:
                # box was not built (e.g. missing reaction in model); skip gracefully
                continue
            if r_id in self.appdata.project.scen_values.keys():
                box.set_value(self.appdata.project.scen_values[r_id])
            elif r_id in self.appdata.project.comp_values.keys():
                box.set_value(self.appdata.project.comp_values[r_id])
            else:
                box.item.setText("")

    def remove_box(self, reaction: str):
        self.delete_box(reaction)
        del self.appdata.project.maps[self.name]["boxes"][reaction]
        del self.reaction_boxes[reaction]
        self.update()
        self.reactionRemoved.emit(reaction)

    def value_changed(self, reaction: str, value: str):
        self.reactionValueChanged.emit(reaction, value)
        self.reaction_boxes[reaction].recolor()

    switchToReactionMask = Signal(str)
    maximizeReaction = Signal(str)
    minimizeReaction = Signal(str)
    setScenValue = Signal(str)
    reactionRemoved = Signal(str)
    reactionValueChanged = Signal(str, str)
    reactionAdded = Signal(str)
    mapChanged = Signal(str)
    broadcastReactionID = Signal(str)


class FluxConstraintSlider(QWidget):
    """Small floating slider to adjust a reaction's flux constraint."""

    def __init__(self, map_view: MapView):
        super().__init__(map_view)
        self.map_view = map_view
        self.reaction_box: ReactionBox = None
        self.scale = 10**self.map_view.appdata.rounding
        self.baseline_value = 0.0
        self.updating = False

        self.setWindowFlags(Qt.Popup)
        self.setLayout(QVBoxLayout())
        self.layout().setContentsMargins(8, 6, 8, 6)
        self.layout().setSpacing(6)

        header = QHBoxLayout()
        self.title = QLabel("")
        header.addWidget(self.title)
        header.addStretch()
        close_btn = QPushButton("X")
        close_btn.setFixedWidth(24)
        close_btn.clicked.connect(self.hide)
        header.addWidget(close_btn)
        self.layout().addLayout(header)

        self.bounds_label = QLabel("")
        self.layout().addWidget(self.bounds_label)

        self.slider = QSlider(Qt.Horizontal)
        self.slider.valueChanged.connect(self._handle_slider_change)
        self.layout().addWidget(self.slider)

        self.percent_toggle = QCheckBox("Adjust by %")
        self.percent_toggle.stateChanged.connect(self._handle_percent_toggle)
        self.layout().addWidget(self.percent_toggle)

        self.percent_slider = QSlider(Qt.Horizontal)
        self.percent_slider.setMinimum(0)
        self.percent_slider.setMaximum(200)
        self.percent_slider.setValue(100)
        self.percent_slider.valueChanged.connect(self._handle_percent_change)
        self.layout().addWidget(self.percent_slider)
        self.percent_slider.setEnabled(False)  # enabled on toggle

        percent_row = QHBoxLayout()
        self.percent_label = QLabel("100% of baseline")
        percent_row.addWidget(self.percent_label)
        self.layout().addLayout(percent_row)

        control_row = QHBoxLayout()
        self.spin = no_scroll(QDoubleSpinBox())
        self.spin.setDecimals(self.map_view.appdata.rounding)
        self.spin.valueChanged.connect(self._handle_spin_change)
        control_row.addWidget(self.spin)

        self.ko_btn = QPushButton("Set 0 (KO)")
        self.ko_btn.clicked.connect(lambda: self._apply_value(0.0))
        control_row.addWidget(self.ko_btn)

        self.clear_btn = QPushButton("Clear")
        self.clear_btn.clicked.connect(self._clear_constraint)
        control_row.addWidget(self.clear_btn)

        self.layout().addLayout(control_row)

    def configure(self, reaction_box: "ReactionBox"):
        self.reaction_box = reaction_box
        self.scale = 10**self.map_view.appdata.rounding

        # Check if reaction exists in model
        if reaction_box.id not in self.map_view.appdata.project.cobra_py_model.reactions:
            QMessageBox.warning(
                self.map_view,
                "Reaction not in model",
                f"Reaction '{reaction_box.id}' is not in the current model. Cannot set flux constraints.",
            )
            return

        reaction = self.map_view.appdata.project.cobra_py_model.reactions.get_by_id(reaction_box.id)
        lb, ub = reaction.bounds

        # Handle infinite bounds
        import math

        if math.isinf(lb) or math.isnan(lb):
            lb = -1000.0  # Default lower bound
        if math.isinf(ub) or math.isnan(ub):
            ub = 1000.0  # Default upper bound

        # make sure slider has a finite usable range
        slider_min = int(lb * self.scale)
        slider_max = int(ub * self.scale)
        if slider_min == slider_max:
            slider_min -= 1
            slider_max += 1

        self.slider.blockSignals(True)
        self.spin.blockSignals(True)
        self.slider.setMinimum(slider_min)
        self.slider.setMaximum(slider_max)
        self.spin.setRange(lb, ub)

        # prefer scenario value, then computed value, otherwise midpoint of bounds
        current = None
        scen_values = self.map_view.appdata.project.scen_values.get(reaction_box.id)
        comp_values = self.map_view.appdata.project.comp_values.get(reaction_box.id)
        if scen_values:
            current = (scen_values[0] + scen_values[1]) / 2
        elif comp_values:
            current = (comp_values[0] + comp_values[1]) / 2
        else:
            current = (lb + ub) / 2

        # If the baseline is 0 the percent control is meaningless, so fall back to abs_tol
        self.baseline_value = current if current is not None else 0.0
        if abs(self.baseline_value) < self.map_view.appdata.abs_tol:
            self.baseline_value = self.map_view.appdata.abs_tol
        self._sync_percent_slider(100)

        self.slider.setValue(int(round(current * self.scale)))
        self.spin.setValue(current)
        self.slider.blockSignals(False)
        self.spin.blockSignals(False)

        self.title.setText(f"{reaction.id} — {reaction.name}")
        self.bounds_label.setText(f"Bounds: [{lb}, {ub}]")
        self.percent_label.setText(f"100% of baseline ({self.baseline_value})")

        # Enable the percent slider by default and sync it
        self.percent_toggle.blockSignals(True)
        self.percent_toggle.setChecked(True)
        self.percent_slider.setEnabled(True)
        self.percent_toggle.blockSignals(False)
        self._sync_percent_from_value(current)

    def _handle_slider_change(self, value: int):
        if self.updating:
            return
        flux_value = value / self.scale
        self.spin.blockSignals(True)
        self.spin.setValue(flux_value)
        self.spin.blockSignals(False)
        self._sync_percent_from_value(flux_value)
        self._apply_value(flux_value)

    def _handle_spin_change(self, value: float):
        if self.updating:
            return
        self.slider.blockSignals(True)
        self.slider.setValue(int(round(value * self.scale)))
        self.slider.blockSignals(False)
        self._sync_percent_from_value(value)
        self._apply_value(value)

    def _handle_percent_toggle(self, state: int):
        use_percent = state == Qt.Checked
        self.percent_slider.setEnabled(use_percent)
        if use_percent:
            self._sync_percent_from_value(self.spin.value())

    def _handle_percent_change(self, value: int):
        if self.updating or self.percent_toggle.checkState() != Qt.Checked:
            return
        percent = value / 100.0
        target = self.baseline_value * percent
        self.updating = True
        self.spin.setValue(target)
        self.slider.setValue(int(round(target * self.scale)))
        self.updating = False
        self.percent_label.setText(f"{value}% of baseline ({self.baseline_value})")
        self._apply_value(target)

    def _sync_percent_slider(self, percent_int: int):
        self.updating = True
        self.percent_slider.setValue(percent_int)
        self.percent_label.setText(f"{percent_int}% of baseline ({self.baseline_value})")
        self.updating = False

    def _sync_percent_from_value(self, value: float):
        if self.baseline_value == 0:
            return
        percent_int = int(round((value / self.baseline_value) * 100))
        percent_int = max(self.percent_slider.minimum(), min(self.percent_slider.maximum(), percent_int))
        self._sync_percent_slider(percent_int)

    def _apply_value(self, value: float):
        if self.reaction_box is None:
            return
        formatted = self.map_view.appdata.format_flux_value(value)
        self.reaction_box.item.setText(formatted)
        self.reaction_box.item.setCursorPosition(0)
        self.reaction_box.value_changed()
        # Recompute immediately for real-time feedback
        self.map_view.central_widget.parent.run_auto_analysis()

    def _clear_constraint(self):
        if self.reaction_box is None:
            return
        self.reaction_box.item.setText("")
        self.reaction_box.value_changed()
        # Run pFBA instead of auto_analysis to restore original optimal flux distribution
        self.map_view.central_widget.parent.pfba()


class CLineEdit(QLineEdit):
    """A special line edit implementation for the use in ReactionBox"""

    def __init__(self, parent):
        self.parent: ReactionBox = parent
        self.accept_next_change_into_history = True
        super().__init__()

    def focusOutEvent(self, event):
        super().focusOutEvent(event)
        self.parent.setSelected(False)
        if self.isModified() and self.parent.map.appdata.auto_fba:
            self.parent.map.central_widget.parent.run_auto_analysis()
        self.parent.update()

    def focusInEvent(self, event):
        # is called before mousePressEvent
        super().focusInEvent(event)
        self.accept_next_change_into_history = True
        self.setModified(False)
        self.parent.setSelected(True)  # in case focus is regained via enterEvent of the map

    def mouseDoubleClickEvent(self, event):
        super().mouseDoubleClickEvent(event)
        self.parent.switch_to_reaction_mask()

    def mousePressEvent(self, event: QMouseEvent):
        # is called after focusInEvent
        super().mousePressEvent(event)
        if event.button() == Qt.MouseButton.LeftButton:
            if not self.parent.map.select:
                for bx in self.parent.map.reaction_boxes.values():
                    bx.setSelected(False)
            self.parent.setSelected(True)
            self.parent.broadcast_reaction_id()
            self.parent.map.show_flux_slider(self.parent)
        event.accept()


class ReactionBox(QGraphicsItem):
    """Handle to the line edits on the map"""

    def __init__(self, parent: MapView, r_id: str, name):
        QGraphicsItem.__init__(self)

        self.map = parent
        self.id = r_id
        self.name = name

        self.setFlag(QGraphicsItem.ItemIsSelectable)
        self.setFlag(QGraphicsItem.ItemIsMovable)
        self.setAcceptHoverEvents(True)
        self.item = CLineEdit(self)
        self.item.setTextMargins(1, -13, 0, -10)  # l t r b
        font = self.item.font()
        point_size = font.pointSize()
        font.setPointSizeF(point_size + 13.0)
        self.item.setFont(font)
        self.item.setAttribute(Qt.WA_TranslucentBackground)

        self.item.setFixedWidth(self.map.appdata.box_width)
        self.item.setMaximumHeight(self.map.appdata.box_height)
        self.item.setMinimumHeight(self.map.appdata.box_height)
        if r_id in self.map.appdata.project.cobra_py_model.reactions:
            r = self.map.appdata.project.cobra_py_model.reactions.get_by_id(r_id)
            text = (
                "Id: "
                + r.id
                + "\nName: "
                + r.name
                + "\nEquation: "
                + r.build_reaction_string()
                + "\nLowerbound: "
                + str(r.lower_bound)
                + "\nUpper bound: "
                + str(r.upper_bound)
                + "\nObjective coefficient: "
                + str(r.objective_coefficient)
            )
        else:
            # Reaction not in model
            text = f"Id: {r_id}\nName: {name}\n(Not in current model)"
        self.item.setToolTip(text)

        self.proxy = None  # proxy is set in add_line_widget after the item has been added

        self.set_default_style()

        self.setCursor(Qt.OpenHandCursor)
        self.setAcceptedMouseButtons(Qt.LeftButton)
        self.item.textEdited.connect(self.value_changed)
        self.item.returnPressed.connect(self.returnPressed)

        self.item.setContextMenuPolicy(Qt.CustomContextMenu)
        self.item.customContextMenuRequested.connect(self.on_context_menu)

        # create context menu
        self.pop_menu = QMenu(parent)
        toggle_knockout_action = QAction("Toggle Knockout", parent)
        self.pop_menu.addAction(toggle_knockout_action)
        toggle_knockout_action.triggered.connect(self.toggle_knockout)
        maximize_action = QAction("maximize flux for this reaction", parent)
        self.pop_menu.addAction(maximize_action)
        maximize_action.triggered.connect(self.emit_maximize_action)
        minimize_action = QAction("minimize flux for this reaction", parent)
        self.pop_menu.addAction(minimize_action)
        set_scen_value_action = QAction("add computed value to scenario", parent)
        set_scen_value_action.triggered.connect(self.emit_set_scen_value_action)
        self.pop_menu.addAction(set_scen_value_action)
        minimize_action.triggered.connect(self.emit_minimize_action)
        switch_action = QAction("switch to reaction mask", parent)
        self.pop_menu.addAction(switch_action)
        switch_action.triggered.connect(self.switch_to_reaction_mask)
        position_action = QAction("set box position...", parent)
        self.pop_menu.addAction(position_action)
        position_action.triggered.connect(self.position)
        remove_action = QAction("remove from map", parent)
        self.pop_menu.addAction(remove_action)
        remove_action.triggered.connect(self.remove)

        self.pop_menu.addSeparator()

    def mousePressEvent(self, event: QGraphicsSceneMouseEvent):
        if (event.button() == Qt.MouseButton.LeftButton) and (event.modifiers() == Qt.AltModifier):
            self.toggle_knockout()
            event.accept()
            return

        super().mousePressEvent(event)
        event.accept()
        if event.button() == Qt.MouseButton.LeftButton:
            if self.map.select:
                self.setSelected(not self.isSelected())
            else:
                self.setSelected(True)
            # Only show flux slider if reaction is in model
            if self.id in self.map.appdata.project.cobra_py_model.reactions:
                self.map.show_flux_slider(self)
        else:
            self.setCursor(Qt.ClosedHandCursor)

    def mouseReleaseEvent(self, event: QGraphicsSceneMouseEvent):
        event.accept()
        self.ungrabMouse()
        if not self.map.select:
            self.setCursor(Qt.OpenHandCursor)
            super().mouseReleaseEvent(event)  # here deselection of the other boxes occurs

    def hoverEnterEvent(self, event):
        if self.map.select:
            self.setCursor(Qt.ArrowCursor)
        else:
            self.setCursor(Qt.OpenHandCursor)
        super().hoverEnterEvent(event)

    def mouseMoveEvent(self, event: QGraphicsSceneMouseEvent):
        # Allow dragging the box if it's being moved
        if event.buttons() == Qt.LeftButton and not self.map.select:
            # Let QGraphicsItem handle the movement
            super().mouseMoveEvent(event)
            # Position will be saved in itemChange
        else:
            # Original drag and drop behavior for other cases
            event.accept()
            drag = QDrag(event.widget())
            mime = QMimeData()
            mime.setText(str(self.id))
            drag.setMimeData(mime)
            drag.exec()

    def itemChange(self, change, value):
        """Called when item properties change, including position."""
        if change == QGraphicsItem.ItemPositionHasChanged:
            # Save the new position to the map
            new_pos = self.pos()
            self.map.appdata.project.maps[self.map.name]["boxes"][self.id] = (new_pos.x(), new_pos.y())
            # Update proxy position to match
            if self.proxy is not None:
                self.proxy.setPos(new_pos)
            self.map.mapChanged.emit(self.id)
            self.map.central_widget.parent.unsaved_changes()
        return super().itemChange(change, value)

    def add_line_widget(self):
        self.proxy = self.map.scene.addWidget(self.item)
        self.proxy.show()

    def returnPressed(self):
        # self.item.clearFocus() # does not yet yield focus...
        self.proxy.clearFocus()  # ...but this does
        self.map.setFocus()
        self.item.accept_next_change_into_history = True  # reset so that next change will be recorded

    def handle_editing_finished(self):
        if self.item.isModified() and self.map.appdata.auto_fba:
            self.map.central_widget.parent.run_auto_analysis()

    # @Slot() # using the decorator gives a connection error?
    def value_changed(self):
        test = self.item.text().replace(" ", "")
        if test == "":
            if not self.item.accept_next_change_into_history:
                if len(self.map.appdata.scenario_past) > 0:
                    self.map.appdata.scenario_past.pop()  # replace previous change
            self.item.accept_next_change_into_history = False
            self.map.value_changed(self.id, test)
            self.set_default_style()
        elif validate_value(self.item.text()):
            if not self.item.accept_next_change_into_history:
                if len(self.map.appdata.scenario_past) > 0:
                    self.map.appdata.scenario_past.pop()  # replace previous change
            self.item.accept_next_change_into_history = False
            self.map.value_changed(self.id, self.item.text())
            if self.id in self.map.appdata.project.scen_values.keys():
                self.set_scen_style()
            else:
                self.set_comp_style()
        else:
            self.set_error_style()

    def set_default_style(self):
        """set the reaction box to error style"""
        palette = self.item.palette()
        role = self.item.backgroundRole()
        color = self.map.appdata.default_color
        color.setAlphaF(0.4)
        palette.setColor(role, color)
        role = self.item.foregroundRole()
        palette.setColor(role, Qt.black)
        self.item.setPalette(palette)

        self.set_font_style(QFont.StyleNormal)

    def set_error_style(self):
        """set the reaction box to error style"""
        self.set_color(Qt.white)
        self.set_fg_color(self.map.appdata.scen_color_bad)
        self.set_font_style(QFont.StyleOblique)

    def set_comp_style(self):
        self.set_color(self.map.appdata.comp_color)
        self.set_font_style(QFont.StyleNormal)

    def set_scen_style(self):
        self.set_color(self.map.appdata.scen_color)
        self.set_font_style(QFont.StyleNormal)

    def set_value(self, value: tuple[float, float]):
        """Sets the text of and reaction box according to the given value"""
        (vl, vu) = value
        if isclose(vl, vu, abs_tol=self.map.appdata.abs_tol):
            self.item.setText(str(round(float(vl), self.map.appdata.rounding)).rstrip("0").rstrip("."))
        else:
            self.item.setText(
                str(round(float(vl), self.map.appdata.rounding)).rstrip("0").rstrip(".")
                + ", "
                + str(round(float(vu), self.map.appdata.rounding)).rstrip("0").rstrip(".")
            )
        self.item.setCursorPosition(0)

    def recolor(self):
        value = self.item.text()
        test = value.replace(" ", "")
        if test == "":
            self.set_default_style()
        elif validate_value(value):
            if self.id in self.map.appdata.project.scen_values.keys():
                self.set_scen_style()
            elif self.id in self.map.appdata.project.comp_values:
                value = self.map.appdata.project.comp_values[self.id]
                (vl, vu) = value
                if math.isclose(vl, vu, abs_tol=self.map.appdata.abs_tol):
                    if self.map.appdata.modes_coloring:
                        if vl == 0:
                            self.set_color(Qt.red)
                        else:
                            self.set_color(Qt.green)
                    else:
                        self.set_comp_style()
                else:
                    if math.isclose(vl, 0.0, abs_tol=self.map.appdata.abs_tol):
                        self.set_color(self.map.appdata.special_color_1)
                    elif math.isclose(vu, 0.0, abs_tol=self.map.appdata.abs_tol):
                        self.set_color(self.map.appdata.special_color_1)
                    elif vl <= 0 and vu >= 0:
                        self.set_color(self.map.appdata.special_color_1)
                    else:
                        self.set_color(self.map.appdata.special_color_2)
            else:
                # Reaction not in model - use default style
                self.set_default_style()
        else:
            self.set_error_style()

    def set_color(self, color: QColor):
        palette = self.item.palette()
        role = self.item.backgroundRole()
        palette.setColor(role, color)
        role = self.item.foregroundRole()
        palette.setColor(role, Qt.black)
        self.item.setPalette(palette)

    def set_font_style(self, style: QFont.Style):
        font = self.item.font()
        font.setStyle(style)
        self.item.setFont(font)

    def set_fg_color(self, color: QColor):
        """set foreground color of the reaction box"""
        palette = self.item.palette()
        role = self.item.foregroundRole()
        palette.setColor(role, color)
        self.item.setPalette(palette)

    def boundingRect(self):
        return QRectF(-15, -15, self.map.appdata.box_width + 15 + 8, self.map.appdata.box_height + 15 + 8)

    def paint(self, painter: QPainter, _option, _widget: QWidget):
        # set color depending on wether the value belongs to the scenario
        if self.isSelected():
            light_blue = QColor(100, 100, 200)
            pen = QPen(light_blue)
            pen.setWidth(6)
            painter.setPen(pen)
            painter.drawRect(0 - 6, 0 - 6, self.map.appdata.box_width + 12, self.map.appdata.box_height + 12)

        if self.id in self.map.appdata.project.scen_values.keys():
            (vl, vu) = self.map.appdata.project.scen_values[self.id]
            ml = self.map.appdata.project.cobra_py_model.reactions.get_by_id(self.id).lower_bound
            mu = self.map.appdata.project.cobra_py_model.reactions.get_by_id(self.id).upper_bound

            if vu < ml or vl > mu:
                pen = QPen(self.map.appdata.scen_color_warn)
                painter.setBrush(self.map.appdata.scen_color_warn)
            else:
                pen = QPen(self.map.appdata.scen_color_good)
                painter.setBrush(self.map.appdata.scen_color_good)

            pen.setWidth(6)
            painter.setPen(pen)
            painter.drawRect(0 - 3, 0 - 3, self.map.appdata.box_width + 6, self.map.appdata.box_height + 6)

            pen.setWidth(1)
            painter.setPen(pen)
            painter.drawEllipse(-15, -15, 20, 20)

        else:
            painter.setPen(Qt.darkGray)
            painter.drawEllipse(-15, -15, 20, 20)

        painter.setPen(Qt.darkGray)
        painter.drawLine(-5, 0, -5, -10)
        painter.drawLine(0, -5, -10, -5)

        self.item.setFixedWidth(self.map.appdata.box_width)

    def setPos(self, x, y):
        self.proxy.setPos(x, y)
        super().setPos(x, y)

    def on_context_menu(self, point):
        # show context menu
        self.pop_menu.exec_(self.item.mapToGlobal(point))

    def position(self):
        position_dialog = BoxPositionDialog(self, self.map)
        position_dialog.exec()

    def remove(self):
        self.map.remove_box(self.id)
        self.map.drag = False

    def switch_to_reaction_mask(self):
        self.map.switchToReactionMask.emit(self.id)
        self.map.drag = False

    def toggle_knockout(self):
        if self.id in self.map.appdata.project.scen_values.keys() and self.map.appdata.project.scen_values[self.id] == (
            0,
            0,
        ):
            self.item.setText("")
        else:
            self.item.setText("0")
        self.value_changed()
        self.map.drag = False

    def emit_maximize_action(self):
        self.map.maximizeReaction.emit(self.id)
        self.map.drag = False

    def emit_set_scen_value_action(self):
        self.map.setScenValue.emit(self.id)
        self.map.drag = False

    def emit_minimize_action(self):
        self.map.minimizeReaction.emit(self.id)
        self.map.drag = False

    def broadcast_reaction_id(self):
        self.map.central_widget.broadcastReactionID.emit(self.id)
        self.map.drag = False


def validate_value(value):
    try:
        _x = float(value)
    except ValueError:
        try:
            (vl, vh) = make_tuple(value)
            if isinstance(vl, (int, float)) and isinstance(vh, (int, float)) and vl <= vh:
                return True
            else:
                return False
        except (ValueError, SyntaxError, TypeError):
            return False
    else:
        return True
