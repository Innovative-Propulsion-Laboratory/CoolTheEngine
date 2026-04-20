import json
import shutil
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any

import matplotlib.backends.backend_pdf

from PyQt6.QtCore import Qt, QThread
from PyQt6.QtGui import QPixmap
from PyQt6.QtWidgets import (
    QApplication,
    QComboBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QGridLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from config_ui_specs import FIELD_SPECS, FieldSpec
from plotter import append_run_log_figure
from ui_file_drop_zone import FileDropZone
from ui_results_page import ResultsPage
from ui_solver_worker import SolverRunResult, SolverWorker


class ConfigWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("CoolTheEngine - Configuration UI")
        self.resize(1366, 768)

        self._spec_map = {spec.key: spec for spec in FIELD_SPECS}
        self._inputs: dict[str, QLineEdit | QComboBox] = {}
        self._enum_list_widgets: dict[str, QListWidget] = {}
        self._run_btn: QPushButton | None = None
        self._save_btn: QPushButton | None = None
        self._json_drop_zone: FileDropZone | None = None
        self._contour_drop_zone: FileDropZone | None = None
        self._contour_path: str | None = None
        self._stack: QStackedWidget | None = None
        self._results_page: ResultsPage | None = None
        self._worker_thread: QThread | None = None
        self._worker: SolverWorker | None = None
        self._last_run_dir: Path | None = None
        self._progress_bar: QProgressBar | None = None

        self._build_ui()

    def _build_ui(self) -> None:
        central = QWidget()
        main_layout = QVBoxLayout(central)
        line_edit_width = 320
        enum_list_height = 58
        enum_list_min_width = line_edit_width
        enum_row_spacing = 6
        combo_width_ratio = 0.7
        channel_image_width = 180
        channel_images_row_spacing = 8
        wall_image_width = 2 * channel_image_width + channel_images_row_spacing

        def build_image_label(image_path: str, width: int | None = None, height: int | None = None) -> QLabel:
            label = QLabel()
            label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            pix = QPixmap(image_path)
            if pix.isNull():
                label.setText(Path(image_path).name)
                return label

            if width is not None and height is not None:
                pix = pix.scaled(width, height,
                                 Qt.AspectRatioMode.KeepAspectRatio,
                                 Qt.TransformationMode.SmoothTransformation)
            elif width is not None:
                pix = pix.scaledToWidth(width, Qt.TransformationMode.SmoothTransformation)
            elif height is not None:
                pix = pix.scaledToHeight(height, Qt.TransformationMode.SmoothTransformation)

            label.setPixmap(pix)
            return label

        header_row = QWidget()
        header_layout = QHBoxLayout(header_row)
        header_layout.setContentsMargins(0, 2, 0, 2)
        header_layout.setSpacing(8)

        title_height = 56
        logo_width = 160

        logo_label = build_image_label("assets/logo_ipl.png", height=title_height)
        logo_label.setFixedHeight(title_height)
        logo_label.setMinimumWidth(logo_width)
        logo_label.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)

        title = QLabel("Cool The Engine - configuration")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title.setMinimumHeight(title_height)
        title.setStyleSheet("font-size: 30px; font-weight: 700; margin: 0;")

        right_spacer = QWidget()
        right_spacer.setFixedWidth(logo_width)

        header_layout.addWidget(logo_label)
        header_layout.addWidget(title, 1)
        header_layout.addWidget(right_spacer)
        main_layout.addWidget(header_row)

        drop_zones_row = QWidget()
        drop_zones_layout = QHBoxLayout(drop_zones_row)
        drop_zones_layout.setContentsMargins(0, 0, 0, 0)
        drop_zones_layout.setSpacing(12)

        self._json_drop_zone = FileDropZone("Drop configuration file (.json)", ".json")
        self._json_drop_zone.clicked.connect(self.open_json)
        self._json_drop_zone.fileDropped.connect(self.open_json_from_path)

        self._contour_drop_zone = FileDropZone("Drop engine contour file (.csv)", ".csv")
        self._contour_drop_zone.clicked.connect(self.open_contour)
        self._contour_drop_zone.fileDropped.connect(self.open_contour_from_path)

        drop_zones_layout.addWidget(self._json_drop_zone, 1)
        drop_zones_layout.addWidget(self._contour_drop_zone, 1)
        main_layout.addWidget(drop_zones_row)

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_content = QWidget()
        groups_grid = QGridLayout(scroll_content)
        groups_grid.setContentsMargins(0, 0, 0, 0)
        groups_grid.setHorizontalSpacing(18)
        groups_grid.setVerticalSpacing(12)
        groups_grid.setColumnStretch(0, 1)
        groups_grid.setColumnStretch(1, 1)

        left_column_widget = QWidget()
        left_column_layout = QVBoxLayout(left_column_widget)
        left_column_layout.setContentsMargins(0, 0, 0, 0)
        left_column_layout.setSpacing(12)

        group_order = ["Engine parameters", "Coolant parameters", "Wall parameters", "Channel sizing"]
        group_boxes: dict[str, QGroupBox] = {}
        group_forms: dict[str, QFormLayout] = {}
        subgroup_seen: dict[str, set[str]] = {}
        for group_name in group_order:
            group_box = QGroupBox(group_name)
            group_box.setStyleSheet("QGroupBox { font-weight: 700; }")

            if group_name == "Channel sizing":
                group_layout = QHBoxLayout(group_box)
                group_layout.setContentsMargins(8, 12, 8, 8)
                group_layout.setSpacing(10)

                channel_form_widget = QWidget()
                form = QFormLayout(channel_form_widget)
                form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
                form.setHorizontalSpacing(12)
                form.setVerticalSpacing(7)
                form.setFormAlignment(Qt.AlignmentFlag.AlignTop)
                group_layout.addWidget(channel_form_widget, 1)

                channel_images_widget = QWidget()
                channel_images_layout = QVBoxLayout(channel_images_widget)
                channel_images_layout.setContentsMargins(0, 0, 0, 0)
                channel_images_layout.setSpacing(0)

                top_images_row = QWidget()
                top_images_layout = QHBoxLayout(top_images_row)
                top_images_layout.setContentsMargins(0, 0, 0, 0)
                top_images_layout.setSpacing(channel_images_row_spacing)

                channel_loc_img = build_image_label("assets/channel_definition_locations.png", width=channel_image_width)
                channel_angle_img = build_image_label("assets/channel_angle_explaination.png", width=channel_image_width)
                channel_loc_img.setFixedWidth(channel_image_width)
                channel_angle_img.setFixedWidth(channel_image_width)

                top_images_layout.addWidget(channel_loc_img)
                top_images_layout.addWidget(channel_angle_img)

                wall_img = build_image_label("assets/wall_explaination.png", width=wall_image_width)
                wall_img.setFixedWidth(wall_image_width)

                channel_images_layout.addWidget(top_images_row, alignment=Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignRight)
                channel_images_layout.addWidget(wall_img, alignment=Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignRight)
                channel_images_layout.addStretch(1)
                group_layout.addWidget(channel_images_widget, 0)
            else:
                form = QFormLayout(group_box)
                form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
                form.setHorizontalSpacing(12)
                form.setVerticalSpacing(7)
                form.setFormAlignment(Qt.AlignmentFlag.AlignTop)

            group_boxes[group_name] = group_box
            group_forms[group_name] = form
            subgroup_seen[group_name] = set()

        left_column_layout.addWidget(group_boxes["Engine parameters"])
        left_column_layout.addWidget(group_boxes["Coolant parameters"])
        left_column_layout.addWidget(group_boxes["Wall parameters"])
        left_column_layout.addStretch(1)

        groups_grid.addWidget(left_column_widget, 0, 0)
        groups_grid.addWidget(group_boxes["Channel sizing"], 0, 1)

        for spec in FIELD_SPECS:
            if not spec.editable:
                continue

            form = group_forms[spec.group]

            if spec.subgroup is not None and spec.subgroup not in subgroup_seen[spec.group]:
                subgroup_title = QLabel(spec.subgroup + ":")
                if spec.group == "Channel sizing":
                    subgroup_title.setStyleSheet("font-size: 15px; font-weight: 700; color: #3F3F46; margin-top: 20px; margin-bottom: 10px;")
                else:
                    subgroup_title.setStyleSheet("font-weight: 600; color: #3F3F46; margin-top: 6px;")
                form.addRow(subgroup_title)
                subgroup_seen[spec.group].add(spec.subgroup)

            if spec.expected == "enum":
                combo = QComboBox()
                combo.setEditable(False)
                for enum_item in spec.enum_cls:  # type: ignore[arg-type]
                    combo.addItem(enum_item.value)
                if not spec.editable:
                    combo.setEnabled(False)
                input_widget = combo

                if spec.allow_list:
                    combo_width = int((line_edit_width - enum_row_spacing) * combo_width_ratio)
                    add_btn_width = (line_edit_width - enum_row_spacing) - combo_width
                    combo.setFixedWidth(combo_width)

                    add_btn = QPushButton("Add to list")
                    add_btn.setFixedWidth(add_btn_width)
                    add_btn.setEnabled(spec.editable)
                    add_btn.clicked.connect(lambda _checked=False, key=spec.key: self._add_enum_list_value(key))

                    selected_list = QListWidget()
                    selected_list.setFixedWidth(enum_list_min_width)
                    selected_list.setFixedHeight(enum_list_height)
                    selected_list.setEnabled(spec.editable)
                    self._enum_list_widgets[spec.key] = selected_list

                    selector_row = QWidget()
                    selector_layout = QHBoxLayout(selector_row)
                    selector_layout.setContentsMargins(0, 0, 0, 0)
                    selector_layout.setSpacing(enum_row_spacing)
                    selector_layout.addWidget(combo)
                    selector_layout.addWidget(add_btn)
                    selector_row.setFixedWidth(line_edit_width)

                    enum_container = QWidget()
                    enum_layout = QVBoxLayout(enum_container)
                    enum_layout.setContentsMargins(0, 0, 0, 0)
                    enum_layout.setSpacing(4)
                    enum_layout.addWidget(selector_row, alignment=Qt.AlignmentFlag.AlignLeft)
                    enum_layout.addWidget(selected_list)
                    enum_container.setFixedWidth(line_edit_width)
                    input_widget = enum_container
                else:
                    combo.setFixedWidth(line_edit_width)
            else:
                line = QLineEdit()
                line.setFixedWidth(line_edit_width)
                line.setReadOnly(not spec.editable)
                if not spec.editable and spec.default is not None:
                    line.setText(str(spec.default))
                line.textChanged.connect(self._update_run_button_state)
                line_row = QWidget()
                line_row_layout = QHBoxLayout(line_row)
                line_row_layout.setContentsMargins(0, 0, 0, 0)
                line_row_layout.setSpacing(0)
                line_row_layout.addStretch(1)
                line_row_layout.addWidget(line)
                input_widget = line_row

            if spec.expected == "enum":
                enum_row = QWidget()
                enum_row_layout = QHBoxLayout(enum_row)
                enum_row_layout.setContentsMargins(0, 0, 0, 0)
                enum_row_layout.setSpacing(0)
                enum_row_layout.addStretch(1)
                enum_row_layout.addWidget(input_widget)
                input_widget = enum_row

            label_widget = QLabel(spec.label)
            if spec.group == "Channel sizing" and spec.subgroup is not None:
                label_widget.setStyleSheet("margin-left: 16px;")

            form.addRow(label_widget, input_widget)
            if spec.expected == "enum":
                self._inputs[spec.key] = combo
                combo.currentIndexChanged.connect(self._update_run_button_state)
            else:
                self._inputs[spec.key] = line

        scroll_area.setWidget(scroll_content)
        main_layout.addWidget(scroll_area, 2)

        bottom_buttons = QHBoxLayout()
        bottom_buttons.addStretch(1)

        btn_save = QPushButton("Save JSON As...")
        btn_save.clicked.connect(self.save_json)
        btn_save.setMinimumWidth(240)
        btn_save.setMinimumHeight(44)
        btn_save.setStyleSheet(
            "QPushButton {"
            "background-color: #2A7FFF;"
            "color: white;"
            "font-weight: 600;"
            "padding: 8px 14px;"
            "border-radius: 6px;"
            "}"
            "QPushButton:disabled {"
            "background-color: #B7B7B7;"
            "color: #666666;"
            "}"
        )
        self._save_btn = btn_save
        bottom_buttons.addWidget(btn_save)

        btn_run = QPushButton("Run current configuration")
        btn_run.clicked.connect(self.run_current_configuration)
        btn_run.setMinimumWidth(240)
        btn_run.setMinimumHeight(44)
        btn_run.setStyleSheet(
            "QPushButton {"
            "background-color: #26A642;"
            "color: white;"
            "font-weight: 600;"
            "padding: 8px 14px;"
            "border-radius: 6px;"
            "}"
            "QPushButton:disabled {"
            "background-color: #B7B7B7;"
            "color: #666666;"
            "}"
        )
        self._run_btn = btn_run
        bottom_buttons.addWidget(btn_run)

        main_layout.addLayout(bottom_buttons)

        self._progress_bar = QProgressBar()
        self._progress_bar.setMinimum(0)
        self._progress_bar.setMaximum(100)
        self._progress_bar.setValue(0)
        self._progress_bar.setTextVisible(True)
        self._progress_bar.setFormat("Bulk run progress: %p%")
        self._progress_bar.setVisible(False)
        main_layout.addWidget(self._progress_bar)

        self._update_run_button_state()

        self._results_page = ResultsPage()
        self._results_page.backRequested.connect(self._show_configuration_page)
        self._results_page.saveRequested.connect(self.save_results)

        self._stack = QStackedWidget()
        self._stack.addWidget(central)
        self._stack.addWidget(self._results_page)
        self._stack.setCurrentIndex(0)

        self.setCentralWidget(self._stack)

    def open_json_from_path(self, path: str) -> None:
        if not path:
            return

        if Path(path).suffix.lower() != ".json":
            QMessageBox.critical(self, "Invalid file", "Only .json files are allowed.")
            if self._json_drop_zone is not None:
                self._json_drop_zone.set_state("error")
            return

        try:
            with open(path, "r", encoding="utf-8") as f:
                data = json.load(f)
        except Exception as exc:
            QMessageBox.critical(self, "Load error", f"Failed to load JSON:\n{exc}")
            if self._json_drop_zone is not None:
                self._json_drop_zone.set_state("error")
            return

        try:
            self._load_into_fields(data)
            if self._json_drop_zone is not None:
                self._json_drop_zone.set_selected_file(path)
        except Exception as exc:
            QMessageBox.critical(self, "Validation error", str(exc))
            if self._json_drop_zone is not None:
                self._json_drop_zone.set_state("error")

        self._update_run_button_state()

    def open_json(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Open configuration JSON", "", "JSON Files (*.json)")
        if not path:
            return
        self.open_json_from_path(path)

    def open_contour_from_path(self, path: str) -> None:
        if not path:
            return

        contour_path = Path(path)
        if contour_path.suffix.lower() != ".csv":
            QMessageBox.critical(self, "Invalid file", "Only .csv files are allowed for contour.")
            if self._contour_drop_zone is not None:
                self._contour_drop_zone.set_state("error")
            return

        if not contour_path.exists() or not contour_path.is_file():
            QMessageBox.critical(self, "Invalid file", "Contour file does not exist.")
            if self._contour_drop_zone is not None:
                self._contour_drop_zone.set_state("error")
            return

        self._contour_path = str(contour_path)
        if self._contour_drop_zone is not None:
            self._contour_drop_zone.set_selected_file(str(contour_path))
        self._update_run_button_state()

    def open_contour(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Open contour CSV", "", "CSV Files (*.csv)")
        if not path:
            return
        self.open_contour_from_path(path)

    def _load_into_fields(self, data: dict[str, Any]) -> None:
        for spec in FIELD_SPECS:
            if spec.key not in self._inputs:
                continue

            widget = self._inputs[spec.key]
            value = data.get(spec.key, spec.default)

            if value is None:
                if isinstance(widget, QLineEdit) and spec.editable:
                    widget.clear()
                continue

            if spec.expected == "enum":
                combo = widget
                enum_values = {e.value for e in spec.enum_cls}  # type: ignore[union-attr]

                if isinstance(value, list):
                    if not spec.allow_list:
                        raise ValueError(f"{spec.key} does not support list values.")
                    if len(value) == 0:
                        raise ValueError(f"{spec.key} list cannot be empty")
                    invalid = [str(v) for v in value if str(v) not in enum_values]
                    if invalid:
                        raise ValueError(f"Invalid enum value(s) for {spec.key}: {', '.join(invalid)}")
                    if spec.key in self._enum_list_widgets:
                        self._clear_enum_list(spec.key)
                        for v in value:
                            self._add_enum_list_item(spec.key, str(v))
                    first = str(value[0])
                    idx = combo.findText(first)
                    if idx >= 0:
                        combo.setCurrentIndex(idx)
                else:
                    if spec.key in self._enum_list_widgets:
                        self._clear_enum_list(spec.key)
                    value_str = str(value)
                    if value_str not in enum_values:
                        raise ValueError(f"Invalid enum value for {spec.key}: {value_str}")
                    idx = combo.findText(value_str)
                    if idx < 0:
                        raise ValueError(f"Invalid enum value for {spec.key}: {value_str}")
                    combo.setCurrentIndex(idx)
            else:
                if isinstance(value, list):
                    display_list = [self._convert_si_to_display_scalar(v, spec) for v in value]
                    text = json.dumps(display_list)
                else:
                    text = str(self._convert_si_to_display_scalar(value, spec))
                widget.setText(text)

    def _convert_scalar(self, raw: Any, spec: FieldSpec) -> int | float:
        if spec.expected == "int":
            if isinstance(raw, bool):
                raise ValueError(f"{spec.key} must be an int, got bool")
            if isinstance(raw, int):
                value = raw
            elif isinstance(raw, float) and raw.is_integer():
                value = int(raw)
            else:
                value = int(str(raw).strip())
        else:
            value = float(raw)

        value_si = value * spec.to_si_factor

        if spec.minimum is not None and value_si < spec.minimum:
            raise ValueError(f"{spec.key} must be >= {spec.minimum}")
        if spec.maximum is not None and value_si > spec.maximum:
            raise ValueError(f"{spec.key} must be <= {spec.maximum}")

        if spec.expected == "int" and spec.to_si_factor == 1.0:
            return int(value_si)

        return float(value_si)

    def _convert_si_to_display_scalar(self, raw: Any, spec: FieldSpec) -> int | float:
        value = float(raw)
        display_value = value / spec.to_si_factor
        if spec.expected == "int" and display_value.is_integer():
            return int(display_value)
        return display_value

    def _parse_field_value(self, spec: FieldSpec) -> Any:
        widget = self._inputs[spec.key]

        if spec.expected == "enum":
            combo = widget
            current = combo.currentText().strip()
            enum_values = {e.value for e in spec.enum_cls}  # type: ignore[union-attr]

            if spec.allow_list and spec.key in self._enum_list_widgets:
                list_widget = self._enum_list_widgets[spec.key]
                if list_widget.count() > 0:
                    values: list[str] = []
                    for i in range(list_widget.count()):
                        item = list_widget.item(i)
                        value = item.data(Qt.ItemDataRole.UserRole)
                        value_str = str(value)
                        if value_str not in enum_values:
                            raise ValueError(f"Invalid value for {spec.key}: {value_str}")
                        values.append(value_str)
                    return values

            if current not in enum_values:
                raise ValueError(f"Invalid value for {spec.key}: {current}")
            return current

        line = widget
        raw_text = line.text().strip()
        if raw_text == "":
            raise ValueError(f"{spec.key} cannot be empty")

        if spec.allow_list and raw_text.startswith("["):
            try:
                arr = json.loads(raw_text)
            except json.JSONDecodeError as exc:
                raise ValueError(f"{spec.key} list must use valid JSON array syntax") from exc

            if not isinstance(arr, list):
                raise ValueError(f"{spec.key} must be a list")
            if len(arr) == 0:
                raise ValueError(f"{spec.key} list cannot be empty")
            return [self._convert_scalar(x, spec) for x in arr]

        if not spec.allow_list and raw_text.startswith("["):
            raise ValueError(f"{spec.key} does not support list values")

        return self._convert_scalar(raw_text, spec)

    def _add_enum_list_value(self, key: str) -> None:
        spec = self._spec_map[key]
        combo = self._inputs[key]
        if not isinstance(combo, QComboBox):
            return
        value = combo.currentText().strip()
        enum_values = {e.value for e in spec.enum_cls}  # type: ignore[union-attr]
        if value not in enum_values:
            QMessageBox.critical(self, "Validation error", f"Invalid value for {key}: {value}")
            return
        self._add_enum_list_item(key, value)

    def _add_enum_list_item(self, key: str, value: str) -> None:
        list_widget = self._enum_list_widgets[key]

        for i in range(list_widget.count()):
            existing = list_widget.item(i).data(Qt.ItemDataRole.UserRole)
            if str(existing) == value:
                return

        item = QListWidgetItem(list_widget)
        item.setData(Qt.ItemDataRole.UserRole, value)

        row_widget = QWidget()
        row_layout = QHBoxLayout(row_widget)
        row_layout.setContentsMargins(4, 2, 4, 2)
        row_layout.setSpacing(6)

        label = QLabel(value)
        remove_btn = QPushButton("Remove")
        remove_btn.setMaximumWidth(90)
        remove_btn.clicked.connect(lambda _checked=False, key=key, it=item: self._remove_enum_list_item(key, it))

        row_layout.addWidget(label)
        row_layout.addStretch(1)
        row_layout.addWidget(remove_btn)

        item.setSizeHint(row_widget.sizeHint())
        list_widget.addItem(item)
        list_widget.setItemWidget(item, row_widget)
        self._update_run_button_state()

    def _remove_enum_list_item(self, key: str, item: QListWidgetItem) -> None:
        list_widget = self._enum_list_widgets[key]
        row = list_widget.row(item)
        if row >= 0:
            list_widget.takeItem(row)
        self._update_run_button_state()

    def _clear_enum_list(self, key: str) -> None:
        if key in self._enum_list_widgets:
            self._enum_list_widgets[key].clear()
        self._update_run_button_state()

    def _is_configuration_complete(self) -> bool:
        try:
            self.build_configuration()
            return True
        except Exception:
            return False

    def _get_run_block_reason(self) -> str | None:
        try:
            self.build_configuration()
        except Exception as exc:
            return f"Configuration is incomplete or invalid: {exc}"

        if self._contour_path is None:
            return "Contour file missing: drop a .csv contour file."

        contour_path = Path(self._contour_path)
        if contour_path.suffix.lower() != ".csv":
            return "Invalid contour file type: contour must be a .csv file."

        if not contour_path.is_file():
            return "Contour file not found: the selected .csv file does not exist anymore."

        return None

    def _update_run_button_state(self) -> None:
        is_complete = self._is_configuration_complete()
        run_block_reason = self._get_run_block_reason()
        is_running = self._worker_thread is not None and self._worker_thread.isRunning()
        can_run = (run_block_reason is None) and (not is_running)

        if self._run_btn is None:
            pass
        else:
            self._run_btn.setEnabled(can_run)
            if is_running:
                self._run_btn.setToolTip("Solver is currently running.")
            else:
                self._run_btn.setToolTip("" if can_run else run_block_reason)

        if self._save_btn is not None:
            self._save_btn.setEnabled(is_complete)

    def build_configuration(self) -> dict[str, Any]:
        config: dict[str, Any] = {}
        for spec in FIELD_SPECS:
            if not spec.editable or spec.key not in self._inputs:
                if spec.default is None:
                    raise ValueError(f"{spec.key} default is missing")
                config[spec.key] = spec.default
                continue

            config[spec.key] = self._parse_field_value(spec)

        return config

    def save_json(self) -> None:
        try:
            config = self.build_configuration()
        except Exception as exc:
            QMessageBox.critical(self, "Validation error", str(exc))
            return

        path, _ = QFileDialog.getSaveFileName(self, "Save configuration JSON", "config.json", "JSON Files (*.json)")
        if not path:
            return

        file_path = Path(path)
        if file_path.suffix.lower() != ".json":
            file_path = file_path.with_suffix(".json")

        try:
            with open(file_path, "w", encoding="utf-8") as f:
                json.dump(config, f, indent=4)
        except Exception as exc:
            QMessageBox.critical(self, "Save error", f"Failed to save JSON:\n{exc}")
            return

        QMessageBox.information(self, "Saved", f"Configuration saved:\n{file_path}")

    def run_current_configuration(self) -> str | None:
        try:
            config = self.build_configuration()
        except Exception as exc:
            QMessageBox.critical(self, "Validation error", str(exc))
            return None

        if self._contour_path is None:
            QMessageBox.critical(self, "Validation error", "Please select a contour CSV file before running.")
            return None

        contour_file = Path(self._contour_path)
        if contour_file.suffix.lower() != ".csv" or not contour_file.is_file():
            QMessageBox.critical(self, "Validation error", "Contour file must be an existing .csv file.")
            return None

        tmp_dir = Path(tempfile.mkdtemp(prefix="cte_run_"))
        config_path = tmp_dir / "input.json"
        contour_output_path = tmp_dir / contour_file.name

        try:
            with open(config_path, "w", encoding="utf-8") as f:
                json.dump(config, f, indent=4)
        except Exception as exc:
            QMessageBox.critical(self, "Run error", f"Failed to write temp configuration:\n{exc}")
            return None

        try:
            shutil.copy2(contour_file, contour_output_path)
        except Exception as exc:
            QMessageBox.critical(self, "Run error", f"Failed to copy contour file:\n{exc}")
            return None

        self._start_solver_worker(str(config_path), str(contour_output_path), str(tmp_dir))
        return str(tmp_dir)

    def _start_solver_worker(self, json_path: str, contour_path: str, run_dir: str) -> None:
        if self._worker_thread is not None and self._worker_thread.isRunning():
            QMessageBox.warning(self, "Solver already running", "Please wait for the current run to finish.")
            return

        if self._run_btn is not None:
            self._run_btn.setEnabled(False)
            self._run_btn.setText("Running...")

        self._worker_thread = QThread(self)
        self._worker = SolverWorker(json_path=json_path, contour_path=contour_path, temp_dir=run_dir)
        self._worker.moveToThread(self._worker_thread)

        self._worker_thread.started.connect(self._worker.run)
        self._worker.progressChanged.connect(self._on_solver_progress)
        self._worker.finished.connect(self._on_solver_finished)
        self._worker.failed.connect(self._on_solver_failed)
        self._worker.finished.connect(self._worker_thread.quit)
        self._worker.failed.connect(self._worker_thread.quit)
        self._worker_thread.finished.connect(self._cleanup_worker)

        self._worker_thread.start()

    def _cleanup_worker(self) -> None:
        if self._run_btn is not None:
            self._run_btn.setText("Run current configuration")
        self._update_run_button_state()

        if self._progress_bar is not None:
            self._progress_bar.setVisible(False)
            self._progress_bar.setValue(0)

        if self._worker is not None:
            self._worker.deleteLater()
            self._worker = None

        if self._worker_thread is not None:
            self._worker_thread.deleteLater()
            self._worker_thread = None

    def _on_solver_finished(self, result_obj: object) -> None:
        result = result_obj
        if not isinstance(result, SolverRunResult):
            QMessageBox.critical(self, "Run error", "Unexpected worker result type.")
            return

        self._last_run_dir = Path(result.temp_dir)

        figures_with_log = append_run_log_figure(result.figures, result.log_text)

        pdf_path = Path(result.output_dir) / "graphs.pdf"
        if len(figures_with_log) > 0:
            pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_path)
            for fig in figures_with_log:
                fig.savefig(pdf, format="pdf")
            pdf.close()

        if self._results_page is not None:
            self._results_page.set_figures(figures_with_log, bulk_run=result.bulk_run)

        if self._stack is not None:
            self._stack.setCurrentIndex(1)

    def _on_solver_failed(self, error_message: str) -> None:
        QMessageBox.critical(self, "Run error", f"Solver failed:\n{error_message}")

    def _on_solver_progress(self, current: int, total: int) -> None:
        if self._progress_bar is None:
            return
        if total <= 0:
            self._progress_bar.setVisible(False)
            return

        self._progress_bar.setVisible(True)
        self._progress_bar.setMinimum(0)
        self._progress_bar.setMaximum(total)
        self._progress_bar.setValue(max(0, min(current, total)))
        self._progress_bar.setFormat(f"Bulk run progress: {current}/{total} (%p%)")

    def _show_configuration_page(self) -> None:
        if self._stack is not None:
            self._stack.setCurrentIndex(0)

    def save_results(self) -> None:
        if self._last_run_dir is None or not self._last_run_dir.is_dir():
            QMessageBox.critical(self, "Save results", "No run results available to save.")
            return

        destination = QFileDialog.getExistingDirectory(self, "Select destination folder")
        if not destination:
            return

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        destination_path = Path(destination) / f"CTE_output_{timestamp}"
        destination_path.mkdir(parents=True, exist_ok=True)

        output_source_dir = self._last_run_dir / "output"
        if not output_source_dir.is_dir():
            QMessageBox.critical(self, "Save results", "No output directory found for the last run.")
            return

        def unique_destination(base_folder: Path, file_name: str) -> Path:
            candidate = base_folder / file_name
            if not candidate.exists():
                return candidate
            stem = Path(file_name).stem
            suffix = Path(file_name).suffix
            idx = 2
            while True:
                candidate = base_folder / f"{stem}_{idx}{suffix}"
                if not candidate.exists():
                    return candidate
                idx += 1

        copied_count = 0
        for ext in (".json", ".csv", ".pdf"):
            for src in output_source_dir.rglob(f"*{ext}"):
                dst = unique_destination(destination_path, src.name)
                shutil.copy2(src, dst)
                copied_count += 1

        QMessageBox.information(self, "Save results", f"Saved {copied_count} file(s) to:\n{destination_path}")


def main() -> None:
    app = QApplication([])
    window = ConfigWindow()
    window.showMaximized()
    app.exec()


if __name__ == "__main__":
    main()
