from __future__ import annotations

from pathlib import Path

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import QFrame, QLabel, QVBoxLayout


class FileDropZone(QFrame):
    clicked = pyqtSignal()
    fileDropped = pyqtSignal(str)

    def __init__(self, prompt_text: str, extension: str) -> None:
        super().__init__()
        self._prompt_text = prompt_text
        self._extension = extension.lower()
        self._selected_file_path: str | None = None

        self.setAcceptDrops(True)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setMinimumHeight(120)
        self._persistent_state = "idle"
        self._is_hover = False
        self._is_drag_active = False

        self._label = QLabel(prompt_text, self)
        self._label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._label.setStyleSheet("font-size: 15px; font-weight: 600;")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 14, 12, 14)
        layout.addWidget(self._label)

        self.setObjectName("jsonDropZone")
        self.set_state("idle")

    def set_selected_file(self, path: str | None) -> None:
        self._selected_file_path = path
        if path:
            self._label.setText(f"Selected: {Path(path).name}")
            self.set_state("success")
        else:
            self._label.setText(self._prompt_text)
            self.set_state("idle")

    def set_state(self, state: str) -> None:
        if state in {"idle", "success", "error"}:
            self._persistent_state = state
            self._is_drag_active = False
        elif state == "active":
            self._is_drag_active = True

        self._refresh_style()

    def _refresh_style(self) -> None:
        if self._is_drag_active:
            self.setStyleSheet(
                "QFrame#jsonDropZone {"
                "background-color: #E8F4FF;"
                "border: 2px dashed #2A7FFF;"
                "border-radius: 10px;"
                "}"
            )
            return

        if self._is_hover:
            self.setStyleSheet(
                "QFrame#jsonDropZone {"
                "background-color: #F7F7F9;"
                "border: 2px solid #8BB8FF;"
                "border-radius: 10px;"
                "}"
            )
            return

        if self._persistent_state == "success":
            self.setStyleSheet(
                "QFrame#jsonDropZone {"
                "background-color: #EAF8EE;"
                "border: 2px solid #2E9E4D;"
                "border-radius: 10px;"
                "}"
            )
        elif self._persistent_state == "error":
            self.setStyleSheet(
                "QFrame#jsonDropZone {"
                "background-color: #FCEBEC;"
                "border: 2px solid #CC3D3D;"
                "border-radius: 10px;"
                "}"
            )
        else:
            self.setStyleSheet(
                "QFrame#jsonDropZone {"
                "background-color: #F7F7F9;"
                "border: 2px dashed #A9ADB6;"
                "border-radius: 10px;"
                "}"
            )

    def enterEvent(self, event) -> None:  # type: ignore[override]
        self._is_hover = True
        self._refresh_style()
        event.accept()

    def leaveEvent(self, event) -> None:  # type: ignore[override]
        self._is_hover = False
        self._refresh_style()
        event.accept()

    def dragEnterEvent(self, event) -> None:  # type: ignore[override]
        if event.mimeData().hasUrls():
            urls = event.mimeData().urls()
            if len(urls) == 1 and urls[0].toLocalFile().lower().endswith(self._extension):
                self.set_state("active")
                event.acceptProposedAction()
                return
        self.set_state("error")
        event.ignore()

    def dragLeaveEvent(self, event) -> None:  # type: ignore[override]
        self.set_state("idle")
        event.accept()

    def dropEvent(self, event) -> None:  # type: ignore[override]
        self.set_state("idle")
        if not event.mimeData().hasUrls():
            event.ignore()
            return

        urls = event.mimeData().urls()
        if len(urls) != 1:
            self.set_state("error")
            event.ignore()
            return

        local_path = urls[0].toLocalFile()
        if not local_path.lower().endswith(self._extension):
            self.set_state("error")
            event.ignore()
            return

        self.fileDropped.emit(local_path)
        event.acceptProposedAction()

    def mousePressEvent(self, event) -> None:  # type: ignore[override]
        if event.button() == Qt.MouseButton.LeftButton:
            self.clicked.emit()
        super().mousePressEvent(event)
