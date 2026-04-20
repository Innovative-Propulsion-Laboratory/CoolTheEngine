from __future__ import annotations

from PyQt6.QtCore import pyqtSignal
from PyQt6.QtGui import QWheelEvent
from PyQt6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QSizePolicy,
    QPushButton,
    QScrollArea,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar


class ScrollFriendlyFigureCanvas(FigureCanvas):
    """Figure canvas that scrolls the parent page when mouse wheel is used."""

    def wheelEvent(self, event: QWheelEvent) -> None:  # type: ignore[override]
        scroll_area = self._find_parent_scroll_area()
        if scroll_area is None:
            super().wheelEvent(event)
            return

        vbar = scroll_area.verticalScrollBar()

        pixel_delta = event.pixelDelta().y()
        if pixel_delta != 0:
            vbar.setValue(vbar.value() - pixel_delta)
            event.accept()
            return

        angle_delta = event.angleDelta().y()
        if angle_delta != 0:
            steps = angle_delta / 120.0
            scroll_amount = int(steps * vbar.singleStep() * 3)
            vbar.setValue(vbar.value() - scroll_amount)
            event.accept()
            return

        super().wheelEvent(event)

    def _find_parent_scroll_area(self) -> QScrollArea | None:
        parent = self.parentWidget()
        while parent is not None:
            if isinstance(parent, QScrollArea):
                return parent
            parent = parent.parentWidget()
        return None


class ResultsPage(QWidget):
    backRequested = pyqtSignal()
    saveRequested = pyqtSignal()

    def __init__(self) -> None:
        super().__init__()

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(10)

        header = QHBoxLayout()
        self._back_btn = QPushButton("Return to configuration")
        self._save_btn = QPushButton("Save results")
        self._title = QLabel("Solver results")
        self._title.setStyleSheet("font-size: 22px; font-weight: 700;")

        button_style = (
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
        self._back_btn.setStyleSheet(button_style)
        self._save_btn.setStyleSheet(button_style)
        self._back_btn.setMinimumHeight(40)
        self._save_btn.setMinimumHeight(40)
        self._back_btn.setMinimumWidth(220)
        self._save_btn.setMinimumWidth(180)

        self._back_btn.clicked.connect(self.backRequested.emit)
        self._save_btn.clicked.connect(self.saveRequested.emit)

        header.addWidget(self._back_btn)
        header.addStretch(1)
        header.addWidget(self._title)
        header.addStretch(1)
        header.addWidget(self._save_btn)

        layout.addLayout(header)

        self._tabs = QTabWidget()
        layout.addWidget(self._tabs, 1)

    def set_figures(self, figures: list, bulk_run: bool = False) -> None:
        self._tabs.clear()
        self._tabs.tabBar().show()

        if len(figures) == 0:
            empty = QWidget()
            empty_layout = QVBoxLayout(empty)
            empty_layout.addWidget(QLabel("No figures were returned by the solver."))
            empty_layout.addStretch(1)
            self._tabs.addTab(empty, "Results")
            return

        if bulk_run:
            # Bulk mode currently returns a single figure: show it directly without tabs.
            self._tabs.tabBar().hide()

            page = QWidget()
            page_layout = QVBoxLayout(page)
            page_layout.setContentsMargins(0, 0, 0, 0)
            page_layout.setSpacing(0)

            fig = figures[0]
            canvas = ScrollFriendlyFigureCanvas(fig)
            toolbar = NavigationToolbar(canvas, page)

            toolbar.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
            canvas.setMinimumSize(0, 0)

            page_layout.addWidget(toolbar)
            page_layout.addWidget(canvas, 1)
            self._tabs.addTab(page, "Results")
            return

        def idx_range(start: int, end: int) -> list[int]:
            return list(range(start, end + 1))

        group_specs: list[tuple[str, list[int]]] = [
            ("Input data", idx_range(1, 2)),
            ("Geometrical data", idx_range(3, 11)),
            ("Hot gas data", idx_range(12, 20) + idx_range(24, 25)),
            ("Heat fluxes", idx_range(21, 23) + idx_range(26, 27)),
            ("Coolant data", idx_range(28, 35)),
            ("Wall data", idx_range(36, 38)),
        ]

        any_tab_added = False
        for tab_name, figure_numbers in group_specs:
            selected = [n for n in figure_numbers if 1 <= n <= len(figures)]
            if len(selected) == 0:
                continue

            scroll = QScrollArea()
            scroll.setWidgetResizable(True)

            content = QWidget()
            content_layout = QVBoxLayout(content)
            content_layout.setContentsMargins(10, 10, 10, 10)
            content_layout.setSpacing(16)

            for fig_number in selected:
                fig = figures[fig_number - 1]

                block = QWidget()
                block_layout = QVBoxLayout(block)
                block_layout.setContentsMargins(0, 0, 0, 0)
                block_layout.setSpacing(0)

                block_title = QLabel(f"Figure {fig_number}")
                block_title.setStyleSheet("font-size: 14px; font-weight: 600; margin-bottom: 4px;")

                canvas = ScrollFriendlyFigureCanvas(fig)
                toolbar = NavigationToolbar(canvas, block)

                fig_w_in, fig_h_in = fig.get_size_inches()
                dpi = float(fig.dpi) if fig.dpi else 100.0
                canvas_width_px = max(int(fig_w_in * dpi), 900)
                canvas_height_px = max(int(fig_h_in * dpi), 600)
                canvas.setMinimumWidth(canvas_width_px)
                canvas.setMinimumHeight(canvas_height_px)
                canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
                toolbar.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
                block.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

                block_layout.addWidget(block_title)
                block_layout.addWidget(toolbar)
                block_layout.addWidget(canvas)
                content_layout.addWidget(block)

            content_layout.addStretch(1)
            scroll.setWidget(content)
            self._tabs.addTab(scroll, tab_name)
            any_tab_added = True

        if not any_tab_added:
            fallback = QWidget()
            fallback_layout = QVBoxLayout(fallback)
            fallback_layout.setContentsMargins(0, 0, 0, 0)
            fallback_layout.setSpacing(0)

            for idx, fig in enumerate(figures, start=1):
                canvas = ScrollFriendlyFigureCanvas(fig)
                toolbar = NavigationToolbar(canvas, fallback)
                fallback_layout.addWidget(QLabel(f"Figure {idx}"))
                fallback_layout.addWidget(toolbar)
                fallback_layout.addWidget(canvas)

            self._tabs.addTab(fallback, "Results")
