from __future__ import annotations

import contextlib
import io
from dataclasses import dataclass
from pathlib import Path

from PyQt6.QtCore import QObject, pyqtSignal
import matplotlib

matplotlib.use("Agg")

from CTE_oop import TaskManager


@dataclass
class SolverRunResult:
    temp_dir: str
    output_dir: str
    bulk_run: bool
    figures: list
    log_text: str


class SolverWorker(QObject):
    finished = pyqtSignal(object)
    failed = pyqtSignal(str)
    progressChanged = pyqtSignal(int, int)

    def __init__(self, json_path: str, contour_path: str, temp_dir: str) -> None:
        super().__init__()
        self._json_path = json_path
        self._contour_path = contour_path
        self._temp_dir = temp_dir

    def run(self) -> None:
        try:
            output_dir = Path(self._temp_dir) / "output"
            output_dir.mkdir(parents=True, exist_ok=True)

            log_buffer = io.StringIO()
            with contextlib.redirect_stdout(log_buffer), contextlib.redirect_stderr(log_buffer):
                manager = TaskManager(
                    self._json_path,
                    self._contour_path,
                    output_dir=str(output_dir),
                    progress_callback=self._on_progress,
                )
                bulk_run, figs = manager.run()

            if isinstance(figs, list):
                figure_list = figs
            elif figs is None:
                figure_list = []
            else:
                figure_list = [figs]

            result = SolverRunResult(
                temp_dir=self._temp_dir,
                output_dir=str(output_dir),
                bulk_run=bulk_run,
                figures=figure_list,
                log_text=log_buffer.getvalue(),
            )
            self.finished.emit(result)
        except Exception as exc:
            self.failed.emit(str(exc))

    def _on_progress(self, current: int, total: int) -> None:
        self.progressChanged.emit(current, total)
