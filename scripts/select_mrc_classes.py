#!/usr/bin/env python3
"""Interactively select class averages from an MRC stack.

Usage:
    python scripts/select_mrc_classes.py /path/to/classes.mrcs

Controls inside the viewer:
    - Click image: toggle selected/unselected
    - n / right arrow / space: next page
    - p / left arrow: previous page
    - a: select all classes on current page
    - u: unselect all classes on current page
    - s: save selection file
    - q: save and quit

Output:
    Writes a text file next to the stack with the same base name and ".txt"
    extension. One row per class average:
        1 = selected
        0 = unselected
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List

import numpy as np

try:
    import mrcfile
except ImportError as exc:
    raise SystemExit(
        "Missing dependency: mrcfile. Install with: pip install mrcfile"
    ) from exc

try:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
except ImportError as exc:
    raise SystemExit(
        "Missing dependency: matplotlib. Install with: pip install matplotlib"
    ) from exc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Interactively select class averages from an MRC stack and save 1/0 rows."
    )
    parser.add_argument("mrc_stack", type=Path, help="Path to MRC/MRCS class-average stack")
    parser.add_argument(
        "--per-page",
        type=int,
        default=100,
        help="Number of class averages shown per page (default: 100)",
    )
    parser.add_argument(
        "--cols",
        type=int,
        default=10,
        help="Number of grid columns per page (default: 10)",
    )
    parser.add_argument(
        "--cmap",
        default="gray",
        help="Matplotlib colormap for display (default: gray)",
    )
    return parser.parse_args()


def load_stack(mrc_path: Path) -> np.ndarray:
    if not mrc_path.exists():
        raise FileNotFoundError(f"MRC stack does not exist: {mrc_path}")

    with mrcfile.open(mrc_path, permissive=True) as mrc:
        data = np.asarray(mrc.data)

    if data.ndim == 2:
        data = data[np.newaxis, :, :]
    elif data.ndim != 3:
        raise ValueError(f"Expected 2D or 3D MRC data, got shape {data.shape}")

    return data.astype(np.float32, copy=False)


def normalize_image(img: np.ndarray) -> np.ndarray:
    lo, hi = np.percentile(img, (1.0, 99.0))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo = float(np.min(img))
        hi = float(np.max(img))

    if hi <= lo:
        return np.zeros_like(img, dtype=np.float32)

    out = (img - lo) / (hi - lo)
    return np.clip(out, 0.0, 1.0)


class SelectorUI:
    def __init__(self, stack: np.ndarray, output_txt: Path, per_page: int, cols: int, cmap: str) -> None:
        self.stack = stack
        self.n_classes = stack.shape[0]
        self.output_txt = output_txt
        self.per_page = max(1, per_page)
        self.cols = max(1, cols)
        self.rows = max(1, math.ceil(self.per_page / self.cols))
        self.cmap = cmap

        self.selection = np.zeros(self.n_classes, dtype=np.int8)
        self.page = 0
        self.total_pages = max(1, math.ceil(self.n_classes / self.per_page))

        self.fig = None
        self.axes = None
        self.ax_to_index: Dict[object, int] = {}
        self.rectangles: Dict[int, Rectangle] = {}

        self._try_load_existing_selection()

    def _try_load_existing_selection(self) -> None:
        if not self.output_txt.exists():
            return

        rows: List[int] = []
        with self.output_txt.open("r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line not in {"0", "1"}:
                    return
                rows.append(int(line))

        if len(rows) != self.n_classes:
            return

        self.selection = np.asarray(rows, dtype=np.int8)
        print(f"Loaded existing selection from {self.output_txt}")

    def _page_indices(self) -> np.ndarray:
        start = self.page * self.per_page
        stop = min(start + self.per_page, self.n_classes)
        return np.arange(start, stop, dtype=int)

    def _selected_count(self) -> int:
        return int(np.count_nonzero(self.selection))

    def _status_title(self) -> str:
        return (
            f"{self.output_txt.name} | page {self.page + 1}/{self.total_pages} | "
            f"selected {self._selected_count()}/{self.n_classes}"
        )

    def _set_border(self, idx: int) -> None:
        rect = self.rectangles.get(idx)
        if rect is None:
            return
        rect.set_edgecolor("lime" if self.selection[idx] == 1 else "red")

    def _draw_page(self) -> None:
        self.ax_to_index.clear()
        self.rectangles.clear()

        indices = self._page_indices()
        flat_axes = self.axes.ravel()

        for ax in flat_axes:
            ax.clear()
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_axis_off()

        for slot, idx in enumerate(indices):
            ax = flat_axes[slot]
            ax.set_axis_on()
            ax.set_xticks([])
            ax.set_yticks([])

            img = normalize_image(self.stack[idx])
            ax.imshow(img, cmap=self.cmap, interpolation="nearest")
            ax.set_title(str(idx + 1), fontsize=9)

            rect = Rectangle(
                (0.01, 0.01),
                0.98,
                0.98,
                transform=ax.transAxes,
                fill=False,
                linewidth=2.0,
                edgecolor="lime" if self.selection[idx] == 1 else "red",
            )
            ax.add_patch(rect)

            self.ax_to_index[ax] = idx
            self.rectangles[idx] = rect

        self.fig.suptitle(self._status_title(), fontsize=11)
        self.fig.canvas.draw_idle()

    def save(self) -> None:
        self.output_txt.parent.mkdir(parents=True, exist_ok=True)
        with self.output_txt.open("w", encoding="utf-8") as fh:
            for val in self.selection:
                fh.write(f"{int(val)}\n")
        print(f"Saved selection to {self.output_txt}")

    def _toggle_index(self, idx: int) -> None:
        self.selection[idx] = 1 - self.selection[idx]
        self._set_border(idx)
        self.fig.suptitle(self._status_title(), fontsize=11)
        self.fig.canvas.draw_idle()

    def _set_current_page(self, value: int) -> None:
        self.page = max(0, min(self.total_pages - 1, value))
        self._draw_page()

    def _bulk_set_current_page(self, selected: int) -> None:
        for idx in self._page_indices():
            self.selection[idx] = selected
            self._set_border(int(idx))
        self.fig.suptitle(self._status_title(), fontsize=11)
        self.fig.canvas.draw_idle()

    def on_click(self, event) -> None:
        ax = event.inaxes
        if ax is None or ax not in self.ax_to_index:
            return
        idx = self.ax_to_index[ax]
        self._toggle_index(idx)

    def on_key(self, event) -> None:
        key = (event.key or "").lower()
        if key in {"n", "right", " "}:
            self._set_current_page(self.page + 1)
        elif key in {"p", "left"}:
            self._set_current_page(self.page - 1)
        elif key == "a":
            self._bulk_set_current_page(1)
        elif key == "u":
            self._bulk_set_current_page(0)
        elif key == "s":
            self.save()
        elif key == "q":
            self.save()
            plt.close(self.fig)

    def on_close(self, _event) -> None:
        self.save()

    def run(self) -> None:
        print("Controls: click=toggle, n/p=page, a=select page, u=unselect page, s=save, q=save+quit")

        self.fig, self.axes = plt.subplots(self.rows, self.cols, figsize=(1.7 * self.cols, 1.7 * self.rows))
        self.fig.canvas.mpl_connect("button_press_event", self.on_click)
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.fig.canvas.mpl_connect("close_event", self.on_close)

        self._draw_page()
        plt.tight_layout()
        plt.show()


def main() -> int:
    args = parse_args()
    mrc_path: Path = args.mrc_stack
    output_txt = mrc_path.with_suffix(".txt")

    stack = load_stack(mrc_path)
    if stack.shape[0] == 0:
        raise SystemExit(f"No images in stack: {mrc_path}")

    ui = SelectorUI(
        stack=stack,
        output_txt=output_txt,
        per_page=args.per_page,
        cols=args.cols,
        cmap=args.cmap,
    )
    ui.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
