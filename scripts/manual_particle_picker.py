#!/usr/bin/env python3
from __future__ import annotations
import sys
from dataclasses import dataclass
from typing import Optional
import mrcfile
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.widgets import Button, Slider

@dataclass
class BoxItem:
    x: float
    y: float
    rect: Rectangle

class ParticlePicker:
    def __init__(self, img: np.ndarray, box_size: int) -> None:
        self.img = np.asarray(img, dtype=float)
        self.box_size = int(box_size)
        self.h, self.w = self.img.shape[:2]

        self.fig, self.ax = plt.subplots(figsize=(10, 10))
        self.fig.subplots_adjust(left=0.06, right=0.98, top=0.95, bottom=0.22)

        # Display range for the contrast slider.
        self.base_vmin, self.base_vmax = np.percentile(self.img[np.isfinite(self.img)], [1, 99])
        if not np.isfinite(self.base_vmin) or not np.isfinite(self.base_vmax) or self.base_vmin == self.base_vmax:
            self.base_vmin = float(np.nanmin(self.img))
            self.base_vmax = float(np.nanmax(self.img))
            if self.base_vmin == self.base_vmax:
                self.base_vmin -= 1.0
                self.base_vmax += 1.0

        self.contrast = 1.0
        self.im = self.ax.imshow(
            self.img,
            cmap="gray",
            origin="lower",
            vmin=self.base_vmin,
            vmax=self.base_vmax,
            interpolation="nearest",
        )
        self.ax.set_xlim(0, self.w)
        self.ax.set_ylim(0, self.h)
        self.ax.set_title(
            "Left click: add box | Left drag: move box | Right drag: pan | Wheel: zoom | Delete: remove selected",
            fontsize=10,
        )

        self.boxes: list[BoxItem] = []
        self.selected: Optional[BoxItem] = None
        self.dragging_box = False
        self.drag_offset = (0.0, 0.0)
        self.panning = False
        self.pan_start = (0.0, 0.0)
        self.pan_xlim = (0.0, 1.0)
        self.pan_ylim = (0.0, 1.0)
        self.closed = False

        # Widgets.
        slider_ax = self.fig.add_axes([0.12, 0.10, 0.62, 0.03])
        self.contrast_slider = Slider(
            slider_ax,
            "Contrast",
            0.5,
            5.0,
            valinit=self.contrast,
            valstep=0.01,
        )
        self.contrast_slider.on_changed(self.on_contrast_change)

        button_ax = self.fig.add_axes([0.79, 0.075, 0.16, 0.08])
        self.delete_button = Button(button_ax, "Delete selected")
        self.delete_button.on_clicked(self.delete_selected)

        self.status_ax = self.fig.add_axes([0.12, 0.045, 0.62, 0.02])
        self.status_ax.axis("off")
        self.status_text = self.status_ax.text(
            0.0,
            0.5,
            "No box selected",
            va="center",
            ha="left",
            fontsize=9,
        )

        # Events.
        self.fig.canvas.mpl_connect("button_press_event", self.on_press)
        self.fig.canvas.mpl_connect("button_release_event", self.on_release)
        self.fig.canvas.mpl_connect("motion_notify_event", self.on_motion)
        self.fig.canvas.mpl_connect("scroll_event", self.on_scroll)
        self.fig.canvas.mpl_connect("key_press_event", self.on_key)
        self.fig.canvas.mpl_connect("close_event", self.on_close)

    def clamp_center(self, x: float, y: float) -> tuple[float, float]:
        half = self.box_size / 2.0
        x = float(np.clip(x, half, max(half, self.w - half)))
        y = float(np.clip(y, half, max(half, self.h - half)))
        return x, y

    def current_bounds(self, x: float, y: float) -> tuple[float, float, float, float]:
        half = self.box_size / 2.0
        return x - half, y - half, self.box_size, self.box_size

    def set_selected(self, item: Optional[BoxItem]) -> None:
        if self.selected is item:
            return

        if self.selected is not None:
            self.selected.rect.set_edgecolor("yellow")
            self.selected.rect.set_linewidth(1.5)

        self.selected = item

        if self.selected is not None:
            self.selected.rect.set_edgecolor("red")
            self.selected.rect.set_linewidth(2.2)
            self.status_text.set_text(
                f"Selected box at x={self.selected.x:.1f}, y={self.selected.y:.1f}"
            )
        else:
            self.status_text.set_text("No box selected")

        self.fig.canvas.draw_idle()

    def add_box(self, x: float, y: float) -> None:
        x, y = self.clamp_center(x, y)
        left, bottom, width, height = self.current_bounds(x, y)
        rect = Rectangle(
            (left, bottom),
            width,
            height,
            linewidth=1.5,
            edgecolor="yellow",
            facecolor="none",
        )
        self.ax.add_patch(rect)
        item = BoxItem(x=x, y=y, rect=rect)
        self.boxes.append(item)
        self.set_selected(item)
        self.fig.canvas.draw_idle()

    def box_under_cursor(self, event) -> Optional[BoxItem]:
        for item in reversed(self.boxes):
            hit, _ = item.rect.contains(event)
            if hit:
                return item
        return None

    def update_box_position(self, item: BoxItem, x: float, y: float) -> None:
        x, y = self.clamp_center(x, y)
        item.x = x
        item.y = y
        left, bottom, width, height = self.current_bounds(x, y)
        item.rect.set_xy((left, bottom))
        item.rect.set_width(width)
        item.rect.set_height(height)
        if self.selected is item:
            self.status_text.set_text(f"Selected box at x={x:.1f}, y={y:.1f}")
        self.fig.canvas.draw_idle()

    def delete_selected(self, event=None) -> None:
        if self.selected is None:
            return

        item = self.selected
        try:
            item.rect.remove()
        except ValueError:
            pass

        self.boxes = [b for b in self.boxes if b is not item]
        self.selected = None
        self.status_text.set_text("No box selected")
        self.fig.canvas.draw_idle()

    def on_contrast_change(self, value: float) -> None:
        self.contrast = float(value)
        center = (self.base_vmin + self.base_vmax) / 2.0
        half_span = (self.base_vmax - self.base_vmin) / 2.0
        half_span = half_span / max(self.contrast, 1e-6)
        self.im.set_clim(center - half_span, center + half_span)
        self.fig.canvas.draw_idle()

    def on_press(self, event) -> None:
        if event.inaxes != self.ax or event.xdata is None or event.ydata is None:
            return

        if event.button == 3:
            self.panning = True
            self.pan_start = (event.xdata, event.ydata)
            self.pan_xlim = self.ax.get_xlim()
            self.pan_ylim = self.ax.get_ylim()
            return

        if event.button != 1:
            return

        item = self.box_under_cursor(event)
        if item is not None:
            self.set_selected(item)
            self.dragging_box = True
            self.drag_offset = (item.x - event.xdata, item.y - event.ydata)
        else:
            self.add_box(event.xdata, event.ydata)
            self.dragging_box = True
            self.drag_offset = (0.0, 0.0)

    def on_motion(self, event) -> None:
        if event.inaxes != self.ax or event.xdata is None or event.ydata is None:
            return

        if self.panning:
            dx = event.xdata - self.pan_start[0]
            dy = event.ydata - self.pan_start[1]
            self.ax.set_xlim(self.pan_xlim[0] - dx, self.pan_xlim[1] - dx)
            self.ax.set_ylim(self.pan_ylim[0] - dy, self.pan_ylim[1] - dy)
            self.fig.canvas.draw_idle()
            return

        if self.dragging_box and self.selected is not None:
            x = event.xdata + self.drag_offset[0]
            y = event.ydata + self.drag_offset[1]
            self.update_box_position(self.selected, x, y)

    def on_release(self, event) -> None:
        if event.button == 3:
            self.panning = False
        if event.button == 1:
            self.dragging_box = False

    def on_scroll(self, event) -> None:
        if event.inaxes != self.ax or event.xdata is None or event.ydata is None:
            return

        scale = 1.2
        if event.button == "up":
            zoom = 1.0 / scale
        elif event.button == "down":
            zoom = scale
        else:
            return

        cur_xlim = self.ax.get_xlim()
        cur_ylim = self.ax.get_ylim()
        xdata = event.xdata
        ydata = event.ydata

        new_xlim = [xdata - (xdata - cur_xlim[0]) * zoom, xdata + (cur_xlim[1] - xdata) * zoom]
        new_ylim = [ydata - (ydata - cur_ylim[0]) * zoom, ydata + (cur_ylim[1] - ydata) * zoom]

        self.ax.set_xlim(new_xlim)
        self.ax.set_ylim(new_ylim)
        self.fig.canvas.draw_idle()

    def on_key(self, event) -> None:
        if event.key in {"delete", "backspace"}:
            self.delete_selected()

    def on_close(self, event) -> None:
        self.closed = True

    def save_boxes(self, output_file: str = "positions_all.box") -> None:
        with open(output_file, "w", encoding="utf-8") as f:
            for item in self.boxes:
                x_bl = int(round(item.x - self.box_size / 2.0))
                y_bl = int(round(item.y - self.box_size / 2.0))
                f.write(f"{x_bl}\t{y_bl}\t{self.box_size}\t{self.box_size}\n")


def pick_particles_with_input() -> None:
    file_path = input("Enter the path to your MRC file: ").strip()
    try:
        box_size = int(input("Enter the box size for picking (e.g., 64, 128, 256): "))
    except ValueError:
        print("Invalid input. Using default box size of 64.")
        box_size = 64

    try:
        with mrcfile.open(file_path, mode="r", permissive=True) as mrc:
            data = mrc.data.copy()
            if data.ndim == 3:
                img = data[data.shape[0] // 2]
            else:
                img = data
    except Exception as e:
        print(f"Error: Could not open file. {e}")
        return

    picker = ParticlePicker(img=img, box_size=box_size)
    plt.show()

    picker.save_boxes("positions_all.box")

    print(f"\nFinished! Total particles picked: {len(picker.boxes)}")
    for i, item in enumerate(picker.boxes, start=1):
        print(f"Particle {i}: X={item.x:.2f}, Y={item.y:.2f}")
    print("Saved boxes to positions_all.box")


if __name__ == "__main__":
    pick_particles_with_input()
