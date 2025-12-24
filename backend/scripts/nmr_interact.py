#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
nmr_interact.py — interactive Bruker pdata viewer (1D/2D) with zoom + mouse wheel
NMR orientation:
    - F2 (horizontal) decreases left → right
    - F1 (vertical) decreases top → bottom
Controls:
    - Left mouse: draw a rectangle to zoom (1D: X zoom, 2D: X/Y zoom)
    - Mouse wheel:
            * 1D: vertical zoom (helps reveal low intensities)
            * 2D: changes the contour threshold (low intensities become visible)
    - Middle click: reset view
    - Right mouse: zoom out
    - R: reset     |  S: save PNG
    - 2D: M: contour/image,  + / - : number of levels
"""

import argparse
import os
import numpy as np
import nmrglue as ng
import matplotlib
# Ensure an interactive backend (Qt). If another one is already active, keep it.
try:
    if matplotlib.get_backend().lower() in ("agg", "template"):
        matplotlib.use("QtAgg", force=True)
except Exception:
    pass
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector, SpanSelector
from matplotlib.colors import LogNorm, SymLogNorm

# ---------------- utils NMR ----------------
def _float(v, default=None):
    try:
        return float(v)
    except Exception:
        return default

def _ppm_axis_from_proc(proc: dict, size: int):
    """Descending ppm axis from OFFSET, SW_p (Hz), SF (MHz)."""
    off   = float(proc["OFFSET"])
    sw_hz = float(proc["SW_p"])
    sf    = float(proc["SF"])  # MHz
    dppm  = (sw_hz / sf) / float(size)   # ppm/point
    idx   = np.arange(size, dtype=float)
    return off - idx * dppm               # descending

def load_1d(pdata_dir):
    dic, data = ng.bruker.read_pdata(pdata_dir)
    if data.ndim != 1:
        data = np.ravel(data)
    procs = dic.get("procs", {}) or dic.get("proc", {})
    y = np.asarray(data, dtype=np.float64)
    xppm = _ppm_axis_from_proc(procs, y.size)
    return xppm, y

def load_2d(pdata_dir):
    dic, data = ng.bruker.read_pdata(pdata_dir)
    procs  = dic.get("procs", {})  or dic.get("proc", {})
    proc2s = dic.get("proc2s", {}) or dic.get("proc2", {})
    Z = np.asarray(data, dtype=np.float64, order="C")
    # target orientation: (F1, F2)
    si_f2 = int(procs.get("SI", Z.shape[-1]))
    si_f1 = int(proc2s.get("SI", Z.shape[0]))
    if Z.shape != (si_f1, si_f2):
        if Z.shape == (si_f2, si_f1):
            Z = Z.T
        # otherwise, keep as-is
    f2 = _ppm_axis_from_proc(procs,  Z.shape[1])   # descending
    f1 = _ppm_axis_from_proc(proc2s, Z.shape[0])   # descending
    return f1, f2, Z

# ---------------- viewer 1D ----------------
class OneDViewer:
    def __init__(self, pdata_dir, title="1D"):
        self.title = title
        self.xppm, self.y = load_1d(pdata_dir)
        self.scale_y = 1.0

        self.fig, self.ax = plt.subplots(figsize=(10, 3.2))
        try:
            self.fig.canvas.manager.set_window_title(f"{title} — {pdata_dir}")
        except Exception:
            pass

        # initial plot
        self.line, = self.ax.plot(self.xppm, self.y, lw=0.9, color="black")
        self.ax.set_xlabel("ppm")
        self.ax.set_ylabel("intensity")
        # F2 descending left→right
        self.ax.set_xlim(float(self.xppm.max()), float(self.xppm.min()))
        # Y symmetric around 0 initially
        yabs = np.max(np.abs(self.y)) or 1.0
        self._y_full = (-yabs, yabs)
        self.ax.set_ylim(*self._y_full)
        self.ax.grid(True, alpha=0.25)
        self._x_full = tuple(self.ax.get_xlim())
        self.ax.set_navigate(True)

        # --- interactions ---
        # SpanSelector = zoom horizontal robuste pour 1D
        self.selector = SpanSelector(
            self.ax, self.on_span, "horizontal",
            useblit=True, minspan=1e-6,
            interactive=True, drag_from_anywhere=True
        )
        self.cid_scroll = self.fig.canvas.mpl_connect("scroll_event", self.on_scroll)
        self.cid_button = self.fig.canvas.mpl_connect("button_press_event", self.on_click)
        self.cid_key    = self.fig.canvas.mpl_connect("key_press_event", self.on_key)

        plt.tight_layout()

    # Span (X uniquement)
    def on_span(self, x_min, x_max):
        if x_min is None or x_max is None:
            return
        lo, hi = (x_min, x_max) if x_min < x_max else (x_max, x_min)
        if abs(hi - lo) < 1e-12:
            return
        # keep descending ppm orientation (max on the left)
        self.ax.set_xlim(hi, lo)
        self.fig.canvas.draw_idle()

    # mouse wheel: vertical zoom (helps reveal low intensities)
    def on_scroll(self, event):
        if event.inaxes != self.ax:
            return
        factor = 1.15 if event.button == "up" else (1/1.15)
        ymin, ymax = self.ax.get_ylim()
        cy = 0.0  # centre autour de 0
        new_ymin = cy + (ymin - cy) / factor
        new_ymax = cy + (ymax - cy) / factor
        # limite min
        eps = 1e-12
        if abs(new_ymax - new_ymin) < eps:
            return
        self.ax.set_ylim(new_ymin, new_ymax)
        self.fig.canvas.draw_idle()

    def on_click(self, event):
        # middle click => reset / right click => zoom out (x2)
        if event.button == 2:
            self.reset_view()
        elif event.button == 3:
            xlo, xhi = self.ax.get_xlim()
            span = (xlo - xhi)   # note: xlo>xhi (inverted axis)
            if abs(span) < 1e-12:
                return
            cx = xhi + span/2.0
            new_span = span * 2.0
            self.ax.set_xlim(cx + new_span/2.0, cx - new_span/2.0)
            self.fig.canvas.draw_idle()

    def on_key(self, event):
        if event.key in ("r", "R"):
            self.reset_view()
        elif event.key in ("s", "S"):
            self.save_png()

    def reset_view(self):
        self.ax.set_xlim(*self._x_full)
        self.ax.set_ylim(*self._y_full)
        self.fig.canvas.draw_idle()

    def save_png(self, path=None):
        path = path or f"{self.title.replace(' ','_')}.png"
        self.fig.savefig(path, dpi=200)
        print(f"[1D] PNG saved: {path}")

# ---------------- viewer 2D ----------------
class TwoDViewer:
    def __init__(self, pdata_dir, title="2D (COSY/HSQC)", mode="contour"):
        self.title = title
        self.f1, self.f2, self.Z = load_2d(pdata_dir)
        self.mode = str(mode).lower()  # "contour" | "image"
        self.n_levels = 20
        self.level_start = 0.10   # fraction du max positif
        self.mult = 1.20

        self.fig, self.ax = plt.subplots(figsize=(9.5, 7.8))
        try:
            self.fig.canvas.manager.set_window_title(f"{title} — {pdata_dir}")
        except Exception:
            pass

        # compute initial levels
        self.update_levels()
        # initial plot
        self.im = None
        self.cont = None
        self.draw()

        # orientation NMR
        self.apply_axes_orientation()

        # interactions
        self.selector = RectangleSelector(
            self.ax, self.on_select, button=[1], useblit=True,
            interactive=True, drag_from_anywhere=True,
            minspanx=0.0001, minspany=0.0001
        )
        self.cid_scroll = self.fig.canvas.mpl_connect("scroll_event", self.on_scroll)
        self.cid_button = self.fig.canvas.mpl_connect("button_press_event", self.on_click)
        self.cid_key    = self.fig.canvas.mpl_connect("key_press_event", self.on_key)

        self.ax.set_xlabel("F2 [ppm]")
        self.ax.set_ylabel("F1 [ppm]")
        self.ax.grid(True, alpha=0.2)
        plt.tight_layout()

    def apply_axes_orientation(self):
        # F2: left = max, right = min
        self.ax.set_xlim(float(self.f2.max()), float(self.f2.min()))
        # F1: top = min, bottom = max
        self.ax.set_ylim(float(self.f1.min()), float(self.f1.max()))
        self.ax.invert_yaxis()
        self._x_full = tuple(self.ax.get_xlim())
        self._y_full = tuple(self.ax.get_ylim())

    def update_levels(self):
        Zpos = np.maximum(self.Z, 0.0)
        pmax = float(np.max(Zpos))
        if pmax <= 0:
            pmax = float(np.max(np.abs(self.Z))) or 1.0
        start = max(1e-12, self.level_start) * pmax
        self.levels = [start * (self.mult ** i) for i in range(self.n_levels)]

    def draw(self):
        # efface
        for art in (self.im, self.cont):
            if art is not None:
                try:
                    art.remove()
                except Exception:
                    pass
        self.im = None
        self.cont = None

        if self.mode == "image":
            # image (useful for quickly exploring intensity dynamics)
            absZ = np.abs(self.Z)
            vmax = np.quantile(absZ, 0.998) or (absZ.max() or 1.0)
            vmin = max(vmax * 1e-4, 1e-12)
            # SymLogNorm: helps reveal low positive/negative intensities
            norm = SymLogNorm(linthresh=vmin, vmin=-vmax, vmax=vmax, base=10)
            self.im = self.ax.imshow(
                self.Z, origin="lower",
                extent=[float(self.f2.min()), float(self.f2.max()),
                        float(self.f1.min()), float(self.f1.max())],
                cmap="gray", aspect="auto", interpolation="nearest", norm=norm
            )
        else:
            # contours (positifs)
            X, Y = np.meshgrid(self.f2, self.f1)
            self.cont = self.ax.contour(
                X, Y, self.Z, levels=self.levels, colors="black", linewidths=0.6
            )

        self.fig.canvas.draw_idle()

    # zoom rectangle (2D)
    def on_select(self, eclick, erelease):
        if None in (eclick.xdata, eclick.ydata, erelease.xdata, erelease.ydata):
            return
        x1, x2 = eclick.xdata, erelease.xdata
        y1, y2 = eclick.ydata, erelease.ydata
        xlo, xhi = min(x1, x2), max(x1, x2)
        ylo, yhi = min(y1, y2), max(y1, y2)
        # keep orientation (F2 descending, F1 descending)
        self.ax.set_xlim(xhi, xlo)
        self.ax.set_ylim(yhi, ylo)  # invert_yaxis already active
        self.fig.canvas.draw_idle()

    # mouse wheel: adjust contour threshold (low intensities visible)
    def on_scroll(self, event):
        if event.inaxes != self.ax:
            return
        # scroll up => lower threshold => more weak signals
        factor = 0.85 if event.button == "up" else (1/0.85)
        self.level_start = max(1e-10, self.level_start * factor)
        self.update_levels()
        self.draw()

    def on_click(self, event):
        if event.button == 2:     # reset
            self.reset_view()
        elif event.button == 3:   # zoom out
            xlo, xhi = self.ax.get_xlim()
            ylo, yhi = self.ax.get_ylim()
            # spans (note inverted axes)
            xspan = (xlo - xhi)
            yspan = (yhi - ylo)
            cx = xhi + xspan/2.0
            cy = ylo + yspan/2.0
            self.ax.set_xlim(cx + xspan, cx - xspan)
            self.ax.set_ylim(cy + yspan, cy - yspan)
            self.fig.canvas.draw_idle()

    def on_key(self, event):
        if event.key in ("r", "R"):
            self.reset_view()
        elif event.key in ("s", "S"):
            self.save_png()
        elif event.key in ("m", "M"):
            self.mode = "image" if self.mode == "contour" else "contour"
            self.draw()
        elif event.key in ("+", "equal"):
            self.n_levels = min(80, self.n_levels + 2)
            self.update_levels(); self.draw()
        elif event.key in ("-",):
            self.n_levels = max(2, self.n_levels - 2)
            self.update_levels(); self.draw()

    def reset_view(self):
        self.ax.set_xlim(*self._x_full)
        self.ax.set_ylim(*self._y_full)
        self.fig.canvas.draw_idle()

    def save_png(self, path=None):
        path = path or f"{self.title.replace(' ','_')}.png"
        self.fig.savefig(path, dpi=200)
        print(f"[2D] PNG saved: {path}")

# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", required=True, choices=["1H", "13C", "COSY", "HSQC"])
    ap.add_argument("--pdata_dir", required=True, help="path to pdata (e.g., .../1/pdata/1)")
    ap.add_argument("--mode2d", default="contour", choices=["contour", "image"], help="initial 2D mode")
    args = ap.parse_args()

    if args.method in ("1H", "13C"):
        v = OneDViewer(args.pdata_dir, title=args.method)
    else:
        v = TwoDViewer(args.pdata_dir, title=args.method, mode=args.mode2d)

    print("Controls: left click=zoom | wheel=scale | middle click=reset | right click=zoom out | R=reset | S=save")
    if args.method in ("COSY", "HSQC"):
        print("2D: M=contour/image,  +/-=levels, wheel=contour threshold")
    plt.show()

if __name__ == "__main__":
    main()
