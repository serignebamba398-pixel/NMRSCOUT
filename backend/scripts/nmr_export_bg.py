#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
nmr_export_bg.py — Export PNG + JSON from Bruker pdata (server/frontend compatible)
- 1D (1H, 13C): line plot
- 2D (COSY, HSQC): contour/image

Required orientation:
    * F2 (horizontal) decreases left→right (ppm max on the left) — like TopSpin
    * F1 (vertical) decreases top→bottom (e.g., -2 at top, 10 at bottom)
"""

import argparse
import json
import os
import sys
import numpy as np
import nmrglue as ng
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------- utils ----------------
def _float(v, default=None):
    try: return float(v)
    except Exception: return default

def _int(v, default=None):
    try: return int(v)
    except Exception: return default

def _ppm_axis_from_proc(proc: dict, size: int):
    """Descending ppm axis from OFFSET, SW_p (Hz), SF (MHz)."""
    off   = float(proc["OFFSET"])
    sw_hz = float(proc["SW_p"])
    sf    = float(proc["SF"])  # MHz
    dppm  = (sw_hz / sf) / float(size)   # ppm / point
    idx   = np.arange(size, dtype=float)
    return off - idx * dppm               # descending

def _ensure_orientation_2d(data, si_f1, si_f2):
    """Returns an (F1, F2) matrix: transpose if needed."""
    ny, nx = data.shape
    if (ny, nx) == (si_f1, si_f2):
        return data
    if (ny, nx) == (si_f2, si_f1):
        return data.T
    d1 = abs(ny - si_f1) + abs(nx - si_f2)
    d2 = abs(nx - si_f1) + abs(ny - si_f2)
    return data if d1 <= d2 else data.T

def _level_array(Z, level_start=0.10, n_levels=20):
    pos = float(np.max(np.maximum(Z, 0.0)))
    if pos <= 0:
        pos = float(np.max(np.abs(Z))) or 1.0
    start = max(1e-9, float(level_start)) * pos
    mult = 1.20
    return [start * (mult ** i) for i in range(int(n_levels))]

def _axes_bbox_px(fig, ax, dpi):
    """
    Returns the axes bbox in pixels within the final image.
    - fig_w_px, fig_h_px: PNG dimensions
    - bbox_px: {x,y,w,h} with (0,0) at the top-left of the image
    """
    fig.canvas.draw()  # compute layout
    bb_in = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig_w_px = int(round(fig.get_figwidth()  * dpi))
    fig_h_px = int(round(fig.get_figheight() * dpi))
    x0 = int(round(bb_in.x0 * dpi))
    y0 = int(round(bb_in.y0 * dpi))   # from bottom in Matplotlib coords
    w  = int(round(bb_in.width  * dpi))
    h  = int(round(bb_in.height * dpi))
    # convert y-bottom -> y-top (image origin is top-left)
    y  = fig_h_px - (y0 + h)
    return fig_w_px, fig_h_px, {"x": x0, "y": y, "w": w, "h": h}


# -------------- pdata reading --------------
def load_1d(pdata_dir):
    dic, data = ng.bruker.read_pdata(pdata_dir)
    if data.ndim != 1:
        data = np.ravel(data)
    procs = dic.get("procs", {}) or dic.get("proc", {})
    return data.astype(np.float64, copy=False), procs

def load_2d(pdata_dir):
    dic, data = ng.bruker.read_pdata(pdata_dir)
    procs  = dic.get("procs", {})  or dic.get("proc", {})
    proc2s = dic.get("proc2s", {}) or dic.get("proc2", {})
    si_f2 = _int(procs.get("SI"),  data.shape[-1])
    si_f1 = _int(proc2s.get("SI"), data.shape[0])
    Z = _ensure_orientation_2d(np.array(data, dtype=np.float64, copy=False), si_f1, si_f2)
    return Z, procs, proc2s


# ---------------- plots ----------------
def plot_1d(pdata_dir, out_png, meta_out,
            ppm_min=None, ppm_max=None, width_px=2000, height_px=360, dpi=200,
            transparent=False, line_color="black", line_width=0.9):
    y, procs = load_1d(pdata_dir)
    N   = y.shape[0]
    xppm = _ppm_axis_from_proc(procs, N)         # descending

    fig_w = float(width_px) / float(dpi)
    fig_h = float(height_px) / float(dpi)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)

    ax.plot(xppm, y, color=str(line_color), linewidth=float(line_width))

    # Determine the final scale used
    if ppm_min is not None and ppm_max is not None:
        lo, hi = sorted([float(ppm_min), float(ppm_max)])
        ax.set_xlim(hi, lo)
        # Use user parameters
        final_ppm_left, final_ppm_right = hi, lo
    else:
        ax.set_xlim(float(xppm.max()), float(xppm.min()))
        # Use data scale
        final_ppm_left, final_ppm_right = float(xppm.max()), float(xppm.min())

    ax.set_facecolor((0.96, 0.96, 0.96))
    ax.set_xlabel("ppm"); ax.set_ylabel("intensity")

    # --- exact bbox of the plotting area
    fig_w_px, fig_h_px, bbox_px = _axes_bbox_px(fig, ax, dpi)

    face = "none" if transparent else "white"
    # >> important: no tight layout, to keep stable geometry
    fig.savefig(out_png, dpi=dpi, facecolor=face, bbox_inches=None, pad_inches=0)
    plt.close(fig)

    # [FIX] store user parameters + final scale
    meta = {
        "kind": "1D", "unit": "ppm",
        "ppm_left":  final_ppm_left,      # final scale used (compat)
        "ppm_right": final_ppm_right,     # final scale used (compat)
        "ppm_min": float(ppm_min) if ppm_min is not None else None,  # user parameters (priority)
        "ppm_max": float(ppm_max) if ppm_max is not None else None,  # user parameters (priority)
        "shape": [int(N)],
        "fig_px": {"w": int(fig_w_px), "h": int(fig_h_px), "dpi": int(dpi)},
        "plot_bbox_px": bbox_px,
        "orientation": { "x_descending": True },
        "procs": {k: (float(v) if isinstance(v, (int,float,np.floating)) else v) for k,v in procs.items()},
    }
    with open(meta_out, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)


def plot_2d(pdata_dir, out_png, meta_out,
            f2_ppm_min=None, f2_ppm_max=None, f1_ppm_min=None, f1_ppm_max=None,
            width_px=1400, height_px=1100, dpi=200, mode="contour",
            transparent=True, n_levels=20, level_start=0.10, cmap="gray"):
    Z, procs, proc2s = load_2d(pdata_dir)

    Ny, Nx = Z.shape               # (F1, F2)
    f2 = _ppm_axis_from_proc(procs,  Nx)   # horizontal, descending
    f1 = _ppm_axis_from_proc(proc2s, Ny)   # vertical, descending

    fig_w = float(width_px) / float(dpi)
    fig_h = float(height_px) / float(dpi)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)

    if str(mode).lower() == "image":
        ax.imshow(
            Z, origin="lower",
            extent=[float(f2.min()), float(f2.max()), float(f1.min()), float(f1.max())],
            cmap=str(cmap), aspect="auto", interpolation="nearest",
        )
    else:
        levels = _level_array(Z, level_start=float(level_start), n_levels=int(n_levels))
        X, Y = np.meshgrid(f2, f1)
        ax.contour(X, Y, Z, levels=levels, colors="black", linewidths=0.6)

    # ---- windows / orientation ----
    # F2 like TopSpin: left=max, right=min
    if f2_ppm_min is not None and f2_ppm_max is not None:
        lo, hi = sorted([float(f2_ppm_min), float(f2_ppm_max)])
        ax.set_xlim(hi, lo)
        f2_left, f2_right = hi, lo
    else:
        ax.set_xlim(float(f2.max()), float(f2.min()))
        f2_left, f2_right = float(f2.max()), float(f2.min())

    # F1 decreases top→bottom ⇒ ylim(top=min, bottom=max) and Y inversion
    if f1_ppm_min is not None and f1_ppm_max is not None:
        lo, hi = sorted([float(f1_ppm_min), float(f1_ppm_max)])
        ax.set_ylim(hi, lo)                # top = lo, bottom = hi
        f1_top, f1_bottom = lo, hi
    else:
        ax.set_ylim(float(f1.min()), float(f1.max()))
        ax.invert_yaxis()
        f1_top, f1_bottom = float(f1.min()), float(f1.max())

    ax.set_facecolor((0.96, 0.96, 0.96))
    ax.set_xlabel("F2 [ppm]"); ax.set_ylabel("F1 [ppm]")

    # --- exact bbox of the plotting area
    fig_w_px, fig_h_px, bbox_px = _axes_bbox_px(fig, ax, dpi)

    face = "none" if transparent else "white"
    fig.savefig(out_png, dpi=dpi, facecolor=face, bbox_inches=None, pad_inches=0)
    plt.close(fig)

    # FIX: store user parameters + final scale
    meta = {
        "kind": "2D", "unit": "ppm", "shape": [int(Ny), int(Nx)],
        "mode": str(mode).lower(),
        "f2_left":  f2_left,  "f2_right": f2_right,       # final scale used (compat)
        "f1_top":   f1_top,   "f1_bottom": f1_bottom,     # final scale used (compat)
        # User parameters (priority for frontend)
        "f2_ppm_min": float(f2_ppm_min) if f2_ppm_min is not None else None,
        "f2_ppm_max": float(f2_ppm_max) if f2_ppm_max is not None else None,
        "f1_ppm_min": float(f1_ppm_min) if f1_ppm_min is not None else None,
        "f1_ppm_max": float(f1_ppm_max) if f1_ppm_max is not None else None,
        "fig_px": {"w": int(fig_w_px), "h": int(fig_h_px), "dpi": int(dpi)},
        "plot_bbox_px": bbox_px,
        "orientation": { "f2_descending": True, "f1_top_is_min": True },
        "procs":  {k: (float(v) if isinstance(v, (int,float,np.floating)) else v) for k,v in procs.items()},
        "proc2s": {k: (float(v) if isinstance(v, (int,float,np.floating)) else v) for k,v in proc2s.items()},
    }
    with open(meta_out, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)


# ---------------- CLI ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", required=True, choices=["1H","13C","COSY","HSQC"])
    ap.add_argument("--pdata_dir", required=True)
    ap.add_argument("--out_png", required=True)
    ap.add_argument("--out_json", required=True)
    # 1D
    ap.add_argument("--ppm_min", type=float)
    ap.add_argument("--ppm_max", type=float)
    ap.add_argument("--line_color", default="black")
    ap.add_argument("--line_width", type=float, default=0.9)
    # 2D
    ap.add_argument("--f2_ppm_min", type=float)
    ap.add_argument("--f2_ppm_max", type=float)
    ap.add_argument("--f1_ppm_min", type=float)
    ap.add_argument("--f1_ppm_max", type=float)
    ap.add_argument("--mode", default="contour", choices=["contour","image"])
    ap.add_argument("--n_levels", type=int, default=20)
    ap.add_argument("--level_start", type=float, default=0.10)
    ap.add_argument("--cmap", default="gray")
    # commun
    ap.add_argument("--width_px", type=int)
    ap.add_argument("--height_px", type=int)
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--transparent", type=int, default=0)
    args = ap.parse_args()

    pdata_dir = os.path.abspath(args.pdata_dir)
    if not os.path.isdir(pdata_dir):
        print(f"[export_bg] ERROR: pdata_dir not found: {pdata_dir}", file=sys.stderr)
        sys.exit(2)

    try:
        if args.method in ("1H","13C"):
            w = args.width_px or 2000
            h = args.height_px or 360
            plot_1d(
                pdata_dir, args.out_png, args.out_json,
                ppm_min=args.ppm_min, ppm_max=args.ppm_max,
                width_px=w, height_px=h, dpi=args.dpi,
                transparent=bool(args.transparent),
                line_color=args.line_color, line_width=args.line_width,
            )
        else:
            w = args.width_px or 1400
            h = args.height_px or 1100
            plot_2d(
                pdata_dir, args.out_png, args.out_json,
                f2_ppm_min=args.f2_ppm_min, f2_ppm_max=args.f2_ppm_max,
                f1_ppm_min=args.f1_ppm_min, f1_ppm_max=args.f1_ppm_max,
                width_px=w, height_px=h, dpi=args.dpi,
                mode=args.mode, transparent=bool(args.transparent),
                n_levels=args.n_levels, level_start=args.level_start, cmap=args.cmap,
            )
        print(f"[export_bg] OK -> {args.out_png}")
        sys.exit(0)
    except KeyError as e:
        print(f"[export_bg] ERROR: missing key {e!s}", file=sys.stderr); sys.exit(3)
    except Exception as e:
        print(f"[export_bg] ERROR: {str(e)}", file=sys.stderr); sys.exit(4)

if __name__ == "__main__":
    main()


