#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
extract_bruker_2d.py

Extract processed 2D Bruker NMR spectrum data (contours or raw) from pdata folder.
Returns contour polylines for direct canvas rendering (no PNG needed).

Usage:
    python extract_bruker_2d.py --pdata_dir "path/to/pdata/1" --output json

Output (JSON on stdout):
    {
      "status": "ok",
      "contours": [
        {"level": 0.1, "paths": [[[x1,y1],[x2,y2],...], ...]},
        ...
      ],
      "f2_range": [min, max],
      "f1_range": [min, max],
      "shape": [ny, nx],
      "sf_f2": 700.13,
      "sf_f1": 700.13
    }
"""

import argparse
import json
import sys
import os
import numpy as np


def _float(v, default=None):
    try:
        return float(v)
    except Exception:
        return default


def _ppm_axis_from_proc(proc: dict, size: int):
    """Descending ppm axis from OFFSET, SW_p (Hz), SF (MHz)."""
    off = float(proc.get("OFFSET", 0))
    sw_hz = float(proc.get("SW_p", 10000))
    sf = float(proc.get("SF", 400))  # MHz
    dppm = (sw_hz / sf) / float(size)   # ppm / point
    idx = np.arange(size, dtype=float)
    return off - idx * dppm  # descending


def _level_array(Z, level_start=0.02, n_levels=30, mult=1.4):
    """Generate contour levels using TopSpin-style geometric progression.
    
    TopSpin default: factor = 1.4 (geometric scale)
    Each contour level = start * mult^i
    
    Args:
        Z: 2D data array
        level_start: Starting level as fraction of max (default 0.02 = 2%)
        n_levels: Number of contour levels (default 30 for fine control)
        mult: Multiplication factor between levels (default 1.4 = TopSpin standard)
    """
    pos = float(np.max(np.maximum(Z, 0.0)))
    if pos <= 0:
        pos = float(np.max(np.abs(Z))) or 1.0
    start = max(1e-9, float(level_start)) * pos
    return [start * (mult ** i) for i in range(int(n_levels))]


def extract_contours(Z, f2_ppm, f1_ppm, level_start=0.15, n_levels=10):
    """
    Extract contour lines from 2D data using matplotlib.
    Returns list of contour objects with paths.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    
    levels = _level_array(Z, level_start=level_start, n_levels=n_levels)
    
    # Create figure just to extract contours
    fig, ax = plt.subplots()
    X, Y = np.meshgrid(f2_ppm, f1_ppm)
    cs = ax.contour(X, Y, Z, levels=levels)
    plt.close(fig)
    
    # Extract contour paths using allsegs (cleaner than collections)
    contours = []
    for i, level in enumerate(cs.levels):
        paths = []
        if i < len(cs.allsegs):
            for seg in cs.allsegs[i]:
                if len(seg) < 3:
                    continue
                # Convert to list of [f2, f1] pairs
                path = [[float(pt[0]), float(pt[1])] for pt in seg]
                # Simplify if too many points
                if len(path) > 150:
                    step = max(1, len(path) // 150)
                    path = path[::step]
                    # Always include last point to close loops properly
                    if path[-1] != [[float(seg[-1][0]), float(seg[-1][1])]]:
                        path.append([float(seg[-1][0]), float(seg[-1][1])])
                paths.append(path)
        
        if paths:
            contours.append({
                "level": float(level),
                "paths": paths
            })
    
    return contours


def extract_bruker_2d(pdata_dir, level_start=0.02, n_levels=30, decimate=4):
    """
    Extract processed 2D Bruker spectrum and compute contours.
    
    Args:
        pdata_dir: Path to pdata directory containing 2rr, procs, proc2s files
        level_start: Starting level as fraction of max (default 0.08 = 8%)
        n_levels: Number of contour levels
        decimate: Decimation factor to reduce data size (default 4)
    
    Returns:
        dict with contours, ranges, metadata
    """
    try:
        import nmrglue as ng
    except ImportError:
        return {
            "status": "error",
            "error": "nmrglue not installed"
        }

    try:
        # Read processed Bruker 2D data
        dic, data = ng.bruker.read_pdata(pdata_dir)
        
        procs = dic.get("procs", {}) or dic.get("proc", {})
        proc2s = dic.get("proc2s", {}) or dic.get("proc2", {})
        
        # Get shape
        if data.ndim != 2:
            return {"status": "error", "error": "Not a 2D spectrum"}
        
        ny, nx = data.shape
        si_f2 = int(procs.get("SI", nx))
        si_f1 = int(proc2s.get("SI", ny))
        
        # Ensure correct orientation
        if (ny, nx) != (si_f1, si_f2):
            if (ny, nx) == (si_f2, si_f1):
                data = data.T
                ny, nx = data.shape
        
        Z = np.array(data, dtype=np.float64)
        
        # Create PPM axes
        f2_ppm = _ppm_axis_from_proc(procs, nx)   # horizontal (F2)
        f1_ppm = _ppm_axis_from_proc(proc2s, ny)  # vertical (F1)
        
        # Decimate data for smaller response
        if decimate > 1:
            Z = Z[::decimate, ::decimate]
            f2_ppm = f2_ppm[::decimate]
            f1_ppm = f1_ppm[::decimate]
        
        # Extract contours
        contours = extract_contours(Z, f2_ppm, f1_ppm, level_start, n_levels)
        
        # Get spectral frequencies
        sf_f2 = _float(procs.get("SF"), 400.0)
        sf_f1 = _float(proc2s.get("SF"), 400.0)
        
        return {
            "status": "ok",
            "contours": contours,
            "f2_range": [float(f2_ppm.min()), float(f2_ppm.max())],
            "f1_range": [float(f1_ppm.min()), float(f1_ppm.max())],
            "shape": [int(ny), int(nx)],
            "sf_f2": sf_f2,
            "sf_f1": sf_f1,
            "n_contours": len(contours)
        }
        
    except FileNotFoundError as e:
        return {
            "status": "error",
            "error": f"Bruker files not found in {pdata_dir}: {str(e)}"
        }
    except Exception as e:
        return {
            "status": "error",
            "error": f"Failed to read Bruker 2D spectrum: {str(e)}"
        }


def main():
    parser = argparse.ArgumentParser(
        description="Extract processed 2D Bruker NMR spectrum contours"
    )
    parser.add_argument(
        "--pdata_dir",
        required=True,
        help="Path to pdata directory (containing 2rr, procs, proc2s files)"
    )
    parser.add_argument(
        "--level_start",
        type=float,
        default=0.08,
        help="Starting contour level as fraction of max (default 0.08)"
    )
    parser.add_argument(
        "--n_levels",
        type=int,
        default=15,
        help="Number of contour levels (default 15)"
    )
    parser.add_argument(
        "--decimate",
        type=int,
        default=4,
        help="Decimation factor to reduce data size (default 4)"
    )
    parser.add_argument(
        "--output",
        default="json",
        choices=["json"],
        help="Output format (json)"
    )

    args = parser.parse_args()

    # Validate path exists
    if not os.path.isdir(args.pdata_dir):
        result = {
            "status": "error",
            "error": f"Directory not found: {args.pdata_dir}"
        }
        print(json.dumps(result), file=sys.stdout)
        sys.exit(1)

    # Extract 2D spectrum contours
    result = extract_bruker_2d(
        args.pdata_dir,
        level_start=args.level_start,
        n_levels=args.n_levels,
        decimate=args.decimate
    )

    # Output
    print(json.dumps(result), file=sys.stdout)
    sys.exit(0 if result["status"] == "ok" else 1)


if __name__ == "__main__":
    main()
