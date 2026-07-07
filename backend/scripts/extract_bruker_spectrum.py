#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
extract_bruker_spectrum.py

Extract processed 1D Bruker NMR spectrum data (ppm, intensity) from pdata folder.

Usage:
    python extract_bruker_spectrum.py --pdata_dir "path/to/pdata/1" --output json

Output (JSON on stdout):
    {
      "ppm": [10.0, 9.99, ...],
      "intensity": [0.01, 0.02, ...],
      "n_points": 65536,
      "sw_ppm": 12.0,
      "sf": 700.13,
      "status": "ok"
    }
"""

import argparse
import json
import sys
import os
import numpy as np

def extract_bruker_1d(pdata_dir):
    """
    Extract processed 1D Bruker spectrum using nmrglue.

    Args:
        pdata_dir: Path to pdata/1 directory containing 1r, 1i, procs files

    Returns:
        dict with ppm, intensity arrays and metadata
    """
    try:
        import nmrglue as ng
    except ImportError:
        return {
            "status": "error",
            "error": "nmrglue not installed"
        }

    try:
        # Read processed Bruker data
        dic, data = ng.bruker.read_pdata(pdata_dir)

        # Get spectral parameters
        procs_dic = dic.get('procs', {})

        # Spectral frequency (SF = spectrometer frequency in MHz)
        sf = float(procs_dic.get('SF', 700.13))

        # Spectral width in ppm
        sw_ppm = float(procs_dic.get('SW_p', 12.0)) / sf

        # Offset (O1P - chemical shift of carrier in ppm)
        offset_ppm = float(procs_dic.get('OFFSET', 0.0))

        # Number of points
        si = int(procs_dic.get('SI', len(data)))

        # Create ppm scale (reversed for NMR convention)
        uc = ng.bruker.guess_udic(dic, data)['ndim'] - 1
        uc_ppm = ng.fileiobase.uc_from_udic(ng.bruker.guess_udic(dic, data), uc)
        ppm_scale = uc_ppm.ppm_scale()

        # Ensure data is real and 1D
        if data.ndim > 1:
            data = data.flatten()
        if np.iscomplexobj(data):
            data = np.real(data)

        # Normalize intensity to [0, 1]
        intensity = np.array(data, dtype=float)
        max_int = np.max(np.abs(intensity))
        if max_int > 0:
            intensity = intensity / max_int

        return {
            "status": "ok",
            "ppm": ppm_scale.tolist(),
            "intensity": intensity.tolist(),
            "n_points": len(ppm_scale),
            "sw_ppm": sw_ppm,
            "sf": sf,
            "offset_ppm": offset_ppm
        }

    except FileNotFoundError as e:
        return {
            "status": "error",
            "error": f"Bruker files not found in {pdata_dir}: {str(e)}"
        }
    except Exception as e:
        return {
            "status": "error",
            "error": f"Failed to read Bruker spectrum: {str(e)}",
            "traceback": str(e.__traceback__)
        }


def main():
    parser = argparse.ArgumentParser(
        description="Extract processed 1D Bruker NMR spectrum data"
    )
    parser.add_argument(
        "--pdata_dir",
        required=True,
        help="Path to pdata/1 directory (containing 1r, 1i, procs files)"
    )
    parser.add_argument(
        "--output",
        default="json",
        choices=["json", "csv"],
        help="Output format (json or csv)"
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

    # Extract spectrum
    result = extract_bruker_1d(args.pdata_dir)

    # Output
    if args.output == "json":
        print(json.dumps(result), file=sys.stdout)
    elif args.output == "csv":
        if result["status"] == "ok":
            print("ppm,intensity")
            for ppm, intensity in zip(result["ppm"], result["intensity"]):
                print(f"{ppm},{intensity}")
        else:
            print(f"ERROR: {result['error']}", file=sys.stderr)
            sys.exit(1)

    sys.exit(0 if result["status"] == "ok" else 1)


if __name__ == "__main__":
    main()
