#!/usr/bin/env python
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Bruker 1D processed-spectrum extraction utilities."""

import argparse
import json
import sys
import os
import numpy as np

def extract_bruker_1d(pdata_dir):
    """Extract processed 1D Bruker spectrum data."""
    try:
        import nmrglue as ng
    except ImportError:
        return {
            "status": "error",
            "error": "nmrglue not installed"
        }

    try:
        dic, data = ng.bruker.read_pdata(pdata_dir)

        procs_dic = dic.get('procs', {})

        sf = float(procs_dic.get('SF', 700.13))

        sw_ppm = float(procs_dic.get('SW_p', 12.0)) / sf

        offset_ppm = float(procs_dic.get('OFFSET', 0.0))

        si = int(procs_dic.get('SI', len(data)))

        uc = ng.bruker.guess_udic(dic, data)['ndim'] - 1
        uc_ppm = ng.fileiobase.uc_from_udic(ng.bruker.guess_udic(dic, data), uc)
        ppm_scale = uc_ppm.ppm_scale()

        if data.ndim > 1:
            data = data.flatten()
        if np.iscomplexobj(data):
            data = np.real(data)

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

    if not os.path.isdir(args.pdata_dir):
        result = {
            "status": "error",
            "error": f"Directory not found: {args.pdata_dir}"
        }
        print(json.dumps(result), file=sys.stdout)
        sys.exit(1)

    result = extract_bruker_1d(args.pdata_dir)

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
