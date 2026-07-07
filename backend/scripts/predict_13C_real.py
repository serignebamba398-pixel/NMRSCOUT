#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
predict_13C_real.py (server version, refactored)
------------------------------------------------
Extracts 13C signals from SDF, with options:
    - non-quaternary carbons (default)
    - all carbons

Server compatibility:
    - --sdf_path      : path to a .sdf file or a folder of .sdf files
    - --output_dir    : output folder
    - --task_id       : task identifier (CSV naming)
    - --recursive     : recursive traversal if --sdf_path is a folder
    - --carbon_mode   : {nonquaternary, all}

Output:
    <output_dir>/prediction_13C_<task_id>.csv

Columns (0-based indices, RDKit-consistent):
    Filename | C_index | Nb_H | H_indices | Delta_ppm | ECFP4
"""

import os
import sys
import json
import argparse
import traceback
from typing import List

import pandas as pd
from rdkit import RDLogger

# Robust utils (SDF reading / shift-vector detection / ECFP4)
from ml_sdf_utils import (
    read_molecule_keepH,
    get_atom_shift_vector,
    ecfp4_bitstring,
)

# Reduce RDKit noise
RDLogger.DisableLog('rdApp.*')


# ------------------------- JSON logging utilities (SSE) -------------------------
def jlog(payload: dict):
    """Writes a JSON object to stdout (for server-side SSE)."""
    print(json.dumps(payload, ensure_ascii=False), flush=True)

def log(msg: str):
    jlog({"type": "log", "message": msg})

def progress(pct: float, message: str = None):
    p = max(0, min(100, int(pct)))  # [NEW] normalize 0-100
    obj = {"type": "progress", "current": p, "total": 100, "percent": p}  # [NEW]
    if message:
        obj["message"] = message
    jlog(obj)

def error(msg: str, code: int = 1):
    jlog({"type": "error", "message": msg})
    jlog({"type": "completed", "success": False, "error": msg})
    sys.exit(code)

def completed(out_csv: str, n_rows: int):
    jlog({
        "type": "completed",
        "success": True,
        "message": f"13C extraction completed ({n_rows} rows)",
        "output": out_csv
    })


# ----------------------------- Input files ------------------------------
def list_sdf_files(sdf_path: str, recursive: bool = False) -> List[str]:
    """Returns the list of SDF files to process (single file or all .sdf files in a folder)."""
    if not os.path.exists(sdf_path):
        return []

    if os.path.isfile(sdf_path):
        return [sdf_path]

    files = []
    if recursive:
        for root, _, fs in os.walk(sdf_path):
            for f in fs:
                if f.lower().endswith('.sdf'):
                    files.append(os.path.join(root, f))
    else:
        for f in os.listdir(sdf_path):
            if f.lower().endswith('.sdf'):
                files.append(os.path.join(sdf_path, f))

    return sorted(files)


def ensure_dir_exists(d: str):
    if d and not os.path.exists(d):
        os.makedirs(d, exist_ok=True)


# -------------------------- 13C core functions --------------------------------
def get_protonated_carbons(mol) -> List[int]:
    """
    Returns indices of carbons bearing at least one H (non-quaternary).
    0-based indexing consistent with RDKit and shift vectors.
    """
    prot = []
    for a in mol.GetAtoms():
        if a.GetSymbol() != 'C':
            continue
        h_count = 0
        for nb in a.GetNeighbors():
            if nb.GetSymbol() == 'H':
                h_count += 1
        if h_count > 0:
            prot.append(a.GetIdx())
    return prot


def get_all_carbons(mol) -> List[int]:
    """Returns 0-based indices of all carbons in the molecule."""
    return [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'C']


def carbon_h_neighbors(mol, c_idx: int) -> List[int]:
    """0-based indices of H atoms directly bonded to carbon c_idx."""
    hs = []
    c_atom = mol.GetAtomWithIdx(c_idx)
    for nb in c_atom.GetNeighbors():
        if nb.GetSymbol() == 'H':
            hs.append(nb.GetIdx())
    return hs


def dataframe_for_sdf(sdf_path: str, carbon_mode: str = "nonquaternary") -> pd.DataFrame:
    """
        Builds a DataFrame (0-based):
      Filename | C_index | Nb_H | H_indices | Delta_ppm | ECFP4
        - Delta_ppm: taken at the carbon atom index in the detected shift vector.
        - H_indices: "i,j,k" (0-based) string for simple CSV compatibility.
    """
    try:
        progress(10, f"Reading {os.path.basename(sdf_path)}...")  # [NEW]
        mol = read_molecule_keepH(sdf_path)
        if mol is None:
            log(f"[ERROR] Unreadable SDF: {os.path.basename(sdf_path)}")
            return pd.DataFrame()

        # Shift vector robustly detected (not strictly '(1)')
        shifts = get_atom_shift_vector(sdf_path)  # 0-based aligned
        if not shifts:
            log(f"[WARN] No shift vector detected: {os.path.basename(sdf_path)} (Delta_ppm will be NA)")

        # ECFP4 (bitstring 1024)
        fp = ecfp4_bitstring(mol)

        # Select carbons based on mode
        if carbon_mode == "all":
            c_indices = get_all_carbons(mol)
        else:
            c_indices = get_protonated_carbons(mol)

        progress(40, f"Extracting carbons ({len(c_indices)} found)...")  # [NEW]

        rows = []
        for c_idx in c_indices:
            # 13C ppm at the carbon index, if available
            delta = None
            if 0 <= c_idx < len(shifts):
                try:
                    delta = round(float(shifts[c_idx]), 3)
                except Exception:
                    delta = None

            h_idx_list = carbon_h_neighbors(mol, c_idx)
            rows.append({
                "Filename":  os.path.basename(sdf_path),
                "C_index":   c_idx,  # 0-based
                "Nb_H":      len(h_idx_list),
                "H_indices": "|".join(map(str, h_idx_list)) if h_idx_list else "",
                "Delta_ppm": delta,
                "ECFP4":     fp,
            })

        progress(70, "Building DataFrame...")  # [NEW]
        return pd.DataFrame(rows)

    except Exception as e:
        log(f"[ERROR] Error on {os.path.basename(sdf_path)}: {e}")
        traceback.print_exc(file=sys.stdout)
        return pd.DataFrame()


# ----------------------------------- Main -------------------------------------
def build_argparser():
    p = argparse.ArgumentParser(
        description="13C extraction + ECFP4 from SDF (file or folder), 0-based indices."
    )
    p.add_argument('--sdf_path', required=True, help='SDF file or folder of SDF files')
    p.add_argument('--output_dir', required=True, help='Output folder')
    p.add_argument('--task_id', default='task', help='Task ID to name the CSV')
    p.add_argument('--recursive', action='store_true', help='Recursive traversal if --sdf_path is a folder')
    p.add_argument('--carbon_mode',
                   choices=['nonquaternary', 'all'],
                   default='nonquaternary',
                   help='Carbon selection: nonquaternary (default) or all')
    return p


def main():
    args = build_argparser().parse_args()

    try:
        ensure_dir_exists(args.output_dir)
        out_csv = os.path.join(args.output_dir, f"prediction_13C_{args.task_id}.csv")

        files = list_sdf_files(args.sdf_path, recursive=args.recursive)
        if not files:
            error(f"No SDF found in {args.sdf_path}", code=2)

        log(f"SDF files detected: {len(files)}")
        log(f"Carbon mode: {args.carbon_mode}")
        all_frames: List[pd.DataFrame] = []

        total = len(files)
        for i, sdf in enumerate(files, 1):
            pct = int(5 + ((i - 1) / max(total, 1)) * 85)  # [NEW] 5% -> 90%
            progress(pct, f"Processing {i}/{total}: {os.path.basename(sdf)}")  # [NEW]
            df = dataframe_for_sdf(sdf, carbon_mode=args.carbon_mode)
            if not df.empty:
                all_frames.append(df)

        if not all_frames:
            error("No 13C data extracted: missing shift vector or invalid SDF(s)?", code=3)

        result = pd.concat(all_frames, ignore_index=True)
        progress(95, "Writing CSV file...")  # [NEW]
        result.to_csv(out_csv, index=False, encoding='utf-8')

        progress(100, "13C extraction completed")  # [NEW]
        completed(out_csv, len(result))

    except SystemExit:
        raise
    except Exception as e:
        jlog({"type": "error", "message": f"Exception: {e}"})
        traceback.print_exc(file=sys.stdout)
        error("13C extraction failed (see trace above)", code=1)


if __name__ == '__main__':
    main()