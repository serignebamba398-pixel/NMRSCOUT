#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
predict_COSY_real.py (server version, refactored)
-------------------------------------------------
Reconstructs COSY from SDF files containing an atom-wise shift vector,
robustly detected via ml_sdf_utils.get_atom_shift_vector().

Arguments:
    --sdf_path         : .sdf file or folder
    --output_dir       : output folder
    --task_id          : task identifier (CSV name)
    --recursive        : traverse folders recursively
    --min_bonds / --max_bonds : bond-distance window for COSY correlations
    --no_diagonal      : do not write the diagonal (H,H)

Output:
    <output_dir>/prediction_COSY_<task_id>.csv

Columns (0-based):
    Filename, H1, Delta1, H2, Delta2, Bonds, Group1, Group2, ECFP4, Method
"""

import os
import sys
import json
import argparse
import traceback
from typing import Dict, List, Tuple
from collections import defaultdict

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger, rdBase

# === robust utils: SDF reading, shift vector, ECFP4 ===
# (ensure ml_sdf_utils.py is on PYTHONPATH or in the same folder)
from ml_sdf_utils import (
    read_molecule_keepH,
    get_atom_shift_vector,
    ecfp4_bitstring,
)

# RDKit: reduce noise
RDLogger.DisableLog('rdApp.*')
rdBase.DisableLog('rdApp.error')
rdBase.DisableLog('rdApp.warning')


# ===================================================================
# 0) JSON LOGGING FOR SSE
# ===================================================================
def jlog(payload: dict):
    print(json.dumps(payload, ensure_ascii=False), flush=True)

def log(message: str):
    jlog({"type": "log", "message": message})

def progress(pct: float, message: str = None):
    p = max(0, min(100, int(pct)))  # [NEW] normalize 0-100
    obj = {"type": "progress", "current": p, "total": 100, "percent": p}  # [NEW]
    if message:
        obj["message"] = message
    jlog(obj)

def error(message: str, code: int = 1):
    jlog({"type": "error", "message": message})
    jlog({"type": "completed", "success": False, "error": message})
    sys.exit(code)

def completed(out_csv: str, n_rows: int):
    jlog({
        "type": "completed",
        "success": True,
        "message": f"COSY reconstructed ({n_rows} rows)",
        "output": out_csv
    })


# ===================================================================
# 1) FILE UTILITIES
# ===================================================================
def list_sdf_files(sdf_path: str, recursive: bool = False) -> List[str]:
    if not os.path.exists(sdf_path):
        return []
    if os.path.isfile(sdf_path):
        return [sdf_path]
    files: List[str] = []
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


# ===================================================================
# 2) PROTON GROUPS + DIASTEREOTOPY
# ===================================================================
def build_parent_to_hydrogens(mol: Chem.Mol) -> Dict[int, List[int]]:
    """
    parent_C_idx -> [H_idx,...] (excludes H on O/N/S).
    """
    parent2Hs: Dict[int, List[int]] = defaultdict(list)
    for h in mol.GetAtoms():
        if h.GetSymbol() != 'H':
            continue
        nbrs = h.GetNeighbors()
        if not nbrs:
            continue
        # Exclude OH/NH/SH
        if any(nb.GetSymbol() in ('O', 'N', 'S') for nb in nbrs):
            continue
        for nb in nbrs:
            if nb.GetSymbol() == 'C':
                parent2Hs[nb.GetIdx()].append(h.GetIdx())
                break
    return parent2Hs

def get_h_groups(mol: Chem.Mol,
                 shifts: List[float] = None) -> Tuple[Dict[int, int], Dict[int, List[int]]]:
    """
    Returns:
      - h2group : dict {H_idx -> group_id}
      - group2Hs : dict {group_id -> [H_idx,...]}
    Groups CH3 / CH2 / CH.

    The SDF is already corrected:
    - Equivalent H have identical shifts → grouped together
    - Diastereotopic H have different shifts → separate groups
    """
    parent2Hs = build_parent_to_hydrogens(mol)
    groups: List[List[int]] = []
    used = set()

    # CH3 - always grouped (equivalent by symmetry)
    for hs in parent2Hs.values():
        if len(hs) == 3:
            groups.append(list(hs))
            used.update(hs)

    # CH2 - split if shifts are different (diastereotopic, already corrected in SDF)
    for hs in parent2Hs.values():
        if len(hs) == 2 and not any(h in used for h in hs):
            h1, h2 = hs
            should_split = False

            if shifts is not None:
                if 0 <= h1 < len(shifts) and 0 <= h2 < len(shifts):
                    try:
                        shift1 = float(shifts[h1])
                        shift2 = float(shifts[h2])
                        # If shifts are different, they are diastereotopic
                        if shift1 != shift2:
                            should_split = True
                    except (ValueError, TypeError):
                        pass

            if should_split:
                groups.append([h1])
                groups.append([h2])
            else:
                groups.append(list(hs))
            used.update(hs)

    # CH
    for hs in parent2Hs.values():
        if len(hs) == 1 and hs[0] not in used:
            groups.append([hs[0]])
            used.add(hs[0])

    h2group: Dict[int, int] = {}
    group2Hs: Dict[int, List[int]] = {}
    for gi, grp in enumerate(groups):
        group2Hs[gi] = grp
        for h in grp:
            h2group[h] = gi
    return h2group, group2Hs


# ===================================================================
# 3) COSY PAIRS + SIGNALS
# ===================================================================
def get_cosy_pairs(mol: Chem.Mol,
                   h2group: Dict[int, int],
                   min_bonds: int = 2,
                   max_bonds: int = 3,
                   include_diagonal: bool = True) -> List[Tuple[int, int, int]]:
    """
    Pairs (h1, h2, bond_distance). Avoids intra-group pairs (except diagonal).
    """
    cosy = []
    h_list = sorted(h2group.keys())
    for i, h1 in enumerate(h_list):
        if include_diagonal:
            cosy.append((h1, h1, 0))
        for h2 in h_list[i+1:]:
            if h2group[h1] == h2group[h2]:
                continue
            path = Chem.GetShortestPath(mol, h1, h2)
            if not path:
                continue
            dist = len(path) - 1
            if min_bonds <= dist <= max_bonds:
                cosy.append((h1, h2, dist))
                cosy.append((h2, h1, dist))
    return cosy

def make_signal_map_per_H(shifts: List[float],
                          h2group: Dict[int, int],
                          group2Hs: Dict[int, List[int]]) -> Dict[int, float]:
    """
    H_idx -> delta (ppm); group average when possible, else individual value.
    """
    group_avg: Dict[int, float] = {}
    for gi, members in group2Hs.items():
        vals = [shifts[h] for h in members if 0 <= h < len(shifts)]
        group_avg[gi] = float(np.mean(vals)) if vals else float('nan')

    signal_map: Dict[int, float] = {}
    for h, gi in h2group.items():
        v = group_avg.get(gi, float('nan'))
        if not np.isnan(v):
            signal_map[h] = round(v, 3)
        elif 0 <= h < len(shifts):
            try:
                signal_map[h] = round(float(shifts[h]), 3)
            except Exception:
                signal_map[h] = float('nan')
        else:
            signal_map[h] = float('nan')
    return signal_map


# ===================================================================
# 4) CORE : UN SDF -> DATAFRAME
# ===================================================================
def dataframe_for_sdf(sdf_path: str,
                      min_bonds: int,
                      max_bonds: int,
                      include_diagonal: bool) -> pd.DataFrame:
    try:
        progress(10, f"Reading {os.path.basename(sdf_path)}...")  # [NEW]
        mol = read_molecule_keepH(sdf_path)
        if mol is None:
            log(f"[ERROR] Unreadable molecule: {os.path.basename(sdf_path)}")
            return pd.DataFrame()

        # robust shift vector (regardless of the property name)
        shifts = get_atom_shift_vector(sdf_path)  # 0-based
        if not shifts:
            log(f"[WARN] Shift vector not found: {os.path.basename(sdf_path)}")

        # H groups - passer les shifts pour auto-split des CH2 diastéréotopiques
        progress(30, "Grouping protons...")  # [NEW]
        h2group, group2Hs = get_h_groups(mol, shifts=shifts)
        if not h2group:
            log(f"[WARN] No eligible protons (excluding OH/NH/SH): {os.path.basename(sdf_path)}")
            return pd.DataFrame()

        # COSY pairs
        pairs = get_cosy_pairs(mol, h2group, min_bonds, max_bonds, include_diagonal)
        progress(50, f"Computing COSY correlations ({len(pairs)} pairs)...")  # [NEW]

        # 1H deltas (group average if >1)
        signal_map = make_signal_map_per_H(shifts, h2group, group2Hs)

        # ECFP4
        ecfp4 = ecfp4_bitstring(mol)

        rows = []
        for h1, h2, dist in pairs:
            d1 = signal_map.get(h1, float('nan'))
            d2 = signal_map.get(h2, float('nan'))
            rows.append({
                'Filename': os.path.basename(sdf_path),
                'H1': h1,
                'Delta1': d1 if np.isfinite(d1) else '',
                'H2': h2,
                'Delta2': d2 if np.isfinite(d2) else '',
                'Bonds': dist,
                'Group1': h2group.get(h1, -1),
                'Group2': h2group.get(h2, -1),
                'ECFP4': ecfp4,
                'Method': 'COSY_from_ML_SDF'
            })

        progress(70, "Building DataFrame...")  # [NEW]
        df = pd.DataFrame(rows)
        if not df.empty:
            df['__diag__'] = (df['H1'] == df['H2']).astype(int)
            df.sort_values(by=['__diag__', 'Bonds', 'Delta1'],
                           ascending=[False, True, True],
                           inplace=True)
            df.drop(columns='__diag__', inplace=True, errors='ignore')
        return df

    except Exception as e:
        log(f"[ERROR] Error on {os.path.basename(sdf_path)}: {e}")
        traceback.print_exc(file=sys.stdout)
        return pd.DataFrame()


# ===================================================================
# 5) MAIN / CLI
# ===================================================================
def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Reconstruct COSY (2–3 bonds) from SDF with detected shift vector (0-based indices)."
    )
    p.add_argument('--sdf_path', required=True, help='SDF file or folder of SDF files')
    p.add_argument('--output_dir', required=True, help='Output folder')
    p.add_argument('--task_id', default='task', help='Task ID to name the CSV')
    p.add_argument('--recursive', action='store_true', help='Recursive traversal if --sdf_path is a folder')

    p.add_argument('--min_bonds', type=int, default=2, help='Minimum distance (in bonds)')
    p.add_argument('--max_bonds', type=int, default=3, help='Maximum distance (in bonds)')

    diag = p.add_mutually_exclusive_group(required=False)
    diag.add_argument('--diagonal', dest='include_diagonal', action='store_true', help='Include diagonal (default)')
    diag.add_argument('--no_diagonal', dest='include_diagonal', action='store_false', help="Do not include diagonal")
    p.set_defaults(include_diagonal=True)
    return p


def main():
    args = build_argparser().parse_args()
    try:
        ensure_dir_exists(args.output_dir)
        out_csv = os.path.join(args.output_dir, f"prediction_COSY_{args.task_id}.csv")

        files = list_sdf_files(args.sdf_path, recursive=args.recursive)
        if not files:
            error(f"No SDF found in {args.sdf_path}", code=2)

        log(f"SDF files detected: {len(files)}")
        frames: List[pd.DataFrame] = []
        total = len(files)

        for i, sdf in enumerate(files, 1):
            pct = int(5 + ((i - 1) / max(total, 1)) * 85)  # [NEW] 5% -> 90%
            progress(pct, f"Processing {i}/{total}: {os.path.basename(sdf)}")  # [NEW]
            df = dataframe_for_sdf(
                sdf,
                min_bonds=args.min_bonds,
                max_bonds=args.max_bonds,
                include_diagonal=args.include_diagonal
            )
            if not df.empty:
                frames.append(df)

        if not frames:
            error("No COSY table generated (missing shift vector?)", code=3)

        result = pd.concat(frames, ignore_index=True)
        progress(95, "Writing CSV file...")  # [NEW]
        result.to_csv(out_csv, index=False, encoding='utf-8')

        progress(100, "COSY reconstruction completed")  # [NEW]
        completed(out_csv, len(result))

    except SystemExit:
        raise
    except Exception as e:
        jlog({"type": "error", "message": f"Exception: {e}"})
        traceback.print_exc(file=sys.stdout)
        error("COSY reconstruction failed (see trace above)", code=1)


if __name__ == '__main__':
    main()