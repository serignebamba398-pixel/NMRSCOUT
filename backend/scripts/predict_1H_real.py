#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""predict_1H_real.py (integrated utils version)
------------------------------------------------
Extract 1H signals from SDF files containing an atom-level shift vector stored
in any SDF property (arbitrary name). The vector is retrieved robustly via
ml_sdf_utils.get_atom_shift_vector().

- Preserves RDKit indexing (0-based) and explicit Hs (important).
- Groups protons into CH3/CH2/CH, with optional splitting of diastereotopic CH2.
- Computes an ECFP4 (bitstring) for each molecule.
- Exports a CSV: Filename | Assignment | H_indices | Parent_C_index | GroupType | Delta_ppm | Nb_H | ECFP4

Usage (batch):
    python predict_1H_real.py --input_folder /path/to/sdfs --output_csv /path/out/all_signals_1H.csv
    # + --recursive for subfolders, --split_diastereotopic_ch2 to split CH2

Usage (single):
    python predict_1H_real.py --sdf_path /path/to/file.sdf --output_csv /path/out/signals_1H.csv
"""

import os
import sys
import json
import argparse
import traceback
from typing import Dict, List, Tuple, Iterable, Optional
from collections import defaultdict

import pandas as pd
from rdkit import Chem, RDLogger

# === robust utils for reading + shifts + ECFP4 ===
# Relative path: backend/scripts/ml_sdf_utils.py (adjust if needed)
from ml_sdf_utils import (
    read_molecule_keepH,
    get_atom_shift_vector,
    ecfp4_bitstring,
)

# Reduce RDKit noise
RDLogger.DisableLog('rdApp.*')


# -------------------------------------------------------------------------
# 0) JSON LOGGING FOR SSE ([NEW] ADDED)
# -------------------------------------------------------------------------
def _print_json(obj: dict) -> None:
    try:
        print(json.dumps(obj, ensure_ascii=False), flush=True)
    except Exception:
        print(str(obj), flush=True)

def send_progress(percent: int, message: str = "") -> None:
    p = max(0, min(100, int(percent)))
    _print_json({"type": "progress", "current": p, "total": 100, "percent": p, "message": message})


# -------------------------------------------------------------------------
# 1) I/O utilities
# -------------------------------------------------------------------------
def find_sdf_files(folder: str, recursive: bool = False) -> List[str]:
    """Return the list of .sdf paths in `folder` (recursive if requested)."""
    sdf_paths: List[str] = []
    if recursive:
        for root, _, files in os.walk(folder):
            for f in files:
                if f.lower().endswith('.sdf'):
                    sdf_paths.append(os.path.join(root, f))
    else:
        for f in os.listdir(folder):
            if f.lower().endswith('.sdf'):
                sdf_paths.append(os.path.join(folder, f))
    return sdf_paths


def safe_mkdir_for_file(file_path: str) -> None:
    """Create the parent folder for a file if needed."""
    out_dir = os.path.dirname(os.path.abspath(file_path))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)


# -------------------------------------------------------------------------
# 2) Proton grouping (CH3/CH2/CH) + optional diastereotopic split for CH2
# -------------------------------------------------------------------------
def build_parent_to_hydrogens(mol: Chem.Mol, keep_exchangeable: bool = False) -> Dict[int, List[int]]:
    """
    Return a dict: parent heavy-atom index -> list of attached H atom indices.
    
    Args:
        mol: RDKit molecule
        keep_exchangeable: If False (default), exclude H on heteroatoms (O, N, S)
            to stay within typical aliphatic/aromatic protons.
            If True, include all protons including OH, NH, SH.
    """
    parent2Hs: Dict[int, List[int]] = defaultdict(list)
    for h in mol.GetAtoms():
        if h.GetSymbol() != 'H':
            continue
        nbrs = h.GetNeighbors()
        if not nbrs:
            continue
        # Exclude H bound to O/N/S (OH, NH, SH...) unless keep_exchangeable=True
        if not keep_exchangeable and any(nb.GetSymbol() in ('O', 'N', 'S') for nb in nbrs):
            continue
        for nb in nbrs:
            # If keep_exchangeable, also accept H on O/N/S
            if nb.GetSymbol() == 'C' or (keep_exchangeable and nb.GetSymbol() in ('O', 'N', 'S', 'C')):
                parent2Hs[nb.GetIdx()].append(h.GetIdx())
                break
    return parent2Hs


def is_diastereotopic_pair(mol: Chem.Mol, h1: int, h2: int) -> bool:
    """
    Heuristic: isotope tagging + chirality detection.
    True => diastereotopic protons (should be treated separately).
    """
    try:
        base = Chem.Mol(mol)
        for a in base.GetAtoms():
            a.SetIsotope(0)
        for heavy, light in ((h1, h2), (h2, h1)):
            m = Chem.Mol(base)
            m.GetAtomWithIdx(heavy).SetIsotope(2)
            m.GetAtomWithIdx(light).SetIsotope(1)
            Chem.AssignStereochemistry(m, force=True, cleanIt=True, flagPossibleStereoCenters=True)
            if Chem.FindMolChiralCenters(m, includeUnassigned=False):
                return True
        return False
    except Exception:
        return False


def group_protons_with_parent(
    mol: Chem.Mol,
    split_diastereotopic_ch2: bool = False,
    keep_exchangeable: bool = False
) -> List[Tuple[int, List[int], str]]:
    """
    Return a list of tuples (parent_idx, [H_idx,...], group_type)
    group_type ∈ {"CH3","CH2","CH","XH"} — H on O/N/S are excluded unless keep_exchangeable=True.
    XH = exchangeable proton (OH, NH, SH)
    """
    parent2Hs = build_parent_to_hydrogens(mol, keep_exchangeable=keep_exchangeable)
    groups: List[Tuple[int, List[int], str]] = []
    used = set()

    # CH3
    for c_idx, hs in parent2Hs.items():
        if len(hs) == 3:
            groups.append((c_idx, list(hs), 'CH3'))
            used.update(hs)

    # CH2
    for c_idx, hs in parent2Hs.items():
        if len(hs) == 2 and not any(h in used for h in hs):
            if split_diastereotopic_ch2 and is_diastereotopic_pair(mol, hs[0], hs[1]):
                groups.append((c_idx, [hs[0]], 'CH2'))
                groups.append((c_idx, [hs[1]], 'CH2'))
            else:
                groups.append((c_idx, list(hs), 'CH2'))
            used.update(hs)

    # CH
    for c_idx, hs in parent2Hs.items():
        if len(hs) == 1 and hs[0] not in used:
            groups.append((c_idx, [hs[0]], 'CH'))
            used.add(hs[0])

    return groups


# -------------------------------------------------------------------------
# 3) 1H table for one SDF
# -------------------------------------------------------------------------
def build_1h_table_for_sdf(
    sdf_path: str,
    split_diastereotopic_ch2: bool = False,
    keep_exchangeable: bool = False
) -> pd.DataFrame:
    """
    Build a DataFrame with the following columns:
      Filename | Assignment | H_indices | Parent_C_index | GroupType | Delta_ppm | Nb_H | ECFP4
    - Assignment == H_indices (format "h1,h2,...", for backward compatibility)
    - Delta_ppm = average of the group's H shifts if values exist in the shift vector
    - keep_exchangeable: if True, include exchangeable protons (OH, NH, SH)
    """
    try:
        send_progress(10, f"Reading {os.path.basename(sdf_path)}...")  # [NEW]
        mol = read_molecule_keepH(sdf_path)
        if mol is None:
            print(f"[ERROR] Unreadable SDF: {sdf_path}")
            return pd.DataFrame()

        # Atom-level shift vector (robust to the property name)
        shifts = get_atom_shift_vector(sdf_path)  # 0-based, may be > or < num atoms
        send_progress(30, "Extracting 1H groups...")  # [NEW]
        # 1H groups per parent atom
        groups = group_protons_with_parent(mol, split_diastereotopic_ch2=split_diastereotopic_ch2, keep_exchangeable=keep_exchangeable)
        # ECFP4 (bitstring)
        send_progress(50, "Computing molecular fingerprints...")  # [NEW]
        fp_bitstr = ecfp4_bitstring(mol)

        rows = []
        for parent_c_idx, h_list, gtype in groups:
            # Average shifts for the Hs in the group, if available
            vals = []
            for h_idx in h_list:
                if 0 <= h_idx < len(shifts):
                    try:
                        vals.append(float(shifts[h_idx]))
                    except Exception:
                        pass
            delta = round(sum(vals) / len(vals), 3) if vals else None

            h_str = ','.join(map(str, sorted(h_list)))
            rows.append({
                'Filename':        os.path.basename(sdf_path),
                'Assignment':      h_str,           # compat
                'H_indices':       h_str,           # explicit
                'Parent_C_index':  int(parent_c_idx),
                'GroupType':       gtype,           # CH3 / CH2 / CH
                'Delta_ppm':       delta,
                'Nb_H':            len(h_list),
                'ECFP4':           fp_bitstr
            })

        send_progress(80, "Building DataFrame...")  # [NEW]
        return pd.DataFrame(rows)

    except Exception as e:
        print(f"[ERROR] Error on {sdf_path}: {e}")
        traceback.print_exc(file=sys.stdout)
        return pd.DataFrame()


# -------------------------------------------------------------------------
# 4) Batch folder -> CSV
# -------------------------------------------------------------------------
def process_folder_to_csv(
    input_folder: str,
    output_csv: str,
    recursive: bool = False,
    split_diastereotopic_ch2: bool = False,
    keep_exchangeable: bool = False
) -> int:
    """
    Process all SDF files from a folder and concatenate the tables.
    Returns 0 on success, 1 otherwise.
    """
    sdf_files = find_sdf_files(input_folder, recursive=recursive)
    if not sdf_files:
        print("[ERROR] No .sdf file found.")
        return 1

    all_dfs: List[pd.DataFrame] = []
    for i, sdf in enumerate(sdf_files, 1):  # [NEW] enumerate added
        pct = int(10 + (i / len(sdf_files)) * 80)  # [NEW] 10% -> 90%
        send_progress(pct, f"Processing {i}/{len(sdf_files)}: {os.path.basename(sdf)}")  # [NEW]
        df = build_1h_table_for_sdf(sdf, split_diastereotopic_ch2=split_diastereotopic_ch2, keep_exchangeable=keep_exchangeable)
        if not df.empty:
            all_dfs.append(df)

    if not all_dfs:
        print("[ERROR] No rows produced from any SDF (missing shift vector property?).")
        return 1

    result = pd.concat(all_dfs, ignore_index=True)
    send_progress(95, "Writing CSV file...")  # [NEW]
    safe_mkdir_for_file(output_csv)
    result.to_csv(output_csv, index=False, encoding='utf-8')
    send_progress(100, "1H extraction finished")
    print(f"[OK] 1H export: {output_csv} ({len(result)} rows)")
    return 0


# -------------------------------------------------------------------------
# 5) Single file -> CSV
# -------------------------------------------------------------------------
def process_single_sdf_to_csv(
    sdf_path: str,
    output_csv: str,
    split_diastereotopic_ch2: bool = False,
    keep_exchangeable: bool = False
) -> int:
    df = build_1h_table_for_sdf(sdf_path, split_diastereotopic_ch2=split_diastereotopic_ch2, keep_exchangeable=keep_exchangeable)
    if df.empty:
        print("[ERROR] No 1H data extracted (missing shift vector?).")
        return 1
    send_progress(95, "Writing CSV file...")  # [NEW]
    safe_mkdir_for_file(output_csv)
    df.to_csv(output_csv, index=False, encoding='utf-8')
    send_progress(100, "1H extraction finished")
    print(f"[OK] 1H export: {output_csv} ({len(df)} rows)")
    return 0


# -------------------------------------------------------------------------
# 6) CLI
# -------------------------------------------------------------------------
def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="1H extraction (CH3/CH2/CH, atom indices) + ECFP4 from SDF (file or folder)."
    )
    g_in = p.add_mutually_exclusive_group(required=True)
    g_in.add_argument('--input_folder', type=str, help='Folder containing .sdf files')
    g_in.add_argument('--sdf_path', type=str, help='Single .sdf file')

    p.add_argument('--output_csv', type=str, required=True, help='Output CSV path')
    p.add_argument('--recursive', action='store_true', help='Recursively traverse subfolders')
    p.add_argument('--split_diastereotopic_ch2', action='store_true',
                   help='Split a CH2 into two protons if diastereotopic')
    p.add_argument('--keep_exchangeable', action='store_true',
                   help='Include exchangeable protons (OH, NH, SH). Excluded by default.')
    return p


def main():
    args = build_argparser().parse_args()

    try:
        if args.input_folder:
            code = process_folder_to_csv(
                input_folder=args.input_folder,
                output_csv=args.output_csv,
                recursive=args.recursive,
                split_diastereotopic_ch2=args.split_diastereotopic_ch2,
                keep_exchangeable=args.keep_exchangeable
            )
        else:
            code = process_single_sdf_to_csv(
                sdf_path=args.sdf_path,
                output_csv=args.output_csv,
                split_diastereotopic_ch2=args.split_diastereotopic_ch2,
                keep_exchangeable=args.keep_exchangeable
            )
        sys.exit(code)
    except Exception as e:
        print(f"[FATAL] Fatal error: {e}")
        traceback.print_exc(file=sys.stdout)
        sys.exit(1)


if __name__ == '__main__':
    main()

