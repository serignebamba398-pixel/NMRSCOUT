#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
predict_HSQC_real.py (server version, optional 13C folding on F1)
-----------------------------------------------------------------
Reconstructs HSQC from SDF containing an atom-wise shift vector.
The vector is recovered robustly via ml_sdf_utils.get_atom_shift_vector().

Arguments:
    --sdf_path     : .sdf file or folder of .sdf files
    --output_dir   : output folder
    --task_id      : task identifier (CSV naming)
    --recursive    : recursive traversal if --sdf_path is a folder

    # NEW:
    --no-fold      : disables folding (otherwise folded by default for compatibility)
    --fold-low     : F1 lower bound (ppm), default 4.8618
    --fold-high    : F1 upper bound (ppm), default 69.8395
                                     mid is derived: mid = 2*fold_high - fold_low

Output:
    <output_dir>/prediction_HSQC_<task_id>.csv

Columns (0-based indices):
    Filename, C_index, F1_ppm, F1_folded, F2_ppm, Direction, Nb_H, H_indices, ECFP4, Method
"""

import os
import sys
import json
import argparse
import collections
import traceback
from typing import List, Dict, Tuple

import numpy as np
import pandas as pd

from rdkit import Chem, rdBase, RDLogger
from rdkit.Chem import AllChem, rdmolfiles
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from collections import defaultdict

# === robust utilities (neighbor file: ml_sdf_utils.py) ===
from ml_sdf_utils import (
    read_molecule_keepH,
    get_atom_shift_vector,
    ecfp4_bitstring,
)

# -----------------------------------------------------------
# RDKit logs: reduce noise
# -----------------------------------------------------------
RDLogger.DisableLog('rdApp.*')
rdBase.DisableLog('rdApp.error')
rdBase.DisableLog('rdApp.warning')

# ===========================================================
# 0) JSON LOG (stdout) for SSE
# ===========================================================
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

def completed(out_csv: str, n_rows: int, folded: bool, low: float, high: float):
    mode = "folded" if folded else "unfolded"
    jlog({
        "type": "completed",
        "success": True,
        "message": f"HSQC {mode} ({n_rows} rows, low={low}, high={high})",
        "output": out_csv
    })

# ===========================================================
# 1) FILE UTILITIES
# ===========================================================
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

# ===========================================================
# 2) HSQC CORE FUNCTIONS (structure, chirality, symmetry)
# ===========================================================
def mol_with_atom_index(mol: Chem.Mol, include_H: bool = True) -> Chem.Mol:
    if include_H:
        mol = Chem.AddHs(mol, addCoords=True)
    for idx in range(mol.GetNumAtoms()):
        mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(idx))
    return mol

def get_molecule_data(compound_path: str):
    """Extracts atom_list, connectivity_list, name, and a temporary RDKit .mol"""
    docline_list = []
    with open(compound_path, "r", encoding='utf-8', errors='ignore') as f:
        for line in f:
            docline_list.append(line)

    atom_list = []
    start = False
    for line in docline_list:
        if start:
            if "M  END" in line:
                break
            fields = line.split()
            if len(fields) < 4:
                continue
            atom = fields[3]
            if atom != "0":
                atom_list.append(atom)
        if "V2000" in line:
            start = True

    start_line = len(atom_list) + 4
    connectivity_list = []
    for idx, line in enumerate(docline_list):
        if idx >= start_line and "M  END" not in line:
            add_list = line.split()
            if len(add_list) == 3:
                # fallback to fixed-width columns
                part_1 = str(int(line[:3])); part_2 = str(int(line[3:6]))
                part_3 = str(int(line[6:9])); part_4 = str(int(line[9:12]))
                add_list = [part_1, part_2, part_3, part_4]
            if add_list and add_list[0] not in ["M", ">", ""]:
                connectivity_list.append(add_list)
        if "M  END" in line:
            break

    temp_mol_path = compound_path[:-4] + ".mol"
    with open(temp_mol_path, "w", encoding='utf-8', errors='ignore') as output:
        for line in docline_list:
            output.write(line)
            if "M  END" in line:
                break

    mol = Chem.MolFromMolFile(temp_mol_path, removeHs=False, sanitize=True)
    name = os.path.basename(compound_path).split(".")[0]
    return atom_list, connectivity_list, name, mol

def get_c_h_connectivity(connectivity_list: List[List[str]], atom_list_: List[str]) -> Dict[int, List[int]]:
    """Dictionary {C_index -> [bonded_H_indices]} (0-based indices)."""
    c_h_connectivity_dict: Dict[int, List[int]] = {}
    for entry in connectivity_list:
        selected_atom_nr = int(entry[0]) - 1
        selected_connection_nr = int(entry[1]) - 1
        atom = atom_list_[selected_atom_nr] if 0 <= selected_atom_nr < len(atom_list_) else ""
        connection = atom_list_[selected_connection_nr] if 0 <= selected_connection_nr < len(atom_list_) else ""
        if atom == "C" and connection == "H":
            c_h_connectivity_dict.setdefault(selected_atom_nr, []).append(selected_connection_nr)
        if atom == "H" and connection == "C":
            c_h_connectivity_dict.setdefault(selected_connection_nr, []).append(selected_atom_nr)
    return c_h_connectivity_dict

def get_HSQC_info_data(mol: Chem.Mol = None, smiles: str = "", file_path: str = ""):
    """Returns: count_list (C->H), m (indexed Mol), stereo_smi, am (adj.), elements."""
    if mol is not None:
        m = mol_with_atom_index(mol, True)
        am = Chem.GetAdjacencyMatrix(m)
    elif smiles:
        mol = Chem.MolFromSmiles(smiles)
        m = mol_with_atom_index(mol, True)
        am = Chem.GetAdjacencyMatrix(m)
    elif file_path:
        mol = Chem.MolFromMolFile(file_path, removeHs=False, sanitize=True)
        m = mol_with_atom_index(mol, True)
        am = Chem.GetAdjacencyMatrix(m)

    elements = [(atom.GetSymbol(), idx) for idx, atom in enumerate(m.GetAtoms())]
    count_list: Dict[int, List[int]] = {}
    for atom, nr in elements:
        if atom == "C":
            for idx, connector in enumerate(am[nr]):
                if connector == 1 and elements[idx][0] == "H":
                    count_list.setdefault(nr, []).append(idx)

    try:
        isomers = tuple(EnumerateStereoisomers(m))
    except RuntimeError:
        Chem.RemoveStereochemistry(m)
        isomers = tuple(EnumerateStereoisomers(m))
    stereo_smi = Chem.MolToSmiles(isomers[0], isomericSmiles=True)
    return count_list, m, stereo_smi, am, elements

def get_chiral_carbons_from_rdkit(mol: Chem.Mol) -> List[int]:
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    return [idx for idx, _ in centers]

def check_terminal_chiral_center(chiral_C: int, am: np.ndarray, elements: List[Tuple[str, int]]) -> bool:
    surrounding_list = []
    for idx, val in enumerate(am[chiral_C]):
        if val == 1 and idx != chiral_C and elements[idx][0] != "C":
            surrounding_list.append(elements[idx][0])
    return len(set(surrounding_list)) == 3

def chiral_center_type(chiral_centers: List[int], elements: List[Tuple[str, int]], am: np.ndarray):
    direct_neighbor_dict: Dict[int, List[int]] = defaultdict(list)
    terminal_neighbor_dict: Dict[int, List[int]] = defaultdict(list)
    direct_neighbor_list: List[int] = []
    terminal_neighbor_list: List[int] = []
    for chiral_center in chiral_centers:
        terminal = check_terminal_chiral_center(chiral_center, am, elements)
        for idx, val in enumerate(am[chiral_center]):
            if val == 1 and elements[idx][0] == "C":
                if terminal:
                    terminal_neighbor_list.append(idx)
                    terminal_neighbor_dict[chiral_center].append(idx)
                else:
                    direct_neighbor_list.append(idx)
                    direct_neighbor_dict[chiral_center].append(idx)
    return direct_neighbor_list, direct_neighbor_dict, terminal_neighbor_list, terminal_neighbor_dict

def get_exclude_second_degree_list(mol: Chem.Mol, chiral_centers: List[int], elements: List[Tuple[str, int]], am: np.ndarray):
    ringinfo = mol.GetRingInfo()
    rings = ringinfo.AtomRings()
    ring_atoms = [i for ring in rings for i in ring]
    exclude_second_degree_list: List[int] = []
    if chiral_centers:
        for i in chiral_centers:
            for idx, j in enumerate(am[i]):
                if j == 1 and elements[idx][0] == "C" and idx in ring_atoms:
                    exclude_second_degree_list.append(i)
    return exclude_second_degree_list, []

def get_second_degree_neighbor(direct_neighbor_dict: Dict[int, List[int]], exclude_second_degree_list: List[int],
                               am: np.ndarray, elements: List[Tuple[str, int]],
                               direct_neighbor_list: List[int], chiral_centers: List[int],
                               terminal_neighbor_list: List[int]):
    import copy
    direct_neighbor_dict_ = copy.deepcopy(direct_neighbor_dict)
    for i in exclude_second_degree_list:
        if i in direct_neighbor_dict_:
            del direct_neighbor_dict_[i]
    sec_direct_neighbor: List[int] = []
    for key, values in direct_neighbor_dict_.items():
        for neighbor in values:
            for idx, val in enumerate(am[neighbor]):
                if val == 1 and elements[idx][0] == "C" and idx not in chiral_centers and idx not in direct_neighbor_list and idx not in terminal_neighbor_list:
                    direct_neighbor_dict[key].append(idx)
                    sec_direct_neighbor.append(idx)
    return sec_direct_neighbor, direct_neighbor_dict

def get_all_chiral_split_positions(direct_neighbor_dict: Dict[int, List[int]], terminal_neighbor_dict: Dict[int, List[int]]) -> List[int]:
    all_split_positions: List[int] = []
    for v in direct_neighbor_dict.values():
        all_split_positions.extend(v)
    for v in terminal_neighbor_dict.values():
        all_split_positions.extend(v)
    return all_split_positions

def check_for_symmetric_substituents(all_split_positions: List[int], chiral_centers: List[int], elements: List[Tuple[str, int]], am: np.ndarray) -> List[int]:
    remove_neighbor_list: List[int] = []
    positions_copy = list(all_split_positions)
    for pos in positions_copy:
        surrounding_list = []
        for idx, val in enumerate(am[pos]):
            if val == 1 and elements[idx][1] not in chiral_centers and elements[idx][1] not in all_split_positions:
                surrounding_list.append(elements[idx][0])
        if len(set(surrounding_list)) == 1 and len(surrounding_list) == 3:
            remove_neighbor_list.append(pos)
            if pos in all_split_positions:
                all_split_positions.remove(pos)
    return remove_neighbor_list

def remove_symmetric_centers(all_split_positions: List[int], remove_neighbor_list: List[int], sec_direct_neighbor: List[int],
                             elements: List[Tuple[str, int]], am: np.ndarray) -> List[int]:
    for pos in remove_neighbor_list:
        if pos in sec_direct_neighbor and pos in all_split_positions:
            for idx, val in enumerate(am[pos]):
                if val == 1 and elements[idx][0] == "C" and idx != pos and idx in all_split_positions:
                    all_split_positions.remove(idx)
    return all_split_positions

def check_for_two_H(all_split_positions: List[int], elements: List[Tuple[str, int]], am: np.ndarray) -> List[int]:
    positions_copy = list(all_split_positions)
    for pos in positions_copy:
        num_h = 0
        for idx, val in enumerate(am[pos]):
            if val == 1 and elements[idx][0] == "H":
                num_h += 1
        if num_h != 2:
            if pos in all_split_positions:
                all_split_positions.remove(pos)
    return all_split_positions

def find_symmetric_positions(stereo_smi: str) -> List[List[int]]:
    mol = Chem.MolFromSmiles(stereo_smi)
    if mol is None:
        return []
    z = list(rdmolfiles.CanonicalRankAtoms(mol, breakTies=False))
    matches = mol.GetSubstructMatches(mol, uniquify=False)
    if len(z) != len(set(z)) and len(matches) > 1:
        symmetric = [item for item, count in collections.Counter(z).items() if count > 1]
        all_duplicates = []
        for j in symmetric:
            indices = [i for i, v in enumerate(z) if v == j]
            all_duplicates.append(indices)
        example_match = matches[0]
        sym_dupl_lists = []
        for sub_list in all_duplicates:
            indices_list = []
            for i in sub_list:
                indices_list.append(example_match[i])
            sym_dupl_lists.append(indices_list)
    else:
        sym_dupl_lists = []
    return sym_dupl_lists

def run_chiral_and_symmetry_finder(smiles: str = "", compound_path: str = "", extra_info: bool = False):
    """Returns: sym_dupl_lists, all_split_positions, mol, compound_path."""
    if compound_path:
        mol = Chem.MolFromMolFile(compound_path, removeHs=False, sanitize=True)
        mol = Chem.AddHs(mol, addCoords=True)
    elif smiles:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol, addCoords=True)
        compound_path = "temp.sdf"
    count_list, m, stereo_smi, am, elements = get_HSQC_info_data(mol=mol)
    chiral_centers = get_chiral_carbons_from_rdkit(mol)
    direct_neighbor_list, direct_neighbor_dict, terminal_neighbor_list, terminal_neighbor_dict = chiral_center_type(chiral_centers, elements, am)
    exclude_second_degree_list, _ = get_exclude_second_degree_list(mol, chiral_centers, elements, am)
    sec_direct_neighbor, direct_neighbor_dict = get_second_degree_neighbor(
        direct_neighbor_dict, exclude_second_degree_list, am, elements, direct_neighbor_list, chiral_centers, terminal_neighbor_list
    )
    all_split_positions = get_all_chiral_split_positions(direct_neighbor_dict, terminal_neighbor_dict)
    remove_neighbor_list = check_for_symmetric_substituents(all_split_positions, chiral_centers, elements, am)
    all_split_positions = remove_symmetric_centers(all_split_positions, remove_neighbor_list, sec_direct_neighbor, elements, am)
    all_split_positions = check_for_two_H(all_split_positions, elements, am)
    sym_dupl_lists = find_symmetric_positions(stereo_smi)
    if extra_info:
        print("Chiral centers:", chiral_centers)
        print("All split positions:", all_split_positions)
        print("Symmetric groups:", sym_dupl_lists)
    return sym_dupl_lists, all_split_positions, mol, compound_path

# ===========================================================
# 3) SELECTION / DEDUPLICATION / DATAFRAME
# ===========================================================
def selecting_shifts(c_h_connectivity_dict: Dict[int, List[int]],
                     ml_shifts: List[float]) -> List[List[object]]:
    """
    Associates the 13C shift of the carbon (ml_shifts[c_pos]) and the 1H shift
    (average or individual if CH2 is diastereotopic).
    Returns entries: [C_shift, H_shift, direction, C_index, H_indices_list]

    The SDF is already corrected:
    - Equivalent H have identical shifts → averaged together
    - Diastereotopic H have different shifts → separate entries
    """
    shifts: List[List[object]] = []
    for c_pos, h_list in c_h_connectivity_dict.items():
        if not h_list:
            continue
        c_shift = ml_shifts[c_pos] if 0 <= c_pos < len(ml_shifts) else None
        direction = 1.0 if (len(h_list) % 2) != 0 else -1.0

        # Pour CH2: split si les shifts sont différents (diastéréotopiques)
        should_split = False
        if len(h_list) == 2:
            h1, h2 = h_list[0], h_list[1]
            if 0 <= h1 < len(ml_shifts) and 0 <= h2 < len(ml_shifts):
                try:
                    shift1 = float(ml_shifts[h1])
                    shift2 = float(ml_shifts[h2])
                    # If shifts are different, they are diastereotopic
                    if shift1 != shift2:
                        should_split = True
                except (ValueError, TypeError):
                    pass

        if should_split and len(h_list) == 2:
            for h in h_list:
                h_shift = ml_shifts[h] if 0 <= h < len(ml_shifts) else None
                shifts.append([c_shift, h_shift, direction, c_pos, h_list])
        else:
            valid_h = [ml_shifts[h] for h in h_list if 0 <= h < len(ml_shifts)]
            h_avg = float(np.mean(valid_h)) if valid_h else None
            shifts.append([c_shift, h_avg, direction, c_pos, h_list])
    return shifts

def perform_deduplication_if_symmetric(shifts: List[List[object]],
                                       sym_dupl_lists: List[List[int]]) -> List[List[object]]:
    """If some carbons are symmetric, average their entries."""
    if not sym_dupl_lists:
        return shifts
    new_shifts: List[List[object]] = []
    c_done: List[int] = []
    for carbons in sym_dupl_lists:
        carb_sum = hyd_sum = 0.0
        count = 0
        any_dir = 1.0
        c_pos_agg: List[int] = []
        h_pos_agg: List[List[int]] = []
        for C_shift, H_shift, direction, c_pos, h_list in shifts:
            if c_pos in carbons:
                if C_shift is not None: carb_sum += C_shift
                if H_shift is not None: hyd_sum += H_shift
                count += 1
                any_dir = direction
                c_pos_agg.append(c_pos)
                h_pos_agg.append(h_list)
                c_done.append(c_pos)
        if count > 0:
            c_avg = carb_sum / count
            h_avg = hyd_sum / count if count > 0 else None
            new_shifts.append([c_avg, h_avg, any_dir, c_pos_agg, h_pos_agg])
    for C_shift, H_shift, direction, c_pos, h_list in shifts:
        if c_pos not in c_done:
            new_shifts.append([C_shift, H_shift, direction, [c_pos], h_list])
    return new_shifts

def generate_ml_dataframe(shifts: List[List[object]]) -> pd.DataFrame:
    """DataFrame (without ECFP4 nor folding)."""
    rows = []
    for C_shift, H_shift, direction, c_pos_or_list, h_list in shifts:
        if isinstance(c_pos_or_list, list):
            c_index = c_pos_or_list[0] if len(c_pos_or_list) == 1 else -1
        else:
            c_index = c_pos_or_list
        if h_list and isinstance(h_list[0], list):
            flat_h = sorted({h for sub in h_list for h in sub})
        else:
            flat_h = sorted(h_list) if h_list else []
        rows.append({
            'C_index': c_index,
            'F1_ppm': C_shift if C_shift is not None else np.nan,   # 13C
            'F2_ppm': H_shift if H_shift is not None else np.nan,   # 1H
            'Direction': float(direction),
            'Nb_H': len(flat_h),
            'H_indices': '|'.join(map(str, flat_h)) if flat_h else ''
        })
    return pd.DataFrame(rows)

# ===========================================================
# 4) 13C FOLDING (F1 AXIS) — configurable
# ===========================================================
def fold_13C_shift(vi: float, low: float, high: float) -> float:
    """
    Parametric folding rules (SYMAPS accordion):

      mid = 2*high - low

      - Vi ≤ high            → Vf = Vi
      - high < Vi ≤ mid      → Vf = high − (Vi − high)
      - Vi > mid             → Vf = low + (Vi − mid)
    """
    mid = (2.0 * float(high)) - float(low)
    vi = float(vi)
    if vi <= high:
        return vi
    elif vi <= mid:
        return float(high) - (vi - float(high))
    else:
        return float(low) + (vi - float(mid))

# ===========================================================
# 5) PIPE POUR UN SDF
# ===========================================================
def build_hsqc_table_from_sdf(sdf_path: str, do_fold: bool, fold_low: float, fold_high: float) -> pd.DataFrame:
    """Builds the HSQC table for a single SDF (optional 13C F1 folding)."""
    try:
        progress(10, f"Reading {os.path.basename(sdf_path)}...")  # [NEW]
        mol = read_molecule_keepH(sdf_path)
        if mol is None:
            log(f"[ERROR] Unreadable molecule: {os.path.basename(sdf_path)}")
            return pd.DataFrame()

        # robust shift vector (regardless of the property name)
        ml_shifts = get_atom_shift_vector(sdf_path)
        if not ml_shifts:
            log(f"[WARN] Shift vector not found: {os.path.basename(sdf_path)}")
            return pd.DataFrame()

        # C–H connectivity from the SDF
        progress(25, "Analyzing C–H connectivity...")  # [NEW]
        atom_list, connectivity_list, _, _ = get_molecule_data(sdf_path)
        c_h_connectivity_dict = get_c_h_connectivity(connectivity_list, atom_list)

        # symmetry detection (for deduplication of symmetric carbons)
        progress(40, "Detecting symmetry...")  # [NEW]
        sym_dupl_lists, _, _, _ = run_chiral_and_symmetry_finder(compound_path=sdf_path)

        # select shifts - auto-split CH2 if shifts different (diastereotopic already corrected in SDF)
        shifts = selecting_shifts(c_h_connectivity_dict, ml_shifts)
        shifts = perform_deduplication_if_symmetric(shifts, sym_dupl_lists)
        progress(60, "Computing HSQC correlations...")  # [NEW]

        df = generate_ml_dataframe(shifts)
        if df.empty:
            return df

        # Conversion propre + application optionnelle du folding sur F1 (13C)
        df['F1_ppm'] = pd.to_numeric(df['F1_ppm'], errors='coerce')

        if do_fold:
            mask = df['F1_ppm'].notna()
            df.loc[mask, 'F1_folded'] = df.loc[mask, 'F1_ppm'].map(lambda v: fold_13C_shift(v, fold_low, fold_high))
        else:
            # Unfolded mode: copy to keep a stable column interface
            df['F1_folded'] = df['F1_ppm']

        # ECFP4 (bitstring) on the molecule with explicit H
        progress(75, "Computing molecular fingerprints...")  # [NEW]
        try:
            mol0 = Chem.MolFromMolFile(sdf_path, removeHs=False, sanitize=True)
            if mol0 is None:
                mol0 = mol
            ecfp4 = ecfp4_bitstring(Chem.AddHs(mol0, addCoords=True))
        except Exception:
            ecfp4 = "0" * 1024

        df['ECFP4'] = ecfp4
        df['Method'] = 'HSQC_from_ML_SDF'
        df['Filename'] = os.path.basename(sdf_path)

        cols = ['Filename', 'C_index', 'F1_ppm', 'F1_folded', 'F2_ppm', 'Direction', 'Nb_H', 'H_indices', 'ECFP4', 'Method']
        df = df[cols]
        return df

    except Exception as e:
        log(f"[ERROR] Error on {os.path.basename(sdf_path)}: {e}")
        traceback.print_exc(file=sys.stdout)
        return pd.DataFrame()

# ===========================================================
# 6) CLI (ARGS SERVEUR)
# ===========================================================
def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Reconstruct HSQC (optional 13C F1 folding + ECFP4) from an SDF file or folder (shift vector is auto-detected)."
    )
    p.add_argument('--sdf_path', required=True, help='SDF file or folder of SDF files')
    p.add_argument('--output_dir', required=True, help='Output folder')
    p.add_argument('--task_id', default='task', help='Task ID (CSV naming)')
    p.add_argument('--recursive', action='store_true', help='Recursive traversal if --sdf_path is a folder')

    # --- NEW: folding control ---
    p.add_argument('--no-fold',
                   dest='fold',
                   action='store_false',
                   help='Disable HSQC folding (F1 stays unfolded).')
    p.set_defaults(fold=True)  # compat: folded by default

    p.add_argument('--fold-low', type=float, default=4.8618,
                    help='F1 lower bound (ppm) for folding (default: 4.8618).')
    p.add_argument('--fold-high', type=float, default=69.8395,
                    help='F1 upper bound (ppm) for folding (default: 69.8395). '
                        'Mid is computed as: mid = 2*high - low.')
    return p

def main():
    args = build_argparser().parse_args()
    try:
        ensure_dir_exists(args.output_dir)
        out_csv = os.path.join(args.output_dir, f"prediction_HSQC_{args.task_id}.csv")

        # basic bounds validation
        fold_low = float(args.fold_low)
        fold_high = float(args.fold_high)
        if not np.isfinite(fold_low) or not np.isfinite(fold_high):
            error("Folding bounds must be numeric.", code=2)
        if fold_high <= fold_low:
            error("Invalid folding constraints: fold-high must be > fold-low.", code=2)

        mode_txt = "folded" if args.fold else "unfolded"
        log(f"HSQC mode: {mode_txt} — low={fold_low:.4f}, high={fold_high:.4f}")

        files = list_sdf_files(args.sdf_path, recursive=args.recursive)
        if not files:
            error(f"No SDF found in {args.sdf_path}", code=3)

        log(f"SDF files detected: {len(files)}")
        all_tables: List[pd.DataFrame] = []
        total = len(files)

        for i, sdf in enumerate(files, 1):
            pct = int(5 + ((i - 1) / max(total, 1)) * 85)  # [NEW] 5% -> 90%
            progress(pct, f"Processing {i}/{total}: {os.path.basename(sdf)}")  # [NEW]
            df = build_hsqc_table_from_sdf(sdf, args.fold, fold_low, fold_high)
            if not df.empty:
                all_tables.append(df)

        if not all_tables:
            error("No HSQC table generated (missing shift vector?)", code=4)

        result = pd.concat(all_tables, ignore_index=True)
        progress(95, "Writing CSV file...")  # [NEW]
        result.to_csv(out_csv, index=False, encoding='utf-8')

        progress(100, "HSQC extraction completed")  # [NEW]
        completed(out_csv, len(result), args.fold, fold_low, fold_high)

    except SystemExit:
        raise
    except Exception as e:
        jlog({"type": "error", "message": f"Exception: {e}"})
        traceback.print_exc(file=sys.stdout)
        error("HSQC reconstruction failed (see trace above)", code=1)

if __name__ == '__main__':
    main()