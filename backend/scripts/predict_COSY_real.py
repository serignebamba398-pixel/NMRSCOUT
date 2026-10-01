#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Topological COSY and TOCSY correlation prediction from SDF files."""
import os
import sys
import json
import argparse
import traceback
from typing import Dict, List, Tuple, Optional
from collections import defaultdict, deque
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger, rdBase
from ml_sdf_utils import (
    read_molecule_keepH,
    get_atom_shift_vector,
    ecfp4_bitstring,
)
RDLogger.DisableLog('rdApp.*')
rdBase.DisableLog('rdApp.error')
rdBase.DisableLog('rdApp.warning')

def jlog(payload: dict):
    print(json.dumps(payload, ensure_ascii=False), flush=True)

def log(message: str):
    jlog({
        "type": "log",
        "message": message
    })

def progress(pct: float, message: str = None):
    p = max(0, min(100, int(pct)))
    obj = {
        "type": "progress",
        "current": p,
        "total": 100,
        "percent": p
    }
    if message:
        obj["message"] = message
    jlog(obj)

def error(message: str, code: int = 1):
    jlog({
        "type": "error",
        "message": message
    })
    jlog({
        "type": "completed",
        "success": False,
        "error": message
    })
    sys.exit(code)

def completed(
    out_csv: str,
    n_rows: int,
    experiment: str
):
    jlog({
        "type": "completed",
        "success": True,
        "message": f"{experiment} reconstructed ({n_rows} rows)",
        "output": out_csv,
        "experiment": experiment
    })

# File utilities
def list_sdf_files(
    sdf_path: str,
    recursive: bool = False
) -> List[str]:
    if not os.path.exists(sdf_path):
        return []
    if os.path.isfile(sdf_path):
        return [sdf_path]
    files: List[str] = []
    if recursive:
        for root, _, fs in os.walk(sdf_path):
            for f in fs:
                if f.lower().endswith('.sdf'):
                    files.append(
                        os.path.join(root, f)
                    )
    else:
        for f in os.listdir(sdf_path):
            if f.lower().endswith('.sdf'):
                files.append(
                    os.path.join(sdf_path, f)
                )
    return sorted(files)

def ensure_dir_exists(d: str):
    if d and not os.path.exists(d):
        os.makedirs(
            d,
            exist_ok=True
        )

# Proton grouping
def build_parent_to_hydrogens(
    mol: Chem.Mol
) -> Dict[int, List[int]]:
    """Map carbon atom indices to attached hydrogen indices."""
    parent2Hs: Dict[int, List[int]] = defaultdict(list)
    for h in mol.GetAtoms():
        if h.GetSymbol() != 'H':
            continue
        nbrs = h.GetNeighbors()
        if not nbrs:
            continue
        if any(
            nb.GetSymbol() in ('O', 'N', 'S')
            for nb in nbrs
        ):
            continue
        for nb in nbrs:
            if nb.GetSymbol() == 'C':
                parent2Hs[
                    nb.GetIdx()
                ].append(
                    h.GetIdx()
                )
                break
    return parent2Hs

def get_h_groups(
    mol: Chem.Mol,
    shifts: List[float] = None
) -> Tuple[
    Dict[int, int],
    Dict[int, List[int]]
]:
    """Build proton signal groups and atom-to-group mappings."""
    parent2Hs = build_parent_to_hydrogens(mol)
    groups: List[List[int]] = []
    used = set()
    for hs in parent2Hs.values():
        if len(hs) == 3:
            groups.append(
                list(hs)
            )
            used.update(hs)
    for hs in parent2Hs.values():
        if (
            len(hs) == 2
            and not any(h in used for h in hs)
        ):
            h1, h2 = hs
            should_split = False
            if shifts is not None:
                if (
                    0 <= h1 < len(shifts)
                    and
                    0 <= h2 < len(shifts)
                ):
                    try:
                        shift1 = float(
                            shifts[h1]
                        )
                        shift2 = float(
                            shifts[h2]
                        )
                        if shift1 != shift2:
                            should_split = True
                    except (
                        ValueError,
                        TypeError
                    ):
                        pass
            if should_split:
                groups.append([h1])
                groups.append([h2])
            else:
                groups.append(
                    list(hs)
                )
            used.update(hs)
    for hs in parent2Hs.values():
        if (
            len(hs) == 1
            and
            hs[0] not in used
        ):
            groups.append(
                [hs[0]]
            )
            used.add(
                hs[0]
            )
    h2group: Dict[int, int] = {}
    group2Hs: Dict[int, List[int]] = {}
    for gi, grp in enumerate(groups):
        group2Hs[gi] = grp
        for h in grp:
            h2group[h] = gi
    return (
        h2group,
        group2Hs
    )

# COSY
def get_cosy_pairs(
    mol: Chem.Mol,
    h2group: Dict[int, int],
    min_bonds: int = 2,
    max_bonds: int = 3,
    include_diagonal: bool = True
) -> List[
    Tuple[int, int, int, int, Optional[int]]
]:
    """Return topological COSY proton pairs."""
    cosy = []
    h_list = sorted(
        h2group.keys()
    )
    for i, h1 in enumerate(h_list):
        if include_diagonal:
            cosy.append((
                h1,
                h1,
                0,
                0,
                None
            ))
        for h2 in h_list[i + 1:]:
            if (
                h2group[h1]
                ==
                h2group[h2]
            ):
                continue
            path = Chem.GetShortestPath(
                mol,
                h1,
                h2
            )
            if not path:
                continue
            dist = len(path) - 1
            if (
                min_bonds
                <= dist
                <= max_bonds
            ):
                cosy.append((
                    h1,
                    h2,
                    dist,
                    1,
                    None
                ))
                cosy.append((
                    h2,
                    h1,
                    dist,
                    1,
                    None
                ))
    return cosy

# TOCSY spin graph
def _group_parent_carbon(
    mol: Chem.Mol,
    group2Hs: Dict[int, List[int]]
) -> Dict[int, Optional[int]]:
    """Map proton groups to their parent carbon atoms."""
    group_parent: Dict[int, Optional[int]] = {}
    for group_id, hydrogens in group2Hs.items():
        parent_carbon = None
        for h_idx in hydrogens:
            atom = mol.GetAtomWithIdx(h_idx)
            for nb in atom.GetNeighbors():
                if nb.GetSymbol() == 'C':
                    parent_carbon = nb.GetIdx()
                    break
            if parent_carbon is not None:
                break
        group_parent[group_id] = parent_carbon
    return group_parent

def build_tocsy_spin_graph(
    mol: Chem.Mol,
    group2Hs: Dict[int, List[int]],
    min_bonds: int = 2,
    max_bonds: int = 3
) -> Dict[int, set]:
    """Build the structural TOCSY spin graph."""
    _ = min_bonds, max_bonds
    groups = sorted(group2Hs.keys())
    adjacency = {g: set() for g in groups}
    group_parent = _group_parent_carbon(
        mol,
        group2Hs
    )
    for i, g1 in enumerate(groups):
        c1 = group_parent.get(g1)
        if c1 is None:
            continue
        for g2 in groups[i + 1:]:
            c2 = group_parent.get(g2)
            if c2 is None:
                continue
            connected = False
            if c1 == c2:
                connected = True
            else:
                if mol.GetBondBetweenAtoms(c1, c2) is not None:
                    connected = True
            if connected:
                adjacency[g1].add(g2)
                adjacency[g2].add(g1)
    return adjacency

# TOCSY spin systems
def find_spin_systems(
    adjacency: Dict[int, set]
) -> Tuple[
    Dict[int, int],
    Dict[int, List[int]]
]:
    """Find connected TOCSY spin systems."""
    visited = set()
    group_to_system: Dict[int, int] = {}
    system_to_groups: Dict[int, List[int]] = {}
    system_id = 0
    for start in sorted(adjacency.keys()):
        if start in visited:
            continue
        stack = [start]
        visited.add(start)
        members = []
        while stack:
            current = stack.pop()
            members.append(current)
            group_to_system[current] = system_id
            for neighbor in adjacency[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        system_to_groups[system_id] = sorted(
            members
        )
        system_id += 1
    return (
        group_to_system,
        system_to_groups
    )

# Relay distance
def shortest_relay_steps(
    adjacency: Dict[int, set],
    start_group: int,
    target_group: int
) -> Optional[int]:
    """Return the shortest relay path between proton groups."""
    if start_group == target_group:
        return 0
    visited = {
        start_group
    }
    queue = deque([
        (
            start_group,
            0
        )
    ])
    while queue:
        current, depth = queue.popleft()
        for neighbor in adjacency[current]:
            if neighbor == target_group:
                return depth + 1
            if neighbor not in visited:
                visited.add(
                    neighbor
                )
                queue.append((
                    neighbor,
                    depth + 1
                ))
    return None

# TOCSY pairs
def get_tocsy_pairs(
    mol: Chem.Mol,
    h2group: Dict[int, int],
    group2Hs: Dict[int, List[int]],
    min_bonds: int = 2,
    max_bonds: int = 3,
    include_diagonal: bool = True
) -> List[
    Tuple[int, int, int, int, int]
]:
    """Return topological TOCSY proton pairs."""
    adjacency = build_tocsy_spin_graph(
        mol,
        group2Hs,
        min_bonds=min_bonds,
        max_bonds=max_bonds
    )
    (
        group_to_system,
        system_to_groups
    ) = find_spin_systems(
        adjacency
    )
    log(f"[TOCSY] Detected {len(system_to_groups)} spin systems")
    for sid, groups in system_to_groups.items():
        log(
            f"[TOCSY] SpinSystem {sid}: "
            f"{len(groups)} proton group(s) -> {groups}"
        )
    log(
        f"[TOCSY] Retained H-H distance: "
        f"{min_bonds}-{max_bonds} bonds"
    )
    tocsy = []
    h_list = sorted(h2group.keys())
    for i, h1 in enumerate(h_list):
        g1 = h2group[h1]
        system1 = group_to_system[g1]
        if include_diagonal:
            tocsy.append((
                h1,
                h1,
                0,
                0,
                system1
            ))
        for h2 in h_list[i + 1:]:
            g2 = h2group[h2]
            if g1 == g2:
                continue
            system2 = group_to_system[g2]
            if system1 != system2:
                continue
            path = Chem.GetShortestPath(
                mol,
                h1,
                h2
            )
            if not path:
                continue
            bond_distance = len(path) - 1
            if not (
                min_bonds
                <= bond_distance
                <= max_bonds
            ):
                continue
            relay_steps = shortest_relay_steps(
                adjacency,
                g1,
                g2
            )
            if relay_steps is None:
                continue
            tocsy.append((
                h1,
                h2,
                bond_distance,
                relay_steps,
                system1
            ))
            tocsy.append((
                h2,
                h1,
                bond_distance,
                relay_steps,
                system1
            ))
    log(f"[TOCSY] retained {sum(a != b for a, b, *_ in tocsy)} cross-peak rows (both directions)")
    return tocsy

# Signal mapping
def make_signal_map_per_H(
    shifts: List[float],
    h2group: Dict[int, int],
    group2Hs: Dict[int, List[int]]
) -> Dict[int, float]:
    """Map hydrogen atom indices to chemical shifts."""
    group_avg: Dict[int, float] = {}
    for gi, members in group2Hs.items():
        vals = [
            shifts[h]
            for h in members
            if 0 <= h < len(shifts)
        ]
        group_avg[gi] = (
            float(np.mean(vals))
            if vals
            else float('nan')
        )
    signal_map: Dict[int, float] = {}
    for h, gi in h2group.items():
        v = group_avg.get(
            gi,
            float('nan')
        )
        if not np.isnan(v):
            signal_map[h] = round(
                v,
                3
            )
        elif 0 <= h < len(shifts):
            try:
                signal_map[h] = round(
                    float(shifts[h]),
                    3
                )
            except Exception:
                signal_map[h] = float(
                    'nan'
                )
        else:
            signal_map[h] = float(
                'nan'
            )
    return signal_map

# SDF processing
def dataframe_for_sdf(
    sdf_path: str,
    experiment: str,
    min_bonds: int,
    max_bonds: int,
    include_diagonal: bool
) -> pd.DataFrame:
    try:
        filename = os.path.basename(
            sdf_path
        )
        progress(
            10,
            f"Reading {filename}..."
        )
        mol = read_molecule_keepH(
            sdf_path
        )
        if mol is None:
            log(
                f"[ERROR] Unreadable molecule: "
                f"{filename}"
            )
            return pd.DataFrame()
        shifts = get_atom_shift_vector(
            sdf_path
        )
        if not shifts:
            log(
                f"[WARN] Shift vector not found: "
                f"{filename}"
            )
        progress(
            30,
            "Grouping protons..."
        )
        (
            h2group,
            group2Hs
        ) = get_h_groups(
            mol,
            shifts=shifts
        )
        if not h2group:
            log(
                f"[WARN] No eligible protons "
                f"(excluding OH/NH/SH): "
                f"{filename}"
            )
            return pd.DataFrame()
        if experiment == "COSY":
            pairs = get_cosy_pairs(
                mol,
                h2group,
                min_bonds=min_bonds,
                max_bonds=max_bonds,
                include_diagonal=include_diagonal
            )
            method = "COSY_from_ML_SDF"
        elif experiment == "TOCSY":
            pairs = get_tocsy_pairs(
                mol,
                h2group,
                group2Hs,
                min_bonds=min_bonds,
                max_bonds=max_bonds,
                include_diagonal=include_diagonal
            )
            method = "TOCSY_topological_from_ML_SDF"
        else:
            raise ValueError(
                f"Unknown experiment: {experiment}"
            )
        progress(
            50,
            f"Computing {experiment} correlations "
            f"({len(pairs)} pairs)..."
        )
        signal_map = make_signal_map_per_H(
            shifts,
            h2group,
            group2Hs
        )
        ecfp4 = ecfp4_bitstring(
            mol
        )
        rows = []
        for (
            h1,
            h2,
            dist,
            relay_steps,
            spin_system
        ) in pairs:
            d1 = signal_map.get(
                h1,
                float('nan')
            )
            d2 = signal_map.get(
                h2,
                float('nan')
            )
            row = {
                'Filename':
                    filename,
                'H1':
                    h1,
                'Delta1':
                    d1
                    if np.isfinite(d1)
                    else '',
                'H2':
                    h2,
                'Delta2':
                    d2
                    if np.isfinite(d2)
                    else '',
                'Bonds':
                    dist,
                'Group1':
                    h2group.get(
                        h1,
                        -1
                    ),
                'Group2':
                    h2group.get(
                        h2,
                        -1
                    ),
                'ECFP4':
                    ecfp4,
                'Method':
                    method,
                'RelaySteps':
                    relay_steps,
                'SpinSystem':
                    (
                        spin_system
                        if spin_system is not None
                        else ''
                    )
            }
            if experiment == 'TOCSY':
                row['MaxBonds'] = max_bonds
            rows.append(row)
        progress(
            70,
            "Building DataFrame..."
        )
        df = pd.DataFrame(
            rows
        )
        if not df.empty:
            df['__diag__'] = (
                df['H1']
                ==
                df['H2']
            ).astype(int)
            if experiment == "TOCSY":
                df.sort_values(
                    by=[
                        '__diag__',
                        'SpinSystem',
                        'RelaySteps',
                        'Delta1',
                        'Delta2'
                    ],
                    ascending=[
                        False,
                        True,
                        True,
                        True,
                        True
                    ],
                    inplace=True
                )
            else:
                df.sort_values(
                    by=[
                        '__diag__',
                        'Bonds',
                        'Delta1'
                    ],
                    ascending=[
                        False,
                        True,
                        True
                    ],
                    inplace=True
                )
            df.drop(
                columns='__diag__',
                inplace=True,
                errors='ignore'
            )
        return df
    except Exception as e:
        log(
            f"[ERROR] Error on "
            f"{os.path.basename(sdf_path)}: "
            f"{e}"
        )
        traceback.print_exc(
            file=sys.stdout
        )
        return pd.DataFrame()

def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Reconstruct topological COSY or TOCSY "
            "from SDF with detected 1H shift vector."
        )
    )
    p.add_argument(
        '--sdf_path',
        required=True,
        help='SDF file or folder of SDF files'
    )
    p.add_argument(
        '--output_dir',
        required=True,
        help='Output folder'
    )
    p.add_argument(
        '--task_id',
        default='task',
        help='Task ID used to name the CSV'
    )
    p.add_argument(
        '--experiment',
        type=str.upper,
        choices=[
            'COSY',
            'TOCSY'
        ],
        default='COSY',
        help=(
            'Experiment type: COSY or TOCSY. '
            'Default: COSY'
        )
    )
    p.add_argument(
        '--recursive',
        action='store_true',
        help=(
            'Recursive traversal if '
            '--sdf_path is a folder'
        )
    )
    p.add_argument(
        '--min_bonds',
        type=int,
        default=2,
        help=(
            'Minimum retained H-H bond distance (diagonal excluded). '
            'Default: 2'
        )
    )
    p.add_argument(
        '--max_bonds',
        type=int,
        default=3,
        help=(
            'Maximum H-H bond distance. '
            'For COSY this is the direct correlation cutoff. '
            'For TOCSY this is the maximum retained H-H distance '
            'inside an independently detected spin system. '
            'Default: 3'
        )
    )
    diag = p.add_mutually_exclusive_group(
        required=False
    )
    diag.add_argument(
        '--diagonal',
        dest='include_diagonal',
        action='store_true',
        help='Include diagonal peaks (default)'
    )
    diag.add_argument(
        '--no_diagonal',
        dest='include_diagonal',
        action='store_false',
        help='Do not include diagonal peaks'
    )
    p.set_defaults(
        include_diagonal=True
    )
    return p

def main():
    parser = build_argparser()
    args = parser.parse_args()
    if args.max_bonds < 1:
        parser.error('--max_bonds must be >= 1')
    if not 1 <= args.min_bonds <= args.max_bonds:
        parser.error('Requires 1 <= min_bonds <= max_bonds')
    experiment = args.experiment.upper()
    try:
        ensure_dir_exists(
            args.output_dir
        )
        out_csv = os.path.join(
            args.output_dir,
            f"prediction_{experiment}_{args.task_id}.csv"
        )
        files = list_sdf_files(
            args.sdf_path,
            recursive=args.recursive
        )
        if not files:
            error(
                f"No SDF found in "
                f"{args.sdf_path}",
                code=2
            )
        log(
            f"Experiment: {experiment}"
        )
        log(
            f"SDF files detected: "
            f"{len(files)}"
        )
        if experiment == "COSY":
            log(
                f"COSY local distance: "
                f"{args.min_bonds}-"
                f"{args.max_bonds} bonds"
            )
        else:
            log(
                f"TOCSY topological mode: "
                f"retained H-H distance = "
                f"{args.min_bonds}-"
                f"{args.max_bonds} bonds"
            )
            log(
                "TOCSY spin systems are determined independently "
                "from max_bonds using structural proton connectivity."
            )
        frames: List[pd.DataFrame] = []
        total = len(files)
        for i, sdf in enumerate(
            files,
            1
        ):
            pct = int(
                5
                +
                (
                    (i - 1)
                    /
                    max(total, 1)
                )
                * 85
            )
            progress(
                pct,
                f"Processing {i}/{total}: "
                f"{os.path.basename(sdf)}"
            )
            df = dataframe_for_sdf(
                sdf,
                experiment=experiment,
                min_bonds=args.min_bonds,
                max_bonds=args.max_bonds,
                include_diagonal=args.include_diagonal
            )
            if not df.empty:
                frames.append(
                    df
                )
        if not frames:
            error(
                f"No {experiment} table generated "
                f"(missing shift vector?)",
                code=3
            )
        result = pd.concat(
            frames,
            ignore_index=True
        )
        progress(
            95,
            "Writing CSV file..."
        )
        result.to_csv(
            out_csv,
            index=False,
            encoding='utf-8'
        )
        progress(
            100,
            f"{experiment} reconstruction completed"
        )
        completed(
            out_csv,
            len(result),
            experiment
        )
    except SystemExit:
        raise
    except Exception as e:
        jlog({
            "type": "error",
            "message": f"Exception: {e}"
        })
        traceback.print_exc(
            file=sys.stdout
        )
        error(
            f"{experiment} reconstruction failed "
            f"(see trace above)",
            code=1
        )

if __name__ == '__main__':
    main()
