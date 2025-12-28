#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""integrated_nmr_prediction.py (server, DGL compatible)
=======================================================
- Generates SDF files enriched with predicted 1H/13C shifts.
- Generates JSON file with full spectral data for 1H NMR (multiplicities, J-couplings).
- CPU-only, Windows compatible.
- Aligns the 1H assignment "per H" with the extractor scripts (predict_1H_real.py / COSY / HSQC):
    * one value per H atom (repeat the CH3/CH2/CH group value)
    * 1H group order = sorted by heavy-atom index (consistent mask -> prediction order)

Outputs:
- SDF files: One per molecule with predicted 1H and 13C shifts
- spectral_data_1H.json: JSON with multiplicity simulation data for spectrum visualization
- sdf_manifest.json: Manifest with all molecule info and paths

Usage (single):
    python integrated_nmr_prediction.py \
        --smiles "CCO" \
        --output_dir ./out \
        --model_1h_path ../models/1H_sparsified_proposed_proposed_set2set_1.pt \
        --model_13c_path ../models/13C_sparsified_proposed_proposed_set2set_1.pt \
        --mol_id ethanol \
        --solvent CDCl3 \
        --spectrometer_mhz 700.0

Usage (CSV):
    python integrated_nmr_prediction.py \
        --csv ./mols.csv \
        --output_dir ./out \
        --model_1h_path ../models/1H_sparsified_proposed_proposed_set2set_1.pt \
        --model_13c_path ../models/13C_sparsified_proposed_proposed_set2set_1.pt \
        --id_column sample-id --smiles_column SMILES \
        --solvent DMSO-d6 --spectrometer_mhz 400.0
"""

from __future__ import annotations

import os
import re
import sys
import csv
import json
import time
import uuid
import shutil
import random
import argparse
import traceback
import importlib.util
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Any

# [OK] Global counter for tracking
_GLOBAL_MOLECULE_COUNTER = {"done": 0, "total": 0}

# [OK] Container for all molecules' spectral JSON data
_MOLECULES_SPECTRAL_DATA: List[Dict[str, Any]] = []

import numpy as np
import torch

# CPU only
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "-1")
os.environ.setdefault("DGLBACKEND", "pytorch")
torch.backends.cudnn.enabled = False

# -----------------------
# RDKit
# -----------------------
from rdkit import Chem, RDConfig, rdBase
from rdkit.Chem import AllChem, ChemicalFeatures, rdchem
try:
    from rdkit.Chem import rdmolfiles as rdMolFiles
except ImportError:
    import rdkit.Chem.rdmolfiles as rdMolFiles

rdBase.DisableLog('rdApp.error')
rdBase.DisableLog('rdApp.warning')

# -----------------------
# DGL (robust import)
# -----------------------
try:
    import dgl
    DGL_VERSION = getattr(dgl, "__version__", "unknown")
    try:
        from dgl import graph as _dgl_graph_fn
        def dgl_graph(edges, num_nodes: int):
            src, dst = edges
            if not torch.is_tensor(src):
                src = torch.as_tensor(src, dtype=torch.int64)
            else:
                src = src.to(torch.int64)
            if not torch.is_tensor(dst):
                dst = torch.as_tensor(dst, dtype=torch.int64)
            else:
                dst = dst.to(torch.int64)
            return _dgl_graph_fn((src, dst), num_nodes=int(num_nodes))
    except Exception:
        try:
            from dgl.convert import graph as _dgl_graph_fn
            def dgl_graph(edges, num_nodes: int):
                src, dst = edges
                src = torch.as_tensor(src, dtype=torch.int64)
                dst = torch.as_tensor(dst, dtype=torch.int64)
                return _dgl_graph_fn((src, dst), num_nodes=int(num_nodes))
        except Exception:
            from dgl import DGLGraph
            def dgl_graph(edges, num_nodes: int):
                g = DGLGraph()
                g.add_nodes(int(num_nodes))
                src, dst = edges
                src = torch.as_tensor(src, dtype=torch.int64)
                dst = torch.as_tensor(dst, dtype=torch.int64)
                g.add_edges(src, dst)
                return g
    print(f"[OK] DGL imported (v{DGL_VERSION})", flush=True)
except Exception as e:
    print(f"DGL import failed: {e}", file=sys.stderr)
    try:
        print(f"PyTorch version: {torch.__version__}", file=sys.stderr)
    except Exception:
        pass
    print("\n[INFO] Install a DGL wheel compatible with your torch version:", file=sys.stderr)
    print("  - torch 2.2.x : pip install dgl==2.2.1 -f https://data.dgl.ai/wheels/torch-2.2/repo.html", file=sys.stderr)
    print("  - torch 2.1.x : pip install 'dgl==2.1.*' -f https://data.dgl.ai/wheels/torch-2.1/repo.html", file=sys.stderr)
    sys.exit(1)

# -----------------------
# Local ML imports
# -----------------------
_THIS = os.path.dirname(__file__)
for p in [
    os.path.join(_THIS, 'ml_dependencies'),
    os.path.join(_THIS, '..', 'ml_dependencies'),
    os.path.join(_THIS, '..')
]:
    if os.path.isdir(p) and p not in sys.path:
        sys.path.insert(0, p)

try:
    from gnn_models.mpnn_proposed import nmr_mpnn_PROPOSED
    from gnn_models.mpnn_baseline import nmr_mpnn_BASELINE
    from model import training, inference
    from util import collate_reaction_graphs
    from dataset import GraphDataset
except ImportError as e:
    print(f"Critical import error: {e}", file=sys.stderr)
    print("Make sure ml_dependencies/ contains gnn_models/, model.py, dataset.py, util.py", file=sys.stderr)
    sys.exit(1)

# -----------------------
# Multiplicity simulator imports
# -----------------------
try:
    from sim_mult_rulebased import (
        simulate_multiplicity_per_proton,
        prepare_molecule_for_nmr,
        SimConfig,
        LABEL2IDX,
        DEFAULT_SOLVENT,
        DEFAULT_SPECTROMETER_FREQ_MHZ,
        DEFAULT_TEMPERATURE_K,
    )
    MULTIPLICITY_SIMULATOR_AVAILABLE = True
    print("[OK] Multiplicity simulator (sim_mult_rulebased) imported", flush=True)
except ImportError as e:
    MULTIPLICITY_SIMULATOR_AVAILABLE = False
    print(f"[WARN] Multiplicity simulator not available: {e}", file=sys.stderr)

# -----------------------
# ML constants
# -----------------------
ATOM_LIST = ['H','Li','B','C','N','O','F','Na','Mg','Al','Si','P','S','Cl','K','Ti','Zn','Ge','As','Se','Br','Pd','Ag','Sn','Sb','Te','I','Hg','Tl','Pb','Bi','Ga']
CHARGE_LIST = [1, 2, 3, -1, -2, -3, 0]
DEGREE_LIST = [1, 2, 3, 4, 5, 6, 0]
VALENCE_LIST = [1, 2, 3, 4, 5, 6, 0]
HYBRID_LIST = ['SP','SP2','SP3','SP3D','SP3D2','S','UNSPECIFIED']
HYDROGEN_LIST = [1, 2, 3, 4, 0]
RINGSIZE_LIST = [3, 4, 5, 6, 7, 8]
BOND_LIST = ['SINGLE', 'DOUBLE', 'TRIPLE', 'AROMATIC']

# -----------------------
# RDKit Feature Factory - portable path
# -----------------------
def _get_fdef_path():
    """Return the BaseFeatures.fdef path in a portable way."""
    # 1. Bundled path (preferred)
    local_fdef = os.path.join(_THIS, 'rdkit_data', 'BaseFeatures.fdef')
    if os.path.isfile(local_fdef):
        return local_fdef
    
    # 2. _internal path (PyInstaller)
    internal_fdef = os.path.join(os.path.dirname(_THIS), 'rdkit_data', 'BaseFeatures.fdef')
    if os.path.isfile(internal_fdef):
        return internal_fdef
    
    # 3. Fallback: system RDKit path (dev only)
    system_fdef = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    if os.path.isfile(system_fdef):
        return system_fdef
    
    raise FileNotFoundError(
        f"BaseFeatures.fdef not found. Paths tried:\n"
        f"  - {local_fdef}\n"
        f"  - {internal_fdef}\n"
        f"  - {system_fdef}"
    )

chem_feature_factory = ChemicalFeatures.BuildFeatureFactory(_get_fdef_path())

READOUT_MODE = "proposed_set2set"
NODE_EMB_DIM = 256
READOUT_N_HID_DIM = 512
BATCH_SIZE = 128

# -----------------------
# Logs JSON SSE
# -----------------------
def _print_json(obj: dict) -> None:
    """Print JSON safely for Windows cp1252 console."""
    try:
        # Use ensure_ascii=True to avoid Unicode encoding issues on Windows
        print(json.dumps(obj, ensure_ascii=True), flush=True)
    except Exception:
        # Fallback: convert to ASCII-safe string
        try:
            safe_str = str(obj).encode('ascii', errors='replace').decode('ascii')
            print(safe_str, flush=True)
        except Exception:
            print('{"type":"error","message":"Output encoding error"}', flush=True)

def _sanitize_message(message: str) -> str:
    """Replace Unicode symbols with ASCII equivalents for Windows compatibility."""
    replacements = {
        '\u2713': '[OK]',   # ✓
        '\u2714': '[OK]',   # ✔
        '\u2192': '->',     # →
        '\u2190': '<-',     # ←
        '\u2022': '*',      # •
        '\u00e9': 'e',      # e-acute
        '\u00e8': 'e',      # e-grave
        '\u00ea': 'e',      # e-circumflex
        '\u00e0': 'a',      # a-grave
        '\u00e2': 'a',      # a-circumflex
        '\u00f9': 'u',      # u-grave
        '\u00fb': 'u',      # u-circumflex
        '\u00e7': 'c',      # c-cedilla
    }
    for unicode_char, ascii_equiv in replacements.items():
        message = message.replace(unicode_char, ascii_equiv)
    return message

def send_log(message: str) -> None:
    _print_json({"type": "log", "message": _sanitize_message(message)})

def send_progress(percent: int, message: str = "", completed: int = None, total: int = None) -> None:
    """
    Send JSON progress.
    - percent: 0-100
    - message: descriptive text
    - completed: number of processed molecules (optional)
    - total: total number of molecules (optional)
    """
    p = max(0, min(100, int(percent)))
    obj = {"type": "progress", "current": p, "total": 100, "percent": p, "message": message}
    if completed is not None:
        obj["completedMolecules"] = int(completed)
    if total is not None:
        obj["totalMolecules"] = int(total)
    _print_json(obj)

def send_error(message: str, **diagnostics) -> None:
    obj = {"type": "error", "message": message}
    if diagnostics:
        obj["diagnostics"] = diagnostics
    _print_json(obj)

def send_completed(success: bool, **payload) -> None:
    obj = {"type": "completed", "success": bool(success)}
    obj.update(payload)
    _print_json(obj)

# -----------------------
# Utils
# -----------------------
def ensure_dir(p: str) -> str:
    os.makedirs(p, exist_ok=True)
    return p

def sanitize_filename(name: str) -> str:
    name = (name or "").strip()
    name = re.sub(r"[^\w.\-]+", "_", name)
    return name or "mol"

def human_bytes(n: int) -> str:
    if n is None:
        return "?"
    units = ['B','KB','MB','GB','TB']
    i = 0
    x = float(n)
    while x >= 1024.0 and i < len(units)-1:
        i += 1
        x /= 1024.0
    return f"{x:.1f} {units[i]}"

# -----------------------
# Sample dataset
# -----------------------
def _sorted_h_groups(hx: Dict[int, List[int]]) -> List[Tuple[int, List[int]]]:
    """Return H groups sorted by heavy-atom index (key), and H sorted within each group."""
    return [(k, sorted(v)) for k, v in sorted(hx.items(), key=lambda kv: kv[0]) if v]

class GraphDataset_sample:
    """
    Inference dataset aligned with the original repo:
    - 1H: no ground-truth y; we only encode the mask (heavy atoms bearing H)
          and each H-group size (CH3/CH2/CH) in a deterministic order (increasing heavy-atom index).
    - 13C: same as training, mask on carbons.
    """
    def __init__(self, target, file_path):
        self.target = target
        self.file_path = file_path
        self.load()

    def load(self):
        [mol_dict] = np.load(self.file_path, allow_pickle=True)['data']

        self.n_node = np.asarray(mol_dict['n_node'], dtype=int)
        self.n_edge = np.asarray(mol_dict['n_edge'], dtype=int)
        self.node_attr = np.asarray(mol_dict['node_attr'])
        self.edge_attr = np.asarray(mol_dict['edge_attr'])
        self.src = np.asarray(mol_dict['src'], dtype=np.int64)
        self.dst = np.asarray(mol_dict['dst'], dtype=np.int64)

        self.shift = mol_dict['shift']
        self.mask  = np.asarray(mol_dict['mask'])
        self.smi   = np.asarray(mol_dict['smi'])
        self.h_x_connectivity = mol_dict.get('h_x_connectivity', [])

        self.n_csum = np.concatenate([[0], np.cumsum(self.n_node)])
        self.e_csum = np.concatenate([[0], np.cumsum(self.n_edge)])

    def __getitem__(self, idx):
        n0, n1 = int(self.n_csum[idx]), int(self.n_csum[idx+1])
        e0, e1 = int(self.e_csum[idx]), int(self.e_csum[idx+1])

        g = dgl_graph((self.src[e0:e1], self.dst[e0:e1]), num_nodes=int(self.n_node[idx]))
        g.ndata['node_attr'] = torch.from_numpy(self.node_attr[n0:n1]).float()
        if e1 - e0 > 0:
            g.edata['edge_attr'] = torch.from_numpy(self.edge_attr[e0:e1]).float()
        else:
            feat_dim = self.edge_attr.shape[1] if self.edge_attr.ndim == 2 else 8
            g.edata['edge_attr'] = torch.zeros((0, feat_dim), dtype=torch.float32)

        n_node = self.n_node[idx:idx+1].astype(int)

        shift = np.asarray(self.shift[n0:n1], dtype=float) if np.ndim(self.shift) > 0 else np.zeros((n_node[0],), dtype=float)
        mask  = self.mask[n0:n1].astype(bool)

        if self.target == '1H':
            shift = np.zeros((n_node[0],), dtype=float)
            hx = self.h_x_connectivity[0] if len(self.h_x_connectivity) > 0 else {}
            ordered_groups = _sorted_h_groups(hx)
            numHshifts = np.array([len(v) for _, v in ordered_groups], dtype=int)
            shift_test = np.zeros(int(numHshifts.sum()) if numHshifts.size else 0, dtype=float)
        else:
            numHshifts = np.zeros(n_node, dtype=float)
            shift_test = shift

        return g, n_node, numHshifts, shift_test, shift, mask

    def __len__(self):
        return int(self.n_node.shape[0])

# -----------------------
# Molecule attributes -> sparse graph
# -----------------------
def _DA(mol: Chem.Mol) -> Tuple[List[int], List[int]]:
    D_list, A_list = [], []
    for feat in chem_feature_factory.GetFeaturesForMol(mol):
        fam = feat.GetFamily()
        if fam == 'Donor':    D_list.append(feat.GetAtomIds()[0])
        elif fam == 'Acceptor': A_list.append(feat.GetAtomIds()[0])
    return D_list, A_list

def _chirality(atom: Chem.Atom) -> List[int]:
    if atom.HasProp('Chirality'):
        return [(atom.GetProp('Chirality') == 'Tet_CW'), (atom.GetProp('Chirality') == 'Tet_CCW')]
    return [0, 0]

def _stereochemistry(bond: Chem.Bond) -> List[int]:
    if bond.HasProp('Stereochemistry'):
        return [(bond.GetProp('Stereochemistry') == 'Bond_Cis'), (bond.GetProp('Stereochemistry') == 'Bond_Trans')]
    return [0, 0]

def add_mol_sparsified_graph(mol_dict: Dict[str, Any], mol: Chem.Mol) -> Dict[str, Any]:
    n_node = mol.GetNumAtoms()
    n_edge = mol.GetNumBonds() * 2
    D_list, A_list = _DA(mol)

    atom_fea1 = np.eye(len(ATOM_LIST), dtype=bool)[[ATOM_LIST.index(a.GetSymbol()) for a in mol.GetAtoms()]]
    atom_fea2 = np.eye(len(CHARGE_LIST), dtype=bool)[[CHARGE_LIST.index(a.GetFormalCharge()) for a in mol.GetAtoms()]][:, :-1]
    atom_fea3 = np.eye(len(DEGREE_LIST), dtype=bool)[[DEGREE_LIST.index(a.GetDegree()) for a in mol.GetAtoms()]][:, :-1]
    atom_fea4 = np.eye(len(HYBRID_LIST), dtype=bool)[[HYBRID_LIST.index(str(a.GetHybridization())) for a in mol.GetAtoms()]][:, :-2]
    atom_fea5 = np.eye(len(HYDROGEN_LIST), dtype=bool)[[HYDROGEN_LIST.index(a.GetTotalNumHs(includeNeighbors=True)) for a in mol.GetAtoms()]][:, :-1]
    atom_fea6 = np.eye(len(VALENCE_LIST), dtype=bool)[[VALENCE_LIST.index(a.GetTotalValence()) for a in mol.GetAtoms()]][:, :-1]
    atom_fea7 = np.array([[(j in D_list), (j in A_list)] for j in range(mol.GetNumAtoms())], dtype=bool)
    atom_fea8 = np.array([_chirality(a) for a in mol.GetAtoms()], dtype=bool)
    atom_fea9 = np.array([[a.IsInRingSize(s) for s in RINGSIZE_LIST] for a in mol.GetAtoms()], dtype=bool)
    atom_fea10 = np.array([[a.GetIsAromatic(), a.IsInRing()] for a in mol.GetAtoms()], dtype=bool)

    node_attr = np.concatenate(
        [atom_fea1, atom_fea2, atom_fea3, atom_fea4, atom_fea5, atom_fea6, atom_fea7, atom_fea8, atom_fea9, atom_fea10],
        axis=1
    )

    shift = np.array([float(atom.GetProp('shift')) if atom.HasProp('shift') else 0.0 for atom in mol.GetAtoms()])
    mask = np.array([bool(atom.GetBoolProp('mask')) if atom.HasProp('mask') else False for atom in mol.GetAtoms()])

    mol_dict['n_node'].append(n_node)
    mol_dict['n_edge'].append(n_edge)
    mol_dict['node_attr'].append(node_attr)
    mol_dict['shift'].append(shift)
    mol_dict['mask'].append(mask)
    mol_dict['smi'].append(Chem.MolToSmiles(mol))

    if n_edge > 0:
        bond_fea1 = np.eye(len(BOND_LIST), dtype=bool)[[BOND_LIST.index(str(b.GetBondType())) for b in mol.GetBonds()]]
        bond_fea2 = np.array([_stereochemistry(b) for b in mol.GetBonds()], dtype=bool)
        bond_fea3 = np.array([[b.IsInRing(), b.GetIsConjugated()] for b in mol.GetBonds()], dtype=bool)

        edge_attr = np.concatenate([bond_fea1, bond_fea2, bond_fea3], axis=1).astype(bool)
        edge_attr = np.vstack([edge_attr, edge_attr])

        bond_loc = np.array([[b.GetBeginAtomIdx(), b.GetEndAtomIdx()] for b in mol.GetBonds()], dtype=int)
        src = np.hstack([bond_loc[:, 0], bond_loc[:, 1]])
        dst = np.hstack([bond_loc[:, 1], bond_loc[:, 0]])

        mol_dict['edge_attr'].append(edge_attr)
        mol_dict['src'].append(src)
        mol_dict['dst'].append(dst)
    else:
        mol_dict['edge_attr'].append(np.zeros((0, len(BOND_LIST) + 2 + 2), dtype=bool))
        mol_dict['src'].append(np.zeros((0,), dtype=int))
        mol_dict['dst'].append(np.zeros((0,), dtype=int))

    return mol_dict

# -----------------------
# NPZ preparation (1H / 13C)
# -----------------------
def mol2SDF(mol: Chem.Mol, folder: Optional[str] = None, name: Optional[str] = None) -> str:
    if name is None:
        name = f"{random.randint(0, 999999):06d}"
    if folder is None:
        folder = os.path.join(os.getcwd(), "trash")
    ensure_dir(folder)
    save_path = os.path.join(folder, f"{name}.sdf")
    writer = rdMolFiles.SDWriter(save_path)
    writer.write(mol)
    writer.close()
    return save_path


def get_molecule_data(compound_path: str):
    docline_list: List[str] = []
    start = False
    name = os.path.basename(compound_path).rsplit(".", 1)[0]

    with open(compound_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            docline_list.append(line)

    atom_list: List[str] = []
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

    connectivity_list: List[List[str]] = []
    for idx, line in enumerate(docline_list):
        if idx >= start_line and "M  END" not in line:
            add_list = line.split()
            if len(add_list) == 3:
                part_1 = str(int(line[:3])); part_2 = str(int(line[3:6])); part_3 = str(int(line[6:9])); part_4 = str(int(line[9:12]))
                add_list = [part_1, part_2, part_3, part_4]
            if add_list and add_list[0] not in ["M", ">", ""]:
                connectivity_list.append(add_list)
        if "M  END" in line:
            break

    tmp_mol_path = compound_path[:-4] + ".mol"
    with open(tmp_mol_path, "w", encoding="utf-8", errors="ignore") as out:
        for line in docline_list:
            out.write(line)
            if "M  END" in line:
                break

    mol = Chem.MolFromMolFile(tmp_mol_path, removeHs=False, sanitize=True)
    return atom_list, connectivity_list, docline_list, name, mol

def get_x_h_connectivity(connectivity_list: List[List[str]], atom_list: List[str]) -> Dict[int, List[int]]:
    x_h: Dict[int, List[int]] = {}
    for entry in connectivity_list:
        try:
            a = int(entry[0]) - 1
            b = int(entry[1]) - 1
        except Exception:
            continue
        atom = atom_list[a] if 0 <= a < len(atom_list) else ""
        conn = atom_list[b] if 0 <= b < len(atom_list) else ""
        if (atom in ["O", "S", "N", "P", "C"]) and conn == "H":
            x_h.setdefault(a, []).append(b)
        if atom == "H" and (conn in ["O", "S", "N", "P", "C"]):
            x_h.setdefault(b, []).append(a)
    return x_h

def save_as_npy_for_1H(mol: Chem.Mol, output_dir: str) -> Tuple[str, Dict[int, List[int]]]:
    mol = Chem.AddHs(mol, addCoords=True)
    sdf_path = mol2SDF(mol, output_dir, name="tmp_1H")

    atom_list, connectivity_list, _, _, mol_read = get_molecule_data(sdf_path)
    x_h_conn = get_x_h_connectivity(connectivity_list, atom_list)

    mol_rd = Chem.RemoveHs(mol_read)
    Chem.SanitizeMol(mol_rd)

    si = Chem.FindPotentialStereo(mol_rd)
    for element in si:
        if str(element.type) == 'Atom_Tetrahedral' and str(element.specified) == 'Specified':
            mol_rd.GetAtomWithIdx(element.centeredOn).SetProp('Chirality', str(element.descriptor))
        elif str(element.type) == 'Bond_Double' and str(element.specified) == 'Specified':
            mol_rd.GetBondWithIdx(element.centeredOn).SetProp('Stereochemistry', str(element.descriptor))

    H_HEAVY_LIST = [k for k, v in x_h_conn.items() if v]
    for atom in mol_rd.GetAtoms():
        if atom.GetIdx() in H_HEAVY_LIST:
            atom.SetProp('shift', str(0.0)); atom.SetBoolProp('mask', True)
        else:
            atom.SetProp('shift', str(0.0)); atom.SetBoolProp('mask', False)

    mol_dict = {
        'n_node': [], 'n_edge': [], 'node_attr': [], 'edge_attr': [],
        'src': [], 'dst': [], 'shift': [], 'mask': [], 'smi': [], 'h_x_connectivity': []
    }
    mol_dict = add_mol_sparsified_graph(mol_dict, mol_rd)

    mol_dict['n_node'] = np.array(mol_dict['n_node']).astype(int)
    mol_dict['n_edge'] = np.array(mol_dict['n_edge']).astype(int)
    mol_dict['node_attr'] = np.vstack(mol_dict['node_attr']).astype(bool)
    mol_dict['edge_attr'] = np.vstack(mol_dict['edge_attr']).astype(bool) if mol_dict['edge_attr'] else np.zeros((0, 8), dtype=bool)
    mol_dict['src'] = np.hstack(mol_dict['src']).astype(int) if mol_dict['src'] else np.zeros((0,), dtype=int)
    mol_dict['dst'] = np.hstack(mol_dict['dst']).astype(int) if mol_dict['dst'] else np.zeros((0,), dtype=int)
    mol_dict['shift'] = np.array(mol_dict['shift'])
    mol_dict['mask'] = np.hstack(mol_dict['mask']).astype(bool)
    mol_dict['smi'] = np.array(mol_dict['smi'])
    mol_dict['h_x_connectivity'].append(x_h_conn)

    save_path = os.path.join(output_dir, "sample_1H.npz")
    np.savez_compressed(save_path, data=[mol_dict])
    return save_path, x_h_conn

def save_as_npy_for_13C(mol: Chem.Mol, output_dir: str) -> str:
    mol_rd = Chem.Mol(mol)
    Chem.SanitizeMol(mol_rd)

    si = Chem.FindPotentialStereo(mol_rd)
    for element in si:
        if str(element.type) == 'Atom_Tetrahedral' and str(element.specified) == 'Specified':
            mol_rd.GetAtomWithIdx(element.centeredOn).SetProp('Chirality', str(element.descriptor))
        elif str(element.type) == 'Bond_Double' and str(element.specified) == 'Specified':
            mol_rd.GetBondWithIdx(element.centeredOn).SetProp('Stereochemistry', str(element.descriptor))

    smiles = Chem.MolToSmiles(mol_rd)
    if '.' in smiles:
        frags = Chem.GetMolFrags(mol_rd, asMols=True, sanitizeFrags=True)
        mol_rd = max(frags, key=lambda m: m.GetNumAtoms())
        send_log(f"Multi-fragment molecule: using largest fragment: {Chem.MolToSmiles(mol_rd)}")

    for atom in mol_rd.GetAtoms():
        if atom.GetSymbol() == "C":
            atom.SetProp('shift', str(0.0)); atom.SetBoolProp('mask', True)
        else:
            atom.SetProp('shift', str(0.0)); atom.SetBoolProp('mask', False)

    mol_rd = Chem.RemoveHs(mol_rd)
    mol_dict = {
        'n_node': [], 'n_edge': [], 'node_attr': [], 'edge_attr': [],
        'src': [], 'dst': [], 'shift': [], 'mask': [], 'smi': [], 'h_c_connectivity': []
    }
    mol_dict = add_mol_sparsified_graph(mol_dict, mol_rd)

    mol_dict['n_node'] = np.array(mol_dict['n_node']).astype(int)
    mol_dict['n_edge'] = np.array(mol_dict['n_edge']).astype(int)
    mol_dict['node_attr'] = np.vstack(mol_dict['node_attr']).astype(bool)
    mol_dict['edge_attr'] = np.vstack(mol_dict['edge_attr']).astype(bool) if mol_dict['edge_attr'] else np.zeros((0, 8), dtype=bool)
    mol_dict['src'] = np.hstack(mol_dict['src']).astype(int) if mol_dict['src'] else np.zeros((0,), dtype=int)
    mol_dict['dst'] = np.hstack(mol_dict['dst']).astype(int) if mol_dict['dst'] else np.zeros((0,), dtype=int)
    mol_dict['shift'] = np.hstack(mol_dict['shift'])
    mol_dict['mask'] = np.hstack(mol_dict['mask']).astype(bool)
    mol_dict['smi'] = np.array(mol_dict['smi'])

    save_path = os.path.join(output_dir, "sample_13C.npz")
    np.savez_compressed(save_path, data=[mol_dict])
    return save_path

# -----------------------
# Model loading & inference
# -----------------------
def load_model(target: str, model_path: str, sample_npz_path: str):
    """Load the ML model used for NMR prediction."""
    send_log(f"Loading {target} model from {model_path}")

    [mol_dict] = np.load(sample_npz_path, allow_pickle=True)['data']
    node_dim = mol_dict['node_attr'].shape[1]
    edge_dim = mol_dict['edge_attr'].shape[1] if mol_dict['edge_attr'].shape[0] > 0 else (len(BOND_LIST) + 2 + 2)

    send_log(f"Dimensions {target}: node_dim={node_dim}, edge_dim={edge_dim}")

    try:
        net = nmr_mpnn_PROPOSED(
            node_dim=node_dim,
            edge_dim=edge_dim,
            readout_mode=READOUT_MODE,
            node_embedding_dim=NODE_EMB_DIM,
            readout_n_hidden_dim=READOUT_N_HID_DIM
        )
    except TypeError:
        net = nmr_mpnn_PROPOSED(
            node_dim, edge_dim, READOUT_MODE, NODE_EMB_DIM, READOUT_N_HID_DIM
        )

    try:
        state = torch.load(model_path, map_location=torch.device("cpu"))

        state_dict = None
        if isinstance(state, dict):
            # Some checkpoints wrap the state dict
            if isinstance(state.get('state_dict'), dict):
                state_dict = state['state_dict']
            else:
                state_dict = state
        elif hasattr(state, 'state_dict'):
            # Some checkpoints were saved as full nn.Module objects
            send_log(f"Checkpoint {target} loaded as module: {type(state).__name__}")
            state_dict = state.state_dict()

        if not isinstance(state_dict, dict):
            raise TypeError(f"Unexpected checkpoint for {target}: {type(state).__name__}")

        net.load_state_dict(state_dict)
        net.to("cpu")
        net.eval()
        send_log(f"Model {target} loaded")
    except Exception as e:
        send_error(f"Error loading model {target}", error=str(e))
        raise
    return net


def inference_model(net, sample_npz_path: str, target: str, train_y_mean: float, train_y_std: float) -> np.ndarray:
    from torch.utils.data import DataLoader
    send_log(f"Running inference {target}...")
    
    data_sample = GraphDataset_sample(target, sample_npz_path)
    test_loader = DataLoader(
        dataset=data_sample,
        batch_size=BATCH_SIZE,
        shuffle=False,
        collate_fn=collate_reaction_graphs
    )
    net.to("cpu"); net.eval()
    test_y_pred, _ = inference(net, test_loader, train_y_mean, train_y_std)
    send_log(f"Inference {target} OK: {len(test_y_pred)} values")
    return test_y_pred

def load_std_mean(target: str, graph_representation: str = "sparsified") -> Tuple[float, float]:
    """Load mean/std from the training dataset, with fallback if unavailable."""
    try:
        data = GraphDataset(target, graph_representation)
        from torch.utils.data import DataLoader
        train_loader = DataLoader(
            dataset=data,
            batch_size=BATCH_SIZE,
            shuffle=True,
            collate_fn=collate_reaction_graphs,
            drop_last=True
        )
        train_y = np.hstack([inst[-2][inst[-1]] for inst in iter(train_loader.dataset)])
        mean, std = np.mean(train_y), np.std(train_y)
        send_log(f"Normalization {target}: mean={mean:.2f}, std={std:.2f}")
        return mean, std
    except Exception as e:
        send_log(f"[WARN] Dataset not available for {target}, using defaults")
        return (2.5, 2.0) if target == "1H" else (70.0, 50.0)

    # -----------------------
    # Reconstruct shift lists & SDF
    # -----------------------
def create_shift_list_1H(
    mol_with_H: Chem.Mol,
    x_h_connectivity: Dict[int, List[int]],
    preds_1h: np.ndarray
) -> List[float]:
    """
    Assign ONE prediction per group (CH3/CH2/CH) to ALL H atoms in the group.
    Group order = increasing heavy-atom index (consistent mask -> prediction order).
    """
    n = mol_with_H.GetNumAtoms()
    out = [0.0] * n
    groups = _sorted_h_groups(x_h_connectivity)

    if len(preds_1h) != len(groups):
        send_log(f"[WARN] 1H mismatch: preds={len(preds_1h)} vs groups={len(groups)} "
                 f"(extra groups will receive 0.0).")

    for gi, (_, hlist) in enumerate(groups):
        val = float(preds_1h[gi]) if gi < len(preds_1h) else 0.0
        for hidx in hlist:
            if 0 <= hidx < n:
                out[hidx] = val
    return out

def create_shift_list_13C(mol_with_H: Chem.Mol, preds_13c: np.ndarray) -> List[float]:
    final_list: List[float] = []
    c_count = 0
    for atom in mol_with_H.GetAtoms():
        if atom.GetSymbol() == "C":
            val = float(preds_13c[c_count]) if c_count < len(preds_13c) else 0.0
            final_list.append(val)
            c_count += 1
        else:
            final_list.append(0.0)
    return final_list

def combine_1H_13C_lists(mol_with_H: Chem.Mol, list_1h: List[float], list_13c: List[float]) -> List[float]:
    """
        Final per-atom vector:
            - H -> 1H value
            - C -> 13C value
            - others -> 0.0
    """
    combined: List[float] = []
    for atom in mol_with_H.GetAtoms():
        i = atom.GetIdx()
        if atom.GetSymbol() == "H":
            combined.append(list_1h[i] if i < len(list_1h) else 0.0)
        elif atom.GetSymbol() == "C":
            combined.append(list_13c[i] if i < len(list_13c) else 0.0)
        else:
            combined.append(0.0)
    return combined

def generate_sdf_with_shifts(
    smiles: str,
    mol_original: Chem.Mol,
    shifts_1h: np.ndarray,
    shifts_13c: np.ndarray,
    x_h_connectivity: Dict[int, List[int]],
    output_dir: str,
    outfile_name: Optional[str] = None,
    sdf_property_key: str = "(1)"
) -> str:
    send_log('Generating enriched SDF...')
    
    molH = Chem.AddHs(mol_original, addCoords=True)

    try:
        if molH.GetNumConformers() == 0:
            AllChem.EmbedMolecule(molH, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(molH)
    except Exception:
        pass

    list_H = create_shift_list_1H(molH, x_h_connectivity, shifts_1h)
    list_C = create_shift_list_13C(molH, shifts_13c)
    all_shifts = combine_1H_13C_lists(molH, list_H, list_C)

    ts = int(time.time())
    if not outfile_name:
        outfile_name = f"nmr_prediction_{ts}.sdf"
    sdf_path = os.path.join(output_dir, outfile_name)

    writer = rdMolFiles.SDWriter(sdf_path)
    molH.SetProp('SMILES', smiles)
    molH.SetProp('Method', 'ML_Prediction_1H_13C_SDF_Only')
    molH.SetProp('Timestamp', str(ts))
    shifts_str = " ".join(f"{x:.3f}" for x in all_shifts)
    molH.SetProp(sdf_property_key, shifts_str)
    writer.write(molH)
    writer.close()

    # [OK] Detectable message for the frontend
    send_log(f"SDF written: {sdf_path} | {len(all_shifts)} positions")
    return sdf_path


# -----------------------
# Helpers spécifiques simulation
# -----------------------
def _build_x_h_connectivity_from_mol_with_H(molH: Chem.Mol) -> Dict[int, List[int]]:
    """Construit {heavy_idx: [H_idx,...]} depuis un Mol *avec* H."""
    xh: Dict[int, List[int]] = {}
    for a in molH.GetAtoms():
        if a.GetSymbol() != "H":
            continue
        hidx = a.GetIdx()
        nbs = [n.GetIdx() for n in a.GetNeighbors()]
        if not nbs:
            continue
        heavy = nbs[0]
        xh.setdefault(heavy, []).append(hidx)
    for k in xh:
        xh[k] = sorted(xh[k])
    return xh

def _proton_shifts_ppm_dict_from_list(molH: Chem.Mol, per_atom_H_ppm: List[float]) -> Dict[int, float]:
    """Retourne {H_idx: δ(ppm)}."""
    if len(per_atom_H_ppm) < molH.GetNumAtoms():
        per_atom_H_ppm = list(per_atom_H_ppm) + [0.0] * (molH.GetNumAtoms() - len(per_atom_H_ppm))
    out: Dict[int, float] = {}
    for a in molH.GetAtoms():
        if a.GetSymbol() == "H":
            out[a.GetIdx()] = float(per_atom_H_ppm[a.GetIdx()])
    return out

# -----------------------
# Multiplicités: sidecar + append props
# -----------------------
def _write_multiplicity_json(out_dir: str, mol_id: str, payload: dict) -> str:
    path = os.path.join(out_dir, f"{sanitize_filename(mol_id)}_multiplicity.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
    send_log(f"Multiplicités: JSON -> {os.path.basename(path)}")
    return path

def _append_props_to_sdf(sdf_path: str, props: Dict[str, str]) -> None:
    try:
        suppl = rdMolFiles.SDMolSupplier(sdf_path, removeHs=False, sanitize=False)
        mols = [m for m in suppl if m is not None]
        if not mols:
            return
        m = mols[0]
        for k, v in props.items():
            m.SetProp(k, v)
        w = rdMolFiles.SDWriter(sdf_path)
        w.write(m); w.close()
    except Exception as e:
        send_log(f"⚠️ Impossible d'ajouter des propriétés SDF: {e}")


def save_spectral_json(
    spectral_data: List[Dict[str, Any]],
    output_dir: str,
    filename: str = "spectral_data_1H.json",
) -> str:
    """
    Save the spectral data for all molecules to a JSON file.
    
    Args:
        spectral_data: List of spectral JSON dictionaries (one per molecule)
        output_dir: Output directory
        filename: Output filename
    
    Returns:
        Path to the saved JSON file
    """
    output_path = os.path.join(output_dir, filename)
    
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(spectral_data, f, indent=2, ensure_ascii=False)
    
    send_log(f"Spectral data saved: {output_path} ({len(spectral_data)} molecules)")
    return output_path


# -----------------------
# Per-molecule orchestration
# -----------------------
@dataclass
class ModelBundle:
    net_1h: Optional[Any]
    net_13c: Optional[Any]
    mean_std_1h: Optional[Tuple[float, float]]
    mean_std_13c: Optional[Tuple[float, float]]

def process_one_smiles(
    smiles: str,
    mol_id: str,
    out_dir: str,
    models: ModelBundle,
    model_1h_path: str,
    model_13c_path: str,
    property_key: str = "(1)",
    multiplicity_mod=None,
    sim_params: Optional[dict] = None,
    write_mult_json: bool = True,
    attach_labels_to_sdf: bool = True,
    return_payload: bool = False,  # NEW: return payload instead of writing file (for batch mode)
) -> Tuple[str, int, int, ModelBundle, Optional[str], Optional[dict]]:
    
    # ✅ CORRECTION : Sauvegarder le SMILES ORIGINAL avant toute transformation
    original_input_smiles = (smiles or "").strip()
    
    smiles = original_input_smiles
    if not smiles:
        raise ValueError("SMILES vide.")
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES invalide : {smiles}")
    
    # ✅ VALIDATION STRUCTURE (si simulateur disponible)
    if multiplicity_mod and hasattr(multiplicity_mod, '_validate_molecular_structure'):
        try:
            multiplicity_mod._validate_molecular_structure(mol, original_input_smiles)
            send_log(f"[{mol_id}] ✓ Structure moléculaire validée")
        except ValueError as ve:
            send_error(f"[{mol_id}] Structure moléculaire invalide", error=str(ve))
            raise

    work_dir = ensure_dir(os.path.join(out_dir, f"tmp_{sanitize_filename(mol_id)}"))

    if models.mean_std_1h is None:
        models.mean_std_1h = load_std_mean("1H")
    if models.mean_std_13c is None:
        models.mean_std_13c = load_std_mean("13C")

    train_y_mean_H, train_y_std_H = models.mean_std_1h
    train_y_mean_C, train_y_std_C = models.mean_std_13c

    # 1H
    send_log(f"[{mol_id}] Préparation 1H...")
    npz_1h, x_h_conn = save_as_npy_for_1H(mol, work_dir)
    try:
        [md] = np.load(npz_1h, allow_pickle=True)['data']
        nb_groups = len(_sorted_h_groups(md.get('h_x_connectivity', [{}])[0])) if md.get('h_x_connectivity', []) else 0
        send_log(f"[{mol_id}] NPZ 1H ok: n_node={md['n_node'].tolist()}, mask_true={int(np.sum(md['mask']))}, groupes_H={nb_groups}")
    except Exception as dbg:
        send_log(f"[{mol_id}] NPZ 1H debug fail: {dbg}")

    if models.net_1h is None:
        send_log(f"[{mol_id}] Chargement modèle 1H...")
        models.net_1h = load_model("1H", model_1h_path, npz_1h)
    send_log(f"[{mol_id}] Inférence 1H...")
    pred_1h = inference_model(models.net_1h, npz_1h, "1H", train_y_mean_H, train_y_std_H)

    expected_groups = len(_sorted_h_groups(x_h_conn))
    send_log(f"[{mol_id}] Groupes 1H attendus={expected_groups} | prédits={len(pred_1h)}")

    # 13C - avec gestion d'erreur pour les molécules problématiques
    pred_13c = np.array([])
    npz_13c = None
    try:
        send_log(f"[{mol_id}] Préparation 13C...")
        npz_13c = save_as_npy_for_13C(mol, work_dir)
        if models.net_13c is None:
            send_log(f"[{mol_id}] Chargement modèle 13C...")
            models.net_13c = load_model("13C", model_13c_path, npz_13c)
        send_log(f"[{mol_id}] Inférence 13C...")
        pred_13c = inference_model(models.net_13c, npz_13c, "13C", train_y_mean_C, train_y_std_C)
    except Exception as e:
        send_log(f"[{mol_id}] ⚠ Prédiction 13C échouée (kékulisation): {e}")
        send_log(f"[{mol_id}] Continuation avec seulement 1H...")
        pred_13c = np.zeros(sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'C'))

    # SDF final
    final_name = f"{sanitize_filename(mol_id)}_predicted.sdf"
    send_log(f"[{mol_id}] Génération SDF final...")
    sdf_path = generate_sdf_with_shifts(
        smiles=original_input_smiles,  # ✅ Utiliser SMILES original
        mol_original=mol,
        shifts_1h=pred_1h,
        shifts_13c=pred_13c,
        x_h_connectivity=x_h_conn,
        output_dir=out_dir,
        outfile_name=final_name,
        sdf_property_key=property_key
    )

    # --- Multiplicités 1H (optionnel) — VERSION ENRICHIE AVEC CONSTANTES J

    sidecar_path = None
    if multiplicity_mod is not None:
        try:
            # 1) Configuration simulation
            sp = sim_params or {}
            SimConfig = getattr(multiplicity_mod, "SimConfig", None)
            cfg = SimConfig() if SimConfig else None
            if cfg:
                cfg.separate_diastereotopic_protons = True
                cfg.spectrometer_mhz = float(sp.get("freq_mhz", 700.0))

            # 2) Reconstruire la molécule avec H pour les indices
            molH_for_sim = Chem.AddHs(mol, addCoords=True)
            
            # Générer géométrie 3D si demandé
            if sp.get("use_3d", True) and molH_for_sim.GetNumConformers() == 0:
                try:
                    AllChem.EmbedMolecule(molH_for_sim, randomSeed=42)
                    AllChem.MMFFOptimizeMolecule(molH_for_sim)
                except Exception:
                    pass
            
            x_h_conn_sim = _build_x_h_connectivity_from_mol_with_H(molH_for_sim)
            
            # 3) Créer le mapping des shifts ML PRÉDITS (pas les shifts heuristiques!)
            list_H = create_shift_list_1H(molH_for_sim, x_h_conn, pred_1h)
            proton_shifts_ppm = _proton_shifts_ppm_dict_from_list(molH_for_sim, list_H)

            # ✅ 4) Utiliser simulate_multiplicity_per_proton() avec les VRAIS shifts ML
            # Cela évite le bug de strong_coupling quand tous les shifts sont identiques
            per_proton = multiplicity_mod.simulate_multiplicity_per_proton(
                mol=mol,
                x_h_connectivity=x_h_conn_sim,
                proton_shifts_ppm=proton_shifts_ppm,  # ← shifts ML prédits !
                solvent=sp.get("solvent", "CDCl3"),
                spectrometer_mhz=float(sp.get("freq_mhz", 700.0)),
                temperature_k=float(sp.get("temperature_k", 298.15)),
                include_heteronuclear=True,
                use_3d=bool(sp.get("use_3d", True)),
                config=cfg,
            )

            # ✅ 5) REGROUPEMENT INTELLIGENT DES PROTONS
            # per_proton contient un représentant par groupe (y compris séparation vinylique)
            # On doit reconstruire les groupes en utilisant per_proton comme source de vérité
            
            # Mapper chaque H de per_proton vers son heavy_idx
            repr_to_heavy: Dict[int, int] = {}
            for h_idx in per_proton.keys():
                nbs = [n.GetIdx() for n in molH_for_sim.GetAtomWithIdx(h_idx).GetNeighbors()]
                if nbs:
                    repr_to_heavy[h_idx] = nbs[0]
            
            # Regrouper les protons par heavy_idx ET par pattern de multiplicité identique
            # Cela permet de séparer les protons vinyliques qui ont le même heavy mais des J différents
            heavy_to_groups: Dict[int, List[List[int]]] = {}
            for h_idx, heavy_idx in repr_to_heavy.items():
                data = per_proton[h_idx]
                # Créer une clé unique basée sur label + J_couplings
                J_key = tuple(sorted([(j.get("J_hz", 0), j.get("n_equiv", 1)) for j in data.get("J_couplings", [])]))
                label_key = (data.get("label", "s"), J_key)
                
                if heavy_idx not in heavy_to_groups:
                    heavy_to_groups[heavy_idx] = []
                
                # Chercher un groupe existant avec le même pattern
                found = False
                for group in heavy_to_groups[heavy_idx]:
                    repr_h = group[0]
                    repr_data = per_proton[repr_h]
                    repr_J_key = tuple(sorted([(j.get("J_hz", 0), j.get("n_equiv", 1)) for j in repr_data.get("J_couplings", [])]))
                    repr_label_key = (repr_data.get("label", "s"), repr_J_key)
                    if label_key == repr_label_key:
                        group.append(h_idx)
                        found = True
                        break
                
                if not found:
                    heavy_to_groups[heavy_idx].append([h_idx])
            
            # Ajouter les protons équivalents qui ne sont pas dans per_proton (CH₃, CH₂ équivalents)
            for heavy_idx, h_list in x_h_conn_sim.items():
                repr_h = None
                for h in h_list:
                    if h in per_proton:
                        repr_h = h
                        break
                if repr_h is not None:
                    # Ajouter les H manquants au groupe du représentant
                    for group in heavy_to_groups.get(heavy_idx, []):
                        if repr_h in group:
                            for h in h_list:
                                if h not in group and h not in per_proton:
                                    group.append(h)
                            break

            n_atoms = molH_for_sim.GetNumAtoms()
            mult_per_atom = [""] * n_atoms

            groups_items: List[dict] = []
            group_map_entries: List[str] = []
            group_id = 0

            # Base sans H pour détection OH/NH/COOH
            base_for_exch = Chem.RemoveHs(molH_for_sim)

            # ✅ BOUCLE CORRIGÉE : Itérer sur les groupes construits (heavy_to_groups)
            # Chaque heavy_idx peut avoir plusieurs groupes si les protons sont non-équivalents
            for heavy_idx in sorted(heavy_to_groups.keys()):
                for h_group in heavy_to_groups[heavy_idx]:
                    h_indices = sorted(h_group)
                    
                    if not h_indices:
                        continue
                    
                    # Le représentant est le premier H du groupe (il est dans per_proton)
                    h_idx_repr = h_indices[0]
                    data = per_proton.get(h_idx_repr, {})
                    
                    label_str = data.get("label", "s")
                    label_idx = int(data.get("label_idx", 0))

                    # Déterminer l'atome lourd
                    if 0 <= heavy_idx < base_for_exch.GetNumAtoms():
                        heavy_atom = base_for_exch.GetAtomWithIdx(heavy_idx)
                        heavy_symbol = heavy_atom.GetSymbol()
                    else:
                        heavy_atom = None
                        heavy_symbol = "?"

                    # ✅ Type de groupe : vérifier si c'est un proton séparé (vinylique, diastéréotopique)
                    note = "unique_or_equiv"
                    
                    # Si ce groupe contient moins de H que x_h_conn_sim, les protons ont été séparés
                    original_h_count = len(x_h_conn_sim.get(heavy_idx, []))
                    if len(h_indices) < original_h_count and original_h_count == 2:
                        # Vérifier si c'est vinylique (sp2) ou diastéréotopique (sp3)
                        if heavy_atom is not None:
                            hyb = heavy_atom.GetHybridization()
                            if hyb == rdchem.HybridizationType.SP2:
                                note = "vinyl_nonequiv"
                            else:
                                note = "diastereotopic"
                    
                    # Vérifier geminal partner dans les données
                    partners = data.get("partners", [])
                    for p in partners:
                        partner_idx = p.get("partner_h_index")
                        J_hz = p.get("J_hz", 0)
                        if partner_idx and partner_idx not in h_indices:
                            partner_nbs = [n.GetIdx() for n in molH_for_sim.GetAtomWithIdx(partner_idx).GetNeighbors()]
                            if partner_nbs and partner_nbs[0] == heavy_idx:
                                # Partenaire geminal sur le même heavy
                                if 0 < abs(J_hz) < 4:
                                    note = "vinyl_nonequiv"
                                elif 10 < abs(J_hz) < 15:
                                    note = "diastereotopic"
                                break

                    # Détection protons échangeables
                    is_exchangeable = False
                    exch_type = None
                    try:
                        if hasattr(multiplicity_mod, "_is_exchangeable_H_on_X") and heavy_atom is not None:
                            exch_type = multiplicity_mod._is_exchangeable_H_on_X(base_for_exch, heavy_idx)
                            if exch_type:
                                is_exchangeable = True
                                note = "exchangeable"
                    except Exception:
                        pass

                    # ✅ ASSIGNER LE LABEL À TOUS LES H DU GROUPE
                    for h in h_indices:
                        if 0 <= h < n_atoms:
                            mult_per_atom[h] = label_str

                    # ✅ Collecter les shifts de tous les H du groupe
                    h_shifts = [float(proton_shifts_ppm.get(h, 0.0)) for h in h_indices]

                    # Construction de l'entrée JSON
                    groups_items.append({
                        "group_id": group_id,
                        "heavy_idx": heavy_idx,
                        "heavy_symbol": heavy_symbol,
                        "H_indices": [int(h) for h in h_indices],
                        "H_shifts_ppm": h_shifts,
                        "label": label_str,
                        "label_index": label_idx,
                        "note": note,
                        "second_order_type": data.get("second_order_type", "first_order"),
                        "system_type": data.get("system_type"),
                        "virtual_coupling": bool(data.get("virtual_coupling", False)),
                        "J_couplings": data.get("J_couplings", []),
                        "partners": data.get("partners", []),
                        "satellites_13C_hz": data.get("satellites_13C_hz", []),
                    })

                    # ✅ Mapper TOUS les H de ce groupe au même group_id
                    for h in h_indices:
                        group_map_entries.append(f"{h}:{group_id}")
                    
                    group_id += 1

            # 6) JSON sidecar enrichi
            payload = None
            if write_mult_json or return_payload:
                payload = {
                    "mol_id": mol_id,
                    "smiles": original_input_smiles,
                    "solvent": sp.get("solvent", "CDCl3"),
                    "freq_MHz": float(sp.get("freq_mhz", 700.0)),
                    "temperature_K": float(sp.get("temperature_k", 298.15)),
                    "groups": groups_items,
                    "nb_groups": len(groups_items),
                    "nb_protons": sum(len(g["H_indices"]) for g in groups_items),  # ✅ Total des H
                    "note": "Protons grouped by heavy atom; diastereotopic protons separated only if truly different",
                    "typical_linewidth_hz": float(cfg.line_width_hz) if cfg else 0.8,
                    "tolerance_grouping_hz": float(getattr(multiplicity_mod, "TOLERANCE_J_HZ", 0.15)),
                    "digital_resolution_hz": float(cfg.digital_resolution_hz) if cfg and cfg.digital_resolution_hz else None
                }
                # Write individual file only if NOT in batch mode (return_payload=False)
                if write_mult_json and not return_payload:
                    sidecar_path = _write_multiplicity_json(out_dir, mol_id, payload)

            # 7) Propriétés SDF
            if attach_labels_to_sdf:
                props = {
                    "H1_MULT_LABELS": " ".join(mult_per_atom),
                    "H1_MULT_GROUPS": str(group_id),
                    "H1_MULT_MAP": " ".join(group_map_entries),
                }
                _append_props_to_sdf(sdf_path, props)

        except Exception as e:
            send_log(f"[{mol_id}] ⚠️ Simulation multiplicités échouée: {e}")
            send_log(traceback.format_exc())
            payload = None  # Ensure payload is None on error
    # Nettoyage
    try:
        shutil.rmtree(work_dir, ignore_errors=True)
    except Exception:
        pass

    nH = int(len(pred_1h))
    nC = int(len(pred_13c))
    return sdf_path, nH, nC, models, sidecar_path, payload

# -----------------------
# Modes single / CSV
# -----------------------
def run_single(
    smiles: str,
    mol_id: Optional[str],
    out_dir: str,
    model_1h_path: str,
    model_13c_path: str,
    property_key: str,
    multiplicity_mod=None,
    sim_params: Optional[dict] = None,
) -> None:
    ensure_dir(out_dir)

    models = ModelBundle(net_1h=None, net_13c=None, mean_std_1h=None, mean_std_13c=None)
    mol_id = mol_id or "mol"

    try:
        sdf_path, nH, nC, _, sidecar, _ = process_one_smiles(
            smiles=smiles,
            mol_id=mol_id,
            out_dir=out_dir,
            models=models,
            model_1h_path=model_1h_path,
            model_13c_path=model_13c_path,
            property_key=property_key,
            multiplicity_mod=multiplicity_mod,
            sim_params=sim_params,
            return_payload=False,  # Single mode: write individual file
        )
        
        send_completed(True, sdf=sdf_path, mol_id=mol_id, nb_predicted_H=nH, nb_predicted_C=nC,
                       multiplicity_json=sidecar or "")
    except Exception as e:
        tb = traceback.format_exc()
        send_error(f"Prediction failed for '{mol_id}'", error=str(e), traceback=tb)
        send_completed(False, mol_id=mol_id, error=str(e))

def run_csv(
    csv_path: str,
    out_dir: str,
    smiles_col: str,
    id_col: str,
    model_1h_path: str,
    model_13c_path: str,
    property_key: str,
    multiplicity_mod=None,
    sim_params: Optional[dict] = None,
) -> None:
    ensure_dir(out_dir)

    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames or []
        if smiles_col not in headers or id_col not in headers:
            raise ValueError(f"CSV must contain columns '{id_col}' and '{smiles_col}' (found: {headers})")

    total = sum(1 for _ in open(csv_path, encoding="utf-8")) - 1
    
    # [OK] Initialize the global counter
    _GLOBAL_MOLECULE_COUNTER["total"] = total
    _GLOBAL_MOLECULE_COUNTER["done"] = 0
    
    done = 0

    models = ModelBundle(net_1h=None, net_13c=None, mean_std_1h=None, mean_std_13c=None)
    manifest_rows: List[dict] = []
    all_multiplicity_payloads: List[dict] = []  # NEW: collect all payloads for consolidated JSON
    ok = 0
    ko = 0

    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            done += 1
            mol_id = (row.get(id_col) or f"mol_{done}").strip()
            smi = (row.get(smiles_col) or "").strip()

            try:
                # NEW: use return_payload=True to get payload without writing individual file
                sdf_path, nH, nC, models, _, payload = process_one_smiles(
                    smiles=smi, mol_id=mol_id, out_dir=out_dir, models=models,
                    model_1h_path=model_1h_path, model_13c_path=model_13c_path,
                    property_key=property_key,
                    multiplicity_mod=multiplicity_mod,
                    sim_params=sim_params,
                    return_payload=True,  # Don't write individual JSON, return payload
                )
                ok += 1
                
                # Collect payload for consolidated JSON
                if payload is not None:
                    all_multiplicity_payloads.append(payload)
                
                # [OK] Increment and emit progress
                _GLOBAL_MOLECULE_COUNTER["done"] += 1
                send_log(f"[OK] Molecule {_GLOBAL_MOLECULE_COUNTER['done']}/{_GLOBAL_MOLECULE_COUNTER['total']} processed: {os.path.basename(sdf_path)}")
                
                progress_pct = (_GLOBAL_MOLECULE_COUNTER['done'] / _GLOBAL_MOLECULE_COUNTER['total']) * 98
                send_progress(
                    progress_pct, 
                    f"Molecule {_GLOBAL_MOLECULE_COUNTER['done']}/{_GLOBAL_MOLECULE_COUNTER['total']} (SDF)",
                    completed=_GLOBAL_MOLECULE_COUNTER['done'],
                    total=_GLOBAL_MOLECULE_COUNTER['total']
                )
                
                send_log(f"[{mol_id}] SDF : {os.path.basename(sdf_path)}")
                manifest_rows.append({
                    "mol_id": mol_id,
                    "smiles": smi,
                    "sdf": sdf_path,
                    "nb_predicted_H": nH,
                    "nb_predicted_C": nC,
                    "error": ""
                })
            except Exception as e:
                ko += 1
                send_log(f"[{mol_id}] Error: {e}")
                manifest_rows.append({
                    "mol_id": mol_id,
                    "smiles": smi,
                    "sdf": "",
                    "nb_predicted_H": 0,
                    "nb_predicted_C": 0,
                    "error": str(e)
                })

    # Write consolidated multiplicity JSON (ONE file for ALL molecules)
    consolidated_multiplicity_path = os.path.join(out_dir, "all_multiplicities.json")
    with open(consolidated_multiplicity_path, "w", encoding="utf-8") as mf:
        json.dump({
            "total_molecules": len(all_multiplicity_payloads),
            "ok": ok,
            "ko": ko,
            "timestamp": time.time(),
            "molecules": all_multiplicity_payloads
        }, mf, indent=2, ensure_ascii=False)
    send_log(f"Multiplicités consolidées: {consolidated_multiplicity_path} ({len(all_multiplicity_payloads)} molécules)")

    # Write manifest (summary)
    manifest_path = os.path.join(out_dir, "sdf_manifest.json")
    with open(manifest_path, "w", encoding="utf-8") as mf:
        json.dump({
            "ok": ok,
            "ko": ko,
            "total": total,
            "output_dir": out_dir,
            "multiplicity_json": consolidated_multiplicity_path,  # Reference to consolidated file
            "items": manifest_rows,
            "timestamp": time.time()
        }, mf, indent=2, ensure_ascii=False)

    send_log(f"CSV finished. Success: {ok} | Failures: {ko} | Manifest: {manifest_path}")
    send_completed(True, ok=ok, ko=ko, output_dir=out_dir, manifest=manifest_path, 
                   multiplicity_json=consolidated_multiplicity_path)

# -----------------------
# Simulateur de multiplicités: import dynamique
# -----------------------
def _load_multiplicity_sim(sim_path: str):
    p = Path(sim_path).resolve()
    if not p.exists():
        raise FileNotFoundError(f"Simulator not found: {p}")

    mod_name = "sim_mult_rulebased"

    spec = importlib.util.spec_from_file_location(mod_name, str(p))
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not create spec for {p}")

    mod = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = mod
    spec.loader.exec_module(mod)
    return mod


# -----------------------
# Argparser & main
# -----------------------
def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Génération SDF (déplacements 1H/13C prédits) à partir d'un SMILES ou d'un CSV.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--smiles", help="SMILES à prédire (mode single)")
    src.add_argument("--csv", help="Fichier CSV avec colonnes 'id' et 'smiles' (ou --id_column/--smiles_column)")

    ap.add_argument("--output_dir", required=True, help="Dossier de sortie (SDF générés)")
    ap.add_argument("--model_1h_path", required=True, help="Chemin du modèle 1H (.pt)")
    ap.add_argument("--model_13c_path", required=True, help="Chemin du modèle 13C (.pt)")

    ap.add_argument("--mol_id", default="mol", help="Identifiant de la molécule quand --smiles est utilisé")
    ap.add_argument("--id_column", default="id", help="Nom de la colonne ID (CSV)")
    ap.add_argument("--smiles_column", default="smiles", help="Nom de la colonne SMILES (CSV)")

    ap.add_argument("--property_key", default="(1)", help="Clé de propriété SDF pour stocker la liste de shifts")
    ap.add_argument("--task_id", default=None, help="ID de tâche externe (optionnel, pour logs)")

    # --- Multiplicités 1H (rule-based)
    ap.add_argument("--simulate_multiplicity", action="store_true",
                    help="Active la simulation de multiplicités 1H.")
    ap.add_argument("--simulator_path", default=None,
                    help="Chemin vers sim_mult_rulebased.py (optionnel).")
    ap.add_argument("--sim_solvent", default="CDCl3",
                    help="Solvant pour la simulation (ex: CDCl3, DMSO-d6, C6D6...).")
    ap.add_argument("--sim_freq_mhz", type=float, default=700.0,
                    help="Fréquence spectromètre (MHz).")
    ap.add_argument("--sim_temperature_k", type=float, default=298.15,
                    help="Température (K).")
    ap.add_argument("--sim_use_3d", action="store_true",
                    help="Tente une géométrie 3D (ETKDG) pour Karplus/vinyliques.")
    ap.add_argument("--sim_debug", action="store_true",
                    help="Enregistre les détails (J/groupes/second ordre) dans le JSON.")

    return ap

def validate_models(paths: List[str]) -> None:
    for p in paths:
        if not os.path.exists(p):
            raise FileNotFoundError(f"Modèle introuvable : {p}")
        size = os.path.getsize(p)
        send_log(f"Modèle détecté : {p} ({human_bytes(size)})")

def main(argv: Optional[List[str]] = None) -> None:
    parser = build_argparser()
    args, unknown = parser.parse_known_args(argv)

    if unknown:
        send_log(f"Arguments CLI inconnus ignorés: {unknown}")
    if getattr(args, "task_id", None):
        send_log(f"Task ID: {args.task_id}")

    try:
        validate_models([args.model_1h_path, args.model_13c_path])
    except Exception as e:
        send_error("Validation des modèles échouée", error=str(e))
        sys.exit(1)

    # --- Chargement simulateur (si demandé)
    multiplicity_mod = None
    sim_params = None
    if getattr(args, "simulate_multiplicity", False):
        if args.simulator_path:
            try:
                multiplicity_mod = _load_multiplicity_sim(args.simulator_path)
                send_log(f"✓ Simulateur de multiplicités chargé depuis {args.simulator_path}")
            except Exception as e:
                send_log(f"⚠️ Échec chargement simulateur: {e}")
        else:
            send_log("⚠️ --simulate_multiplicity activé mais --simulator_path non fourni")
        
        if multiplicity_mod:
            sim_params = {
                "solvent": args.sim_solvent,
                "freq_mhz": args.sim_freq_mhz,
                "temperature_k": args.sim_temperature_k,
                "use_3d": bool(args.sim_use_3d),
                "debug": bool(args.sim_debug),
            }

    try:
        if args.smiles:
            run_single(
                smiles=args.smiles,
                mol_id=args.mol_id,
                out_dir=args.output_dir,
                model_1h_path=args.model_1h_path,
                model_13c_path=args.model_13c_path,
                property_key=args.property_key,
                multiplicity_mod=multiplicity_mod,
                sim_params=sim_params
            )
        else:
            run_csv(
                csv_path=args.csv,
                out_dir=args.output_dir,
                smiles_col=args.smiles_column,
                id_col=args.id_column,
                model_1h_path=args.model_1h_path,
                model_13c_path=args.model_13c_path,
                property_key=args.property_key,
                multiplicity_mod=multiplicity_mod,
                sim_params=sim_params
            )
    except Exception as e:
        send_error("Erreur fatale", error=str(e), traceback=traceback.format_exc())
        send_completed(False)
        sys.exit(1)

if __name__ == "__main__":
    main()