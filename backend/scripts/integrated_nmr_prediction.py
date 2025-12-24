#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""integrated_nmr_prediction.py (server, DGL compatible)
=======================================================
- Generates SDF files enriched with predicted 1H/13C shifts.
- CPU-only, Windows compatible.
- Aligns the 1H assignment "per H" with the extractor scripts (predict_1H_real.py / COSY / HSQC):
    * one value per H atom (repeat the CH3/CH2/CH group value)
    * 1H group order = sorted by heavy-atom index (consistent mask -> prediction order)

Usage (single):
    python integrated_nmr_prediction.py \
        --smiles "CCO" \
        --output_dir ./out \
        --model_1h_path ../models/1H_sparsified_proposed_proposed_set2set_1.pt \
        --model_13c_path ../models/13C_sparsified_proposed_proposed_set2set_1.pt \
        --mol_id ethanol

Usage (CSV):
    python integrated_nmr_prediction.py \
        --csv ./mols.csv \
        --output_dir ./out \
        --model_1h_path ../models/1H_sparsified_proposed_proposed_set2set_1.pt \
        --model_13c_path ../models/13C_sparsified_proposed_proposed_set2set_1.pt \
        --id_column sample-id --smiles_column SMILES
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
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Any

# [OK] Global counter for tracking
_GLOBAL_MOLECULE_COUNTER = {"done": 0, "total": 0}

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
from rdkit.Chem import AllChem, ChemicalFeatures
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
    property_key: str = "(1)"
) -> Tuple[str, int, int, ModelBundle]:
    smiles = (smiles or "").strip()
    if not smiles:
        raise ValueError('Empty SMILES.')
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    work_dir = ensure_dir(os.path.join(out_dir, f"tmp_{sanitize_filename(mol_id)}"))

    if models.mean_std_1h is None:
        models.mean_std_1h = load_std_mean("1H")
    if models.mean_std_13c is None:
        models.mean_std_13c = load_std_mean("13C")

    train_y_mean_H, train_y_std_H = models.mean_std_1h
    train_y_mean_C, train_y_std_C = models.mean_std_13c

    # 1H
    send_log(f"[{mol_id}] Preparing 1H...")
    npz_1h, x_h_conn = save_as_npy_for_1H(mol, work_dir)
    
    try:
        [md] = np.load(npz_1h, allow_pickle=True)['data']
        nb_groups = len(_sorted_h_groups(md.get('h_x_connectivity', [{}])[0])) if md.get('h_x_connectivity', []) else 0
        send_log(f"[{mol_id}] NPZ 1H ok: n_node={md['n_node'].tolist()}, mask_true={int(np.sum(md['mask']))}, h_groups={nb_groups}")
    except Exception as dbg:
        send_log(f"[{mol_id}] NPZ 1H debug fail: {dbg}")

    if models.net_1h is None:
        send_log(f"[{mol_id}] Loading 1H model...")
        models.net_1h = load_model("1H", model_1h_path, npz_1h)
    send_log(f"[{mol_id}] Running 1H inference...")
    pred_1h = inference_model(models.net_1h, npz_1h, "1H", train_y_mean_H, train_y_std_H)

    expected_groups = len(_sorted_h_groups(x_h_conn))
    send_log(f"[{mol_id}] Expected 1H groups={expected_groups} | predicted={len(pred_1h)}")

    # 13C
    send_log(f"[{mol_id}] Preparing 13C...")
    npz_13c = save_as_npy_for_13C(mol, work_dir)
    if models.net_13c is None:
        send_log(f"[{mol_id}] Loading 13C model...")
        models.net_13c = load_model("13C", model_13c_path, npz_13c)
    send_log(f"[{mol_id}] Running 13C inference...")
    pred_13c = inference_model(models.net_13c, npz_13c, "13C", train_y_mean_C, train_y_std_C)

    # SDF final
    final_name = f"{sanitize_filename(mol_id)}_predicted.sdf"
    send_log(f"[{mol_id}] Generating final SDF...")
    sdf_path = generate_sdf_with_shifts(
        smiles=smiles,
        mol_original=mol,
        shifts_1h=pred_1h,
        shifts_13c=pred_13c,
        x_h_connectivity=x_h_conn,
        output_dir=out_dir,
        outfile_name=final_name,
        sdf_property_key=property_key
    )

    # [OK] Detectable completion message
    send_log(f"[OK] Molecule {mol_id} : SDF generated successfully -> {os.path.basename(sdf_path)}")

    try:
        shutil.rmtree(work_dir, ignore_errors=True)
    except Exception:
        pass

    nH = int(len(pred_1h))
    nC = int(len(pred_13c))
    return sdf_path, nH, nC, models

# -----------------------
# Modes single / CSV
# -----------------------
def run_single(
    smiles: str,
    mol_id: Optional[str],
    out_dir: str,
    model_1h_path: str,
    model_13c_path: str,
    property_key: str
) -> None:
    ensure_dir(out_dir)

    models = ModelBundle(net_1h=None, net_13c=None, mean_std_1h=None, mean_std_13c=None)
    mol_id = mol_id or "mol"

    try:
        sdf_path, nH, nC, _ = process_one_smiles(
            smiles=smiles,
            mol_id=mol_id,
            out_dir=out_dir,
            models=models,
            model_1h_path=model_1h_path,
            model_13c_path=model_13c_path,
            property_key=property_key
        )
        send_completed(True, sdf=sdf_path, mol_id=mol_id, nb_predicted_H=nH, nb_predicted_C=nC)
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
    property_key: str
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
    ok = 0
    ko = 0

    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            done += 1
            mol_id = (row.get(id_col) or f"mol_{done}").strip()
            smi = (row.get(smiles_col) or "").strip()

            try:
                sdf_path, nH, nC, models = process_one_smiles(
                    smiles=smi, mol_id=mol_id, out_dir=out_dir, models=models,
                    model_1h_path=model_1h_path, model_13c_path=model_13c_path,
                    property_key=property_key
                )
                ok += 1
                
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

    manifest_path = os.path.join(out_dir, "sdf_manifest.json")
    with open(manifest_path, "w", encoding="utf-8") as mf:
        json.dump({
            "ok": ok,
            "ko": ko,
            "total": total,
            "output_dir": out_dir,
            "items": manifest_rows,
            "timestamp": time.time()
        }, mf, indent=2, ensure_ascii=False)

    send_log(f"CSV finished. Success: {ok} | Failures: {ko} | Manifest: {manifest_path}")
    send_completed(True, ok=ok, ko=ko, output_dir=out_dir, manifest=manifest_path)

# -----------------------
# Argparser & main
# -----------------------
def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Generate SDF (predicted 1H/13C shifts) from a SMILES or a CSV.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--smiles", help="SMILES to predict (single mode)")
    src.add_argument("--csv", help="CSV file with 'id' and 'smiles' columns (or --id_column/--smiles_column)")

    ap.add_argument("--output_dir", required=True, help="Output folder (generated SDFs)")
    ap.add_argument("--model_1h_path", required=True, help="Path to 1H model (.pt)")
    ap.add_argument("--model_13c_path", required=True, help="Path to 13C model (.pt)")

    ap.add_argument("--mol_id", default="mol", help="Molecule identifier when using --smiles")
    ap.add_argument("--id_column", default="id", help="ID column name (CSV)")
    ap.add_argument("--smiles_column", default="smiles", help="SMILES column name (CSV)")

    ap.add_argument("--property_key", default="(1)", help="SDF property key used to store the shift list")
    ap.add_argument("--task_id", default=None, help="External task ID (optional, for logs)")

    return ap

def validate_models(paths: List[str]) -> None:
    for p in paths:
        if not os.path.exists(p):
            raise FileNotFoundError(f"Model not found: {p}")
        size = os.path.getsize(p)
        send_log(f"Model detected: {p} ({human_bytes(size)})")

def main(argv: Optional[List[str]] = None) -> None:
    parser = build_argparser()
    args, unknown = parser.parse_known_args(argv)

    if unknown:
        send_log(f"Unknown CLI arguments ignored: {unknown}")
    if getattr(args, "task_id", None):
        send_log(f"Task ID: {args.task_id}")

    try:
        validate_models([args.model_1h_path, args.model_13c_path])
    except Exception as e:
        send_error('Model validation failed', error=str(e))
        sys.exit(1)

    try:
        if args.smiles:
            run_single(
                smiles=args.smiles,
                mol_id=args.mol_id,
                out_dir=args.output_dir,
                model_1h_path=args.model_1h_path,
                model_13c_path=args.model_13c_path,
                property_key=args.property_key
            )
        else:
            run_csv(
                csv_path=args.csv,
                out_dir=args.output_dir,
                smiles_col=args.smiles_column,
                id_col=args.id_column,
                model_1h_path=args.model_1h_path,
                model_13c_path=args.model_13c_path,
                property_key=args.property_key
            )
    except Exception as e:
        send_error('Fatal error', error=str(e), traceback=traceback.format_exc())
        send_completed(False)
        sys.exit(1)

if __name__ == "__main__":
    main()