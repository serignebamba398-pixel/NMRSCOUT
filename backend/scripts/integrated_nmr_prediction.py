#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Integrated NMR prediction pipeline and supporting utilities."""

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
import threading
import importlib.util
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Any
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

_GLOBAL_MOLECULE_COUNTER = {"done": 0, "total": 0}
_PRINT_LOCK = threading.Lock()
_MODEL_LOCK = threading.Lock()

_MOLECULES_SPECTRAL_DATA: List[Dict[str, Any]] = []

import numpy as np
import torch

os.environ.setdefault("CUDA_VISIBLE_DEVICES", "-1")
os.environ.setdefault("DGLBACKEND", "pytorch")
torch.backends.cudnn.enabled = False

# RDKit
from rdkit import Chem, RDConfig, rdBase
from rdkit.Chem import AllChem, ChemicalFeatures, rdchem
try:
    from rdkit.Chem import rdmolfiles as rdMolFiles
except ImportError:
    import rdkit.Chem.rdmolfiles as rdMolFiles

rdBase.DisableLog('rdApp.error')
rdBase.DisableLog('rdApp.warning')

AROMATIC_13C_CORRECTIONS_AVAILABLE = False

# DGL compatibility
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

# Local machine-learning dependencies
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

# Multiplicity simulation
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

# Fingerprint database predictor
try:
    from fingerprint_db_predictor import FingerprintPredictor
    _FINGERPRINT_PREDICTOR: Optional["FingerprintPredictor"] = None
    FINGERPRINT_DB_AVAILABLE = True
    print("[OK] Fingerprint database predictor imported", flush=True)
except ImportError as e:
    FINGERPRINT_DB_AVAILABLE = False
    _FINGERPRINT_PREDICTOR = None
    print(f"[WARN] Fingerprint database predictor not available: {e}", file=sys.stderr)

def get_fingerprint_predictor() -> Optional["FingerprintPredictor"]:
    """Return the shared fingerprint predictor instance."""
    global _FINGERPRINT_PREDICTOR
    if not FINGERPRINT_DB_AVAILABLE:
        return None
    if _FINGERPRINT_PREDICTOR is None:
        try:
            _FINGERPRINT_PREDICTOR = FingerprintPredictor()

            r6_count = len(_FINGERPRINT_PREDICTOR._db.get(6, {})) if _FINGERPRINT_PREDICTOR._db else 0
            r5_count = len(_FINGERPRINT_PREDICTOR._db.get(5, {})) if _FINGERPRINT_PREDICTOR._db else 0
            r4_count = len(_FINGERPRINT_PREDICTOR._db.get(4, {})) if _FINGERPRINT_PREDICTOR._db else 0
            r3_count = len(_FINGERPRINT_PREDICTOR._db.get(3, {})) if _FINGERPRINT_PREDICTOR._db else 0
            r2_count = len(_FINGERPRINT_PREDICTOR._db.get(2, {})) if _FINGERPRINT_PREDICTOR._db else 0
            print(f"[OK] Fingerprint predictor initialized (r6={r6_count}, r5={r5_count}, r4={r4_count}, r3={r3_count}, r2={r2_count})", flush=True)
        except Exception as e:
            print(f"[WARN] Failed to initialize fingerprint predictor: {e}", file=sys.stderr)
            return None
    return _FINGERPRINT_PREDICTOR

try:
    from diastereotopic_delta_table import (
        get_solvent_multiplier,
        is_aromatic_solvent,
        DELTA_CYCLOHEXANE,
        DELTA_CHIRAL_ADJ,
        DELTA_DEFAULT,
        SOLVENT_MULTIPLIER,
        AMINO_ACID_BETA_DELTA,
        ANISOTROPIC_SYSTEMS,
        CYCLIC_DELTA,
    )
    DIASTEREOTOPIC_TABLE_AVAILABLE = True
    print("[OK] Diastereotopic delta table imported", flush=True)
except ImportError as e:
    DIASTEREOTOPIC_TABLE_AVAILABLE = False
    print(f"[WARN] Diastereotopic delta table not available: {e}", file=sys.stderr)

    DELTA_CYCLOHEXANE = 0.48
    DELTA_CHIRAL_ADJ = 0.12
    DELTA_DEFAULT = 0.10
    SOLVENT_MULTIPLIER = {"CDCl3": 1.0, "DMSO-d6": 0.95, "C6D6": 1.6, "default": 1.0}
    AMINO_ACID_BETA_DELTA = {"default": (0.04, 0.01, "fallback")}
    ANISOTROPIC_SYSTEMS = {"calixarene_cone": (1.43, 1.70, "fallback")}
    CYCLIC_DELTA = {"cyclohexane": (0.48, 0.50, "fallback")}

    def get_solvent_multiplier(s): return SOLVENT_MULTIPLIER.get(s, 1.0)
    def is_aromatic_solvent(s): return s in ("C6D6", "pyridine-d5")

# Prediction constants
DELTA_GEMINAL_CH3_DOUBLE_BOND = 0.06

# Graph feature constants
ATOM_LIST = ['H','Li','B','C','N','O','F','Na','Mg','Al','Si','P','S','Cl','K','Ti','Zn','Ge','As','Se','Br','Pd','Ag','Sn','Sb','Te','I','Hg','Tl','Pb','Bi','Ga']
CHARGE_LIST = [1, 2, 3, -1, -2, -3, 0]
DEGREE_LIST = [1, 2, 3, 4, 5, 6, 0]
VALENCE_LIST = [1, 2, 3, 4, 5, 6, 0]
HYBRID_LIST = ['SP','SP2','SP3','SP3D','SP3D2','S','UNSPECIFIED']
HYDROGEN_LIST = [1, 2, 3, 4, 0]
RINGSIZE_LIST = [3, 4, 5, 6, 7, 8]
BOND_LIST = ['SINGLE', 'DOUBLE', 'TRIPLE', 'AROMATIC']

def _get_fdef_path():
    """Return a portable path to BaseFeatures.fdef."""

    local_fdef = os.path.join(_THIS, 'rdkit_data', 'BaseFeatures.fdef')
    if os.path.isfile(local_fdef):
        return local_fdef

    internal_fdef = os.path.join(os.path.dirname(_THIS), 'rdkit_data', 'BaseFeatures.fdef')
    if os.path.isfile(internal_fdef):
        return internal_fdef

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

# JSON logging
def _print_json(obj: dict) -> None:
    """Print JSON safely and serially."""
    with _PRINT_LOCK:
        try:

            print(json.dumps(obj, ensure_ascii=True), flush=True)
        except Exception:

            try:
                safe_str = str(obj).encode('ascii', errors='replace').decode('ascii')
                print(safe_str, flush=True)
            except Exception:
                print('{"type":"error","message":"Output encoding error"}', flush=True)

def _sanitize_message(message: str) -> str:
    """Convert selected Unicode symbols to ASCII-safe text."""
    replacements = {
        '\u2713': '[OK]',
        '\u2714': '[OK]',
        '\u2192': '->',
        '\u2190': '<-',
        '\u2022': '*',
        '\u00e9': 'e',
        '\u00e8': 'e',
        '\u00ea': 'e',
        '\u00e0': 'a',
        '\u00e2': 'a',
        '\u00f9': 'u',
        '\u00fb': 'u',
        '\u00e7': 'c',
    }
    for unicode_char, ascii_equiv in replacements.items():
        message = message.replace(unicode_char, ascii_equiv)
    return message

def send_log(message: str) -> None:
    _print_json({"type": "log", "message": _sanitize_message(message)})

def send_progress(percent: int, message: str = "", completed: int = None, total: int = None) -> None:
    """Emit a JSON progress event."""

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

# General utilities
def ensure_dir(p: str) -> str:
    os.makedirs(p, exist_ok=True)
    return p

def sanitize_filename(name: str, max_len: int = 80) -> str:
    import unicodedata, hashlib
    name = (name or "").strip()

    name = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode("ascii")
    name = re.sub(r"[^\w.\-]+", "_", name)
    if not name:
        return "mol"
    if len(name) > max_len:
        h = hashlib.md5(name.encode()).hexdigest()[:8]
        name = name[:max_len] + "_" + h
    return name

# Incremental resume utilities
def _load_existing_predictions(resume_dir: str) -> Tuple[Dict[str, dict], Dict[str, Any]]:
    """Load reusable predictions from a previous output directory."""

    mol_data = {}
    params = {}

    checkpoint_path = os.path.join(resume_dir, "checkpoint.json")
    if os.path.exists(checkpoint_path):
        try:
            with open(checkpoint_path, "r", encoding="utf-8") as f:
                checkpoint = json.load(f)
            opts = checkpoint.get("options", {})
            params = {
                "freq_mhz": opts.get("freq_mhz", opts.get("sim_freq_mhz", 700.0)),
                "solvent": opts.get("solvent", opts.get("sim_solvent", "CDCl3")),
                "temperature_k": opts.get("temperature_k", opts.get("sim_temperature_k", 298.15)),
                "ph": opts.get("ph", opts.get("sim_ph", 7.4)),
            }
        except Exception as e:
            send_log(f"[RESUME] [WARN] Impossible de lire checkpoint.json: {e}")

    manifest_path = os.path.join(resume_dir, "sdf_manifest.json")
    if os.path.exists(manifest_path):
        try:
            with open(manifest_path, "r", encoding="utf-8") as f:
                manifest = json.load(f)
            for item in manifest.get("items", []):
                mol_id = item.get("mol_id", "")
                if mol_id and not item.get("error"):
                    sdf_file = item.get("sdf", "")

                    if sdf_file and not os.path.isabs(sdf_file):
                        sdf_file = os.path.join(resume_dir, os.path.basename(sdf_file))
                    mol_data[mol_id] = {
                        "smiles": item.get("smiles", ""),
                        "sdf_path": sdf_file,
                        "nb_predicted_H": item.get("nb_predicted_H", 0),
                        "nb_predicted_C": item.get("nb_predicted_C", 0),
                    }
        except Exception as e:
            send_log(f"[RESUME] [WARN] Impossible de lire sdf_manifest.json: {e}")

    mult_path = os.path.join(resume_dir, "all_multiplicities.json")
    if os.path.exists(mult_path):
        try:
            with open(mult_path, "r", encoding="utf-8") as f:
                mult_data = json.load(f)
            for mol in mult_data.get("molecules", []):
                mol_id = mol.get("mol_id", "")
                if mol_id and mol_id in mol_data:
                    mol_data[mol_id]["payload"] = mol
        except Exception as e:
            send_log(f"[RESUME] [WARN] Impossible de lire all_multiplicities.json: {e}")

    return mol_data, params

def _params_match(old_params: Dict[str, Any], new_params: Dict[str, Any]) -> bool:
    """Check whether essential simulation parameters match."""

    keys = ["freq_mhz", "solvent", "temperature_k", "ph"]
    for k in keys:
        old_val = old_params.get(k)
        new_val = new_params.get(k)

        if isinstance(old_val, float) and isinstance(new_val, float):
            if abs(old_val - new_val) > 0.01:
                return False
        elif old_val != new_val:
            return False
    return True

def _copy_sdf_to_output(src_path: str, dest_dir: str) -> str:
    """Copy an SDF file and upgrade legacy 2D coordinates for the 3D viewer."""
    if not os.path.exists(src_path):
        return ""
    dest_path = os.path.join(dest_dir, os.path.basename(src_path))
    try:
        import shutil
        shutil.copy2(src_path, dest_path)
        _upgrade_sdf_coordinates_to_3d(dest_path)
        return dest_path
    except Exception as e:
        send_log(f"[RESUME] [WARN] Erreur copie SDF {src_path}: {e}")
        return ""

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

# Molecular pattern helpers
def _detect_geminal_ch3_on_double_bond(mol: Chem.Mol, heavy_to_H: Dict[int, List[int]]) -> Dict[int, Tuple[int, str]]:
    """Detect geminal methyl pairs attached to a substituted double bond."""

    result = {}
    base = Chem.RemoveHs(mol)

    pattern = Chem.MolFromSmarts('[CH3;X4]-[C;X3]([CH3;X4])=[C;X3]-[#6]')
    if pattern is None:
        return result

    matches = mol.GetSubstructMatches(pattern)

    for match in matches:

        ch3_1_idx = match[0]
        c_quat_idx = match[1]
        ch3_2_idx = match[2]
        c_vinyl_idx = match[3]

        if ch3_1_idx not in heavy_to_H or ch3_2_idx not in heavy_to_H:
            continue
        if len(heavy_to_H.get(ch3_1_idx, [])) != 3 or len(heavy_to_H.get(ch3_2_idx, [])) != 3:
            continue

        if ch3_1_idx < ch3_2_idx:
            result[ch3_1_idx] = (ch3_2_idx, 'cis')
            result[ch3_2_idx] = (ch3_1_idx, 'trans')
        else:
            result[ch3_2_idx] = (ch3_1_idx, 'cis')
            result[ch3_1_idx] = (ch3_2_idx, 'trans')

    return result

def _detect_geminal_dimethyl_13c(mol: Chem.Mol) -> List[Tuple[int, int, int]]:
    """Detect geminal dimethyl carbon pairs."""

    result = []

    for c_idx in range(mol.GetNumAtoms()):
        c_atom = mol.GetAtomWithIdx(c_idx)

        if c_atom.GetSymbol() != 'C':
            continue

        if c_atom.GetHybridization() not in [rdchem.HybridizationType.SP2, rdchem.HybridizationType.SP3]:
            continue

        ch3_neighbors = []
        for neighbor in c_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and neighbor.GetDegree() == 1:
                ch3_neighbors.append(neighbor.GetIdx())

        if len(ch3_neighbors) == 2:
            result.append((c_idx, ch3_neighbors[0], ch3_neighbors[1]))

    return result

def _sorted_h_groups(hx: Dict[int, List[int]]) -> List[Tuple[int, List[int]]]:
    """Return hydrogen groups sorted by heavy-atom index."""
    return [(k, sorted(v)) for k, v in sorted(hx.items(), key=lambda kv: kv[0]) if v]

# Graph dataset
class GraphDataset_sample:
    """Inference dataset compatible with the original graph pipeline."""

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

def _safe_mol_to_smiles(mol: Chem.Mol) -> str:
    """Convert a molecule to SMILES with fallbacks."""
    try:
        return Chem.MolToSmiles(mol)
    except Exception:
        try:

            return Chem.MolToSmiles(mol, kekuleSmiles=False)
        except Exception:
            try:

                return Chem.MolToSmiles(mol, canonical=True, kekuleSmiles=False)
            except Exception:

                return "UNKNOWN"

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
    mol_dict['smi'].append(_safe_mol_to_smiles(mol))

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

# Dataset preparation
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

def save_as_npy_for_1H(mol: Chem.Mol, output_dir: str) -> Tuple[str, Dict[int, List[int]], Chem.Mol]:
    """Prepare the 1H inference dataset and preserve atom-index consistency."""

    mol_with_h = Chem.AddHs(mol, addCoords=True)
    sdf_path = mol2SDF(mol_with_h, output_dir, name="tmp_1H")

    atom_list, connectivity_list, _, _, mol_read = get_molecule_data(sdf_path)
    x_h_conn = get_x_h_connectivity(connectivity_list, atom_list)

    mol_rd = Chem.RemoveHs(mol_read)

    sanitized = False
    try:
        Chem.SanitizeMol(mol_rd)
        sanitized = True
    except Exception:
        pass

    if not sanitized:
        try:

            Chem.SanitizeMol(mol_rd, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE)
            Chem.SetAromaticity(mol_rd, Chem.AromaticityModel.AROMATICITY_SIMPLE)
            sanitized = True
        except Exception:
            pass

    if not sanitized:
        try:

            Chem.SanitizeMol(mol_rd, sanitizeOps=(
                Chem.SanitizeFlags.SANITIZE_ALL ^
                Chem.SanitizeFlags.SANITIZE_KEKULIZE ^
                Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
            ))
            mol_rd.UpdatePropertyCache(strict=False)
        except Exception:
            mol_rd.UpdatePropertyCache(strict=False)

    try:
        si = Chem.FindPotentialStereo(mol_rd)
        for element in si:
            if str(element.type) == 'Atom_Tetrahedral' and str(element.specified) == 'Specified':
                mol_rd.GetAtomWithIdx(element.centeredOn).SetProp('Chirality', str(element.descriptor))
            elif str(element.type) == 'Bond_Double' and str(element.specified) == 'Specified':
                mol_rd.GetBondWithIdx(element.centeredOn).SetProp('Stereochemistry', str(element.descriptor))
    except Exception:

        pass

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

    return save_path, x_h_conn, mol_read

def save_as_npy_for_13C(mol: Chem.Mol, output_dir: str) -> str:
    mol_rd = Chem.Mol(mol)

    sanitized = False
    try:
        Chem.SanitizeMol(mol_rd)
        sanitized = True
    except Exception:
        pass

    if not sanitized:
        try:

            Chem.SanitizeMol(mol_rd, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE)
            Chem.SetAromaticity(mol_rd, Chem.AromaticityModel.AROMATICITY_SIMPLE)
            sanitized = True
        except Exception:
            pass

    if not sanitized:
        try:

            Chem.SanitizeMol(mol_rd, sanitizeOps=(
                Chem.SanitizeFlags.SANITIZE_ALL ^
                Chem.SanitizeFlags.SANITIZE_KEKULIZE ^
                Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
            ))
            mol_rd.UpdatePropertyCache(strict=False)
            sanitized = True
        except Exception:
            pass

    if not sanitized:

        try:
            mol_rd.UpdatePropertyCache(strict=False)
        except Exception:
            pass

    try:
        si = Chem.FindPotentialStereo(mol_rd)
        for element in si:
            if str(element.type) == 'Atom_Tetrahedral' and str(element.specified) == 'Specified':
                mol_rd.GetAtomWithIdx(element.centeredOn).SetProp('Chirality', str(element.descriptor))
            elif str(element.type) == 'Bond_Double' and str(element.specified) == 'Specified':
                mol_rd.GetBondWithIdx(element.centeredOn).SetProp('Stereochemistry', str(element.descriptor))
    except Exception:

        pass

    smiles = _safe_mol_to_smiles(mol_rd)
    if '.' in smiles:

        frags = Chem.GetMolFrags(mol_rd, asMols=True, sanitizeFrags=False)
        mol_rd = max(frags, key=lambda m: m.GetNumAtoms())

        try:
            Chem.SanitizeMol(mol_rd)
        except Exception as e:
            send_log(f"Kekulization échoué pour le fragment, utilisation d'une approche alternative: {e}")
            try:
                Chem.SanitizeMol(mol_rd, sanitizeOps=Chem.SANITIZE_ALL^Chem.SANITIZE_KEKULIZE)
                Chem.SetAromaticity(mol_rd, Chem.AromaticityModel.AROMATICITY_SIMPLE)
            except Exception as e2:
                send_log(f"Sanitization alternative échoué: {e2}")
        send_log(f"Multi-fragment molecule: using largest fragment: {_safe_mol_to_smiles(mol_rd)}")

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

# Model loading and inference
def load_model(target: str, model_path: str, sample_npz_path: str):
    """Load an NMR prediction model."""
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

            if isinstance(state.get('state_dict'), dict):
                state_dict = state['state_dict']
            else:
                state_dict = state
        elif hasattr(state, 'state_dict'):

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
    """Load training mean and standard deviation with fallbacks."""
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

# Shift reconstruction
def create_shift_list_1H(
    mol_with_H: Chem.Mol,
    x_h_connectivity: Dict[int, List[int]],
    preds_1h: np.ndarray
) -> List[float]:
    """Assign one predicted shift to each hydrogen group."""

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
    """Combine per-atom 1H and 13C shift vectors."""

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
    sdf_property_key: str = "(1)",
    mol_with_h: Optional[Chem.Mol] = None,
) -> str:
    # Reuse the prediction molecule so the simulation can use this same
    # conformer. Coordinates change; atom order, graph and shifts do not.
    molH = mol_with_h if mol_with_h is not None else Chem.AddHs(mol_original)
    _ensure_3d_conformer(molH)

    list_H = create_shift_list_1H(molH, x_h_connectivity, shifts_1h)
    list_C = create_shift_list_13C(molH, shifts_13c)
    all_shifts = combine_1H_13C_lists(molH, list_H, list_C)

    ts = int(time.time())
    if not outfile_name:
        outfile_name = f"nmr_prediction_{ts}.sdf"
    sdf_path = os.path.join(output_dir, outfile_name)

    molH.SetProp('SMILES', smiles)
    molH.SetProp('Method', 'ML_Prediction_1H_13C_SDF_Only')
    molH.SetProp('Timestamp', str(ts))
    shifts_str = " ".join(f"{x:.3f}" for x in all_shifts)
    molH.SetProp(sdf_property_key, shifts_str)
    _write_2d_mol_sidecar(molH, sdf_path)
    writer = rdMolFiles.SDWriter(sdf_path)
    writer.write(molH)
    writer.close()

    return sdf_path


def _has_3d_conformer(mol: Chem.Mol) -> bool:
    if mol is None or mol.GetNumConformers() == 0:
        return False
    conf = mol.GetConformer()
    if conf.Is3D():
        return True
    z_values = [conf.GetAtomPosition(i).z for i in range(mol.GetNumAtoms())]
    return bool(z_values) and (max(z_values) - min(z_values) > 1e-3)


def _ensure_3d_conformer(mol: Chem.Mol) -> bool:
    """Generate a reproducible 3D conformer; keep a 2D fallback on failure."""
    if _has_3d_conformer(mol):
        return True

    mol.RemoveAllConformers()
    embedded = False
    try:
        # Legacy SDFs are read with sanitize=False so their stored properties
        # remain untouched; ETKDG still needs a populated valence cache.
        mol.UpdatePropertyCache(strict=False)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.maxIterations = 200
        params.numThreads = 1
        embedded = AllChem.EmbedMolecule(mol, params) == 0
    except Exception as exc:
        send_log(f"[3D] ETKDG conformer generation failed: {exc}")

    if embedded:
        try:
            if AllChem.MMFFHasAllMoleculeParams(mol):
                AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            elif AllChem.UFFHasAllMoleculeParams(mol):
                AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        except Exception as exc:
            # ETKDG coordinates remain usable if force-field optimization fails.
            send_log(f"[3D] Force-field optimization skipped: {exc}")
        return _has_3d_conformer(mol)

    send_log("[3D] Conformer generation failed; retaining a 2D layout as fallback")
    try:
        from rdkit.Chem import rdCoordGen
        rdCoordGen.AddCoords(mol)
    except Exception:
        try:
            AllChem.Compute2DCoords(mol)
        except Exception:
            pass
    return False


def _make_2d_molblock(mol: Chem.Mol) -> str:
    """Create an independent 2D depiction with the exact same atom ordering."""
    depiction = Chem.Mol(mol)
    depiction.RemoveAllConformers()
    depiction.UpdatePropertyCache(strict=False)
    try:
        from rdkit.Chem import rdCoordGen
        rdCoordGen.AddCoords(depiction)
    except Exception:
        AllChem.Compute2DCoords(depiction)
    return rdMolFiles.MolToMolBlock(depiction)


def _write_2d_mol_sidecar(mol: Chem.Mol, sdf_path: str) -> str:
    """Write a same-index 2D MOL beside the 3D prediction SDF."""
    mol_path = os.path.splitext(sdf_path)[0] + ".mol"
    with open(mol_path, "w", encoding="ascii", newline="\n") as handle:
        handle.write(_make_2d_molblock(mol))
    return mol_path


def _upgrade_sdf_coordinates_to_3d(sdf_path: str) -> bool:
    """Upgrade a legacy prediction SDF and regenerate its 2D MOL sidecar."""
    suppl = rdMolFiles.SDMolSupplier(sdf_path, removeHs=False, sanitize=False)
    mols = [mol for mol in suppl if mol is not None]
    if len(mols) != 1:
        send_log(f"[RESUME] [WARN] SDF 3D upgrade skipped: expected one molecule in {sdf_path}")
        return False
    mol = mols[0]
    has_3d = _has_3d_conformer(mol)
    if not has_3d and not _ensure_3d_conformer(mol):
        return False
    needs_write = not has_3d
    if mol.HasProp("NMRScout_2D_MolBlock_Base64"):
        mol.ClearProp("NMRScout_2D_MolBlock_Base64")
        needs_write = True
    try:
        _write_2d_mol_sidecar(mol, sdf_path)
    except Exception as exc:
        send_log(f"[RESUME] [WARN] 2D MOL sidecar generation failed: {exc}")
        return False
    if not needs_write:
        return True
    writer = rdMolFiles.SDWriter(sdf_path)
    try:
        writer.write(mol)
    finally:
        writer.close()
    send_log(f"[RESUME] Legacy SDF upgraded with 3D coordinates: {os.path.basename(sdf_path)}")
    return True

# Simulation helpers
def _build_x_h_connectivity_from_mol_with_H(molH: Chem.Mol, solvent: str = "CDCl3") -> Dict[int, List[int]]:
    """Build heavy-atom to hydrogen connectivity."""

    protic_solvents = {'D2O', 'd2o', 'CD3OD', 'cd3od', 'MeOD', 'meod', 'CD3COCD3', 'cd3cocd3'}
    exclude_exchangeable = solvent in protic_solvents

    xh: Dict[int, List[int]] = {}
    for a in molH.GetAtoms():
        if a.GetSymbol() != "H":
            continue
        hidx = a.GetIdx()
        nbs = [n.GetIdx() for n in a.GetNeighbors()]
        if not nbs:
            continue
        heavy = nbs[0]
        heavy_atom = molH.GetAtomWithIdx(heavy)
        heavy_symbol = heavy_atom.GetSymbol()

        if exclude_exchangeable and heavy_symbol in ('O', 'N', 'S'):

            continue

        xh.setdefault(heavy, []).append(hidx)
    for k in xh:
        xh[k] = sorted(xh[k])
    return xh

def _proton_shifts_ppm_dict_from_list(molH: Chem.Mol, per_atom_H_ppm: List[float]) -> Dict[int, float]:
    """Build a proton-index to shift mapping."""
    if len(per_atom_H_ppm) < molH.GetNumAtoms():
        per_atom_H_ppm = list(per_atom_H_ppm) + [0.0] * (molH.GetNumAtoms() - len(per_atom_H_ppm))
    out: Dict[int, float] = {}
    for a in molH.GetAtoms():
        if a.GetSymbol() == "H":
            out[a.GetIdx()] = float(per_atom_H_ppm[a.GetIdx()])
    return out

def _write_multiplicity_json(out_dir: str, mol_id: str, payload: dict) -> str:
    path = os.path.join(out_dir, f"{sanitize_filename(mol_id)}_multiplicity.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
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
        send_log(f"[WARN] Impossible d'ajouter des propriétés SDF: {e}")

def _update_sdf_with_corrected_shifts(
    sdf_path: str,
    groups_items: List[Dict[str, Any]],
    property_key: str = "(1)"
) -> None:
    """Update 1H shifts stored in an SDF file."""

    try:
        suppl = rdMolFiles.SDMolSupplier(sdf_path, removeHs=False, sanitize=False)
        mols = [m for m in suppl if m is not None]
        if not mols:
            send_log("[WARN] Impossible de lire le SDF pour mise à jour des shifts")
            return

        m = mols[0]

        current_shifts_str = m.GetProp(property_key) if m.HasProp(property_key) else ""
        if not current_shifts_str:
            send_log("[WARN] Pas de shifts trouvés dans le SDF")
            return

        current_shifts = [float(x) for x in current_shifts_str.split()]

        corrected_h_shifts: Dict[int, float] = {}
        for group in groups_items:
            h_indices = group.get("H_indices", [])
            h_shifts = group.get("H_shifts_ppm", [])
            heavy_symbol = group.get("heavy_symbol", "C")

            if heavy_symbol not in ("N", "O", "S"):
                for i, h_idx in enumerate(h_indices):
                    if i < len(h_shifts):
                        corrected_h_shifts[h_idx] = h_shifts[i]

        updated = False
        for h_idx, corr_shift in corrected_h_shifts.items():
            if h_idx < len(current_shifts):
                old_val = current_shifts[h_idx]
                if abs(old_val - corr_shift) > 0.001:
                    current_shifts[h_idx] = corr_shift
                    updated = True

        if updated:

            new_shifts_str = " ".join(f"{x:.3f}" for x in current_shifts)
            m.SetProp(property_key, new_shifts_str)
            m.SetProp("Shifts_Corrected", "diastereotopic_delta_applied")

            w = rdMolFiles.SDWriter(sdf_path)
            w.write(m)
            w.close()

    except Exception as e:
        send_log(f"[WARN] Erreur mise à jour SDF shifts: {e}")

def save_spectral_json(
    spectral_data: List[Dict[str, Any]],
    output_dir: str,
    filename: str = "spectral_data_1H.json",
) -> str:
    """Save spectral data to JSON."""

    output_path = os.path.join(output_dir, filename)

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(spectral_data, f, indent=2, ensure_ascii=False)

    send_log(f"Spectral data saved: {output_path} ({len(spectral_data)} molecules)")
    return output_path

# Per-molecule orchestration
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
    return_payload: bool = False,
) -> Tuple[str, int, int, ModelBundle, Optional[str], Optional[dict]]:

    original_input_smiles = (smiles or "").strip()

    smiles = original_input_smiles
    if not smiles:
        raise ValueError("SMILES vide.")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES invalide : {smiles}")

    if multiplicity_mod and hasattr(multiplicity_mod, '_validate_molecular_structure'):
        try:
            multiplicity_mod._validate_molecular_structure(mol, original_input_smiles)
        except ValueError as ve:
            send_error(f"[{mol_id}] Structure moléculaire invalide", error=str(ve))
            raise

    work_dir = ensure_dir(os.path.join(out_dir, f"tmp_{sanitize_filename(mol_id)}"))

    with _MODEL_LOCK:
        if models.mean_std_1h is None:
            models.mean_std_1h = load_std_mean("1H")
        if models.mean_std_13c is None:
            models.mean_std_13c = load_std_mean("13C")

    train_y_mean_H, train_y_std_H = models.mean_std_1h
    train_y_mean_C, train_y_std_C = models.mean_std_13c

    npz_1h, x_h_conn, mol_with_h = save_as_npy_for_1H(mol, work_dir)

    with _MODEL_LOCK:
        if models.net_1h is None:
            models.net_1h = load_model("1H", model_1h_path, npz_1h)
    pred_1h = inference_model(models.net_1h, npz_1h, "1H", train_y_mean_H, train_y_std_H)

    expected_groups = len(_sorted_h_groups(x_h_conn))

    fp_predictor = get_fingerprint_predictor()

    fp_1h_used = 0
    ml_1h_used = 0
    fp_1h_exact = 0
    fp_1h_r6 = 0
    fp_1h_r5 = 0
    fp_1h_r4 = 0
    fp_1h_r3 = 0
    fp_1h_r2 = 0

    solvent_for_db = (sim_params or {}).get('solvent', 'CDCl3')

    sorted_groups = _sorted_h_groups(x_h_conn)
    new_pred_1h = list(pred_1h) if len(pred_1h) else []
    prediction_source_1h = ['ML'] * len(new_pred_1h)

    individual_h_shifts: Dict[int, float] = {}
    exact_group_shifts = (
        fp_predictor.predict_exact_proton_groups(mol_with_h, solvent_for_db)
        if fp_predictor is not None and fp_predictor.is_available()
        else {}
    )

    for gi, (heavy_idx, h_list) in enumerate(sorted_groups):
        if gi >= len(new_pred_1h):
            break

        if not h_list:
            ml_1h_used += 1
            continue

        group_fp_found = False
        group_fp_results = []
        if exact_group_shifts and all(h_idx in exact_group_shifts for h_idx in h_list):
            group_fp_results = [
                (h_idx, float(exact_group_shifts[h_idx]), 99, 0)
                for h_idx in h_list
            ]
        else:
            for h_idx in h_list:
                if fp_predictor is not None and fp_predictor.is_available():
                    result_fp = fp_predictor.predict_atom(mol_with_h, h_idx, solvent_for_db, '1H')
                    if result_fp is not None:
                        shift_fp, std_fp, count_fp, radius_used = result_fp
                        group_fp_results.append((h_idx, float(shift_fp), count_fp, radius_used))

        # The three protons of a methyl group are equivalent on the NMR time
        # scale.  Per-hydrogen fingerprint keys can become unstable when atom
        # ordering or unspecified stereochemistry differs from the database.
        # Reject an incoherent fingerprint assignment for the complete methyl
        # group instead of mixing unrelated shifts (for example 10.14, 3.55,
        # and 3.46 ppm in one CH3 group).
        heavy_atom = mol_with_h.GetAtomWithIdx(heavy_idx)
        is_methyl = heavy_atom.GetSymbol() == 'C' and len(h_list) == 3
        if is_methyl and group_fp_results:
            methyl_values = [item[1] for item in group_fp_results]
            if len(group_fp_results) != 3 or max(methyl_values) - min(methyl_values) > 0.15:
                send_log(
                    f"[{mol_id}] [WARN] Rejected incoherent fingerprint CH3 assignment "
                    f"at atom {heavy_idx}: {methyl_values}; using the ML group prediction"
                )
                group_fp_results = []

        if group_fp_results:
            group_shift = float(np.median([item[1] for item in group_fp_results]))
            new_pred_1h[gi] = group_shift
            group_fp_found = True

            for h_idx, shift_fp, count_fp, radius_used in group_fp_results:
                # All methyl protons receive the same robust group value.
                individual_h_shifts[h_idx] = group_shift if is_methyl else shift_fp
                fp_1h_used += 1
                is_exact = (radius_used == 0 and count_fp == 99)
                if is_exact:
                    fp_1h_exact += 1
                elif radius_used == 6:
                    fp_1h_r6 += 1
                elif radius_used == 5:
                    fp_1h_r5 += 1
                elif radius_used == 4:
                    fp_1h_r4 += 1
                elif radius_used == 3:
                    fp_1h_r3 += 1
                else:
                    fp_1h_r2 += 1

            if gi < len(prediction_source_1h):
                radii = {item[3] for item in group_fp_results}
                exact_group = all(item[3] == 0 and item[2] == 99 for item in group_fp_results)
                prediction_source_1h[gi] = 'FP-EXACT' if exact_group else f"FP-R{min(radii)}"

        if not group_fp_found:
            ml_1h_used += 1

    if fp_1h_used > 0:
        pred_1h = np.array(new_pred_1h)
        send_log(f"[{mol_id}] [OK] 1H sources: FP={fp_1h_used} (EXACT={fp_1h_exact}, R6={fp_1h_r6}, R5={fp_1h_r5}, R4={fp_1h_r4}, R3={fp_1h_r3}, R2={fp_1h_r2}), ML={ml_1h_used}")

    corrections_ch2_indole = 0
    pred_1h_mutable = list(pred_1h) if len(pred_1h) else []

    for gi, (heavy_idx, h_list) in enumerate(sorted_groups):
        if gi >= len(pred_1h_mutable) or not h_list:
            continue

        if len(h_list) == 2:
            h_idx = h_list[0]
            atom_h = mol_with_h.GetAtomWithIdx(h_idx)

            neighbors_of_h = list(atom_h.GetNeighbors())
            if not neighbors_of_h:
                continue
            carbon = neighbors_of_h[0]

            if carbon.GetSymbol() != 'C' or carbon.GetHybridization() != Chem.HybridizationType.SP3:
                continue

            c_neighbors = [n for n in carbon.GetNeighbors() if n.GetSymbol() == 'C']
            if len(c_neighbors) != 2:
                continue

            has_indole = False
            has_vinyl = False

            for neighbor_c in c_neighbors:

                if neighbor_c.GetHybridization() == Chem.HybridizationType.SP2:

                    if neighbor_c.GetIsAromatic():

                        for n2 in neighbor_c.GetNeighbors():
                            if n2.GetSymbol() == 'N' and n2.GetIsAromatic():
                                has_indole = True
                                break

                        if not has_indole:
                            for n2 in neighbor_c.GetNeighbors():
                                if n2.GetSymbol() == 'C' and n2.GetIsAromatic():
                                    for n3 in n2.GetNeighbors():
                                        if n3.GetSymbol() == 'N' and n3.GetIsAromatic():
                                            has_indole = True
                                            break
                                if has_indole:
                                    break

                    for bond in neighbor_c.GetBonds():
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            other_idx = bond.GetOtherAtomIdx(neighbor_c.GetIdx())
                            other_atom = mol.GetAtomWithIdx(other_idx)
                            if other_atom.GetSymbol() == 'C' and other_idx != carbon.GetIdx():
                                has_vinyl = True
                                break

            if has_indole and has_vinyl:
                old_shift = pred_1h_mutable[gi]

                if 3.0 <= old_shift <= 3.6:
                    correction = 0.30
                    pred_1h_mutable[gi] = old_shift + correction
                    corrections_ch2_indole += 1

    if corrections_ch2_indole > 0:
        pred_1h = np.array(pred_1h_mutable)

    pred_13c = np.array([])
    npz_13c = None
    try:
        npz_13c = save_as_npy_for_13C(mol, work_dir)
        with _MODEL_LOCK:
            if models.net_13c is None:
                models.net_13c = load_model("13C", model_13c_path, npz_13c)
        pred_13c = inference_model(models.net_13c, npz_13c, "13C", train_y_mean_C, train_y_std_C)
    except Exception as e:
        send_log(f"[{mol_id}] [WARN] Prédiction 13C échouée: {e}")

        try:
            mol_clean = Chem.Mol(mol)
            Chem.SanitizeMol(mol_clean, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE)
            Chem.SetAromaticity(mol_clean, Chem.AromaticityModel.AROMATICITY_SIMPLE)
            npz_13c = save_as_npy_for_13C(mol_clean, work_dir)
            with _MODEL_LOCK:
                if models.net_13c is None:
                    models.net_13c = load_model("13C", model_13c_path, npz_13c)
            pred_13c = inference_model(models.net_13c, npz_13c, "13C", train_y_mean_C, train_y_std_C)
        except Exception as e2:
            send_log(f"[{mol_id}] [WARN] Récupération 13C échouée: {e2}")
            pred_13c = np.zeros(sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'C'))

    hose_13c_used = 0
    fp_13c_used = 0
    ml_13c_used = 0
    hose_13c_s4 = 0
    hose_13c_s3 = 0
    fp_13c_exact = 0
    fp_13c_r6 = 0
    fp_13c_r5 = 0
    fp_13c_r4 = 0
    fp_13c_r3 = 0
    fp_13c_r2 = 0

    if len(pred_13c) > 0:

        carbon_indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'C']

        if len(carbon_indices) == len(pred_13c):
            new_pred_13c = list(pred_13c)
            mol_with_h_13c = Chem.AddHs(mol)
            prediction_source_13c = ['ML'] * len(new_pred_13c)

            for ci, c_idx in enumerate(carbon_indices):
                ml_value = new_pred_13c[ci]
                found = False

                if fp_predictor is not None and fp_predictor.is_available():
                    result_fp = fp_predictor.predict_atom(mol_with_h_13c, c_idx, solvent_for_db, '13C')
                    if result_fp is not None:
                        shift_fp, std_fp, count_fp, radius_used = result_fp
                        new_pred_13c[ci] = shift_fp
                        fp_13c_used += 1

                        is_exact = (radius_used == 0 and count_fp == 99)
                        if is_exact:
                            fp_13c_exact += 1
                        elif radius_used == 6:
                            fp_13c_r6 += 1
                        elif radius_used == 5:
                            fp_13c_r5 += 1
                        elif radius_used == 4:
                            fp_13c_r4 += 1
                        elif radius_used == 3:
                            fp_13c_r3 += 1
                        else:
                            fp_13c_r2 += 1
                        prediction_source_13c[ci] = 'FP-EXACT' if is_exact else f'FP-R{radius_used}'
                        found = True

                if not found:
                    ml_13c_used += 1

            if fp_13c_used > 0:
                pred_13c = np.array(new_pred_13c)
                send_log(f"[{mol_id}] [OK] 13C sources: FP={fp_13c_used} (EXACT={fp_13c_exact}, R6={fp_13c_r6}, R5={fp_13c_r5}, R4={fp_13c_r4}, R3={fp_13c_r3}, R2={fp_13c_r2}), ML={ml_13c_used}")
            else:
                ml_13c_used = len(carbon_indices)
        else:
            send_log(f"[{mol_id}] [WARN] Mismatch carbones ({len(carbon_indices)}) vs prédictions ({len(pred_13c)})")

    if len(pred_13c) > 0:
        try:
            carbon_indices_corr = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'C']

            geminal_dimethyl_pairs = _detect_geminal_dimethyl_13c(mol)

            corrections_applied_13c = 0
            for c_central_idx, ch3_1_idx, ch3_2_idx in geminal_dimethyl_pairs:
                try:

                    ci_1 = carbon_indices_corr.index(ch3_1_idx)
                    ci_2 = carbon_indices_corr.index(ch3_2_idx)

                    shift_1 = pred_13c[ci_1]
                    shift_2 = pred_13c[ci_2]

                    if abs(shift_1 - shift_2) < 1.0:

                        pred_13c[ci_1] = shift_1 + 3.5
                        pred_13c[ci_2] = shift_2 - 3.5
                        corrections_applied_13c += 1
                except (ValueError, IndexError):
                    pass
        except Exception as e:
            send_log(f"[{mol_id}] [WARN] Erreur correction 13C diastéréotopique: {e}")

    final_name = f"{sanitize_filename(mol_id)}_predicted.sdf"
    sdf_path = generate_sdf_with_shifts(
        smiles=original_input_smiles,
        mol_original=mol,
        shifts_1h=pred_1h,
        shifts_13c=pred_13c,
        x_h_connectivity=x_h_conn,
        output_dir=out_dir,
        outfile_name=final_name,
        sdf_property_key=property_key,
        mol_with_h=mol_with_h
    )

    sidecar_path = None
    payload = None
    if multiplicity_mod is not None:
        try:

            sp = sim_params or {}
            SimConfig = getattr(multiplicity_mod, "SimConfig", None)
            cfg = SimConfig() if SimConfig else None
            if cfg:
                cfg.separate_diastereotopic_protons = True
                cfg.spectrometer_mhz = float(sp.get("freq_mhz", 700.0))

                if "use_conformational_averaging" in sp:
                    cfg.use_conformational_averaging = bool(sp["use_conformational_averaging"])
                if "smart_conformer_mode" in sp:
                    cfg.smart_conformer_mode = bool(sp["smart_conformer_mode"])
                if "n_conformers" in sp:
                    cfg.n_conformers = int(sp["n_conformers"])
                if "fast_mode_for_large_molecules" in sp:
                    cfg.fast_mode_for_large_molecules = bool(sp["fast_mode_for_large_molecules"])

            molH_for_sim = mol_with_h

            use_3d = sp.get("use_3d", True)
            n_heavy_atoms_sim = molH_for_sim.GetNumHeavyAtoms()

            if use_3d and molH_for_sim.GetNumConformers() == 0 and n_heavy_atoms_sim <= 50:
                try:
                    AllChem.EmbedMolecule(molH_for_sim, randomSeed=42)
                    AllChem.MMFFOptimizeMolecule(molH_for_sim)
                except Exception:
                    pass

            solvent_for_sim = sp.get("solvent", "CDCl3")
            x_h_conn_sim = _build_x_h_connectivity_from_mol_with_H(molH_for_sim, solvent=solvent_for_sim)

            list_H = create_shift_list_1H(molH_for_sim, x_h_conn, pred_1h)
            proton_shifts_ppm = _proton_shifts_ppm_dict_from_list(molH_for_sim, list_H)

            geminal_ch3_info = _detect_geminal_ch3_on_double_bond(molH_for_sim, x_h_conn_sim)
            if geminal_ch3_info:
                half_delta = DELTA_GEMINAL_CH3_DOUBLE_BOND / 2.0

                processed_pairs = set()
                corrections_applied = 0
                for heavy_idx, (partner_idx, position) in geminal_ch3_info.items():

                    pair_key = tuple(sorted([heavy_idx, partner_idx]))
                    if pair_key in processed_pairs:
                        continue
                    processed_pairs.add(pair_key)

                    h_list1 = x_h_conn_sim.get(heavy_idx, [])
                    h_list2 = x_h_conn_sim.get(partner_idx, [])
                    if h_list1 and h_list2:
                        avg_shift1 = sum(proton_shifts_ppm.get(h, 0) for h in h_list1) / len(h_list1)
                        avg_shift2 = sum(proton_shifts_ppm.get(h, 0) for h in h_list2) / len(h_list2)

                        if abs(avg_shift1 - avg_shift2) < 0.01:

                            for h_idx in h_list1:
                                if h_idx in proton_shifts_ppm:
                                    if position == 'cis':
                                        proton_shifts_ppm[h_idx] -= half_delta
                                    else:
                                        proton_shifts_ppm[h_idx] += half_delta

                            partner_position = geminal_ch3_info[partner_idx][1]
                            for h_idx in h_list2:
                                if h_idx in proton_shifts_ppm:
                                    if partner_position == 'cis':
                                        proton_shifts_ppm[h_idx] -= half_delta
                                    else:
                                        proton_shifts_ppm[h_idx] += half_delta
                            corrections_applied += 1

            for h_idx, fp_shift in individual_h_shifts.items():
                if h_idx in proton_shifts_ppm:
                    proton_shifts_ppm[h_idx] = fp_shift

            per_proton = multiplicity_mod.simulate_multiplicity_per_proton(
                mol=molH_for_sim,
                x_h_connectivity=x_h_conn_sim,
                proton_shifts_ppm=proton_shifts_ppm,
                solvent=sp.get("solvent", "CDCl3"),
                spectrometer_mhz=float(sp.get("freq_mhz", 700.0)),
                temperature_k=float(sp.get("temperature_k", 298.15)),
                include_heteronuclear=True,
                use_3d=bool(sp.get("use_3d", True)),
                config=cfg,

                enable_dynamics=bool(sp.get("enable_dynamics", False)),
                ph=float(sp.get("ph", 7.4)),
            )

            repr_to_heavy: Dict[int, int] = {}
            for h_idx in per_proton.keys():
                nbs = [n.GetIdx() for n in molH_for_sim.GetAtomWithIdx(h_idx).GetNeighbors()]
                if nbs:
                    repr_to_heavy[h_idx] = nbs[0]

            heavy_to_groups: Dict[int, List[List[int]]] = {}

            diastereotopic_heavy_set = set()
            heavy_proton_count_in_per_proton: Dict[int, int] = {}
            for h_idx, heavy_idx in repr_to_heavy.items():
                heavy_proton_count_in_per_proton[heavy_idx] = heavy_proton_count_in_per_proton.get(heavy_idx, 0) + 1
            for heavy_idx, count in heavy_proton_count_in_per_proton.items():
                if count == 2:

                    diastereotopic_heavy_set.add(heavy_idx)

            for h_idx, heavy_idx in repr_to_heavy.items():
                data = per_proton[h_idx]

                J_key = tuple(sorted([(j.get("J_hz", 0), j.get("n_equiv", 1)) for j in data.get("J_couplings", [])]))
                label_key = (data.get("label", "s"), J_key)

                if heavy_idx not in heavy_to_groups:
                    heavy_to_groups[heavy_idx] = []

                if heavy_idx in diastereotopic_heavy_set:

                    heavy_to_groups[heavy_idx].append([h_idx])
                else:

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

            for heavy_idx, h_list in x_h_conn_sim.items():
                repr_h = None
                for h in h_list:
                    if h in per_proton:
                        repr_h = h
                        break
                if repr_h is not None:

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

            base_for_exch = Chem.RemoveHs(molH_for_sim)

            for heavy_idx in sorted(heavy_to_groups.keys()):
                for h_group in heavy_to_groups[heavy_idx]:
                    h_indices = sorted(h_group)

                    if not h_indices:
                        continue

                    h_idx_repr = h_indices[0]
                    data = per_proton.get(h_idx_repr, {})

                    label_str = data.get("label", "s")
                    label_idx = int(data.get("label_idx", 0))

                    if 0 <= heavy_idx < base_for_exch.GetNumAtoms():
                        heavy_atom = base_for_exch.GetAtomWithIdx(heavy_idx)
                        heavy_symbol = heavy_atom.GetSymbol()
                    else:
                        heavy_atom = None
                        heavy_symbol = "?"

                    note = "unique_or_equiv"

                    original_h_count = len(x_h_conn_sim.get(heavy_idx, []))
                    if len(h_indices) < original_h_count and original_h_count == 2:

                        if heavy_atom is not None:
                            hyb = heavy_atom.GetHybridization()
                            if hyb == rdchem.HybridizationType.SP2:
                                note = "vinyl_nonequiv"
                            else:
                                note = "diastereotopic"

                    partners = data.get("partners", [])
                    for p in partners:
                        partner_idx = p.get("partner_h_index")
                        J_hz = p.get("J_hz", 0)
                        if partner_idx and partner_idx not in h_indices:
                            partner_nbs = [n.GetIdx() for n in molH_for_sim.GetAtomWithIdx(partner_idx).GetNeighbors()]
                            if partner_nbs and partner_nbs[0] == heavy_idx:

                                if 0 < abs(J_hz) < 4:
                                    note = "vinyl_nonequiv"
                                elif 10 < abs(J_hz) < 15:
                                    note = "diastereotopic"
                                break

                    is_exchangeable = False
                    exch_type = None
                    try:

                        if hasattr(multiplicity_mod, "_is_exchangeable_H_on_X") and heavy_atom is not None:
                            exch_type = multiplicity_mod._is_exchangeable_H_on_X(base_for_exch, heavy_idx)
                            if exch_type:
                                is_exchangeable = True
                                note = "exchangeable"

                        if not is_exchangeable and heavy_atom is not None:
                            heavy_sym = heavy_atom.GetSymbol()

                            if heavy_sym in ('O', 'N', 'S'):
                                is_exchangeable = True
                                exch_type = f"{heavy_sym}H"
                                note = "exchangeable"
                    except Exception:
                        pass

                    for h in h_indices:
                        if 0 <= h < n_atoms:
                            mult_per_atom[h] = label_str

                    h_shifts_raw = []
                    for h in h_indices:
                        h_data = per_proton.get(h, {})
                        corrected_shift = h_data.get("corrected_shift_ppm")
                        h_note = h_data.get("note", "")

                        NONEQUIV_NOTES = ("diastereotopic", "vinyl_nonequiv", "allyl_nonequiv")
                        is_nonequiv = (note in NONEQUIV_NOTES or h_note in NONEQUIV_NOTES)

                        if is_nonequiv and corrected_shift is not None:

                            h_shifts_raw.append(float(corrected_shift))
                        elif h in individual_h_shifts:

                            h_shifts_raw.append(float(individual_h_shifts[h]))
                        elif corrected_shift is not None:

                            h_shifts_raw.append(float(corrected_shift))
                        else:

                            h_shifts_raw.append(float(proton_shifts_ppm.get(h, 0.0)))

                    h_shifts = []
                    if label_str == "s" and len(h_shifts_raw) > 1:

                        shifts_unique = set(round(s, 3) for s in h_shifts_raw)
                        if len(shifts_unique) == 1:

                            avg_shift = sum(h_shifts_raw) / len(h_shifts_raw)
                            h_shifts = [avg_shift]
                        else:

                            h_shifts = h_shifts_raw
                    else:

                        h_shifts = h_shifts_raw

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

                    for h in h_indices:
                        group_map_entries.append(f"{h}:{group_id}")

                    group_id += 1

            if write_mult_json or return_payload:
                payload = {
                    "mol_id": mol_id,
                    "smiles": original_input_smiles,
                    "solvent": sp.get("solvent", "CDCl3"),
                    "freq_MHz": float(sp.get("freq_mhz", 700.0)),
                    "temperature_K": float(sp.get("temperature_k", 298.15)),
                    "groups": groups_items,
                    "nb_groups": len(groups_items),
                    "nb_protons": (nb_protons := sum(len(g["H_indices"]) for g in groups_items)),
                    "note": "Protons grouped by heavy atom; diastereotopic protons separated only if truly different",
                    "typical_linewidth_hz": float(cfg.line_width_hz) if cfg else 1.2,
                    "tolerance_grouping_hz": float(getattr(multiplicity_mod, "TOLERANCE_J_HZ", 0.15)),
                    "digital_resolution_hz": float(cfg.digital_resolution_hz) if cfg and cfg.digital_resolution_hz else 0.05,

                    "include_connected_groups": True,
                    "max_hamiltonian_spins": 8,
                    "force_all_hamiltonian": nb_protons <= 8
                }

                if write_mult_json and not return_payload:
                    sidecar_path = _write_multiplicity_json(out_dir, mol_id, payload)

            if attach_labels_to_sdf:
                props = {
                    "H1_MULT_LABELS": " ".join(mult_per_atom),
                    "H1_MULT_GROUPS": str(group_id),
                    "H1_MULT_MAP": " ".join(group_map_entries),
                }
                _append_props_to_sdf(sdf_path, props)

            if groups_items:
                _update_sdf_with_corrected_shifts(sdf_path, groups_items, property_key)

        except Exception as e:
            send_log(f"[{mol_id}] [WARN] Simulation multiplicités échouée: {e}")
            send_log(traceback.format_exc())
            payload = None

    try:
        shutil.rmtree(work_dir, ignore_errors=True)
    except Exception:
        pass

    nH = int(len(pred_1h))
    nC = int(len(pred_13c))
    return sdf_path, nH, nC, models, sidecar_path, payload

# Single and CSV execution modes
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
            return_payload=False,
        )

        send_completed(True, sdf=sdf_path, mol_id=mol_id, nb_predicted_H=nH, nb_predicted_C=nC,
                       multiplicity_json=sidecar or "")
    except Exception as e:
        tb = traceback.format_exc()
        send_error(f"Prediction failed for '{mol_id}'", error=str(e), traceback=tb)
        send_completed(False, mol_id=mol_id, error=str(e))

def _run_csv_parallel(
    molecules: List[dict],
    out_dir: str,
    model_1h_path: str,
    model_13c_path: str,
    property_key: str,
    sim_params: Optional[dict],
    num_workers: int,
    manifest_rows: List[dict],
    all_multiplicity_payloads: List[dict],
) -> None:
    """Process CSV molecules with a process pool."""

    try:
        from parallel_prediction_worker import worker_process_molecule
    except ImportError:
        send_log("[WARN] Could not import parallel worker - falling back to sequential mode")

        raise

    tasks = []
    for mol in molecules:
        task = {
            'mol_data': {'mol_id': mol['mol_id'], 'smiles': mol['smiles']},
            'out_dir': out_dir,
            'model_1h_path': model_1h_path,
            'model_13c_path': model_13c_path,
            'property_key': property_key,
            'sim_params': sim_params,
        }
        tasks.append((mol, task))

    completed_count = 0
    with ProcessPoolExecutor(max_workers=num_workers) as executor:

        future_to_mol = {}
        for mol, task in tasks:
            future = executor.submit(
                worker_process_molecule,
                task['mol_data'],
                task['out_dir'],
                task['model_1h_path'],
                task['model_13c_path'],
                task['property_key'],
                task['sim_params'],
            )
            future_to_mol[future] = mol

        for future in as_completed(future_to_mol):
            mol = future_to_mol[future]
            try:
                success, result, mol_id = future.result()

                completed_count += 1
                _GLOBAL_MOLECULE_COUNTER["done"] = completed_count

                progress_pct = (completed_count / _GLOBAL_MOLECULE_COUNTER['total']) * 98
                status = "SDF" if success else "Error"
                send_progress(
                    progress_pct,
                    f"Molecule {completed_count}/{_GLOBAL_MOLECULE_COUNTER['total']} ({status})",
                    completed=completed_count,
                    total=_GLOBAL_MOLECULE_COUNTER['total']
                )

                if success:
                    send_log(f"[OK] Molecule {completed_count}/{_GLOBAL_MOLECULE_COUNTER['total']} processed: {result['mol_id']}")
                    if result.get('payload'):
                        all_multiplicity_payloads.append(result['payload'])
                else:
                    send_log(f"[{mol_id}] Error: {result.get('error', 'Unknown error')}")

                manifest_rows.append(result)

            except Exception as e:

                send_log(f"[{mol['mol_id']}] Worker error: {e}")
                completed_count += 1
                _GLOBAL_MOLECULE_COUNTER["done"] = completed_count

                manifest_rows.append({
                    'mol_id': mol['mol_id'],
                    'smiles': mol['smiles'],
                    'sdf': '',
                    'nb_predicted_H': 0,
                    'nb_predicted_C': 0,
                    'payload': None,
                    'error': f"Worker crashed: {str(e)}"
                })

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
    num_workers: int = 1,
    resume_from_dir: Optional[str] = None,
) -> None:
    ensure_dir(out_dir)

    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    existing_mol_data = {}
    existing_params = {}
    current_params = {}
    skipped_count = 0
    reused_payloads = []
    reused_manifest_rows = []

    if resume_from_dir and os.path.isdir(resume_from_dir):
        send_log(f"[RESUME] Chargement des prédictions existantes depuis: {resume_from_dir}")
        existing_mol_data, existing_params = _load_existing_predictions(resume_from_dir)
        send_log(f"[RESUME] {len(existing_mol_data)} molécules trouvées dans le dossier précédent")
        if existing_params:
            send_log(f"[RESUME] Paramètres précédents: freq={existing_params.get('freq_mhz')}MHz, "
                     f"solvent={existing_params.get('solvent')}, T={existing_params.get('temperature_k')}K, "
                     f"pH={existing_params.get('ph')}")

    if sim_params:
        current_params = {
            "freq_mhz": sim_params.get("freq_mhz", 700.0),
            "solvent": sim_params.get("solvent", "CDCl3"),
            "temperature_k": sim_params.get("temperature_k", 298.15),
            "ph": sim_params.get("ph", 7.4),
        }
    else:
        current_params = {"freq_mhz": 700.0, "solvent": "CDCl3", "temperature_k": 298.15, "ph": 7.4}

    params_match = _params_match(existing_params, current_params) if existing_params else False
    if resume_from_dir and existing_params:
        if params_match:
            send_log(f"[RESUME] [OK] Paramètres identiques - réutilisation des molécules déjà prédites")
        else:
            send_log(f"[RESUME] [WARN] Paramètres différents - toutes les molécules seront re-prédites")
            send_log(f"[RESUME] Nouveaux paramètres: freq={current_params.get('freq_mhz')}MHz, "
                     f"solvent={current_params.get('solvent')}, T={current_params.get('temperature_k')}K, "
                     f"pH={current_params.get('ph')}")

    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames or []
        if smiles_col not in headers or id_col not in headers:
            raise ValueError(f"CSV must contain columns '{id_col}' and '{smiles_col}' (found: {headers})")

    done = 0

    models = ModelBundle(net_1h=None, net_13c=None, mean_std_1h=None, mean_std_13c=None)
    manifest_rows: List[dict] = []
    all_multiplicity_payloads: List[dict] = []
    ok = 0
    ko = 0

    all_molecules = []
    seen_ids = set()
    duplicates = 0
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for idx, row in enumerate(reader, start=1):
            mol_id = (row.get(id_col) or f"mol_{idx}").strip()
            smi = (row.get(smiles_col) or "").strip()
            if mol_id in seen_ids:
                duplicates += 1
                continue
            seen_ids.add(mol_id)
            all_molecules.append({'mol_id': mol_id, 'smiles': smi, 'index': idx})
    if duplicates > 0:
        send_log(f"[DEDUP] {duplicates} doublons ignores ({len(all_molecules)} molecules uniques sur {len(all_molecules) + duplicates})")
    total = len(all_molecules)

    _GLOBAL_MOLECULE_COUNTER["total"] = total
    _GLOBAL_MOLECULE_COUNTER["done"] = 0

    molecules_to_process = []
    for mol_data in all_molecules:
        mol_id = mol_data['mol_id']

        if params_match and mol_id in existing_mol_data:
            existing = existing_mol_data[mol_id]
            existing_sdf = existing.get("sdf_path", "")

            if existing_sdf and os.path.exists(existing_sdf):

                new_sdf_path = _copy_sdf_to_output(existing_sdf, out_dir)
                if new_sdf_path:

                    skipped_count += 1
                    ok += 1
                    _GLOBAL_MOLECULE_COUNTER["done"] += 1

                    reused_manifest_rows.append({
                        "mol_id": mol_id,
                        "smiles": mol_data['smiles'],
                        "sdf": new_sdf_path,
                        "nb_predicted_H": existing.get("nb_predicted_H", 0),
                        "nb_predicted_C": existing.get("nb_predicted_C", 0),
                        "error": "",
                        "reused": True,
                    })

                    if "payload" in existing and existing["payload"]:
                        reused_payloads.append(existing["payload"])

                    send_log(f"[RESUME] [OK] {mol_id}: réutilisé (paramètres identiques)")
                    continue

        molecules_to_process.append(mol_data)

    if skipped_count > 0:
        send_log(f"[RESUME] {skipped_count} molécules réutilisées, {len(molecules_to_process)} à (re)prédire")

    manifest_rows.extend(reused_manifest_rows)
    all_multiplicity_payloads.extend(reused_payloads)

    for mol_data in molecules_to_process:
        mol_id = mol_data['mol_id']
        smi = mol_data['smiles']

        send_log(f"[{mol_id}] Starting prediction...")
        try:
            sdf_path, nH, nC, models, _, payload = process_one_smiles(
                smiles=smi, mol_id=mol_id, out_dir=out_dir, models=models,
                model_1h_path=model_1h_path, model_13c_path=model_13c_path,
                property_key=property_key,
                multiplicity_mod=multiplicity_mod,
                sim_params=sim_params,
                return_payload=True,
            )
            ok += 1

            if payload is not None:
                all_multiplicity_payloads.append(payload)

            _GLOBAL_MOLECULE_COUNTER["done"] += 1
            current_done = _GLOBAL_MOLECULE_COUNTER["done"]
            send_log(f"[OK] Molecule {current_done}/{total} processed: {os.path.basename(sdf_path)}")
            progress_pct = (current_done / total) * 98
            send_progress(
                progress_pct,
                f"Molecule {current_done}/{total} (SDF)",
                completed=current_done,
                total=total
            )

            manifest_rows.append({
                "mol_id": mol_id, "smiles": smi, "sdf": sdf_path,
                "nb_predicted_H": nH, "nb_predicted_C": nC, "error": ""
            })
        except Exception as e:
            ko += 1
            _GLOBAL_MOLECULE_COUNTER["done"] += 1
            current_done = _GLOBAL_MOLECULE_COUNTER["done"]

            send_log(f"[{mol_id}] Error: {e}")
            progress_pct = (current_done / total) * 98
            send_progress(
                progress_pct,
                f"Molecule {current_done}/{total} (Error)",
                completed=current_done,
                total=total
            )

            manifest_rows.append({
                "mol_id": mol_id, "smiles": smi, "sdf": "",
                "nb_predicted_H": 0, "nb_predicted_C": 0, "error": str(e)
            })

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

    manifest_path = os.path.join(out_dir, "sdf_manifest.json")
    with open(manifest_path, "w", encoding="utf-8") as mf:
        json.dump({
            "ok": ok,
            "ko": ko,
            "total": total,
            "output_dir": out_dir,
            "multiplicity_json": consolidated_multiplicity_path,
            "items": manifest_rows,
            "timestamp": time.time()
        }, mf, indent=2, ensure_ascii=False)

    send_log(f"CSV finished. Success: {ok} | Failures: {ko} | Manifest: {manifest_path}")
    send_completed(True, ok=ok, ko=ko, output_dir=out_dir, manifest=manifest_path,
                   multiplicity_json=consolidated_multiplicity_path)

# Dynamic simulator loading
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

# Command-line interface
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

    ap.add_argument("--resume-from-dir", default=None, dest="resume_from_dir",
                    help="Dossier d'une prédiction précédente pour reprise incrémentale. "
                         "Les molécules déjà prédites avec les mêmes paramètres seront réutilisées.")

    ap.add_argument("--solvent", default="CDCl3", help="Solvant (CDCl3, DMSO-d6, D2O, etc.)")
    ap.add_argument("--spectrometer_mhz", type=float, default=700.0, help="Fréquence du spectromètre (MHz)")
    ap.add_argument("--temperature_k", type=float, default=298.15, help="Température (K)")
    ap.add_argument("--ph", type=float, default=7.4, help="pH de la solution")

    ap.add_argument("--workers", type=int, default=0,
                    help="Nombre de threads parallèles pour CSV (0=auto, 1=séquentiel)")

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

    ap.add_argument("--sim_enable_dynamics", action="store_true",
                    help="Active le post-processeur dynamique (pH, échange, ring flipping).")
    ap.add_argument("--sim_ph", type=float, default=7.4,
                    help="pH de la solution pour le post-processeur dynamique (défaut: 7.4).")

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

    # Always retain the conditions selected in Prediction Parameters.  The
    # fingerprint lookup needs the selected solvent even if the optional
    # multiplicity simulator cannot be loaded.
    ph_value = float(getattr(args, 'sim_ph', getattr(args, 'ph', 7.4)))
    if ph_value < 0.0 or ph_value > 14.0:
        send_log(f"[WARNING] pH invalide: {ph_value}, fallback vers 7.4")
        ph_value = 7.4

    sim_params = {
        "solvent": args.sim_solvent,
        "freq_mhz": args.sim_freq_mhz,
        "temperature_k": args.sim_temperature_k,
        "use_3d": bool(args.sim_use_3d),
        "debug": bool(args.sim_debug),
        "enable_dynamics": bool(getattr(args, 'sim_enable_dynamics', False)),
        "ph": ph_value,
    }
    send_log(
        f"[PARAMETERS] solvent={sim_params['solvent']}, "
        f"frequency={sim_params['freq_mhz']} MHz, "
        f"temperature={sim_params['temperature_k']} K, pH={sim_params['ph']}"
    )

    multiplicity_mod = None
    if getattr(args, "simulate_multiplicity", False):
        if args.simulator_path:
            try:
                multiplicity_mod = _load_multiplicity_sim(args.simulator_path)
                send_log(f"[OK] Simulateur de multiplicités chargé depuis {args.simulator_path}")
            except Exception as e:
                send_log(f"[WARN] Échec chargement simulateur: {e}")
        else:
            send_log("[WARN] --simulate_multiplicity activé mais --simulator_path non fourni")

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
                sim_params=sim_params,
                num_workers=getattr(args, 'workers', 0),
                resume_from_dir=getattr(args, 'resume_from_dir', None)
            )
    except Exception as e:
        send_error("Erreur fatale", error=str(e), traceback=traceback.format_exc())
        send_completed(False)
        sys.exit(1)

if __name__ == "__main__":
    main()
