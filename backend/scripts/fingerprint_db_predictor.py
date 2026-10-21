#!/usr/bin/env python
# -*- coding: utf-8 -*-

import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


import os
import re
import pickle
import json
from collections import defaultdict
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdMolDescriptors import CalcNumHeteroatoms

REVISION = 1
KNOWN_IDS = {'bmse000064', 'bmse000937'}
INCHIKEY = 'KPGXRSRHYNQIFU-N-UHFFFAOYSA-N'

def reviewed_assignments():
    return [
        {'atom_index': i, 'nucleus': nucleus, 'ambiguity_code': 4,
         'candidates': candidates, 'legacy_shift': candidates[-1]}
        for indices, nucleus, candidates in [
            ([0, 1], '13C', [38.573, 33.325]),
            ([10, 11, 12, 13], '1H', [2.996, 2.428]),
        ] for i in indices
    ]

def ambiguous_atoms(mol):
    """Return unresolved atom assignments."""
    result = set()
    if mol.HasProp('BMRB_ID') and mol.GetProp('BMRB_ID').strip() in KNOWN_IDS:
        result.update((a['atom_index'], a['nucleus']) for a in reviewed_assignments())
    if mol.HasProp('NMR_Shift_Ambiguities'):
        for assignment in json.loads(mol.GetProp('NMR_Shift_Ambiguities')):
            result.add((int(assignment['atom_index']), assignment['nucleus']))
    values = {}
    if mol.HasProp('NMR_Shifts'):
        for line in mol.GetProp('NMR_Shifts').splitlines():
            parts = line.split()
            if len(parts) != 4:
                continue
            i, symbol, h, c = parts
            nucleus = {'H': '1H', 'C': '13C'}.get(symbol)
            if nucleus:
                value = float(h if symbol == 'H' else c)
                if value != 0:
                    values.setdefault((int(i), nucleus), set()).add(value)
    result.update(key for key, candidates in values.items() if len(candidates) > 1)
    if mol.HasProp('NMR_Shift_Resolutions'):
        resolutions = json.loads(mol.GetProp('NMR_Shift_Resolutions'))
        for assignment in resolutions:
            index = int(assignment['atom_index'])
            key = (index, assignment.get('nucleus'))
            if (assignment.get('method') == 'GISSMO_graph_verified'
                    and key[1] == '1H'
                    and assignment.get('source_url', '').startswith('https://gissmo.bmrb.io/entry/')
                    and assignment.get('source_atom_numbers')
                    and values.get(key) == {float(assignment['shift_ppm'])}
                    and mol.GetAtomWithIdx(index).GetSymbol() == 'H'):
                result.discard(key)
            if (assignment.get('method') == 'BMRB_HSQC_linked_to_GISSMO'
                    and key[1] == '13C'
                    and assignment.get('source_url', '').startswith('https://bmrb.io/ftp/')
                    and assignment.get('additional_source_hashes')
                    and values.get(key) == {float(assignment['shift_ppm'])}
                    and mol.GetAtomWithIdx(index).GetSymbol() == 'C'):
                linked = assignment.get('linked_proton_shifts', {})
                verified_h = {int(r['atom_index']): float(r['shift_ppm']) for r in resolutions
                              if r.get('method') == 'GISSMO_graph_verified' and r.get('nucleus') == '1H'}
                if linked and all(values.get((int(h), '1H')) == {float(v)}
                                  and verified_h.get(int(h)) == float(v)
                                  and mol.GetBondBetweenAtoms(index, int(h)) is not None
                                  for h, v in linked.items()):
                    result.discard(key)
    return result

try:
    from beyond_features_stereo import (
        compute_beyond_with_stereo_safe as compute_beyond_discrete_safe,
        ensure_gasteiger_charges
    )
    BEYOND_AVAILABLE = True
except ImportError:
    BEYOND_AVAILABLE = False
    def compute_beyond_discrete_safe(mol, atom_idx, include_chiral=True):

        return (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    def ensure_gasteiger_charges(mol):
        return False

def is_exchangeable_proton(mol, atom_idx):
    """Check whether a proton is exchangeable."""
    atom = mol.GetAtomWithIdx(atom_idx)

    if atom.GetSymbol() != 'H':
        return False

    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() in ['O', 'N', 'S']:
            return True

    return False

def classify_carbon_radius_requirement(mol, atom_idx):
    """Return the minimum fingerprint radius for an atom."""
    atom = mol.GetAtomWithIdx(atom_idx)
    symbol = atom.GetSymbol()

    if symbol == 'H':
        neighbors = list(atom.GetNeighbors())
        if not neighbors:
            return 2
        atom = neighbors[0]
        atom_idx = atom.GetIdx()
        symbol = atom.GetSymbol()

    if symbol != 'C':
        return 2

    ring_info = mol.GetRingInfo()

    def has_double_bond_to(a_idx, target_symbol=None):
        a = mol.GetAtomWithIdx(a_idx)
        for n in a.GetNeighbors():
            if target_symbol and n.GetSymbol() != target_symbol:
                continue
            bond = mol.GetBondBetweenAtoms(a_idx, n.GetIdx())
            if bond and bond.GetBondTypeAsDouble() >= 2.0:
                return True
        return False

    def count_rings(a_idx):
        return sum(1 for r in ring_info.AtomRings() if a_idx in r)

    def is_conjugated_sp2(a_idx):
        """Run is conjugated sp2."""
        a = mol.GetAtomWithIdx(a_idx)
        if 'SP2' not in str(a.GetHybridization()):
            return False
        for n in a.GetNeighbors():
            if n.GetSymbol() == 'C' and 'SP2' in str(n.GetHybridization()):
                return True
            if n.GetIsAromatic():
                return True
        return False

    def neighbor_symbols():
        return {n.GetSymbol() for n in atom.GetNeighbors()}

    def neighbor_has_hetero_aromatic():
        """Run neighbor has hetero aromatic."""
        for n in atom.GetNeighbors():
            if n.GetIsAromatic():
                for r in ring_info.AtomRings():
                    if n.GetIdx() in r:
                        ring_atoms = [mol.GetAtomWithIdx(i) for i in r]
                        if any(a.GetSymbol() in ('N', 'O', 'S') for a in ring_atoms):
                            return True
        return False

    n_H = sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == 'H')

    degree = sum(1 for n in atom.GetNeighbors() if n.GetSymbol() != 'H')
    is_arom = atom.GetIsAromatic()
    is_in_ring = atom.IsInRing()
    nsymbols = neighbor_symbols()

    if n_H == 3 and degree == 1:
        heavy_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() != 'H']
        if heavy_neighbors:
            parent = heavy_neighbors[0]
            if parent.GetSymbol() == 'C' and parent.GetTotalNumHs() == 0:

                ch3_count = sum(
                    1 for n in parent.GetNeighbors()
                    if n.GetSymbol() == 'C' and n.GetTotalNumHs() == 3 and n.GetDegree() == 1
                )
                if ch3_count >= 2:

                    for n in parent.GetNeighbors():
                        if n.GetSymbol() == 'C' and has_double_bond_to(parent.GetIdx(), 'C'):
                            return 99

    if n_H == 3 and degree == 1:
        heavy_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() != 'H']
        if heavy_neighbors:
            parent = heavy_neighbors[0]
            if parent.GetSymbol() == 'C':
                if not has_double_bond_to(parent.GetIdx()):
                    if not parent.GetIsAromatic():
                        return 2

    if n_H == 2 and not is_arom and not is_in_ring:
        has_hetero = any(s not in ('C', 'H') for s in nsymbols)
        has_db = any(
            mol.GetBondBetweenAtoms(atom_idx, n.GetIdx()).GetBondTypeAsDouble() >= 2.0
            for n in atom.GetNeighbors()
            if mol.GetBondBetweenAtoms(atom_idx, n.GetIdx()) is not None
        )

        has_adjacent_db = any(
            has_double_bond_to(n.GetIdx(), 'C')
            for n in atom.GetNeighbors()
            if n.GetSymbol() == 'C'
        )
        if not has_hetero and not has_db and not has_adjacent_db:
            return 2

    if n_H == 3 and degree == 1:
        heavy_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() != 'H']
        if heavy_neighbors:
            parent = heavy_neighbors[0]
            if parent.GetSymbol() in ('O', 'N', 'S'):
                return 2

    if n_H >= 1 and not is_arom:
        if has_double_bond_to(atom_idx, 'C'):
            if not is_conjugated_sp2(atom_idx):
                return 2

    if has_double_bond_to(atom_idx, 'O'):
        return 3

    if is_arom and count_rings(atom_idx) >= 2:
        return 3

    if is_arom:
        for n in atom.GetNeighbors():
            if n.GetSymbol() in ('O', 'N') and not n.GetIsAromatic():
                return 3

    if not is_arom and is_conjugated_sp2(atom_idx):
        return 3

    if n_H >= 1 and not is_arom:
        if neighbor_has_hetero_aromatic():
            return 3

    if n_H >= 2 and not is_arom and not is_in_ring:
        for n in atom.GetNeighbors():
            if n.GetSymbol() == 'C' and has_double_bond_to(n.GetIdx(), 'C'):
                return 3

    if is_in_ring and not is_arom:
        for r in ring_info.AtomRings():
            if atom_idx in r:
                ring_atoms = [mol.GetAtomWithIdx(i) for i in r]
                if any(a.GetSymbol() in ('N', 'O', 'S') for a in ring_atoms):
                    return 3

    if is_arom:
        return 3

    if n_H == 0 and count_rings(atom_idx) >= 3:
        return 4

    if n_H == 1 and has_double_bond_to(atom_idx, 'O'):

        for n in atom.GetNeighbors():
            if n.GetSymbol() == 'C':
                if n.GetIsAromatic():
                    return 4

                if is_conjugated_sp2(n.GetIdx()):
                    return 4

    return 2

def parse_mol_block_robust(mol_block: str):
    """Parse a MOL block with tolerant validation."""
    lines = mol_block.strip().split('\n')

    header_idx = None
    for i, line in enumerate(lines):
        if 'V2000' in line or 'V3000' in line:
            header_idx = i
            break

    if header_idx is None:
        return None, None

    mol = Chem.RWMol()

    atom_lines = []
    idx = header_idx + 1
    while idx < len(lines):
        line = lines[idx]

        if len(line.strip()) == 0 or line.strip().startswith('M '):
            break

        try:

            if len(line) < 34:

                break

            x_str = line[0:10].strip()
            y_str = line[10:20].strip()
            z_str = line[20:30].strip()
            symbol = line[30:34].strip()

            float(x_str)
            float(y_str)
            float(z_str)

            if len(symbol) > 0 and len(symbol) <= 3 and symbol[0].isalpha():
                atom_lines.append((symbol, line))
                idx += 1
            else:
                break
        except:

            break

    for symbol, line in atom_lines:
        try:
            atom = Chem.Atom(symbol)
            mol.AddAtom(atom)
        except:
            return None, None

    n_atoms = len(atom_lines)
    if n_atoms == 0:
        return None, None

    bond_start = header_idx + 1 + n_atoms
    for i in range(bond_start, len(lines)):
        line = lines[i].strip()
        if line.startswith('M '):
            break

        try:
            parts = line.split()
            if len(parts) < 3:
                continue

            atom1 = int(parts[0]) - 1
            atom2 = int(parts[1]) - 1
            bond_type = int(parts[2])

            if atom1 < 0 or atom1 >= n_atoms or atom2 < 0 or atom2 >= n_atoms:
                continue

            if bond_type == 1:
                bt = Chem.BondType.SINGLE
            elif bond_type == 2:
                bt = Chem.BondType.DOUBLE
            elif bond_type == 3:
                bt = Chem.BondType.TRIPLE
            else:
                bt = Chem.BondType.SINGLE

            mol.AddBond(atom1, atom2, bt)
        except:
            pass

    try:
        mol_clean = mol.GetMol()

        Chem.SanitizeMol(mol_clean)
        return mol_clean, n_atoms
    except:

        try:
            mol_clean = mol.GetMol()

            mol_clean.UpdatePropertyCache(strict=False)
            return mol_clean, n_atoms
        except:
            return None, None

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_DATA_DIR = os.path.join(_THIS_DIR, '..', 'ml_dependencies', 'data')

BMRB_DIR = os.path.join(_DATA_DIR, 'sdf_bmrb_rdkit_v2')
NMRSHIFTDB_FILE = os.path.join(_DATA_DIR, 'nmrshiftdb2withsignals.sd')
NMR_SPECTRA_FILE = os.path.join(_DATA_DIR, 'all_nmr_spectra_nmrshiftdb_format.sdf')
NPMRD_SDF_DIR = os.path.join(_DATA_DIR, 'NPMRD')
EXPERIMENTAL_DATA_DIR = os.path.join(_DATA_DIR, 'with_experimental_data')
DB_CACHE_FILE = os.path.join(_DATA_DIR, 'fingerprint_database.pkl')

RADIUS_ULTRA_EXACT = 6
RADIUS_HYPER_EXACT = 5
RADIUS_VERY_EXACT = 4
RADIUS_EXACT = 3
RADIUS_FALLBACK = 2
SUPPORTED_RADII = [RADIUS_ULTRA_EXACT, RADIUS_HYPER_EXACT, RADIUS_VERY_EXACT, RADIUS_EXACT, RADIUS_FALLBACK]

from rdkit.Chem import SubstructMatchParameters as _SMP
_SUBSTRUCT_PARAMS = _SMP()
_SUBSTRUCT_PARAMS.maxMatches = 1

def _substruct_match_with_timeout(mol, mol_canon, timeout=1.0):
    """Run substructure matching with a timeout."""
    import threading
    result = []
    def _run():
        try:
            result.append(mol.GetSubstructMatch(mol_canon, _SUBSTRUCT_PARAMS))
        except Exception:
            pass
    t = threading.Thread(target=_run, daemon=True)
    t.start()
    t.join(timeout)
    if t.is_alive():

        return ()
    return result[0] if result else ()

MAX_STD_13C = 3.0
MAX_STD_1H = 0.35

FALLBACK_MAX_STD_13C = 5.0
FALLBACK_MAX_STD_1H = 0.6

MIN_SAMPLES_R6 = 1
MIN_SAMPLES_R5 = 1
MIN_SAMPLES_R4 = 1
MIN_SAMPLES_R3 = 3
MIN_SAMPLES_R2 = 3

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)

SOLVENT_MAP = {

    'cdcl3': 'CDCl3',
    'chloroform': 'CDCl3',
    'chloroform-d': 'CDCl3',
    'chloroform-d1': 'CDCl3',
    'deuterochloroform': 'CDCl3',
    'chloroform-d1 (cdcl3)': 'CDCl3',
    'chcl3': 'CDCl3',
    'cdcl3,tms': 'CDCl3',

    'd2o': 'D2O',
    'deuterium oxide': 'D2O',
    'water': 'D2O',
    'h2o': 'D2O',
    'deuteriumoxide (d2o)': 'D2O',
    'water (d2o)': 'D2O',
    'water (h2o)': 'D2O',

    'dmso': 'DMSO-d6',
    'dmso-d6': 'DMSO-d6',
    'deuterated dmso': 'DMSO-d6',
    'dimethylsulphoxide-d6 (dmso-d6, c2d6os)': 'DMSO-d6',
    'dimethyl sulfoxide-d6': 'DMSO-d6',
    'dimco-d6': 'DMSO-d6',
    '(cd3)2so': 'DMSO-d6',
    'dmco-d6': 'DMSO-d6',
    'dimethylsulphoxide-d6 (dmso-d6, c2d6so)': 'DMSO-d6',

    'meod': 'MeOD',
    'meod-d4': 'MeOD',
    'methanol-d4': 'MeOD',
    'cd3od': 'MeOD',
    'meoh-d4': 'MeOD',
    'methanol-d4 (cd3od)': 'MeOD',
    'methanol-d4 (meod-d4)': 'MeOD',
    'methanol': 'MeOD',
    'meoh': 'MeOD',
    'methanol-d1 (ch3od)': 'MeOD',
    'methanol-d3(cd3oh)': 'MeOD',
    'methanol-d4/d2o (cd3od/d4)': 'MeOD',

    'acetone-d6': 'Acetone-d6',
    'acetone': 'Acetone-d6',
    'acetone-d6 ((cd3)2co)': 'Acetone-d6',
    'cd3cocd': 'Acetone-d6',
    'acetone-d6 + 1 drop d2o': 'Acetone-d6',

    'pyridine-d5': 'Pyridine-d5',
    'pyridin-d5 (c5d5n)': 'Pyridine-d5',

    'benzene-d6': 'Benzene-d6',
    'c6d6': 'Benzene-d6',
    'benzene-d6 (c6d6)': 'Benzene-d6',
    'benzene': 'Benzene-d6',

    'cd2cl2': 'CD2Cl2',
    'methylenchloride-d2 (cd2cl2)': 'CD2Cl2',
    'thf-d8': 'THF-d8',
    'tetrahydrofuran-d8 (thf-d8, c4d4o)': 'THF-d8',
    'toluene-d8': 'Toluene-d8',
    'toluol-d8': 'Toluene-d8',
    'toluene': 'Toluene-d8',
    'cd3cn': 'CD3CN',
    'acetonitrile-d3(cd3cn)': 'CD3CN',
    'cdcn': 'CD3CN',
    'cd3cn or cdcl3': 'CD3CN',
    'ccl4': 'CCl4',
    'tetrachloro-methane (ccl4)': 'CCl4',
    'ethanol': 'Ethanol',
}

def normalize_solvent(raw: str) -> str:
    """Normalize a solvent name."""
    if not raw:
        return 'Unreported'

    solvent = raw.strip()
    if ':' in solvent and solvent.split(':')[0].isdigit():
        solvent = solvent.split(':', 1)[1].strip()

    if ';' in solvent:
        solvent = solvent.split(';')[0].strip()

    key = solvent.lower().strip()

    unknown_variants = ['unknown', 'unreported', 'not reported', 'not specified',
                        'n/a', 'na', 'none', 'undefined', '']
    if key in unknown_variants:
        return 'Unreported'

    return SOLVENT_MAP.get(key, solvent)

def is_valid_solvent(solvent: str) -> bool:
    """Check whether a solvent label is valid."""
    if not solvent:
        return False

    return True

SOLVENT_VARIANTS = {}
for raw_key, normalized in SOLVENT_MAP.items():
    if normalized not in SOLVENT_VARIANTS:
        SOLVENT_VARIANTS[normalized] = set()

    if raw_key != normalized.lower():

        SOLVENT_VARIANTS[normalized].add(raw_key.title().replace('-D', '-d').replace(' D', ' d'))

SOLVENT_VARIANTS.setdefault('MeOD', set()).add('Methanol-D4 (MeOD-d4)')

def get_solvent_variants(normalized_solvent: str) -> list:
    """Return known variants of a normalized solvent."""
    variants = [normalized_solvent]
    if normalized_solvent in SOLVENT_VARIANTS:
        variants.extend(SOLVENT_VARIANTS[normalized_solvent])
    return variants

def get_atom_fingerprint_raw(mol, atom_idx: int, radius: int = RADIUS_VERY_EXACT, heavy_only: bool = False) -> int:
    """Compute a raw atom-centered Morgan fingerprint."""

    if heavy_only:
        mol_for_fp = Chem.RemoveHs(mol)

        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() == 'H':

            neighbors = list(atom.GetNeighbors())
            if neighbors:
                parent_idx = neighbors[0].GetIdx()

                heavy_count = 0
                for i in range(parent_idx):
                    if mol.GetAtomWithIdx(i).GetSymbol() != 'H':
                        heavy_count += 1
                atom_idx_in_no_h = heavy_count
            else:
                return 0
        else:

            heavy_count = 0
            for i in range(atom_idx):
                if mol.GetAtomWithIdx(i).GetSymbol() != 'H':
                    heavy_count += 1
            atom_idx_in_no_h = heavy_count

        if atom_idx_in_no_h >= mol_for_fp.GetNumAtoms():
            return 0
    else:
        mol_for_fp = mol
        atom_idx_in_no_h = atom_idx

    info = {}
    AllChem.GetMorganFingerprint(mol_for_fp, radius, bitInfo=info)

    atom_bits_with_radius = []
    for bit_id, atom_list in info.items():
        for a_idx, rad in atom_list:
            if a_idx == atom_idx_in_no_h:

                atom_bits_with_radius.append((bit_id, rad))

    if not atom_bits_with_radius:
        return 0

    return hash(tuple(sorted(atom_bits_with_radius)))

def get_atom_fingerprint(mol, atom_idx: int, radius: int = RADIUS_VERY_EXACT, heavy_only: bool = False) -> int:
    """Compute an atom-centered Morgan fingerprint."""
    atom = mol.GetAtomWithIdx(atom_idx)

    if atom.GetSymbol() == 'H':
        neighbors = list(atom.GetNeighbors())
        if neighbors:
            parent_idx = neighbors[0].GetIdx()

            parent_fp = get_atom_fingerprint_raw(mol, parent_idx, radius, heavy_only=heavy_only)

            h_fp = get_atom_fingerprint_raw(mol, atom_idx, radius, heavy_only=heavy_only)
            return (radius << 62) | (parent_fp << 32) | h_fp

    raw_fp = get_atom_fingerprint_raw(mol, atom_idx, radius, heavy_only=heavy_only)
    return (radius << 62) | raw_fp


def get_atom_env_string(mol, atom_idx: int) -> str:
    """Return a compact atom-environment description."""
    atom = mol.GetAtomWithIdx(atom_idx)
    env_parts = [atom.GetSymbol()]
    hyb = str(atom.GetHybridization()).replace('SP', 'sp')
    env_parts.append(hyb)
    env_parts.append(f"H{atom.GetTotalNumHs()}")
    neighbors = []
    for n in atom.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(atom_idx, n.GetIdx())
        bond_type = str(bond.GetBondType())[0]
        neighbors.append(f"{n.GetSymbol()}{bond_type}")
    neighbors.sort()
    env_parts.append("_".join(neighbors))

    fast_2d = get_fast_2d_descriptors(mol, atom_idx)
    env_parts.append(fast_2d)

    return "|".join(env_parts)

def get_fast_2d_descriptors(mol, atom_idx: int) -> str:
    """Return compact 2D atom descriptors."""
    atom = mol.GetAtomWithIdx(atom_idx)
    descriptors = []

    ring_info = "acyclic"
    if atom.IsInRing():
        ring_info_obj = mol.GetRingInfo()

        min_ring_size = 999
        for ring_size in [3, 4, 5, 6, 7, 8, 9, 10]:
            if ring_info_obj.IsAtomInRingOfSize(atom_idx, ring_size):
                min_ring_size = min(min_ring_size, ring_size)
        if min_ring_size < 999:
            ring_info = f"ring{min_ring_size}"
        else:
            ring_info = "ring"
    descriptors.append(ring_info)

    aniso_groups = 0
    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() in ['O', 'N'] and neighbor.GetTotalNumHs() == 0:
            aniso_groups += 1
    descriptors.append(f"aniso{aniso_groups}")

    heavy_neighbors = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() > 1)
    descriptors.append(f"heavy{heavy_neighbors}")

    return "|".join(descriptors)

def _process_mol_block_worker(args):
    """Process one molecular block."""
    from rdkit.Chem import inchi as rdkit_inchi

    mol_block, mol_count = args

    try:
        mol_name = mol_block.split('\n')[0].strip() or f'mol#{mol_count}'

        m_end_idx = mol_block.find('M  END')
        if m_end_idx < 0:
            return None

        mol_text = mol_block[:m_end_idx + 6]

        mol, num_heavy = parse_mol_block_robust(mol_text)
        if mol is None:
            return None

        mol_with_h = Chem.AddHs(mol)

        beyond_available = False
        try:
            if BEYOND_AVAILABLE:
                ensure_gasteiger_charges(mol_with_h)
                beyond_available = True
        except:
            pass

        inchikey = None
        canon_map = {}
        try:
            inchikey = rdkit_inchi.MolToInchiKey(mol)
            if inchikey:
                canon_smiles = Chem.MolToSmiles(mol)
                mol_canon_no_h = Chem.MolFromSmiles(canon_smiles)
                if mol_canon_no_h is not None:
                    match_no_h = _substruct_match_with_timeout(mol, mol_canon_no_h, timeout=1.0)
                    for canon_idx, sdf_idx in enumerate(match_no_h):
                        canon_map[sdf_idx] = canon_idx

                    mol_canon_with_h = Chem.AddHs(mol_canon_no_h)
                    match_with_h = _substruct_match_with_timeout(mol_with_h, mol_canon_with_h, timeout=1.0)
                    if match_with_h:
                        for canon_idx, sdf_idx in enumerate(match_with_h):
                            canon_map[sdf_idx] = canon_idx
        except Exception:
            inchikey = None

        solvent_match = re.search(r'<Solvent>\s*\n([^\n]+)', mol_block)
        solvent_map = {}
        if solvent_match:
            solvent_line = solvent_match.group(1).strip()
            parts = re.findall(r'(\d+):(.+?)(?=\s+\d+:|$)', solvent_line)
            for idx, solv in parts:
                solvent_map[idx] = normalize_solvent(solv.strip())

        fp_records = []
        exact_records = []
        composite_records = []
        count = 0

        for nucleus, pattern in [('13C', r'<Spectrum 13C (\d+)>\s*\n([^\n<]+)'),
                                 ('1H', r'<Spectrum 1H (\d+)>\s*\n([^\n<]+)')]:
            for spec_match in re.finditer(pattern, mol_block):
                spec_idx = spec_match.group(1)
                spec_data = spec_match.group(2).strip()
                solvent = solvent_map.get(spec_idx, 'Unknown')
                solvent_normalized = normalize_solvent(solvent)

                if not is_valid_solvent(solvent_normalized):
                    continue

                for entry in spec_data.split('|'):
                    parts = entry.split(';')
                    if len(parts) >= 3:
                        try:
                            shift = float(parts[0])
                            atom_idx_no_h = int(parts[2])

                            if atom_idx_no_h >= mol.GetNumAtoms():
                                continue

                            if nucleus == '1H':
                                if shift < -5.0 or shift > 14.0:
                                    continue
                            elif nucleus == '13C':
                                if shift < -5.0 or shift > 250.0:
                                    continue

                            if nucleus == '1H':

                                atom = mol_with_h.GetAtomWithIdx(atom_idx_no_h)
                                if atom.GetSymbol() == 'H':
                                    h_indices = [atom_idx_no_h]
                                else:
                                    h_indices = [n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == 'H']

                                for h_idx in h_indices:
                                    env_str = get_atom_env_string(mol_with_h, h_idx)
                                    for radius in SUPPORTED_RADII:
                                        fp = get_atom_fingerprint(mol_with_h, h_idx, radius)
                                        fp_records.append((radius, fp, solvent_normalized, nucleus, shift, env_str))

                                        if beyond_available:
                                            try:
                                                beyond = compute_beyond_discrete_safe(mol_with_h, h_idx, include_chiral=False)
                                                composite_fp = (fp, beyond)
                                                composite_records.append((radius, composite_fp, solvent_normalized, nucleus, shift))
                                            except Exception as e:

                                                logger.debug(f"Beyond features failed for H{h_idx}: {e}")
                                                pass

                                    if inchikey and h_idx in canon_map:
                                        exact_records.append((inchikey, solvent_normalized, nucleus, canon_map[h_idx], shift))
                                    count += 1
                            else:

                                env_str = get_atom_env_string(mol_with_h, atom_idx_no_h)
                                for radius in SUPPORTED_RADII:
                                    fp = get_atom_fingerprint(mol_with_h, atom_idx_no_h, radius)
                                    fp_records.append((radius, fp, solvent_normalized, nucleus, shift, env_str))

                                    if beyond_available:
                                        try:
                                            beyond = compute_beyond_discrete_safe(mol_with_h, atom_idx_no_h, include_chiral=False)
                                            composite_fp = (fp, beyond)
                                            composite_records.append((radius, composite_fp, solvent_normalized, nucleus, shift))
                                        except:
                                            pass

                                if inchikey and atom_idx_no_h in canon_map:
                                    exact_records.append((inchikey, solvent_normalized, nucleus, canon_map[atom_idx_no_h], shift))
                                count += 1
                        except (ValueError, IndexError):
                            pass

        return (fp_records, exact_records, composite_records, count)

    except Exception:
        return None

def _process_nmr_spectra_worker(args):
    """Process one NMR spectra record."""
    from rdkit.Chem import inchi as rdkit_inchi

    mol_block, mol_count = args

    try:
        m_end_idx = mol_block.find('M  END')
        if m_end_idx < 0:
            return None

        mol_text = mol_block[:m_end_idx + 6]
        mol, num_heavy = parse_mol_block_robust(mol_text)
        if mol is None:
            return None

        mol_with_h = Chem.AddHs(mol)

        beyond_available = False
        try:
            if BEYOND_AVAILABLE:
                ensure_gasteiger_charges(mol_with_h)
                beyond_available = True
        except:
            pass

        inchikey = None
        canon_map = {}
        try:
            inchikey = rdkit_inchi.MolToInchiKey(mol)
            if inchikey:
                canon_smiles = Chem.MolToSmiles(mol)
                mol_canon_no_h = Chem.MolFromSmiles(canon_smiles)
                if mol_canon_no_h is not None:
                    match_no_h = _substruct_match_with_timeout(mol, mol_canon_no_h, timeout=1.0)
                    for canon_idx, sdf_idx in enumerate(match_no_h):
                        canon_map[sdf_idx] = canon_idx
                    mol_canon_with_h = Chem.AddHs(mol_canon_no_h)
                    match_with_h = _substruct_match_with_timeout(mol_with_h, mol_canon_with_h, timeout=1.0)
                    if match_with_h:
                        for canon_idx, sdf_idx in enumerate(match_with_h):
                            canon_map[sdf_idx] = canon_idx
        except Exception:
            inchikey = None

        solvent_match = re.search(r'<Solvent>\s*\(\d+\)\s*\n([^\n]+)', mol_block)
        solvent_map = {}
        if solvent_match:
            solvent_line = solvent_match.group(1).strip()
            parts = re.findall(r'(\d+):(.+?)(?=\s*\d+:|$)', solvent_line)
            for idx, solv in parts:
                solvent_map[idx] = normalize_solvent(solv.strip())
            if not solvent_map:
                solvent_map['0'] = normalize_solvent(solvent_line)

        fp_records = []
        exact_records = []
        composite_records = []
        count = 0

        for nucleus, pattern in [('13C', r'<Spectrum 13C (\d+)>\s*\(\d+\)\s*\n([^\n<]+)'),
                                 ('1H', r'<Spectrum 1H (\d+)>\s*\(\d+\)\s*\n([^\n<]+)')]:
            for spec_match in re.finditer(pattern, mol_block):
                spec_idx = spec_match.group(1)
                spec_data = spec_match.group(2).strip()
                solvent = solvent_map.get(spec_idx, solvent_map.get('0', 'Unknown'))
                solvent_normalized = normalize_solvent(solvent)

                if not is_valid_solvent(solvent_normalized):
                    continue

                for entry in spec_data.split('|'):
                    if not entry.strip():
                        continue
                    parts = entry.split(';')
                    if len(parts) >= 3:
                        try:
                            shift = float(parts[0])
                            atom_idx = int(parts[2])

                            if atom_idx >= mol_with_h.GetNumAtoms():
                                continue

                            if nucleus == '1H' and (shift < -5.0 or shift > 14.0):
                                continue
                            if nucleus == '13C' and (shift < -5.0 or shift > 250.0):
                                continue

                            env_str = get_atom_env_string(mol_with_h, atom_idx)
                            for radius in SUPPORTED_RADII:
                                fp = get_atom_fingerprint(mol_with_h, atom_idx, radius)
                                fp_records.append((radius, fp, solvent_normalized, nucleus, shift, env_str))

                                if beyond_available:
                                    try:
                                        beyond = compute_beyond_discrete_safe(mol_with_h, atom_idx, include_chiral=False)
                                        composite_fp = (fp, beyond)
                                        composite_records.append((radius, composite_fp, solvent_normalized, nucleus, shift))
                                    except:
                                        pass

                            if inchikey and atom_idx in canon_map:
                                exact_records.append((inchikey, solvent_normalized, nucleus, canon_map[atom_idx], shift))
                            count += 1
                        except (ValueError, IndexError):
                            pass

        return (fp_records, exact_records, composite_records, count)

    except Exception:
        return None

def _process_npmrd_file_worker(args):
    """Process one NPMRD file."""
    from rdkit.Chem import inchi as rdkit_inchi

    filepath, file_idx = args

    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        m_end_idx = content.find('M  END')
        if m_end_idx < 0:
            return None

        mol_text = content[:m_end_idx + 6]
        mol, num_heavy = parse_mol_block_robust(mol_text)
        if mol is None:
            return None

        mol_with_h = Chem.AddHs(mol)

        beyond_available = False
        try:
            if BEYOND_AVAILABLE:
                ensure_gasteiger_charges(mol_with_h)
                beyond_available = True
        except:
            pass

        inchikey = None
        canon_map = {}
        try:
            inchikey = rdkit_inchi.MolToInchiKey(mol)
            if inchikey:
                canon_smiles = Chem.MolToSmiles(mol)
                mol_canon_no_h = Chem.MolFromSmiles(canon_smiles)
                if mol_canon_no_h is not None:
                    match_no_h = _substruct_match_with_timeout(mol, mol_canon_no_h, timeout=1.0)
                    for canon_idx, sdf_idx in enumerate(match_no_h):
                        canon_map[sdf_idx] = canon_idx
                    mol_canon_with_h = Chem.AddHs(mol_canon_no_h)
                    match_with_h = _substruct_match_with_timeout(mol_with_h, mol_canon_with_h, timeout=1.0)
                    if match_with_h:
                        for canon_idx, sdf_idx in enumerate(match_with_h):
                            canon_map[sdf_idx] = canon_idx
        except Exception:
            inchikey = None

        solvent_match = re.search(r'> <Solvent>\s*\n([^\n]+)', content)
        solvent = normalize_solvent(solvent_match.group(1).strip()) if solvent_match else 'Unknown'

        if not is_valid_solvent(solvent):
            solvent = 'Unknown'

        fp_records = []
        exact_records = []
        composite_records = []
        count = 0

        spec_13c_match = re.search(r'> <Spectrum_13C>\s*\n([^\n]+)', content)
        if spec_13c_match:
            spec_data = spec_13c_match.group(1).strip()
            for entry in spec_data.split('|'):
                if not entry.strip():
                    continue
                parts = entry.split(';')
                if len(parts) >= 3:
                    try:
                        shift = float(parts[0])

                        atom_idx = int(parts[2]) - 1
                        if atom_idx < 0 or atom_idx >= mol_with_h.GetNumAtoms():
                            continue
                        if shift < -5.0 or shift > 250.0:
                            continue

                        env_str = get_atom_env_string(mol_with_h, atom_idx)
                        for radius in SUPPORTED_RADII:
                            fp = get_atom_fingerprint(mol_with_h, atom_idx, radius)
                            fp_records.append((radius, fp, solvent, '13C', shift, env_str))

                            if beyond_available:
                                try:
                                    beyond = compute_beyond_discrete_safe(mol_with_h, atom_idx, include_chiral=False)
                                    composite_fp = (fp, beyond)
                                    composite_records.append((radius, composite_fp, solvent, '13C', shift))
                                except:
                                    pass

                        if inchikey and atom_idx in canon_map:
                            exact_records.append((inchikey, solvent, '13C', canon_map[atom_idx], shift))
                        count += 1
                    except (ValueError, IndexError):
                        pass

        spec_1h_match = re.search(r'> <Spectrum_1H>\s*\n([^\n]+)', content)
        if not spec_1h_match:
            spec_1h_match = re.search(r'> <Spectrum_Unknown>\s*\n([^\n]+)', content)

        if spec_1h_match:
            spec_data = spec_1h_match.group(1).strip()
            for entry in spec_data.split('|'):
                if not entry.strip():
                    continue
                parts = entry.split(';')
                if len(parts) >= 3:
                    try:
                        shift = float(parts[0])

                        for idx_str in parts[2:]:
                            atom_idx = int(idx_str) - 1
                            if atom_idx < 0 or atom_idx >= mol_with_h.GetNumAtoms():
                                continue
                            if shift < -5.0 or shift > 14.0:
                                continue

                            env_str = get_atom_env_string(mol_with_h, atom_idx)
                            for radius in SUPPORTED_RADII:
                                fp = get_atom_fingerprint(mol_with_h, atom_idx, radius)
                                fp_records.append((radius, fp, solvent, '1H', shift, env_str))

                                if beyond_available:
                                    try:
                                        beyond = compute_beyond_discrete_safe(mol_with_h, atom_idx, include_chiral=False)
                                        composite_fp = (fp, beyond)
                                        composite_records.append((radius, composite_fp, solvent, '1H', shift))
                                    except:
                                        pass

                            if inchikey and atom_idx in canon_map:
                                exact_records.append((inchikey, solvent, '1H', canon_map[atom_idx], shift))
                            count += 1
                    except (ValueError, IndexError):
                        pass

        return (fp_records, exact_records, composite_records, count)

    except Exception:
        return None

def _process_experimental_file_worker(args):
    """Process one experimental data file."""
    from rdkit.Chem import inchi as rdkit_inchi

    filepath, file_idx = args

    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        m_end_idx = content.find('M  END')
        if m_end_idx < 0:
            return None

        mol_text = content[:m_end_idx + 6]
        mol = Chem.MolFromMolBlock(mol_text)
        if mol is None:
            return None

        num_heavy = mol.GetNumAtoms()
        mol_with_h = Chem.AddHs(mol)

        beyond_available = False
        try:
            if BEYOND_AVAILABLE:
                ensure_gasteiger_charges(mol_with_h)
                beyond_available = True
        except:
            pass

        inchikey = None
        canon_map = {}
        try:
            inchikey = rdkit_inchi.MolToInchiKey(mol)
            if inchikey:
                canon_smiles = Chem.MolToSmiles(mol)
                mol_canon_no_h = Chem.MolFromSmiles(canon_smiles)
                if mol_canon_no_h is not None:
                    match_no_h = _substruct_match_with_timeout(mol, mol_canon_no_h, timeout=1.0)
                    for canon_idx, sdf_idx in enumerate(match_no_h):
                        canon_map[sdf_idx] = canon_idx
                    mol_canon_with_h = Chem.AddHs(mol_canon_no_h)
                    match_with_h = _substruct_match_with_timeout(mol_with_h, mol_canon_with_h, timeout=1.0)
                    if match_with_h:
                        for canon_idx, sdf_idx in enumerate(match_with_h):
                            canon_map[sdf_idx] = canon_idx
        except Exception:
            inchikey = None

        solvent_match = re.search(r'<Solvent>\s*\(\d+\)\s*\n([^\n]+)', content)
        solvent = normalize_solvent(solvent_match.group(1).strip()) if solvent_match else 'Unknown'

        fp_records = []
        exact_records = []
        composite_records = []
        count = 0

        spec_match = re.search(r'<Spectrum 1H>\s*\(\d+\)\s*\n([^\n]+)', content)
        if spec_match:
            spec_data = spec_match.group(1).strip()
            for entry in spec_data.split('|'):
                if not entry.strip():
                    continue
                parts = entry.split(';')
                if len(parts) >= 3:
                    try:
                        shift = float(parts[0])
                        atom_idx = int(parts[2])
                        if atom_idx >= mol_with_h.GetNumAtoms():
                            continue
                        if shift < -5.0 or shift > 14.0:
                            continue

                        atom = mol_with_h.GetAtomWithIdx(atom_idx)
                        if atom.GetSymbol() == 'H':
                            h_indices = [atom_idx]
                        else:
                            h_indices = [n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == 'H']

                        for h_idx in h_indices:
                            env_str = get_atom_env_string(mol_with_h, h_idx)
                            for radius in SUPPORTED_RADII:
                                fp = get_atom_fingerprint(mol_with_h, h_idx, radius)
                                fp_records.append((radius, fp, solvent, '1H', shift, env_str))

                                if beyond_available:
                                    try:
                                        beyond = compute_beyond_discrete_safe(mol_with_h, h_idx, include_chiral=False)
                                        composite_fp = (fp, beyond)
                                        composite_records.append((radius, composite_fp, solvent, '1H', shift))
                                    except:
                                        pass

                            if inchikey and h_idx in canon_map:
                                exact_records.append((inchikey, solvent, '1H', canon_map[h_idx], shift))
                            count += 1
                    except (ValueError, IndexError):
                        pass

        spec_13c_match = re.search(r'<Spectrum 13C>\s*\(\d+\)\s*\n([^\n]+)', content)
        if spec_13c_match:
            spec_data_13c = spec_13c_match.group(1).strip()
            for entry in spec_data_13c.split('|'):
                if not entry.strip():
                    continue
                parts = entry.split(';')
                if len(parts) >= 3:
                    try:
                        shift = float(parts[0])
                        atom_idx = int(parts[2])
                        if atom_idx >= mol_with_h.GetNumAtoms():
                            continue
                        if shift < -5.0 or shift > 250.0:
                            continue

                        env_str = get_atom_env_string(mol_with_h, atom_idx)
                        for radius in SUPPORTED_RADII:
                            fp = get_atom_fingerprint(mol_with_h, atom_idx, radius)
                            fp_records.append((radius, fp, solvent, '13C', shift, env_str))

                            if beyond_available:
                                try:
                                    beyond = compute_beyond_discrete_safe(mol_with_h, atom_idx, include_chiral=False)
                                    composite_fp = (fp, beyond)
                                    composite_records.append((radius, composite_fp, solvent, '13C', shift))
                                except:
                                    pass

                        if inchikey and atom_idx in canon_map:
                            exact_records.append((inchikey, solvent, '13C', canon_map[atom_idx], shift))
                        count += 1
                    except (ValueError, IndexError):
                        pass

        return (fp_records, exact_records, composite_records, count)

    except Exception:
        return None

class FingerprintPredictor:
    """Morgan fingerprint-based NMR shift predictor."""

    _instance = None
    _db = None
    _db_composite = None
    _db_mtime = None
    _exact_db = None
    _exact_1h_groups = None

    def __new__(cls):
        """Return the singleton instance."""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):

        needs_reload = False

        if FingerprintPredictor._db is None:
            needs_reload = True
        elif os.path.exists(DB_CACHE_FILE):

            current_mtime = os.path.getmtime(DB_CACHE_FILE)
            if FingerprintPredictor._db_mtime is None or current_mtime > FingerprintPredictor._db_mtime:
                print(f"[FingerprintPredictor] Database file updated, reloading...", flush=True)
                needs_reload = True

        if needs_reload:
            self._load_database()

    def _load_database(self):
        """Load the fingerprint database cache."""
        if not os.path.exists(DB_CACHE_FILE):
            print(f"[FingerprintPredictor] Cache not found: {DB_CACHE_FILE}")
            print("[FingerprintPredictor] Building database...")
            self._build_database()
            return

        try:
            with open(DB_CACHE_FILE, 'rb') as f:
                data = pickle.load(f)

            if data.get('version') not in ('14.0',):
                print(f"[FingerprintPredictor] Cache version mismatch ({data.get('version')} != 14.0), rebuilding...")
                self._build_database()
                return

            repaired = self._repair_bmrb_ambiguities(data)

            FingerprintPredictor._db = {}
            for radius in SUPPORTED_RADII:
                FingerprintPredictor._db[radius] = data['indexes'].get(radius, {})
            FingerprintPredictor._fp_to_env = data.get('fp_to_env', {})
            FingerprintPredictor._exact_db = data.get('exact_db', {})
            FingerprintPredictor._exact_1h_groups = data.get('exact_1h_groups', {})

            FingerprintPredictor._db_composite = {}
            composite_indexes = data.get('indexes_composite', {})
            for radius in SUPPORTED_RADII:
                FingerprintPredictor._db_composite[radius] = composite_indexes.get(radius, {})

            if RADIUS_ULTRA_EXACT not in data['indexes'] or len(data['indexes'].get(RADIUS_ULTRA_EXACT, {})) == 0:
                print(f"[FingerprintPredictor] Radius {RADIUS_ULTRA_EXACT} not in cache, rebuilding...")
                self._build_database()
                return

            for r in SUPPORTED_RADII:
                count = len(FingerprintPredictor._db.get(r, {}))
                print(f"[FingerprintPredictor] Loaded {count} fingerprints (radius {r})")

            comp_total = sum(len(FingerprintPredictor._db_composite.get(r, {})) for r in SUPPORTED_RADII)
            if comp_total > 0:
                print(f"[FingerprintPredictor] Loaded {comp_total} composite fingerprints total")

            FingerprintPredictor._db_mtime = os.path.getmtime(DB_CACHE_FILE)
            if repaired:
                self._save_cache()
                FingerprintPredictor._db_mtime = os.path.getmtime(DB_CACHE_FILE)

        except Exception as e:
            print(f"[FingerprintPredictor] Error loading cache: {e}")
            self._build_database()

    def _build_database(self):
        """Build the fingerprint database."""
        import time
        start = time.time()

        db_by_radius = {r: defaultdict(lambda: defaultdict(lambda: defaultdict(list))) for r in SUPPORTED_RADII}
        fp_to_env = {}

        db_composite_by_radius = {r: defaultdict(lambda: defaultdict(lambda: defaultdict(list))) for r in SUPPORTED_RADII}

        exact_db = {}
        FingerprintPredictor._exact_db = exact_db
        exact_1h_groups = {}
        FingerprintPredictor._exact_1h_groups = exact_1h_groups

        count_bmrb = self._load_bmrb(
            db_by_radius, fp_to_env, exact_db,
            db_composite=db_composite_by_radius,
            exact_1h_groups=exact_1h_groups,
        )

        count_nmrdb = self._load_nmrshiftdb(db_by_radius, fp_to_env, exact_db=exact_db, db_composite=db_composite_by_radius)

        count_spectra = self._load_nmr_spectra(db_by_radius, fp_to_env, exact_db=exact_db, db_composite=db_composite_by_radius)

        count_npmrd = self._load_npmrd(db_by_radius, fp_to_env, exact_db=exact_db, db_composite=db_composite_by_radius)

        count_experimental = self._load_experimental_data(db_by_radius, fp_to_env, exact_db, db_composite=db_composite_by_radius)

        elapsed = time.time() - start
        total = count_bmrb + count_nmrdb + count_spectra + count_npmrd + count_experimental
        print(f"[FingerprintPredictor] Built database: {total} shifts in {elapsed:.1f}s")
        for r in SUPPORTED_RADII:
            print(f"[FingerprintPredictor] Radius {r}: {len(db_by_radius[r])} unique fingerprints")

        if BEYOND_AVAILABLE:
            for r in SUPPORTED_RADII:
                comp_count = len(db_composite_by_radius[r])
                if comp_count > 0:
                    print(f"[FingerprintPredictor] Composite R{r}: {comp_count} unique fingerprints")

        FingerprintPredictor._db = {r: dict(db_by_radius[r]) for r in SUPPORTED_RADII}
        FingerprintPredictor._db_composite = {r: dict(db_composite_by_radius[r]) for r in SUPPORTED_RADII}
        FingerprintPredictor._fp_to_env = fp_to_env

        self._save_cache()

    def _repair_bmrb_ambiguities(self, data):
        """Repair reviewed legacy BMRB ambiguity records."""
        if data.get('bmrb_ambiguity_revision', 0) >= REVISION:
            return False
        from collections import Counter
        contributions = {r: defaultdict(lambda: defaultdict(lambda: defaultdict(list))) for r in SUPPORTED_RADII}
        composite = {r: defaultdict(lambda: defaultdict(lambda: defaultdict(list))) for r in SUPPORTED_RADII}
        found = set()
        for filename in sorted(os.listdir(BMRB_DIR)):
            if not filename.endswith('.sdf'):
                continue
            for mol in Chem.SDMolSupplier(os.path.join(BMRB_DIR, filename), removeHs=False):
                if mol is None or not mol.HasProp('BMRB_ID'):
                    continue
                entry = mol.GetProp('BMRB_ID').strip()
                if entry not in KNOWN_IDS:
                    continue
                if entry in found or Chem.MolToInchiKey(Chem.RemoveHs(mol)) != INCHIKEY:
                    raise ValueError('Unexpected BMRB source for ambiguity repair')
                found.add(entry)
                ensure_gasteiger_charges(mol)
                for a in reviewed_assignments():
                    self._add_record(contributions, {}, mol, a['atom_index'], a['legacy_shift'],
                                     'D2O', a['nucleus'], db_composite_by_radius=composite)
        if found != KNOWN_IDS:
            raise ValueError('Reviewed BMRB sources missing; cache needs rebuilding')
        changes = []
        for name, removed in [('indexes', contributions), ('indexes_composite', composite)]:
            if name == 'indexes_composite' and data.get(name) and not BEYOND_AVAILABLE:
                raise ValueError('Beyond features required to repair composite cache')
            for radius, fps in removed.items():
                for fp, solvents in fps.items():
                    for solvent, nuclei in solvents.items():
                        for nucleus, shifts in nuclei.items():
                            bucket = data.get(name, {}).get(radius, {}).get(fp, {}).get(solvent, {})
                            existing = bucket.get(nucleus, [])
                            counts = Counter(shifts)
                            if Counter(existing) & counts != counts:
                                raise ValueError('Legacy BMRB contributions not found; rebuild required')
                            kept = []
                            for value in existing:
                                if counts[value]:
                                    counts[value] -= 1
                                else:
                                    kept.append(value)
                            changes.append((bucket, nucleus, kept))

        exact = data.get('exact_db', {}).get(INCHIKEY, {}).get('D2O', {})
        for bucket, nucleus, kept in changes:
            bucket[nucleus] = kept
        for nucleus, invalid in [('1H', 2.428), ('13C', 33.325)]:
            values = exact.get(nucleus, {})
            for index in list(values):
                if values[index] == invalid:
                    del values[index]
        data['bmrb_ambiguity_revision'] = REVISION
        print('[FingerprintPredictor] Removed ambiguous bmse000064/bmse000937 assignments', flush=True)
        return True

    def _save_cache(self):
        """Save the fingerprint database cache."""
        try:

            indexes_serializable = {}
            stats = {}
            for radius in SUPPORTED_RADII:
                indexes_serializable[radius] = {}
                db = FingerprintPredictor._db.get(radius, {})
                for fp, solvent_dict in db.items():
                    indexes_serializable[radius][fp] = {}
                    for solvent, nuc_dict in solvent_dict.items():
                        indexes_serializable[radius][fp][solvent] = dict(nuc_dict)
                stats[radius] = len(db)

            composite_serializable = {}
            composite_stats = {}
            for radius in SUPPORTED_RADII:
                composite_serializable[radius] = {}
                db_comp = (FingerprintPredictor._db_composite or {}).get(radius, {})
                for fp, solvent_dict in db_comp.items():
                    composite_serializable[radius][fp] = {}
                    for solvent, nuc_dict in solvent_dict.items():
                        composite_serializable[radius][fp][solvent] = dict(nuc_dict)
                composite_stats[radius] = len(db_comp)

            data = {
                'bmrb_ambiguity_revision': REVISION,
                'version': '14.0',
                'indexes': indexes_serializable,
                'indexes_composite': composite_serializable,
                'fp_to_env': FingerprintPredictor._fp_to_env,
                'exact_db': FingerprintPredictor._exact_db or {},
                'exact_1h_groups': FingerprintPredictor._exact_1h_groups or {},
                'total_records': sum(stats.values()),
                'stats': {'fps_by_radius': stats, 'composite_by_radius': composite_stats}
            }

            os.makedirs(os.path.dirname(DB_CACHE_FILE), exist_ok=True)
            temporary_cache = DB_CACHE_FILE + '.tmp'
            with open(temporary_cache, 'wb') as f:
                pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
            os.replace(temporary_cache, DB_CACHE_FILE)

            print(f"[FingerprintPredictor] Cache saved to {DB_CACHE_FILE}")
        except Exception as e:
            print(f"[FingerprintPredictor] Could not save cache: {e}")

    def _add_record(self, db_by_radius, fp_to_env, mol, atom_idx, shift, solvent, nucleus, env_str=None, db_composite_by_radius=None):
        """Add one fingerprint record."""
        if env_str is None:
            env_str = get_atom_env_string(mol, atom_idx)

        solvent_normalized = normalize_solvent(solvent)

        if not is_valid_solvent(solvent_normalized):
            return

        if nucleus == '1H':
            if shift < -5.0 or shift > 14.0:
                return
        elif nucleus == '13C':
            if shift < -5.0 or shift > 250.0:
                return

        beyond = None
        if db_composite_by_radius is not None and BEYOND_AVAILABLE:
            try:

                beyond = compute_beyond_discrete_safe(mol, atom_idx, include_chiral=False)
            except:
                pass

        for radius in SUPPORTED_RADII:
            fp = get_atom_fingerprint(mol, atom_idx, radius)
            db_by_radius[radius][fp][solvent_normalized][nucleus].append(shift)
            if fp not in fp_to_env:
                fp_to_env[fp] = f"r{radius}|{env_str}"

            if beyond is not None and db_composite_by_radius is not None:
                composite_fp = (fp, beyond)
                db_composite_by_radius[radius][composite_fp][solvent_normalized][nucleus].append(shift)

    def _store_1h_shift(self, mol_with_h, atom_idx, shift, solvent, inchikey,
                        canon_map, exact_db, db, fp_to_env, db_composite_by_radius=None):
        """Store one proton shift assignment."""
        atom = mol_with_h.GetAtomWithIdx(atom_idx)
        count = 0

        if atom.GetSymbol() == 'H':

            self._add_record(db, fp_to_env, mol_with_h, atom_idx, shift, solvent, '1H',
                           db_composite_by_radius=db_composite_by_radius)
            count += 1

            if inchikey and exact_db is not None and atom_idx in canon_map:
                h_canon_idx = canon_map[atom_idx]
                exact_db.setdefault(inchikey, {}).setdefault(solvent, {}).setdefault('1H', {})[h_canon_idx] = shift

        else:

            h_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'H']

            for h_atom in h_neighbors:
                h_idx = h_atom.GetIdx()

                self._add_record(db, fp_to_env, mol_with_h, h_idx, shift, solvent, '1H',
                               db_composite_by_radius=db_composite_by_radius)
                count += 1

                if inchikey and exact_db is not None and h_idx in canon_map:
                    h_canon_idx = canon_map[h_idx]
                    exact_db.setdefault(inchikey, {}).setdefault(solvent, {}).setdefault('1H', {})[h_canon_idx] = shift

        return count

    def _load_bmrb(self, db, fp_to_env, exact_db=None, db_composite=None,
                   exact_1h_groups=None) -> int:
        """Load BMRB records."""
        from rdkit.Chem import inchi as rdkit_inchi

        if not os.path.exists(BMRB_DIR):
            return 0

        sdf_files = [f for f in os.listdir(BMRB_DIR) if f.endswith('.sdf')]
        if not sdf_files:
            print("[FingerprintPredictor] BMRB: No SDF file found")
            return 0

        bmrb_file = os.path.join(BMRB_DIR, sdf_files[0])
        count = 0
        exact_count = 0

        try:

            suppl = Chem.SDMolSupplier(bmrb_file, removeHs=False, sanitize=True)

            for mol_idx, mol in enumerate(suppl):
                if mol is None:
                    continue

                if mol_idx % 100 == 0:
                    mol_name = mol.GetProp('_Name') if mol.HasProp('_Name') else f'mol#{mol_idx}'
                    print(f"  [BMRB] mol {mol_idx}: {mol_name}", flush=True)

                try:

                    raw_solvent = mol.GetProp('Solvent') if mol.HasProp('Solvent') else 'Unknown'
                    solvent = normalize_solvent(raw_solvent)

                    if not mol.HasProp('NMR_Shifts'):
                        continue

                    nmr_shifts_text = mol.GetProp('NMR_Shifts').strip()
                    shift_lines = [line.strip() for line in nmr_shifts_text.split('\n') if line.strip()]
                    unresolved = ambiguous_atoms(mol)

                    if mol.GetNumAtoms() > 0 and all(atom.GetSymbol() != 'H' for atom in mol.GetAtoms()):
                        mol = Chem.AddHs(mol)

                    if db_composite is not None and BEYOND_AVAILABLE:
                        ensure_gasteiger_charges(mol)

                    inchikey = None
                    canon_map = {}
                    if exact_db is not None:
                        try:
                            mol_no_h = Chem.RemoveHs(mol)
                            inchikey = rdkit_inchi.MolToInchiKey(mol_no_h)
                            if inchikey:

                                canon_smiles = Chem.MolToSmiles(mol_no_h)
                                mol_canon_no_h = Chem.MolFromSmiles(canon_smiles)
                                if mol_canon_no_h is not None:
                                    match_no_h = _substruct_match_with_timeout(mol_no_h, mol_canon_no_h, timeout=1.0)

                                    noh_to_canon = {sdf_i: canon_i for canon_i, sdf_i in enumerate(match_no_h)}

                                    for noh_idx, canon_idx in noh_to_canon.items():
                                        canon_map[noh_idx] = canon_idx

                                    mol_canon_with_h = Chem.AddHs(mol_canon_no_h)
                                    mol_with_h_for_match = Chem.AddHs(mol)
                                    match_with_h = _substruct_match_with_timeout(mol_with_h_for_match, mol_canon_with_h, timeout=1.0)
                                    if match_with_h:

                                        for canon_idx, sdf_idx in enumerate(match_with_h):
                                            canon_map[sdf_idx] = canon_idx
                        except Exception:
                            inchikey = None

                    proton_groups = defaultdict(lambda: defaultdict(list))

                    for line in shift_lines:
                        parts = line.split()
                        if len(parts) != 4:
                            continue

                        try:
                            atom_idx = int(parts[0])
                            atom_symbol = parts[1]
                            shift_1h = float(parts[2])
                            shift_13c = float(parts[3])
                        except (ValueError, IndexError):
                            continue

                        if atom_idx >= mol.GetNumAtoms():
                            continue

                        if shift_1h != 0.0 and atom_symbol == 'H' and (atom_idx, '1H') not in unresolved:

                            h_atom = mol.GetAtomWithIdx(atom_idx)
                            heavy_neighbors = [n for n in h_atom.GetNeighbors() if n.GetAtomicNum() > 1]
                            if inchikey and heavy_neighbors:
                                heavy_idx = heavy_neighbors[0].GetIdx()
                                heavy_key = str(get_atom_fingerprint(
                                    mol, heavy_idx, RADIUS_ULTRA_EXACT, heavy_only=True
                                ))
                                proton_groups[heavy_key][heavy_idx].append(float(shift_1h))

                            added = self._store_1h_shift(mol, atom_idx, shift_1h, solvent,
                                                         inchikey, canon_map, exact_db, db, fp_to_env,
                                                         db_composite_by_radius=db_composite)
                            count += added
                            if added > 0:
                                exact_count += added

                        if shift_13c != 0.0 and atom_symbol == 'C' and (atom_idx, '13C') not in unresolved:
                            self._add_record(db, fp_to_env, mol, atom_idx, shift_13c, solvent, '13C',
                                           db_composite_by_radius=db_composite)
                            count += 1
                            if inchikey and exact_db is not None and atom_idx in canon_map:
                                canon_idx = canon_map[atom_idx]
                                exact_db.setdefault(inchikey, {}).setdefault(solvent, {}).setdefault('13C', {})[canon_idx] = shift_13c
                                exact_count += 1

                    if inchikey and exact_1h_groups is not None and proton_groups:
                        solvent_groups = exact_1h_groups.setdefault(inchikey, {}).setdefault(solvent, {})
                        for heavy_key, instances in proton_groups.items():
                            records = solvent_groups.setdefault(heavy_key, [])
                            for heavy_idx, shifts in instances.items():
                                attached_h = sum(
                                    1 for n in mol.GetAtomWithIdx(heavy_idx).GetNeighbors()
                                    if n.GetAtomicNum() == 1
                                )
                                record = {
                                    'n_h': attached_h,
                                    'shifts': sorted(float(v) for v in shifts),
                                }
                                # Preserve multiplicity of structurally identical
                                # groups (for example the two methyl groups of a
                                # symmetric molecule).  They are matched as a
                                # multiset at prediction time.
                                records.append(record)

                except Exception:

                    continue

        except Exception:

            pass

        print(f"[FingerprintPredictor] BMRB: {count} shifts, {exact_count} exact-match entries ({len(exact_db) if exact_db else 0} molecules total in exact_db)")
        return count

    def _load_nmrshiftdb(self, db, fp_to_env, exact_db=None, db_composite=None) -> int:
        """Load NMRShiftDB records."""
        if not os.path.exists(NMRSHIFTDB_FILE):
            return 0

        def _iter_mol_blocks(filepath):
            current_mol = []
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if line.strip() == '$$$$':
                        if current_mol:
                            yield '\n'.join(current_mol)
                            current_mol = []
                    else:
                        current_mol.append(line.rstrip())

        import time as _time
        _t0 = _time.time()

        print("  [NMRShiftDB] Loading mol blocks...", flush=True)
        mol_blocks = [(block, i+1) for i, block in enumerate(_iter_mol_blocks(NMRSHIFTDB_FILE))]
        total_mols = len(mol_blocks)
        print(f"  [NMRShiftDB] {total_mols} molecules loaded in {_time.time() - _t0:.1f}s", flush=True)

        n_workers = max(1, multiprocessing.cpu_count() - 1)
        print(f"  [NMRShiftDB] Starting parallel processing with {n_workers} workers...", flush=True)

        count = 0
        exact_count = 0
        processed = 0
        batch_size = 1000

        _t_parallel = _time.time()
        with ProcessPoolExecutor(max_workers=n_workers) as executor:

            futures = {executor.submit(_process_mol_block_worker, args): args[1] for args in mol_blocks}

            for future in as_completed(futures):
                processed += 1

                if processed % batch_size == 0:
                    elapsed = _time.time() - _t_parallel
                    rate = processed / elapsed if elapsed > 0 else 0
                    eta = (total_mols - processed) / rate if rate > 0 else 0
                    print(f"  [NMRShiftDB] {processed}/{total_mols} | {count} shifts | {elapsed:.0f}s | {rate:.0f} mol/s | ETA {eta/60:.1f}min", flush=True)

                result = future.result()
                if result is None:
                    continue

                if len(result) == 4:
                    fp_records, exact_records, composite_records, rec_count = result
                else:

                    fp_records, exact_records, rec_count = result
                    composite_records = []

                count += rec_count

                for radius, fp, solvent, nucleus, shift, env_str in fp_records:
                    db[radius][fp][solvent][nucleus].append(shift)
                    if fp not in fp_to_env:
                        fp_to_env[fp] = f"r{radius}|{env_str}"

                if db_composite is not None and composite_records:
                    for radius, composite_fp, solvent, nucleus, shift in composite_records:
                        db_composite[radius][composite_fp][solvent][nucleus].append(shift)

                if exact_db is not None:
                    for inchikey, solvent, nucleus, canon_idx, shift in exact_records:
                        exact_db.setdefault(inchikey, {}).setdefault(solvent, {}).setdefault(nucleus, {})[canon_idx] = shift
                        exact_count += 1

        elapsed_total = _time.time() - _t0
        print(f"[FingerprintPredictor] NMRShiftDB: {count} shifts, {exact_count} exact-match entries ({len(exact_db) if exact_db else 0} molecules) in {elapsed_total:.1f}s ({n_workers} workers)")
        return count

    def _load_nmr_spectra(self, db, fp_to_env, exact_db=None, db_composite=None) -> int:
        """Load NMR spectra records."""
        if not os.path.exists(NMR_SPECTRA_FILE):
            return 0

        def _iter_mol_blocks(filepath):
            current_mol = []
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if line.strip() == '$$$$':
                        if current_mol:
                            yield '\n'.join(current_mol)
                            current_mol = []
                    else:
                        current_mol.append(line.rstrip())

        import time as _time
        _t0 = _time.time()

        print("  [NMRSpectra] Loading mol blocks...", flush=True)
        mol_blocks = [(block, i+1) for i, block in enumerate(_iter_mol_blocks(NMR_SPECTRA_FILE))]
        total_mols = len(mol_blocks)
        print(f"  [NMRSpectra] {total_mols} molecules loaded in {_time.time() - _t0:.1f}s", flush=True)

        n_workers = max(1, multiprocessing.cpu_count() - 1)
        print(f"  [NMRSpectra] Starting parallel processing with {n_workers} workers...", flush=True)

        count = 0
        exact_count = 0
        processed = 0
        batch_size = 1000

        _t_parallel = _time.time()
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(_process_nmr_spectra_worker, args): args[1] for args in mol_blocks}

            for future in as_completed(futures):
                processed += 1
                if processed % batch_size == 0:
                    elapsed = _time.time() - _t_parallel
                    rate = processed / elapsed if elapsed > 0 else 0
                    eta = (total_mols - processed) / rate if rate > 0 else 0
                    print(f"  [NMRSpectra] {processed}/{total_mols} | {count} shifts | {elapsed:.0f}s | {rate:.0f} mol/s | ETA {eta/60:.1f}min", flush=True)

                result = future.result()
                if result is None:
                    continue

                if len(result) == 4:
                    fp_records, exact_records, composite_records, rec_count = result
                else:
                    fp_records, exact_records, rec_count = result
                    composite_records = []

                count += rec_count

                for radius, fp, solvent, nucleus, shift, env_str in fp_records:
                    db[radius][fp][solvent][nucleus].append(shift)
                    if fp not in fp_to_env:
                        fp_to_env[fp] = f"r{radius}|{env_str}"

                if db_composite is not None and composite_records:
                    for radius, composite_fp, solvent, nucleus, shift in composite_records:
                        db_composite[radius][composite_fp][solvent][nucleus].append(shift)

                if exact_db is not None:
                    for inchikey, solvent, nucleus, canon_idx, shift in exact_records:
                        exact_db.setdefault(inchikey, {}).setdefault(solvent, {}).setdefault(nucleus, {})[canon_idx] = shift
                        exact_count += 1

        elapsed_total = _time.time() - _t0
        print(f"[FingerprintPredictor] NMRSpectra: {count} shifts, {exact_count} exact-match entries ({len(exact_db) if exact_db else 0} molecules) in {elapsed_total:.1f}s ({n_workers} workers)")
        return count

    def _load_npmrd(self, db, fp_to_env, exact_db=None, db_composite=None) -> int:
        """Load NPMRD records."""
        if not os.path.exists(NPMRD_SDF_DIR):
            return 0

        import glob
        import time as _time
        _t0 = _time.time()

        sdf_files = glob.glob(os.path.join(NPMRD_SDF_DIR, '*.sdf'))
        total_files = len(sdf_files)
        print(f"  [NPMRD] {total_files} files found", flush=True)

        n_workers = max(1, multiprocessing.cpu_count() - 1)
        print(f"  [NPMRD] Starting parallel processing with {n_workers} workers...", flush=True)

        count = 0
        exact_count = 0
        composite_count = 0
        processed = 0
        batch_size = 100

        file_args = [(filepath, i) for i, filepath in enumerate(sdf_files)]

        _t_parallel = _time.time()
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(_process_npmrd_file_worker, args): args[1] for args in file_args}

            for future in as_completed(futures):
                processed += 1
                if processed % batch_size == 0:
                    elapsed = _time.time() - _t_parallel
                    rate = processed / elapsed if elapsed > 0 else 0
                    eta = (total_files - processed) / rate if rate > 0 else 0
                    print(f"  [NPMRD] {processed}/{total_files} | {count} shifts | {elapsed:.0f}s | {rate:.0f} files/s | ETA {eta:.1f}s", flush=True)

                result = future.result()
                if result is None:
                    continue

                if len(result) == 4:
                    fp_records, exact_records, composite_records, rec_count = result
                else:
                    fp_records, exact_records, rec_count = result
                    composite_records = []
                count += rec_count

                for radius, fp, solvent, nucleus, shift, env_str in fp_records:
                    db[radius][fp][solvent][nucleus].append(shift)
                    if fp not in fp_to_env:
                        fp_to_env[fp] = f"r{radius}|{env_str}"

                if exact_db is not None:
                    for inchikey, solvent, nucleus, canon_idx, shift in exact_records:
                        exact_db.setdefault(inchikey, {}).setdefault(solvent, {}).setdefault(nucleus, {})[canon_idx] = shift
                        exact_count += 1

                if db_composite is not None and composite_records:
                    for radius, composite_fp, solvent, nucleus, shift in composite_records:
                        db_composite[radius][composite_fp][solvent][nucleus].append(shift)
                        composite_count += 1

        elapsed_total = _time.time() - _t0
        print(f"[FingerprintPredictor] NPMRD: {count} shifts from {total_files} files, {exact_count} exact-match entries, {composite_count} composite ({len(exact_db) if exact_db else 0} molecules) in {elapsed_total:.1f}s ({n_workers} workers)")
        return count

    def _parse_npmrd_spectrum(self, mol_with_h, spec_data, solvent, nucleus, db, fp_to_env, exact_db=None, inchikey=None, canon_map=None) -> tuple:
        """Parse an NPMRD spectrum record."""
        if canon_map is None:
            canon_map = {}
        count = 0
        exact_count = 0
        entries = spec_data.split('|')

        for entry in entries:
            if not entry.strip():
                continue

            parts = entry.split(';')
            if len(parts) < 3:
                continue

            try:
                shift = float(parts[0])

                atom_indices = []
                for idx_str in parts[2:]:
                    try:
                        sdf_idx = int(idx_str)

                        atom_idx = sdf_idx - 1
                        atom_indices.append(atom_idx)
                    except ValueError:
                        pass

                if not atom_indices:
                    continue

                atom_idx = atom_indices[0]

                if atom_idx >= mol_with_h.GetNumAtoms():
                    continue

                atom = mol_with_h.GetAtomWithIdx(atom_idx)
                atom_symbol = atom.GetSymbol()

                if nucleus == '13C' and atom_symbol != 'C':

                    if atom_symbol == 'H':
                        neighbors = list(atom.GetNeighbors())
                        if neighbors and neighbors[0].GetSymbol() == 'C':
                            atom_idx = neighbors[0].GetIdx()
                        else:
                            continue
                    else:
                        continue

                elif nucleus == '1H':

                    added = self._store_1h_shift(mol_with_h, atom_idx, shift, solvent,
                                                 inchikey, canon_map, exact_db, db, fp_to_env)
                    count += added
                    if inchikey and exact_db is not None:
                        exact_count += added
                    continue

                self._add_record(db, fp_to_env, mol_with_h, atom_idx, shift, solvent, nucleus)
                count += 1
                if inchikey and exact_db is not None and atom_idx in canon_map:
                    canon_idx = canon_map[atom_idx]
                    exact_db.setdefault(inchikey, {}).setdefault(solvent, {}).setdefault(nucleus, {})[canon_idx] = shift
                    exact_count += 1

            except (ValueError, IndexError) as e:

                pass

        return count, exact_count

    def _load_experimental_data(self, db, fp_to_env, exact_db=None, db_composite=None) -> int:
        """Load experimental spectrum records."""
        if not os.path.exists(EXPERIMENTAL_DATA_DIR):
            return 0

        import glob
        import time as _time
        _t0 = _time.time()

        sdf_files = glob.glob(os.path.join(EXPERIMENTAL_DATA_DIR, '*.sdf'))
        total_files = len(sdf_files)
        print(f"  [Experimental] {total_files} files found", flush=True)

        if total_files == 0:
            return 0

        n_workers = max(1, multiprocessing.cpu_count() - 1)
        print(f"  [Experimental] Starting parallel processing with {n_workers} workers...", flush=True)

        count = 0
        exact_count = 0
        composite_count = 0
        processed = 0
        batch_size = max(10, total_files // 10)

        file_args = [(filepath, i) for i, filepath in enumerate(sdf_files)]

        _t_parallel = _time.time()
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(_process_experimental_file_worker, args): args[1] for args in file_args}

            for future in as_completed(futures):
                processed += 1
                if processed % batch_size == 0:
                    elapsed = _time.time() - _t_parallel
                    rate = processed / elapsed if elapsed > 0 else 0
                    print(f"  [Experimental] {processed}/{total_files} | {count} shifts | {elapsed:.0f}s | {rate:.0f} files/s", flush=True)

                result = future.result()
                if result is None:
                    continue

                if len(result) == 4:
                    fp_records, exact_records, composite_records, rec_count = result
                else:
                    fp_records, exact_records, rec_count = result
                    composite_records = []
                count += rec_count

                for radius, fp, solvent, nucleus, shift, env_str in fp_records:
                    db[radius][fp][solvent][nucleus].append(shift)
                    if fp not in fp_to_env:
                        fp_to_env[fp] = f"r{radius}|{env_str}"

                if exact_db is not None:
                    for inchikey, solvent, nucleus, canon_idx, shift in exact_records:
                        exact_db.setdefault(inchikey, {}).setdefault(solvent, {}).setdefault(nucleus, {})[canon_idx] = shift
                        exact_count += 1

                if db_composite is not None and composite_records:
                    for radius, composite_fp, solvent, nucleus, shift in composite_records:
                        db_composite[radius][composite_fp][solvent][nucleus].append(shift)
                        composite_count += 1

        elapsed_total = _time.time() - _t0
        print(f"[FingerprintPredictor] Experimental data: {count} shifts, {exact_count} exact-match entries, {composite_count} composite in {elapsed_total:.1f}s ({n_workers} workers)")
        return count

    def _select_best_n(self, shifts: list, n: int = 10) -> list:
        """Select a compact subset of shift values."""
        import numpy as np

        if len(shifts) <= n:
            return shifts

        shifts_array = np.array(shifts)

        q1, q3 = np.percentile(shifts_array, [25, 75])
        iqr = q3 - q1

        if iqr > 0:

            lower_bound = q1 - 1.5 * iqr
            upper_bound = q3 + 1.5 * iqr
            filtered = shifts_array[(shifts_array >= lower_bound) & (shifts_array <= upper_bound)]
        else:
            filtered = shifts_array

        if len(filtered) <= n:
            return filtered.tolist()

        sorted_filtered = np.sort(filtered)
        densities = np.array([np.sum(np.abs(sorted_filtered - x) < np.std(sorted_filtered) * 0.5)
                             for x in sorted_filtered])
        center_idx = np.argmax(densities)
        center_value = sorted_filtered[center_idx]

        distances = np.abs(sorted_filtered - center_value)
        nearest_indices = np.argsort(distances)[:n]

        return sorted_filtered[nearest_indices].tolist()

    def _is_bimodal(self, shifts_list: list, nucleus: str) -> bool:
        """Check whether a shift distribution is bimodal."""
        import numpy as np

        return False

        shifts_array = np.array(shifts_list)

        try:
            from scipy import stats as scipy_stats

            skewness = scipy_stats.skew(shifts_array)

            kurtosis_excess = scipy_stats.kurtosis(shifts_array, fisher=True)

            kurtosis_pearson = kurtosis_excess + 3

            if kurtosis_pearson <= 0:
                kurtosis_pearson = 0.01

            bimodality_coeff = (skewness**2 + 1) / kurtosis_pearson

            if bimodality_coeff > 0.59:
                logger.debug(f"Bimodal detected: Sarle coefficient = {bimodality_coeff:.3f} > 0.59")
                return True
        except Exception as e:
            logger.debug(f"Sarle coefficient calculation failed (scipy may not be available): {e}")

        sorted_shifts = np.sort(shifts_array)
        overall_std = np.std(shifts_array)

        min_gap_absolute = 0.5 if nucleus == '1H' else 5.0
        min_gap_relative = 2.0 * overall_std
        min_gap = max(min_gap_absolute, min_gap_relative)

        gaps = np.diff(sorted_shifts)
        max_gap_idx = np.argmax(gaps)
        max_gap = gaps[max_gap_idx]

        if max_gap > min_gap:

            group1_size = max_gap_idx + 1
            group2_size = len(sorted_shifts) - group1_size

            if group1_size >= 3 and group2_size >= 3:
                logger.debug(f"Bimodal detected: Gap = {max_gap:.2f} ppm > {min_gap:.2f} ppm (groups: {group1_size}+{group2_size})")
                return True

        return False

    def predict_exact_proton_groups(self, mol_with_h, solvent: str) -> Dict[int, float]:
        """Return exact 1H shifts mapped by heavy-atom groups.

        Values are matched as a multiset per invariant heavy-atom fingerprint.
        This preserves two distinct shifts for a diastereotopic CH2 without
        relying on explicit-H indices or SMILES atom order.
        """
        if not FingerprintPredictor._exact_1h_groups:
            return {}

        try:
            from rdkit.Chem import inchi as rdkit_inchi
            mol_no_h = Chem.RemoveHs(mol_with_h)
            inchikey = rdkit_inchi.MolToInchiKey(mol_no_h)
            molecule_entry = FingerprintPredictor._exact_1h_groups.get(inchikey, {})
            solvent_entry = None
            for sol_key in get_solvent_variants(normalize_solvent(solvent)):
                if sol_key in molecule_entry:
                    solvent_entry = molecule_entry[sol_key]
                    break
            if not solvent_entry:
                return {}

            query_groups = defaultdict(list)
            for heavy in mol_with_h.GetAtoms():
                if heavy.GetAtomicNum() == 1:
                    continue
                h_indices = sorted(n.GetIdx() for n in heavy.GetNeighbors() if n.GetAtomicNum() == 1)
                if not h_indices:
                    continue
                heavy_key = str(get_atom_fingerprint(
                    mol_with_h, heavy.GetIdx(), RADIUS_ULTRA_EXACT, heavy_only=True
                ))
                query_groups[heavy_key].append((heavy.GetIdx(), h_indices))

            assignments = {}
            for heavy_key, groups in query_groups.items():
                records = solvent_entry.get(heavy_key, [])
                compatible = [r for r in records if int(r.get('n_h', 0)) == len(groups[0][1])]
                if len(compatible) != len(groups):
                    continue

                groups_sorted = sorted(groups, key=lambda item: item[0])
                records_sorted = sorted(
                    compatible,
                    key=lambda r: tuple(float(v) for v in r.get('shifts', [])),
                )
                for (_, h_indices), record in zip(groups_sorted, records_sorted):
                    shifts = [float(v) for v in record.get('shifts', [])]
                    if len(shifts) == 1 and len(h_indices) > 1:
                        shifts = shifts * len(h_indices)
                    if len(shifts) != len(h_indices):
                        continue
                    for h_idx, shift in zip(h_indices, sorted(shifts)):
                        assignments[h_idx] = shift
            return assignments
        except Exception as exc:
            logger.debug(f"Exact 1H group-match error: {exc}")
            return {}

    def predict_atom(self, mol_with_h, atom_idx: int, solvent: str, nucleus: str) -> Optional[Tuple[float, float, int, int]]:
        """Predict the shift for one atom."""
        if FingerprintPredictor._db is None:
            return None

        if nucleus == '1H' and is_exchangeable_proton(mol_with_h, atom_idx):
            logger.debug(f"Atom {atom_idx}: Exchangeable proton → ML")
            return None

        solvent_norm = normalize_solvent(solvent)
        import numpy as np

        # Legacy exact 1H entries are indexed by explicit-H position and are
        # unsafe until the group-preserving v14 cache is rebuilt.
        if FingerprintPredictor._exact_db and nucleus != '1H':
            try:
                from rdkit.Chem import inchi as rdkit_inchi
                mol_no_h = Chem.RemoveHs(mol_with_h)
                inchikey = rdkit_inchi.MolToInchiKey(mol_no_h)
                if inchikey and inchikey in FingerprintPredictor._exact_db:
                    exact_entry = FingerprintPredictor._exact_db[inchikey]

                    canon_smiles = Chem.MolToSmiles(mol_no_h)
                    mol_canon_no_h = Chem.MolFromSmiles(canon_smiles)
                    if mol_canon_no_h is not None:
                        match_no_h = _substruct_match_with_timeout(mol_no_h, mol_canon_no_h, timeout=1.0)
                        sdf_to_canon = {sdf_i: canon_i for canon_i, sdf_i in enumerate(match_no_h)}
                        heavy_idx_in_noh = sum(
                            1 for i in range(atom_idx)
                            if mol_with_h.GetAtomWithIdx(i).GetAtomicNum() > 1
                        )
                        canon_idx = sdf_to_canon.get(heavy_idx_in_noh)
                        if canon_idx is not None:
                            for sol_key in get_solvent_variants(solvent_norm):
                                nuc_dict = exact_entry.get(sol_key, {}).get(nucleus, {})
                                shift = nuc_dict.get(canon_idx)
                                if shift is None:
                                    shift = nuc_dict.get(str(canon_idx))
                                if shift is not None:
                                    logger.debug(f"Atom {atom_idx}: Exact-match InChIKey {inchikey}, shift={shift}")
                                    return shift, 0.0, 99, 0
            except Exception as e:
                logger.debug(f"Exact-match error: {e}")

        min_radius_required = classify_carbon_radius_requirement(mol_with_h, atom_idx)

        if min_radius_required == 99:
            logger.debug(f"Atom {atom_idx}: Requires 3D/ML (gem-dimethyl E/Z) → skip fingerprint")
            return None

        if BEYOND_AVAILABLE and FingerprintPredictor._db_composite:
            try:
                ensure_gasteiger_charges(mol_with_h)

                beyond = compute_beyond_discrete_safe(mol_with_h, atom_idx, include_chiral=False)

                for radius in [RADIUS_ULTRA_EXACT, RADIUS_HYPER_EXACT, RADIUS_VERY_EXACT, RADIUS_EXACT, RADIUS_FALLBACK]:
                    if radius < min_radius_required:
                        continue

                    db_comp = FingerprintPredictor._db_composite.get(radius, {})
                    if not db_comp:
                        continue

                    fp = get_atom_fingerprint(mol_with_h, atom_idx, radius)

                    composite_fp = (fp, beyond)

                    if composite_fp in db_comp:

                        COMPOSITE_SOLVENT_FALLBACK = {
                            'CDCl3': ['Acetone-d6', 'MeOD', 'DMSO-d6', 'Unreported'],
                            'Acetone-d6': ['CDCl3', 'MeOD', 'DMSO-d6', 'Unreported'],
                            'MeOD': ['Acetone-d6', 'CDCl3', 'DMSO-d6', 'Unreported'],
                            'DMSO-d6': ['MeOD', 'Acetone-d6', 'CDCl3', 'Unreported'],
                            'Unreported': ['CDCl3', 'DMSO-d6', 'MeOD', 'Acetone-d6'],
                            'D2O': []
                        }
                        composite_fallbacks = COMPOSITE_SOLVENT_FALLBACK.get(solvent_norm, [])

                        MIN_SAMPLES_COMPOSITE = {
                            RADIUS_ULTRA_EXACT: MIN_SAMPLES_R6,
                            RADIUS_HYPER_EXACT: MIN_SAMPLES_R5,
                            RADIUS_VERY_EXACT: MIN_SAMPLES_R4,
                            RADIUS_EXACT: MIN_SAMPLES_R3,
                            RADIUS_FALLBACK: MIN_SAMPLES_R2
                        }
                        min_samples_comp = MIN_SAMPLES_COMPOSITE.get(radius, 3)

                        solvent_variants = get_solvent_variants(solvent_norm)
                        exact_shifts = []
                        for sv in solvent_variants:
                            if sv in db_comp[composite_fp]:
                                exact_shifts.extend(db_comp[composite_fp][sv].get(nucleus, []))

                        if len(exact_shifts) >= min_samples_comp:

                            if self._is_bimodal(exact_shifts, nucleus):
                                logger.debug(f"Atom {atom_idx}: COMPOSITE R{radius} {solvent_norm} BIMODAL → fallback ML")
                                continue

                            import numpy as np
                            mean = float(np.mean(exact_shifts))
                            std = float(np.std(exact_shifts)) if len(exact_shifts) > 1 else 0.0
                            logger.debug(f"Atom {atom_idx}: COMPOSITE R{radius} {solvent_norm}, n={len(exact_shifts)}")
                            return mean, std, len(exact_shifts), radius

                        for try_solvent in composite_fallbacks:

                            fallback_variants = get_solvent_variants(try_solvent)
                            fallback_shifts = []
                            for fv in fallback_variants:
                                if fv in db_comp[composite_fp]:
                                    fallback_shifts.extend(db_comp[composite_fp][fv].get(nucleus, []))

                            if len(fallback_shifts) >= min_samples_comp:

                                if self._is_bimodal(fallback_shifts, nucleus):
                                    logger.debug(f"Atom {atom_idx}: COMPOSITE R{radius} {try_solvent} BIMODAL → essayer fallback suivant")
                                    continue

                                import numpy as np
                                mean = float(np.mean(fallback_shifts))
                                std = float(np.std(fallback_shifts)) if len(fallback_shifts) > 1 else 0.0
                                logger.debug(f"Atom {atom_idx}: COMPOSITE R{radius} {try_solvent} (fallback n={len(fallback_shifts)})")
                                return mean, std, len(fallback_shifts), radius

                        if exact_shifts:

                            if len(exact_shifts) >= 10 and self._is_bimodal(exact_shifts, nucleus):
                                logger.debug(f"Atom {atom_idx}: COMPOSITE R{radius} {solvent_norm} (low n) BIMODAL → fallback ML")
                                continue

                            import numpy as np
                            mean = float(np.mean(exact_shifts))
                            std = float(np.std(exact_shifts)) if len(exact_shifts) > 1 else 0.0
                            logger.debug(f"Atom {atom_idx}: COMPOSITE R{radius} {solvent_norm} (low n={len(exact_shifts)})")
                            return mean, std, len(exact_shifts), radius

            except Exception as e:
                logger.debug(f"Composite lookup error: {e}")

        MIN_COUNT_BY_RADIUS = {
            RADIUS_ULTRA_EXACT: MIN_SAMPLES_R6,
            RADIUS_HYPER_EXACT: MIN_SAMPLES_R5,
            RADIUS_VERY_EXACT: MIN_SAMPLES_R4,
            RADIUS_EXACT: MIN_SAMPLES_R3,
            RADIUS_FALLBACK: MIN_SAMPLES_R2
        }

        MAX_STD = MAX_STD_13C if nucleus == '13C' else MAX_STD_1H
        FALLBACK_STD = FALLBACK_MAX_STD_13C if nucleus == '13C' else FALLBACK_MAX_STD_1H
        MAX_VALUES_FOR_MEAN = 999999

        SOLVENTLESS_FALLBACK = ['Unreported', 'unknown', 'Unknown']

        SOLVENT_FALLBACK_HIERARCHY = {
            'CDCl3': ['Acetone-d6', 'MeOD', 'DMSO-d6', 'Unreported'],
            'Acetone-d6': ['CDCl3', 'MeOD', 'DMSO-d6', 'Unreported'],
            'MeOD': ['Acetone-d6', 'CDCl3', 'DMSO-d6', 'Unreported'],
            'DMSO-d6': ['MeOD', 'Acetone-d6', 'CDCl3', 'Unreported'],
            'Unreported': ['CDCl3', 'DMSO-d6', 'MeOD', 'Acetone-d6'],
            'D2O': []
        }

        fp_by_radius = {r: get_atom_fingerprint(mol_with_h, atom_idx, r)
                       for r in SUPPORTED_RADII}

        best_result = None
        best_radius = None
        best_solvent = None

        fallback_solvents = SOLVENT_FALLBACK_HIERARCHY.get(solvent_norm, [])

        radii_to_try = [RADIUS_ULTRA_EXACT, RADIUS_HYPER_EXACT, RADIUS_VERY_EXACT, RADIUS_EXACT, RADIUS_FALLBACK]

        for radius in radii_to_try:

            if radius < min_radius_required:
                continue

            db = FingerprintPredictor._db.get(radius, {})
            if not db:
                continue

            fp = fp_by_radius[radius]
            min_count = MIN_COUNT_BY_RADIUS.get(radius, 1)

            if fp not in db:
                continue

            solvent_variants = get_solvent_variants(solvent_norm)
            exact_shifts = []
            for sv in solvent_variants:
                if sv in db[fp]:
                    exact_shifts.extend(db[fp][sv].get(nucleus, []))

            if len(exact_shifts) >= min_count:
                mean, std, count = self._compute_stats_no_check(exact_shifts, MAX_VALUES_FOR_MEAN, nucleus)

                if self._is_bimodal(exact_shifts, nucleus):
                    logger.debug(f"Atom {atom_idx}: Morgan R{radius} {solvent_norm} BIMODAL → skip")

                else:

                    if std <= MAX_STD:
                        logger.debug(f"Atom {atom_idx}: {solvent_norm} R{radius}, n={count}, std={std:.2f}")
                        return mean, std, count, radius

                    if std <= FALLBACK_STD and (best_result is None or std < best_result[1]):
                        best_result = (mean, std, count, exact_shifts)
                        best_radius = radius
                        best_solvent = solvent_norm

            if len(exact_shifts) < min_count or (len(exact_shifts) >= min_count and self._is_bimodal(exact_shifts, nucleus)):

                for solvent_candidate in fallback_solvents:
                    candidate_variants = get_solvent_variants(solvent_candidate)
                    candidate_shifts = []
                    for cv in candidate_variants:
                        if cv in db[fp]:
                            candidate_shifts.extend(db[fp][cv].get(nucleus, []))

                    if len(candidate_shifts) >= min_count:
                        mean, std, count = self._compute_stats_no_check(candidate_shifts, MAX_VALUES_FOR_MEAN, nucleus)

                        if self._is_bimodal(candidate_shifts, nucleus):
                            logger.debug(f"Atom {atom_idx}: Morgan R{radius} {solvent_candidate} BIMODAL → essayer fallback suivant")
                            continue

                        if std <= MAX_STD:
                            logger.debug(f"Atom {atom_idx}: {solvent_candidate} R{radius} (fallback), n={count}, std={std:.2f}")
                            return mean, std, count, radius

                        if std <= FALLBACK_STD and (best_result is None or std < best_result[1]):
                            best_result = (mean, std, count, candidate_shifts)
                            best_radius = radius
                            best_solvent = solvent_candidate

        if best_result is not None:
            mean, std, count, shifts_list = best_result

            if len(shifts_list) >= 10 and self._is_bimodal(shifts_list, nucleus):
                logger.debug(f"Atom {atom_idx}: Best result R{best_radius} from {best_solvent} is BIMODAL → ML")
                return None

            logger.debug(f"Atom {atom_idx}: Using best result R{best_radius} from {best_solvent}")
            return (mean, std, count, best_radius)

        logger.debug(f"Atom {atom_idx}: No match in {solvent_norm} (R6/R5/R4/R3/R2) → ML")
        return None

    def _search_any_solvent(self, mol_with_h, atom_idx: int, nucleus: str,
                          min_count_by_radius: dict, max_values: int, exclude_solventless: bool = False) -> Optional[Tuple[float, float, int, int]]:
        """Search matching records across solvents."""
        import numpy as np

        SOLVENTLESS_FALLBACK = ['Unreported', 'unknown', 'Unknown']

        for radius in SUPPORTED_RADII:
            db = FingerprintPredictor._db.get(radius, {})
            if not db:
                continue

            fp = get_atom_fingerprint(mol_with_h, atom_idx, radius)

            if fp in db:

                all_shifts = []
                for solvent, nuc_dict in db[fp].items():

                    if exclude_solventless and solvent in SOLVENTLESS_FALLBACK:
                        continue

                    shifts = nuc_dict.get(nucleus, [])
                    if shifts:
                        all_shifts.extend(shifts)

                if len(all_shifts) >= min_count_by_radius.get(radius, 1):
                    mean, std, count = self._compute_stats(all_shifts, max_values, nucleus)
                    logger.debug(f"Atom {atom_idx}: Any solvent R{radius}, n={count}, std={std:.2f}")
                    return mean, std, count, radius

        return None

    def _select_best_solventless(self, shifts: list, max_count: int) -> list:
        """Select representative solventless values."""
        import numpy as np

        if len(shifts) <= max_count:
            return shifts

        shifts_array = np.array(shifts)

        q1, q3 = np.percentile(shifts_array, [25, 75])
        iqr = q3 - q1

        if iqr > 0:
            lower_bound = q1 - 2.0 * iqr
            upper_bound = q3 + 2.0 * iqr
            filtered = shifts_array[(shifts_array >= lower_bound) & (shifts_array <= upper_bound)]
        else:
            filtered = shifts_array

        if len(filtered) <= max_count:
            return filtered.tolist()

        n_bins = min(20, max(5, len(filtered) // 10))
        hist, bin_edges = np.histogram(filtered, bins=n_bins)

        densest_bin_idx = np.argmax(hist)
        densest_bin_center = (bin_edges[densest_bin_idx] + bin_edges[densest_bin_idx + 1]) / 2

        distances = np.abs(filtered - densest_bin_center)
        nearest_indices = np.argsort(distances)[:max_count]

        return filtered[nearest_indices].tolist()

    def _remove_outliers_iqr(self, shifts_list: list, iqr_factor: float = 1.5) -> list:
        """Remove IQR-based outliers."""
        import numpy as np

        if len(shifts_list) < 4:

            return shifts_list

        shifts_array = np.array(shifts_list)

        q1 = np.percentile(shifts_array, 25)
        q3 = np.percentile(shifts_array, 75)
        iqr = q3 - q1

        lower_bound = q1 - iqr_factor * iqr
        upper_bound = q3 + iqr_factor * iqr

        filtered = shifts_array[(shifts_array >= lower_bound) & (shifts_array <= upper_bound)]

        return filtered.tolist()

    def _filter_coherent_shifts(self, shifts_list: list, nucleus: str) -> list:
        """Filter shifts by statistical coherence."""
        import numpy as np

        max_std = MAX_STD_13C if nucleus == '13C' else MAX_STD_1H

        if len(shifts_list) <= 3:
            return shifts_list

        shifts_array = np.array(shifts_list)

        current_std = np.std(shifts_array)
        if current_std <= max_std:
            return shifts_list

        best_subset = []
        best_size = 0

        percentile_ranges = [
            (0, 100),
            (5, 95),
            (10, 90),
            (15, 85),
            (20, 80),
            (25, 75),
            (30, 70),
            (35, 65),
            (40, 60)
        ]

        for pct_low, pct_high in percentile_ranges:
            low_val = np.percentile(shifts_array, pct_low)
            high_val = np.percentile(shifts_array, pct_high)

            mask = (shifts_array >= low_val) & (shifts_array <= high_val)
            subset = shifts_array[mask]

            if len(subset) < 2:
                continue

            subset_std = np.std(subset)

            if subset_std <= max_std and len(subset) > best_size:
                best_subset = subset.tolist()
                best_size = len(subset)
                break

        if not best_subset:
            median = np.median(shifts_array)
            distances = np.abs(shifts_array - median)

            n_keep = max(3, len(shifts_array) // 4)
            closest_indices = np.argsort(distances)[:n_keep]
            best_subset = shifts_array[closest_indices].tolist()

        return best_subset

    def _compute_stats(self, shifts_list: list, max_values: int, nucleus: str = '13C') -> Tuple[float, float, int]:
        """Compute filtered shift statistics."""
        import numpy as np

        original_count = len(shifts_list)

        cleaned_shifts = self._remove_outliers_iqr(shifts_list)

        coherent_shifts = self._filter_coherent_shifts(cleaned_shifts, nucleus)

        if len(coherent_shifts) > max_values:
            selected = self._select_best_n(coherent_shifts, max_values)
        else:
            selected = coherent_shifts

        return float(np.mean(selected)), float(np.std(selected)), original_count

    def _compute_stats_no_check(self, shifts_list: list, max_values: int, nucleus: str = '13C') -> Tuple[float, float, int]:
        """Compute shift statistics without threshold validation."""
        import numpy as np

        original_count = len(shifts_list)

        cleaned_shifts = self._remove_outliers_iqr(shifts_list)

        coherent_shifts = self._filter_coherent_shifts(cleaned_shifts, nucleus)

        if len(coherent_shifts) > max_values:
            selected = self._select_best_n(coherent_shifts, max_values)
        else:
            selected = coherent_shifts

        return float(np.mean(selected)), float(np.std(selected)), original_count

    def _compute_stats_with_original_count(self, selected_shifts: list, original_count: int) -> Tuple[float, float, int]:
        """Compute statistics while preserving the original count."""
        import numpy as np
        return float(np.mean(selected_shifts)), float(np.std(selected_shifts)), original_count

    def _select_best_n(self, shifts_list: list, n: int) -> list:
        """Select a compact subset of shift values."""
        import numpy as np

        if len(shifts_list) <= n:
            return shifts_list

        sorted_shifts = sorted(shifts_list)

        best_std = float('inf')
        best_start = 0

        for i in range(len(sorted_shifts) - n + 1):
            window = sorted_shifts[i:i+n]
            window_std = np.std(window)
            if window_std < best_std:
                best_std = window_std
                best_start = i

        return sorted_shifts[best_start:best_start+n]

    def get_environment(self, mol_with_h, atom_idx: int) -> str:
        """Return the atom-environment description."""
        return get_atom_env_string(mol_with_h, atom_idx)

    def is_available(self) -> bool:
        """Return whether the database is loaded."""
        return FingerprintPredictor._db is not None and len(FingerprintPredictor._db) > 0

def get_predictor() -> FingerprintPredictor:
    """Return the singleton predictor instance."""
    return FingerprintPredictor()

if __name__ == '__main__':
    import sys

    smiles = "CCO"
    solvent = "CDCl3"

    if len(sys.argv) >= 2:
        smiles = sys.argv[1]
    if len(sys.argv) >= 3:
        solvent = sys.argv[2]

    print(f"Testing FingerprintPredictor with: {smiles} in {solvent}")
    print("=" * 60)

    predictor = FingerprintPredictor()

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Invalid SMILES: {smiles}")
        sys.exit(1)

    mol_with_h = Chem.AddHs(mol)

    print(f"\nPredictions (multi-radius hierarchy):")
    print("-" * 60)

    found_count = 0
    not_found_count = 0

    for atom_idx in range(mol_with_h.GetNumAtoms()):
        atom = mol_with_h.GetAtomWithIdx(atom_idx)
        symbol = atom.GetSymbol()

        if symbol == 'C':
            nucleus = '13C'
        elif symbol == 'H':
            nucleus = '1H'
        else:
            continue

        env = predictor.get_environment(mol_with_h, atom_idx)
        result = predictor.predict_atom(mol_with_h, atom_idx, solvent, nucleus)

        if result:
            shift, std, count, radius_used = result
            print(f"  {atom_idx:3d} | {nucleus:3s} | {env:30s} | {shift:7.2f} ± {std:5.2f} (n={count}, R{radius_used})")
            found_count += 1
        else:
            print(f"  {atom_idx:3d} | {nucleus:3s} | {env:30s} | NOT FOUND (use ML fallback)")
            not_found_count += 1

    print("-" * 60)
    print(f"Found: {found_count}, Not found (need ML): {not_found_count}")
