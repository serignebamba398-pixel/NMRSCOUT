# backend/scripts/ml_sdf_utils.py
# -*- coding: utf-8 -*-
import os, re, json
from typing import Dict, List, Tuple, Optional
from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect

RDLogger.DisableLog('rdApp.*')

# --------- molecule loading: keep H + indexing ----------
def read_molecule_keepH(sdf_path: str) -> Optional[Chem.Mol]:
    try:
        suppl = Chem.SDMolSupplier(sdf_path, removeHs=False, sanitize=False)
        if not suppl or suppl[0] is None:
            return None
        mol = suppl[0]
        try: mol.UpdatePropertyCache(strict=False)
        except: pass
        # Sanitize "best effort" (if it fails, keep the raw mol)
        try: Chem.SanitizeMol(mol)
        except: pass
        try: Chem.AssignStereochemistry(mol, force=True, cleanIt=True, flagPossibleStereoCenters=True)
        except: pass
        return mol
    except Exception:
        return None

# --------- parse all SDF properties -----------
def _strip_bom(s: str) -> str:
    return s.replace('\ufeff', '')

def _norm_key(k: str) -> str:
    k = _strip_bom(k or '').strip()
    # normalize: lowercase, remove diacritics, spaces -> _
    k = k.encode('ascii', 'ignore').decode('ascii')
    k = re.sub(r'[\s\-./\\|;:,]+', '_', k.lower())
    k = re.sub(r'[^a-z0-9_]+', '', k)
    k = k.strip('_')
    return k

def scan_sdf_properties(sdf_path: str) -> Dict[str,str]:
    """
    Returns {normalized_key: multiline_text} for all properties.
    Supports: '>  <NAME>' ... (following lines until next property or $$$$)
    """
    props: Dict[str, str] = {}
    cur_key = None
    buf = []
    with open(sdf_path, 'r', encoding='utf-8', errors='ignore') as f:
        for raw in f:
            line = raw.rstrip('\n')
            if line.startswith('>'):
                # flush previous
                if cur_key is not None:
                    props[_norm_key(cur_key)] = '\n'.join(buf).strip()
                    buf = []
                # extract name between <...>
                m = re.search(r'<([^>]+)>', line)
                cur_key = m.group(1) if m else '(unnamed)'
            elif line.strip() == '$$$$':
                if cur_key is not None:
                    props[_norm_key(cur_key)] = '\n'.join(buf).strip()
                    cur_key, buf = None, []
                break
            else:
                if cur_key is not None:
                    buf.append(_strip_bom(line))
    if cur_key is not None and buf:
        props[_norm_key(cur_key)] = '\n'.join(buf).strip()
    return props

def _parse_float_vector(s: str) -> List[float]:
    out: List[float] = []
    if not s: return out
    # replace commas with spaces, then split
    for tok in re.split(r'[,\s]+', s.strip()):
        if not tok: continue
        try:
            out.append(float(tok))
        except:
            # ignore non-numeric tokens
            pass
    return out

# --------- select the RIGHT vector (per-atom) ----------
def get_atom_shift_vector(sdf_path: str) -> List[float]:
    """
        Heuristic:
            - read all properties; try common keys first: (1), pred_shifts...
            - pick a vector whose length == number of atoms (ideal); otherwise the longest one with length ≥ n_atoms.
    """
    mol = read_molecule_keepH(sdf_path)
    if mol is None:
        return []
    n_atoms = mol.GetNumAtoms()

    props = scan_sdf_properties(sdf_path)

    # candidates ordered by name priority
    priority_keys = [
        '1', '_1', '(1)', 'ml_shifts', 'pred_shifts', 'predicted_shifts',
        'shifts', 'shifts_ppm', 'ml_shift_vector', 'shift_vector', 'ppm'
    ]
    # map potential normalized names to present ones
    def find_in_props(key: str) -> Optional[str]:
        k = _norm_key(key)
        if k in props: return k
        # also search contains
        for kk in props.keys():
            if k == kk or k in kk:
                return kk
        return None

    # 1) try priority names
    for k in priority_keys:
        kk = find_in_props(k)
        if kk:
            vec = _parse_float_vector(props[kk])
            if len(vec) >= n_atoms:
                return vec

    # 2) otherwise try all: choose the one that fits best
    best = []
    for kk, txt in props.items():
        vec = _parse_float_vector(txt)
        if len(vec) == n_atoms:
            return vec
        if len(vec) > len(best):
            best = vec

    return best

# --------- ECFP4 (bitstring) ----------
def ecfp4_bitstring(mol: Chem.Mol, radius: int = 2, nBits: int = 1024) -> str:
    try:
        fp = GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
        return fp.ToBitString()
    except Exception:
        try:
            fp = GetMorganFingerprintAsBitVect(Chem.RemoveHs(mol), radius, nBits=nBits)
            return fp.ToBitString()
        except Exception:
            return "0" * nBits
