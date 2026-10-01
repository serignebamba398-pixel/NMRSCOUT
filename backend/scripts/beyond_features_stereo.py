#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""Stereo-aware molecular feature utilities."""


from rdkit import Chem

from rdkit.Chem import AllChem

from rdkit.Chem.rdMolDescriptors import CalcNumAtomStereoCenters

import hashlib


# Gasteiger charge helper

def ensure_gasteiger_charges(mol):

    """Compute Gasteiger charges once and report whether they are available."""


    try:

        # Reuse existing charges when available

        mol.GetAtomWithIdx(0).GetProp('_GasteigerCharge')

        return True

    except:

        pass


    try:

        AllChem.ComputeGasteigerCharges(mol)

        return True

    except:

        return False


# Internal 2D feature helpers


BEYOND_2D_AVAILABLE = True


def _get_gasteiger_charge(atom):

    """Return the Gasteiger charge for an atom, or 0.0 when unavailable."""

    try:

        import numpy as np

        charge = float(atom.GetProp('_GasteigerCharge'))

        if np.isnan(charge) or np.isinf(charge):

            return 0.0

        return charge

    except:

        return 0.0


def _detect_nearest_carbonyl_type(mol, atom_idx):

    """Return the nearest carbonyl category within three bonds."""


    visited = set()

    queue = [(atom_idx, 0)]

    best_type = 0

    best_dist = 99


    while queue:

        current, dist = queue.pop(0)

        if dist > 3 or current in visited:

            continue

        visited.add(current)

        c_atom = mol.GetAtomWithIdx(current)


        if c_atom.GetSymbol() == 'C':

            has_double_O = False

            has_single_O_ester = False

            has_single_OH = False

            has_single_N = False


            for nb in c_atom.GetNeighbors():

                bond = mol.GetBondBetweenAtoms(current, nb.GetIdx())

                if bond is None:

                    continue


                bond_order = bond.GetBondTypeAsDouble()


                if nb.GetSymbol() == 'O' and bond_order >= 2.0:

                    has_double_O = True

                elif nb.GetSymbol() == 'O' and bond_order < 2.0:

                    # Distinguish OH from O-R

                    h_count = sum(1 for n in nb.GetNeighbors() if n.GetSymbol() == 'H')

                    if h_count > 0:

                        has_single_OH = True

                    else:

                        has_single_O_ester = True

                elif nb.GetSymbol() == 'N' and bond_order < 2.0:

                    has_single_N = True


            if has_double_O and dist < best_dist:

                best_dist = dist

                if has_single_OH:

                    best_type = 4  # Carboxylic acid

                elif has_single_O_ester:

                    best_type = 2  # Ester

                elif has_single_N:

                    best_type = 3  # Amide

                else:

                    best_type = 1  # Ketone/aldehyde


        # Continue BFS through heavy atoms

        if dist < 3:

            for nb in c_atom.GetNeighbors():

                if nb.GetSymbol() != 'H':

                    queue.append((nb.GetIdx(), dist + 1))


    return best_type


def _detect_heterocycle_type(mol, atom_idx):

    """Return the aromatic heterocycle ring-size category for an atom."""


    atom = mol.GetAtomWithIdx(atom_idx)


    if not atom.GetIsAromatic():

        return 0


    ring_info = mol.GetRingInfo()

    aromatic_ring_sizes = []


    for ring in ring_info.AtomRings():

        if atom_idx in ring:

            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]

            # Require an aromatic ring containing a heteroatom

            if all(a.GetIsAromatic() for a in ring_atoms):

                has_hetero = any(a.GetSymbol() in ('N', 'O', 'S') for a in ring_atoms)

                if has_hetero:

                    aromatic_ring_sizes.append(len(ring))


    if not aromatic_ring_sizes:

        return 0


    # Return the smallest aromatic heterocycle size

    return min(aromatic_ring_sizes)


def _detect_donor_acceptor_environment(mol, atom_idx):

    """Return the local donor/acceptor environment category."""


    visited = set()

    queue = [(atom_idx, 0)]

    donor_score = 0

    acceptor_score = 0


    while queue:

        current, dist = queue.pop(0)

        if dist > 3 or current in visited:

            continue

        visited.add(current)

        c_atom = mol.GetAtomWithIdx(current)


        if current != atom_idx:

            symbol = c_atom.GetSymbol()


            # Detect donor and acceptor groups

            if symbol in ('N', 'O') and not c_atom.GetIsAromatic():

                # Check whether the atom belongs to an acceptor group

                is_in_acceptor = False

                for nb in c_atom.GetNeighbors():

                    bond = mol.GetBondBetweenAtoms(current, nb.GetIdx())

                    if bond and bond.GetBondTypeAsDouble() >= 2.0 and nb.GetSymbol() in ('C', 'N', 'S'):

                        is_in_acceptor = True

                        break


                if is_in_acceptor:

                    acceptor_score += 2

                elif c_atom.GetTotalNumHs() > 0:

                    donor_score += 2  # Strong donor

                else:

                    donor_score += 1  # Weak donor


            # Weak acceptor

            if symbol in ('F', 'Cl', 'Br', 'I'):

                acceptor_score += 1


        # Continue BFS

        if dist < 3:

            for nb in c_atom.GetNeighbors():

                if nb.GetSymbol() != 'H':

                    queue.append((nb.GetIdx(), dist + 1))


    # Classify the environment

    if acceptor_score > donor_score:

        return 4 if acceptor_score >= 3 else 3

    elif donor_score > acceptor_score:

        return 1 if donor_score >= 2 else 2

    return 0


def _detect_alpha_to_heteroatom(mol, atom_idx):

    """Return the category of an adjacent heteroatom, if present."""


    atom = mol.GetAtomWithIdx(atom_idx)


    # Use the parent atom for hydrogen inputs

    if atom.GetSymbol() == 'H':

        parents = [n for n in atom.GetNeighbors()]

        if not parents:

            return 0

        atom = parents[0]


    for neighbor in atom.GetNeighbors():

        symbol = neighbor.GetSymbol()

        if symbol == 'O':

            return 1

        elif symbol == 'N':

            return 2

        elif symbol == 'S':

            return 3

        elif symbol in ('F', 'Cl', 'Br', 'I'):

            return 4


    return 0


def _detect_aromatic_system_type(mol, atom_idx):

    """Return the aromatic-system category associated with an atom."""


    atom = mol.GetAtomWithIdx(atom_idx)


    # Use the parent atom for hydrogen inputs

    if atom.GetSymbol() == 'H':

        parents = [n for n in atom.GetNeighbors()]

        if not parents:

            return 0

        atom = parents[0]


    # Find an aromatic neighbor

    aromatic_neighbor = None

    for nb in atom.GetNeighbors():

        if nb.GetIsAromatic():

            aromatic_neighbor = nb

            break


    if aromatic_neighbor is None:

        # Fall back to the atom itself when aromatic

        if not atom.GetIsAromatic():

            return 0

        aromatic_neighbor = atom


    # Inspect the aromatic ring system

    ring_info = mol.GetRingInfo()


    for ring in ring_info.AtomRings():

        if aromatic_neighbor.GetIdx() not in ring:

            continue


        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]

        if not all(a.GetIsAromatic() for a in ring_atoms):

            continue


        ring_size = len(ring)

        n_count = sum(1 for a in ring_atoms if a.GetSymbol() == 'N')

        o_count = sum(1 for a in ring_atoms if a.GetSymbol() == 'O')

        s_count = sum(1 for a in ring_atoms if a.GetSymbol() == 'S')


        # Check whether the ring is fused

        is_fused = False

        for other_ring in ring_info.AtomRings():

            if other_ring != ring:

                shared = set(ring) & set(other_ring)

                if len(shared) >= 2:

                    is_fused = True

                    break


        # Classify the ring

        if ring_size == 5:

            if is_fused and n_count >= 1:

                # Check for an NH-containing fused ring

                for a in ring_atoms:

                    if a.GetSymbol() == 'N' and a.GetTotalNumHs() > 0:

                        return 5  # Indole

            if n_count >= 2:

                return 6  # Imidazole

            if n_count == 1:

                return 4  # Pyrrole

            if o_count >= 1 or s_count >= 1:

                return 7  # Furan/thiophene


        elif ring_size == 6:

            if n_count >= 2:

                return 3  # Pyrimidine

            if n_count == 1:

                return 2  # Pyridine

            if n_count == 0 and o_count == 0 and s_count == 0:

                return 1  # Benzene


    return 0


def _detect_allylic_position(mol, atom_idx):

    """Return the allylic-environment category for an atom."""


    atom = mol.GetAtomWithIdx(atom_idx)


    # Use the parent atom for hydrogen inputs

    if atom.GetSymbol() == 'H':

        parents = [n for n in atom.GetNeighbors()]

        if not parents:

            return 0

        atom = parents[0]


    # Require an sp3 carbon

    if atom.GetSymbol() != 'C':

        return 0

    if str(atom.GetHybridization()) != 'SP3':

        return 0


    # Find non-aromatic sp2 carbon neighbors

    sp2_neighbors = []

    for nb in atom.GetNeighbors():

        if nb.GetSymbol() == 'C' and str(nb.GetHybridization()) == 'SP2':


            if not nb.GetIsAromatic():

                sp2_neighbors.append(nb)


    if len(sp2_neighbors) == 0:

        return 0


    if len(sp2_neighbors) >= 2:

        return 3  # Bis-allylic


    # Check for a prenyl environment

    sp2_nb = sp2_neighbors[0]

    methyl_count = 0

    for nb2 in sp2_nb.GetNeighbors():

        if nb2.GetIdx() != atom.GetIdx():

            if nb2.GetSymbol() == 'C' and nb2.GetTotalNumHs() == 3:

                methyl_count += 1


    if methyl_count >= 2:

        return 4  # Prenyl


    # Check for conjugation

    for nb2 in sp2_nb.GetNeighbors():

        if nb2.GetIdx() != atom.GetIdx():

            if nb2.GetSymbol() == 'C' and str(nb2.GetHybridization()) == 'SP2':

                if not nb2.GetIsAromatic():

                    return 2  # Conjugated


    return 1  # Simple allylic


def _detect_vinylic_conjugation(mol, atom_idx):

    """Return the vinylic-conjugation category for an atom."""


    atom = mol.GetAtomWithIdx(atom_idx)


    # Use the parent atom for hydrogen inputs

    if atom.GetSymbol() == 'H':

        parents = [n for n in atom.GetNeighbors()]

        if not parents:

            return 0

        parent = parents[0]

    else:

        parent = atom


    # Require a carbon parent

    if parent.GetSymbol() != 'C':

        return 0


    # Search within four bonds


    visited = set()

    queue = [(parent.GetIdx(), 0)]


    has_vinyl = False

    has_aromatic_conj = False

    is_exocyclic = False


    while queue:

        current_idx, dist = queue.pop(0)

        if dist > 4 or current_idx in visited:

            continue

        visited.add(current_idx)


        current = mol.GetAtomWithIdx(current_idx)


        # Find C=C bonds

        for nb in current.GetNeighbors():

            bond = mol.GetBondBetweenAtoms(current_idx, nb.GetIdx())

            if bond and bond.GetBondTypeAsDouble() >= 2.0:

                if nb.GetSymbol() == 'C' and not nb.GetIsAromatic():

                    has_vinyl = True


                    # Check whether the double bond is exocyclic

                    ring_info = mol.GetRingInfo()

                    current_in_ring = any(current_idx in r for r in ring_info.AtomRings())

                    nb_in_ring = any(nb.GetIdx() in r for r in ring_info.AtomRings())

                    if current_in_ring != nb_in_ring:

                        is_exocyclic = True


                    # Check aromatic conjugation

                    for nb2 in nb.GetNeighbors():

                        if nb2.GetIsAromatic():

                            has_aromatic_conj = True


        if dist < 4:

            for nb in current.GetNeighbors():

                if nb.GetSymbol() != 'H':

                    queue.append((nb.GetIdx(), dist + 1))


    if not has_vinyl:

        return 0


    if is_exocyclic and has_aromatic_conj:

        return 3  # Exocyclic vinylidene

    if has_aromatic_conj:

        return 4 if parent.GetTotalNumHs() >= 1 else 2  # Styrene-like or vinyl-aromatic


    return 1  # Simple vinylic


# Local stereochemistry


def get_stereo_environment(mol, atom_idx, radius=3):

    """Return local chiral and double-bond stereochemistry signatures."""


    try:

        # Assign stereochemistry

        Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)

    except:

        pass


    # Build the local environment

    env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom_idx)


    atoms_in_env = set([atom_idx])

    bonds_in_env = set()


    for bond_idx in env:

        bond = mol.GetBondWithIdx(bond_idx)

        atoms_in_env.add(bond.GetBeginAtomIdx())

        atoms_in_env.add(bond.GetEndAtomIdx())

        bonds_in_env.add(bond_idx)


    # Collect chiral centers

    chiral_list = []

    for idx in sorted(atoms_in_env):

        atom = mol.GetAtomWithIdx(idx)


        # Use ChiralTag when available

        chiral_tag = atom.GetChiralTag()

        if chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:

            chiral_list.append('R')

        elif chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:

            chiral_list.append('S')

        elif atom.HasProp('_CIPCode'):

            # Fall back to the CIP code

            cip = atom.GetProp('_CIPCode')

            if cip in ['R', 'S']:

                chiral_list.append(cip)


    # Collect E/Z double-bond assignments

    ez_list = []

    for bond_idx in sorted(bonds_in_env):

        bond = mol.GetBondWithIdx(bond_idx)


        if bond.GetBondType() != Chem.BondType.DOUBLE:

            continue


        # Read the stereo assignment

        stereo = bond.GetStereo()


        if stereo == Chem.BondStereo.STEREOE:

            ez_list.append('E')

        elif stereo == Chem.BondStereo.STEREOZ:

            ez_list.append('Z')

        elif stereo == Chem.BondStereo.STEREOTRANS:

            ez_list.append('E')

        elif stereo == Chem.BondStereo.STEREOCIS:

            ez_list.append('Z')

        else:


            begin_atom = bond.GetBeginAtom()

            end_atom = bond.GetEndAtom()


            begin_neighbors = [n for n in begin_atom.GetNeighbors() if n.GetIdx() != end_atom.GetIdx()]

            end_neighbors = [n for n in end_atom.GetNeighbors() if n.GetIdx() != begin_atom.GetIdx()]


            if len(begin_neighbors) >= 1 and len(end_neighbors) >= 1:


                pass


    # Build compact signatures

    chiral_sig = ''.join(chiral_list) if chiral_list else ''

    ez_sig = ''.join(ez_list) if ez_list else ''


    return (chiral_sig, ez_sig, len(chiral_list), len(ez_list))


def encode_relative_chirality(chiral_sig):

    """Encode a relative chirality signature as a deterministic integer category."""


    if not chiral_sig:

        return 0


    chiral_list = list(chiral_sig)


    if len(chiral_list) == 1:

        return 1  # A single center maps to one category


    # Encode pairwise relative relationships

    relations = []

    for i in range(len(chiral_list) - 1):

        if chiral_list[i] == chiral_list[i + 1]:

            relations.append("like")

        else:

            relations.append("unlike")


    rel_string = "-".join(relations)

    rel_hash = int(hashlib.md5(rel_string.encode()).hexdigest(), 16)

    return (rel_hash % 997) + 2  # Range: 2-999


def _distance_to_nearest_heteroatom(mol, atom_idx, max_dist=8):

    """Return the BFS distance to the nearest heteroatom."""


    from collections import deque


    hetero_symbols = {'O', 'N', 'S', 'F', 'Cl', 'Br', 'I'}


    atom = mol.GetAtomWithIdx(atom_idx)

    if atom.GetSymbol() in hetero_symbols:

        return 0


    # Breadth-first search

    visited = {atom_idx}

    queue = deque([(atom_idx, 0)])


    while queue:

        current_idx, dist = queue.popleft()


        if dist >= max_dist:

            break


        current_atom = mol.GetAtomWithIdx(current_idx)

        for neighbor in current_atom.GetNeighbors():

            neighbor_idx = neighbor.GetIdx()

            if neighbor_idx in visited:

                continue


            visited.add(neighbor_idx)


            # Check for a heteroatom

            if neighbor.GetSymbol() in hetero_symbols:

                return dist + 1


            queue.append((neighbor_idx, dist + 1))


    return max_dist


def _distance_to_nearest_aromatic(mol, atom_idx, max_dist=8):

    """Return the BFS distance to the nearest aromatic atom."""


    from collections import deque


    atom = mol.GetAtomWithIdx(atom_idx)

    if atom.GetIsAromatic():

        return 0


    # Breadth-first search

    visited = {atom_idx}

    queue = deque([(atom_idx, 0)])


    while queue:

        current_idx, dist = queue.popleft()


        if dist >= max_dist:

            break


        current_atom = mol.GetAtomWithIdx(current_idx)

        for neighbor in current_atom.GetNeighbors():

            neighbor_idx = neighbor.GetIdx()

            if neighbor_idx in visited:

                continue


            visited.add(neighbor_idx)


            # Check for an aromatic atom

            if neighbor.GetIsAromatic():

                return dist + 1


            queue.append((neighbor_idx, dist + 1))


    return max_dist


def encode_stereo_environment(chiral_sig, ez_sig):

    """Encode stereochemistry signatures as deterministic numeric codes."""


    if not chiral_sig:

        chiral_code = 0

    else:

        # Use deterministic MD5 hashing

        chiral_hash = int(hashlib.md5(chiral_sig.encode()).hexdigest(), 16)

        chiral_code = (chiral_hash % 999) + 1


    if not ez_sig:

        ez_code = 0

    else:

        # Use deterministic MD5 hashing

        ez_hash = int(hashlib.md5(ez_sig.encode()).hexdigest(), 16)

        ez_code = (ez_hash % 999) + 1


    return (chiral_code, ez_code)


# Extended feature computation


def compute_beyond_with_stereo(mol, atom_idx, radius_stereo=6, include_chiral=True):

    """Compute the extended feature tuple for the selected atom."""


    atom = mol.GetAtomWithIdx(atom_idx)


    # Use the parent atom for hydrogen inputs

    if atom.GetSymbol() == 'H':

        parents = list(atom.GetNeighbors())

        if parents:

            parent_idx = parents[0].GetIdx()

        else:

            return (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

    else:

        parent_idx = atom_idx


    # Base 2D features


    if BEYOND_2D_AVAILABLE:

        # Compute base features

        ensure_gasteiger_charges(mol)


        carbonyl_type = _detect_nearest_carbonyl_type(mol, parent_idx)

        hetero_ring = _detect_heterocycle_type(mol, parent_idx)

        charge = _get_gasteiger_charge(atom)

        charge_bin = max(-2, min(2, int(round(charge * 5))))

        donor_acceptor = _detect_donor_acceptor_environment(mol, parent_idx)

        alpha_hetero = _detect_alpha_to_heteroatom(mol, atom_idx)

        aromatic_type = _detect_aromatic_system_type(mol, atom_idx)

        allylic_pos = _detect_allylic_position(mol, atom_idx)

        vinyl_conj = _detect_vinylic_conjugation(mol, atom_idx)

    else:

        # Fallback values

        carbonyl_type = 0

        hetero_ring = 0

        charge_bin = 0

        donor_acceptor = 0

        alpha_hetero = 0

        aromatic_type = 0

        allylic_pos = 0

        vinyl_conj = 0


    # Stereochemistry features


    try:

        chiral_sig, ez_sig, num_chiral, num_ez = get_stereo_environment(mol, atom_idx, radius_stereo)


        # Encode relative configuration


        relative_config_code = encode_relative_chirality(chiral_sig)


        # Disable chirality when requested

        if not include_chiral:

            relative_config_code = 0


        # Encode E/Z assignments deterministically

        if not ez_sig:

            ez_code = 0

        else:

            ez_hash = int(hashlib.md5(ez_sig.encode()).hexdigest(), 16)

            ez_code = (ez_hash % 999) + 1


    except Exception:

        relative_config_code = 0

        ez_code = 0


    # Atom-state features


    # Feature 10: attached hydrogen count

    atom_for_features = mol.GetAtomWithIdx(parent_idx)

    nH = sum(1 for n in atom_for_features.GetNeighbors() if n.GetSymbol() == 'H')


    # Feature 11: hybridization code

    hyb_str = str(atom_for_features.GetHybridization())

    if 'SP3' in hyb_str:

        hybridization_code = 3

    elif 'SP2' in hyb_str:

        hybridization_code = 2

    elif 'SP' in hyb_str:

        hybridization_code = 1

    else:

        hybridization_code = 0


    # Distance features


    # Feature 12: distance to the nearest heteroatom

    dist_heteroatom = _distance_to_nearest_heteroatom(mol, parent_idx, max_dist=8)


    # Feature 13: distance to the nearest aromatic atom

    dist_aromatic = _distance_to_nearest_aromatic(mol, parent_idx, max_dist=8)


    # Return the 14-element feature tuple


    return (

        carbonyl_type,          # 0

        hetero_ring,            # 1

        charge_bin,             # 2

        donor_acceptor,         # 3

        alpha_hetero,           # 4

        aromatic_type,          # 5

        allylic_pos,            # 6

        vinyl_conj,             # 7

        relative_config_code,   # 8

        ez_code,                # 9

        nH,                     # 10

        hybridization_code,     # 11

        dist_heteroatom,        # 12

        dist_aromatic           # 13

    )


def compute_beyond_with_stereo_safe(mol, atom_idx, radius_stereo=6, include_chiral=True):

    """Return the extended feature tuple, or zeros on failure."""


    try:

        return compute_beyond_with_stereo(mol, atom_idx, radius_stereo, include_chiral)

    except Exception:

        return (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)


# Utilities


def beyond_stereo_to_string(beyond_tuple):

    """Convert a 14-element feature tuple to a readable string."""


    if len(beyond_tuple) != 14:

        return str(beyond_tuple)


    carbonyl_names = {0: 'none', 1: 'ketone', 2: 'ester', 3: 'amide', 4: 'acid'}

    alpha_names = {0: 'none', 1: 'O', 2: 'N', 3: 'S', 4: 'X'}

    aromatic_names = {0: 'none', 1: 'benzene', 2: 'pyridine', 3: 'pyrimidine',

                      4: 'pyrrole', 5: 'indole', 6: 'imidazole', 7: 'furan/thio'}

    hyb_names = {0: 'unknown', 1: 'sp', 2: 'sp2', 3: 'sp3'}


    carb, ring, chg, da, alpha, arom, ally, vinyl, rel_chiral, ez, nH, hyb, d_hetero, d_arom = beyond_tuple


    parts = [

        f"carb={carbonyl_names.get(carb, carb)}",

        f"ring={ring}",

        f"chg={chg}",

        f"da={da}",

        f"alpha={alpha_names.get(alpha, alpha)}",

        f"arom={aromatic_names.get(arom, arom)}",

        f"ally={ally}",

        f"vinyl={vinyl}",

        f"rel_chiral={rel_chiral}",  # Relative configuration

        f"ez={ez}",

        f"nH={nH}",

        f"hyb={hyb_names.get(hyb, hyb)}",

        f"d_hetero={d_hetero}",  # Heteroatom distance

        f"d_arom={d_arom}"       # Aromatic distance

    ]


    return '|'.join(parts)


# Self-test


if __name__ == '__main__':

    import sys


    test_cases = [

        ("SDF original", "C[C@H]1C(=O)N/C(=C\\\C2=C(NC3=CC=CC=C32)C(C)(C)C=C)/C(=O)N1"),

        ("User @@", "C[C@@H]1C(=O)N/C(=C\\\C2=C(NC3=CC=CC=C32)C(C)(C)C=C)/C(=O)N1"),

        ("User /", "C[C@H]1C(=O)N/C(=C/C2=C(NC3=CC=CC=C32)C(C)(C)C=C)/C(=O)N1"),

    ]


    print("=" * 80)

    print("TEST DES FEATURES BEYOND + STÉRÉO")

    print("=" * 80)


    for label, smiles in test_cases:

        print(f"\n{label}:")

        print(f"  SMILES: {smiles}")


        mol = Chem.MolFromSmiles(smiles)

        if mol:

            mol_h = Chem.AddHs(mol)


            # Compute features for atom 24

            beyond = compute_beyond_with_stereo(mol_h, 24, radius_stereo=3)

            print(f"  Atome 24: {beyond_stereo_to_string(beyond)}")


            # Display stereochemistry-related features

            chiral, ez, total = beyond[8], beyond[9], beyond[10]

            print(f"  Stéréo: chiral={chiral}, ez={ez}, total={total}")


    print("\n" + "=" * 80)

    print("VÉRIFICATION: Les codes sont-ils DIFFÉRENTS?")

    print("=" * 80)


    results = {}

    for label, smiles in test_cases:

        mol = Chem.MolFromSmiles(smiles)

        if mol:

            mol_h = Chem.AddHs(mol)

            beyond = compute_beyond_with_stereo(mol_h, 24)

            results[label] = beyond


    # Compare results

    unique_tuples = set(results.values())

    print(f"\nNombre de tuples uniques: {len(unique_tuples)} / {len(results)}")


    if len(unique_tuples) == len(results):

        print("PARFAIT! Chaque stéréoisomère a un tuple unique!")

        print("  Le composite fingerprint séparera automatiquement les stéréoisomères")

    else:

        print("PROBLÈME: Certains stéréoisomères ont le même tuple")

        print("  Ils seront fusionnés (moyenne biaisée)")


    print("\nDétail:")

    for label, beyond in results.items():

        print(f"  {label}: {beyond[8:]}")  # Stereochemistry-related slice
