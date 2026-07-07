#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
beyond_features_stereo.py
=========================
Extension de beyond_features.py avec support de la stéréochimie.

PROBLÈME RÉSOLU:
Le fingerprint Morgan 2D ne distingue pas les stéréoisomères.
Résultat: fusion de molécules E/Z différentes → moyenne biaisée.

SOLUTION:
Ajouter 3 features stéréochimiques au tuple beyond:
  beyond_standard = (carbonyl, hetero, charge, donor, alpha, arom, allyl, vinyl)  # 8 features
  beyond_stereo = beyond_standard + (chiral_env, ez_env, total_stereo)            # 11 features

Composite fingerprint robuste:
  composite_fp = hash((morgan_fp, beyond_stereo))

Maintenant les stéréoisomères ont des fingerprints DIFFÉRENTS!
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcNumAtomStereoCenters
import hashlib

# Fonction utilitaire pour Gasteiger (autonome, copié de beyond_features.py)
def ensure_gasteiger_charges(mol):
    """
    Calcule les charges Gasteiger une seule fois par molécule.
    Idempotent: ne recalcule pas si déjà présent.

    Returns:
        True si les charges sont disponibles, False sinon
    """
    try:
        # Vérifier si déjà calculé
        mol.GetAtomWithIdx(0).GetProp('_GasteigerCharge')
        return True
    except:
        pass

    try:
        AllChem.ComputeGasteigerCharges(mol)
        return True
    except:
        return False

# =============================================================================
# FONCTIONS INTERNES POUR FEATURES 2D (copiées de beyond_features.py)
# =============================================================================
# Ces fonctions sont maintenant autonomes pour éviter le code mort

BEYOND_2D_AVAILABLE = True

def _get_gasteiger_charge(atom):
    """Recupere la charge Gasteiger d'un atome, 0.0 si non disponible."""
    try:
        import numpy as np
        charge = float(atom.GetProp('_GasteigerCharge'))
        if np.isnan(charge) or np.isinf(charge):
            return 0.0
        return charge
    except:
        return 0.0


def _detect_nearest_carbonyl_type(mol, atom_idx):
    """
    Detecte le type de C=O le plus proche (jusqu'a 3 liaisons).

    Returns:
        0 = aucun carbonyle proche
        1 = cetone/aldehyde (C=O sans O-R ni N-R)
        2 = ester (C=O avec O-R)
        3 = amide (C=O avec N-R)
        4 = acide carboxylique (C=O avec O-H)

    Impact sur shift CH3:
        - Ester:    ~2.02 ppm
        - Cetone:   ~2.15 ppm
        - Amide:    ~2.00 ppm
        - Acide:    ~2.10 ppm
    """
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
                    # Verifier si c'est OH ou O-R
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
                    best_type = 4  # Acide carboxylique
                elif has_single_O_ester:
                    best_type = 2  # Ester
                elif has_single_N:
                    best_type = 3  # Amide
                else:
                    best_type = 1  # Cetone/aldehyde

        # Continuer BFS vers les voisins lourds
        if dist < 3:
            for nb in c_atom.GetNeighbors():
                if nb.GetSymbol() != 'H':
                    queue.append((nb.GetIdx(), dist + 1))

    return best_type


def _detect_heterocycle_type(mol, atom_idx):
    """
    Detecte le type d'heterocycle aromatique.

    Returns:
        0 = non aromatique ou pas dans un heterocycle
        5 = cycle aromatique 5-membre (imidazole, pyrrole, furan, thiophene)
        6 = cycle aromatique 6-membre (pyridine, pyrimidine)

    Impact sur shift:
        - Imidazole H2: ~7.5-8.5 ppm
        - Pyrimidine H: ~8.5-9.5 ppm
    """
    atom = mol.GetAtomWithIdx(atom_idx)

    if not atom.GetIsAromatic():
        return 0

    ring_info = mol.GetRingInfo()
    aromatic_ring_sizes = []

    for ring in ring_info.AtomRings():
        if atom_idx in ring:
            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
            # Verifier si le cycle est aromatique ET contient un heteroatome
            if all(a.GetIsAromatic() for a in ring_atoms):
                has_hetero = any(a.GetSymbol() in ('N', 'O', 'S') for a in ring_atoms)
                if has_hetero:
                    aromatic_ring_sizes.append(len(ring))

    if not aromatic_ring_sizes:
        return 0

    # Retourner la taille du plus petit cycle heterocyclique aromatique
    return min(aromatic_ring_sizes)


def _detect_donor_acceptor_environment(mol, atom_idx):
    """
    Detecte l'environnement donneur/accepteur dans un rayon de 3 liaisons.

    Returns:
        0 = neutre
        1 = donneur fort (NH2, OH avec H)
        2 = donneur faible (NR2, OR)
        3 = accepteur faible (halogene)
        4 = accepteur fort (C=O, NO2, CN)

    Impact: effet mesomere sur systemes aromatiques
    """
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

            # Detecter groupes donneurs/accepteurs
            if symbol in ('N', 'O') and not c_atom.GetIsAromatic():
                # Verifier si dans un groupe accepteur (C=O, C=N)
                is_in_acceptor = False
                for nb in c_atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(current, nb.GetIdx())
                    if bond and bond.GetBondTypeAsDouble() >= 2.0 and nb.GetSymbol() in ('C', 'N', 'S'):
                        is_in_acceptor = True
                        break

                if is_in_acceptor:
                    acceptor_score += 2
                elif c_atom.GetTotalNumHs() > 0:
                    donor_score += 2  # NH, OH = donneur fort
                else:
                    donor_score += 1  # NR2, OR = donneur faible

            # Halogenes = accepteurs faibles
            if symbol in ('F', 'Cl', 'Br', 'I'):
                acceptor_score += 1

        # Continuer BFS
        if dist < 3:
            for nb in c_atom.GetNeighbors():
                if nb.GetSymbol() != 'H':
                    queue.append((nb.GetIdx(), dist + 1))

    # Classifier
    if acceptor_score > donor_score:
        return 4 if acceptor_score >= 3 else 3
    elif donor_score > acceptor_score:
        return 1 if donor_score >= 2 else 2
    return 0


def _detect_alpha_to_heteroatom(mol, atom_idx):
    """
    Detecte si l'atome est en position alpha d'un heteroatome.

    Returns:
        0 = pas en alpha
        1 = alpha-O (ether, ester)
        2 = alpha-N (amine, amide)
        3 = alpha-S (thioether)
        4 = alpha-halogen

    Impact majeur sur shift 1H et 13C:
        - Alpha-O: +3.5 ppm (13C), +0.5 ppm (1H)
        - Alpha-N: +2.5 ppm (13C), +0.3 ppm (1H)
    """
    atom = mol.GetAtomWithIdx(atom_idx)

    # Pour les H, regarder le parent
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
    """
    Detecte le type de systeme aromatique auquel l'atome est attache.

    Returns:
        0 = pas attache a un aromatique
        1 = benzene (6-membre, C seulement)
        2 = pyridine (6-membre, 1 N)
        3 = pyrimidine (6-membre, 2+ N)
        4 = pyrrole (5-membre, 1 N, NH)
        5 = indole (5+6 fusionne, NH)
        6 = imidazole (5-membre, 2 N)
        7 = furan/thiophene (5-membre, O ou S)

    Impact sur shift:
        - Indole C3-CH2: plus blinde que benzene-CH2
        - Pyridine-CH2: plus deblinde
    """
    atom = mol.GetAtomWithIdx(atom_idx)

    # Pour les H, regarder le parent
    if atom.GetSymbol() == 'H':
        parents = [n for n in atom.GetNeighbors()]
        if not parents:
            return 0
        atom = parents[0]

    # Chercher un voisin aromatique
    aromatic_neighbor = None
    for nb in atom.GetNeighbors():
        if nb.GetIsAromatic():
            aromatic_neighbor = nb
            break

    if aromatic_neighbor is None:
        # Verifier si l'atome lui-meme est aromatique
        if not atom.GetIsAromatic():
            return 0
        aromatic_neighbor = atom

    # Analyser le systeme aromatique
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

        # Verifier si fusionne (indole = 5+6)
        is_fused = False
        for other_ring in ring_info.AtomRings():
            if other_ring != ring:
                shared = set(ring) & set(other_ring)
                if len(shared) >= 2:
                    is_fused = True
                    break

        # Classification
        if ring_size == 5:
            if is_fused and n_count >= 1:
                # Verifier si c'est un indole (NH present)
                for a in ring_atoms:
                    if a.GetSymbol() == 'N' and a.GetTotalNumHs() > 0:
                        return 5  # Indole
            if n_count >= 2:
                return 6  # Imidazole
            if n_count == 1:
                return 4  # Pyrrole
            if o_count >= 1 or s_count >= 1:
                return 7  # Furan/Thiophene

        elif ring_size == 6:
            if n_count >= 2:
                return 3  # Pyrimidine
            if n_count == 1:
                return 2  # Pyridine
            if n_count == 0 and o_count == 0 and s_count == 0:
                return 1  # Benzene

    return 0


def _detect_allylic_position(mol, atom_idx):
    """
    Detecte si l'atome est en position allylique (a cote d'une C=C).

    Returns:
        0 = pas allylique
        1 = allylique simple (CH2-C=C ou CH3-C=C)
        2 = allylique conjugue (CH2-C=C-C=C)
        3 = bis-allylique (C=C-CH2-C=C)
        4 = prenyle (CH2-C=C(CH3)2)

    Impact sur shift:
        - Allylique: ~0.3-0.5 ppm vers champs faibles
        - Prenyle Ar-CH2: ~3.3-3.5 ppm
    """
    atom = mol.GetAtomWithIdx(atom_idx)

    # Pour les H, regarder le parent
    if atom.GetSymbol() == 'H':
        parents = [n for n in atom.GetNeighbors()]
        if not parents:
            return 0
        atom = parents[0]

    # L'atome doit etre un C sp3
    if atom.GetSymbol() != 'C':
        return 0
    if str(atom.GetHybridization()) != 'SP3':
        return 0

    # Chercher des voisins sp2 (C=C)
    sp2_neighbors = []
    for nb in atom.GetNeighbors():
        if nb.GetSymbol() == 'C' and str(nb.GetHybridization()) == 'SP2':
            # Verifier que c'est une vraie C=C (pas un aromatique)
            if not nb.GetIsAromatic():
                sp2_neighbors.append(nb)

    if len(sp2_neighbors) == 0:
        return 0

    if len(sp2_neighbors) >= 2:
        return 3  # Bis-allylique

    # Verifier si prenyle (C=C(CH3)2)
    sp2_nb = sp2_neighbors[0]
    methyl_count = 0
    for nb2 in sp2_nb.GetNeighbors():
        if nb2.GetIdx() != atom.GetIdx():
            if nb2.GetSymbol() == 'C' and nb2.GetTotalNumHs() == 3:
                methyl_count += 1

    if methyl_count >= 2:
        return 4  # Prenyle

    # Verifier conjugaison (C=C-C=C)
    for nb2 in sp2_nb.GetNeighbors():
        if nb2.GetIdx() != atom.GetIdx():
            if nb2.GetSymbol() == 'C' and str(nb2.GetHybridization()) == 'SP2':
                if not nb2.GetIsAromatic():
                    return 2  # Conjugue

    return 1  # Allylique simple


def _detect_vinylic_conjugation(mol, atom_idx):
    """
    Detecte si l'atome est sur une double liaison conjuguee a un aromatique.

    Returns:
        0 = pas de conjugaison vinylique
        1 = vinylique simple (C=CH2)
        2 = vinylique-aromatique (Ar-C=C)
        3 = vinylidene exocyclique (cycle-C=CH-)
        4 = styrene-like (Ar-CH=CH2)

    Impact sur shift:
        - Vinylidene sur DKP-indole: ~0.2 ppm deblindage du CH3 adjacent
    """
    atom = mol.GetAtomWithIdx(atom_idx)

    # Pour les H, regarder le parent
    if atom.GetSymbol() == 'H':
        parents = [n for n in atom.GetNeighbors()]
        if not parents:
            return 0
        parent = parents[0]
    else:
        parent = atom

    # Le parent doit etre un C
    if parent.GetSymbol() != 'C':
        return 0

    # Chercher dans un rayon de 4 liaisons
    # pour capturer conjugaison a travers amide (CH3-CH-C=O-N-C=C)
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

        # Chercher C=C
        for nb in current.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(current_idx, nb.GetIdx())
            if bond and bond.GetBondTypeAsDouble() >= 2.0:
                if nb.GetSymbol() == 'C' and not nb.GetIsAromatic():
                    has_vinyl = True

                    # Verifier si exocyclique (une extremite dans un cycle)
                    ring_info = mol.GetRingInfo()
                    current_in_ring = any(current_idx in r for r in ring_info.AtomRings())
                    nb_in_ring = any(nb.GetIdx() in r for r in ring_info.AtomRings())
                    if current_in_ring != nb_in_ring:
                        is_exocyclic = True

                    # Verifier conjugaison aromatique
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
        return 3  # Vinylidene exocyclique
    if has_aromatic_conj:
        return 4 if parent.GetTotalNumHs() >= 1 else 2  # Styrene-like ou vinyl-aromatique

    return 1  # Vinylique simple


# =============================================================================
# DÉTECTION DE STÉRÉOCHIMIE LOCALE AMÉLIORÉE
# =============================================================================

def get_stereo_environment(mol, atom_idx, radius=3):
    """
    Détecte TOUTE la stéréochimie dans l'environnement d'un atome.

    Returns:
        tuple: (chiral_signature, ez_signature, num_chiral, num_ez)
    """
    try:
        # Assigner la stéréochimie (crucial!)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    except:
        pass

    # Trouver l'environnement
    env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom_idx)

    atoms_in_env = set([atom_idx])
    bonds_in_env = set()

    for bond_idx in env:
        bond = mol.GetBondWithIdx(bond_idx)
        atoms_in_env.add(bond.GetBeginAtomIdx())
        atoms_in_env.add(bond.GetEndAtomIdx())
        bonds_in_env.add(bond_idx)

    # Collecter les centres chiraux
    chiral_list = []
    for idx in sorted(atoms_in_env):
        atom = mol.GetAtomWithIdx(idx)

        # Méthode 1: ChiralTag
        chiral_tag = atom.GetChiralTag()
        if chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            chiral_list.append('R')
        elif chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            chiral_list.append('S')
        elif atom.HasProp('_CIPCode'):
            # Méthode 2: CIP code
            cip = atom.GetProp('_CIPCode')
            if cip in ['R', 'S']:
                chiral_list.append(cip)

    # Collecter les doubles liaisons E/Z
    ez_list = []
    for bond_idx in sorted(bonds_in_env):
        bond = mol.GetBondWithIdx(bond_idx)

        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue

        # Vérifier stéréo
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
            # Vérifier si la double liaison PEUT avoir une stéréo
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()

            # Double liaison peut être E/Z si chaque C a 2 substituants différents
            begin_neighbors = [n for n in begin_atom.GetNeighbors() if n.GetIdx() != end_atom.GetIdx()]
            end_neighbors = [n for n in end_atom.GetNeighbors() if n.GetIdx() != begin_atom.GetIdx()]

            if len(begin_neighbors) >= 1 and len(end_neighbors) >= 1:
                # Potentiellement E/Z mais non spécifié → marquer comme '?'
                # (peut causer des problèmes de fusion)
                pass

    # Créer signatures
    chiral_sig = ''.join(chiral_list) if chiral_list else ''
    ez_sig = ''.join(ez_list) if ez_list else ''

    return (chiral_sig, ez_sig, len(chiral_list), len(ez_list))


def encode_relative_chirality(chiral_sig):
    """
    Encode la configuration RELATIVE entre centres chiraux.

    CHIMIE IMPORTANTE:
    - Énantiomères (R,R) et (S,S) → MÊME code (mêmes shifts NMR)
    - Énantiomères (R,S) et (S,R) → MÊME code (mêmes shifts NMR)
    - Diastéréomères (R,R) vs (R,S) → codes DIFFÉRENTS (shifts différents)

    Args:
        chiral_sig: String comme "RRS" ou "SS" ou "R" etc.

    Returns:
        int: 0 si pas de chiral
             1 si 1 seul centre (pas de distinction R vs S)
             2-999 selon le pattern relatif like/unlike

    Exemples:
        "R" ou "S" → 1 (un seul centre, indistinguable)
        "RR" ou "SS" → même code (like-like, énantiomères)
        "RS" ou "SR" → même code (like-unlike, énantiomères)
        "RRS" → code A (like-unlike)
        "RSS" → code B (unlike-like) ≠ code A (diastéréomères!)
    """
    if not chiral_sig:
        return 0

    chiral_list = list(chiral_sig)  # "RRS" → ['R', 'R', 'S']

    if len(chiral_list) == 1:
        return 1  # Un seul centre : R et S indistinguables en NMR achiral

    # 2+ centres : encoder les relations like/unlike
    relations = []
    for i in range(len(chiral_list) - 1):
        if chiral_list[i] == chiral_list[i + 1]:
            relations.append("like")
        else:
            relations.append("unlike")

    rel_string = "-".join(relations)
    rel_hash = int(hashlib.md5(rel_string.encode()).hexdigest(), 16)
    return (rel_hash % 997) + 2  # 2-999


def _distance_to_nearest_heteroatom(mol, atom_idx, max_dist=8):
    """
    Calcule la distance BFS au plus proche hétéroatome (O, N, S, halogène).

    USAGE:
    - Discrimine CH3-O-C (distance=1) vs CH3-C-C-C (distance=8+)
    - Discrimine CH2-NH2 (distance=1) vs CH2-alkyle (distance=8+)

    Args:
        mol: Molécule RDKit (doit avoir des H explicites)
        atom_idx: Index de l'atome (typiquement un H)
        max_dist: Distance maximale à explorer (8 bonds)

    Returns:
        int: Distance en nombre de bonds (0 si l'atome lui-même est hétéro,
             max_dist si aucun trouvé)
    """
    from collections import deque

    hetero_symbols = {'O', 'N', 'S', 'F', 'Cl', 'Br', 'I'}

    atom = mol.GetAtomWithIdx(atom_idx)
    if atom.GetSymbol() in hetero_symbols:
        return 0

    # BFS
    visited = {atom_idx}
    queue = deque([(atom_idx, 0)])  # (index, distance)

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

            # Vérifier si c'est un hétéroatome
            if neighbor.GetSymbol() in hetero_symbols:
                return dist + 1

            queue.append((neighbor_idx, dist + 1))

    return max_dist


def _distance_to_nearest_aromatic(mol, atom_idx, max_dist=8):
    """
    Calcule la distance BFS au plus proche atome aromatique.

    USAGE:
    - Discrimine CH2-Ph (benzylique, distance=1) vs CH2-alkyle (distance=8+)
    - Discrimine CH3 sur cycle aromatique (distance=1) vs CH3 alkyle (distance=8+)

    Args:
        mol: Molécule RDKit (doit avoir des H explicites)
        atom_idx: Index de l'atome (typiquement un H)
        max_dist: Distance maximale à explorer (8 bonds)

    Returns:
        int: Distance en nombre de bonds (0 si l'atome lui-même est aromatique,
             max_dist si aucun trouvé)
    """
    from collections import deque

    atom = mol.GetAtomWithIdx(atom_idx)
    if atom.GetIsAromatic():
        return 0

    # BFS
    visited = {atom_idx}
    queue = deque([(atom_idx, 0)])  # (index, distance)

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

            # Vérifier si c'est un atome aromatique
            if neighbor.GetIsAromatic():
                return dist + 1

            queue.append((neighbor_idx, dist + 1))

    return max_dist


def encode_stereo_environment(chiral_sig, ez_sig):
    """
    Encode les signatures stéréo en codes numériques.

    IMPORTANT: Utilise hashlib.md5() au lieu de hash() pour être DÉTERMINISTE
    entre différents processus Python (hash() utilise un seed aléatoire).

    Returns:
        tuple: (chiral_code, ez_code)
            0 = pas de stéréo
            1-999 = différentes configurations

    NOTE: Cette fonction garde l'ancien comportement (config absolue).
    Pour config relative, utiliser encode_relative_chirality() directement.
    """
    if not chiral_sig:
        chiral_code = 0
    else:
        # Hash déterministe avec MD5
        chiral_hash = int(hashlib.md5(chiral_sig.encode()).hexdigest(), 16)
        chiral_code = (chiral_hash % 999) + 1

    if not ez_sig:
        ez_code = 0
    else:
        # Hash déterministe avec MD5
        ez_hash = int(hashlib.md5(ez_sig.encode()).hexdigest(), 16)
        ez_code = (ez_hash % 999) + 1

    return (chiral_code, ez_code)


# =============================================================================
# FONCTION PRINCIPALE ÉTENDUE
# =============================================================================

def compute_beyond_with_stereo(mol, atom_idx, radius_stereo=6, include_chiral=True):
    """
    Calcule les features beyond + features stéréochimiques + features discriminantes.

    Args:
        radius_stereo: Rayon pour détecter la stéréo (défaut: 6)
        include_chiral: Si False, relative_config_code=0 (fusionne énantiomères R/S)
                       True = exact-match, False = composite fingerprint

    Returns:
        tuple (12 éléments):
            [0-7]: Features 2D classiques (beyond_features.py)
                0. carbonyl_type
                1. hetero_ring
                2. charge_bin
                3. donor_acceptor
                4. alpha_hetero
                5. aromatic_type
                6. allylic_pos
                7. vinyl_conj

            [8-9]: Features stéréochimiques (ÉTAPE 3)
                8. relative_config_code: Config relative (like/unlike entre centres)
                   → (R,R) et (S,S) = même code (énantiomères)
                   → (R,R) et (R,S) = codes différents (diastéréomères)
                9. ez_code: Config E/Z (0 ou 1-999) ← TOUJOURS gardé

            [10-11]: Features discriminantes (ÉTAPE 1)
               10. nH: Nombre d'hydrogènes (0, 1, 2, 3)
               11. hybridization_code: SP=1, SP2=2, SP3=3

    CHIMIE IMPORTANTE:
        - Énantiomères (R vs S seul) → MÊMES shifts en solvant achiral
        - Diastéréomères (E vs Z) → shifts DIFFÉRENTS
        → Pour composite: include_chiral=False (fusionne R/S, sépare E/Z)

    Usage:
        # Composite fingerprint (fusionne énantiomères, sépare E/Z)
        beyond = compute_beyond_with_stereo(mol, idx, include_chiral=False)
        composite_fp = (morgan_fp, beyond)
    """
    atom = mol.GetAtomWithIdx(atom_idx)

    # Pour les H, utiliser le parent
    if atom.GetSymbol() == 'H':
        parents = list(atom.GetNeighbors())
        if parents:
            parent_idx = parents[0].GetIdx()
        else:
            return (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    else:
        parent_idx = atom_idx

    # =========================================================================
    # PARTIE 1: Features 2D classiques (8 features)
    # =========================================================================

    if BEYOND_2D_AVAILABLE:
        # Utiliser les fonctions existantes
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
        # Fallback si beyond_features pas disponible
        carbonyl_type = 0
        hetero_ring = 0
        charge_bin = 0
        donor_acceptor = 0
        alpha_hetero = 0
        aromatic_type = 0
        allylic_pos = 0
        vinyl_conj = 0

    # =========================================================================
    # PARTIE 2: Features stéréochimiques (2 features) - ÉTAPE 3
    # =========================================================================

    try:
        chiral_sig, ez_sig, num_chiral, num_ez = get_stereo_environment(mol, atom_idx, radius_stereo)

        # ÉTAPE 3: Configuration RELATIVE (pas absolue)
        # → Énantiomères (R,R) et (S,S) ont le même code
        # → Diastéréomères (R,R) vs (R,S) ont des codes différents
        relative_config_code = encode_relative_chirality(chiral_sig)

        # CHIMIE: Ignorer R/S pour composite (énantiomères = mêmes shifts)
        if not include_chiral:
            relative_config_code = 0

        # E/Z : encoder normalement (diastéréomères → shifts différents)
        if not ez_sig:
            ez_code = 0
        else:
            ez_hash = int(hashlib.md5(ez_sig.encode()).hexdigest(), 16)
            ez_code = (ez_hash % 999) + 1

    except Exception:
        relative_config_code = 0
        ez_code = 0

    # =========================================================================
    # PARTIE 3: Features discriminantes (2 features) - ÉTAPE 1
    # =========================================================================
    # Ces features garantissent que CH3, CH2, CH et C quaternaire ne fusionnent jamais,
    # et qu'un CH2 sp3 (alkyle) ne fusionne pas avec un CH2 sp2 (vinyl)

    # Feature 10: Nombre d'hydrogènes attachés (0, 1, 2, 3)
    atom_for_features = mol.GetAtomWithIdx(parent_idx)
    nH = sum(1 for n in atom_for_features.GetNeighbors() if n.GetSymbol() == 'H')

    # Feature 11: Code d'hybridation (SP=1, SP2=2, SP3=3)
    hyb_str = str(atom_for_features.GetHybridization())
    if 'SP3' in hyb_str:
        hybridization_code = 3
    elif 'SP2' in hyb_str:
        hybridization_code = 2
    elif 'SP' in hyb_str:
        hybridization_code = 1
    else:
        hybridization_code = 0

    # =========================================================================
    # PARTIE 4: Features de distance (2 features) - ÉTAPE 4
    # =========================================================================
    # Ces features discriminent les environnements selon la proximité d'hétéroatomes
    # et d'atomes aromatiques, permettant de séparer :
    # - CH3-O-C vs CH3-C-C-C (acétate vs alkyle)
    # - CH2-Ph vs CH2-alkyle (benzylique vs alkyle)
    # - CH2-NH2 vs CH2-alkyle (amine vs alkyle)

    # Feature 12: Distance au plus proche hétéroatome (O, N, S, halogène)
    dist_heteroatom = _distance_to_nearest_heteroatom(mol, parent_idx, max_dist=8)

    # Feature 13: Distance au plus proche atome aromatique
    dist_aromatic = _distance_to_nearest_aromatic(mol, parent_idx, max_dist=8)

    # =========================================================================
    # RETOUR: 14 features au total (étapes 1+3+4)
    # =========================================================================
    # Étape 1: +nH +hybridization (11→13 features)
    # Étape 3: relative_config_code remplace chiral_code, -total_stereo (13→12)
    # Étape 4: +dist_heteroatom +dist_aromatic (12→14 features)

    return (
        carbonyl_type,          # 0
        hetero_ring,            # 1
        charge_bin,             # 2
        donor_acceptor,         # 3
        alpha_hetero,           # 4
        aromatic_type,          # 5
        allylic_pos,            # 6
        vinyl_conj,             # 7
        relative_config_code,   # 8 (MODIFIÉ étape 3: config relative, pas absolue)
        ez_code,                # 9
        nH,                     # 10 (NOUVEAU étape 1)
        hybridization_code,     # 11 (NOUVEAU étape 1)
        dist_heteroatom,        # 12 (NOUVEAU étape 4)
        dist_aromatic           # 13 (NOUVEAU étape 4)
    )


def compute_beyond_with_stereo_safe(mol, atom_idx, radius_stereo=6, include_chiral=True):
    """
    Version safe qui retourne (0,)*14 en cas d'erreur.
    """
    try:
        return compute_beyond_with_stereo(mol, atom_idx, radius_stereo, include_chiral)
    except Exception:
        return (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)


# =============================================================================
# UTILITIES
# =============================================================================

def beyond_stereo_to_string(beyond_tuple):
    """
    Convertit le tuple (14 features) en string lisible.

    Example:
        (2, 0, 0, 0, 1, 5, 4, 0, 180, 0, 3, 3, 1, 8) →
        "carb=ester|ring=0|chg=0|da=0|alpha=O|arom=indole|ally=prenyl|vinyl=0|rel_chiral=180|ez=0|nH=3|hyb=sp3|d_hetero=1|d_arom=8"
    """
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
        f"rel_chiral={rel_chiral}",  # Config relative
        f"ez={ez}",
        f"nH={nH}",
        f"hyb={hyb_names.get(hyb, hyb)}",
        f"d_hetero={d_hetero}",  # Distance hétéroatome (étape 4)
        f"d_arom={d_arom}"       # Distance aromatique (étape 4)
    ]

    return '|'.join(parts)


# =============================================================================
# TEST
# =============================================================================

if __name__ == '__main__':
    import sys

    test_cases = [
        ("SDF original", "C[C@H]1C(=O)N/C(=C\\C2=C(NC3=CC=CC=C32)C(C)(C)C=C)/C(=O)N1"),
        ("User @@", "C[C@@H]1C(=O)N/C(=C\\C2=C(NC3=CC=CC=C32)C(C)(C)C=C)/C(=O)N1"),
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

            # Calculer pour atome 24 (CH3)
            beyond = compute_beyond_with_stereo(mol_h, 24, radius_stereo=3)
            print(f"  Atome 24: {beyond_stereo_to_string(beyond)}")

            # Afficher juste les features stéréo
            chiral, ez, total = beyond[8], beyond[9], beyond[10]
            print(f"  → Stéréo: chiral={chiral}, ez={ez}, total={total}")

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

    # Comparer
    unique_tuples = set(results.values())
    print(f"\nNombre de tuples uniques: {len(unique_tuples)} / {len(results)}")

    if len(unique_tuples) == len(results):
        print("✓ PARFAIT! Chaque stéréoisomère a un tuple unique!")
        print("  → Le composite fingerprint séparera automatiquement les stéréoisomères")
    else:
        print("✗ PROBLÈME: Certains stéréoisomères ont le même tuple")
        print("  → Ils seront fusionnés (moyenne biaisée)")

    print("\nDétail:")
    for label, beyond in results.items():
        print(f"  {label}: {beyond[8:]}")  # Juste les 3 features stéréo
