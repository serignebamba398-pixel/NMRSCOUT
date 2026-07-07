# nmr_plugins/dynamics.py - VERSION POST-PROCESSEUR PRODUCTION
# ============================================================
# POST-PROCESSEUR DYNAMIQUE pour RMN métabolomique
# Reçoit J_static, applique corrections dynamiques, retourne J_effectif
# ============================================================

"""
POST-PROCESSEUR PUR pour effets dynamiques RMN en métabolomique.

ARCHITECTURE:
- Reçoit J_static calculé par d'autres modules (extensions, core, etc.)
- Applique corrections dynamiques (échange, T, pH, solvant)
- Retourne J_effectif ou None
- NE CALCULE PAS de J théorique (pas de duplication)

UTILISATION:
    from dynamics import apply_dynamics_to_coupling
    J_static = coupling_extensions.get_coupling(...)  # ou autres modules
    J_final = apply_dynamics_to_coupling(mol, h1_idx, h2_idx, J_static, ph=7.4, solvent="D2O")
"""

from typing import Optional, Dict, Any, Tuple
from functools import lru_cache
import warnings
import math

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    raise ImportError("RDKit requis pour dynamics.py")

# Cache pour les parents d'hydrogènes (optimisation)
_parents_cache: Dict[int, Optional[int]] = {}

# =============================================================================
# CONSTANTES UNIQUES ET VALIDÉES
# =============================================================================

# Facteurs d'échange selon solvant (données expérimentales métabolomique)
SOLVENT_EXCHANGE_FACTORS: Dict[str, Dict[str, float]] = {
    'D2O': {  # Métabolomique standard (plasma, urine)
        'OH_exchange': 10.0,    # Échange ultra-rapide (< 1 ms)
        'NH2_exchange': 5.0,    # Amines primaires (rapide)
        'COOH_exchange': 8.0,   # Acides carboxyliques (rapide)
        'SH_exchange': 0.1,     # Thiols (lents, ~1 s)
        'imide_exchange': 0.05, # NH amide (très lent)
    },
    'CD3OD': {  # Extraction lipides, métabolites hydrophobes
        'OH_exchange': 1.0,     # Échange ralenti
        'NH2_exchange': 0.5,
        'COOH_exchange': 2.0,
        'SH_exchange': 0.05,
        'imide_exchange': 0.01,
    },
    'CDCl3': {  # Lipides, composés neutres
        'OH_exchange': 0.1,     # Échange très lent (> 1 s)
        'NH2_exchange': 0.1,
        'COOH_exchange': 0.5,
        'SH_exchange': 0.01,
        'imide_exchange': 0.005,
    },
    'DMSO-d6': {  # Composés peu solubles
        'OH_exchange': 0.2,
        'NH2_exchange': 0.3,
        'COOH_exchange': 0.8,
        'SH_exchange': 0.02,
        'imide_exchange': 0.01,
    },
}

# Effets de pH sur état de protonation (pKa métabolomiques)
PH_PROTONATION_EFFECTS: Dict[str, Dict[str, Any]] = {
    'carboxyl': {
        'pKa_range': (3.5, 5.5),  # Acides organiques
        'effect_below_pKa': 0.2,   # COOH (échange lent)
        'effect_above_pKa': 10.0,  # COO⁻ (échange ultra-rapide)
    },
    'amino': {
        'pKa_range': (9.0, 11.0),   # Amines aliphatiques
        'effect_below_pKa': 0.1,    # NH₃⁺ (lent)
        'effect_above_pKa': 5.0,    # NH₂ (rapide)
    },
    'imidazole': {
        'pKa_range': (6.0, 7.2),    # Histidine
        'effect_below_pKa': 0.5,
        'effect_above_pKa': 2.0,
    },
    'phosphate': {
        'pKa_range': (6.0, 7.4),    # ATP, ADP
        'effect_below_pKa': 0.3,
        'effect_above_pKa': 8.0,
    },
    'phenol': {
        'pKa_range': (9.5, 11.0),   # Tyrosine
        'effect_below_pKa': 0.2,
        'effect_above_pKa': 10.0,
    },
}

# Énergies d'activation pour processus dynamiques (kJ/mol)
DYNAMIC_ACTIVATION_ENERGIES: Dict[str, float] = {
    'ring_flip_pyranose': 45.0,      # Glucose, galactose, mannose
    'ring_flip_furanose': 35.0,      # Ribose, fructose
    'ring_flip_cyclohexane': 45.0,   # Stéroïdes, cyclitols
    'ring_flip_piperidine': 46.0,    # Pipéridine, morpholine
    'amide_rotation_peptide': 80.0,  # Liaisons peptidiques (N-CO)
    'amide_rotation_simple': 75.0,   # Amides secondaires
    'amide_rotation_formamide': 75.0,
    'amide_rotation_dmf': 18.0,      # Amides tertiaires (rare métabolomique)
    'keto_enol_tautomerism': 60.0,   # Pyruvate, acétoacétate
    'imine_enamine_tautomerism': 55.0,
    'nucleobase_tautomerism': 50.0,  # Purines, pyrimidines
}

# Facteurs de correction température pour différents processus dynamiques
# {température_K: facteur_correction}
TEMPERATURE_CORRECTION_FACTORS: Dict[str, Dict[int, float]] = {
    'sugar_ring_flipping': {
        273: 0.8,   # 0°C - Ring flip ralenti
        283: 0.7,   # 10°C
        298: 0.6,   # 25°C - Standard NMR
        310: 0.55,  # 37°C - Corps humain
        323: 0.5,   # 50°C
    },
    'amide_rotation': {
        273: 1.0,   # 0°C - Rotation très lente
        283: 1.0,   # 10°C
        298: 1.0,   # 25°C - Standard NMR (rotation lente)
        310: 0.95,  # 37°C
        323: 0.85,  # 50°C - Rotation plus rapide
        353: 0.6,   # 80°C - Rotation moyennée
    },
}

# =============================================================================
# SYSTÈMES RIGIDES - NE PAS RÉDUIRE LES COUPLAGES J
# =============================================================================
# Ces systèmes sont structurellement rigides et ne subissent pas de moyennage
# conformationnel significatif → J doit rester CONSTANT

RIGID_RING_SIZES = {3, 4}  # Cyclopropane, cyclobutane (contraintes de tension)

# Taille minimale pour considérer un système comme bicyclique ponté rigide
# NOTE: Norbornane a 2 cycles SSSR mais avec 3 atomes partagés = pont rigide
MIN_RINGS_FOR_BRIDGED_RIGID = 2  # ≥2 cycles avec pont = rigide
MIN_SHARED_ATOMS_FOR_BRIDGE = 3   # ≥3 atomes partagés entre cycles = pont

# SMARTS patterns pour détection rapide de systèmes rigides
RIGID_SYSTEM_SMARTS = {
    # Aromatiques
    'benzene': 'c1ccccc1',
    'pyridine': 'c1ccncc1',
    'pyrrole': 'c1cc[nH]c1',
    'furan': 'c1ccoc1',
    'thiophene': 'c1ccsc1',
    'imidazole': 'c1c[nH]cn1',
    'pyrazole': 'c1cc[nH]n1',
    
    # Bicycliques aromatiques
    'naphthalene': 'c1ccc2ccccc2c1',
    'indole': 'c1ccc2[nH]ccc2c1',
    'quinoline': 'c1ccc2ncccc2c1',
    
    # Petits cycles saturés rigides
    'cyclopropane': 'C1CC1',
    'cyclobutane': 'C1CCC1',
    'epoxide': 'C1OC1',
    'aziridine': 'C1NC1',
    
    # Bicycliques pontés
    'norbornane': 'C1CC2CCC1C2',
    'adamantane': 'C1C2CC3CC1CC(C2)C3',
    
    # Alcènes/alcynes conjugués rigides
    'allene': 'C=C=C',
    'butadiene_trans': 'C=CC=C',
    
    # Carbonyles conjugués
    'enone': 'C=CC=O',
    'acrolein': 'C=CC=O',
}

# Systèmes où la rigidité est TOUJOURS appliquée
ALWAYS_RIGID_PATTERNS = [
    'cyclopropane', 'cyclobutane', 'epoxide', 'aziridine',
    'norbornane', 'adamantane', 'allene',
]

# Températures de référence pour corrections
TEMPERATURE_REFERENCE = 298.15  # 25°C
TEMPERATURE_HUMAN_BODY = 310.15  # 37°C

# =============================================================================
# POST-PROCESSEUR DYNAMIQUE
# =============================================================================

class DynamicsPostProcessor:
    """
    POST-PROCESSEUR PUR: Applique les effets dynamiques INSTRUMENTAUX à J_static.
    
    Architecture:
    - Input: J_static (calculé par modules statiques)
    - Output: J_effectif (corrigé pour dynamique)
    - Pas de calcul de J théorique
    """
    
    def __init__(
        self,
        ph: float = 7.4,
        solvent: str = "D2O",
        temperature_k: float = 298.15,
        validate_params: bool = True
    ):
        """
        Initialise les conditions métabolomiques.
        
        Args:
            ph: pH de l'échantillon (7.4 = plasma humain)
            solvent: Solvant RMN (D2O, CD3OD, CDCl3, DMSO-d6)
            temperature_k: Température en Kelvin (298.15 = 25°C)
            validate_params: Valider en production (True)
        """
        if validate_params:
            self._validate_parameters(ph, solvent, temperature_k)
        
        self.ph = ph
        self.solvent = solvent
        self.temperature_k = temperature_k
        
        # Cache local pour les parents d'hydrogènes (optimisation)
        self._parents_cache: Dict[int, Optional[int]] = {}
    
    def _validate_parameters(self, ph: float, solvent: str, temperature_k: float):
        """Valide les paramètres de simulation."""
        # Validation basique du pH (plage chimique)
        if not (0.0 <= ph <= 14.0):
            raise ValueError(f"pH {ph} hors plage valide [0.0, 14.0]")
        
        if solvent not in SOLVENT_EXCHANGE_FACTORS:
            raise ValueError(f"Solvant '{solvent}' non supporté. Valides: {list(SOLVENT_EXCHANGE_FACTORS.keys())}")
        
        temp_c = temperature_k - 273.15
        if not (-20 <= temp_c <= 100):
            import warnings
            warnings.warn(f"Température {temp_c:.1f}°C hors plage métabolomique [-20, 100]°C", UserWarning)
    
    def apply_dynamics(
        self,
        mol: Chem.Mol,
        h1_idx: int,
        h2_idx: int,
        J_static: float
    ) -> Optional[float]:
        """
        Applique les corrections dynamiques INSTRUMENTALES à un couplage statique.
        
        Args:
            mol: Molécule RDKit (AVEC hydrogènes explicites)
            h1_idx, h2_idx: Indices des protons
            J_static: Couplage calculé par modules statiques (Hz)
        
        Returns:
            J effectif avec corrections dynamiques (Hz)
            None si J_static invalide (< 0.001 Hz)
        """
        # Validation d'entrée
        if J_static is None or abs(J_static) < 0.001:
            return None
        
        try:
            # ÉTAPE 0: VÉRIFICATION SYSTÈME RIGIDE → PAS DE RÉDUCTION
            # ============================================================
            # Les systèmes rigides (aromatiques, cycles tendus, bicycliques pontés,
            # alcènes/alcynes, conjugués) ne subissent PAS de moyennage dynamique
            if self._should_not_reduce_coupling(mol, h1_idx, h2_idx):
                return J_static  # ✅ SYSTÈME RIGIDE → J INCHANGÉ
            
            # ÉTAPE 1: Échange rapide → Suppression coupante
            if self._is_rapid_exchange(mol, h1_idx) or self._is_rapid_exchange(mol, h2_idx):
                # Le couplage n'est pas observable → signal singulet
                return 0.0
            
            # ÉTAPE 2: Calcul du facteur de correction cumulatif
            correction_factor = 1.0
            
            # Ring flipping (cycles 5-6 membres)
            ring_factor = self._get_ring_flipping_correction(mol, h1_idx, h2_idx)
            if ring_factor < 1.0:
                correction_factor *= ring_factor
            
            # Tautomérisme
            tautomer_factor = self._get_tautomer_correction(mol, h1_idx, h2_idx)
            if tautomer_factor < 1.0:
                correction_factor *= tautomer_factor
            
            # Rotation amide
            amide_factor = self._get_amide_rotation_correction(mol, h1_idx, h2_idx)
            if amide_factor < 1.0:
                correction_factor *= amide_factor
            
            # ÉTAPE 3: Application et validation
            J_effective = J_static * correction_factor
            
            # Arrondir à 3 décimales pour éviter artefacts numériques
            J_effective = round(J_effective, 3)
            
            return J_effective if abs(J_effective) > 0.001 else 0.0
            
        except Exception as e:
            # En cas d'erreur, retourner None pour éviter valeurs fausses
            warnings.warn(f"Erreur critique dans apply_dynamics: {e}. Retour None.", UserWarning)
            return None
    
    def _is_rapid_exchange(self, mol: Chem.Mol, h_idx: int) -> bool:
        """
        Détection AVANCÉE de protons échangeables rapides.
        
        Critères:
        - OH ionisé (pH > pKa)
        - NH₂ amine (pH > pKa)
        - COOH déprotéonné (pH > pKa)
        - SH toujours lent (conservé)
        """
        h_atom = mol.GetAtomWithIdx(h_idx)
        if h_atom.GetSymbol() != 'H':
            return False
        
        for neighbor in h_atom.GetNeighbors():
            parent = neighbor
            parent_sym = parent.GetSymbol()
            
            # === OXYGÈNE (OH, COOH) ===
            if parent_sym == 'O':
                factor = SOLVENT_EXCHANGE_FACTORS[self.solvent].get('OH_exchange', 1.0)
                
                # Vérifier si c'est un COOH
                is_carboxyl = False
                for nn in parent.GetNeighbors():
                    if nn.GetSymbol() == 'C':
                        # Chercher C=O
                        for nnn in nn.GetNeighbors():
                            if nnn.GetSymbol() == 'O' and nnn.GetIdx() != parent.GetIdx():
                                bond = mol.GetBondBetweenAtoms(nn.GetIdx(), nnn.GetIdx())
                                if bond and bond.GetBondTypeAsDouble() == 2.0:
                                    is_carboxyl = True
                                    factor = SOLVENT_EXCHANGE_FACTORS[self.solvent].get('COOH_exchange', 1.0)
                
                # Correction pH pour COOH
                if is_carboxyl:
                    if self.ph > 5.5:  # Au-dessus pKa → COO⁻
                        factor = PH_PROTONATION_EFFECTS['carboxyl']['effect_above_pKa']
                    else:
                        factor = PH_PROTONATION_EFFECTS['carboxyl']['effect_below_pKa']
                
                return factor > 5.0
            
            # === AZOTE (NH₂ amine, pas NH amide) ===
            elif parent_sym == 'N':
                if self._is_amide_nitrogen(mol, parent.GetIdx()):
                    return False  # NH amide → échange très lent
                
                factor = SOLVENT_EXCHANGE_FACTORS[self.solvent].get('NH2_exchange', 1.0)
                
                # Correction pH
                if self.ph < 7.0:  # pH acide → NH₃⁺
                    factor *= 0.1
                else:  # pH basique → NH₂
                    factor = PH_PROTONATION_EFFECTS['amino']['effect_above_pKa']
                
                return factor > 2.0
            
            # === SOUFRE (SH) ===
            elif parent_sym == 'S':
                # Toujours lent en métabolomique
                return False
        
        return False
    
    def _is_amide_nitrogen(self, mol: Chem.Mol, n_idx: int) -> bool:
        """Vérifie si N est dans liaison amide/peptide (échange LENT)."""
        n_atom = mol.GetAtomWithIdx(n_idx)
        if n_atom.GetSymbol() != 'N':
            return False
        
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                # Vérifier C=O
                for nn in neighbor.GetNeighbors():
                    if nn.GetSymbol() == 'O':
                        bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nn.GetIdx())
                        if bond and bond.GetBondTypeAsDouble() == 2.0:
                            return True
        
        return False
    
    def _get_ring_flipping_correction(self, mol: Chem.Mol, h1_idx: int, h2_idx: int) -> float:
        """
        Correction pour ring flipping (cycles 5-6 membres SATURÉS uniquement).
        
        IMPORTANT: Les cycles AROMATIQUES sont RIGIDES et ne subissent PAS de ring flipping!
        Seuls les cycles saturés (cyclohexane, pyranose, furanose) sont concernés.
        
        Returns: facteur (1.0 = pas d'effet, 0.3 = fort moyennage)
        """
        # Identifie atomes parents (avec cache local)
        h1_parent = self._get_hydrogen_parent_cached(mol, h1_idx)
        h2_parent = self._get_hydrogen_parent_cached(mol, h2_idx)
        
        if h1_parent is None or h2_parent is None:
            return 1.0
        
        # ✅ CORRECTION CRITIQUE: Vérifier si les atomes sont aromatiques
        atom1 = mol.GetAtomWithIdx(h1_parent)
        atom2 = mol.GetAtomWithIdx(h2_parent)
        
        # Si l'un des atomes est aromatique → PAS de ring flipping
        if atom1.GetIsAromatic() or atom2.GetIsAromatic():
            return 1.0  # ✅ Cycles aromatiques RIGIDES - pas de correction
        
        # Trouver cycles communs
        ring_info = mol.GetRingInfo()
        
        for ring in ring_info.AtomRings():
            if len(ring) in [5, 6]:
                if (h1_parent in ring) or (h2_parent in ring):
                    ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                    symbols = [a.GetSymbol() for a in ring_atoms]
                    
                    # ✅ CORRECTION: Vérifier que le cycle n'est PAS aromatique
                    is_aromatic_ring = any(a.GetIsAromatic() for a in ring_atoms)
                    if is_aromatic_ring:
                        return 1.0  # Cycle aromatique → rigide → pas de correction
                    
                    # Cycles SATURÉS avec O (sucres) ou C uniquement (cyclohexane)
                    if all(s in ['C', 'O', 'N'] for s in symbols):
                        # Interpolation température
                        temp_rounded = round(self.temperature_k)
                        factors = TEMPERATURE_CORRECTION_FACTORS.get('sugar_ring_flipping', {298: 0.6})
                        
                        # Trouver température la plus proche
                        closest_temp = min(factors.keys(), key=lambda t: abs(t - temp_rounded))
                        
                        return factors[closest_temp]
        
        return 1.0
    
    def _get_tautomer_correction(self, mol: Chem.Mol, h1_idx: int, h2_idx: int) -> float:
        """
        Correction pour tautomérisme (céto-énol, imine-énamine).
        
        IMPORTANT: Les cycles aromatiques stables (benzène, indole, pyridine, etc.)
        ne subissent PAS de tautomérisme rapide → pas de correction.
        
        Cas concernés:
        - β-dicarbonyles (acétoacétate)
        - Imines/énamines non aromatiques
        - Nucléobases (cas particulier)
        
        MVP: Facteur fixe 0.7 si système tautomère RÉEL détecté
        """
        h1_parent = self._get_hydrogen_parent_cached(mol, h1_idx)
        h2_parent = self._get_hydrogen_parent_cached(mol, h2_idx)
        
        if h1_parent is None or h2_parent is None:
            return 1.0
        
        atom1 = mol.GetAtomWithIdx(h1_parent)
        atom2 = mol.GetAtomWithIdx(h2_parent)
        
        # ✅ CORRECTION CRITIQUE: Les aromatiques dans des cycles STABLES ne tautomerisent pas
        # Indole, benzène, pyridine, etc. sont des systèmes aromatiques STABLES
        if atom1.GetIsAromatic() and atom2.GetIsAromatic():
            # Vérifier si les deux sont dans le MÊME cycle aromatique stable
            ring_info = mol.GetRingInfo()
            rings1 = [set(r) for r in ring_info.AtomRings() if h1_parent in r]
            rings2 = [set(r) for r in ring_info.AtomRings() if h2_parent in r]
            
            # S'ils partagent un cycle aromatique → système stable → pas de tautomérisme
            for r1 in rings1:
                for r2 in rings2:
                    if r1 & r2:  # Intersection non vide
                        shared_ring = r1 & r2
                        # Vérifier si tous les atomes du cycle partagé sont aromatiques
                        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in shared_ring):
                            return 1.0  # ✅ Cycle aromatique stable - pas de correction
        
        # Tautomérisme réel: β-dicarbonyles, imines non aromatiques
        # Détecter motif C=O...C=C-OH ou C=N...C=C-NH
        if str(atom1.GetHybridization()) == 'SP2' and str(atom2.GetHybridization()) == 'SP2':
            # Vérifier si c'est un vrai système tautomère (pas aromatique stable)
            if not atom1.GetIsAromatic() and not atom2.GetIsAromatic():
                # Vérifier hétéroatomes adjacents (O, N pour tautomérisme)
                if any(n.GetSymbol() in ['O', 'N'] for n in atom1.GetNeighbors()) or \
                   any(n.GetSymbol() in ['O', 'N'] for n in atom2.GetNeighbors()):
                    return 0.7  # Facteur pour tautomérisme réel
        
        return 1.0
    
    def _get_amide_rotation_correction(self, mol: Chem.Mol, h1_idx: int, h2_idx: int) -> float:
        """
        Correction pour rotation amide lente.
        
        Returns: 1.0 (lent) ou 0.8 (partiel à haute T)
        """
        for h_idx in [h1_idx, h2_idx]:
            h_atom = mol.GetAtomWithIdx(h_idx)
            if h_atom.GetSymbol() != 'H':
                continue
            
            for neighbor in h_atom.GetNeighbors():
                if neighbor.GetSymbol() == 'N' and self._is_amide_nitrogen(mol, neighbor.GetIdx()):
                    # NH amide détecté → rotation lente
                    temp_rounded = round(self.temperature_k)
                    factors = TEMPERATURE_CORRECTION_FACTORS.get('amide_rotation', {298: 1.0})
                    
                    closest_temp = min(factors.keys(), key=lambda t: abs(t - temp_rounded))
                    return factors[closest_temp]
        
        return 1.0
    
    def _get_hydrogen_parent_cached(self, mol: Chem.Mol, h_idx: int) -> Optional[int]:
        """
        Récupère parent d'hydrogène avec cache local (optimisation).
        
        Cache évite 4 appels répétés aux mêmes H dans apply_dynamics.
        """
        if h_idx in self._parents_cache:
            return self._parents_cache[h_idx]
        
        h_atom = mol.GetAtomWithIdx(h_idx)
        if h_atom.GetSymbol() != 'H':
            self._parents_cache[h_idx] = None
            return None
        
        parent = None
        for neighbor in h_atom.GetNeighbors():
            parent = neighbor.GetIdx()
            break
        
        self._parents_cache[h_idx] = parent
        return parent

    # =========================================================================
    # DÉTECTION DES SYSTÈMES RIGIDES - NE PAS RÉDUIRE LES COUPLAGES J
    # =========================================================================
    
    def _should_not_reduce_coupling(self, mol: Chem.Mol, h1_idx: int, h2_idx: int) -> bool:
        """
        Détecte si la paire H-H est dans un système RIGIDE qui NE DOIT PAS être réduit.
        
        Systèmes rigides:
        1. Cycles aromatiques (benzène, pyridine, pyrrole, indole, etc.)
        2. Petits cycles saturés tendus (cyclopropane, cyclobutane, époxyde)
        3. Systèmes bicycliques pontés (norbornane, adamantane)
        4. Alcènes/alcynes rigides (butadiène trans, allène, styrène)
        5. Carbonyles conjugués avec aromatiques
        6. Systèmes sp2 conjugués étendus
        
        Returns:
            True si système RIGIDE → NE PAS réduire J
            False si système FLEXIBLE → réduction possible
        """
        # Obtenir atomes parents
        h1_parent = self._get_hydrogen_parent_cached(mol, h1_idx)
        h2_parent = self._get_hydrogen_parent_cached(mol, h2_idx)
        
        if h1_parent is None or h2_parent is None:
            return False  # Par défaut, permettre réduction
        
        atom1 = mol.GetAtomWithIdx(h1_parent)
        atom2 = mol.GetAtomWithIdx(h2_parent)
        
        # === 1. VÉRIFIER AROMATIQUE (même cycle aromatique) ===
        if self._is_in_same_aromatic_ring(mol, h1_parent, h2_parent, atom1, atom2):
            return True  # RIGIDE - pas de réduction
        
        # === 2. VÉRIFIER CYCLE SATURÉ RIGIDE (3-4 membres) ===
        if self._is_in_strained_ring(mol, h1_parent, h2_parent, atom1, atom2):
            return True  # RIGIDE - cycles tendus
        
        # === 3. VÉRIFIER SYSTÈME BICYCLIQUE PONTÉ ===
        if self._is_bridged_rigid_system(mol, h1_parent, h2_parent):
            return True  # RIGIDE - bicyclique ponté
        
        # === 4. VÉRIFIER ALCÈNE/ALCYNE RIGIDE ===
        if self._is_rigid_unsaturated_system(mol, h1_parent, h2_parent, atom1, atom2):
            return True  # RIGIDE - liaison double/triple
        
        # === 5. VÉRIFIER CONJUGAISON AVEC AROMATIQUE ===
        if self._is_conjugated_with_aromatic(mol, h1_parent, h2_parent):
            return True  # RIGIDE - conjugué avec aromatique
        
        return False  # Système flexible - réduction possible
    
    def _is_in_same_aromatic_ring(self, mol: Chem.Mol, h1_parent: int, h2_parent: int,
                                   atom1: Chem.Atom, atom2: Chem.Atom) -> bool:
        """
        Vérifie si les deux atomes sont dans le MÊME cycle aromatique.
        
        Couvre: benzène, naphthalène, pyridine, pyrrole, furane, thiophène,
                imidazole, indole, quinoléine, etc.
        """
        # Si l'un n'est pas aromatique, pas de cycle aromatique commun
        if not (atom1.GetIsAromatic() and atom2.GetIsAromatic()):
            # Cas spécial: un seul est aromatique mais ils sont adjacents
            # (ex: H sur carbone sp3 adjacent à benzène - NOT rigide pour ce couplage)
            return False
        
        ring_info = mol.GetRingInfo()
        
        for ring in ring_info.AtomRings():
            if h1_parent in ring and h2_parent in ring:
                # Vérifier que TOUS les atomes du cycle sont aromatiques
                all_aromatic = all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
                if all_aromatic:
                    return True  # Dans même cycle aromatique complet
        
        # Vérifier aussi les cycles fusionnés (naphthalène, indole)
        # Si les deux sont aromatiques et partagent un système de cycles fusionnés
        rings_containing_1 = [set(r) for r in ring_info.AtomRings() if h1_parent in r]
        rings_containing_2 = [set(r) for r in ring_info.AtomRings() if h2_parent in r]
        
        for r1 in rings_containing_1:
            for r2 in rings_containing_2:
                # Cycles fusionnés partagent des atomes
                shared = r1 & r2
                if shared:
                    # Vérifier si le système fusionné est entièrement aromatique
                    combined = r1 | r2
                    if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in combined):
                        return True
        
        return False
    
    def _is_in_strained_ring(self, mol: Chem.Mol, h1_parent: int, h2_parent: int,
                              atom1: Chem.Atom, atom2: Chem.Atom) -> bool:
        """
        Vérifie si les atomes sont dans un petit cycle SATURÉ TENDU (3-4 membres).
        
        Couvre: cyclopropane, cyclobutane, époxyde, aziridine, β-lactame, β-lactone
        """
        if not (atom1.IsInRing() and atom2.IsInRing()):
            return False
        
        ring_info = mol.GetRingInfo()
        
        for ring in ring_info.AtomRings():
            ring_size = len(ring)
            
            # Seulement 3-4 membres (cycles tendus)
            if ring_size not in RIGID_RING_SIZES:
                continue
            
            if h1_parent in ring and h2_parent in ring:
                # Vérifier que c'est SATURÉ (pas de doubles liaisons dans le cycle)
                is_saturated = True
                for i in range(ring_size):
                    bond = mol.GetBondBetweenAtoms(ring[i], ring[(i+1) % ring_size])
                    if bond and bond.GetBondTypeAsDouble() >= 2.0:
                        is_saturated = False
                        break
                
                if is_saturated:
                    return True  # Cycle saturé tendu → RIGIDE
        
        return False
    
    def _is_bridged_rigid_system(self, mol: Chem.Mol, atom1_idx: int, atom2_idx: int) -> bool:
        """
        Détecte systèmes bicycliques PONTÉS rigides.
        
        Couvre: norbornane (bicyclo[2.2.1]), adamantane (tricycle), cubane, prismane,
                décaline trans, stéroïdes (partie rigide)
        
        Critères:
        1. Atome dans ≥2 cycles avec pont (≥3 atomes partagés entre cycles)
        2. Système avec ≥2 cycles qui partagent ≥3 atomes (pont)
        """
        ring_info = mol.GetRingInfo()
        atom_rings = list(ring_info.AtomRings())
        
        if len(atom_rings) < 2:
            return False  # Besoin d'au moins 2 cycles pour un pont
        
        # Compter combien de cycles contiennent chaque atome
        def count_rings(atom_idx):
            return sum(1 for ring in atom_rings if atom_idx in ring)
        
        rings_atom1 = count_rings(atom1_idx)
        rings_atom2 = count_rings(atom2_idx)
        
        # Si un atome est dans ≥2 cycles → potentiellement ponté
        if rings_atom1 >= MIN_RINGS_FOR_BRIDGED_RIGID or rings_atom2 >= MIN_RINGS_FOR_BRIDGED_RIGID:
            # Vérifier si c'est vraiment un pont (≥3 atomes partagés entre 2 cycles)
            for i, r1 in enumerate(atom_rings):
                for r2 in atom_rings[i+1:]:
                    shared = set(r1) & set(r2)
                    # Pont = ≥3 atomes partagés (norbornane-like)
                    if len(shared) >= MIN_SHARED_ATOMS_FOR_BRIDGE:
                        # Vérifier si nos atomes sont dans ce système
                        combined = set(r1) | set(r2)
                        if atom1_idx in combined or atom2_idx in combined:
                            return True
        
        # Vérifier aussi atomes voisins (tête de pont adjacente)
        atom1 = mol.GetAtomWithIdx(atom1_idx)
        atom2 = mol.GetAtomWithIdx(atom2_idx)
        
        for neighbor in atom1.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if count_rings(n_idx) >= MIN_RINGS_FOR_BRIDGED_RIGID:
                # Vérifier pont via voisin
                for i, r1 in enumerate(atom_rings):
                    if n_idx in r1:
                        for r2 in atom_rings[i+1:]:
                            if n_idx in r2:
                                shared = set(r1) & set(r2)
                                if len(shared) >= MIN_SHARED_ATOMS_FOR_BRIDGE:
                                    return True
        
        for neighbor in atom2.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if count_rings(n_idx) >= MIN_RINGS_FOR_BRIDGED_RIGID:
                for i, r1 in enumerate(atom_rings):
                    if n_idx in r1:
                        for r2 in atom_rings[i+1:]:
                            if n_idx in r2:
                                shared = set(r1) & set(r2)
                                if len(shared) >= MIN_SHARED_ATOMS_FOR_BRIDGE:
                                    return True
        
        # Décaline trans: 2 cycles fusionnés avec géométrie trans
        # Détection simplifiée: si 2 cycles de 6 partagent exactement 2 atomes adjacents
        six_rings = [set(r) for r in atom_rings if len(r) == 6]
        if len(six_rings) >= 2:
            for i, r1 in enumerate(six_rings):
                for r2 in six_rings[i+1:]:
                    shared = r1 & r2
                    if len(shared) == 2:
                        # Fusion trans → rigide
                        shared_list = list(shared)
                        bond = mol.GetBondBetweenAtoms(shared_list[0], shared_list[1])
                        if bond:  # Atomes de jonction adjacents
                            if atom1_idx in (r1 | r2) and atom2_idx in (r1 | r2):
                                return True
        
        return False
    
    def _is_rigid_unsaturated_system(self, mol: Chem.Mol, h1_parent: int, h2_parent: int,
                                      atom1: Chem.Atom, atom2: Chem.Atom) -> bool:
        """
        Vérifie si la paire est liée par une liaison double/triple RIGIDE.
        
        Couvre: alcènes (éthylène, styrène), alcynes, allène, cumulènes
        """
        # Liaison directe entre les atomes parents?
        bond = mol.GetBondBetweenAtoms(h1_parent, h2_parent)
        
        if bond:
            bond_order = bond.GetBondTypeAsDouble()
            
            # Liaison double ou triple → RIGIDE
            if bond_order >= 2.0:
                return True
            
            # Simple mais les deux sont sp2 conjugués → potentiellement rigide
            if str(atom1.GetHybridization()) == 'SP2' and str(atom2.GetHybridization()) == 'SP2':
                # Vérifier si fait partie d'un système conjugué étendu
                if self._is_part_of_conjugated_system(mol, h1_parent, h2_parent):
                    return True
        
        # Pas de liaison directe - vérifier si séparés par liaison double
        # Ex: vinyle où H-C=C-H (couplage trans/cis vicinal)
        try:
            path = Chem.GetShortestPath(mol, h1_parent, h2_parent)
        except:
            return False
        
        if path and len(path) <= 4:  # Couplage à courte distance
            for i in range(len(path) - 1):
                bond = mol.GetBondBetweenAtoms(path[i], path[i+1])
                if bond and bond.GetBondTypeAsDouble() >= 2.0:
                    # Présence de liaison double dans le chemin
                    # Vérifier que ce n'est pas un carbonyle isolé
                    atom_i = mol.GetAtomWithIdx(path[i])
                    atom_j = mol.GetAtomWithIdx(path[i+1])
                    
                    # C=C (alcène) → rigide
                    if atom_i.GetSymbol() == 'C' and atom_j.GetSymbol() == 'C':
                        return True
                    
                    # C≡C (alcyne) → rigide
                    if bond.GetBondTypeAsDouble() == 3.0:
                        return True
        
        return False
    
    def _is_part_of_conjugated_system(self, mol: Chem.Mol, atom1_idx: int, atom2_idx: int) -> bool:
        """
        Vérifie si les atomes font partie d'un système π conjugué étendu.
        
        Ex: butadiène (C=C-C=C), styrène (Ph-C=C), énone (C=C-C=O)
        """
        atom1 = mol.GetAtomWithIdx(atom1_idx)
        atom2 = mol.GetAtomWithIdx(atom2_idx)
        
        # Les deux doivent être sp2
        if not (str(atom1.GetHybridization()) == 'SP2' and str(atom2.GetHybridization()) == 'SP2'):
            return False
        
        # Vérifier alternance simple-double autour des atomes
        def has_adjacent_double_bond(atom):
            for bond in atom.GetBonds():
                if bond.GetBondTypeAsDouble() == 2.0:
                    return True
            return False
        
        return has_adjacent_double_bond(atom1) or has_adjacent_double_bond(atom2)
    
    def _is_conjugated_with_aromatic(self, mol: Chem.Mol, atom1_idx: int, atom2_idx: int) -> bool:
        """
        Vérifie si les atomes sont conjugués avec un système aromatique.
        
        Ex: styrène (Ph-CH=CH2), benzaldéhyde (Ph-CHO), acide benzoïque
        
        Ces systèmes sont RIGIDES car la conjugaison restreint la rotation.
        """
        atom1 = mol.GetAtomWithIdx(atom1_idx)
        atom2 = mol.GetAtomWithIdx(atom2_idx)
        
        # Si l'un est directement aromatique
        if atom1.GetIsAromatic() or atom2.GetIsAromatic():
            # Vérifier si l'autre est sp2 (conjugué)
            non_aromatic = atom2 if atom1.GetIsAromatic() else atom1
            
            if str(non_aromatic.GetHybridization()) == 'SP2':
                return True  # Conjugué avec aromatique → rigide
        
        # Vérifier via chemin court (max 3 liaisons vers aromatique)
        try:
            path = Chem.GetShortestPath(mol, atom1_idx, atom2_idx)
        except:
            return False
        
        if not path:
            return False
        
        # Chercher atome aromatique dans le voisinage du chemin
        for atom_idx in path:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetIsAromatic():
                return True
            
            # Vérifier voisins directs du chemin
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIsAromatic():
                    # Voisin aromatique + chemin de sp2 → conjugué
                    if str(atom.GetHybridization()) == 'SP2':
                        return True
        
        return False


# =============================================================================
# DIASTÉRÉOTOPES: FUSION OU SÉPARATION SELON DYNAMIQUE CONFORMATIONNELLE
# =============================================================================
# Les diastéréotopes (CH₂ près d'un centre chiral) peuvent apparaître en:
# - DEUX patterns distincts (régime lent: k << πΔν)
# - UN pattern moyenné (régime rapide: k >> πΔν)
#
# Critères de décision:
# 1. Température: T > 320K favorise moyennage
# 2. Rigidité: cycles, centres chiraux adjacents → séparation
# 3. Flexibilité: chaînes acycliques → possible moyennage
# 4. Taux d'échange: k_exchange vs seuil de coalescence
# =============================================================================

# Énergies d'activation pour rotation conformationnelle (kJ/mol)
CONFORMATIONAL_ACTIVATION_ENERGIES = {
    'acyclic_CH2_CH2': 12.0,       # Rotation libre éthane-like (~3 kcal/mol)
    'acyclic_CH2_CHR': 14.0,       # Légèrement plus contraint
    'acyclic_CH2_CR2': 16.0,       # Plus encombré
    'cyclopentane_ring': 25.0,     # Pseudorotation cyclopentane
    'cyclohexane_ring': 45.0,      # Ring flip cyclohexane
    'cyclohexane_ring_chair': 45.0,
    'cycloheptane_ring': 30.0,     # Plus flexible que 6
    'sugar_pyranose': 45.0,        # Comparable cyclohexane
    'sugar_furanose': 25.0,        # Plus flexible (5 membres)
    'piperidine': 46.0,            # N-inversion + ring flip
    'adjacent_to_chiral': 18.0,    # Contraint par chiralité
    'benzylic_CH2': 15.0,          # Adjacent aromatique
    'alpha_carbonyl': 20.0,        # Adjacent C=O
}

# Seuil de coalescence: k = π × Δν × √2
# Δν typique pour diastéréotopes: 10-100 Hz selon rigidité
TYPICAL_DELTA_NU_HZ = {
    'rigid_cycle': 50.0,           # Cyclohexane, sucres
    'semi_rigid': 30.0,            # Cyclopentane
    'flexible_acyclic': 15.0,      # Chaîne acyclique
    'adjacent_chiral': 40.0,       # Centre chiral à distance 1
    'distant_chiral': 20.0,        # Centre chiral à distance 2-3
}

# Constante de Boltzmann
R_GAS = 8.314  # J/(mol·K)


def _calculate_exchange_rate(
    activation_energy_kj: float,
    temperature_k: float,
    prefactor: float = 1e13  # Facteur pré-exponentiel typique (s⁻¹)
) -> float:
    """
    Calcule le taux d'échange conformationnel selon Arrhenius.
    
    k = A × exp(-Ea / RT)
    
    Args:
        activation_energy_kj: Énergie d'activation en kJ/mol
        temperature_k: Température en Kelvin
        prefactor: Facteur pré-exponentiel A (typiquement 10¹³ s⁻¹)
    
    Returns:
        Taux d'échange k en s⁻¹
    """
    if temperature_k <= 0:
        return 0.0
    
    Ea_J = activation_energy_kj * 1000  # kJ → J
    exponent = -Ea_J / (R_GAS * temperature_k)
    
    # Éviter overflow
    if exponent < -100:
        return 0.0
    if exponent > 100:
        return prefactor  # Saturation
    
    return prefactor * math.exp(exponent)


def _get_coalescence_threshold(delta_nu_hz: float) -> float:
    """
    Calcule le seuil de coalescence pour l'échange rapide.
    
    Critère de coalescence: k_c = π × Δν × √2
    
    Si k >> k_c → échange rapide → UN pattern moyenné
    Si k << k_c → échange lent → DEUX patterns distincts
    
    Args:
        delta_nu_hz: Différence de fréquence en Hz entre les deux signaux
    
    Returns:
        Seuil de coalescence k_c en s⁻¹
    """
    return math.pi * abs(delta_nu_hz) * math.sqrt(2)


def _classify_ch2_environment(
    mol: Chem.Mol,
    heavy_idx: int
) -> Tuple[str, float, float]:
    """
    Classifie l'environnement d'un CH₂ pour déterminer sa dynamique.
    
    Args:
        mol: Molécule RDKit avec H explicites
        heavy_idx: Index du carbone portant les deux H
    
    Returns:
        Tuple (type_env, activation_energy_kj, delta_nu_hz)
    """
    try:
        atom = mol.GetAtomWithIdx(heavy_idx)
        
        if atom.GetSymbol() != 'C':
            return ('unknown', 15.0, 20.0)
        
        # Vérifier si dans un cycle
        if atom.IsInRing():
            ring_info = mol.GetRingInfo()
            ring_sizes = [len(r) for r in ring_info.AtomRings() if heavy_idx in r]
            
            if ring_sizes:
                min_ring = min(ring_sizes)
                
                # Cycles rigides (3-4 membres)
                if min_ring in (3, 4):
                    return ('rigid_small_ring', 80.0, 80.0)  # Très rigide
                
                # Cyclopentane (pseudorotation rapide)
                if min_ring == 5:
                    # Vérifier si c'est un furanose
                    has_O_in_ring = any(
                        mol.GetAtomWithIdx(i).GetSymbol() == 'O'
                        for r in ring_info.AtomRings() if heavy_idx in r
                        for i in r
                    )
                    if has_O_in_ring:
                        return ('sugar_furanose', 25.0, 35.0)
                    return ('cyclopentane', 25.0, 25.0)
                
                # Cyclohexane (ring flip lent)
                if min_ring == 6:
                    # Vérifier si pyranose
                    has_O_in_ring = any(
                        mol.GetAtomWithIdx(i).GetSymbol() == 'O'
                        for r in ring_info.AtomRings() if heavy_idx in r
                        for i in r
                    )
                    if has_O_in_ring:
                        return ('sugar_pyranose', 45.0, 50.0)
                    return ('cyclohexane', 45.0, 50.0)
                
                # Cycles plus grands (plus flexibles)
                if min_ring >= 7:
                    return ('large_ring', 30.0, 30.0)
        
        # CH₂ acyclique - vérifier voisinage
        neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() != 'H']

        # Adjacent à un centre chiral?
        for nb in neighbors:
            if nb.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                return ('adjacent_chiral', 18.0, 40.0)

        # ✅ NOUVEAU: Adjacent à un centre prochiral quaternaire?
        # Ex: citrate CH2-C(OH)(COOH)(CH2COOH) - le C quaternaire est prochiral
        # Les protons restent diastéréotopiques même avec rotation libre
        for nb in neighbors:
            if nb.GetSymbol() == 'C':
                # Compter les voisins lourds du voisin
                nb_heavy_neighbors = [n for n in nb.GetNeighbors() if n.GetSymbol() != 'H']
                # C'est un carbone quaternaire sp3? (4 voisins lourds, pas de double liaison)
                if len(nb_heavy_neighbors) == 4:
                    has_double_bond = any(
                        bond.GetBondTypeAsDouble() >= 2
                        for bond in nb.GetBonds()
                    )
                    if not has_double_bond:
                        # Vérifier si les substituants sont asymétriques (prochiral)
                        try:
                            ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
                            sub_ranks = [ranks[n.GetIdx()] for n in nb_heavy_neighbors]
                            # ≥3 types différents = prochiral (ex: OH, COOH, CH2, CH2)
                            if len(set(sub_ranks)) >= 3:
                                return ('adjacent_prochiral_quaternary', 50.0, 80.0)
                            # ≥2 types différents et l'un n'est pas le CH2 lui-même
                            if len(set(sub_ranks)) >= 2:
                                ch2_rank = ranks[heavy_idx]
                                other_ranks = [r for r in sub_ranks if r != ch2_rank]
                                if len(set(other_ranks)) >= 1:
                                    return ('adjacent_prochiral_quaternary', 50.0, 80.0)
                        except Exception:
                            pass

        # Adjacent à un aromatique?
        for nb in neighbors:
            if nb.GetIsAromatic():
                return ('benzylic', 15.0, 25.0)

        # Adjacent à un carbonyle?
        for nb in neighbors:
            if nb.GetSymbol() == 'C':
                for nn in nb.GetNeighbors():
                    if nn.GetSymbol() == 'O':
                        bond = mol.GetBondBetweenAtoms(nb.GetIdx(), nn.GetIdx())
                        if bond and bond.GetBondTypeAsDouble() >= 1.5:
                            return ('alpha_carbonyl', 20.0, 30.0)
        
        # CH₂ acyclique standard (rotation libre)
        return ('acyclic_flexible', 12.0, 15.0)
        
    except Exception:
        return ('unknown', 15.0, 20.0)


def should_merge_diastereotopic_CH2(
    mol: Chem.Mol,
    heavy_idx: int,
    h_a_idx: int,
    h_b_idx: int,
    temperature_k: float = 298.15,
    spectrometer_mhz: float = 700.0
) -> Tuple[bool, str, Dict[str, Any]]:
    """
    Détermine si deux protons diastéréotopiques doivent être fusionnés (UN pattern)
    ou séparés (DEUX patterns) selon la dynamique conformationnelle.
    
    RÈGLE FONDAMENTALE:
    - Si k_exchange >> k_coalescence → FUSIONNER (échange rapide)
    - Si k_exchange << k_coalescence → SÉPARER (échange lent)
    
    Args:
        mol: Molécule RDKit avec H explicites
        heavy_idx: Index du carbone portant les deux H
        h_a_idx, h_b_idx: Indices des deux protons
        temperature_k: Température en Kelvin (défaut 298.15 = 25°C)
        spectrometer_mhz: Fréquence du spectromètre (affecte Δν en Hz)
    
    Returns:
        Tuple:
        - should_merge: True si UN pattern (échange rapide), False si DEUX patterns
        - reason: Explication textuelle
        - details: Dict avec k_exchange, k_coalescence, ratio, etc.
    """
    details = {
        'temperature_k': temperature_k,
        'spectrometer_mhz': spectrometer_mhz,
        'environment': 'unknown',
        'activation_energy_kj': 0.0,
        'k_exchange': 0.0,
        'k_coalescence': 0.0,
        'exchange_ratio': 0.0,
    }
    
    try:
        # 1. Classifier l'environnement du CH₂
        env_type, Ea_kj, delta_nu_base = _classify_ch2_environment(mol, heavy_idx)
        details['environment'] = env_type
        details['activation_energy_kj'] = Ea_kj
        
        # 2. Ajuster Δν selon la fréquence du spectromètre
        # Δν (Hz) = Δδ (ppm) × freq (MHz)
        # Pour diastéréotopes: Δδ typique 0.01-0.1 ppm
        delta_nu_hz = delta_nu_base * (spectrometer_mhz / 700.0)
        details['delta_nu_hz'] = delta_nu_hz
        
        # 3. Calculer le taux d'échange conformationnel
        k_exchange = _calculate_exchange_rate(Ea_kj, temperature_k)
        details['k_exchange'] = k_exchange
        
        # 4. Calculer le seuil de coalescence
        k_coalescence = _get_coalescence_threshold(delta_nu_hz)
        details['k_coalescence'] = k_coalescence
        
        # 5. Calculer le ratio d'échange
        if k_coalescence > 0:
            exchange_ratio = k_exchange / k_coalescence
        else:
            exchange_ratio = float('inf')
        details['exchange_ratio'] = exchange_ratio
        
        # =====================================================================
        # RÈGLES DE DÉCISION
        # =====================================================================
        
        # RÈGLE 0: Adjacent à un centre PROCHIRAL QUATERNAIRE → TOUJOURS SÉPARER
        # Ex: citrate CH2-C(OH)(COOH)(R) - le C quaternaire est prochiral
        # Les protons sont sur des faces différentes du centre prochiral
        # Même avec rotation libre, ils restent diastéréotopiques
        if env_type == 'adjacent_prochiral_quaternary':
            return (False, "Adjacent à centre prochiral quaternaire → DEUX patterns (diastéréotopes)", details)

        # RÈGLE 1: Systèmes RIGIDES → TOUJOURS SÉPARER
        # Petits cycles (3-4), bicycliques pontés, centres chiraux adjacents directs
        if env_type in ('rigid_small_ring',):
            return (False, "Cycle tendu rigide (3-4 membres) → DEUX patterns", details)
        
        # RÈGLE 2: Cycles 6 membres (cyclohexane, pyranose) → SÉPARATION à T < 350K
        # Le ring flip est lent à température ambiante
        if env_type in ('cyclohexane', 'sugar_pyranose'):
            if temperature_k < 350:
                return (False, f"Cycle 6 membres rigide à {temperature_k:.0f}K → DEUX patterns", details)
            else:
                return (True, f"Cycle 6 membres avec ring flip rapide à {temperature_k:.0f}K → UN pattern", details)
        
        # RÈGLE 3: Cyclopentane/furanose → FUSION possible dès 298K
        # Pseudorotation rapide (~10⁹ s⁻¹)
        if env_type in ('cyclopentane', 'sugar_furanose'):
            if exchange_ratio > 10:
                return (True, "Cycle 5 membres avec pseudorotation rapide → UN pattern", details)
            elif temperature_k > 310:
                return (True, f"Cycle 5 membres à T={temperature_k:.0f}K → UN pattern", details)
            else:
                return (False, "Cycle 5 membres mais Δν élevé → DEUX patterns", details)
        
        # RÈGLE 4: Adjacent à centre chiral → SÉPARATION (sauf T très élevée)
        if env_type == 'adjacent_chiral':
            if temperature_k < 340:
                return (False, "Adjacent à centre chiral → DEUX patterns", details)
            else:
                return (True, f"Adjacent à chiral mais T={temperature_k:.0f}K → moyennage partiel", details)
        
        # RÈGLE 5: CH₂ benzylique → dépend de l'asymétrie
        if env_type == 'benzylic':
            if exchange_ratio > 50:
                return (True, "CH₂ benzylique avec rotation rapide → UN pattern", details)
            else:
                return (False, "CH₂ benzylique avec environnement asymétrique → DEUX patterns", details)
        
        # RÈGLE 6: CH₂ acyclique flexible → FUSION si rotation libre
        if env_type == 'acyclic_flexible':
            if exchange_ratio > 100:
                return (True, "CH₂ acyclique avec rotation libre → UN pattern moyenné", details)
            elif temperature_k > 320:
                return (True, f"CH₂ acyclique à T={temperature_k:.0f}K → moyennage", details)
            else:
                # Même acyclique, si proche d'un centre chiral → séparation
                # Vérifier si un centre chiral est à distance 1-3
                return (False, "CH₂ acyclique mais environnement asymétrique → DEUX patterns", details)
        
        # RÈGLE 7: Grands cycles (≥7 membres) → plus flexibles
        if env_type == 'large_ring':
            if exchange_ratio > 20:
                return (True, "Grand cycle flexible → UN pattern", details)
            else:
                return (False, "Grand cycle mais contraint → DEUX patterns", details)
        
        # RÈGLE 8: Alpha-carbonyle → dépend de la température
        if env_type == 'alpha_carbonyl':
            if temperature_k > 330 or exchange_ratio > 30:
                return (True, "CH₂ α-carbonyle avec rotation → UN pattern", details)
            else:
                return (False, "CH₂ α-carbonyle contraint → DEUX patterns", details)
        
        # DÉFAUT: Comparer k_exchange vs k_coalescence
        if exchange_ratio > 10:
            return (True, f"Échange rapide (ratio={exchange_ratio:.1f}) → UN pattern", details)
        elif exchange_ratio < 0.1:
            return (False, f"Échange lent (ratio={exchange_ratio:.1f}) → DEUX patterns", details)
        else:
            # Zone intermédiaire → élargissement des raies (intermediate exchange)
            # Considérer comme séparé avec note
            return (False, f"Échange intermédiaire (ratio={exchange_ratio:.1f}) → DEUX patterns élargis", details)
        
    except Exception as e:
        warnings.warn(f"Erreur should_merge_diastereotopic_CH2: {e}", UserWarning)
        # Par défaut: séparer (plus conservateur)
        return (False, f"Erreur: {e} → défaut DEUX patterns", details)


def get_averaged_diastereotopic_shift(
    shift_a_ppm: float,
    shift_b_ppm: float,
    population_a: float = 0.5,
    population_b: float = 0.5
) -> float:
    """
    Calcule le shift moyenné pour deux protons diastéréotopiques en échange rapide.
    
    δ_obs = P_a × δ_a + P_b × δ_b
    
    Args:
        shift_a_ppm, shift_b_ppm: Shifts individuels des deux protons
        population_a, population_b: Populations des deux conformères (somme = 1)
    
    Returns:
        Shift moyenné en ppm
    """
    total_pop = population_a + population_b
    if total_pop <= 0:
        return (shift_a_ppm + shift_b_ppm) / 2.0
    
    return (population_a * shift_a_ppm + population_b * shift_b_ppm) / total_pop


def get_averaged_diastereotopic_J(
    J_list_a: list,
    J_list_b: list,
    population_a: float = 0.5,
    population_b: float = 0.5
) -> list:
    """
    Calcule les couplages moyennés pour deux protons diastéréotopiques en échange rapide.
    
    J_obs = P_a × J_a + P_b × J_b (pour chaque couplage)
    
    Args:
        J_list_a, J_list_b: Listes de couplages pour chaque proton
        population_a, population_b: Populations des deux conformères
    
    Returns:
        Liste de J moyennés
    """
    total_pop = population_a + population_b
    if total_pop <= 0:
        total_pop = 1.0
    
    # Aligner les listes par taille
    max_len = max(len(J_list_a), len(J_list_b))
    J_a = list(J_list_a) + [0.0] * (max_len - len(J_list_a))
    J_b = list(J_list_b) + [0.0] * (max_len - len(J_list_b))
    
    averaged = []
    for ja, jb in zip(J_a, J_b):
        j_avg = (population_a * ja + population_b * jb) / total_pop
        if abs(j_avg) > 0.5:  # Seuil de détection
            averaged.append(round(j_avg, 1))
    
    return averaged


# =============================================================================
# API PRINCIPALE (à utiliser dans sim_mult_rulebased.py)
# =============================================================================

def apply_dynamics_to_coupling(
    mol: Chem.Mol,
    h1_idx: int,
    h2_idx: int,
    J_static: float,
    ph: float = 7.4,
    solvent: str = "D2O",
    temperature_k: float = 298.15,
    use_cache: bool = True
) -> Optional[float]:
    """
    Fonction PRINCIPALE: applique la dynamique à un J statique.
    
    Usage recommandé dans orchestrateur:
        J_static = coupling_extensions.get_coupling(...)  # ou autres modules
        J_final = apply_dynamics_to_coupling(mol, h1_idx, h2_idx, J_static, ph=7.4)
    
    Args:
        mol: Molécule RDKit (AVEC hydrogènes explicites)
        h1_idx, h2_idx: Indices des protons
        J_static: Couplage calculé par modules statiques (Hz)
        ph: pH de la solution [1.0, 14.0]
        solvent: Solvant RMN (D2O, CD3OD, CDCl3, DMSO-d6)
        temperature_k: Température en Kelvin [273.15, 373.15]
        use_cache: Utiliser cache global (True en production, False en test)
    
    Returns:
        J effectif avec corrections dynamiques (Hz)
        None si J_static invalide ou erreur critique
    """
    if J_static is None or abs(J_static) < 0.001:
        return None
    
    if not use_cache:
        # Calcul direct (debug, tests)
        processor = DynamicsPostProcessor(ph, solvent, temperature_k)
        return processor.apply_dynamics(mol, h1_idx, h2_idx, J_static)
    
    # Calcul avec cache (production)
    try:
        mol_smiles = Chem.MolToSmiles(mol, canonical=True)
        cache_key = f"{mol_smiles}_{h1_idx}_{h2_idx}_{J_static:.3f}_{ph:.2f}_{solvent}_{temperature_k:.2f}"
    except Exception as e:
        warnings.warn(f"Erreur génération clé cache: {e}. Calcul sans cache.", UserWarning)
        processor = DynamicsPostProcessor(ph, solvent, temperature_k)
        return processor.apply_dynamics(mol, h1_idx, h2_idx, J_static)
    
    # Cache global simplifié (plus rapide que pickle)
    if not hasattr(apply_dynamics_to_coupling, '_cache'):
        apply_dynamics_to_coupling._cache = {}
    
    if cache_key in apply_dynamics_to_coupling._cache:
        return apply_dynamics_to_coupling._cache[cache_key]
    
    # Calcul et stockage
    processor = DynamicsPostProcessor(ph, solvent, temperature_k)
    J_final = processor.apply_dynamics(mol, h1_idx, h2_idx, J_static)
    
    # Limiter taille cache
    if len(apply_dynamics_to_coupling._cache) < 50000:
        apply_dynamics_to_coupling._cache[cache_key] = J_final
    
    return J_final


# =============================================================================
# TESTS UNITAIRES (déplacer dans tests/test_dynamics.py en production)
# =============================================================================

def _run_tests():
    """Tests rapides validation. Ne pas exécuter en production."""
    print("=" * 70)
    print("DYNAMICS POST-PROCESSOR - TESTS VALIDATION")
    print("=" * 70)
    
    # Test 1: Échange rapide (éthanol)
    print("\n🧪 Test 1: Éthanol (OH échangeable, pH=7.4, D2O)")
    mol = Chem.MolFromSmiles("CCO")
    mol = Chem.AddHs(mol)
    h_atoms = [i for i, atom in enumerate(mol.GetAtoms()) if atom.GetSymbol() == 'H']
    
    if len(h_atoms) >= 2:
        # Couplage OH-CH2 (devrait être supprimé)
        J_test = apply_dynamics_to_coupling(
            mol, h_atoms[4], h_atoms[2], J_static=6.0, ph=7.4, solvent="D2O", use_cache=False
        )
        print(f"  J statique: 6.0 Hz → J effectif: {J_test} Hz (attendu: 0.0)")
        assert J_test == 0.0, "❌ OH devrait être échangeable rapidement"
        print("  ✅ PASS")
    
    # Test 2: Ring flipping (glucose)
    print("\n🧪 Test 2: Glucose β-D (ring flipping, 25°C)")
    mol = Chem.MolFromSmiles("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O")
    mol = Chem.AddHs(mol)
    h_atoms = [i for i, atom in enumerate(mol.GetAtoms()) if atom.GetSymbol() == 'H']
    
    if len(h_atoms) >= 2:
        J_test = apply_dynamics_to_coupling(
            mol, h_atoms[0], h_atoms[1], J_static=8.0, ph=7.4, solvent="D2O", use_cache=False
        )
        print(f"  J statique: 8.0 Hz → J effectif: {J_test} Hz (attendu: < 8.0)")
        assert J_test < 8.0 and J_test > 0, "❌ Ring flipping devrait réduire J"
        print("  ✅ PASS")
    
    # Test 3: Cache fonctionnel
    print("\n🧪 Test 3: Cache global (50k entrées)")
    mol = Chem.MolFromSmiles("CC(N)C(O)=O")  # Alanine
    mol = Chem.AddHs(mol)
    h_atoms = [i for i, atom in enumerate(mol.GetAtoms()) if atom.GetSymbol() == 'H']
    
    if len(h_atoms) >= 2:
        J1 = apply_dynamics_to_coupling(mol, h_atoms[0], h_atoms[1], J_static=7.0, ph=7.4)
        J2 = apply_dynamics_to_coupling(mol, h_atoms[0], h_atoms[1], J_static=7.0, ph=7.4)
        assert J1 == J2, "❌ Cache ne retourne pas la même valeur"
        print("  ✅ Cache fonctionnel")
        
        # Stats cache
        if hasattr(apply_dynamics_to_coupling, '_cache'):
            print(f"  Cache size: {len(apply_dynamics_to_coupling._cache)}")
    
    print("\n" + "=" * 70)
    print("🎉 TOUS LES TESTS PASSÉS - Module prêt pour production")
    print("=" * 70)


def clear_dynamics_cache():
    """Nettoie le cache global (utilitaire pour tests/benchmark)."""
    if hasattr(apply_dynamics_to_coupling, '_cache'):
        apply_dynamics_to_coupling._cache.clear()
        print("✅ Cache dynamics nettoyé")


# =============================================================================
# MAIN (exécuter tests uniquement si script lancé directement)
# =============================================================================

if __name__ == "__main__":
    # Désactiver les tests en production
    import os
    if os.environ.get('RUN_DYNAMICS_TESTS', '0') == '1':
        _run_tests()
    else:
        print("=" * 70)
        print("DYNAMICS POST-PROCESSOR - MODULE PRÊT POUR PRODUCTION")
        print("=" * 70)
        print("📋 Pour lancer les tests: RUN_DYNAMICS_TESTS=1 python dynamics.py")
        print("\n⚡ Usage:")
        print("  J_static = coupling_extensions.get_coupling(...)")
        print("  J_final = apply_dynamics_to_coupling(mol, h1_idx, h2_idx, J_static, ph=7.4)")
        print("\n💡 Optimisations:")
        print("  - Cache 50k entrées")
        print("  - Gestion pH/solvant/température")
        print("  - Post-traitement seulement (pas de calcul J)")
        print("=" * 70)
