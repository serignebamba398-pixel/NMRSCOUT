# nmr_plugins/dynamics.py - VERSION POST-PROCESSEUR PRODUCTION
# ============================================================
# Dynamic post-processing for metabolomic NMR
# Apply dynamic corrections to J_static and return the effective coupling.
# ============================================================

"""Apply dynamic corrections to static NMR couplings.

Static couplings come from other modules. This module applies exchange,
temperature, pH, and solvent corrections and returns an effective J or None.
It does not calculate theoretical static couplings.

Usage:
    from dynamics import apply_dynamics_to_coupling
    J_static = coupling_extensions.get_coupling(...)
    J_final = apply_dynamics_to_coupling(mol, h1_idx, h2_idx, J_static, ph=7.4, solvent="D2O")"""

from typing import Optional, Dict, Any, Tuple
from functools import lru_cache
import warnings
import math

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    raise ImportError("RDKit requis pour dynamics.py")

# Cache hydrogen parent atoms.
_parents_cache: Dict[int, Optional[int]] = {}

# =============================================================================
# Shared constants
# =============================================================================

# Solvent exchange factors from experimental metabolomics data
SOLVENT_EXCHANGE_FACTORS: Dict[str, Dict[str, float]] = {
    'D2O': {  # Standard metabolomics (plasma, urine)
        'OH_exchange': 10.0,    # Very fast exchange (< 1 ms)
        'NH2_exchange': 5.0,    # Amines primaires (rapide)
        'COOH_exchange': 8.0,   # Acides carboxyliques (rapide)
        'SH_exchange': 0.1,     # Thiols (lents, ~1 s)
        'imide_exchange': 0.05, # Amide NH (very slow)
    },
    'CD3OD': {  # Lipid extraction and hydrophobic metabolites
        'OH_exchange': 1.0,     # Slower exchange
        'NH2_exchange': 0.5,
        'COOH_exchange': 2.0,
        'SH_exchange': 0.05,
        'imide_exchange': 0.01,
    },
    'CDCl3': {  # Lipids and neutral compounds
        'OH_exchange': 0.1,     # Very slow exchange (> 1 s)
        'NH2_exchange': 0.1,
        'COOH_exchange': 0.5,
        'SH_exchange': 0.01,
        'imide_exchange': 0.005,
    },
    'DMSO-d6': {  # Poorly soluble compounds
        'OH_exchange': 0.2,
        'NH2_exchange': 0.3,
        'COOH_exchange': 0.8,
        'SH_exchange': 0.02,
        'imide_exchange': 0.01,
    },
}

# pH effects on protonation state (metabolite pKa values)
PH_PROTONATION_EFFECTS: Dict[str, Dict[str, Any]] = {
    'carboxyl': {
        'pKa_range': (3.5, 5.5),  # Acides organiques
        'effect_below_pKa': 0.2,   # COOH (slow exchange)
        'effect_above_pKa': 10.0,  # COO⁻ (very fast exchange)
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

# Activation energies for dynamic processes (kJ/mol)
DYNAMIC_ACTIVATION_ENERGIES: Dict[str, float] = {
    'ring_flip_pyranose': 45.0,      # Glucose, galactose, mannose
    'ring_flip_furanose': 35.0,      # Ribose, fructose
    'ring_flip_cyclohexane': 45.0,   # Steroids and cyclitols
    'ring_flip_piperidine': 46.0,    # Piperidine and morpholine
    'amide_rotation_peptide': 80.0,  # Liaisons peptidiques (N-CO)
    'amide_rotation_simple': 75.0,   # Amides secondaires
    'amide_rotation_formamide': 75.0,
    'amide_rotation_dmf': 18.0,      # Tertiary amides (rare in metabolomics)
    'keto_enol_tautomerism': 60.0,   # Pyruvate and acetoacetate
    'imine_enamine_tautomerism': 55.0,
    'nucleobase_tautomerism': 50.0,  # Purines, pyrimidines
}

# Temperature correction factors for dynamic processes
# {temperature_K: correction_factor}
TEMPERATURE_CORRECTION_FACTORS: Dict[str, Dict[int, float]] = {
    'sugar_ring_flipping': {
        273: 0.8,   # 0°C - Ring flip ralenti
        283: 0.7,   # 10°C
        298: 0.6,   # 25°C - Standard NMR
        310: 0.55,  # 37°C - Corps humain
        323: 0.5,   # 50°C
    },
    'amide_rotation': {
        273: 1.0,   # 0°C: very slow rotation
        283: 1.0,   # 10°C
        298: 1.0,   # 25°C - Standard NMR (rotation lente)
        310: 0.95,  # 37°C
        323: 0.85,  # 50°C - Rotation plus rapide
        353: 0.6,   # 80°C: averaged rotation
    },
}

# =============================================================================
# Preserve J couplings in rigid systems.
# =============================================================================
# These structurally rigid systems do not undergo significant
# conformational averaging, so J remains constant.

RIGID_RING_SIZES = {3, 4}  # Cyclopropane and cyclobutane (ring strain)

# Minimum size for classification as a rigid bridged bicyclic system
# Norbornane has two SSSR rings with three shared atoms forming a rigid bridge.
MIN_RINGS_FOR_BRIDGED_RIGID = 2  # At least two bridged rings indicate rigidity.
MIN_SHARED_ATOMS_FOR_BRIDGE = 3   # At least three shared ring atoms indicate a bridge.

# SMARTS patterns for rigid systems
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
    
    # Small rigid saturated rings
    'cyclopropane': 'C1CC1',
    'cyclobutane': 'C1CCC1',
    'epoxide': 'C1OC1',
    'aziridine': 'C1NC1',
    
    # Bridged bicyclic systems
    'norbornane': 'C1CC2CCC1C2',
    'adamantane': 'C1C2CC3CC1CC(C2)C3',
    
    # Rigid conjugated alkenes and alkynes
    'allene': 'C=C=C',
    'butadiene_trans': 'C=CC=C',
    
    # Conjugated carbonyls
    'enone': 'C=CC=O',
    'acrolein': 'C=CC=O',
}

# Systems always treated as rigid
ALWAYS_RIGID_PATTERNS = [
    'cyclopropane', 'cyclobutane', 'epoxide', 'aziridine',
    'norbornane', 'adamantane', 'allene',
]

# Reference temperatures for corrections
TEMPERATURE_REFERENCE = 298.15  # 25°C
TEMPERATURE_HUMAN_BODY = 310.15  # 37°C

# =============================================================================
# POST-PROCESSEUR DYNAMIQUE
# =============================================================================

class DynamicsPostProcessor:
    """Apply instrumental dynamic effects to J_static.

    Static modules supply J_static; this class returns a dynamically corrected
    coupling without recalculating the theoretical static value."""
    
    def __init__(
        self,
        ph: float = 7.4,
        solvent: str = "D2O",
        temperature_k: float = 298.15,
        validate_params: bool = True
    ):
        """Initialize metabolomic conditions.

        Args:
            ph: Sample pH (7.4 for human plasma).
            solvent: NMR solvent (D2O, CD3OD, CDCl3, DMSO-d6).
            temperature_k: Temperature in kelvin (298.15 = 25°C).
            validate_params: Enable parameter validation for production use."""
        if validate_params:
            self._validate_parameters(ph, solvent, temperature_k)
        
        self.ph = ph
        self.solvent = solvent
        self.temperature_k = temperature_k
        
        # Local cache of hydrogen parent atoms
        self._parents_cache: Dict[int, Optional[int]] = {}
    
    def _validate_parameters(self, ph: float, solvent: str, temperature_k: float):
        """Validate simulation parameters."""
        # Check pH against the supported chemical range.
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
        """Apply instrumental dynamic corrections to a static coupling.

        Args:
            mol: RDKit molecule with explicit hydrogens.
            h1_idx, h2_idx: Proton indices.
            J_static: Coupling from static modules, in Hz.

        Returns:
            Effective coupling in Hz, or None for invalid J_static (< 0.001 Hz)."""
        # Validate inputs.
        if J_static is None or abs(J_static) < 0.001:
            return None
        
        try:
            # Preserve couplings in rigid systems.
            # ============================================================
            # Rigid aromatic, strained-ring, bridged bicyclic, and conjugated
            # alkene/alkyne systems are not subject to dynamic averaging here.
            if self._should_not_reduce_coupling(mol, h1_idx, h2_idx):
                return J_static  # Rigid system: leave J unchanged.
            
            # Suppress couplings under fast exchange.
            if self._is_rapid_exchange(mol, h1_idx) or self._is_rapid_exchange(mol, h2_idx):
                # The coupling is not observable; the signal is a singlet.
                return 0.0
            
            # Calculate the cumulative correction factor.
            correction_factor = 1.0
            
            # Ring flipping (cycles 5-6 membres)
            ring_factor = self._get_ring_flipping_correction(mol, h1_idx, h2_idx)
            if ring_factor < 1.0:
                correction_factor *= ring_factor
            
            # Tautomerism
            tautomer_factor = self._get_tautomer_correction(mol, h1_idx, h2_idx)
            if tautomer_factor < 1.0:
                correction_factor *= tautomer_factor
            
            # Rotation amide
            amide_factor = self._get_amide_rotation_correction(mol, h1_idx, h2_idx)
            if amide_factor < 1.0:
                correction_factor *= amide_factor
            
            # Apply and validate the correction.
            J_effective = J_static * correction_factor
            
            # Round to three decimals to avoid numerical artifacts.
            J_effective = round(J_effective, 3)
            
            return J_effective if abs(J_effective) > 0.001 else 0.0
            
        except Exception as e:
            # Return None on error rather than an unreliable value.
            warnings.warn(f"Erreur critique dans apply_dynamics: {e}. Retour None.", UserWarning)
            return None
    
    def _is_rapid_exchange(self, mol: Chem.Mol, h_idx: int) -> bool:
        """Detect rapidly exchanging protons.

        Criteria used here:
        - Ionized OH at pH > pKa.
        - Amine NH₂ at pH > pKa.
        - Deprotonated COOH at pH > pKa.
        - SH is treated as slow and retained."""
        h_atom = mol.GetAtomWithIdx(h_idx)
        if h_atom.GetSymbol() != 'H':
            return False
        
        for neighbor in h_atom.GetNeighbors():
            parent = neighbor
            parent_sym = parent.GetSymbol()
            
            # Oxygen (OH, COOH)
            if parent_sym == 'O':
                factor = SOLVENT_EXCHANGE_FACTORS[self.solvent].get('OH_exchange', 1.0)
                
                # Check for COOH.
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
                
                # pH correction for COOH
                if is_carboxyl:
                    if self.ph > 5.5:  # Au-dessus pKa → COO⁻
                        factor = PH_PROTONATION_EFFECTS['carboxyl']['effect_above_pKa']
                    else:
                        factor = PH_PROTONATION_EFFECTS['carboxyl']['effect_below_pKa']
                
                return factor > 5.0
            
            # === AZOTE (NH₂ amine, pas NH amide) ===
            elif parent_sym == 'N':
                if self._is_amide_nitrogen(mol, parent.GetIdx()):
                    return False  # Amide NH: very slow exchange
                
                factor = SOLVENT_EXCHANGE_FACTORS[self.solvent].get('NH2_exchange', 1.0)
                
                # Correction pH
                if self.ph < 7.0:  # pH acide → NH₃⁺
                    factor *= 0.1
                else:  # pH basique → NH₂
                    factor = PH_PROTONATION_EFFECTS['amino']['effect_above_pKa']
                
                return factor > 2.0
            
            # === SOUFRE (SH) ===
            elif parent_sym == 'S':
                # Treated as slow in metabolomic conditions.
                return False
        
        return False
    
    def _is_amide_nitrogen(self, mol: Chem.Mol, n_idx: int) -> bool:
        """Check whether N belongs to an amide or peptide bond (slow exchange)."""
        n_atom = mol.GetAtomWithIdx(n_idx)
        if n_atom.GetSymbol() != 'N':
            return False
        
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                # Check for C=O.
                for nn in neighbor.GetNeighbors():
                    if nn.GetSymbol() == 'O':
                        bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nn.GetIdx())
                        if bond and bond.GetBondTypeAsDouble() == 2.0:
                            return True
        
        return False
    
    def _get_ring_flipping_correction(self, mol: Chem.Mol, h1_idx: int, h2_idx: int) -> float:
        """Correct ring flipping in saturated five- and six-membered rings.

        Aromatic rings are treated as rigid. This correction applies to
        cyclohexane, pyranose, and furanose.

        Returns:
            A factor from 1.0 (no effect) to 0.3 (strong averaging)."""
        # Find parent atoms using the local cache.
        h1_parent = self._get_hydrogen_parent_cached(mol, h1_idx)
        h2_parent = self._get_hydrogen_parent_cached(mol, h2_idx)
        
        if h1_parent is None or h2_parent is None:
            return 1.0
        
        # Check atom aromaticity.
        atom1 = mol.GetAtomWithIdx(h1_parent)
        atom2 = mol.GetAtomWithIdx(h2_parent)
        
        # Skip ring flipping when either atom is aromatic.
        if atom1.GetIsAromatic() or atom2.GetIsAromatic():
            return 1.0  # Aromatic rings: no correction.
        
        # Trouver cycles communs
        ring_info = mol.GetRingInfo()
        
        for ring in ring_info.AtomRings():
            if len(ring) in [5, 6]:
                if (h1_parent in ring) or (h2_parent in ring):
                    ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                    symbols = [a.GetSymbol() for a in ring_atoms]
                    
                    # Exclude aromatic rings.
                    is_aromatic_ring = any(a.GetIsAromatic() for a in ring_atoms)
                    if is_aromatic_ring:
                        return 1.0  # Aromatic rings are treated as rigid.
                    
                    # Saturated rings containing O (sugars) or only C (cyclohexane)
                    if all(s in ['C', 'O', 'N'] for s in symbols):
                        # Temperature interpolation
                        temp_rounded = round(self.temperature_k)
                        factors = TEMPERATURE_CORRECTION_FACTORS.get('sugar_ring_flipping', {298: 0.6})
                        
                        # Find the nearest temperature.
                        closest_temp = min(factors.keys(), key=lambda t: abs(t - temp_rounded))
                        
                        return factors[closest_temp]
        
        return 1.0
    
    def _get_tautomer_correction(self, mol: Chem.Mol, h1_idx: int, h2_idx: int) -> float:
        """Correct keto-enol and imine-enamine tautomerism.

        Stable aromatic rings (benzene, indole, pyridine) receive no correction.
        Candidate systems include beta-dicarbonyls, nonaromatic imines/enamines,
        and nucleobases. A detected tautomeric system receives a fixed 0.7 factor."""
        h1_parent = self._get_hydrogen_parent_cached(mol, h1_idx)
        h2_parent = self._get_hydrogen_parent_cached(mol, h2_idx)
        
        if h1_parent is None or h2_parent is None:
            return 1.0
        
        atom1 = mol.GetAtomWithIdx(h1_parent)
        atom2 = mol.GetAtomWithIdx(h2_parent)
        
        # Exclude stable aromatic rings from the tautomerism correction.
        # Treat indole, benzene, and pyridine as stable aromatic systems.
        if atom1.GetIsAromatic() and atom2.GetIsAromatic():
            # Check whether both atoms share a stable aromatic ring.
            ring_info = mol.GetRingInfo()
            rings1 = [set(r) for r in ring_info.AtomRings() if h1_parent in r]
            rings2 = [set(r) for r in ring_info.AtomRings() if h2_parent in r]
            
            # A shared aromatic ring indicates a stable system.
            for r1 in rings1:
                for r2 in rings2:
                    if r1 & r2:  # Intersection non vide
                        shared_ring = r1 & r2
                        # Check aromaticity of every atom in the shared ring.
                        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in shared_ring):
                            return 1.0  # Stable aromatic ring: no correction.
        
        # Tautomeric candidates: beta-dicarbonyls and nonaromatic imines
        # Detect C=O...C=C-OH or C=N...C=C-NH motifs.
        if str(atom1.GetHybridization()) == 'SP2' and str(atom2.GetHybridization()) == 'SP2':
            # Exclude stable aromatic systems from tautomeric candidates.
            if not atom1.GetIsAromatic() and not atom2.GetIsAromatic():
                # Check adjacent O and N atoms for tautomerism.
                if any(n.GetSymbol() in ['O', 'N'] for n in atom1.GetNeighbors()) or \
                   any(n.GetSymbol() in ['O', 'N'] for n in atom2.GetNeighbors()):
                    return 0.7  # Tautomerism correction factor
        
        return 1.0
    
    def _get_amide_rotation_correction(self, mol: Chem.Mol, h1_idx: int, h2_idx: int) -> float:
        """Correct slow amide rotation. Return 1.0 for slow rotation or 0.8 for partial averaging at high temperature."""
        for h_idx in [h1_idx, h2_idx]:
            h_atom = mol.GetAtomWithIdx(h_idx)
            if h_atom.GetSymbol() != 'H':
                continue
            
            for neighbor in h_atom.GetNeighbors():
                if neighbor.GetSymbol() == 'N' and self._is_amide_nitrogen(mol, neighbor.GetIdx()):
                    # Amide NH indicates slow rotation.
                    temp_rounded = round(self.temperature_k)
                    factors = TEMPERATURE_CORRECTION_FACTORS.get('amide_rotation', {298: 1.0})
                    
                    closest_temp = min(factors.keys(), key=lambda t: abs(t - temp_rounded))
                    return factors[closest_temp]
        
        return 1.0
    
    def _get_hydrogen_parent_cached(self, mol: Chem.Mol, h_idx: int) -> Optional[int]:
        """Retrieve a hydrogen parent from the local cache, avoiding four repeated lookups in apply_dynamics."""
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
    # Detect rigid systems to preserve their J couplings.
    # =========================================================================
    
    def _should_not_reduce_coupling(self, mol: Chem.Mol, h1_idx: int, h2_idx: int) -> bool:
        """Check whether an H-H pair belongs to a rigid system.

        Covers aromatic rings; strained saturated three- and four-membered rings;
        bridged bicyclic systems; rigid alkenes, alkynes and cumulenes;
        carbonyls conjugated with aromatics; and extended conjugated sp2 systems.

        Returns:
            True to preserve J in a rigid system, False to allow reduction."""
        # Obtenir atomes parents
        h1_parent = self._get_hydrogen_parent_cached(mol, h1_idx)
        h2_parent = self._get_hydrogen_parent_cached(mol, h2_idx)
        
        if h1_parent is None or h2_parent is None:
            return False  # Allow reduction by default.
        
        atom1 = mol.GetAtomWithIdx(h1_parent)
        atom2 = mol.GetAtomWithIdx(h2_parent)
        
        # Check for a shared aromatic ring.
        if self._is_in_same_aromatic_ring(mol, h1_parent, h2_parent, atom1, atom2):
            return True  # Rigid: no reduction.
        
        # Check for a rigid saturated three- or four-membered ring.
        if self._is_in_strained_ring(mol, h1_parent, h2_parent, atom1, atom2):
            return True  # RIGIDE - cycles tendus
        
        # Check for a bridged bicyclic system.
        if self._is_bridged_rigid_system(mol, h1_parent, h2_parent):
            return True  # Rigid bridged bicyclic system
        
        # Check for a rigid alkene or alkyne.
        if self._is_rigid_unsaturated_system(mol, h1_parent, h2_parent, atom1, atom2):
            return True  # RIGIDE - liaison double/triple
        
        # Check conjugation with an aromatic system.
        if self._is_conjugated_with_aromatic(mol, h1_parent, h2_parent):
            return True  # Rigid through aromatic conjugation
        
        return False  # Flexible system: reduction allowed.
    
    def _is_in_same_aromatic_ring(self, mol: Chem.Mol, h1_parent: int, h2_parent: int,
                                   atom1: Chem.Atom, atom2: Chem.Atom) -> bool:
        """Check whether both atoms share an aromatic ring, including benzene, naphthalene, pyridine, pyrrole, furan, thiophene, imidazole, indole, and quinoline."""
        # Both atoms must be aromatic to share an aromatic ring.
        if not (atom1.GetIsAromatic() and atom2.GetIsAromatic()):
            # Handle adjacent atoms when only one is aromatic.
            # For example, H on an sp3 carbon next to benzene does not imply a rigid coupling.
            return False
        
        ring_info = mol.GetRingInfo()
        
        for ring in ring_info.AtomRings():
            if h1_parent in ring and h2_parent in ring:
                # Check that every ring atom is aromatic.
                all_aromatic = all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
                if all_aromatic:
                    return True  # Both atoms share a fully aromatic ring.
        
        # Also check fused rings such as naphthalene and indole.
        # Both aromatic atoms may share a fused ring system.
        rings_containing_1 = [set(r) for r in ring_info.AtomRings() if h1_parent in r]
        rings_containing_2 = [set(r) for r in ring_info.AtomRings() if h2_parent in r]
        
        for r1 in rings_containing_1:
            for r2 in rings_containing_2:
                # Fused rings share atoms.
                shared = r1 & r2
                if shared:
                    # Check whether the entire fused system is aromatic.
                    combined = r1 | r2
                    if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in combined):
                        return True
        
        return False
    
    def _is_in_strained_ring(self, mol: Chem.Mol, h1_parent: int, h2_parent: int,
                              atom1: Chem.Atom, atom2: Chem.Atom) -> bool:
        """Check for a strained saturated three- or four-membered ring, including cyclopropane, cyclobutane, epoxide, aziridine, beta-lactam, and beta-lactone."""
        if not (atom1.IsInRing() and atom2.IsInRing()):
            return False
        
        ring_info = mol.GetRingInfo()
        
        for ring in ring_info.AtomRings():
            ring_size = len(ring)
            
            # Seulement 3-4 membres (cycles tendus)
            if ring_size not in RIGID_RING_SIZES:
                continue
            
            if h1_parent in ring and h2_parent in ring:
                # Require saturation: no double bonds in the ring.
                is_saturated = True
                for i in range(ring_size):
                    bond = mol.GetBondBetweenAtoms(ring[i], ring[(i+1) % ring_size])
                    if bond and bond.GetBondTypeAsDouble() >= 2.0:
                        is_saturated = False
                        break
                
                if is_saturated:
                    return True  # Treat strained saturated rings as rigid.
        
        return False
    
    def _is_bridged_rigid_system(self, mol: Chem.Mol, atom1_idx: int, atom2_idx: int) -> bool:
        """Detect rigid bridged bicyclic systems.

        Candidates include norbornane, adamantane, cubane, prismane, trans-decalin,
        and rigid steroid fragments. Look for atoms in at least two rings
        sharing at least three atoms, or a ring pair with such a bridge."""
        ring_info = mol.GetRingInfo()
        atom_rings = list(ring_info.AtomRings())
        
        if len(atom_rings) < 2:
            return False  # A bridge requires at least two rings.
        
        # Count rings containing each atom.
        def count_rings(atom_idx):
            return sum(1 for ring in atom_rings if atom_idx in ring)
        
        rings_atom1 = count_rings(atom1_idx)
        rings_atom2 = count_rings(atom2_idx)
        
        # An atom in at least two rings may belong to a bridge.
        if rings_atom1 >= MIN_RINGS_FOR_BRIDGED_RIGID or rings_atom2 >= MIN_RINGS_FOR_BRIDGED_RIGID:
            # Check for at least three shared atoms between two rings.
            for i, r1 in enumerate(atom_rings):
                for r2 in atom_rings[i+1:]:
                    shared = set(r1) & set(r2)
                    # A norbornane-like bridge has at least three shared atoms.
                    if len(shared) >= MIN_SHARED_ATOMS_FOR_BRIDGE:
                        # Check whether the target atoms belong to this system.
                        combined = set(r1) | set(r2)
                        if atom1_idx in combined or atom2_idx in combined:
                            return True
        
        # Also check neighboring bridgehead atoms.
        atom1 = mol.GetAtomWithIdx(atom1_idx)
        atom2 = mol.GetAtomWithIdx(atom2_idx)
        
        for neighbor in atom1.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if count_rings(n_idx) >= MIN_RINGS_FOR_BRIDGED_RIGID:
                # Check for a bridge through a neighbor.
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
        
        # Trans-decalin: two fused rings with trans geometry.
        # Approximate detection: two six-membered rings sharing two adjacent atoms.
        six_rings = [set(r) for r in atom_rings if len(r) == 6]
        if len(six_rings) >= 2:
            for i, r1 in enumerate(six_rings):
                for r2 in six_rings[i+1:]:
                    shared = r1 & r2
                    if len(shared) == 2:
                        # Fusion trans → rigide
                        shared_list = list(shared)
                        bond = mol.GetBondBetweenAtoms(shared_list[0], shared_list[1])
                        if bond:  # Adjacent ring-junction atoms
                            if atom1_idx in (r1 | r2) and atom2_idx in (r1 | r2):
                                return True
        
        return False
    
    def _is_rigid_unsaturated_system(self, mol: Chem.Mol, h1_parent: int, h2_parent: int,
                                      atom1: Chem.Atom, atom2: Chem.Atom) -> bool:
        """Check whether the pair is linked through a rigid double or triple bond, including alkenes, alkynes, allene, and cumulenes."""
        # Check for a direct bond between the parent atoms.
        bond = mol.GetBondBetweenAtoms(h1_parent, h2_parent)
        
        if bond:
            bond_order = bond.GetBondTypeAsDouble()
            
            # Liaison double ou triple → RIGIDE
            if bond_order >= 2.0:
                return True
            
            # A single bond between conjugated sp2 atoms may be rigid.
            if str(atom1.GetHybridization()) == 'SP2' and str(atom2.GetHybridization()) == 'SP2':
                # Check for an extended conjugated system.
                if self._is_part_of_conjugated_system(mol, h1_parent, h2_parent):
                    return True
        
        # Without a direct bond, look for an intervening double bond.
        # For example, H-C=C-H gives a cis/trans vicinal coupling.
        try:
            path = Chem.GetShortestPath(mol, h1_parent, h2_parent)
        except:
            return False
        
        if path and len(path) <= 4:  # Short-range coupling
            for i in range(len(path) - 1):
                bond = mol.GetBondBetweenAtoms(path[i], path[i+1])
                if bond and bond.GetBondTypeAsDouble() >= 2.0:
                    # Check for a double bond along the path.
                    # Exclude isolated carbonyls.
                    atom_i = mol.GetAtomWithIdx(path[i])
                    atom_j = mol.GetAtomWithIdx(path[i+1])
                    
                    # Treat C=C as rigid.
                    if atom_i.GetSymbol() == 'C' and atom_j.GetSymbol() == 'C':
                        return True
                    
                    # C≡C (alcyne) → rigide
                    if bond.GetBondTypeAsDouble() == 3.0:
                        return True
        
        return False
    
    def _is_part_of_conjugated_system(self, mol: Chem.Mol, atom1_idx: int, atom2_idx: int) -> bool:
        """Check for an extended conjugated pi system, such as butadiene (C=C-C=C), styrene (Ph-C=C), or an enone (C=C-C=O)."""
        atom1 = mol.GetAtomWithIdx(atom1_idx)
        atom2 = mol.GetAtomWithIdx(atom2_idx)
        
        # Both atoms must be sp2.
        if not (str(atom1.GetHybridization()) == 'SP2' and str(atom2.GetHybridization()) == 'SP2'):
            return False
        
        # Check for alternating single and double bonds around the atoms.
        def has_adjacent_double_bond(atom):
            for bond in atom.GetBonds():
                if bond.GetBondTypeAsDouble() == 2.0:
                    return True
            return False
        
        return has_adjacent_double_bond(atom1) or has_adjacent_double_bond(atom2)
    
    def _is_conjugated_with_aromatic(self, mol: Chem.Mol, atom1_idx: int, atom2_idx: int) -> bool:
        """Check conjugation with an aromatic system.

        Examples include styrene, benzaldehyde, and benzoic acid. These systems
        are treated as rigid because conjugation restricts rotation."""
        atom1 = mol.GetAtomWithIdx(atom1_idx)
        atom2 = mol.GetAtomWithIdx(atom2_idx)
        
        # Check whether either atom is aromatic.
        if atom1.GetIsAromatic() or atom2.GetIsAromatic():
            # Check whether the other atom is sp2 and conjugated.
            non_aromatic = atom2 if atom1.GetIsAromatic() else atom1
            
            if str(non_aromatic.GetHybridization()) == 'SP2':
                return True  # Aromatic conjugation indicates rigidity.
        
        # Check paths within three bonds of an aromatic atom.
        try:
            path = Chem.GetShortestPath(mol, atom1_idx, atom2_idx)
        except:
            return False
        
        if not path:
            return False
        
        # Look for aromatic atoms near the path.
        for atom_idx in path:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetIsAromatic():
                return True
            
            # Check direct neighbors of the path.
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIsAromatic():
                    # An aromatic neighbor and an sp2 path indicate conjugation.
                    if str(atom.GetHybridization()) == 'SP2':
                        return True
        
        return False


# =============================================================================
# Merge or separate diastereotopic signals according to conformational dynamics.
# =============================================================================
# Diastereotopic CH₂ protons near a chiral center can produce:
# - Two distinct patterns in slow exchange (k << πΔν).
# - One averaged pattern in fast exchange (k >> πΔν).
#
# Decision criteria:
# 1. Temperature: T > 320 K favors averaging.
# 2. Rigidity: rings and adjacent chiral centers favor separation.
# 3. Flexibility: acyclic chains may allow averaging.
# 4. Exchange rate relative to the coalescence threshold.
# =============================================================================

# Activation energies for conformational rotation (kJ/mol)
CONFORMATIONAL_ACTIVATION_ENERGIES = {
    'acyclic_CH2_CH2': 12.0,       # Free ethane-like rotation (~3 kcal/mol)
    'acyclic_CH2_CHR': 14.0,       # Slightly more constrained
    'acyclic_CH2_CR2': 16.0,       # More sterically hindered
    'cyclopentane_ring': 25.0,     # Pseudorotation cyclopentane
    'cyclohexane_ring': 45.0,      # Ring flip cyclohexane
    'cyclohexane_ring_chair': 45.0,
    'cycloheptane_ring': 30.0,     # Plus flexible que 6
    'sugar_pyranose': 45.0,        # Comparable cyclohexane
    'sugar_furanose': 25.0,        # Plus flexible (5 membres)
    'piperidine': 46.0,            # N-inversion + ring flip
    'adjacent_to_chiral': 18.0,    # Constrained by chirality
    'benzylic_CH2': 15.0,          # Adjacent aromatique
    'alpha_carbonyl': 20.0,        # Adjacent C=O
}

# Coalescence threshold: k = π × Δν × √2
# Typical diastereotopic Δν: 10-100 Hz, depending on rigidity
TYPICAL_DELTA_NU_HZ = {
    'rigid_cycle': 50.0,           # Cyclohexane, sucres
    'semi_rigid': 30.0,            # Cyclopentane
    'flexible_acyclic': 15.0,      # Acyclic chain
    'adjacent_chiral': 40.0,       # Chiral center one bond away
    'distant_chiral': 20.0,        # Chiral center two or three bonds away
}

# Boltzmann constant
R_GAS = 8.314  # J/(mol·K)


def _calculate_exchange_rate(
    activation_energy_kj: float,
    temperature_k: float,
    prefactor: float = 1e13  # Typical pre-exponential factor (s⁻¹)
) -> float:
    """Calculate the Arrhenius conformational exchange rate.

    k = A × exp(-Ea / RT). activation_energy_kj is in kJ/mol, temperature_k
    in kelvin, and prefactor is A (typically 10¹³ s⁻¹). Return k in s⁻¹."""
    if temperature_k <= 0:
        return 0.0
    
    Ea_J = activation_energy_kj * 1000  # kJ → J
    exponent = -Ea_J / (R_GAS * temperature_k)
    
    # Avoid overflow.
    if exponent < -100:
        return 0.0
    if exponent > 100:
        return prefactor  # Saturation
    
    return prefactor * math.exp(exponent)


def _get_coalescence_threshold(delta_nu_hz: float) -> float:
    """Calculate the coalescence threshold k_c = π × Δν × √2 in s⁻¹.

    delta_nu_hz is the separation of the two signals in Hz.
    k >> k_c gives fast exchange and one averaged pattern;
    k << k_c gives slow exchange and two distinct patterns."""
    return math.pi * abs(delta_nu_hz) * math.sqrt(2)


def _classify_ch2_environment(
    mol: Chem.Mol,
    heavy_idx: int
) -> Tuple[str, float, float]:
    """Classify a CH₂ environment for dynamic corrections.

    mol must have explicit hydrogens; heavy_idx identifies their carbon.
    Return (environment_type, activation_energy_kj, delta_nu_hz)."""
    try:
        atom = mol.GetAtomWithIdx(heavy_idx)
        
        if atom.GetSymbol() != 'C':
            return ('unknown', 15.0, 20.0)
        
        # Check ring membership.
        if atom.IsInRing():
            ring_info = mol.GetRingInfo()
            ring_sizes = [len(r) for r in ring_info.AtomRings() if heavy_idx in r]
            
            if ring_sizes:
                min_ring = min(ring_sizes)
                
                # Cycles rigides (3-4 membres)
                if min_ring in (3, 4):
                    return ('rigid_small_ring', 80.0, 80.0)  # Very rigid
                
                # Cyclopentane (pseudorotation rapide)
                if min_ring == 5:
                    # Check for furanose.
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
                    # Check for pyranose.
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
        
        # Inspect the neighborhood of an acyclic CH₂.
        neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() != 'H']

        # Check for an adjacent chiral center.
        for nb in neighbors:
            if nb.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                return ('adjacent_chiral', 18.0, 40.0)

        # Check for an adjacent quaternary prochiral center.
        # For example, the quaternary C in citrate CH2-C(OH)(COOH)(CH2COOH) is prochiral.
        # The protons remain diastereotopic despite free rotation.
        for nb in neighbors:
            if nb.GetSymbol() == 'C':
                # Count the neighbor's heavy-atom neighbors.
                nb_heavy_neighbors = [n for n in nb.GetNeighbors() if n.GetSymbol() != 'H']
                # Check for a quaternary sp3 carbon: four heavy neighbors and no double bonds.
                if len(nb_heavy_neighbors) == 4:
                    has_double_bond = any(
                        bond.GetBondTypeAsDouble() >= 2
                        for bond in nb.GetBonds()
                    )
                    if not has_double_bond:
                        # Check substituent asymmetry for prochirality.
                        try:
                            ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
                            sub_ranks = [ranks[n.GetIdx()] for n in nb_heavy_neighbors]
                            # At least three distinct types indicate prochirality (OH, COOH, CH2, CH2).
                            if len(set(sub_ranks)) >= 3:
                                return ('adjacent_prochiral_quaternary', 50.0, 80.0)
                            # At least two distinct types, with one different from the CH2 itself.
                            if len(set(sub_ranks)) >= 2:
                                ch2_rank = ranks[heavy_idx]
                                other_ranks = [r for r in sub_ranks if r != ch2_rank]
                                if len(set(other_ranks)) >= 1:
                                    return ('adjacent_prochiral_quaternary', 50.0, 80.0)
                        except Exception:
                            pass

        # Check for an adjacent aromatic atom.
        for nb in neighbors:
            if nb.GetIsAromatic():
                return ('benzylic', 15.0, 25.0)

        # Check for an adjacent carbonyl.
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
    """Decide whether conformational dynamics merge two diastereotopic signals.

    Fast exchange (k_exchange >> k_coalescence) favors one pattern;
    slow exchange favors two patterns.

    Args:
        mol: RDKit molecule with explicit hydrogens.
        heavy_idx: Carbon bearing the two protons.
        h_a_idx, h_b_idx: Proton indices.
        temperature_k: Temperature in kelvin (default 298.15 = 25°C).
        spectrometer_mhz: Spectrometer frequency, used to scale Δν in Hz.

    Returns:
        (should_merge, reason, details), where details includes exchange and
        coalescence rates and their ratio."""
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
        # Classify the CH₂ environment.
        env_type, Ea_kj, delta_nu_base = _classify_ch2_environment(mol, heavy_idx)
        details['environment'] = env_type
        details['activation_energy_kj'] = Ea_kj
        
        # Scale Δν to the spectrometer frequency.
        # Δν (Hz) = Δδ (ppm) × freq (MHz)
        # Typical diastereotopic Δδ: 0.01-0.1 ppm
        delta_nu_hz = delta_nu_base * (spectrometer_mhz / 700.0)
        details['delta_nu_hz'] = delta_nu_hz
        
        # Calculate the conformational exchange rate.
        k_exchange = _calculate_exchange_rate(Ea_kj, temperature_k)
        details['k_exchange'] = k_exchange
        
        # Calculate the coalescence threshold.
        k_coalescence = _get_coalescence_threshold(delta_nu_hz)
        details['k_coalescence'] = k_coalescence
        
        # Calculate the exchange ratio.
        if k_coalescence > 0:
            exchange_ratio = k_exchange / k_coalescence
        else:
            exchange_ratio = float('inf')
        details['exchange_ratio'] = exchange_ratio
        
        # =====================================================================
        # Decision rules
        # =====================================================================
        
        # Always separate protons adjacent to a quaternary prochiral center.
        # For example, citrate CH2-C(OH)(COOH)(R) has a quaternary prochiral carbon.
        # The protons lie on different faces of the prochiral center.
        # They remain diastereotopic despite free rotation.
        if env_type == 'adjacent_prochiral_quaternary':
            return (False, "Adjacent à centre prochiral quaternaire → DEUX patterns (diastéréotopes)", details)

        # Always separate signals in rigid systems.
        # Small rings (3-4), bridged bicyclic systems, and directly adjacent chiral centers
        if env_type in ('rigid_small_ring',):
            return (False, "Cycle tendu rigide (3-4 membres) → DEUX patterns", details)
        
        # Separate six-membered ring signals (cyclohexane, pyranose) below 350 K.
        # Ring flipping is treated as slow at room temperature.
        if env_type in ('cyclohexane', 'sugar_pyranose'):
            if temperature_k < 350:
                return (False, f"Cycle 6 membres rigide à {temperature_k:.0f}K → DEUX patterns", details)
            else:
                return (True, f"Cycle 6 membres avec ring flip rapide à {temperature_k:.0f}K → UN pattern", details)
        
        # Allow cyclopentane/furanose merging from 298 K.
        # Pseudorotation rapide (~10⁹ s⁻¹)
        if env_type in ('cyclopentane', 'sugar_furanose'):
            if exchange_ratio > 10:
                return (True, "Cycle 5 membres avec pseudorotation rapide → UN pattern", details)
            elif temperature_k > 310:
                return (True, f"Cycle 5 membres à T={temperature_k:.0f}K → UN pattern", details)
            else:
                return (False, "Cycle 5 membres mais Δν élevé → DEUX patterns", details)
        
        # Separate signals next to chiral centers except at very high temperatures.
        if env_type == 'adjacent_chiral':
            if temperature_k < 340:
                return (False, "Adjacent à centre chiral → DEUX patterns", details)
            else:
                return (True, f"Adjacent à chiral mais T={temperature_k:.0f}K → moyennage partiel", details)
        
        # Benzylic CH₂ behavior depends on asymmetry.
        if env_type == 'benzylic':
            if exchange_ratio > 50:
                return (True, "CH₂ benzylique avec rotation rapide → UN pattern", details)
            else:
                return (False, "CH₂ benzylique avec environnement asymétrique → DEUX patterns", details)
        
        # Merge flexible acyclic CH₂ signals under free rotation.
        if env_type == 'acyclic_flexible':
            if exchange_ratio > 100:
                return (True, "CH₂ acyclique avec rotation libre → UN pattern moyenné", details)
            elif temperature_k > 320:
                return (True, f"CH₂ acyclique à T={temperature_k:.0f}K → moyennage", details)
            else:
                # Keep acyclic signals separate near a chiral center.
                # Check for a chiral center one to three bonds away.
                return (False, "CH₂ acyclique mais environnement asymétrique → DEUX patterns", details)
        
        # Large rings (at least seven members) are more flexible.
        if env_type == 'large_ring':
            if exchange_ratio > 20:
                return (True, "Grand cycle flexible → UN pattern", details)
            else:
                return (False, "Grand cycle mais contraint → DEUX patterns", details)
        
        # Alpha-carbonyl behavior depends on temperature.
        if env_type == 'alpha_carbonyl':
            if temperature_k > 330 or exchange_ratio > 30:
                return (True, "CH₂ α-carbonyle avec rotation → UN pattern", details)
            else:
                return (False, "CH₂ α-carbonyle contraint → DEUX patterns", details)
        
        # Default: compare exchange and coalescence rates.
        if exchange_ratio > 10:
            return (True, f"Échange rapide (ratio={exchange_ratio:.1f}) → UN pattern", details)
        elif exchange_ratio < 0.1:
            return (False, f"Échange lent (ratio={exchange_ratio:.1f}) → DEUX patterns", details)
        else:
            # Intermediate exchange broadens lines.
            # Keep signals separate and add an explanatory note.
            return (False, f"Échange intermédiaire (ratio={exchange_ratio:.1f}) → DEUX patterns élargis", details)
        
    except Exception as e:
        warnings.warn(f"Erreur should_merge_diastereotopic_CH2: {e}", UserWarning)
        # Default to separate signals.
        return (False, f"Erreur: {e} → défaut DEUX patterns", details)


def get_averaged_diastereotopic_shift(
    shift_a_ppm: float,
    shift_b_ppm: float,
    population_a: float = 0.5,
    population_b: float = 0.5
) -> float:
    """Calculate the population-weighted shift under fast exchange.

    δ_obs = P_a × δ_a + P_b × δ_b. Individual shifts are in ppm;
    conformer populations must sum to one. Return the averaged shift in ppm."""
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
    """Calculate population-weighted couplings under fast exchange.

    For each coupling, J_obs = P_a × J_a + P_b × J_b.
    J_list_a and J_list_b contain the proton couplings; population_a and
    population_b are the conformer populations. Return the averaged J list."""
    total_pop = population_a + population_b
    if total_pop <= 0:
        total_pop = 1.0
    
    # Align the list lengths.
    max_len = max(len(J_list_a), len(J_list_b))
    J_a = list(J_list_a) + [0.0] * (max_len - len(J_list_a))
    J_b = list(J_list_b) + [0.0] * (max_len - len(J_list_b))
    
    averaged = []
    for ja, jb in zip(J_a, J_b):
        j_avg = (population_a * ja + population_b * jb) / total_pop
        if abs(j_avg) > 0.5:  # Detection threshold
            averaged.append(round(j_avg, 1))
    
    return averaged


# =============================================================================
# Public API used by sim_mult_rulebased.py
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
    """Apply dynamic corrections to a static J coupling.

    Example:
        J_static = coupling_extensions.get_coupling(...)
        J_final = apply_dynamics_to_coupling(mol, h1_idx, h2_idx, J_static, ph=7.4)

    Args:
        mol: RDKit molecule with explicit hydrogens.
        h1_idx, h2_idx: Proton indices.
        J_static: Static coupling in Hz.
        ph: Solution pH in [1.0, 14.0].
        solvent: NMR solvent (D2O, CD3OD, CDCl3, DMSO-d6).
        temperature_k: Temperature in kelvin, in [273.15, 373.15].
        use_cache: Enable the global cache for production; disable for tests.

    Returns:
        Effective J in Hz, or None for invalid input or a critical error."""
    if J_static is None or abs(J_static) < 0.001:
        return None
    
    if not use_cache:
        # Uncached calculation for debugging and tests
        processor = DynamicsPostProcessor(ph, solvent, temperature_k)
        return processor.apply_dynamics(mol, h1_idx, h2_idx, J_static)
    
    # Cached calculation for production
    try:
        mol_smiles = Chem.MolToSmiles(mol, canonical=True)
        cache_key = f"{mol_smiles}_{h1_idx}_{h2_idx}_{J_static:.3f}_{ph:.2f}_{solvent}_{temperature_k:.2f}"
    except Exception as e:
        warnings.warn(f"Erreur génération clé cache: {e}. Calcul sans cache.", UserWarning)
        processor = DynamicsPostProcessor(ph, solvent, temperature_k)
        return processor.apply_dynamics(mol, h1_idx, h2_idx, J_static)
    
    # Use a dictionary cache to avoid pickle overhead.
    if not hasattr(apply_dynamics_to_coupling, '_cache'):
        apply_dynamics_to_coupling._cache = {}
    
    if cache_key in apply_dynamics_to_coupling._cache:
        return apply_dynamics_to_coupling._cache[cache_key]
    
    # Calculate and cache the result.
    processor = DynamicsPostProcessor(ph, solvent, temperature_k)
    J_final = processor.apply_dynamics(mol, h1_idx, h2_idx, J_static)
    
    # Limiter taille cache
    if len(apply_dynamics_to_coupling._cache) < 50000:
        apply_dynamics_to_coupling._cache[cache_key] = J_final
    
    return J_final


# =============================================================================
# Unit tests; intended for tests/test_dynamics.py in production
# =============================================================================

def _run_tests():
    """Run validation checks outside production."""
    print("=" * 70)
    print("DYNAMICS POST-PROCESSOR - TESTS VALIDATION")
    print("=" * 70)
    
    # Test fast exchange in ethanol.
    print("\nTest 1: Ethanol (exchangeable OH, pH=7.4, D2O)")
    mol = Chem.MolFromSmiles("CCO")
    mol = Chem.AddHs(mol)
    h_atoms = [i for i, atom in enumerate(mol.GetAtoms()) if atom.GetSymbol() == 'H']
    
    if len(h_atoms) >= 2:
        # The OH-CH2 coupling should be suppressed.
        J_test = apply_dynamics_to_coupling(
            mol, h_atoms[4], h_atoms[2], J_static=6.0, ph=7.4, solvent="D2O", use_cache=False
        )
        print(f"  J statique: 6.0 Hz → J effectif: {J_test} Hz (attendu: 0.0)")
        assert J_test == 0.0, "OH should exchange rapidly"
        print("  PASS")
    
    # Test 2: Ring flipping (glucose)
    print("\nTest 2: Glucose β-D (ring flipping, 25°C)")
    mol = Chem.MolFromSmiles("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O")
    mol = Chem.AddHs(mol)
    h_atoms = [i for i, atom in enumerate(mol.GetAtoms()) if atom.GetSymbol() == 'H']
    
    if len(h_atoms) >= 2:
        J_test = apply_dynamics_to_coupling(
            mol, h_atoms[0], h_atoms[1], J_static=8.0, ph=7.4, solvent="D2O", use_cache=False
        )
        print(f"  J statique: 8.0 Hz → J effectif: {J_test} Hz (attendu: < 8.0)")
        assert J_test < 8.0 and J_test > 0, "Ring flipping should reduce J"
        print("  PASS")
    
    # Test 3: Cache fonctionnel
    print("\nTest 3: Global cache (50k entries)")
    mol = Chem.MolFromSmiles("CC(N)C(O)=O")  # Alanine
    mol = Chem.AddHs(mol)
    h_atoms = [i for i, atom in enumerate(mol.GetAtoms()) if atom.GetSymbol() == 'H']
    
    if len(h_atoms) >= 2:
        J1 = apply_dynamics_to_coupling(mol, h_atoms[0], h_atoms[1], J_static=7.0, ph=7.4)
        J2 = apply_dynamics_to_coupling(mol, h_atoms[0], h_atoms[1], J_static=7.0, ph=7.4)
        assert J1 == J2, "Cache returned a different value"
        print("  PASS: cache lookup")
        
        # Stats cache
        if hasattr(apply_dynamics_to_coupling, '_cache'):
            print(f"  Cache size: {len(apply_dynamics_to_coupling._cache)}")
    
    print("\n" + "=" * 70)
    print("All tests passed.")
    print("=" * 70)


def clear_dynamics_cache():
    """Clear the global cache for tests and benchmarks."""
    if hasattr(apply_dynamics_to_coupling, '_cache'):
        apply_dynamics_to_coupling._cache.clear()
        print("Dynamics cache cleared")


# =============================================================================
# Run tests only when invoked as a script.
# =============================================================================

if __name__ == "__main__":
    # Disable tests in production.
    import os
    if os.environ.get('RUN_DYNAMICS_TESTS', '0') == '1':
        _run_tests()
    else:
        print("=" * 70)
        print("DYNAMICS POST-PROCESSOR - MODULE PRÊT POUR PRODUCTION")
        print("=" * 70)
        print("Run tests with: RUN_DYNAMICS_TESTS=1 python dynamics.py")
        print("\nUsage:")
        print("  J_static = coupling_extensions.get_coupling(...)")
        print("  J_final = apply_dynamics_to_coupling(mol, h1_idx, h2_idx, J_static, ph=7.4)")
        print("\nOptimizations:")
        print("  - Cache 50k entrées")
        print("  - Gestion pH/solvant/température")
        print("  - Post-traitement seulement (pas de calcul J)")
        print("=" * 70)
