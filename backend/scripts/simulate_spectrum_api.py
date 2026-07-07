#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
simulate_spectrum_api.py

Script CLI pour simuler un spectre RMN 1H complet avec formes de raie.
Utilise DIRECTEMENT les données de multiplicité (JSON) générées par le simulateur.

DEUX MODES :
1. --json_path : Utilise un fichier JSON de multiplicité existant (recommandé, plus précis)
2. --smiles : Calcule les multiplicités à la volée (fallback)

Usage:
    # Mode 1: Utiliser le JSON existant (RECOMMANDÉ)
    python simulate_spectrum_api.py --json_path "output/mol_multiplicity.json" --freq_mhz 700
    
    # Mode 2: Calculer depuis SMILES (fallback)
    python simulate_spectrum_api.py --smiles "CCO" --freq_mhz 700
    
Output (JSON sur stdout):
    {
      "mol_id": "mol",
      "smiles": "CCO",
      "freq_mhz": 700,
      "ppm": [0.0, 0.001, ...],
      "intensity": [0.0, 0.1, ...],
      "peaks": [{"ppm": 1.2, "intensity": 1.0, "label": "t", "J_hz": [7.0]}, ...],
      "n_points": 120000,
      "status": "ok"
    }
"""

import argparse
import json
import sys
import os
import traceback

# Add parent directory to path for imports
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)


def main():
    parser = argparse.ArgumentParser(description="Simulate 1H NMR spectrum from multiplicity JSON")
    
    # Input options (mutually exclusive in practice)
    parser.add_argument("--json_path", help="Path to multiplicity JSON file (preferred)")
    parser.add_argument("--json_data", help="Inline JSON data for a single molecule")
    parser.add_argument("--smiles", help="SMILES string (fallback, will compute multiplicities)")
    
    # Common options
    parser.add_argument("--mol_id", default="mol", help="Molecule ID")
    parser.add_argument("--freq_mhz", type=float, default=700.0, help="Spectrometer frequency in MHz")
    parser.add_argument("--solvent", default="CDCl3", help="Solvent")
    parser.add_argument("--line_shape", default="voigt", choices=["lorentzian", "gaussian", "voigt"], help="Line shape")
    parser.add_argument("--digital_resolution_hz", type=float, default=0.05, help="Digital resolution in Hz")
    parser.add_argument("--ppm_min", type=float, default=-0.5, help="Minimum ppm")
    parser.add_argument("--ppm_max", type=float, default=12.0, help="Maximum ppm")
    parser.add_argument("--line_width_hz", type=float, default=1.2, help="Default line width in Hz")
    # Paramètres dynamics post-processeur
    parser.add_argument("--enable_dynamics", action="store_true", help="Enable dynamics post-processor")
    parser.add_argument("--temperature_k", type=float, default=298.0, help="Temperature in Kelvin for dynamics")
    parser.add_argument("--ph", type=float, default=7.4, help="pH for dynamics post-processor (1.0-14.0)")
    
    # Mode batch pour traiter toutes les molécules d'un JSON consolidé
    parser.add_argument("--all", action="store_true", help="Process ALL molecules in a consolidated JSON (batch mode)")
    
    args = parser.parse_args()
    
    try:
        import numpy as np
        
        # Determine the source of multiplicity data
        mult_data = None
        
        # Priority 1: JSON file path
        if args.json_path and os.path.exists(args.json_path):
            with open(args.json_path, "r", encoding="utf-8") as f:
                mult_data = json.load(f)
            # 🔍 DIAGNOSTIC: Vérifier la présence de 'partners' dans le JSON lu du fichier
            chain_groups = [g for g in mult_data.get('groups', []) if 0.8 <= g.get('H_shifts_ppm', [0])[0] <= 1.6]
            total_partners = sum(len(g.get('partners', [])) for g in chain_groups)
            print(f"🔍 [DIAGNOSTIC] JSON lu depuis {args.json_path}: {len(chain_groups)} groupes chaîne, {total_partners} partners totaux", file=sys.stderr)

        # Priority 2: Inline JSON data
        elif args.json_data:
            mult_data = json.loads(args.json_data)
            # 🔍 DIAGNOSTIC: Vérifier la présence de 'partners' dans le JSON reçu
            chain_groups = [g for g in mult_data.get('groups', []) if 0.8 <= g.get('H_shifts_ppm', [0])[0] <= 1.6]
            total_partners = sum(len(g.get('partners', [])) for g in chain_groups)
            print(f"🔍 [DIAGNOSTIC] JSON reçu via --json_data: {len(chain_groups)} groupes chaîne, {total_partners} partners totaux", file=sys.stderr)
            
        # Priority 3: Compute from SMILES
        elif args.smiles:
            mult_data = compute_multiplicities_from_smiles(
                args.smiles, 
                args.mol_id,
                args.freq_mhz,
                args.solvent
            )
        else:
            output_error("Either --json_path, --json_data, or --smiles must be provided")
            return
        
        if mult_data is None:
            output_error("Failed to get multiplicity data")
            return
        
        # ============================================
        # MODE BATCH: Traiter toutes les molécules
        # ============================================
        if getattr(args, 'all', False) and "molecules" in mult_data:
            all_spectra = []
            molecules = mult_data.get("molecules", [])
            
            for mol in molecules:
                try:
                    mol_id = mol.get("mol_id", "unknown")
                    result = generate_spectrum_from_multiplicity_json(
                        mol,
                        freq_mhz=args.freq_mhz,
                        line_shape=args.line_shape,
                        ppm_min=args.ppm_min,
                        ppm_max=args.ppm_max,
                        digital_resolution_hz=args.digital_resolution_hz,
                        line_width_hz=args.line_width_hz,
                        enable_dynamics=args.enable_dynamics,
                        temperature_k=args.temperature_k,
                        ph=args.ph,
                        solvent=args.solvent,
                    )
                    result["mol_id"] = mol_id
                    result["smiles"] = mol.get("smiles", "")
                    all_spectra.append(result)
                except Exception as e:
                    # En cas d'erreur sur une molécule, l'ajouter avec statut error
                    all_spectra.append({
                        "status": "error",
                        "mol_id": mol.get("mol_id", "unknown"),
                        "smiles": mol.get("smiles", ""),
                        "error": str(e)
                    })
            
            # Output all spectra as JSON array
            batch_result = {
                "status": "ok",
                "mode": "batch",
                "total_molecules": len(molecules),
                "successful": sum(1 for s in all_spectra if s.get("status") == "ok"),
                "failed": sum(1 for s in all_spectra if s.get("status") == "error"),
                "spectra": all_spectra
            }
            print(json.dumps(batch_result, ensure_ascii=False))
            return
        
        # ============================================
        # MODE SINGLE: Une seule molécule
        # ============================================
        # Handle consolidated JSON (multiple molecules) vs single molecule
        if "molecules" in mult_data:
            # Consolidated JSON - find the right molecule
            mol_id = args.mol_id
            molecule_data = None
            for mol in mult_data["molecules"]:
                if mol.get("mol_id") == mol_id:
                    molecule_data = mol
                    break
            if molecule_data is None and mult_data["molecules"]:
                molecule_data = mult_data["molecules"][0]  # Use first if not found
        else:
            # Single molecule JSON
            molecule_data = mult_data
        
        if molecule_data is None:
            output_error("No molecule data found in JSON")
            return
        
        # Generate spectrum from multiplicity data
        result = generate_spectrum_from_multiplicity_json(
            molecule_data,
            freq_mhz=args.freq_mhz,
            line_shape=args.line_shape,
            ppm_min=args.ppm_min,
            ppm_max=args.ppm_max,
            digital_resolution_hz=args.digital_resolution_hz,
            line_width_hz=args.line_width_hz,
            enable_dynamics=args.enable_dynamics,
            temperature_k=args.temperature_k,
            ph=args.ph,
            solvent=args.solvent,
        )
        
        # Output JSON to stdout
        print(json.dumps(result, ensure_ascii=False))
        
    except Exception as e:
        output_error(str(e), traceback.format_exc())


def output_error(message: str, tb: str = ""):
    """Output error as JSON"""
    result = {
        "status": "error",
        "error": message,
        "traceback": tb[:2000] if tb else "",
    }
    print(json.dumps(result, ensure_ascii=False))
    sys.exit(1)


def compute_multiplicities_from_smiles(smiles: str, mol_id: str, freq_mhz: float, solvent: str) -> dict:
    """
    Compute multiplicities from SMILES using the full simulation pipeline.
    This is the fallback when no pre-computed JSON is available.
    """
    try:
        import sim_mult_rulebased as sim_mod
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        molH = Chem.AddHs(mol)
        AllChem.EmbedMolecule(molH, AllChem.ETKDG())
        
        # Create config
        if hasattr(sim_mod, 'SimulationConfig'):
            cfg = sim_mod.SimulationConfig(
                freq_mhz=freq_mhz,
                solvent=solvent,
                line_width_hz=1.0,
            )
        else:
            cfg = None
        
        # Try to use compute_all_H if available
        if hasattr(sim_mod, 'compute_all_H'):
            result = sim_mod.compute_all_H(molH, cfg)
            if result and hasattr(result, 'per_proton'):
                # Convert to our expected format
                groups = []
                for h_idx, data in result.per_proton.items():
                    groups.append({
                        "group_id": len(groups),
                        "H_indices": [h_idx],
                        "H_shifts_ppm": [data.get("shift_ppm", 0.0)],
                        "label": data.get("label", "s"),
                        "J_couplings": data.get("J_couplings", []),
                        "partners": data.get("partners", []),
                    })
                return {
                    "mol_id": mol_id,
                    "smiles": smiles,
                    "freq_MHz": freq_mhz,
                    "solvent": solvent,
                    "groups": groups,
                }
        
        # Fallback: estimate shifts from structure
        groups = []
        for atom in molH.GetAtoms():
            if atom.GetAtomicNum() != 1:
                continue
            h_idx = atom.GetIdx()
            shift = _estimate_shift(molH, h_idx)
            groups.append({
                "group_id": len(groups),
                "H_indices": [h_idx],
                "H_shifts_ppm": [shift],
                "label": "s",  # Default to singlet
                "J_couplings": [],
                "partners": [],
            })
        
        return {
            "mol_id": mol_id,
            "smiles": smiles,
            "freq_MHz": freq_mhz,
            "solvent": solvent,
            "groups": groups,
        }
        
    except Exception as e:
        print(f"[WARN] compute_multiplicities_from_smiles error: {e}", file=sys.stderr)
        return None


def merge_groups_by_heavy_idx(groups: list) -> list:
    """
    Les protons diastéréotopiques/vinyliques sur le même carbone sont 
    chimiquement non-équivalents et doivent rester SÉPARÉS.
    
    Chaque proton garde:
    - Son propre shift (prédit par ML)
    - Ses propres J_couplings individuels
    - Son propre pattern de multiplicité
    
    Cette fonction ne fait PAS de fusion - elle retourne les groupes tels quels.
    Les protons vinyliques =CH₂ apparaîtront comme deux dd distincts.
    
    Returns:
        Liste de groupes inchangée
    """
    # Pas de fusion - retourner les groupes tels quels
    # Les protons diastéréotopiques sont non-équivalents et doivent rester séparés
    return groups
    return result_groups


def _superpose_groups(group_list: list, heavy_idx: int) -> list:
    """
    Superpose les groupes au même shift moyen, en gardant leurs J individuels.
    
    Chaque proton garde son pattern de couplage (label, J_couplings, partners)
    mais tous sont déplacés vers le shift moyen.
    
    Returns:
        Liste de groupes modifiés (même taille que l'entrée)
    """
    import sys
    
    # Calculer le shift moyen de tous les protons
    all_shifts = []
    for g in group_list:
        all_shifts.extend(g.get("H_shifts_ppm", []))
    
    if not all_shifts:
        return group_list
    
    avg_shift = sum(all_shifts) / len(all_shifts)
    
    # Log
    original_shifts = [g.get("H_shifts_ppm", [0])[0] for g in group_list]
    print(f"[INFO] Superposition heavy_idx={heavy_idx}: "
          f"shifts {[f'{s:.3f}' for s in original_shifts]} -> {avg_shift:.4f} ppm", 
          file=sys.stderr)
    
    # Modifier chaque groupe pour utiliser le shift moyen
    result = []
    for g in group_list:
        # Copier le groupe pour ne pas modifier l'original
        new_g = g.copy()
        
        # Remplacer tous les shifts par le shift moyen
        n_h = len(g.get("H_indices", [1]))
        new_g["H_shifts_ppm"] = [avg_shift] * n_h
        
        # Marquer comme superposé
        new_g["note"] = "superposed_equivalent"
        
        # Garder label, J_couplings, partners INCHANGÉS
        # (déjà copiés par g.copy())
        
        result.append(new_g)
    
    return result

def _merge_group_list(group_list: list, heavy_idx: int) -> dict:
    """
    Fusionne une liste de groupes en un seul groupe.
    
    - Moyenne des shifts
    - Union des H_indices
    - Combinaison des J_couplings (moyenne si différents pour même partenaire)
    - Label recalculé selon le nombre de protons et couplings
    """
    # Collecter tous les H_indices et shifts
    all_h_indices = []
    all_shifts = []
    all_j_couplings = []
    all_partners = []
    
    for g in group_list:
        all_h_indices.extend(g.get("H_indices", []))
        all_shifts.extend(g.get("H_shifts_ppm", []))
        all_j_couplings.extend(g.get("J_couplings", []))
        all_partners.extend(g.get("partners", []))
    
    # Calculer le shift moyen
    avg_shift = sum(all_shifts) / len(all_shifts) if all_shifts else 0.0
    n_protons = len(all_h_indices)
    
    # =========================================================================
    # FUSION DES J_COUPLINGS
    # =========================================================================
    # Pour les protons fusionnés (ex: vinyl =CH₂), les couplages vers le même
    # partenaire doivent être moyennés car ils apparaissent comme un seul signal.
    # Ex: J_trans=16 Hz et J_cis=10 Hz vers H33 → J_moyen=13 Hz
    # 
    # IMPORTANT: Les couplages ENTRE les protons fusionnés eux-mêmes (geminal)
    # doivent être EXCLUS car on ne couple pas avec soi-même après fusion.
    # =========================================================================
    
    # 1. Construire un dictionnaire partner_h_index -> liste de J
    #    EXCLURE les partenaires qui font partie des protons fusionnés
    partner_j_values = {}
    for p in all_partners:
        ph_idx = p.get("partner_h_index")
        j_hz = p.get("J_hz", 0)
        if ph_idx is not None and abs(j_hz) > 0.1:
            # EXCLURE les couplages entre protons fusionnés (geminal interne)
            if ph_idx in all_h_indices:
                continue  # Ce partenaire est un des protons fusionnés, ignorer
            if ph_idx not in partner_j_values:
                partner_j_values[ph_idx] = []
            partner_j_values[ph_idx].append(j_hz)
    
    # 2. Moyenner les J vers chaque partenaire externe
    merged_j = []
    for ph_idx, j_list in partner_j_values.items():
        avg_j = sum(j_list) / len(j_list)
        # n_equiv = 1 car c'est un seul partenaire
        merged_j.append({"J_hz": round(avg_j, 2), "n_equiv": 1})
    
    # Log de debug
    print(f"[DEBUG] Fusion heavy_idx={heavy_idx}: H_indices={all_h_indices}, "
          f"J_fusionnés={merged_j}", file=__import__('sys').stderr)
    
    # 3. Si pas de partners, utiliser les J_couplings directement (dédupliqués)
    if not merged_j and all_j_couplings:
        j_values = {}
        for jc in all_j_couplings:
            j_hz = jc.get("J_hz", 0)
            if abs(j_hz) > 0.1:
                key = round(j_hz * 2) / 2
                if key not in j_values:
                    j_values[key] = {"J_hz": j_hz, "n_equiv": jc.get("n_equiv", 1), "count": 1}
                else:
                    old = j_values[key]
                    old["J_hz"] = (old["J_hz"] * old["count"] + j_hz) / (old["count"] + 1)
                    old["count"] += 1
        merged_j = [{"J_hz": v["J_hz"], "n_equiv": v["n_equiv"]} for v in j_values.values()]
    
    # Dédupliquer les partners (par partner_h_index) avec J moyenné
    # EXCLURE les partenaires internes (protons fusionnés)
    partner_dict = {}
    for p in all_partners:
        ph_idx = p.get("partner_h_index")
        if ph_idx is not None and ph_idx not in all_h_indices:  # Exclure internes
            if ph_idx not in partner_dict:
                partner_dict[ph_idx] = p.copy()
            else:
                # Moyenner les J si présent
                old_j = partner_dict[ph_idx].get("J_hz", 0)
                new_j = p.get("J_hz", 0)
                if old_j and new_j:
                    partner_dict[ph_idx]["J_hz"] = (old_j + new_j) / 2
    
    merged_partners = list(partner_dict.values())
    
    # Déterminer le label basé sur les J restants
    label = _determine_label_from_j(merged_j)
    
    # Prendre les métadonnées du premier groupe
    first = group_list[0]
    
    return {
        "group_id": 0,  # Sera réindexé
        "heavy_idx": heavy_idx,
        "heavy_symbol": first.get("heavy_symbol", "C"),
        "H_indices": all_h_indices,
        "H_shifts_ppm": [avg_shift] * n_protons,
        "label": label,
        "label_index": first.get("label_index", 0),
        "note": "merged_equivalent",
        "second_order_type": first.get("second_order_type", "first_order"),
        "system_type": first.get("system_type"),
        "virtual_coupling": first.get("virtual_coupling", False),
        "J_couplings": merged_j,
        "partners": merged_partners,
        "satellites_13C_hz": first.get("satellites_13C_hz", []),
    }


def _determine_label_from_j(j_couplings: list) -> str:
    """
    Détermine le label de multiplicité à partir des J_couplings.
    """
    if not j_couplings:
        return "s"
    
    # Compter le nombre total de couplages
    total_couplings = sum(jc.get("n_equiv", 1) for jc in j_couplings)
    
    # Labels standards
    labels = {
        1: "d",
        2: "t",
        3: "q",
        4: "quint",
        5: "sext",
        6: "sept",
    }
    
    if len(j_couplings) == 1:
        return labels.get(total_couplings, "m")
    elif len(j_couplings) == 2:
        # dd, dt, dq, etc.
        n1 = j_couplings[0].get("n_equiv", 1)
        n2 = j_couplings[1].get("n_equiv", 1)
        l1 = labels.get(n1, "m")
        l2 = labels.get(n2, "m")
        return f"{l1}{l2[0] if l2 else ''}"
    else:
        return "m"


def find_connected_group_fragments(groups: list) -> list:
    """
    Trouve les fragments de groupes connectés via le champ 'partners'.

    Utilise Union-Find pour identifier les systèmes de spins isolés.
    Deux groupes sont connectés s'ils partagent un couplage J (via partners).

    Args:
        groups: Liste de groupes avec leurs H_indices et partners

    Returns:
        Liste de fragments, chaque fragment étant une liste de groupes
        Ex: [[g0, g1, g2], [g3, g4], [g5]] = 3 systèmes isolés
    """
    if not groups:
        return []

    # Créer un mapping H_index -> group_id
    h_to_group = {}
    group_by_id = {}
    for g in groups:
        gid = g.get("group_id", id(g))
        group_by_id[gid] = g
        for h_idx in g.get("H_indices", []):
            h_to_group[h_idx] = gid

    # Union-Find
    parent = {g.get("group_id", id(g)): g.get("group_id", id(g)) for g in groups}

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    # Connecter les groupes via partners
    for g in groups:
        gid = g.get("group_id", id(g))
        partners = g.get("partners", [])
        for p in partners:
            partner_h_idx = p.get("partner_h_index")
            if partner_h_idx is not None and partner_h_idx in h_to_group:
                partner_gid = h_to_group[partner_h_idx]
                union(gid, partner_gid)

    # Regrouper par composante connexe
    components = {}
    for g in groups:
        gid = g.get("group_id", id(g))
        root = find(gid)
        if root not in components:
            components[root] = []
        components[root].append(g)

    fragments = list(components.values())

    # Log si fragmentation utile
    if len(fragments) > 1:
        total_spins = sum(len(g.get("H_indices", [])) for g in groups)
        print(f"[INFO] GROUP FRAGMENTATION: {len(groups)} groups → {len(fragments)} isolated systems", file=sys.stderr)
        for i, frag in enumerate(fragments):
            n_spins = sum(len(g.get("H_indices", [])) for g in frag)
            gids = [g.get("group_id") for g in frag]
            print(f"[INFO]   Fragment {i+1}: {len(frag)} groups, {n_spins} spins, IDs={gids}", file=sys.stderr)

    return fragments


def _estimate_shift(mol, h_idx: int) -> float:
    """Estimate chemical shift for a hydrogen atom."""
    from rdkit import Chem
    
    atom = mol.GetAtomWithIdx(h_idx)
    neighbors = atom.GetNeighbors()
    
    if not neighbors:
        return 2.0
    
    heavy = neighbors[0]
    symbol = heavy.GetSymbol()
    hyb = heavy.GetHybridization()
    
    if symbol == 'C':
        if hyb == Chem.rdchem.HybridizationType.SP3:
            shift = 1.5
            for n in heavy.GetNeighbors():
                if n.GetAtomicNum() in [7, 8, 9, 17, 35]:
                    shift += 1.5
            return min(shift, 5.0)
        elif hyb == Chem.rdchem.HybridizationType.SP2:
            if heavy.GetIsAromatic():
                return 7.2
            return 5.5
        elif hyb == Chem.rdchem.HybridizationType.SP:
            return 2.5
    elif symbol == 'N':
        return 5.0
    elif symbol == 'O':
        return 3.5
    
    return 2.0


def merge_strongly_coupled_groups(groups: list, freq_mhz: float, max_spins: int = 8, threshold: float = 3.0) -> list:
    """
    Pré-traitement: détecte et marque les systèmes fortement couplés.

    1) Trouve le noyau fortement couplé (Δν/J < threshold)
    2) Si total ≤ max_spins : inclut tout dans le Hamiltonien
    3) Si total > max_spins : noyau seul + élargir les voisins (couplage virtuel)

    Args:
        groups: liste des groupes du JSON
        freq_mhz: fréquence du spectromètre
        max_spins: limite de spins pour le Hamiltonien (8 = 256×256, rapide)
        threshold: seuil Δν/J en dessous duquel considérer comme fortement couplé (3.0 = standard)

    Returns:
        groups modifiés avec:
        - hamiltonian_system_id: ID du système à simuler ensemble
        - hamiltonian_role: 'full' ou 'core'
        - virtual_coupling: True si couplé au noyau mais exclu (pour élargissement)
        - virtual_linewidth_extra_hz: élargissement additionnel en Hz

    Référence: Claridge, High-Resolution NMR Techniques (20-40% du J)
    """
    n = len(groups)
    if n == 0:
        return groups

    # Étape 1: Construire les matrices de couplage
    J_matrix = [[0.0] * n for _ in range(n)]
    strongly_coupled = [[False] * n for _ in range(n)]

    # Indexer les groupes par H_indices pour recherche rapide
    h_idx_to_group = {}
    for gi, group in enumerate(groups):
        for h_idx in group.get("H_indices", []):
            h_idx_to_group[h_idx] = gi

    for i in range(n):
        for j in range(i + 1, n):
            gi = groups[i]
            gj = groups[j]

            # Chercher si les groupes sont couplés (vérifier les deux sens)
            J_ij = 0.0
            # Direction i → j
            for partner in gi.get("partners", []):
                partner_h_idx = partner.get("partner_h_index")
                if partner_h_idx in gj.get("H_indices", []):
                    J_ij = abs(partner.get("J_hz", 0.0))
                    break
            # Direction j → i (si pas trouvé)
            if J_ij == 0:
                for partner in gj.get("partners", []):
                    partner_h_idx = partner.get("partner_h_index")
                    if partner_h_idx in gi.get("H_indices", []):
                        J_ij = abs(partner.get("J_hz", 0.0))
                        break

            if J_ij > 0:
                J_matrix[i][j] = J_matrix[j][i] = J_ij

                # Calculer Δν en Hz
                shift_i = gi.get("H_shifts_ppm", [0])[0] if gi.get("H_shifts_ppm") else 0
                shift_j = gj.get("H_shifts_ppm", [0])[0] if gj.get("H_shifts_ppm") else 0
                delta_hz = abs(shift_i - shift_j) * freq_mhz

                # Couplage fort si Δν/J < threshold
                if J_ij > 0 and delta_hz / J_ij < threshold:
                    strongly_coupled[i][j] = True
                    strongly_coupled[j][i] = True

    # Étape 2: Trouver les noyaux fortement couplés (BFS)
    visited = [False] * n
    systems = []

    for i in range(n):
        if visited[i] or not any(strongly_coupled[i]):
            continue

        # BFS pour trouver le cluster
        core = set()
        queue = [i]
        while queue:
            current = queue.pop(0)
            if current in core:
                continue
            core.add(current)
            for j in range(n):
                if j not in core and strongly_coupled[current][j]:
                    queue.append(j)

        for c in core:
            visited[c] = True

        # Étape 3: Trouver les voisins couplés au noyau
        neighbors = set()
        for member in core:
            for j in range(n):
                if j not in core and J_matrix[member][j] > 0:
                    neighbors.add(j)

        # Compter les spins
        core_spins = sum(len(groups[k].get("H_indices", [])) for k in core)
        neighbor_spins = sum(len(groups[k].get("H_indices", [])) for k in neighbors)
        total_spins = core_spins + neighbor_spins

        systems.append({
            "core": core,
            "neighbors": neighbors,
            "core_spins": core_spins,
            "total_spins": total_spins,
        })

    # Étape 4: Assigner les rôles
    system_id = 0
    for system in systems:
        core = system["core"]
        neighbors = system["neighbors"]
        total = system["total_spins"]

        if total <= max_spins:
            # CAS A: Tout tient → Hamiltonien complet
            all_members = core | neighbors
            for k in all_members:
                groups[k]["hamiltonian_system_id"] = system_id
                groups[k]["hamiltonian_role"] = "full"
                groups[k]["hamiltonian_total_spins"] = total
                groups[k]["second_order_type"] = "strong_coupling"
            print(f"[INFO] Virtual coupling system {system_id}: {len(all_members)} groups, {total} spins → full Hamiltonian", file=sys.stderr)
        else:
            # CAS B: Trop de spins → noyau seul + élargir voisins
            for k in core:
                groups[k]["hamiltonian_system_id"] = system_id
                groups[k]["hamiltonian_role"] = "core"
                groups[k]["hamiltonian_total_spins"] = system["core_spins"]
                groups[k]["second_order_type"] = "strong_coupling"

            for k in neighbors:
                # NE PAS inclure dans le Hamiltonien mais marquer pour élargissement
                groups[k]["virtual_coupling"] = True
                groups[k]["virtual_coupling_source"] = system_id
                # Facteur d'élargissement: ~30% du J vers le noyau (Claridge)
                J_to_core = max(J_matrix[k][m] for m in core) if core else 0.0
                groups[k]["virtual_linewidth_extra_hz"] = J_to_core * 0.3

            print(f"[INFO] Virtual coupling system {system_id}: core {len(core)} groups ({system['core_spins']} spins), "
                  f"neighbors {len(neighbors)} groups → broadened", file=sys.stderr)

        system_id += 1

    return groups


def get_effective_linewidth(group: dict, base_linewidth_hz: float) -> float:
    """
    Calcule la largeur de raie effective incluant le couplage virtuel.

    Args:
        group: groupe du JSON
        base_linewidth_hz: largeur de raie de base en Hz

    Returns:
        largeur de raie effective en Hz
    """
    extra = group.get("virtual_linewidth_extra_hz", 0.0)
    return base_linewidth_hz + extra


def simulate_spectrum_hamiltonian(
    mult_data: dict,
    freq_mhz: float,
    line_shape: str,
    ppm_min: float,
    ppm_max: float,
    digital_resolution_hz: float,
    default_lw: float,
    sim_mod,
    enable_dynamics: bool = False,
    temperature_k: float = 298.0,
    ph: float = 7.4,
    solvent: str = None,
) -> dict:
    """
    Use advanced Hamiltonian simulation from sim_mult_rulebased.py.

    This provides exact quantum mechanical simulation for:
    - AB systems (two coupled protons with similar shifts)
    - ABX systems (three coupled protons)
    - AA'BB' systems (symmetrically coupled systems)
    - Larger systems via fragmentation

    Returns None if simulation fails (triggers fallback to 1st order).
    """
    import numpy as np

    try:
        mol_id = mult_data.get("mol_id", "mol")
        smiles = mult_data.get("smiles", "")
        # Priorité: solvant passé en paramètre > solvant dans JSON > défaut CDCl3
        if solvent is None:
            solvent = mult_data.get("solvent", "CDCl3")
        groups = mult_data.get("groups", [])

        # Filter out exchangeable protons based on solvent (D2O exchanges OH/NH)
        # Le module dynamics gère automatiquement la visibilité selon le solvant
        hide_exchangeable = solvent.upper() in ('D2O', 'D20', 'DEUTERIUM_OXIDE')
        if hide_exchangeable:
            original_count = len(groups)
            # Filtrer par "note" OU par heavy_symbol (O, N, S = protons échangeables)
            def is_exchangeable_h(g):
                if g.get("note") == "exchangeable":
                    return True
                # Fallback: si heavy_symbol est O, N, ou S → proton échangeable (OH, NH, SH)
                heavy_sym = g.get("heavy_symbol", "")
                if heavy_sym in ("O", "N", "S"):
                    return True
                return False
            
            groups = [g for g in groups if not is_exchangeable_h(g)]
            filtered_count = original_count - len(groups)
            if filtered_count > 0:
                print(f"[INFO] Hamiltonian: Filtered out {filtered_count} exchangeable protons (OH, NH, SH) - solvent={solvent}", file=sys.stderr)

        # ═══════════════════════════════════════════════════════════════════════════
        # DÉTECTION DU COUPLAGE VIRTUEL (avant séparation des groupes)
        # Fusionne les groupes fortement couplés et marque les voisins pour élargissement
        # ═══════════════════════════════════════════════════════════════════════════
        max_hamiltonian_spins = mult_data.get("max_hamiltonian_spins", 8)
        virtual_coupling_threshold = mult_data.get("virtual_coupling_threshold", 3.0)
        groups = merge_strongly_coupled_groups(groups, freq_mhz, max_hamiltonian_spins, virtual_coupling_threshold)

        # Calculate spectrum dimensions
        ppm_range = ppm_max - ppm_min
        hz_range = ppm_range * freq_mhz
        n_points = int(hz_range / digital_resolution_hz) + 1
        n_points = min(n_points, 200000)

        ppm_axis = np.linspace(ppm_max, ppm_min, n_points)
        spectrum = np.zeros(n_points, dtype=np.float64)
        peaks_info = []

        # Separate groups by simulation method
        second_order_groups = []  # AA'BB', AB, ABX, etc.
        first_order_groups = []   # s, d, t, q, dd, dt, etc.

        # Compter le nombre total de protons
        total_protons = sum(len(g.get("H_indices", [])) for g in groups)

        # Option: Forcer TOUS les groupes dans l'Hamiltonien (pour systèmes couplés complets)
        force_all_hamiltonian = mult_data.get("force_all_hamiltonian", False)

        # DÉTECTION AUTOMATIQUE: Si ≤8 protons, simuler tout ensemble pour un spectre plus réaliste
        auto_force_threshold = mult_data.get("auto_force_threshold", 8)
        if not force_all_hamiltonian and total_protons <= auto_force_threshold:
            force_all_hamiltonian = True
            print(f"[INFO] Auto-Hamiltonian: {total_protons} protons ≤ {auto_force_threshold} → simulation complète activée", file=sys.stderr)

        for group in groups:
            second_order_type = group.get("second_order_type", "first_order")

            # Use Hamiltonian for second-order systems OR if force_all_hamiltonian is True
            # ⚠️ IMPORTANT: Ne pas utiliser system_type seul pour décider du mode de simulation!
            # system_type est juste un tag descriptif (ex: "saturated_diketopiperazine"),
            # pas une indication de couplage fort. Un quartet couplé à un CH3 externe
            # doit être traité en first-order, pas en Hamiltonien sans ses partenaires!
            is_truly_second_order = second_order_type in ["strong_coupling", "weak_second_order", "moderate_second_order"]

            if force_all_hamiltonian or is_truly_second_order:
                second_order_groups.append(group)
            else:
                first_order_groups.append(group)

        # Option: Inclure les groupes first-order CONNECTÉS aux second-order
        # pour simuler la chaîne alkyle complète (plus lent mais plus réaliste)
        include_connected = mult_data.get("include_connected_groups", False)
        max_hamiltonian_spins = mult_data.get("max_hamiltonian_spins", 8)  # Limite par défaut (8 = 256x256, rapide)

        # AUTO-ACTIVATION: Si des groupes fortement couplés existent, activer include_connected
        # automatiquement pour simuler le sous-système couplé complet (ex: ABX diastéréotopique)
        # Cela permet un multiplet réaliste sans simuler TOUS les protons
        # ✅ FIX: Inclure moderate_second_order (Δν/J ~ 3-5) pour systèmes ABX diastéréotopiques
        has_strong_or_moderate_coupling = any(
            g.get("second_order_type") in ["strong_coupling", "moderate_second_order"]
            for g in second_order_groups
        )
        if has_strong_or_moderate_coupling and not include_connected:
            include_connected = True
            print(f"[INFO] Auto-include connected groups: strong/moderate coupling détecté → activation", file=sys.stderr)

        if include_connected and second_order_groups:
            # Compter les spins actuels dans second_order_groups
            current_spins = sum(len(g.get("H_indices", [])) for g in second_order_groups)

            # Collecter tous les H_indices des groupes second-order
            second_order_h_indices = set()
            for g in second_order_groups:
                second_order_h_indices.update(g.get("H_indices", []))

            # Trouver les groupes first-order qui sont couplés aux second-order
            # Prioriser les groupes avec shift proche des second-order (pour le massif)
            candidates = []
            for group in first_order_groups:
                partners = group.get("partners", [])
                for p in partners:
                    partner_h_idx = p.get("partner_h_index")
                    if partner_h_idx in second_order_h_indices:
                        avg_shift = sum(group.get("H_shifts_ppm", [0])) / max(1, len(group.get("H_shifts_ppm", [1])))
                        n_protons = len(group.get("H_indices", []))
                        candidates.append((group, avg_shift, n_protons))
                        break

            # Trier par shift décroissant (groupes plus proches de la région 1.2-1.5 ppm d'abord)
            candidates.sort(key=lambda x: -x[1])

            # Ajouter jusqu'à la limite de spins
            for group, avg_shift, n_protons in candidates:
                if current_spins + n_protons <= max_hamiltonian_spins:
                    first_order_groups.remove(group)
                    second_order_groups.append(group)
                    current_spins += n_protons
                    print(f"[INFO] Including connected group {group.get('group_id')} ({group.get('label')}, {avg_shift:.2f} ppm, {n_protons}H) → total {current_spins} spins", file=sys.stderr)

        # LIMITE STRICTE: Si les groupes second-ordre dépassent max_hamiltonian_spins,
        # déplacer les groupes excédentaires vers first-order (traités analytiquement)
        total_second_order_spins = sum(len(g.get("H_indices", [])) for g in second_order_groups)
        if total_second_order_spins > max_hamiltonian_spins:
            print(f"[WARN] {total_second_order_spins} spins > max {max_hamiltonian_spins} → réduction", file=sys.stderr)
            # Trier par priorité: garder les vrais systèmes AA'BB' d'abord
            def priority(g):
                st = g.get("system_type")
                sot = g.get("second_order_type", "")
                if st in ["AA'BB'", "AB", "ABX"]:
                    return 0
                elif sot == "strong_coupling":
                    return 1
                elif sot == "weak_second_order":
                    return 2
                else:
                    return 3
            second_order_groups.sort(key=priority)

            # Garder seulement jusqu'à max_hamiltonian_spins
            kept = []
            current = 0
            for g in second_order_groups:
                n = len(g.get("H_indices", []))
                if current + n <= max_hamiltonian_spins:
                    kept.append(g)
                    current += n
                else:
                    first_order_groups.append(g)
                    print(f"[INFO] Groupe {g.get('group_id')} ({g.get('label')}) → first-order (limite spins)", file=sys.stderr)
            second_order_groups = kept

        print(f"[INFO] {len(second_order_groups)} second-order groups, {len(first_order_groups)} first-order groups", file=sys.stderr)
        print(f"[DEBUG] Second-order groups: {len(second_order_groups)}, spins: {sum(len(g.get('H_indices', [])) for g in second_order_groups)}", file=sys.stderr)

        # ═══════════════════════════════════════════════════════════════════════════
        # FRAGMENTATION AU NIVEAU GROUPE: Identifier les systèmes de spins isolés
        # Chaque fragment sera simulé indépendamment avec son propre Hamiltonian
        # ═══════════════════════════════════════════════════════════════════════════
        use_group_fragmentation = mult_data.get("use_group_fragmentation", True)
        group_fragments = []

        if use_group_fragmentation and len(second_order_groups) > 1:
            # Fragmenter les groupes second-order basé sur la connectivité partners
            all_fragments = find_connected_group_fragments(second_order_groups)

            # Traiter chaque fragment: vérifier la limite de spins
            for frag in all_fragments:
                frag_spins = sum(len(g.get("H_indices", [])) for g in frag)
                if frag_spins <= max_hamiltonian_spins:
                    # Fragment OK, garder tel quel
                    group_fragments.append(frag)
                else:
                    # Fragment trop grand, appliquer la réduction
                    print(f"[WARN] Fragment {len(frag)} groups, {frag_spins} spins > max {max_hamiltonian_spins} → réduction", file=sys.stderr)

                    # Trier par priorité
                    def priority(g):
                        st = g.get("system_type")
                        sot = g.get("second_order_type", "")
                        if st in ["AA'BB'", "AB", "ABX"]:
                            return 0
                        elif sot == "strong_coupling":
                            return 1
                        elif sot == "weak_second_order":
                            return 2
                        else:
                            return 3
                    frag.sort(key=priority)

                    # Garder jusqu'à max_hamiltonian_spins
                    kept = []
                    current = 0
                    for g in frag:
                        n = len(g.get("H_indices", []))
                        if current + n <= max_hamiltonian_spins:
                            kept.append(g)
                            current += n
                        else:
                            first_order_groups.append(g)
                            print(f"[INFO] Groupe {g.get('group_id')} ({g.get('label')}) → first-order (fragment trop grand)", file=sys.stderr)

                    if kept:
                        group_fragments.append(kept)
        else:
            # Pas de fragmentation: un seul fragment avec tous les groupes
            if second_order_groups:
                group_fragments.append(second_order_groups)

        print(f"[INFO] {len(group_fragments)} fragment(s) Hamiltonian à simuler", file=sys.stderr)
        for i, frag in enumerate(group_fragments):
            n_spins = sum(len(g.get("H_indices", [])) for g in frag)
            gids = [g.get("group_id") for g in frag]
            print(f"[INFO]   Fragment {i+1}: {len(frag)} groups, {n_spins} spins, IDs={gids}", file=sys.stderr)

        # Define default linewidth for both parts (default_lw is already passed as parameter)
        if default_lw is None:
            default_lw = mult_data.get("typical_linewidth_hz", 1.0)

        # PART 1: Handle first-order systems with analytical multiplets
        for group in first_order_groups:
            shifts = group.get("H_shifts_ppm", [])
            if not shifts:
                continue

            ppm_center = float(shifts[0])
            label = group.get("label", "s")
            j_couplings = group.get("J_couplings", [])
            j_hz_list = []

            for jc in j_couplings:
                j_val = jc.get("J_hz", 0)
                n_equiv = jc.get("n_equiv", 1)
                if abs(j_val) > 0.1:
                    for _ in range(n_equiv):
                        j_hz_list.append(j_val)

            # Generate analytical multiplet pattern
            pattern = generate_multiplet_from_label_and_j(
                center_ppm=ppm_center,
                label=label,
                j_couplings=j_hz_list,
                freq_mhz=freq_mhz,
                n_h=len(group.get("H_indices", [1]))
            )

            # Add to spectrum (apply virtual coupling broadening if marked)
            effective_lw = get_effective_linewidth(group, default_lw)
            lw_ppm = effective_lw / freq_mhz
            for peak_ppm, rel_intensity in pattern:
                add_peak_to_spectrum(
                    spectrum, ppm_axis,
                    peak_ppm,
                    rel_intensity * len(group.get("H_indices", [1])),
                    lw_ppm,
                    line_shape
                )

            # Build peak info with virtual coupling metadata
            peak_info = {
                "ppm": ppm_center,
                "intensity": float(len(group.get("H_indices", [1]))),
                "label": label,
                "multiplicity": label_to_name(label),
                "J_hz": j_hz_list,
                "n_protons": len(group.get("H_indices", [1])),
                "group_id": group.get("group_id", -1),
                "H_indices": group.get("H_indices", []),
                "note": group.get("note", "First-order analytical"),
            }
            # Add virtual coupling metadata if present
            if group.get("virtual_coupling"):
                peak_info["virtual_coupling"] = True
                peak_info["virtual_linewidth_extra_hz"] = group.get("virtual_linewidth_extra_hz", 0.0)
                peak_info["effective_linewidth_hz"] = effective_lw
            peaks_info.append(peak_info)

        # PART 2: Handle second-order systems with Hamiltonian simulation
        # Simuler chaque fragment indépendamment
        if group_fragments and hasattr(sim_mod, 'simulate_with_fragmentation'):
            for frag_idx, fragment_groups in enumerate(group_fragments):
                print(f"[INFO] Simulating fragment {frag_idx + 1}/{len(group_fragments)}...", file=sys.stderr)

                # Build spin system from this fragment's groups only
                all_shifts_hz = []
                group_to_spins = {}

                spin_idx = 0
                for group in fragment_groups:
                    shifts = group.get("H_shifts_ppm", [])
                    h_indices = group.get("H_indices", [])
                    if not shifts:
                        continue

                    # ✅ CORRECTION: Créer TOUS les spins même s'ils sont équivalents
                    # Pour un CH3 avec H_indices:[22,23,24] et shifts:[0.88,0.88,0.88],
                    # il faut créer 3 spins à 0.88 ppm pour que le simulateur Hamiltonien
                    # puisse calculer correctement le système A3X2 (triplet)
                    n_protons = len(h_indices)
                    n_shifts = len(shifts)

                    group_spins = []
                    if n_shifts == 1 and n_protons > 1:
                        # Un seul shift pour plusieurs protons équivalents → créer n_protons spins
                        shift_hz = shifts[0] * freq_mhz
                        for _ in range(n_protons):
                            all_shifts_hz.append(shift_hz)
                            group_spins.append(spin_idx)
                            spin_idx += 1
                    else:
                        # Cas normal: un shift par proton (ou diastéréotopique)
                        for shift_ppm in shifts:
                            shift_hz = shift_ppm * freq_mhz
                            all_shifts_hz.append(shift_hz)
                            group_spins.append(spin_idx)
                            spin_idx += 1

                    group_id = group.get("group_id", len(group_to_spins))
                    group_to_spins[group_id] = (group, group_spins)

                n_spins = len(all_shifts_hz)
                print(f"[DEBUG] Fragment {frag_idx + 1}: {len(fragment_groups)} groups, {n_spins} spins", file=sys.stderr)

                # Only simulate if there are spins
                if n_spins > 0:
                    # Build J-coupling matrix
                    j_matrix = np.zeros((n_spins, n_spins))

                    # Build H_index to spin_idx mapping
                    h_idx_to_spin = {}
                    for group_id, (group, spins) in group_to_spins.items():
                        h_indices = group.get("H_indices", [])
                        for i, h_idx in enumerate(h_indices):
                            if i < len(spins):
                                h_idx_to_spin[h_idx] = spins[i]

                    # Add couplings from partners (using partner_h_index)
                    for group_id, (group, spins) in group_to_spins.items():
                        partners = group.get("partners", [])
                        for partner_info in partners:
                            partner_h_idx = partner_info.get("partner_h_index")
                            j_hz = partner_info.get("J_hz", 0)

                            if partner_h_idx is not None and abs(j_hz) > 0.1:
                                partner_spin_idx = h_idx_to_spin.get(partner_h_idx)
                                if partner_spin_idx is not None:
                                    # Add coupling between all spins in this group and the partner
                                    for i in spins:
                                        if j_matrix[i, partner_spin_idx] == 0:  # Avoid overwriting
                                            j_matrix[i, partner_spin_idx] = j_hz
                                            j_matrix[partner_spin_idx, i] = j_hz

                    # Also add geminal couplings between diastereotopic pairs (same heavy_idx)
                    heavy_to_spins = {}
                    for group_id, (group, spins) in group_to_spins.items():
                        heavy_idx = group.get("heavy_idx")
                        if heavy_idx is not None:
                            if heavy_idx not in heavy_to_spins:
                                heavy_to_spins[heavy_idx] = []
                            heavy_to_spins[heavy_idx].extend(spins)

                    # Add geminal J for CH2 diastereotopic pairs
                    for heavy_idx, spins_list in heavy_to_spins.items():
                        if len(spins_list) == 2:
                            # Diastereotopic pair - add geminal coupling
                            j_gem = 12.5  # Typical geminal J for CH2
                            i, j = spins_list[0], spins_list[1]
                            if j_matrix[i, j] == 0:  # Don't overwrite if already set
                                j_matrix[i, j] = -j_gem  # Geminal J is typically negative
                                j_matrix[j, i] = -j_gem

                    # Call the advanced simulator
                    sim_result = sim_mod.simulate_with_fragmentation(
                        shifts_hz=np.array(all_shifts_hz),
                        j_matrix=j_matrix,
                        freq_mhz=freq_mhz,
                        linewidth_hz=default_lw,
                        ppm_axis=ppm_axis,
                        line_shape=line_shape,
                        enable_dynamics=enable_dynamics,
                        ph=ph,
                    )

                    if sim_result is not None and hasattr(sim_result, 'spectrum'):
                        # COMBINE first-order and second-order spectra (don't overwrite!)
                        spectrum += sim_result.spectrum

                    # Extract peaks from simulation result
                    if hasattr(sim_result, 'peaks') and sim_result.peaks:
                        for peak in sim_result.peaks:
                            peaks_info.append({
                                "ppm": peak.get("ppm", 0),
                                "intensity": peak.get("intensity", 1.0),
                                "label": peak.get("label", "m"),
                                "multiplicity": peak.get("multiplicity", "multiplet"),
                                "J_hz": peak.get("J_hz", []),
                                "n_protons": peak.get("n_protons", 1),
                                "group_id": peak.get("group_id", -1),
                                "H_indices": peak.get("H_indices", []),
                                "note": peak.get("note", ""),
                            })
                    else:
                        # Generate peak info from groups in this fragment
                        for group in fragment_groups:
                            shifts = group.get("H_shifts_ppm", [])
                            if shifts:
                                # ✅ CORRIGÉ: Extraire les J_couplings du groupe JSON
                                j_couplings = group.get("J_couplings", [])
                                j_hz_list = []
                                for jc in j_couplings:
                                    j_val = jc.get("J_hz", 0)
                                    n_equiv = jc.get("n_equiv", 1)
                                    if abs(j_val) > 0.1:
                                        # Ajouter le J autant de fois que n_equiv
                                        for _ in range(n_equiv):
                                            j_hz_list.append(j_val)

                                peaks_info.append({
                                    "ppm": float(shifts[0]),
                                    "intensity": float(len(group.get("H_indices", [1]))),
                                    "label": group.get("label", "m"),
                                    "multiplicity": label_to_name(group.get("label", "m")),
                                    "J_hz": j_hz_list,
                                    "n_protons": len(group.get("H_indices", [1])),
                                    "group_id": group.get("group_id", -1),
                                    "H_indices": group.get("H_indices", []),
                                    "note": group.get("note", "Hamiltonian simulation"),
                                })

        # Normalize spectrum (whether from first-order, second-order, or both)
        max_intensity = np.max(spectrum)
        if max_intensity > 0:
            spectrum = spectrum / max_intensity

        # Return combined spectrum
        return {
            "status": "ok",
            "mol_id": mol_id,
            "smiles": smiles,
            "freq_mhz": freq_mhz,
            "solvent": solvent,
            "line_shape": line_shape,
            "ppm": ppm_axis.tolist(),
            "intensity": spectrum.tolist(),
            "peaks": sorted(peaks_info, key=lambda p: -p["ppm"]),
            "n_points": n_points,
            "n_groups": len(groups),
            "total_protons": sum(p.get("n_protons", 1) for p in peaks_info),
            "ppm_range": [ppm_min, ppm_max],
            "simulation_method": "hybrid_first_and_second_order",
        }

    except Exception as e:
        print(f"[ERROR] Hamiltonian simulation failed: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        return None


def generate_spectrum_from_multiplicity_json(
    mult_data: dict,
    freq_mhz: float = 700.0,
    line_shape: str = "lorentzian",
    ppm_min: float = -0.5,
    ppm_max: float = 12.0,
    digital_resolution_hz: float = 0.1,
    line_width_hz: float = 1.0,
    enable_dynamics: bool = False,
    temperature_k: float = 298.0,
    ph: float = 7.4,
    solvent: str = None,
) -> dict:
    """
    Generate a complete 1H NMR spectrum from multiplicity JSON data.

    Uses ADVANCED HAMILTONIAN SIMULATION for realistic spectra:
    - Exact quantum mechanical simulation for strongly coupled systems (AB, ABX, AA'BB')
    - Adaptive line widths based on proton type and molecular environment
    - First-order approximation only for simple uncoupled systems

    Uses the groups data with:
    - H_shifts_ppm: chemical shifts
    - label: multiplicity (s, d, t, dd, etc.)
    - J_couplings: coupling constants in Hz
    - partners: coupling partners info
    - note: "exchangeable" for OH, NH protons (auto-filtered based on solvent, e.g., D2O)
    """
    import numpy as np

    # Try to import the advanced simulator
    try:
        import sim_mult_rulebased as sim_mod
        use_advanced_sim = True
        print("[INFO] Using advanced Hamiltonian simulator", file=sys.stderr)
    except ImportError:
        sim_mod = None
        use_advanced_sim = False
        print("[WARN] Advanced simulator not available, using first-order approximation", file=sys.stderr)

    # Extract metadata
    mol_id = mult_data.get("mol_id", "mol")
    smiles = mult_data.get("smiles", "")
    # Priorité: solvant passé en paramètre > solvant dans JSON > défaut CDCl3
    if solvent is None:
        solvent = mult_data.get("solvent", "CDCl3")
    groups = mult_data.get("groups", [])

    # DEBUG: Log all notes in groups
    notes_found = [(g.get("group_id"), g.get("note"), g.get("heavy_symbol"), g.get("H_shifts_ppm", [])) for g in groups]
    print(f"[DEBUG] Solvent={solvent}, Groups with notes: {notes_found}", file=sys.stderr)

    # Filter out exchangeable protons based on solvent (D2O exchanges OH/NH)
    # Le module dynamics gère automatiquement la visibilité selon le solvant
    hide_exchangeable = solvent.upper() in ('D2O', 'D20', 'DEUTERIUM_OXIDE')
    print(f"[DEBUG] hide_exchangeable={hide_exchangeable} (solvent.upper()={solvent.upper()})", file=sys.stderr)
    if hide_exchangeable:
        original_count = len(groups)
        # Filtrer par "note" OU par heavy_symbol (O, N, S = protons échangeables)
        def is_exchangeable(g):
            if g.get("note") == "exchangeable":
                return True
            # Fallback: si heavy_symbol est O, N, ou S → proton échangeable (OH, NH, SH)
            heavy_sym = g.get("heavy_symbol", "")
            if heavy_sym in ("O", "N", "S"):
                return True
            return False
        
        groups = [g for g in groups if not is_exchangeable(g)]
        filtered_count = original_count - len(groups)
        print(f"[INFO] Filtered out {filtered_count} exchangeable protons (OH, NH, SH) - solvent={solvent}", file=sys.stderr)

    # =========================================================================
    # PROTONS DIASTÉRÉOTOPIQUES / VINYLIQUES
    # =========================================================================
    # Les protons sur le même carbone (même heavy_idx) mais chimiquement 
    # non-équivalents (vinyl_nonequiv, diastereotopic) restent SÉPARÉS.
    # Chaque proton garde son shift et ses J individuels.
    # Pas de fusion - c'est scientifiquement correct.
    # =========================================================================

    # Override freq if specified in JSON
    json_freq = mult_data.get("freq_MHz", mult_data.get("freq_mhz"))
    if json_freq:
        freq_mhz = float(json_freq)

    # Auto-calculate PPM range from shifts if groups exist
    if groups and len(groups) > 0:
        all_shifts = []
        for group in groups:
            shifts = group.get("H_shifts_ppm", [])
            if shifts:
                all_shifts.extend([float(s) for s in shifts])

        if all_shifts:
            min_shift = min(all_shifts)
            max_shift = max(all_shifts)
            # Add margin (10% on each side, minimum 0.5 ppm)
            margin = max(0.5, (max_shift - min_shift) * 0.1)
            ppm_min = max(-1.0, min_shift - margin)
            ppm_max = min(15.0, max_shift + margin)
            print(f"[INFO] Auto PPM range: {ppm_min:.2f} - {ppm_max:.2f} ppm", file=sys.stderr)

    # Calculate spectrum dimensions
    ppm_range = ppm_max - ppm_min
    hz_range = ppm_range * freq_mhz
    n_points = int(hz_range / digital_resolution_hz) + 1
    n_points = min(n_points, 200000)  # Cap at 200k points

    ppm_axis = np.linspace(ppm_max, ppm_min, n_points)  # Reversed (NMR convention)
    spectrum = np.zeros(n_points, dtype=np.float64)
    peaks_info = []

    # Use linewidth from parameter (user override), fallback to JSON value
    # PRIORITÉ: line_width_hz (paramètre utilisateur) > typical_linewidth_hz (JSON)
    default_lw = line_width_hz if line_width_hz is not None else mult_data.get("typical_linewidth_hz", 1.0)

    # TRY ADVANCED SIMULATION IF AVAILABLE
    if use_advanced_sim and hasattr(sim_mod, 'simulate_with_fragmentation'):
        try:
            # Build spin system from groups
            result = simulate_spectrum_hamiltonian(
                mult_data, freq_mhz, line_shape,
                ppm_min, ppm_max, digital_resolution_hz,
                default_lw, sim_mod, enable_dynamics, temperature_k, ph, solvent
            )
            if result is not None:
                return result
            # If failed, fall back to simple simulation below
            print("[WARN] Hamiltonian simulation failed, falling back to 1st order", file=sys.stderr)
        except Exception as e:
            print(f"[WARN] Hamiltonian simulation error: {e}, falling back to 1st order", file=sys.stderr)

    # FALLBACK: Simple first-order simulation
    # Process each group
    for group in groups:
        # Get shift(s) for this group
        shifts = group.get("H_shifts_ppm", [])
        if not shifts:
            continue

        # Use the first shift as the center (all H in group are equivalent)
        center_ppm = float(shifts[0])
        if not (ppm_min <= center_ppm <= ppm_max):
            continue

        # Get multiplicity info
        label = group.get("label", "s")

        # Get J couplings - try multiple sources
        j_couplings = []

        # Source 1: J_couplings array
        if group.get("J_couplings"):
            for jc in group["J_couplings"]:
                if isinstance(jc, dict):
                    j_hz = jc.get("J_hz", jc.get("j_hz", jc.get("J", 0)))
                    n_equiv = jc.get("n_equiv", 1)  # Nombre de protons équivalents couplés
                    print(f"[DEBUG] Group {group.get('group_id')}: label={label}, J_hz={j_hz}, n_equiv={n_equiv}", file=sys.stderr)
                    if j_hz and abs(float(j_hz)) > 0.1:
                        # Ajouter J_hz autant de fois que n_equiv pour générer le bon multiplet
                        # triplet (t) = couplage avec 2 H équivalents = [J, J]
                        # quartet (q) = couplage avec 3 H équivalents = [J, J, J]
                        for _ in range(int(n_equiv)):
                            j_couplings.append(float(j_hz))
                        print(f"[DEBUG] Generated j_couplings={j_couplings}", file=sys.stderr)
                elif isinstance(jc, (int, float)):
                    if abs(float(jc)) > 0.1:
                        j_couplings.append(float(jc))

        # Source 2: partners array
        if not j_couplings and group.get("partners"):
            for partner in group["partners"]:
                j_hz = partner.get("J_hz", partner.get("j_hz", 0))
                if j_hz and abs(float(j_hz)) > 0.1:
                    j_couplings.append(float(j_hz))

        # Get number of equivalent protons
        n_h = len(group.get("H_indices", [1]))

        # Get linewidth (adaptive if available)
        lw_hz = default_lw

        # Generate multiplet pattern based on label and J couplings
        pattern = generate_multiplet_from_label_and_j(
            center_ppm,
            label,
            j_couplings,
            freq_mhz,
            n_h
        )

        # Add each line of the multiplet to the spectrum
        lw_ppm = lw_hz / freq_mhz
        for line_ppm, line_intensity in pattern:
            add_peak_to_spectrum(
                spectrum, ppm_axis,
                line_ppm, line_intensity * n_h,  # Scale by number of protons
                lw_ppm,
                line_shape
            )

        # Store peak info
        peaks_info.append({
            "ppm": center_ppm,
            "intensity": float(n_h),
            "label": label,
            "multiplicity": label_to_name(label),
            "J_hz": j_couplings,
            "n_protons": n_h,
            "group_id": group.get("group_id", -1),
            "H_indices": group.get("H_indices", []),
            "note": group.get("note", ""),
        })

    # Normalize spectrum
    max_intensity = np.max(spectrum)
    if max_intensity > 0:
        spectrum = spectrum / max_intensity

    # Spectre propre sans bruit ni baseline artificielle

    return {
        "status": "ok",
        "mol_id": mol_id,
        "smiles": smiles,
        "freq_mhz": freq_mhz,
        "solvent": solvent,
        "line_shape": line_shape,
        "ppm": ppm_axis.tolist(),
        "intensity": spectrum.tolist(),
        "peaks": sorted(peaks_info, key=lambda p: -p["ppm"]),  # Sort by ppm descending
        "n_points": n_points,
        "n_groups": len(groups),
        "total_protons": sum(p.get("n_protons", 1) for p in peaks_info),
        "ppm_range": [ppm_min, ppm_max],
        "simulation_method": "first_order_fallback",
    }


def generate_multiplet_from_label_and_j(
    center_ppm: float, 
    label: str, 
    j_couplings: list, 
    freq_mhz: float,
    n_h: int = 1
) -> list:
    """
    Generate multiplet pattern based on multiplicity label and J couplings.
    
    Returns list of (ppm, relative_intensity) tuples.
    """
    import numpy as np
    
    # If no J couplings, generate pattern from label
    if not j_couplings:
        j_couplings = estimate_j_from_label(label)
    
    if not j_couplings:
        return [(center_ppm, 1.0)]
    
    # Start with single line at center
    lines = [(0.0, 1.0)]  # (offset_hz, intensity)
    
    # Apply each coupling to split the pattern
    for j_hz in j_couplings:
        if abs(j_hz) < 0.1:
            continue
        new_lines = []
        for offset, intensity in lines:
            # Split into two lines
            new_lines.append((offset - j_hz / 2, intensity / 2))
            new_lines.append((offset + j_hz / 2, intensity / 2))
        lines = new_lines
    
    # Merge overlapping lines (within 1 Hz to handle floating point errors)
    lines = merge_close_lines(lines, tolerance_hz=1.0)
    
    # Convert to ppm
    pattern = []
    for offset_hz, intensity in lines:
        ppm = center_ppm + offset_hz / freq_mhz
        pattern.append((ppm, intensity))
    
    return pattern


def estimate_j_from_label(label: str) -> list:
    """
    Estimate J couplings from multiplicity label when not provided.
    Uses typical J values.
    """
    label_lower = label.lower().strip()
    
    # Common patterns
    patterns = {
        "s": [],
        "d": [7.0],
        "t": [7.0, 7.0],
        "q": [7.0, 7.0, 7.0],
        "quint": [7.0] * 4,
        "sext": [7.0] * 5,
        "sept": [7.0] * 6,
        "dd": [10.0, 5.0],
        "dt": [10.0, 7.0, 7.0],
        "dq": [10.0, 7.0, 7.0, 7.0],
        "td": [7.0, 7.0, 3.0],
        "tt": [7.0, 7.0, 3.0, 3.0],
        "ddd": [12.0, 8.0, 4.0],
        "dddd": [12.0, 8.0, 5.0, 2.0],
        "m": [],  # Multiplet - no simple pattern
        "br s": [],
        "br d": [7.0],
    }
    
    return patterns.get(label_lower, [])


def merge_close_lines(lines: list, tolerance_hz: float = 0.3) -> list:
    """Merge lines that are very close together."""
    if not lines:
        return lines
    
    # Sort by offset
    sorted_lines = sorted(lines, key=lambda x: x[0])
    merged = [sorted_lines[0]]
    
    for offset, intensity in sorted_lines[1:]:
        last_offset, last_intensity = merged[-1]
        if abs(offset - last_offset) < tolerance_hz:
            # Merge with previous - weighted average position
            total_int = last_intensity + intensity
            if total_int > 0:
                new_offset = (last_offset * last_intensity + offset * intensity) / total_int
            else:
                new_offset = (last_offset + offset) / 2
            merged[-1] = (new_offset, total_int)
        else:
            merged.append((offset, intensity))
    
    return merged


def add_peak_to_spectrum(spectrum, ppm_axis, center_ppm, intensity, width_ppm, shape_type):
    """
    Add a peak shape to the spectrum at the given position.
    """
    import numpy as np
    
    if shape_type == "lorentzian":
        # Lorentzian: L(x) = (gamma/2)^2 / ((x-x0)^2 + (gamma/2)^2)
        gamma = width_ppm
        peak = intensity * (gamma/2)**2 / ((ppm_axis - center_ppm)**2 + (gamma/2)**2)
    elif shape_type == "gaussian":
        # Gaussian: G(x) = exp(-0.5 * ((x-x0)/sigma)^2)
        sigma = width_ppm / 2.355  # FWHM to sigma
        peak = intensity * np.exp(-0.5 * ((ppm_axis - center_ppm) / sigma)**2)
    else:  # voigt (pseudo-Voigt approximation optimized for overlapping signals)
        # 70% Lorentzian + 30% Gaussian for better overlap visualization
        gamma = width_ppm
        sigma = width_ppm / 2.355
        lorentz = (gamma/2)**2 / ((ppm_axis - center_ppm)**2 + (gamma/2)**2)
        gauss = np.exp(-0.5 * ((ppm_axis - center_ppm) / sigma)**2)
        peak = intensity * (0.7 * lorentz + 0.3 * gauss)
    
    spectrum += peak


def add_realistic_baseline(spectrum: 'np.ndarray', noise_level: float = 0.002,
                          baseline_amplitude: float = 0.005, n_waves: int = 3) -> 'np.ndarray':
    """
    Add realistic baseline artifacts to NMR spectrum.

    Args:
        spectrum: Input spectrum array
        noise_level: Random noise amplitude (default 0.002 = 0.2%)
        baseline_amplitude: Baseline wave amplitude (default 0.005 = 0.5%)
        n_waves: Number of baseline waves (default 3)

    Returns:
        Spectrum with added baseline and noise
    """
    import numpy as np

    n_points = len(spectrum)

    # 1. White noise (high frequency)
    white_noise = np.random.normal(0, noise_level, n_points)

    # 2. Baseline drift (low frequency waves)
    baseline = np.zeros(n_points)
    x = np.linspace(0, 2 * np.pi * n_waves, n_points)

    # Add multiple sine/cosine waves with random phases
    for i in range(n_waves):
        phase = np.random.uniform(0, 2 * np.pi)
        amplitude = baseline_amplitude * np.random.uniform(0.5, 1.0)
        baseline += amplitude * np.sin(x * (i + 1) / n_waves + phase)

    # 3. Combine spectrum + baseline + noise
    result = spectrum + baseline + white_noise

    # Ensure no negative values
    result = np.maximum(result, 0)

    return result


def label_to_name(label: str) -> str:
    """Convert multiplicity label to full name."""
    mapping = {
        "s": "singlet",
        "d": "doublet",
        "t": "triplet",
        "q": "quartet",
        "quint": "quintet",
        "sext": "sextet",
        "sept": "septet",
        "oct": "octet",
        "dd": "doublet of doublets",
        "dt": "doublet of triplets",
        "dq": "doublet of quartets",
        "td": "triplet of doublets",
        "tt": "triplet of triplets",
        "ddd": "doublet of doublet of doublets",
        "dddd": "doublet of doublet of doublet of doublets",
        "ddt": "doublet of doublet of triplets",
        "dtd": "doublet of triplet of doublets",
        "m": "multiplet",
        "br s": "broad singlet",
        "br d": "broad doublet",
        "br t": "broad triplet",
        "app t": "apparent triplet",
        "app q": "apparent quartet",
    }
    return mapping.get(label.lower().strip(), label)


if __name__ == "__main__":
    main()
