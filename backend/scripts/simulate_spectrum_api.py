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
    parser.add_argument("--line_shape", default="lorentzian", choices=["lorentzian", "gaussian", "voigt"], help="Line shape")
    parser.add_argument("--digital_resolution_hz", type=float, default=0.1, help="Digital resolution in Hz")
    parser.add_argument("--ppm_min", type=float, default=-0.5, help="Minimum ppm")
    parser.add_argument("--ppm_max", type=float, default=12.0, help="Maximum ppm")
    parser.add_argument("--line_width_hz", type=float, default=1.0, help="Default line width in Hz")
    
    args = parser.parse_args()
    
    try:
        import numpy as np
        
        # Determine the source of multiplicity data
        mult_data = None
        
        # Priority 1: JSON file path
        if args.json_path and os.path.exists(args.json_path):
            with open(args.json_path, "r", encoding="utf-8") as f:
                mult_data = json.load(f)
                
        # Priority 2: Inline JSON data
        elif args.json_data:
            mult_data = json.loads(args.json_data)
            
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


def simulate_spectrum_hamiltonian(
    mult_data: dict,
    freq_mhz: float,
    line_shape: str,
    ppm_min: float,
    ppm_max: float,
    digital_resolution_hz: float,
    default_lw: float,
    sim_mod
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
        solvent = mult_data.get("solvent", "CDCl3")
        groups = mult_data.get("groups", [])

        # Calculate spectrum dimensions
        ppm_range = ppm_max - ppm_min
        hz_range = ppm_range * freq_mhz
        n_points = int(hz_range / digital_resolution_hz) + 1
        n_points = min(n_points, 200000)

        ppm_axis = np.linspace(ppm_max, ppm_min, n_points)
        spectrum = np.zeros(n_points, dtype=np.float64)
        peaks_info = []

        # Use simulate_with_fragmentation if available (handles large systems)
        if hasattr(sim_mod, 'simulate_with_fragmentation'):
            # Build spin system from groups
            all_shifts_hz = []
            group_to_spins = {}

            spin_idx = 0
            for group in groups:
                shifts = group.get("H_shifts_ppm", [])
                if not shifts:
                    continue

                group_spins = []
                for shift_ppm in shifts:
                    shift_hz = shift_ppm * freq_mhz
                    all_shifts_hz.append(shift_hz)
                    group_spins.append(spin_idx)
                    spin_idx += 1

                group_id = group.get("group_id", len(group_to_spins))
                group_to_spins[group_id] = (group, group_spins)

            n_spins = len(all_shifts_hz)
            if n_spins == 0:
                return None

            # Build J-coupling matrix
            j_matrix = np.zeros((n_spins, n_spins))

            # Add couplings from partners
            for group_id, (group, spins) in group_to_spins.items():
                partners = group.get("partners", [])
                for partner_info in partners:
                    partner_group_id = partner_info.get("group_id")
                    j_hz = partner_info.get("J_hz", 0)

                    if partner_group_id in group_to_spins and abs(j_hz) > 0.1:
                        partner_spins = group_to_spins[partner_group_id][1]
                        # Add coupling between all spins in these groups
                        for i in spins:
                            for j in partner_spins:
                                j_matrix[i, j] = j_hz
                                j_matrix[j, i] = j_hz

            # Call the advanced simulator
            sim_result = sim_mod.simulate_with_fragmentation(
                shifts_hz=np.array(all_shifts_hz),
                j_matrix=j_matrix,
                freq_mhz=freq_mhz,
                linewidth_hz=default_lw,
                ppm_axis=ppm_axis,
                line_shape=line_shape
            )

            if sim_result is not None and hasattr(sim_result, 'spectrum'):
                spectrum = sim_result.spectrum

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
                    # Generate peak info from groups
                    for group in groups:
                        shifts = group.get("H_shifts_ppm", [])
                        if shifts:
                            peaks_info.append({
                                "ppm": float(shifts[0]),
                                "intensity": float(len(group.get("H_indices", [1]))),
                                "label": group.get("label", "m"),
                                "multiplicity": label_to_name(group.get("label", "m")),
                                "J_hz": [],
                                "n_protons": len(group.get("H_indices", [1])),
                                "group_id": group.get("group_id", -1),
                                "H_indices": group.get("H_indices", []),
                                "note": group.get("note", "Hamiltonian simulation"),
                            })

                # Normalize
                max_intensity = np.max(spectrum)
                if max_intensity > 0:
                    spectrum = spectrum / max_intensity

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
                    "simulation_method": "hamiltonian",
                }

        return None  # Fallback to 1st order

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
    solvent = mult_data.get("solvent", "CDCl3")
    groups = mult_data.get("groups", [])

    # Override freq if specified in JSON
    json_freq = mult_data.get("freq_MHz", mult_data.get("freq_mhz"))
    if json_freq:
        freq_mhz = float(json_freq)

    # Calculate spectrum dimensions
    ppm_range = ppm_max - ppm_min
    hz_range = ppm_range * freq_mhz
    n_points = int(hz_range / digital_resolution_hz) + 1
    n_points = min(n_points, 200000)  # Cap at 200k points

    ppm_axis = np.linspace(ppm_max, ppm_min, n_points)  # Reversed (NMR convention)
    spectrum = np.zeros(n_points, dtype=np.float64)
    peaks_info = []

    # Get linewidth from JSON if available
    default_lw = mult_data.get("typical_linewidth_hz", line_width_hz)

    # TRY ADVANCED SIMULATION IF AVAILABLE
    if use_advanced_sim and hasattr(sim_mod, 'simulate_with_fragmentation'):
        try:
            # Build spin system from groups
            result = simulate_spectrum_hamiltonian(
                mult_data, freq_mhz, line_shape,
                ppm_min, ppm_max, digital_resolution_hz,
                default_lw, sim_mod
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
                    if j_hz and abs(float(j_hz)) > 0.1:
                        j_couplings.append(float(j_hz))
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
    
    # Merge overlapping lines (within 0.5 Hz)
    lines = merge_close_lines(lines, tolerance_hz=0.3)
    
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
    else:  # voigt (pseudo-Voigt approximation)
        # 50% Lorentzian + 50% Gaussian
        gamma = width_ppm
        sigma = width_ppm / 2.355
        lorentz = (gamma/2)**2 / ((ppm_axis - center_ppm)**2 + (gamma/2)**2)
        gauss = np.exp(-0.5 * ((ppm_axis - center_ppm) / sigma)**2)
        peak = intensity * 0.5 * (lorentz + gauss)
    
    spectrum += peak


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
