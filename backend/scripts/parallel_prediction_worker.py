#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
parallel_prediction_worker.py
==============================
Worker function for parallel molecule prediction using multiprocessing.
Each worker loads its own copy of the models and processes molecules independently.

IMPORTANT: This module must be importable as a separate module for multiprocessing.
Do not use relative imports - use absolute imports or add parent to sys.path.
"""

import os
import sys
import traceback
from typing import Dict, Any, Optional, Tuple

# Global worker-local model cache to avoid reloading models for each molecule
_WORKER_MODELS = None

def _initialize_worker_models(model_1h_path: str, model_13c_path: str):
    """Initialize models once per worker process."""
    global _WORKER_MODELS
    if _WORKER_MODELS is not None:
        return _WORKER_MODELS

    # Set CPU-only mode
    os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
    os.environ["DGLBACKEND"] = "pytorch"

    # Import torch and disable CUDA
    import torch
    torch.backends.cudnn.enabled = False

    # Import the main module's classes
    # We need to be careful about imports in multiprocessing context
    from integrated_nmr_prediction import ModelBundle

    # Create a fresh model bundle for this worker
    _WORKER_MODELS = {
        'bundle': ModelBundle(net_1h=None, net_13c=None, mean_std_1h=None, mean_std_13c=None),
        'model_1h_path': model_1h_path,
        'model_13c_path': model_13c_path,
    }

    return _WORKER_MODELS


def worker_process_molecule(
    mol_data: Dict[str, Any],
    out_dir: str,
    model_1h_path: str,
    model_13c_path: str,
    property_key: str = "(1)",
    sim_params: Optional[dict] = None,
) -> Tuple[bool, Dict[str, Any], str]:
    """
    Worker function to process a single molecule in a separate process.

    Args:
        mol_data: Dict with 'mol_id' and 'smiles'
        out_dir: Output directory
        model_1h_path: Path to 1H model
        model_13c_path: Path to 13C model
        property_key: SDF property key
        sim_params: Simulation parameters

    Returns:
        (success: bool, result: dict, mol_id: str)
    """
    mol_id = mol_data['mol_id']
    smiles = mol_data['smiles']

    try:
        # Initialize worker models (cached after first call)
        worker_models = _initialize_worker_models(model_1h_path, model_13c_path)

        # Import process_one_smiles
        from integrated_nmr_prediction import process_one_smiles

        # Try to import multiplicity module
        multiplicity_mod = None
        try:
            import sim_mult_rulebased as multiplicity_mod
        except ImportError:
            pass

        # Process the molecule using worker's model bundle
        sdf_path, nH, nC, models_updated, _, payload = process_one_smiles(
            smiles=smiles,
            mol_id=mol_id,
            out_dir=out_dir,
            models=worker_models['bundle'],
            model_1h_path=model_1h_path,
            model_13c_path=model_13c_path,
            property_key=property_key,
            multiplicity_mod=multiplicity_mod,
            sim_params=sim_params,
            return_payload=True,
        )

        # Update the cached models for next molecule
        worker_models['bundle'] = models_updated

        return True, {
            'mol_id': mol_id,
            'smiles': smiles,
            'sdf': sdf_path,
            'nb_predicted_H': nH,
            'nb_predicted_C': nC,
            'payload': payload,
            'error': ''
        }, mol_id

    except Exception as e:
        error_msg = str(e)
        # Don't include full traceback in return value to avoid serialization issues
        return False, {
            'mol_id': mol_id,
            'smiles': smiles,
            'sdf': '',
            'nb_predicted_H': 0,
            'nb_predicted_C': 0,
            'payload': None,
            'error': error_msg
        }, mol_id
