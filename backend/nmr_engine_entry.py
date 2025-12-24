#!/usr/bin/env python
"""
Main entry point for bundled NMR executable.
Dispatches to different scripts according to command.
"""
import sys
import os
import traceback

# Add scripts directory to path
if getattr(sys, 'frozen', False):
    # Running as compiled
    BASE_DIR = sys._MEIPASS
    os.environ['DGLBACKEND'] = 'pytorch'
    # Default model paths in bundle
    MODELS_DIR = os.path.join(BASE_DIR, 'models')
else:
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    MODELS_DIR = os.path.join(BASE_DIR, 'models')

sys.path.insert(0, os.path.join(BASE_DIR, 'scripts'))
sys.path.insert(0, os.path.join(BASE_DIR, 'ml_dependencies'))

def inject_model_paths():
    """Inject default model paths if not provided and running as frozen."""
    if not getattr(sys, 'frozen', False):
        return
    
    # Check if model paths are already in argv
    has_model_1h = any('--model_1h_path' in arg for arg in sys.argv)
    has_model_13c = any('--model_13c_path' in arg for arg in sys.argv)
    
    if not has_model_1h:
        model_1h = os.path.join(MODELS_DIR, '1H_sparsified_proposed_proposed_set2set_1.pt')
        if os.path.exists(model_1h):
            sys.argv.extend(['--model_1h_path', model_1h])
    
    if not has_model_13c:
        model_13c = os.path.join(MODELS_DIR, '13C_sparsified_proposed_proposed_set2set_1.pt')
        if os.path.exists(model_13c):
            sys.argv.extend(['--model_13c_path', model_13c])

def main():
    if len(sys.argv) < 2:
        print("Usage: nmr_engine <command> [args]")
        print("Commands:")
        print("  test_imports           - Test if all ML imports work")
        print("  integrated_prediction  - Run NMR prediction with ML")
        print("  extract_smiles         - Extract SMILES from SDF")
        print("  predict_1H             - Predict 1H NMR")
        print("  predict_13C            - Predict 13C NMR")
        sys.exit(1)
    
    command = sys.argv[1]
    # Remove command from argv so scripts see correct args
    sys.argv = [sys.argv[0]] + sys.argv[2:]
    
    if command == "test_imports":
        print("Testing critical imports...")
        
        # Test numpy first (simplest)
        print("Step 1: Testing NumPy...")
        try:
            import numpy as np
            print(f"NumPy {np.__version__}")
        except Exception as e:
            print(f"NumPy ERROR: {e}")
            traceback.print_exc()
            sys.exit(1)
            
        # Test PyTorch - with debug
        print("Step 2: Testing PyTorch...", flush=True)
        try:
            print("  2a: Setting DLL directories...", flush=True)
            import os
            import ctypes
            
            base = sys._MEIPASS if getattr(sys, 'frozen', False) else os.path.dirname(__file__)
            torch_lib = os.path.join(base, 'torch', 'lib')
            
            if os.path.exists(torch_lib):
                # Add to DLL directory if Windows 10+
                if hasattr(os, 'add_dll_directory'):
                    os.add_dll_directory(torch_lib)
                    print(f"  Added DLL directory: {torch_lib}", flush=True)
                
                # Also try loading critical DLLs manually
                print("  2b: Pre-loading torch DLLs...", flush=True)
                for dll_name in ['c10.dll', 'torch_cpu.dll', 'torch.dll']:
                    dll_path = os.path.join(torch_lib, dll_name)
                    if os.path.exists(dll_path):
                        try:
                            ctypes.CDLL(dll_path)
                            print(f"    Loaded {dll_name}", flush=True)
                        except Exception as e:
                            print(f"    Failed to load {dll_name}: {e}", flush=True)
            
            print("  2c: Importing torch...", flush=True)
            sys.stdout.flush()
            sys.stderr.flush()
            
            # Use faulthandler to get Python traceback on crash
            import faulthandler
            faulthandler.enable()
            
            import torch
            print(f"  PyTorch {torch.__version__}", flush=True)
        except Exception as e:
            print(f"PyTorch ERROR: {e}", flush=True)
            traceback.print_exc()
            sys.exit(1)
        # Test DGL
        print("Step 3: Testing DGL...", flush=True)
        try:
            import dgl
            print(f"DGL {dgl.__version__}", flush=True)
            g = dgl.graph(([0, 1], [1, 2]))
            print(f"DGL graph test: OK ({g.num_nodes()} nodes)", flush=True)
        except Exception as e:
            print(f"DGL ERROR: {e}", flush=True)
            traceback.print_exc()
        # Test RDKit
        print("Step 4: Testing RDKit...", flush=True)
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles("CCO")
            print(f"RDKit OK: {Chem.MolToSmiles(mol)}", flush=True)
        except Exception as e:
            print(f"RDKit ERROR: {e}", flush=True)
            traceback.print_exc()
            
        print("Import test complete!", flush=True)
        sys.exit(0)
    elif command == "integrated_prediction":
        print("[ENTRY] integrated_prediction command", file=sys.stderr)
        sys.stderr.flush()
        
        try:
            print("[ENTRY] Calling inject_model_paths...", file=sys.stderr)
            sys.stderr.flush()
            inject_model_paths()
            
            print("[ENTRY] Importing integrated_nmr_prediction...", file=sys.stderr)
            sys.stderr.flush()
            
            try:
                from integrated_nmr_prediction import main as run_prediction
                print("[ENTRY] Import successful!", file=sys.stderr)
            except ImportError as ie:
                print(f"[ENTRY] Import ERROR: {ie}", file=sys.stderr)
                traceback.print_exc(file=sys.stderr)
                sys.stderr.flush()
                sys.exit(1)
            
            print("[ENTRY] Running prediction main...", file=sys.stderr)
            sys.stderr.flush()
            run_prediction()
        except Exception as e:
            print(f"[ENTRY] ERROR: {e}", file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
            sys.stderr.flush()
            sys.exit(1)
    elif command == "extract_smiles":
        from extract_smiles import main as run_extract
        run_extract()
    elif command == "predict_1H":
        from predict_1H_real import main as run_1h
        run_1h()
    elif command == "predict_13C":
        from predict_13C_real import main as run_13c
        run_13c()
    elif command == "predict_COSY":
        from predict_COSY_real import main as run_cosy
        run_cosy()
    elif command == "predict_HSQC":
        from predict_HSQC_real import main as run_hsqc
        run_hsqc()
    elif command == "nmr_export_bg":
        # Export spectral background (PNG + JSON) from Bruker pdata
        print("[ENTRY] nmr_export_bg command", file=sys.stderr)
        sys.stderr.flush()
        try:
            # nmr_export_bg.py uses argparse, so we just need to import and run it
            import argparse
            # The script is in the scripts directory which was added to path
            from nmr_export_bg import main as run_export_bg
            run_export_bg()
        except Exception as e:
            print(f"[ENTRY] nmr_export_bg ERROR: {e}", file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
            sys.stderr.flush()
            sys.exit(1)
    else:
        print(f"Unknown command: {command}")
        sys.exit(1)

if __name__ == "__main__":
    main()
