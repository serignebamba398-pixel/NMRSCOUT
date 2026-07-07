#!/usr/bin/env python
"""
Simple launcher for nmr_engine bundle
"""
import sys
import os

# Force CPU mode
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
os.environ["DGLBACKEND"] = "pytorch"
os.environ["KMP_DUPLICATE_LIB_OK"] = "True"

# Add DLL directories BEFORE any imports (critical for Windows PyTorch)
if getattr(sys, 'frozen', False) and sys.platform == 'win32':
    base_dir = sys._MEIPASS
    dll_dirs = [
        os.path.join(base_dir, 'torch', 'lib'),
        os.path.join(base_dir, 'dgl'),
        base_dir,
    ]
    for dll_dir in dll_dirs:
        if os.path.isdir(dll_dir) and hasattr(os, 'add_dll_directory'):
            try:
                os.add_dll_directory(dll_dir)
            except Exception:
                pass
    # Also add to PATH
    os.environ['PATH'] = ';'.join(dll_dirs) + ';' + os.environ.get('PATH', '')

# Run main directly
if __name__ == "__main__":
    # Check if there's an external override for nmr_engine_entry
    # This allows hot-patching without recompiling the bundle
    override_path = os.environ.get('NMR_SCRIPTS_OVERRIDE', '')
    override_entry = os.path.join(override_path, 'nmr_engine_entry.py') if override_path else ''
    
    if override_path and os.path.isfile(override_entry):
        print(f"[launcher] Loading nmr_engine_entry from override: {override_entry}", file=sys.stderr)
        # Insert override path at beginning so imports find corrected scripts
        sys.path.insert(0, override_path)
        # Load from override
        import importlib.util
        spec = importlib.util.spec_from_file_location("nmr_engine_entry", override_entry)
        nmr_engine_entry = importlib.util.module_from_spec(spec)
        sys.modules['nmr_engine_entry'] = nmr_engine_entry
        spec.loader.exec_module(nmr_engine_entry)
    else:
        import nmr_engine_entry
    
    nmr_engine_entry.main()
