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
    import nmr_engine_entry
    nmr_engine_entry.main()
