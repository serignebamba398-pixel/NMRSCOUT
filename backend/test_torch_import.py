#!/usr/bin/env python
"""Test minimal torch import"""
import sys
import os
import faulthandler

# Enable faulthandler immediately
faulthandler.enable()
print("[TEST] Faulthandler enabled", flush=True)

# Force CPU mode
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["USE_CUDA"] = "0"
os.environ["DGLBACKEND"] = "pytorch"
os.environ["KMP_DUPLICATE_LIB_OK"] = "True"
print("[TEST] Environment set", flush=True)

# Check if frozen
if getattr(sys, 'frozen', False):
    base = sys._MEIPASS
    print(f"[TEST] Running frozen from: {base}", flush=True)
    
    # Setup DLL paths
    torch_lib = os.path.join(base, 'torch', 'lib')
    if os.path.exists(torch_lib):
        if hasattr(os, 'add_dll_directory'):
            os.add_dll_directory(torch_lib)
            print(f"[TEST] Added DLL dir: {torch_lib}", flush=True)
        os.environ['PATH'] = torch_lib + ';' + os.environ.get('PATH', '')
else:
    print("[TEST] Running from source", flush=True)

# Try imports step by step
print("\n[TEST] Step 1: Import ctypes...", flush=True)
import ctypes
print("[TEST] ctypes OK", flush=True)

print("\n[TEST] Step 2: Import numpy...", flush=True)
import numpy as np
print(f"[TEST] numpy OK: {np.__version__}", flush=True)

print("\n[TEST] Step 3: Try to load torch DLLs manually...", flush=True)
if getattr(sys, 'frozen', False):
    torch_lib = os.path.join(sys._MEIPASS, 'torch', 'lib')
    for dll in ['c10.dll', 'torch_cpu.dll', 'torch_global_deps.dll', 'torch_python.dll']:
        dll_path = os.path.join(torch_lib, dll)
        if os.path.exists(dll_path):
            try:
                handle = ctypes.CDLL(dll_path)
                print(f"[TEST] Loaded {dll}: handle={handle}", flush=True)
            except Exception as e:
                print(f"[TEST] FAILED to load {dll}: {e}", flush=True)

print("\n[TEST] Step 4: Import torch._C (core)...", flush=True)
sys.stdout.flush()
sys.stderr.flush()

try:
    import torch._C
    print(f"[TEST] torch._C OK", flush=True)
except Exception as e:
    print(f"[TEST] torch._C FAILED: {e}", flush=True)
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n[TEST] Step 5: Import torch...", flush=True)
sys.stdout.flush()
sys.stderr.flush()

try:
    import torch
    print(f"[TEST] torch OK: {torch.__version__}", flush=True)
    print(f"[TEST] CUDA available: {torch.cuda.is_available()}", flush=True)
except Exception as e:
    print(f"[TEST] torch FAILED: {e}", flush=True)
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n[TEST] All imports successful!", flush=True)
