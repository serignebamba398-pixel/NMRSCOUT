# -*- mode: python ; coding: utf-8 -*-
# MINIMAL test spec - just torch import test

import sys
sys.setrecursionlimit(10000)

import os
from pathlib import Path
from PyInstaller.utils.hooks import collect_dynamic_libs

block_cipher = None
BACKEND_DIR = Path(r"C:\Users\seck-s\Documents\_UNC\VERSIONDOCKER00Last\annotation-nmr\backend")

# Force CPU mode  
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["DGLBACKEND"] = "pytorch"

binaries = []

# Collect ONLY dynamic libs for torch
try:
    binaries.extend(collect_dynamic_libs('torch'))
    print(f"Collected torch binaries: {len([b for b in binaries])}")
except Exception as e:
    print(f"Failed to collect torch libs: {e}")

# Minimal hidden imports
hiddenimports = [
    'torch',
    'torch._C',
    'torch.cuda',
    'numpy',
]

# Exclude CUDA libs
excludes = [
    'torch.backends.cudnn',
    'torch.distributed',
    'nvidia', 
    'triton',
    'matplotlib',
    'tkinter',
]

a = Analysis(
    [str(BACKEND_DIR / "test_torch_import.py")],
    pathex=[str(BACKEND_DIR)],
    binaries=binaries,
    datas=[],
    hiddenimports=hiddenimports,
    excludes=excludes,
    runtime_hooks=[str(BACKEND_DIR / "hooks" / "runtimehook-torch.py")],
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="test_torch", 
          strip=False, upx=False, console=True)
coll = COLLECT(exe, a.binaries, a.zipfiles, a.datas, strip=False, upx=False, name="test_torch")
