# -*- mode: python ; coding: utf-8 -*-
# MINIMAL spec - avoid collect_all completely
# Only add what's strictly necessary

import sys
sys.setrecursionlimit(10000)

import os
from pathlib import Path
from PyInstaller.utils.hooks import collect_dynamic_libs, collect_data_files, collect_submodules, collect_all

block_cipher = None
APP_NAME = "nmr_engine_production"
try:
    BACKEND_DIR = Path(__file__).resolve().parent
except NameError:
    # PyInstaller executes spec via exec() and may not define __file__.
    # SPECPATH is usually injected by PyInstaller; fall back to current working directory.
    BACKEND_DIR = Path(globals().get('SPECPATH', os.getcwd())).resolve()

# Force CPU mode  
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["DGLBACKEND"] = "pytorch"

datas = []
binaries = []
hiddenimports = []

# Custom application data ONLY
for subdir in ["scripts", "ml_dependencies", "models"]:
    p = BACKEND_DIR / subdir
    if p.exists():
        datas.append((str(p), subdir))

# Collect NumPy (required by runtime + test_imports)
try:
    import numpy as _np
    from PyInstaller.utils.hooks import collect_dynamic_libs, collect_data_files, collect_submodules

    # Include numpy package directory explicitly (robust even when import graph misses it)
    _numpy_dir = Path(_np.__file__).resolve().parent
    datas.append((str(_numpy_dir), 'numpy'))

    # Wheels often ship extra DLLs in a sibling folder (numpy.libs)
    _numpy_libs = _numpy_dir.parent / 'numpy.libs'
    if _numpy_libs.exists():
        datas.append((str(_numpy_libs), 'numpy.libs'))

    # Also collect using hooks to capture dynamic libs + hidden imports
    datas.extend(collect_data_files('numpy'))
    binaries.extend(collect_dynamic_libs('numpy'))
    hiddenimports.extend(collect_submodules('numpy'))

    print("Collected numpy (package dir + hooks)")
except Exception as e:
    print(f"Warning: numpy collection failed: {e}")

# Collect torch COMPLETELY using collect_all
try:
    tmp_ret = collect_all('torch')
    datas.extend(tmp_ret[0])
    binaries.extend(tmp_ret[1])
    hiddenimports.extend(tmp_ret[2])
    print(f"Collected torch via collect_all")
except Exception as e:
    print(f"Warning: torch collection failed: {e}")

# Collect DGL (add to binaries before filtering)
try:
    dgl_libs = collect_dynamic_libs('dgl')
    binaries.extend(dgl_libs)
    hiddenimports.extend(collect_submodules('dgl'))
    print(f"Collected {len(dgl_libs)} DGL binaries")
except Exception as e:
    print(f"Warning: DGL collection failed: {e}")

# Collect RDKit (add to binaries before filtering)
try:
    rdkit_libs = collect_dynamic_libs('rdkit')
    binaries.extend(rdkit_libs)
    print(f"Collected {len(rdkit_libs)} RDKit binaries")
except Exception as e:
    print(f"Warning: RDKit collection failed: {e}")

# Filter binaries to remove problematic MKL/OpenMP libraries
# AND normalize to TOC format (dest_name, src_path, 'BINARY')
filtered_binaries = []
for item in binaries:
    if len(item) >= 3:
        dest, src, kind = item[:3]
    elif len(item) == 2:
        # Handle (src, dest_dir) from collect_dynamic_libs
        val1, val2 = item
        if os.path.isabs(val1):
            src = val1
            dest_dir = val2
            dest = os.path.join(dest_dir, os.path.basename(src))
        else:
            dest = val1
            src = val2
        kind = 'BINARY'
    else:
        continue
        
    lower_src = src.lower()
    if 'mkl_blacs' in lower_src or 'mkl_pgi' in lower_src or 'omptarget' in lower_src:
        print(f"Excluding binary: {src}")
        continue
    filtered_binaries.append((dest, src, kind))

binaries = filtered_binaries

# Manually add missing MKL DLLs (version 1) from Conda environment
# PyTorch may need version 1 DLLs, not version 2
try:
    # Try multiple potential paths for MKL v1 DLLs
    potential_paths = [
        os.path.join(sys.prefix, 'Library', 'bin'),
    ]
    
    mkl_v1_dlls = [
        'mkl_core.1.dll',
        'mkl_intel_thread.1.dll',
        'mkl_avx512.1.dll',
        'mkl_avx2.1.dll',
        'mkl_avx.1.dll',
        'mkl_def.1.dll',
        'mkl_rt.1.dll',
        'mkl_sequential.1.dll',
        'mkl_mc.1.dll',
        'mkl_mc3.1.dll',
        'mkl_vml_avx512.1.dll',
        'mkl_vml_avx2.1.dll',
        'mkl_vml_avx.1.dll',
        'mkl_vml_def.1.dll',
        'mkl_vml_cmpt.1.dll',
        'mkl_vml_mc.1.dll',
        'mkl_vml_mc3.1.dll',
    ]
    
    for conda_bin in potential_paths:
        if os.path.exists(conda_bin):
            for dll in mkl_v1_dlls:
                dll_path = os.path.join(conda_bin, dll)
                if os.path.exists(dll_path):
                    print(f"Adding missing DLL: {dll_path}")
                    binaries.append((dll, dll_path, 'BINARY'))
            break  # Only use the first valid path
except Exception as e:
    print(f"Warning: Failed to add extra MKL DLLs: {e}")

# Minimal hidden imports
hiddenimports.extend([
    'dgl.backend.pytorch',
    'dgl.nn.pytorch',
    'dgl.optim.pytorch',
    'dgl._ffi',
    'dgl._ffi.base',
    'torch._C',
    'torch.cuda',
    'torch.distributed',
    'torch.distributed.rpc',
    'torch.distributed.rpc.api',
    'torch.distributed.rpc.backend_registry',
    'numpy',
    'scipy',
    'pandas',
    'sklearn',
    'rdkit',
    'rdkit.Chem',
    # For nmr_export_bg.py spectral background export
    'nmrglue',
    'nmrglue.bruker',
    'nmrglue.fileio',
    'nmrglue.fileio.bruker',
    'matplotlib',
    'matplotlib.pyplot',
    'matplotlib.backends',
    'matplotlib.backends.backend_agg',
])

excludes = [
    'torch.backends.cudnn',
    'torchvision',
    'nvidia', 
    'triton',
    'tkinter',
    'jupyter',
    'jupyterlab',
    'notebook',
    'IPython',
    'ipykernel',
]

# Convert all datas to TOC format (dest, src, type) manually
final_datas = []
for item in datas:
    if len(item) == 2:
        src, dest = item
        if os.path.isfile(src):
            # Fix: dest is a directory, append filename to get full target path
            target_dest = os.path.join(dest, os.path.basename(src))
            final_datas.append((target_dest, src, 'DATA'))
        elif os.path.isdir(src):
            for root, dirs, files in os.walk(src):
                for f in files:
                    abs_file = os.path.join(root, f)
                    rel_path = os.path.relpath(abs_file, src)
                    target_dest = os.path.join(dest, rel_path)
                    final_datas.append((target_dest, abs_file, 'DATA'))
    elif len(item) == 3:
        final_datas.append(item)

a = Analysis(
    [str(BACKEND_DIR / "launcher.py")],
    pathex=[str(BACKEND_DIR)],
    binaries=[], 
    datas=[],    
    hiddenimports=hiddenimports,
    hookspath=[str(BACKEND_DIR / "hooks")],
    excludes=excludes,
    runtime_hooks=[str(BACKEND_DIR / "hooks" / "runtimehook-torch.py")],
    noarchive=False,
)

a.binaries.extend(binaries)
a.datas.extend(final_datas)

# Filter out datas with paths > 200 chars
filtered_datas = []
for dest, src, dtype in a.datas:
    full_path = os.path.join('dist', APP_NAME, '_internal', dest)
    if len(full_path) < 200 and 'jupyter' not in dest.lower() and 'labextension' not in dest.lower():
        filtered_datas.append((dest, src, dtype))
a.datas = filtered_datas

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name=APP_NAME, 
          strip=False, upx=False, console=True)
coll = COLLECT(exe, a.binaries, a.zipfiles, a.datas, strip=False, upx=False, name=APP_NAME)
