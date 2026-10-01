"""PyInstaller hook for DGL (Deep Graph Library). Collect binaries and hidden imports; based on PyInstaller issue #6916."""
from PyInstaller.utils.hooks import collect_all, collect_dynamic_libs, copy_metadata
import os
import glob

datas, binaries, hiddenimports = [], [], []

# Collect all DGL data and Python files
tmp = collect_all('dgl', include_py_files=True)
datas, binaries, hiddenimports = tmp[0], tmp[1], tmp[2]

# Collect DGL binaries explicitly
try:
    binaries.extend(collect_dynamic_libs('dgl'))
except Exception:
    pass

# Metadata
try:
    datas.extend(copy_metadata('dgl'))
except Exception:
    pass

# CRITICAL: DGL hidden imports for PyTorch backend
# These are not auto-detected by PyInstaller
hiddenimports.extend([
    'dgl.backend.pytorch',
    'dgl.distributed.nn.pytorch',
    'dgl.distributed.optim.pytorch',
    'dgl.optim.pytorch',
    'dgl._ffi',
    'dgl._ffi.base',
    'dgl.nn.pytorch',
])

# Windows: manually collect DGL DLLs if needed
try:
    import dgl
    dgl_path = os.path.dirname(dgl.__file__)
    for ext in ['*.dll', '*.pyd']:
        for f in glob.glob(os.path.join(dgl_path, '**', ext), recursive=True):
            binaries.append((f, 'dgl'))
except Exception:
    pass

    
    # Collect from dgl root and subdirectories
    for pattern in ['*.dll', '*.pyd', '**/*.dll', '**/*.pyd']:
        for f in glob.glob(os.path.join(dgl_path, pattern), recursive=True):
            if f not in [b[0] for b in binaries]:  # Avoid duplicates
                binaries.append((f, 'dgl'))
except ImportError:
    pass

# Critical DGL imports (dynamically loaded)
critical = [
    'dgl.backend.pytorch',
    'dgl.nn.pytorch',
    'dgl.dataloading',
    'dgl.optim.pytorch',
    'dgl._ffi',
    'dgl._ffi._cy3',
]

for imp in critical:
    if imp not in hiddenimports:
        hiddenimports.append(imp)
