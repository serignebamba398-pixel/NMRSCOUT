"""
PyInstaller hook for torch (production-ready)
Windows: manually collect all DLLs from torch/lib
Based on official PyTorch packaging guide
"""
from PyInstaller.utils.hooks import (
    collect_data_files, 
    collect_dynamic_libs, 
    collect_submodules,
    copy_metadata
)
import os
import glob

# Standard collection with source files (required for JIT)
datas = collect_data_files('torch', include_py_files=True)
datas += copy_metadata('torch')
binaries = collect_dynamic_libs('torch')
hiddenimports = collect_submodules('torch')

# CRITICAL Windows: manually collect EVERY DLL from torch/lib
# PyInstaller may miss dynamically loaded native libraries
try:
    import torch
    torch_lib = os.path.join(os.path.dirname(torch.__file__), 'lib')
    if os.path.isdir(torch_lib):
        # Collect .dll files
        for dll in glob.glob(os.path.join(torch_lib, '*.dll')):
            binaries.append((dll, 'torch/lib'))
        # Collect .pyd files (Python extensions)
        for pyd in glob.glob(os.path.join(torch_lib, '*.pyd')):
            binaries.append((pyd, 'torch/lib'))
except ImportError:
    pass

# Critical submodules not auto-detected
hiddenimports += [
    'torch._C',
    'torch.nn.modules.lazy',
    'torch.backends.cudnn',
    'torch.utils._device_count',
]
