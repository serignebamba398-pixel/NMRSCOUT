"""
Runtime hook for PyTorch CPU on Windows
Sets environment variables and DLL paths for frozen application
Based on community guide for CPU-only PyTorch + DGL
"""
import os
import sys

# Disable CUDA entirely - for CPU-only PyTorch build
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["USE_CUDA"] = "0"
os.environ["DGLBACKEND"] = "pytorch"
os.environ["KMP_DUPLICATE_LIB_OK"] = "True"
os.environ["TORCH_CUDA_ARCH_LIST"] = ""

# Note: No CUDA stub needed when using torch+cpu version
# The CPU-only build of PyTorch already handles this correctly

# Windows 3.8+: Register DLL paths explicitly
if sys.platform == 'win32' and hasattr(os, 'add_dll_directory') and getattr(sys, 'frozen', False):
    base_path = sys._MEIPASS
    dll_paths = [
        base_path,
        os.path.join(base_path, 'torch'),
        os.path.join(base_path, 'torch', 'lib'),
        os.path.join(base_path, 'dgl'),
    ]
    
    for path in dll_paths:
        if os.path.isdir(path):
            try:
                os.add_dll_directory(path)
            except OSError:
                pass
    
    # Also update PATH for compatibility
    os.environ['PATH'] = os.pathsep.join(dll_paths) + os.pathsep + os.environ.get('PATH', '')
