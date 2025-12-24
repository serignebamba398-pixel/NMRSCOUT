"""
Custom hook for PyInstaller to collect torch/lib
Ensures torch native libraries are bundled
"""
from PyInstaller.utils.hooks import collect_all, get_package_paths
from pathlib import Path
import os

def hook(mod):
    """Hook to collect torch/lib"""
    try:
        torch_base, torch_pkg = get_package_paths("torch")
        torch_lib = Path(torch_pkg) / "lib"
        
        if not torch_lib.exists():
            return {}
        
        # Collect torch/lib files
        datas = []
        binaries = []
        
        # Add torch/lib as a data dir to preserve structure
        datas.append((str(torch_lib), "torch/lib"))
        
        # Also add each DLL individually to binaries
        for dll in torch_lib.glob("*.dll"):
            binaries.append((str(dll), "torch/lib"))
        
        return {"datas": datas, "binaries": binaries}
    except:
        return {}
