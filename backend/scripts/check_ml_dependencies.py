#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
check_ml_dependencies.py - Check ML dependencies
"""

import sys
import subprocess
import os
import platform

def check_module(module_name):
    """Check whether a module is installed."""
    try:
        __import__(module_name)
        return True
    except ImportError:
        return False

def check_all_dependencies():
    """Check all dependencies."""
    
    print("=== Checking ML dependencies ===")
    
    # 1. PyTorch
    if not check_module("torch"):
        print("PyTorch missing")
        return False
    else:
        import torch
        print(f"PyTorch {torch.__version__}")
    
    # 2. DGL
    if not check_module("dgl"):
        print("DGL missing")
        return False
    else:
        import dgl
        print(f"DGL {dgl.__version__}")
    
    # 3. Other dependencies
    other_deps = ["numpy", "pandas", "rdkit"]
    
    for module in other_deps:
        if not check_module(module):
            print(f"{module} missing")
            return False
        else:
            m = __import__(module)
            version = getattr(m, "__version__", "?")
            print(f"{module} {version}")
    
    # 4. Test final DGL
    print("\n=== DGL test ===")
    try:
        import torch
        import dgl
        import numpy as np
        
        # DGL graph creation test
        g = dgl.graph(([0, 1], [1, 2]))
        print(f"DGL test: graph created with {g.num_nodes()} nodes")
        
        # Test PyTorch
        x = torch.randn(3, 3)
        print(f"PyTorch test: tensor created {x.shape}")
        
        print("\nAll dependencies are working!")
        return True
        
    except Exception as e:
        print(f"\nError during test: {e}")
        return False

if __name__ == "__main__":
    success = check_all_dependencies()
    sys.exit(0 if success else 1)
