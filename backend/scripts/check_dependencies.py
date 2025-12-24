#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
check_dependencies.py - Checks and installs ML dependencies
"""

import sys
import subprocess
import os

def check_and_install_dgl():
    try:
        import dgl
        print(f"[OK] DGL {dgl.__version__} already installed")
        return True
    except ImportError:
        print("[X] DGL not found. Installing...")
        
        # Detect PyTorch version
        try:
            import torch
            torch_version = torch.__version__.split('+')[0]
            print(f"PyTorch {torch_version} detected")
        except:
            print("PyTorch not found. Installing...")
            subprocess.check_call([
                sys.executable, "-m", "pip", "install", 
                "torch==2.1.0+cpu", "-f", 
                "https://download.pytorch.org/whl/torch_stable.html"
            ])
        
        # Install DGL
        try:
            if os.name == 'nt':  # Windows
                subprocess.check_call([
                    sys.executable, "-m", "pip", "install", "dgl",
                    "-f", "https://data.dgl.ai/wheels/repo.html"
                ])
            else:  # Linux/Mac
                subprocess.check_call([
                    sys.executable, "-m", "pip", "install", "dgl"
                ])
            
            # Verification
            import dgl
            print(f"[OK] DGL {dgl.__version__} installed successfully")
            return True
            
        except Exception as e:
            print(f"[X] Error installing DGL: {e}")
            return False

def check_all_dependencies():
    """Checks all required dependencies"""
    
    dependencies = {
        'torch': 'torch==2.1.0+cpu -f https://download.pytorch.org/whl/torch_stable.html',
        'numpy': 'numpy',
        'pandas': 'pandas',
        'rdkit': 'rdkit',
    }
    
    missing = []
    
    for module, install_cmd in dependencies.items():
        try:
            __import__(module)
            print(f"[OK] {module} OK")
        except ImportError:
            print(f"[X] {module} missing")
            missing.append(install_cmd)
    
    # Install missing dependencies
    if missing:
        print("\nInstalling missing dependencies...")
        for pkg in missing:
            try:
                cmd = [sys.executable, "-m", "pip", "install"] + pkg.split()
                subprocess.check_call(cmd)
            except Exception as e:
                print(f"Error installing {pkg}: {e}")
    
    # DGL last (requires PyTorch)
    check_and_install_dgl()
    
    # Final verification
    print("\n=== Final verification ===")
    all_ok = True
    for module in ['torch', 'dgl', 'numpy', 'pandas', 'rdkit']:
        try:
            m = __import__(module)
            version = getattr(m, '__version__', 'unknown')
            print(f"[OK] {module} {version}")
        except:
            print(f"[X] {module} ERROR")
            all_ok = False
    
    return all_ok

if __name__ == "__main__":
    success = check_all_dependencies()
    sys.exit(0 if success else 1)