# -*- mode: python ; coding: utf-8 -*-

block_cipher = None

import sys
import os

# Add the scripts directory to the path so analysis can find modules like ml_sdf_utils
sys.path.append(os.path.abspath('scripts'))
sys.path.append(os.path.abspath('ml_dependencies'))

a = Analysis(
    ['py_entry.py'],
    pathex=[],
    binaries=[],
    datas=[
        ('scripts/*.py', 'scripts'),
        ('ml_dependencies/*.py', 'ml_dependencies'),
        ('ml_dependencies/data', 'ml_dependencies/data'),
        ('ml_dependencies/gnn_models', 'ml_dependencies/gnn_models'),
        ('models/*.pt', 'models'),  # Include the pytorch models
    ],
    hiddenimports=[
        'rdkit', 'rdkit.Chem', 'rdkit.Chem.AllChem', 
        'torch', 'dgl', 'dgllife', 'pandas', 'numpy', 'scipy',
        'ml_sdf_utils', # Ensure this is picked up
        'sklearn.metrics', 'sklearn.ensemble'
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='nmr_engine',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='nmr_engine',
)
