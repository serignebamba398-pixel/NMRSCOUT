"""
Hook for Python standard library encodings
CRITICAL for frozen applications
"""
from PyInstaller.utils.hooks import collect_data_files, collect_submodules

datas = collect_data_files('encodings')
hiddenimports = collect_submodules('encodings')
