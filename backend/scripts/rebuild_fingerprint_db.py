#!/usr/bin/env python
"""Rebuild the fingerprint database."""

import os
import sys
import multiprocessing

def main():
    # Suppress RDKit warnings
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')

    # Add script directory to path
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

    from fingerprint_db_predictor import FingerprintPredictor

    print("Starting fingerprint database rebuild...")
    print("=" * 60)

    predictor = FingerprintPredictor()

    print("=" * 60)
    print("Rebuild complete!")


if __name__ == '__main__':
    # Required for Windows multiprocessing
    multiprocessing.freeze_support()
    main()
