import sys
import json
from rdkit import Chem

def main():
    if len(sys.argv) < 2:
        print(json.dumps({"ok": False, "error": "Missing SDF path argument"}))
        sys.exit(1)

    sdf_path = sys.argv[1]
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False, sanitize=False)
    
    try:
        m = suppl[0]
    except Exception:
        m = None

    if m is None:
        print(json.dumps({"ok": False, "error": "SDF illisible"}))
        sys.exit(2)

    try:
        mol = Chem.RemoveHs(m)
    except:
        mol = m
        
    smi = Chem.MolToSmiles(mol)
    print(json.dumps({"ok": True, "smiles": smi}))

if __name__ == "__main__":
    main()
