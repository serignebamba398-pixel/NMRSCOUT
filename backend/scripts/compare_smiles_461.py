from rdkit import Chem
from rdkit.Chem import Draw, rdMolDescriptors, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import os

s1 = "C[C@H]1C(=O)N[C@H](C(=O)N1)CC2=C(NC3=C(C=C(C=C23)CC=C(C)C)CC=C(C)C)C(C)(C)C=C"
s2 = "C[C@H]1C(=O)N[C@H](C(=O)N1)CC2=C(NC3=C2C(=C(C=C3)CC=C(C)C)CC=C(C)C)C(C)(C)C=C"

m1 = Chem.MolFromSmiles(s1)
m2 = Chem.MolFromSmiles(s2)

print("DNP000042 :", Chem.MolToSmiles(m1))
print("NPA007877 :", Chem.MolToSmiles(m2))
print()
print("Formule DNP000042 :", rdMolDescriptors.CalcMolFormula(m1))
print("Formule NPA007877 :", rdMolDescriptors.CalcMolFormula(m2))
print()

# Comparer les InChI pour identifier la différence structurale
from rdkit.Chem.inchi import MolToInchi
i1 = MolToInchi(m1)
i2 = MolToInchi(m2)
print("InChI DNP000042 :", i1)
print()
print("InChI NPA007877 :", i2)
print()
print("Identiques ?", i1 == i2)
print()

# Compter les atomes par type
def atom_env(mol, label):
    print(f"\n{label} — environnements des CH2 :")
    mol_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_h, randomSeed=42)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            hcount = atom.GetTotalNumHs()
            if hcount == 2:
                neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
                idx = atom.GetIdx()
                print(f"  C{idx} (CH2) — voisins: {neighbors}")

atom_env(m1, "DNP000042")
atom_env(m2, "NPA007877")
