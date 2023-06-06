from rdkit import Chem
import os

# folder path
dir_path = "C:\\Users\\ifige\\Documents\\GitHub\\WATERS\\data\\molecules"

# list to store files
molecules = []

# Iterate directory
for path in os.listdir(dir_path):
    # check if current path is a file
    if os.path.isfile(os.path.join(dir_path, path)):
        molecules.append(path)

compound_lst = []

for file in molecules:
    path = dir_path + "/" + file
    suppl = Chem.SDMolSupplier(path, removeHs=False)
    for mol in suppl:
        compound_lst.append(mol)

for mol in compound_lst:
    print(Chem.MolToMolBlock(mol))

# Alternative
for mol_2 in compound_lst:
    print(mol_2.GetProp('s_user_Name'))
    for i, atom in enumerate(mol.GetAtoms()):
        positions = mol.GetConformer().GetAtomPosition(i)
        print(atom.GetSymbol(), positions.x, positions.y, positions.z)


