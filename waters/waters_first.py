from rdkit import Chem
import os

# folder path
dir_path = "/Users/pantelispanka/Jaqpot/WATERS/data/sdf_multiple_entries"

# list to store files
sdfs = []

# Iterate directory
for path in os.listdir(dir_path):
    # check if current path is a file
    if os.path.isfile(os.path.join(dir_path, path)):
        sdfs.append(path)

water_mols = []

for file in sdfs:
    path = dir_path + "/" + file
    suppl = Chem.SDMolSupplier(path)
    for mol in suppl:
        water_mols.append(mol)

for mol in water_mols:
    print(mol.GetProp('szmap_neut_diff_free_energy'))
    for i, atom in enumerate(mol.GetAtoms()):
        positions = mol.GetConformer().GetAtomPosition(i)
        print(atom.GetSymbol(), positions.x, positions.y, positions.z)