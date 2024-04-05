from rdkit import Chem
import os
import numpy as np

hydrophobic_waters = "magenta > 0"
hydrophilic_waters = "yellow < 0"

hydrophobic_atoms = ['C', 'S', 'F']
# hydrophilic_atoms = ['N', 'O']
hydrophilic_atoms = ['N', 'O', 'Br', 'Cl']

"hydrophobic_atoms <-> hydrophobic_waters > 0"
"hydrophilic_waters <-> hydrophilic_waters < 0"

# folder path
dir_path = "/Users/pantelispanka/Jaqpot/WATERS/data/sdf_multiple_entries"


# water_data = "/Users/pantelispanka/Jaqpot/WATERS/data/Perampanel_analogue_data/Water_data"
water_data = "/Users/pantelispanka/Jaqpot/WATERS/data/Water_data"

# ligand_data = "/Users/pantelispanka/Jaqpot/WATERS/data/Perampanel_analogue_data/Ligand_data"
ligand_data = "/Users/pantelispanka/Jaqpot/WATERS/data/Ligand_data"

# list to store files
sdfs = []

# Iterate directory
for path in os.listdir(water_data):
    # check if current path is a file
    if os.path.isfile(os.path.join(water_data, path)):
        sdfs.append(path)

water_mols = []

for file in sdfs:
    path = water_data + "/" + file
    suppl = Chem.SDMolSupplier(path)
    for mol in suppl:
        water_mols.append(mol)

min_x = None
min_y = None
min_z = None
max_x = None
max_y = None
max_z = None

for mol in water_mols:
    # print(mol.GetProp('szmap_neut_diff_free_energy'))
    for i, atom in enumerate(mol.GetAtoms()):
        positions = mol.GetConformer().GetAtomPosition(i)
        if min_x is None:
            min_x = positions.x
        if min_y is None:
            min_y = positions.y
        if min_z is None:
            min_z = positions.z
        if max_x is None:
            max_x = positions.x
        if max_y is None:
            max_y = positions.y
        if max_z is None:
            max_z = positions.z
        if positions.x < min_x:
            min_x = positions.x
        if positions.y < min_y:
            min_y = positions.y
        if positions.z < min_z:
            min_z = positions.z
        if positions.x > max_x:
            max_x = positions.x
        if positions.y > max_y:
            max_y = positions.y
        if positions.z > max_z:
            max_z = positions.z
        # print(atom.GetSymbol(), positions.x, positions.y, positions.z)
print("START")
print("max_x: " + str(max_x))
print("max_y: " + str(max_y))
print("max_z: " + str(max_z))
print("min_x: " + str(min_x))
print("min_y: " + str(min_y))
print("min_z: " + str(min_z))


if min_x < 0:
    min_x = min_x - 1
if min_y < 0:
    min_y = min_y - 1
if min_z < 0:
    min_z = min_z - 1
if min_x > 0:
    min_x = min_x + 1
if min_y > 0:
    min_y = min_y + 1
if min_z > 0:
    min_z = min_z + 1

if max_x < 0:
    max_x = max_x - 1
if max_y < 0:
    max_y = max_y - 1
if max_z < 0:
    max_z = max_z - 1
if max_x > 0:
    max_x = max_x + 1
if max_y > 0:
    max_y = max_y + 1
if max_z > 0:
    max_z = max_z + 1

print("END")
print("max_x: " + str(max_x))
print("max_y: " + str(max_y))
print("max_z: " + str(max_z))
print("min_x: " + str(min_x))
print("min_y: " + str(min_y))
print("min_z: " + str(min_z))

min_x_int = int(min_x)
min_y_int = int(min_y)
min_z_int = int(min_z)
max_x_int = int(max_x)
max_y_int = int(max_y)
max_z_int = int(max_z)

print("INTS")
print(min_x_int)
print(min_y_int)
print(min_z_int)
print(max_x_int)
print(max_y_int)
print(max_z_int)


grid = {}
for mol in water_mols:
    # print("Free energy: " + str(mol.GetProp('szmap_neut_diff_free_energy')))
    for i, atom in enumerate(mol.GetAtoms()):
        positions = mol.GetConformer().GetAtomPosition(i)
        # print("CENTER: x -> " + str(int(positions.x)), " , y -> " + str(int(positions.y)), " , z -> " + str(int(positions.z)))
        # x = int(round(positions.x, 1))
        # y = int(round(positions.y, 1))
        # z = int(round(positions.z, 1))
        x = round(positions.x, 2)
        y = round(positions.y, 2)
        z = round(positions.z, 2)


        # x = abs(int(positions.x) + min_x_int)
        # y = abs(int(positions.y) + min_y_int)
        # z = abs(int(positions.z) + min_z_int)
        # print("CENTER TRANSLATED: x -> " + str(x), " , y -> " + str(y), " , z -> " +str(z))
        key = str(x) + "_" + str(y) + "_" + str(z)
        try:
            grid[key]['x'] = x
            grid[key]['y'] = y
            grid[key]['z'] = z
            grid[key]['energies'].append(float(mol.GetProp('szmap_neut_diff_free_energy')))
        except KeyError as e:
            grid[key] = {}
            grid[key]['x'] = x
            grid[key]['y'] = y
            grid[key]['z'] = z
            grid[key]['energies'] = []
            grid[key]['energies'].append(float(mol.GetProp('szmap_neut_diff_free_energy')))

import statistics

for key in grid:
    # print(key)
    # print(grid[key]['energies'])
    # print(statistics.mean(grid[key]['energies']))
    grid[key]['mean'] = statistics.mean(grid[key]['energies'])
    try:
        # print(statistics.variance(grid[key]['energies']))
        grid[key]['std'] = statistics.variance(grid[key]['energies'])
    except statistics.StatisticsError as e:
        grid[key]['std'] = 0.0


for key in grid:
    if grid[key]['mean'] > 0.0:
        print(grid[key]['x'])
        print(grid[key]['y'])
        print(grid[key]['z'])
        print(grid[key]['energies'])
        print(grid[key]['std'])

mol_sdfs = []

# Iterate directory
for path in os.listdir(ligand_data):
    # check if current path is a file
    if os.path.isfile(os.path.join(ligand_data, path)):
        mol_sdfs.append(path)


ligand_mols = []

for file in mol_sdfs[:40]:
    path = ligand_data + "/" + file
    suppl = Chem.SDMolSupplier(path)
    for mol in suppl:
        ligand_mols.append(mol)



print(len(ligand_mols))

import math


def distance(x1, y1, z1, x2, y2, z2):
    d = 0.0
    d = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    return d

atoms_set = set()

mol_maps = []
mol_map = {}
for ind, mol in enumerate(ligand_mols):
    # print(mol.GetProp('szmap_neut_diff_free_energy'))
    print(ind)
    print(mol_sdfs[ind])

    mol_map[mol_sdfs[ind]] = {}
    mol_map[mol_sdfs[ind]]['phil_phil'] = {}
    mol_map[mol_sdfs[ind]]['phil_phob'] = {}
    mol_map[mol_sdfs[ind]]['phob_phob'] = {}
    mol_map[mol_sdfs[ind]]['phob_phil'] = {}
    mol_map[mol_sdfs[ind]]['phil_phil']['atoms'] = []
    mol_map[mol_sdfs[ind]]['phil_phil']['grid_energies'] = []
    mol_map[mol_sdfs[ind]]['phil_phil']['std'] = []
    mol_map[mol_sdfs[ind]]['phil_phob']['atoms'] = []
    mol_map[mol_sdfs[ind]]['phil_phob']['grid_energies'] = []
    mol_map[mol_sdfs[ind]]['phil_phob']['std'] = []
    mol_map[mol_sdfs[ind]]['phob_phob']['atoms'] = []
    mol_map[mol_sdfs[ind]]['phob_phob']['grid_energies'] = []
    mol_map[mol_sdfs[ind]]['phob_phob']['std'] = []
    mol_map[mol_sdfs[ind]]['phob_phil']['atoms'] = []
    mol_map[mol_sdfs[ind]]['phob_phil']['grid_energies'] = []
    mol_map[mol_sdfs[ind]]['phob_phil']['std'] = []
    for i, atom in enumerate(mol.GetAtoms()):
        positions = mol.GetConformer().GetAtomPosition(i)
        # print(atom.GetSymbol())
        # print(positions.x)
        # print(positions.y)
        # print(positions.z)
        atom_symbol = atom.GetSymbol()
        atoms_set.add(atom_symbol)
        atom_x = positions.x
        atom_y = positions.y
        atom_z = positions.z
        for key in grid:
            grid_x = grid[key]['x']
            grid_y = grid[key]['y']
            grid_z = grid[key]['z']
            energy = grid[key]['energies']
            std = grid[key]['std']
            dist = distance(atom_x, atom_y, atom_z, grid_x, grid_y, grid_z)
            if dist < 0.8:
                if atom_symbol in hydrophilic_atoms:
                    energy_mean = sum(energy) / len(energy)
                    if -0 > energy_mean > -12:
                        atom_in = f"{atom_symbol}_in_phil"
                        try:
                            mol_map[mol_sdfs[ind]]['phil_phil'][atom_in].append(1)
                        except KeyError as e:
                            mol_map[mol_sdfs[ind]]['phil_phil'][atom_in] = [1]
                        mol_map[mol_sdfs[ind]]['phil_phil']['atoms'].append(atom_symbol)
                        mol_map[mol_sdfs[ind]]['phil_phil']['grid_energies'].append(energy_mean)
                        mol_map[mol_sdfs[ind]]['phil_phil']['std'].append(std)
                    if energy_mean > 0:
                        atom_in = f"{atom_symbol}_in_phob"
                        try:
                            mol_map[mol_sdfs[ind]]['phil_phob'][atom_in].append(1)
                        except KeyError as e:
                            mol_map[mol_sdfs[ind]]['phil_phob'][atom_in] = [1]
                        mol_map[mol_sdfs[ind]]['phil_phob']['atoms'].append(atom_symbol)
                        mol_map[mol_sdfs[ind]]['phil_phob']['grid_energies'].append(energy_mean)
                        mol_map[mol_sdfs[ind]]['phil_phob']['std'].append(std)
                if atom_symbol in hydrophobic_atoms:
                    energy_mean = sum(energy) / len(energy)
                    if -0 > energy_mean > -12:
                        atom_in = f"{atom_symbol}_in_phil"
                        try:
                            mol_map[mol_sdfs[ind]]['phob_phil'][atom_in].append(1)
                        except KeyError as e:
                            mol_map[mol_sdfs[ind]]['phob_phil'][atom_in] = [1]
                        mol_map[mol_sdfs[ind]]['phob_phil']['atoms'].append(atom_symbol)
                        mol_map[mol_sdfs[ind]]['phob_phil']['grid_energies'].append(energy_mean)
                        mol_map[mol_sdfs[ind]]['phob_phil']['std'].append(std)
                    if energy_mean > 0:
                        atom_in = f"{atom_symbol}_in_phob"
                        try:
                            mol_map[mol_sdfs[ind]]['phob_phob'][atom_in].append(1)
                        except KeyError as e:
                            mol_map[mol_sdfs[ind]]['phob_phob'][atom_in] = [1]
                        mol_map[mol_sdfs[ind]]['phob_phob']['atoms'].append(atom_symbol)
                        mol_map[mol_sdfs[ind]]['phob_phob']['grid_energies'].append(energy_mean)
                        mol_map[mol_sdfs[ind]]['phob_phob']['std'].append(std)


                    # print("-----FOUND CLOSE HYDROPHIL----")
                    # print(atom_symbol)
                    # print("----THE ENERGY----")
                    # energy_mean = sum(energy) / len(energy)
                    # print(energy_mean)
            # if grid[key]['mean'] > 0.0:
            #     grid_x = grid[key]['x']
            #     grid_y = grid[key]['y']
            #     grid_z = grid[key]['z']
            #     energy = grid[key]['energies']
            #     std = grid[key]['std']

print(atoms_set)
print("-------MOL MAP-----")
print(len(mol_map))
for mol in mol_map:
    print(mol_map[mol])
    print("\n")



dataset = {}

# dataset["7L10"] = 4.02
dataset["7L11"] = 0.140
dataset["7M8X"] = 0.470
dataset["7M8M"] = 0.12
dataset["7L12"] = 0.128
dataset["7M8Y"] = 0.110
dataset["7M8N"] = 0.1
dataset["7M8O"] = 0.037
dataset["7L13"] = 0.018
dataset["7M8P"] = 0.02
dataset["7M91"] = 0.025
dataset["7L14"] = 0.170
dataset["7M8Z"] = 0.25
dataset["7M90"] = 0.25


dataset_class = {}
# dataset["7L10"] = 0
dataset_class["7L11"] = 0
dataset_class["7M8X"] = 0
dataset_class["7M8M"] = 0
dataset_class["7L12"] = 0
dataset_class["7M8Y"] = 0
dataset_class["7M8N"] = 1
dataset_class["7M8O"] = 1
dataset_class["7L13"] = 1
dataset_class["7M8P"] = 1
dataset_class["7M91"] = 1
dataset_class["7L14"] = 0
dataset_class["7M8Z"] = 0
dataset_class["7M90"] = 0


model_datas = []

print("---PRINTING---")
for key in mol_map:
    model_data = {}
    print("MOLECULE")
    print(key)
    molecule = key.split("-")[1].split(".")[0]
    if molecule != "7L10":
        print(molecule)
        model_data['Id'] = molecule
        print("IC50")
        print(dataset[molecule])
        model_data['IC50'] = dataset[molecule]
        model_data['CLASS'] = dataset_class[molecule]
        print(key.split("-")[1].split(".")[0])
        print(mol_map[key])
        print('phil_phil')
        print(len(mol_map[key]['phil_phil']['atoms']))
        model_data['phil_phil'] = len(mol_map[key]['phil_phil']['atoms'])
        print('phil_phob')
        print(len(mol_map[key]['phil_phob']['atoms']))
        model_data['phil_phob'] = len(mol_map[key]['phil_phob']['atoms'])
        print('phob_phob')
        print(len(mol_map[key]['phob_phob']['atoms']))
        model_data['phob_phob'] = len(mol_map[key]['phob_phob']['atoms'])
        print('phob_phil')
        print(len(mol_map[key]['phob_phil']['atoms']))
        model_data['phob_phil'] = len(mol_map[key]['phob_phil']['atoms'])

        energy = mol_map[key]['phob_phil']['grid_energies']
        std = mol_map[key]['phob_phil']['std']
        energies_phob_phil = mol_map[key]['phob_phil']['grid_energies']



        try:
            energies_phob_phil_mean = sum(energies_phob_phil) / len(energies_phob_phil)
        except ZeroDivisionError as e:
            energies_phob_phil_mean = 0
            print(str(e))
        print("energies_phob_phil")
        print(energies_phob_phil_mean)

        model_data['energies_phob_phil_mean'] = energies_phob_phil_mean

        energy = mol_map[key]['phob_phob']['grid_energies']
        std = mol_map[key]['phob_phob']['std']
        energies_phob_phob = mol_map[key]['phob_phob']['grid_energies']
        try:
            energies_phob_phob_mean = sum(energies_phob_phob) / len(energies_phob_phob)
        except ZeroDivisionError as e:
            energies_phob_phob_mean = 0
            print(str(e))
        print("energies_phob_phob")
        print(energies_phob_phob_mean)

        model_data['energies_phob_phob_mean'] = energies_phob_phob_mean

        energy = mol_map[key]['phil_phil']['grid_energies']
        std = mol_map[key]['phil_phil']['std']
        energies_phil_phil = mol_map[key]['phil_phil']['grid_energies']
        try:
            energies_phil_phil_mean = sum(energies_phil_phil) / len(energies_phil_phil)
        except ZeroDivisionError as e:
            energies_phil_phil_mean = 0
            print(str(e))
        print("energies_phil_phil")
        print(energies_phil_phil_mean)

        model_data['energies_phil_phil_mean'] = energies_phil_phil_mean

        energy = mol_map[key]['phil_phob']['grid_energies']
        std = mol_map[key]['phil_phob']['std']
        energies_phil_phob = mol_map[key]['phil_phob']['grid_energies']
        try:
            energies_phil_phob_mean = sum(energies_phil_phob) / len(energies_phil_phob)
        except ZeroDivisionError  as e:
            energies_phil_phob_mean = 0
            print(str(e))
        print("energies_phil_phil")
        print(energies_phil_phob_mean)

        model_data['energies_phil_phob_mean'] = energies_phil_phob_mean

        for atom in hydrophobic_atoms:
            map = f"{atom}_in_phil"
            map2 = f"{atom}_in_phob"
            try:
                val = len(mol_map[key]['phob_phob'][map])
                model_data[map] = val
            except KeyError as e:
                model_data[map] = 0
            try:
                val = len(mol_map[key]['phob_phil'][map2])
                model_data[map2] = val
            except KeyError as e:
                model_data[map2] = 0
            try:
                val = len(mol_map[key]['phil_phob'][map])
                model_data[map] = val
            except KeyError as e:
                model_data[map] = 0
            try:
                val = len(mol_map[key]['phil_phil'][map2])
                model_data[map2] = val
            except KeyError as e:
                model_data[map2] = 0

        for atom in hydrophilic_atoms:
            map = f"{atom}_in_phil"
            map2 = f"{atom}_in_phob"
            try:
                val = len(mol_map[key]['phil_phob'][map])
                model_data[map] = val
            except KeyError as e:
                model_data[map] = 0
            try:
                val = len(mol_map[key]['phil_phil'][map2])
                model_data[map2] = val
            except KeyError as e:
                model_data[map2] = 0
            try:
                val = len(mol_map[key]['phil_phob'][map])
                model_data[map] = val
            except KeyError as e:
                model_data[map] = 0
            try:
                val = len(mol_map[key]['phil_phil'][map2])
                model_data[map2] = val
            except KeyError as e:
                model_data[map2] = 0




        model_datas.append(model_data)

print(model_datas)

import pandas as pd


df = pd.DataFrame.from_records(model_datas, index="Id")

print(df)

from sklearn.model_selection import train_test_split

train, test = train_test_split(df, test_size=0.2)

print(train)
print(test)

print(list(train))

X_train = train[['phil_phil', 'phil_phob', 'phob_phob', 'phob_phil', 'energies_phob_phil_mean', 'energies_phob_phob_mean', 'energies_phil_phil_mean', 'energies_phil_phob_mean']]
Y_train = train[['IC50']]
Y_train_class = train[['CLASS']]


X_test = test[['phil_phil', 'phil_phob', 'phob_phob', 'phob_phil', 'energies_phob_phil_mean', 'energies_phob_phob_mean', 'energies_phil_phil_mean', 'energies_phil_phob_mean']]
Y_test = test[['IC50']]
Y_test_class = test[['CLASS']]

# X = df[['phil_phil', 'phil_phob', 'phob_phob', 'phob_phil', 'energies_phob_phil_mean', 'energies_phob_phob_mean', 'energies_phil_phil_mean', 'energies_phil_phob_mean', 'C_in_phil', 'C_in_phob', 'S_in_phil', 'S_in_phob', 'F_in_phil', 'F_in_phob', 'N_in_phil', 'N_in_phob', 'O_in_phil', 'O_in_phob', 'Br_in_phil', 'Br_in_phob', 'Cl_in_phil', 'Cl_in_phob']]
X = df[['phil_phil', 'phil_phob', 'phob_phob', 'phob_phil', 'energies_phob_phil_mean', 'energies_phob_phob_mean', 'energies_phil_phil_mean', 'energies_phil_phob_mean']]
y = df['IC50']
from sklearn.linear_model import LinearRegression
from sklearn.svm import LinearSVR
from sklearn.tree import ExtraTreeRegressor
from sklearn.tree import ExtraTreeClassifier
from sklearn.ensemble import RandomForestRegressor

# linear = LinearRegression()
# model = linear.fit(X_train, Y_train)
#
svr = LinearSVR()
# model = svr.fit(X_train, Y_train)

tree = ExtraTreeRegressor()
model = tree.fit(X_train, Y_train)

from sklearn.tree import DecisionTreeRegressor

regressor = DecisionTreeRegressor(random_state=0)

# tree_c = ExtraTreeClassifier()
# model_c = tree_c.fit(X_train, Y_train_class)
from sklearn.model_selection import cross_val_score

print("---CROSS---")
print(cross_val_score(tree, X, y, cv=4))

# rfr = RandomForestRegressor(max_depth=6, random_state=0)
# model = rfr.fit(X_train, Y_train)


preds_train = model.predict(X_train)
print("---TRAIN---")
print(Y_train.values)
print(preds_train)

print("---TEST---")
preds_test = model.predict(X_test)
print(Y_test.values)
print(preds_test)

# preds_train_c = model_c.predict(X_train)
# print("---TRAIN---")
# print(Y_train_class.values)
# print(preds_train_c)

# print("---TEST---")
# preds_test_c = model_c.predict(X_test)
# print(Y_test_class.values)
# print(preds_test_c)


from sklearn.metrics import r2_score

print(r2_score(Y_train, preds_train))
print(r2_score(Y_test, preds_test))



# import matplotlib.pyplot as plt
#
# def explode(data):
#     size = np.array(data.shape)*2
#     data_e = np.zeros(size - 1, dtype=data.dtype)
#     data_e[::2, ::2, ::2] = data
#     return data_e
#
# # build up the numpy logo
# n_voxels = np.zeros((4, 3, 4), dtype=bool)
# n_voxels[0, 0, :] = True
# n_voxels[-1, 0, :] = True
# n_voxels[1, 0, 2] = True
# n_voxels[2, 0, 1] = True
# facecolors = np.where(n_voxels, '#FFD65DC0', '#7A88CCC0')
# edgecolors = np.where(n_voxels, '#BFAB6E', '#7D84A6')
# filled = np.ones(n_voxels.shape)
#
# # upscale the above voxel image, leaving gaps
# filled_2 = explode(filled)
# fcolors_2 = explode(facecolors)
# ecolors_2 = explode(edgecolors)
#
# # Shrink the gaps
# x, y, z = np.indices(np.array(filled_2.shape) + 1).astype(float) // 2
# x[0::2, :, :] += 0.05
# y[:, 0::2, :] += 0.05
# z[:, :, 0::2] += 0.05
# x[1::2, :, :] += 0.95
# y[:, 1::2, :] += 0.95
# z[:, :, 1::2] += 0.95
#
# ax = plt.figure().add_subplot(projection='3d')
# ax.voxels(x, y, z, filled_2, facecolors=fcolors_2, edgecolors=ecolors_2)
# ax.set_aspect('equal')
#
# plt.show()