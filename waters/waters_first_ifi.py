from rdkit import Chem
import os
import regex
import numpy as np


# folder path
#dir_path = "/Users/pantelispanka/Jaqpot/WATERS/data/sdf_multiple_entries"
dir_path = "C:\Users\ifige\Documents\GitHub\SZMAP\szmap_unique_waters"

#dir_path = "C:\\Users\\ifige\\Documents\\GitHub\\WATERS\\data\\estrogen_receptor_water_data\\"
# list to store files
sdfs = []

round_at = 3
maps_number = 180
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
        #print(atom.GetSymbol(), mol.GetProp('szmap_neut_diff_free_energy'), mol.GetProp('s_user_Crystal_Structure_ID'), positions.x, positions.y, positions.z)
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

import numpy
grid = {}
for mol in water_mols:
    #print("Free energy: " + str(mol.GetProp('szmap_neut_diff_free_energy')))
    for i, atom in enumerate(mol.GetAtoms()):
        positions = mol.GetConformer().GetAtomPosition(i)
        #print("CENTER: x -> " + str(int(positions.x)), " , y -> " + str(int(positions.y)), " , z -> " + str(int(positions.z)))
        x = float(round(positions.x, round_at))
        y = float(round(positions.y, round_at))
        z = float(round(positions.z, round_at))

        # x = abs(int(positions.x) + min_x_int)
        # y = abs(int(positions.y) + min_y_int)
        # z = abs(int(positions.z) + min_z_int)
        #print("CENTER TRANSLATED: x -> " + str(x), " , y -> " + str(y), " , z -> " +str(z))
        key = str(x) + "_" + str(y) + "_" + str(z)
        try:
            grid[key]['x'] = x
            grid[key]['y'] = y
            grid[key]['z'] = z
            #grid[key]['IDs'].append(str(mol.GetProp('s_user_Crystal_Structure_ID')))
            grid[key]['energies'].append(float(mol.GetProp('szmap_neut_diff_free_energy')))

        except KeyError as e:
            grid[key] = {}
            grid[key]['x'] = x
            grid[key]['y'] = y
            grid[key]['z'] = z
            #grid[key]['IDs'] = []
            #grid[key]['IDs'].append(str(mol.GetProp('s_user_Crystal_Structure_ID')))
            grid[key]['energies'] = []
            grid[key]['energies'].append(float(mol.GetProp('szmap_neut_diff_free_energy')))
import statistics


for key in grid:
    #print(key)
    # print(grid[key]['energies'])
    # print(statistics.mean(grid[key]['energies']))
    grid[key]['mean'] = statistics.mean(grid[key]['energies'])
   # try:
        # print(statistics.variance(grid[key]['energies']))
        #grid[key]['std'] = statistics.variance(grid[key]['energies'])
   # except statistics.StatisticsError as e:
       # grid[key]['std'] = 0.0


#for key in grid:
 #   if grid[key]['mean'] > -10:
  #      print(grid[key]['x'])
  #      print(grid[key]['y'])
   #     print(grid[key]['z'])
    #    print(grid[key]['IDs'])
     #   print(numpy.unique(grid[key]['IDs']))
      #  print(grid[key]['energies'])
       # print(grid[key]['mean'])


for key in grid:
    if len(numpy.unique(grid[key]['IDs'])) > 2:
        print(grid[key]['x'])
        print(grid[key]['y'])
        print(grid[key]['z'])
        #print(grid[key]['IDs'])
        #print(numpy.unique(grid[key]['IDs']))
        print(grid[key]['energies'])
        print(grid[key]['mean'])

x_1 = []
y_1 = []
z_1 = []
E_1 = []
atom_C = []
l_1 = []
ID_1 = []

x_2 = []
y_2 = []
z_2 = []
E_2 = []
l_2 = []
ID_2 = []
atom_O = []

x_3 = []
y_3 = []
z_3 = []
E_3 = []
l_3 = []
ID_3 = []
atom_N = []

names_C = []
names_O = []
names_N = []

#for key in grid:
#    if len(numpy.unique(grid[key]['IDs'])) > maps_number:
#        if grid[key]['mean'] > 0:
#            x_1.append(grid[key]['x'])
#            y_1.append(grid[key]['y'])
#           z_1.append(grid[key]['z'])
#            E_1.append(grid[key]['mean'])
#            atom_C.append(str('C'))
#            l_1.append(len(numpy.unique(grid[key]['IDs'])))
#            ID_1.append(numpy.unique(grid[key]['IDs']))

for key in grid:
     if len(numpy.unique(grid[key]['IDs'])) > maps_number:
        if grid[key]['mean'] < -3:
            x_2.append(grid[key]['x'])
            y_2.append(grid[key]['y'])
            z_2.append(grid[key]['z'])
            E_2.append(grid[key]['mean'])
            atom_O.append(str('O'))
            l_2.append(len(numpy.unique(grid[key]['IDs'])))
            ID_2.append(numpy.unique(grid[key]['IDs']))

for key in grid:
    if len(numpy.unique(grid[key]['IDs'])) > maps_number:
        if (grid[key]['mean'] > -3) and (grid[key]['mean'] < 0):
            x_3.append(grid[key]['x'])
            y_3.append(grid[key]['y'])
            z_3.append(grid[key]['z'])
            E_3.append(grid[key]['mean'])
            atom_N.append(str('N'))
            l_3.append(len(numpy.unique(grid[key]['IDs'])))
            ID_3.append(numpy.unique(grid[key]['IDs']))

coord = []

#for a1, b1, c1, d1 in zip(atom_C, x_1, y_1, z_1):
#    coord.append([a1, b1, c1, d1])

#for a2, b2, c2, d2 in zip(atom_O, x_1, y_1, z_1):
#    coord.append([a2, b2, c2, d2])

#for a3, b3, c3, d3 in zip(atom_N, x_1, y_1, z_1):
#    coord.append([a3, b3, c3, d3])

#print(coord)

#length = len(l)
#for item in range(1, length):
    #names.append("Oxygen_" + str(item))
#print(names)

#length = len(l)
#print(E)
for item in E_1:
    names_C.append("Carbon_" + str(item))
print(names_C)
print(len(names_C))

for item in E_2:
    names_O.append("Oxygen_" + str(item))
print(names_O)
print(len(names_O))

for item in E_3:
    names_N.append("Nitrogen_" + str(item))
print(names_N)
print(len(names_N))


for name, ha, ia, ja, ka, ea, ids in zip(names_C, atom_C, x_1, y_1, z_1, E_1, ID_1):
    with open("C:\\Users\\ifige\\Documents\\GitHub\\WATERS\\data\\results_ER\\results_maps_number_5_round_at_3\\" + name + ".xyz", 'w') as xyz_file:
                   xyz_file.write("%d\n%s%f%s%s\n" % (1, "Carbon ", ea, " ", ids))
                   xyz_file.write("%s\t%f\t%f\t%f\n" % (ha, ia, ja, ka))

for name, hb, ib, jb, kb, eb, ids in zip(names_O, atom_O, x_2, y_2, z_2, E_2, ID_2):
    with open("C:\\Users\\ifige\\Documents\\GitHub\\WATERS\\data\\results_ER\\results_maps_number_5_round_at_3\\" + name + ".xyz", 'w') as xyz_file:
                   xyz_file.write("%d\n%s%f%s%s\n" % (1, "Oxygen ", eb, " ", ids))
                   xyz_file.write("%s\t%f\t%f\t%f\n" % (hb, ib, jb, kb))

for name, hc, ic, jc, kc, ec, ids in zip(names_N, atom_N, x_3, y_3, z_3, E_3, ID_3):
    with open("C:\\Users\\ifige\\Documents\\GitHub\\WATERS\\data\\results_ER\\results_maps_number_5_round_at_3\\" + name + ".xyz", 'w') as xyz_file:
                  xyz_file.write("%d\n%s%f%s%s\n" % (1, "Nitrogen ", ec, " ", ids))
                  xyz_file.write("%s\t%f\t%f\t%f\n" % (hc, ic, jc, kc))


#print(ID)
#writer = Chem.SDWriter('C:\\Users\\ifige\\Documents\\GitHub\\WATERS\\data\\results.sdf')



#for key in grid:
    #if len(numpy.unique(grid[key]['IDs'])) > 6:
                #writer.write(mol, )
                #positions = mol.GetConformer().GetAtomPosition(i)
                #x = int(grid[key]['x'])
                #y = int(grid[key]['y'])
                #z = int(grid[key]['z'])


    #print(grid[key]['IDs'])
        #print(numpy.unique(grid[key]['IDs']))
        #print(grid[key]['energies'])
        #print(grid[key]['mean'])

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



