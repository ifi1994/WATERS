from rdkit import Chem
import os
import numpy as np
import statistics

# === INPUTS ===
dir_path = r"C:\Users\ifige\Documents\GitHub\SZMAP\szmap_unique_waters"
out_path = r"C:\Users\ifige\Documents\GitHub\SZMAP\in_house_script_results"

round_at = 1
maps_number = 180  # (kept for compatibility but not used for selection anymore)

# === Collect input SDF files ===
sdfs = []
for path in os.listdir(dir_path):
    if os.path.isfile(os.path.join(dir_path, path)):
        sdfs.append(path)

water_mols = []
for file in sdfs:
    path = os.path.join(dir_path, file)
    suppl = Chem.SDMolSupplier(path)
    for mol in suppl:
        if mol is not None:
            water_mols.append(mol)

# === Build grid dictionary ===
grid = {}
for mol in water_mols:
    for i, atom in enumerate(mol.GetAtoms()):
        pos = mol.GetConformer().GetAtomPosition(i)
        x = float(round(pos.x, round_at))
        y = float(round(pos.y, round_at))
        z = float(round(pos.z, round_at))
        key = f"{x}_{y}_{z}"
        try:
            grid[key]['energies'].append(float(mol.GetProp('szmap_neut_diff_free_energy')))
            grid[key]['count'] += 1
        except KeyError:
            grid[key] = {
                'x': x,
                'y': y,
                'z': z,
                'energies': [float(mol.GetProp('szmap_neut_diff_free_energy'))],
                'count': 1
            }

# === Compute mean energies ===
for key in grid:
    grid[key]['mean'] = statistics.mean(grid[key]['energies'])

# === Define thresholds as fractions of number of input files ===
n_files = len(sdfs)
fractions = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.75, 0.90]
thresholds = {int(f * 100): f * n_files for f in fractions}

# === Ensure output directory exists ===
os.makedirs(out_path, exist_ok=True)

# === Write out molecules for each threshold ===
for perc, cutoff in thresholds.items():
    outfile = os.path.join(out_path, f"selected_{perc}.sdf")
    w = Chem.SDWriter(outfile)
    for key, data in grid.items():
        if data['count'] >= cutoff:
            # Create a single-atom molecule
            mol = Chem.RWMol()

            # Assign atom type depending on mean energy
            if data['mean'] < -3:
                atom = Chem.Atom('O')  # Oxygen
            elif -3 <= data['mean'] < 0:
                atom = Chem.Atom('N')  # Nitrogen
            else:
                atom = Chem.Atom('C')  # Carbon

            idx = mol.AddAtom(atom)
            conf = Chem.Conformer(1)
            conf.SetAtomPosition(idx, (data['x'], data['y'], data['z']))
            mol.AddConformer(conf)

            # Attach properties
            mol.SetProp('mean_energy', str(data['mean']))
            mol.SetProp('count', str(data['count']))

            w.write(mol)
    w.close()
    print(f"✅ Wrote {outfile}")
