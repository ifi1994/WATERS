import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from minisom import MiniSom

# --- Configuration ---
# You can change these paths if needed
INPUT_DIR = r"C:\Users\ifige\Documents\GitHub\SZMAP\szmap_unique_waters"
OUTPUT_DIR = r"C:\Users\ifige\Documents\GitHub\SZMAP\Kohonen_results"
MIN_CLUSTER_SIZE = 250


# --- Main Functions ---

def parse_sdf_files(input_dir):
    """
    Parses all SDF files in the specified directory to extract oxygen atom coordinates,
    szmap_neut_diff_free_energy, and szmap_probability.
    Returns a list of dictionaries with data for each water molecule.
    """
    print(f"Parsing SDF files from: {input_dir}")
    water_molecules = []

    # Use glob to find all files ending with .sdf
    file_paths = glob.glob(os.path.join(input_dir, "*.sdf"))
    if not file_paths:
        print(f"No .sdf files found in {input_dir}. Please check the path and file extensions.")
        return []

    for file_path in file_paths:
        with open(file_path, 'r') as f:
            content = f.read()

        # Split the content into individual water molecule blocks
        water_blocks = content.split('$$$$\n')

        for block in water_blocks:
            if block.strip() == "":
                continue

            lines = block.split('\n')
            mol_name = lines[0].strip()

            # Find the line containing the oxygen atom's coordinates
            oxygen_coords = None
            for line in lines[4:]:
                parts = line.split()
                if len(parts) >= 4 and parts[3].strip() == 'O':
                    try:
                        x = float(parts[0])
                        y = float(parts[1])
                        z = float(parts[2])
                        oxygen_coords = (x, y, z)
                        break
                    except (ValueError, IndexError):
                        continue

            if oxygen_coords:
                # Find the szmap_neut_diff_free_energy and szmap_probability
                free_energy = None
                probability = None
                for line in lines:
                    if '<szmap_neut_diff_free_energy>' in line:
                        try:
                            free_energy = float(lines[lines.index(line) + 1].strip())
                        except (ValueError, IndexError):
                            pass
                    if '<szmap_probability>' in line:
                        try:
                            probability = float(lines[lines.index(line) + 1].strip())
                        except (ValueError, IndexError):
                            pass

                water_molecules.append({
                    "id": f"{os.path.basename(file_path)}--{mol_name}",
                    "coords": oxygen_coords,
                    "block": block + '$$$$\n',
                    "free_energy": free_energy,
                    "probability": probability
                })

    print(f"Successfully parsed {len(water_molecules)} water molecules.")
    return water_molecules


def run_som_clustering(data, som_x=10, som_y=10, som_sigma=1.0, som_learning_rate=0.5):
    """
    Initializes and trains a Self-Organizing Map (SOM) on the input data,
    which now includes coordinates and free energy.
    Returns the trained SOM model and the best matching unit (BMU) for each data point.
    """
    print("Initializing and training the Self-Organizing Map...")

    # SOM is 2D and now needs to be trained on 4D data (x, y, z, free_energy)
    # The input dimension of the SOM is now 4
    som = MiniSom(som_x, som_y, 4, sigma=som_sigma, learning_rate=som_learning_rate, activation_distance='euclidean',
                  neighborhood_function='gaussian')
    som.random_weights_init(data)

    # Train the SOM
    som.train_random(data, 1000)

    # Get the best matching unit (BMU) for each data point
    cluster_assignments = np.array([som.winner(x) for x in data])

    print("SOM training complete.")
    return som, cluster_assignments


def filter_clusters(water_molecules, som, cluster_assignments, min_size):
    """
    Filters out clusters that contain fewer than the specified minimum number of molecules.
    Returns the filtered list of molecules and their corresponding cluster assignments.
    """
    print(f"Filtering clusters with fewer than {min_size} molecules...")

    som_x, som_y, _ = som.get_weights().shape
    cluster_labels = np.array([(y * som_x + x) for (x, y) in cluster_assignments])
    unique_clusters, counts = np.unique(cluster_labels, return_counts=True)

    valid_cluster_ids = unique_clusters[counts >= min_size]

    filtered_indices = []
    for cluster_id in valid_cluster_ids:
        indices_for_cluster = np.where(cluster_labels == cluster_id)[0]
        filtered_indices.extend(indices_for_cluster)

    filtered_molecules = [water_molecules[i] for i in filtered_indices]
    filtered_assignments = cluster_assignments[filtered_indices]

    print(f"Retained {len(valid_cluster_ids)} clusters with {len(filtered_molecules)} total molecules.")

    return filtered_molecules, filtered_assignments, valid_cluster_ids


def select_key_waters(water_molecules, cluster_assignments):
    """
    Selects a single representative water molecule for each cluster based on
    free energy and probability.
    """
    print("Selecting key waters for each cluster...")

    som_x, som_y, _ = som.get_weights().shape
    cluster_labels = np.array([(y * som_x + x) for (x, y) in cluster_assignments])
    unique_clusters = np.unique(cluster_labels)

    key_waters = []

    for cluster_id in unique_clusters:
        # Get molecules and their data for the current cluster
        cluster_indices = np.where(cluster_labels == cluster_id)[0]
        cluster_mols = [water_molecules[i] for i in cluster_indices]

        # Use a mask to filter out None values and create a list of valid molecules
        valid_mols = [mol for mol in cluster_mols if mol['free_energy'] is not None and mol['probability'] is not None]

        if not valid_mols:
            print(f"Cluster {cluster_id} has no valid probability or energy data. Skipping key water selection.")
            continue

        # Find the highest probability in the cluster
        max_prob = max([mol['probability'] for mol in valid_mols])

        # Find all molecules with this highest probability
        max_prob_mols = [mol for mol in valid_mols if mol['probability'] == max_prob]

        # Calculate mean energy to determine selection logic
        mean_energy = np.mean([mol['free_energy'] for mol in valid_mols])

        selected_mol = None
        if mean_energy >= 0:
            # Positive mean energy: find highest energy among max probability molecules
            selected_mol = max(max_prob_mols, key=lambda mol: mol['free_energy'])
        else:
            # Negative mean energy: find lowest energy among max probability molecules
            selected_mol = min(max_prob_mols, key=lambda mol: mol['free_energy'])

        if selected_mol:
            key_waters.append(selected_mol)

    print(f"Selected {len(key_waters)} key waters.")
    return key_waters


def generate_outputs(water_molecules, som, cluster_assignments, prefix=""):
    """
    Generates all requested output files (SDF, CSV, PNG).
    """
    print(f"Generating outputs in directory: {OUTPUT_DIR}")

    # Create the output directory if it doesn't exist
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        print(f"Created directory: {OUTPUT_DIR}")

    # --- Prepare Data for Output ---
    data_points = [mol['coords'] for mol in water_molecules]
    som_x, som_y, _ = som.get_weights().shape
    cluster_labels = [(y * som_x + x) for (x, y) in cluster_assignments]

    # --- Generate CSV File ---
    csv_path = os.path.join(OUTPUT_DIR, f"kohonen_clustering_results{prefix}.csv")
    df = pd.DataFrame({
        "molecule_id": [mol['id'] for mol in water_molecules],
        "x": [p[0] for p in data_points],
        "y": [p[1] for p in data_points],
        "z": [p[2] for p in data_points],
        "cluster_id": cluster_labels
    })
    df.to_csv(csv_path, index=False)
    print(f"CSV file saved to: {csv_path}")

    # --- Generate Clustered SDF File ---
    sdf_path = os.path.join(OUTPUT_DIR, f"clustered_waters{prefix}.sdf")
    with open(sdf_path, 'w') as f:
        for i, mol in enumerate(water_molecules):
            block = mol['block']
            v2000_end_line = block.find("M  END") + len("M  END")
            modified_block = (block[:v2000_end_line] +
                              f'\n> <kohonen_cluster_id>\n{cluster_labels[i]}' +
                              block[v2000_end_line:])
            f.write(modified_block)

    print(f"SDF file saved to: {sdf_path}")

    # --- Generate PNG Plot (U-Matrix) ---
    plt.figure(figsize=(10, 10))
    distance_map = som.distance_map()
    plt.pcolor(distance_map.T, cmap='bone_r')
    plt.colorbar()
    markers = ['o', 's', 'D', 'v', '^', '<', '>', 'p', 'h']
    colors = plt.cm.tab10(range(10))
    for i, coords in enumerate(data_points):
        win_x, win_y = cluster_assignments[i]
        cluster_id = cluster_labels[i] % len(colors)
        plt.plot(win_x + 0.5, win_y + 0.5,
                 marker=markers[cluster_id % len(markers)],
                 markeredgecolor=colors[cluster_id],
                 markerfacecolor='None',
                 markersize=12,
                 markeredgewidth=2)

    plt.title('Kohonen Self-Organizing Map Clustering (U-Matrix)')
    plt.axis([0, som_x, 0, som_y])
    png_path = os.path.join(OUTPUT_DIR, f"kohonen_cluster_map{prefix}.png")
    plt.savefig(png_path)
    plt.close()

    print(f"PNG plot saved to: {png_path}")

    # --- Generate 3D PNG Plot of the Clusters ---
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    data_points_np = np.array(data_points)
    cluster_labels_np = np.array(cluster_labels)
    unique_clusters = np.unique(cluster_labels_np)
    colors = plt.cm.nipy_spectral(np.linspace(0, 1, len(unique_clusters)))

    for i, cluster_id in enumerate(unique_clusters):
        cluster_points = data_points_np[cluster_labels_np == cluster_id]
        ax.scatter(cluster_points[:, 0], cluster_points[:, 1], cluster_points[:, 2],
                   c=[colors[i]], marker='o')

    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')
    ax.set_title('3D Visualization of Water Molecule Clusters')
    png_3d_path = os.path.join(OUTPUT_DIR, f"3d_water_clusters{prefix}.png")
    plt.savefig(png_3d_path)
    plt.close()

    print(f"3D PNG plot saved to: {png_3d_path}")
    print("All output files have been generated.")


def generate_centroids_sdf(water_molecules, som, cluster_assignments, prefix=""):
    """
    Calculates the centroid of each cluster and exports their coordinates and mean free energy
    to a new SDF file.
    """
    print("Calculating cluster centroids and generating centroids SDF file...")
    data_points = np.array([mol['coords'] for mol in water_molecules])
    free_energies = np.array([mol['free_energy'] for mol in water_molecules])
    som_x, som_y, _ = som.get_weights().shape
    cluster_labels = np.array([(y * som_x + x) for (x, y) in cluster_assignments])
    unique_clusters = np.unique(cluster_labels)
    centroids = []

    for cluster_id in unique_clusters:
        cluster_indices = np.where(cluster_labels == cluster_id)
        cluster_points = data_points[cluster_indices]
        cluster_energies = free_energies[cluster_indices]
        centroid = np.mean(cluster_points, axis=0)
        mean_energy = np.mean(cluster_energies)

        centroids.append({
            "id": f"Centroid_Cluster_{cluster_id}",
            "coords": centroid,
            "mean_free_energy": mean_energy
        })

    sdf_path = os.path.join(OUTPUT_DIR, f"cluster_centroids{prefix}.sdf")
    with open(sdf_path, 'w') as f:
        for centroid in centroids:
            # Determine atom type based on mean free energy
            atom_type = 'O' if centroid['mean_free_energy'] < 0 else 'C'

            f.write(f"{centroid['id']}\n")
            f.write("\n")
            f.write("\n")
            f.write("  1  0  0  0  0  0  0  0  0  0999 V2000\n")
            f.write(
                f"{centroid['coords'][0]:10.4f}{centroid['coords'][1]:10.4f}{centroid['coords'][2]:10.4f} {atom_type:<3s} 0  0  0  0  0  0  0  0  0  0  0  0\n")
            f.write("M  END\n")
            f.write(f"> <mean_free_energy>\n{centroid['mean_free_energy']:.4f}\n")
            f.write("$$$$\n")

    print(f"Centroids SDF file saved to: {sdf_path}")


def generate_centroids_png(water_molecules, som, cluster_assignments, prefix=""):
    """
    Generates a 3D PNG plot of the cluster centroids.
    """
    print("Generating centroids PNG plot...")
    data_points = np.array([mol['coords'] for mol in water_molecules])
    som_x, som_y, _ = som.get_weights().shape
    cluster_labels = np.array([(y * som_x + x) for (x, y) in cluster_assignments])
    unique_clusters = np.unique(cluster_labels)
    centroids_data = []

    for cluster_id in unique_clusters:
        cluster_points = data_points[np.where(cluster_labels == cluster_id)]
        centroid = np.mean(cluster_points, axis=0)
        centroids_data.append({"coords": centroid, "id": cluster_id})

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    centroids_np = np.array([c['coords'] for c in centroids_data])
    colors = plt.cm.nipy_spectral(np.linspace(0, 1, len(centroids_data)))

    ax.scatter(centroids_np[:, 0], centroids_np[:, 1], centroids_np[:, 2],
               c=colors, marker='X', s=200, label='Centroids')

    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')
    ax.set_title(f"3D Visualization of Cluster Centroids ({'All' if not prefix else 'Filtered'})")
    png_path = os.path.join(OUTPUT_DIR, f"cluster_centroids_plot{prefix}.png")
    plt.savefig(png_path)
    plt.close()
    print(f"Centroids PNG plot saved to: {png_path}")


def generate_key_waters_sdf(key_waters, prefix=""):
    """
    Generates an SDF file for the selected key water molecules.
    """
    print("Generating key waters SDF file...")

    sdf_path = os.path.join(OUTPUT_DIR, f"key_waters{prefix}.sdf")
    with open(sdf_path, 'w') as f:
        for mol in key_waters:
            # Find the line with the O atom and replace it with C if energy is positive
            lines = mol['block'].split('\n')
            new_block_lines = []

            for line in lines:
                parts = line.split()
                if len(parts) >= 4 and parts[3].strip() == 'O':
                    if mol['free_energy'] is not None and mol['free_energy'] >= 0:
                        parts[3] = 'C'
                    new_block_lines.append("     " + "     ".join(parts))
                else:
                    new_block_lines.append(line)

            modified_block = '\n'.join(new_block_lines)
            f.write(modified_block)

    print(f"Key waters SDF file saved to: {sdf_path}")


def generate_key_waters_png(key_waters, prefix=""):
    """
    Generates a 3D PNG plot for the selected key water molecules.
    """
    print("Generating key waters PNG plot...")

    data_points = np.array([mol['coords'] for mol in key_waters])
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(data_points[:, 0], data_points[:, 1], data_points[:, 2],
               marker='o', s=100, color='blue')

    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')
    ax.set_title(f"3D Visualization of Key Water Molecules ({'All' if not prefix else 'Filtered'})")
    png_path = os.path.join(OUTPUT_DIR, f"key_waters_plot{prefix}.png")
    plt.savefig(png_path)
    plt.close()

    print(f"Key waters PNG plot saved to: {png_path}")


def generate_statistics_txt(water_molecules, som, cluster_assignments, prefix=""):
    """
    Generates a text file with statistics for each cluster.
    """
    print("Generating statistics text file...")
    som_x, som_y, _ = som.get_weights().shape
    cluster_labels = np.array([(y * som_x + x) for (x, y) in cluster_assignments])
    unique_clusters = np.unique(cluster_labels)
    stats_path = os.path.join(OUTPUT_DIR, f"statistics{prefix}.txt")

    with open(stats_path, 'w') as f:
        f.write(f"--- Cluster Statistics Report ({'All' if not prefix else 'Filtered'}) ---\n\n")

        for cluster_id in unique_clusters:
            cluster_indices = np.where(cluster_labels == cluster_id)
            cluster_molecules = [water_molecules[i] for i in cluster_indices[0]]
            num_molecules = len(cluster_molecules)

            free_energies = [mol['free_energy'] for mol in cluster_molecules if mol['free_energy'] is not None]
            probabilities = [mol['probability'] for mol in cluster_molecules if mol['probability'] is not None]

            mean_energy = np.mean(free_energies) if free_energies else 0
            num_positive = sum(1 for e in free_energies if e > 0)
            num_negative = sum(1 for e in free_energies if e < 0)

            f.write(f"Cluster ID: {cluster_id}\n")
            f.write(f"  Number of Molecules: {num_molecules}\n")
            f.write(f"  Mean Free Energy: {mean_energy:.4f}\n")
            f.write(f"  Number of Positive Free Energies: {num_positive}\n")
            f.write(f"  Number of Negative Free Energies: {num_negative}\n")
            f.write(f"  All Free Energies: {', '.join([f'{e:.2f}' for e in free_energies])}\n\n")
            f.write(f"  All Probabilities: {', '.join([f'{p:.2f}' for p in probabilities])}\n\n")

    print(f"Statistics text file saved to: {stats_path}")


def generate_key_waters_statistics_txt(key_waters, prefix=""):
    """
    Generates a text file with statistics for the selected key water molecules.
    """
    print("Generating key waters statistics text file...")
    stats_path = os.path.join(OUTPUT_DIR, f"key_waters_statistics{prefix}.txt")

    with open(stats_path, 'w') as f:
        f.write(f"--- Key Water Statistics Report ({'All' if not prefix else 'Filtered'}) ---\n\n")

        for mol in key_waters:
            f.write(f"Molecule ID: {mol['id']}\n")
            f.write(f"  Coordinates: {mol['coords']}\n")
            f.write(f"  szmap_neut_diff_free_energy: {mol['free_energy']:.4f}\n")
            f.write(f"  szmap_probability: {mol['probability']:.4f}\n\n")

    print(f"Key waters statistics text file saved to: {stats_path}")


# --- Script Execution ---
if __name__ == "__main__":
    try:
        water_molecules = parse_sdf_files(INPUT_DIR)

        if water_molecules:
            molecules_with_free_energy = [mol for mol in water_molecules if
                                          mol['free_energy'] is not None and mol['probability'] is not None]

            if not molecules_with_free_energy:
                print(
                    "No molecules with 'szmap_neut_diff_free_energy' and 'szmap_probability' data found. Cannot perform 4D clustering.")
            else:
                # Prepare a 4D array with coordinates and free energy for clustering
                coords_energy_array = np.array(
                    [[mol['coords'][0], mol['coords'][1], mol['coords'][2], mol['free_energy']] for mol in
                     molecules_with_free_energy])

                som, cluster_assignments = run_som_clustering(coords_energy_array, som_x=20, som_y=20)

                print("\n--- Generating files for ALL clusters ---")
                generate_outputs(molecules_with_free_energy, som, cluster_assignments, prefix="_all")
                generate_centroids_sdf(molecules_with_free_energy, som, cluster_assignments, prefix="_all")
                generate_centroids_png(molecules_with_free_energy, som, cluster_assignments, prefix="_all")
                generate_statistics_txt(molecules_with_free_energy, som, cluster_assignments, prefix="_all")
                key_waters_all = select_key_waters(molecules_with_free_energy, cluster_assignments)
                if key_waters_all:
                    generate_key_waters_sdf(key_waters_all, prefix="_key_all")
                    generate_key_waters_png(key_waters_all, prefix="_key_all")
                    generate_key_waters_statistics_txt(key_waters_all, prefix="_key_all")

                filtered_molecules, filtered_assignments, filtered_cluster_ids = filter_clusters(
                    molecules_with_free_energy, som, cluster_assignments, MIN_CLUSTER_SIZE)

                if filtered_molecules:
                    print("\n--- Generating files for FILTERED clusters ---")
                    generate_outputs(filtered_molecules, som, filtered_assignments, prefix="_filtered")
                    generate_centroids_sdf(filtered_molecules, som, filtered_assignments, prefix="_filtered")
                    generate_centroids_png(filtered_molecules, som, filtered_assignments, prefix="_filtered")
                    generate_statistics_txt(filtered_molecules, som, filtered_assignments, prefix="_filtered")
                    key_waters_filtered = select_key_waters(filtered_molecules, filtered_assignments)
                    if key_waters_filtered:
                        generate_key_waters_sdf(key_waters_filtered, prefix="_key_filtered")
                        generate_key_waters_png(key_waters_filtered, prefix="_key_filtered")
                        generate_key_waters_statistics_txt(key_waters_filtered, prefix="_key_filtered")
                else:
                    print("No clusters met the minimum size requirement. No filtered output files were generated.")

    except FileNotFoundError as e:
        print(f"Error: {e}. Please ensure the input directory path is correct.")
    except ImportError as e:
        print(f"Error: {e}. Please install the required library (MiniSom) using 'pip install MiniSom'.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
