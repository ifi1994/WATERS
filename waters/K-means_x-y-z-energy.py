import os
import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from mpl_toolkits.mplot3d import Axes3D

# --- User-defined Parameters ---
# The minimum number of clusters to check for the Elbow Method.
MIN_CLUSTERS = 100
# The maximum number of clusters to check for the Elbow Method.
MAX_CLUSTERS = 150
# The minimum number of atoms a cluster must contain to be included in the final output.
MIN_CLUSTER_SIZE = 100

# --- File Paths ---
# The directory containing your SDF files.
INPUT_DIR = r"C:\Users\ifige\Documents\GitHub\SZMAP\szmap_unique_waters"
# The directory where the results will be saved.
OUTPUT_DIR = r"C:\Users\ifige\Documents\GitHub\SZMAP\K-means_results"


def find_optimal_k(sse_list):
    """
    Finds the optimal number of clusters (K) using the Elbow Method.
    This function finds the point on the curve with the maximum distance
    from a straight line connecting the first and last points.

    Args:
        sse_list (list): A list of SSE values for each K.

    Returns:
        int: The optimal K value.
    """
    if len(sse_list) < 3:
        # Not enough points to find an elbow.
        return len(sse_list)

    # Coordinates of the first and last points of the curve
    start_point = np.array([MIN_CLUSTERS, sse_list[0]])
    end_point = np.array([MAX_CLUSTERS, sse_list[-1]])
    line_vector = end_point - start_point

    distances = []
    for i, sse in enumerate(sse_list):
        point = np.array([MIN_CLUSTERS + i, sse])
        # Vector from the start point to the current point
        point_vector = point - start_point
        # Calculate the cross product of the vectors
        cross_product = np.abs(line_vector[0] * point_vector[1] - line_vector[1] * point_vector[0])
        # Calculate the length of the line vector
        line_length = np.linalg.norm(line_vector)
        # Calculate the distance
        distance = cross_product / line_length
        distances.append(distance)

    optimal_index = np.argmax(distances)
    return MIN_CLUSTERS + optimal_index


def parse_sdf_file(filepath):
    """
    Parses a single SDF file to extract oxygen coordinates, molecule block,
    and szmap_neut_diff_free_energy for each molecule.

    Args:
        filepath (str): The path to the SDF file.

    Returns:
        list: A list of dictionaries, where each dictionary contains
              the molecule's data, including its block, oxygen coordinates,
              free energy value, and probability.
    """
    molecules_data = []
    current_mol_data = {'block': '', 'oxygen_coord': None, 'free_energy': None, 'probability': None}
    is_mol_block = False
    reading_free_energy = False
    reading_probability = False

    with open(filepath, 'r') as f:
        for line in f:
            current_mol_data['block'] += line

            if '$$$$' in line:
                if current_mol_data and current_mol_data['oxygen_coord'] is not None:
                    molecules_data.append(current_mol_data)
                current_mol_data = {'block': '', 'oxygen_coord': None, 'free_energy': None, 'probability': None}
                reading_free_energy = False
                reading_probability = False
                continue

            if 'V2000' in line:
                is_mol_block = True
            elif 'M END' in line:
                is_mol_block = False

            if is_mol_block and ' O ' in line:
                parts = line.split()
                try:
                    x = float(parts[0])
                    y = float(parts[1])
                    z = float(parts[2])
                    current_mol_data['oxygen_coord'] = [x, y, z]
                except (ValueError, IndexError):
                    print(f"Warning: Could not parse coordinates from line: {line.strip()}")

            # Use flags to read the value from the next line
            if reading_free_energy:
                try:
                    current_mol_data['free_energy'] = float(line.strip())
                except (ValueError, IndexError):
                    current_mol_data['free_energy'] = None
                reading_free_energy = False

            if reading_probability:
                try:
                    current_mol_data['probability'] = float(line.strip())
                except (ValueError, IndexError):
                    current_mol_data['probability'] = None
                reading_probability = False

            if '> <szmap_neut_diff_free_energy>' in line:
                reading_free_energy = True

            if '> <szmap_probability>' in line:
                reading_probability = True

    if current_mol_data['block'].strip() and current_mol_data['oxygen_coord'] is not None:
        molecules_data.append(current_mol_data)

    return molecules_data


def main():
    """
    Main function to orchestrate the clustering process.
    """
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        print(f"Created output directory: {OUTPUT_DIR}")

    print("Starting data parsing from SDF files...")

    all_molecules = []

    sdf_files = glob.glob(os.path.join(INPUT_DIR, "*.sdf"))
    if not sdf_files:
        print(f"No .sdf files found in {INPUT_DIR}. Exiting.")
        return

    for sdf_file in sdf_files:
        mols = parse_sdf_file(sdf_file)
        all_molecules.extend(mols)

    all_features = []
    clustered_molecules = []

    for mol in all_molecules:
        if mol['oxygen_coord'] is not None and mol['free_energy'] is not None:
            # Create a 4-dimensional feature vector: x, y, z, and free energy
            all_features.append(mol['oxygen_coord'] + [mol['free_energy']])
            clustered_molecules.append(mol)

    if not all_features:
        print("No oxygen atoms with free energy values found in the SDF files. Exiting.")
        return

    print(f"Found {len(all_features)} oxygen atoms with energy data for clustering.")

    X = np.array(all_features)

    print("Generating Elbow Plot to find optimal K...")
    sse = []
    k_range = range(MIN_CLUSTERS, MAX_CLUSTERS + 1)
    for k in k_range:
        km = KMeans(n_clusters=k, random_state=42, n_init=10)
        km.fit(X)
        sse.append(km.inertia_)

    plt.figure(figsize=(10, 6))
    plt.plot(k_range, sse, marker='o')
    plt.title('Elbow Method for Optimal K (4D Clustering)')
    plt.xlabel('Number of Clusters (K)')
    plt.ylabel('Sum of Squared Errors (SSE)')
    plt.xticks(k_range)
    elbow_plot_path = os.path.join(OUTPUT_DIR, 'elbow_method_plot.png')
    plt.savefig(elbow_plot_path)
    plt.close()
    print(f"Elbow plot saved to {elbow_plot_path}")

    optimal_k = find_optimal_k(sse)

    if optimal_k == 1 and MAX_CLUSTERS > 1:
        print("Warning: The optimal K was found to be 1, which may not be a useful clustering.")
        print("Falling back to the second best option if available.")
        distances = []
        start_point = np.array([MIN_CLUSTERS, sse[0]])
        end_point = np.array([MAX_CLUSTERS, sse[-1]])
        line_vector = end_point - start_point
        for i, sse_val in enumerate(sse):
            point = np.array([MIN_CLUSTERS + i, sse_val])
            point_vector = point - start_point
            distance = np.abs(line_vector[0] * point_vector[1] - line_vector[1] * point_vector[0]) / np.linalg.norm(
                line_vector)
            distances.append(distance)

        sorted_distances = sorted(distances, reverse=True)
        if len(sorted_distances) > 1:
            second_best_distance = sorted_distances[1]
            optimal_k = MIN_CLUSTERS + distances.index(second_best_distance)
        else:
            print("No second best option available. Using K=1.")

    print(f"\nOptimal number of clusters (K) calculated as: {optimal_k}")

    print(f"Performing final K-means clustering with {optimal_k} clusters...")
    kmeans = KMeans(n_clusters=optimal_k, random_state=42, n_init=10)
    kmeans.fit(X)
    labels = kmeans.labels_
    centroids = kmeans.cluster_centers_

    print("\nProcessing and filtering clusters...")

    clusters = {i: [] for i in range(optimal_k)}
    for i, mol in enumerate(clustered_molecules):
        mol['cluster_id'] = labels[i]
        clusters[labels[i]].append(mol)

    filtered_clusters_data = {k: v for k, v in clusters.items() if len(v) >= MIN_CLUSTER_SIZE}
    print(f"Filtered down to {len(filtered_clusters_data)} clusters with at least {MIN_CLUSTER_SIZE} members.")

    if not filtered_clusters_data:
        print("No clusters met the minimum size requirement. Exiting.")
        return

    filtered_coords = []
    filtered_labels = []
    filtered_free_energies = []
    filtered_molecules = []

    for cluster_id, cluster_mols in filtered_clusters_data.items():
        for mol in cluster_mols:
            filtered_coords.append(mol['oxygen_coord'])
            filtered_labels.append(cluster_id)
            filtered_free_energies.append(mol['free_energy'])
            filtered_molecules.append(mol)

    filtered_X = np.array(filtered_coords)
    filtered_centroids = np.array([centroids[c] for c in filtered_clusters_data.keys()])

    print("\nGenerating output files for filtered clusters...")

    print("Generating 3D cluster plot...")
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(filtered_X[:, 0], filtered_X[:, 1], filtered_X[:, 2], c=filtered_labels, cmap='viridis', s=50,
                         alpha=0.8)
    ax.scatter(filtered_centroids[:, 0], filtered_centroids[:, 1], filtered_centroids[:, 2], marker='X', s=200, c='red',
               label='Centroids')
    ax.set_title(f'4D Clustering (K={optimal_k}) with Size Filter')
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')
    ax.legend()
    plt.colorbar(scatter, label='Cluster ID')
    cluster_plot_path = os.path.join(OUTPUT_DIR, 'filtered_cluster_plot.png')
    plt.savefig(cluster_plot_path)
    plt.close()
    print(f"3D cluster plot saved to {cluster_plot_path}")

    print("Generating CSV file...")
    data = {
        'x': filtered_X[:, 0],
        'y': filtered_X[:, 1],
        'z': filtered_X[:, 2],
        'cluster_id': filtered_labels,
        'szmap_neut_diff_free_energy': filtered_free_energies
    }
    df = pd.DataFrame(data)
    csv_path = os.path.join(OUTPUT_DIR, 'filtered_clustering_results.csv')
    df.to_csv(csv_path, index=False)
    print(f"CSV file saved to {csv_path}")

    print("Generating new SDF file with cluster data...")
    new_sdf_path = os.path.join(OUTPUT_DIR, 'filtered_clustered_waters.sdf')
    with open(new_sdf_path, 'w') as f_out:
        for mol_data in filtered_molecules:
            block = mol_data['block']
            cluster_id_tag = f"> <cluster_id>\n{mol_data['cluster_id']}\n\n"
            block_lines = block.split('\n')
            updated_block = ""
            for line in block_lines:
                updated_block += line + '\n'
                if line.strip() == 'M END':
                    updated_block += cluster_id_tag
            f_out.write(updated_block)
            f_out.write('$$$$\n')

    print(f"New SDF file saved to {new_sdf_path}")

    print("Generating cluster statistics file...")
    stats_path = os.path.join(OUTPUT_DIR, 'cluster_stats.txt')
    with open(stats_path, 'w') as f_out:
        f_out.write("Cluster Statistics (for clusters with at least {} atoms)\n".format(MIN_CLUSTER_SIZE))
        f_out.write("=" * 60 + "\n\n")

        for cluster_id, cluster_mols in filtered_clusters_data.items():
            energies = [mol['free_energy'] for mol in cluster_mols if mol['free_energy'] is not None]
            probabilities = [mol['probability'] for mol in cluster_mols if mol['probability'] is not None]

            f_out.write(f"Cluster ID: {cluster_id}\n")
            f_out.write(f"Number of oxygen atoms: {len(cluster_mols)}\n")

            energies_str = ', '.join([f"{e:.4f}" for e in energies]) if energies else "No energy values found"
            f_out.write(f"List of szmap_neut_diff_free_energy: [{energies_str}]\n")

            probabilities_str = ', '.join(
                [f"{p:.4f}" for p in probabilities]) if probabilities else "No probability values found"
            f_out.write(f"List of szmap_probability: [{probabilities_str}]\n")

            pos_energies = [e for e in energies if e > 0]
            neg_energies = [e for e in energies if e < 0]
            f_out.write(f"Number of positive energies: {len(pos_energies)}\n")
            f_out.write(f"Number of negative energies: {len(neg_energies)}\n")

            avg_energy = np.mean(energies) if energies else np.nan
            avg_energy_str = f"{avg_energy:.4f}" if not np.isnan(avg_energy) else "N/A"
            f_out.write(f"Average energy of the cluster: {avg_energy_str}\n\n")

            # Find and write the Best Spot information
            best_spot = None
            if not np.isnan(avg_energy) and avg_energy >= 0:
                best_spot = max(cluster_mols, key=lambda mol: (
                mol['free_energy'] if mol['free_energy'] is not None else -np.inf,
                mol['probability'] if mol['probability'] is not None else -np.inf))
            elif not np.isnan(avg_energy):
                best_spot = min(cluster_mols, key=lambda mol: (
                mol['free_energy'] if mol['free_energy'] is not None else np.inf,
                mol['probability'] if mol['probability'] is not None else -np.inf))

            f_out.write("Best Spot Information:\n")
            if best_spot:
                energy_str = f"{best_spot['free_energy']:.4f}" if best_spot['free_energy'] is not None else "N/A"
                probability_str = f"{best_spot['probability']:.4f}" if best_spot['probability'] is not None else "N/A"
                coords_str = f"({best_spot['oxygen_coord'][0]:.4f}, {best_spot['oxygen_coord'][1]:.4f}, {best_spot['oxygen_coord'][2]:.4f})"
                f_out.write(f"  Energy: {energy_str}\n")
                f_out.write(f"  Probability: {probability_str}\n")
                f_out.write(f"  Coordinates: {coords_str}\n")
            else:
                f_out.write("  No best spot found for this cluster.\n")
            f_out.write("\n" + "-" * 60 + "\n\n")

    print(f"Cluster statistics file saved to {stats_path}")

    print("Generating SDF file with cluster centroids...")
    centroids_sdf_path = os.path.join(OUTPUT_DIR, 'cluster_centroids.sdf')
    with open(centroids_sdf_path, 'w') as f_out:
        for cluster_id, cluster_mols in filtered_clusters_data.items():
            energies = [mol['free_energy'] for mol in cluster_mols if mol['free_energy'] is not None]
            avg_energy = np.mean(energies) if energies else np.nan
            centroid_coords = centroids[cluster_id]

            avg_energy_str = f"{avg_energy:.4f}" if not np.isnan(avg_energy) else "N/A"

            f_out.write(f"Centroid for Cluster {cluster_id}\n\n\n")
            f_out.write(f" 1 0 0 0 0 0 0 0 0 0999 V2000\n")
            atom_symbol = 'C' if avg_energy >= 0 else 'O'
            f_out.write(
                f"{centroid_coords[0]:10.4f}{centroid_coords[1]:10.4f}{centroid_coords[2]:10.4f} {atom_symbol}   0 0 0 0 0 0 0 0 0 0 0 0\n")
            f_out.write("M  END\n")
            f_out.write(f"> <cluster_id>\n{cluster_id}\n\n")
            f_out.write(f"> <average_szmap_neut_diff_free_energy>\n{avg_energy_str}\n\n")
            f_out.write("$$$$\n")

    print(f"Centroids SDF file saved to {centroids_sdf_path}")

    print("Generating 3D plot of centroids...")
    fig_centroids = plt.figure(figsize=(12, 10))
    ax_centroids = fig_centroids.add_subplot(111, projection='3d')
    centroid_coords = np.array([centroids[c] for c in filtered_clusters_data.keys()])
    centroid_labels = np.array(list(filtered_clusters_data.keys()))
    scatter_centroids = ax_centroids.scatter(centroid_coords[:, 0], centroid_coords[:, 1], centroid_coords[:, 2],
                                             c=centroid_labels, cmap='viridis', s=100, alpha=1.0)
    ax_centroids.set_title('3D Plot of Cluster Centroids (4D Clustering)')
    ax_centroids.set_xlabel('X Coordinate')
    ax_centroids.set_ylabel('Y Coordinate')
    ax_centroids.set_zlabel('Z Coordinate')
    plt.colorbar(scatter_centroids, label='Cluster ID')
    centroids_plot_path = os.path.join(OUTPUT_DIR, 'centroids_plot.png')
    plt.savefig(centroids_plot_path)
    plt.close()
    print(f"3D centroids plot saved to {centroids_plot_path}")

    print("Generating SDF file with best spots...")
    best_spots_sdf_path = os.path.join(OUTPUT_DIR, 'best_spots.sdf')
    best_spots = []
    with open(best_spots_sdf_path, 'w') as f_out:
        for cluster_id, cluster_mols in filtered_clusters_data.items():
            energies = [mol['free_energy'] for mol in cluster_mols if mol['free_energy'] is not None]
            avg_energy = np.mean(energies) if energies else np.nan

            best_spot = None
            if not np.isnan(avg_energy) and avg_energy >= 0:
                best_spot = max(cluster_mols, key=lambda mol: (
                mol['free_energy'] if mol['free_energy'] is not None else -np.inf,
                mol['probability'] if mol['probability'] is not None else -np.inf))
            elif not np.isnan(avg_energy):
                best_spot = min(cluster_mols, key=lambda mol: (
                mol['free_energy'] if mol['free_energy'] is not None else np.inf,
                mol['probability'] if mol['probability'] is not None else -np.inf))

            if best_spot:
                best_spots.append(best_spot)
                best_spot_block = best_spot['block']

                # Replace ' O ' with ' C ' if average energy is positive
                if avg_energy >= 0:
                    best_spot_block = best_spot_block.replace(' O   ', ' C   ', 1)

                energy_str = f"{best_spot['free_energy']:.4f}" if best_spot['free_energy'] is not None else "N/A"
                probability_str = f"{best_spot['probability']:.4f}" if best_spot['probability'] is not None else "N/A"
                avg_energy_str = f"{avg_energy:.4f}" if not np.isnan(avg_energy) else "N/A"

                best_spot_block += f"> <cluster_id>\n{best_spot['cluster_id']}\n\n"
                best_spot_block += f"> <szmap_neut_diff_free_energy>\n{energy_str}\n\n"
                best_spot_block += f"> <szmap_probability>\n{probability_str}\n\n"
                best_spot_block += f"> <average_cluster_energy>\n{avg_energy_str}\n\n"
                f_out.write(best_spot_block)
                f_out.write("$$$$\n")

    print(f"Best spots SDF file saved to {best_spots_sdf_path}")

    print("Generating 3D plot of best spots...")
    fig_best_spots = plt.figure(figsize=(12, 10))
    ax_best_spots = fig_best_spots.add_subplot(111, projection='3d')
    best_spots_coords = np.array([spot['oxygen_coord'] for spot in best_spots])
    best_spots_labels = np.array([spot['cluster_id'] for spot in best_spots])
    scatter_best_spots = ax_best_spots.scatter(best_spots_coords[:, 0], best_spots_coords[:, 1],
                                               best_spots_coords[:, 2], c=best_spots_labels, cmap='viridis', s=100,
                                               alpha=1.0)
    ax_best_spots.set_title('3D Plot of Best Spots (4D Clustering)')
    ax_best_spots.set_xlabel('X Coordinate')
    ax_best_spots.set_ylabel('Y Coordinate')
    ax_best_spots.set_zlabel('Z Coordinate')
    plt.colorbar(scatter_best_spots, label='Cluster ID')
    best_spots_plot_path = os.path.join(OUTPUT_DIR, 'best_spots_plot.png')
    plt.savefig(best_spots_plot_path)
    plt.close()
    print(f"3D best spots plot saved to {best_spots_plot_path}")

    print("\nClustering process completed successfully!")


if __name__ == "__main__":
    main()
