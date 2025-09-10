import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.cluster import AgglomerativeClustering
from collections import Counter, defaultdict
from mpl_toolkits.mplot3d import Axes3D
from sklearn.preprocessing import StandardScaler
import csv

# --- Configuration ---
# Set the input and output directories as specified in the prompt.
INPUT_DIR = r"C:\Users\ifige\Documents\GitHub\SZMAP\szmap_unique_waters"
OUTPUT_DIR = r"C:\Users\ifige\Documents\GitHub\SZMAP\Hierarchical_cl_results"
# A large number of initial clusters to capture all small groupings.
N_CLUSTERS = 500
# The minimum number of molecules required to form a valid cluster.
MIN_CLUSTER_SIZE = 200


def parse_sdf_file(filepath):
    """
    Parses a single SDF file to extract coordinates, energy, and probability
    for each water molecule.

    Args:
        filepath (str): The path to the SDF file.

    Returns:
        A list of dictionaries, where each dictionary contains the parsed data
        for a single water molecule.
    """
    water_data = []
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()

        # A simple state machine to parse molecule blocks
        current_data = {}

        for line in lines:
            if 'V2000' in line:
                current_data = {}
            elif '$$$$' in line:
                if 'coords' in current_data and 'energy' in current_data and 'probability' in current_data:
                    water_data.append(current_data)
                current_data = {}

            parts = line.strip().split()
            if len(parts) >= 4 and parts[3] == 'O':
                current_data['coords'] = (float(parts[0]), float(parts[1]), float(parts[2]))
            elif "> <szmap_neut_diff_free_energy>" in line:
                energy_line_index = lines.index(line) + 1
                if energy_line_index < len(lines):
                    current_data['energy'] = float(lines[energy_line_index].strip())
            elif "> <szmap_probability>" in line:
                probability_line_index = lines.index(line) + 1
                if probability_line_index < len(lines):
                    current_data['probability'] = float(lines[probability_line_index].strip())

    except FileNotFoundError:
        print(f"File not found: {filepath}")
    except Exception as e:
        print(f"Error parsing file {filepath}: {e}")

    return water_data


def write_molecule_to_sdf(out_f, coords, energy, atom_type='O'):
    """
    Writes a single molecule block to an SDF file with a specific atom type.
    """
    out_f.write("water\n")
    out_f.write(" -OEChem-09022514583D\n")
    out_f.write("\n")
    out_f.write("  3  2  0  0  0  0  0  0  0  0999 V2000\n")
    out_f.write(f"{coords[0]:10.4f}{coords[1]:10.4f}{coords[2]:10.4f} {atom_type} 0 0 0 0 0 0 0 0 0 0 0 0\n")
    out_f.write(f"{-12.5495:10.4f}{-23.3919:10.4f}{-5.8184:10.4f} H 0 0 0 0 0 0 0 0 0 0 0 0\n")
    out_f.write(f"{-12.7006:10.4f}{-22.5587:10.4f}{-4.5635:10.4f} H 0 0 0 0 0 0 0 0 0 0 0 0\n")
    out_f.write("  1  2  1  0  0  0  0\n")
    out_f.write("  1  3  1  0  0  0  0\n")
    out_f.write("M  END\n")
    out_f.write(f"> <szmap_neut_diff_free_energy>\n")
    out_f.write(f"{energy}\n")
    out_f.write("$$$$\n")


def main():
    """
    Main function to orchestrate the Hierarchical clustering process.
    """
    print(f"Reading SDF files from: {INPUT_DIR}")

    if not os.path.exists(INPUT_DIR):
        print(f"Error: Input directory does not exist at {INPUT_DIR}")
        return

    # Create the output directory if it doesn't exist
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        print(f"Created output directory: {OUTPUT_DIR}")

    all_sdf_files = glob.glob(os.path.join(INPUT_DIR, '*.sdf'))
    if not all_sdf_files:
        print("No SDF files found. Please check your input directory and file extension.")
        return

    all_water_data = []

    # Process all SDF files and collect data
    for file in all_sdf_files:
        data = parse_sdf_file(file)
        all_water_data.extend(data)

    if not all_water_data:
        print("No water data was found in the SDF files.")
        return

    print(f"Found {len(all_water_data)} water molecules for clustering.")

    # Convert all data to a NumPy array for later use (plotting, CSV, etc.)
    data_array = np.array([list(d['coords']) + [d['energy'], d['probability']] for d in all_water_data])

    # Create a separate array with only the coordinates for clustering
    clustering_data = data_array[:, :3]

    # Scale the coordinate data to ensure all dimensions contribute equally to the clustering
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(clustering_data)

    # Perform Hierarchical Clustering with a large number of clusters to get granular results
    print(f"Performing Hierarchical Clustering with {N_CLUSTERS} initial clusters on coordinates only...")

    clustering = AgglomerativeClustering(n_clusters=N_CLUSTERS, linkage='ward')
    raw_cluster_labels = clustering.fit_predict(scaled_data)

    # Filter clusters based on the minimum size
    cluster_counts = Counter(raw_cluster_labels)
    final_cluster_labels = np.array(
        [label if cluster_counts[label] >= MIN_CLUSTER_SIZE else -1 for label in raw_cluster_labels])

    print(
        f"Clustering complete. Identified {len(set(final_cluster_labels)) - 1} valid clusters and labeled the rest as noise.")
    print("Generating output files...")

    # --- Generate Outputs ---

    # 1. Dendrogram plot
    plt.figure(figsize=(15, 7))
    linkage_matrix = linkage(scaled_data, method='ward')
    dendrogram(linkage_matrix)
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('Sample Index')
    plt.ylabel('Distance')
    dendrogram_path = os.path.join(OUTPUT_DIR, 'dendrogram.png')
    plt.savefig(dendrogram_path)
    plt.close()
    print(f"Dendrogram saved to: {dendrogram_path}")

    # 2. 3D scatter plot of the clusters and noise, colored by energy
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot noise points in a distinct color
    noise_mask = (final_cluster_labels == -1)
    ax.scatter(data_array[noise_mask, 0], data_array[noise_mask, 1], data_array[noise_mask, 2], c='k', marker='x', s=50,
               label='Noise')

    # Plot clustered points and color by energy
    scatter = ax.scatter(data_array[~noise_mask, 0], data_array[~noise_mask, 1], data_array[~noise_mask, 2],
                         c=data_array[~noise_mask, 3], cmap='plasma', s=50)

    ax.set_title(f'Hierarchical Clustering Results (3D Plot, {len(set(final_cluster_labels)) - 1} Clusters)')
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')

    # Add a color bar for the energy values
    cbar = fig.colorbar(scatter, ax=ax, pad=0.1)
    cbar.set_label('SZMAP Neut Diff Free Energy')

    clusters_3d_path = os.path.join(OUTPUT_DIR, 'clusters_3d_plot_all.png')
    plt.legend()
    plt.savefig(clusters_3d_path)
    plt.close()
    print(f"3D cluster plot (all points) saved to: {clusters_3d_path}")

    # 3. New 3D scatter plot showing ONLY the clustered waters
    fig_clusters = plt.figure(figsize=(10, 8))
    ax_clusters = fig_clusters.add_subplot(111, projection='3d')

    valid_mask = (final_cluster_labels != -1)
    # Get a unique list of valid cluster IDs to use for coloring
    unique_clusters = np.unique(final_cluster_labels[valid_mask])

    # Get a colormap with enough colors
    cmap = plt.cm.get_cmap('gist_rainbow', len(unique_clusters))

    for i, cluster_id in enumerate(unique_clusters):
        cluster_mask = (final_cluster_labels == cluster_id)
        ax_clusters.scatter(data_array[cluster_mask, 0], data_array[cluster_mask, 1], data_array[cluster_mask, 2],
                            color=cmap(i), s=50, label=f'Cluster {cluster_id}')

    ax_clusters.set_title(f'Valid Clusters Only ({len(unique_clusters)} Clusters)')
    ax_clusters.set_xlabel('X Coordinate')
    ax_clusters.set_ylabel('Y Coordinate')
    ax_clusters.set_zlabel('Z Coordinate')

    clusters_only_path = os.path.join(OUTPUT_DIR, 'clusters_3d_plot_valid_clusters_only.png')
    plt.legend()
    plt.savefig(clusters_only_path)
    plt.close()
    print(f"3D plot for valid clusters saved to: {clusters_only_path}")

    # 4. New 3D scatter plot showing ONLY the noise waters
    fig_noise = plt.figure(figsize=(10, 8))
    ax_noise = fig_noise.add_subplot(111, projection='3d')

    ax_noise.scatter(data_array[noise_mask, 0], data_array[noise_mask, 1], data_array[noise_mask, 2],
                     c='k', marker='x', s=50, label='Noise')

    ax_noise.set_title('Noise Waters Only')
    ax_noise.set_xlabel('X Coordinate')
    ax_noise.set_ylabel('Y Coordinate')
    ax_noise.set_zlabel('Z Coordinate')

    noise_only_path = os.path.join(OUTPUT_DIR, 'clusters_3d_plot_noise_only.png')
    plt.legend()
    plt.savefig(noise_only_path)
    plt.close()
    print(f"3D plot for noise waters saved to: {noise_only_path}")

    # 5. CSV file with coordinates and cluster labels
    csv_path = os.path.join(OUTPUT_DIR, 'cluster_data.csv')
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['x', 'y', 'z', 'szmap_neut_diff_free_energy', 'szmap_probability', 'cluster_id'])
        for i, data in enumerate(all_water_data):
            writer.writerow(
                [data['coords'][0], data['coords'][1], data['coords'][2], data['energy'], data['probability'],
                 final_cluster_labels[i]])
    print(f"CSV data saved to: {csv_path}")

    # Group data by final cluster labels for the new outputs
    clustered_data = defaultdict(list)
    for i, label in enumerate(final_cluster_labels):
        if label != -1:
            clustered_data[label].append(all_water_data[i])

    # 6. New SDF file with cluster means
    mean_sdf_path = os.path.join(OUTPUT_DIR, 'mean_cluster_waters.sdf')
    with open(mean_sdf_path, 'w') as out_f:
        print("Writing mean-value SDF file...")
        for cluster_id, data_list in clustered_data.items():
            coords = np.mean([d['coords'] for d in data_list], axis=0)
            energy = np.mean([d['energy'] for d in data_list])
            atom_type = 'O' if energy < 0 else 'C'
            write_molecule_to_sdf(out_f, coords, energy, atom_type)
        print(f"Mean-value SDF file saved to: {mean_sdf_path}")

    # 7. New SDF file with representative molecules
    representative_sdf_path = os.path.join(OUTPUT_DIR, 'representative_cluster_waters.sdf')
    with open(representative_sdf_path, 'w') as out_f:
        print("Writing representative-molecule SDF file...")
        for cluster_id, data_list in clustered_data.items():
            # Find the molecule(s) with the highest probability
            max_prob = max([d['probability'] for d in data_list])
            highest_prob_mols = [d for d in data_list if d['probability'] == max_prob]

            # Determine if we should pick lowest or highest energy
            negative_energy_count = sum(1 for d in data_list if d['energy'] < 0)
            positive_energy_count = len(data_list) - negative_energy_count

            if negative_energy_count >= positive_energy_count:
                # Pick the one with the lowest energy
                representative_mol = min(highest_prob_mols, key=lambda x: x['energy'])
            else:
                # Pick the one with the highest energy
                representative_mol = max(highest_prob_mols, key=lambda x: x['energy'])

            coords = representative_mol['coords']
            energy = representative_mol['energy']
            atom_type = 'O' if energy < 0 else 'C'
            write_molecule_to_sdf(out_f, coords, energy, atom_type)
        print(f"Representative-molecule SDF file saved to: {representative_sdf_path}")

    # 8. Text file with cluster statistics
    stats_path = os.path.join(OUTPUT_DIR, 'cluster_statistics.txt')
    with open(stats_path, 'w') as out_f:
        print("Writing cluster statistics file...")
        out_f.write(f"Cluster Statistics (Valid Clusters Only)\n\n")

        for cluster_id, data_list in clustered_data.items():
            energies = [d['energy'] for d in data_list]
            avg_energy = np.mean(energies)

            negative_count = sum(1 for e in energies if e < 0)
            positive_count = len(energies) - negative_count

            out_f.write(f"--- Cluster {cluster_id} ---\n")
            out_f.write(f"Number of molecules: {len(data_list)}\n")
            out_f.write(f"Average Energy: {avg_energy:.4f}\n")
            out_f.write(f"Number of negative energies: {negative_count}\n")
            out_f.write(f"Number of positive energies: {positive_count}\n")
            out_f.write(f"All Energies: {', '.join([f'{e:.4f}' for e in energies])}\n\n")
        print(f"Cluster statistics saved to: {stats_path}")

    # 9. New SDF file with cluster labels
    print("Writing new SDF file with cluster IDs...")

    # Read all input files to get all molecule data in order
    all_molecule_blocks = []
    for file in all_sdf_files:
        with open(file, 'r') as in_f:
            lines = in_f.read().strip().split('$$$$\n')
            for block in lines:
                if block.strip():
                    all_molecule_blocks.append(block)

    # Write out the new SDF file with cluster IDs
    output_sdf_path = os.path.join(OUTPUT_DIR, 'hierarchical_clustered_waters.sdf')
    with open(output_sdf_path, 'w') as out_f:
        for i, block in enumerate(all_molecule_blocks):
            if i < len(final_cluster_labels):
                out_f.write(block)
                out_f.write(f"\n> <szmap_cluster_id>\n")
                out_f.write(f"{final_cluster_labels[i]}\n")
                out_f.write("$$$$\n")
    print(f"New SDF file saved to: {output_sdf_path}")
    print("Hierarchical clustering process finished successfully.")


if __name__ == "__main__":
    main()