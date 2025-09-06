import os

def filter_sdf(input_file, output_file):
    seen_oxygens = set()
    keep_molecules = []

    with open(input_file, "r") as f:
        molecule_block = []
        oxygen_coords = None

        for line in f:
            molecule_block.append(line)

            # Detect end of molecule
            if line.strip() == "$$$$":
                # Extract oxygen coordinates (first O line)
                for atom_line in molecule_block:
                    parts = atom_line.split()
                    if len(parts) >= 4 and parts[3] == "O":
                        oxygen_coords = tuple(parts[0:3])  # (x, y, z) as strings
                        break

                # Keep only first occurrence of each oxygen
                if oxygen_coords not in seen_oxygens:
                    seen_oxygens.add(oxygen_coords)
                    keep_molecules.extend(molecule_block)

                # Reset for next molecule
                molecule_block = []

    # Write filtered molecules to output file
    with open(output_file, "w") as f:
        f.writelines(keep_molecules)


def process_all_files(input_folder, output_folder):
    # Make sure output folder exists
    os.makedirs(output_folder, exist_ok=True)

    for filename in os.listdir(input_folder):
        if filename.lower().endswith(".sdf"):
            input_file = os.path.join(input_folder, filename)

            # Build output filename with "_cleaned"
            base, ext = os.path.splitext(filename)
            output_file = os.path.join(output_folder, f"{base}_cleaned{ext}")

            print(f"Processing {input_file} -> {output_file}")
            filter_sdf(input_file, output_file)


if __name__ == "__main__":
    input_folder = r"C:\Users\ifige\Documents\GitHub\SZMAP\szmap_data"
    output_folder = r"C:\Users\ifige\Documents\GitHub\SZMAP\szmap_unique_waters"

    process_all_files(input_folder, output_folder)
    print("✅ All files processed successfully!")
