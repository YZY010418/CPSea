import os
import numpy as np
from Bio.PDB import PDBParser
from tqdm import tqdm
import argparse

def add_cb(input_array):
    # from protein mpnn
    # The virtual Cβ coordinates were calculated using ideal angle and bond length definitions: b = Cα - N, c = C - Cα, a = cross(b, c), Cβ = -0.58273431*a + 0.56802827*b - 0.54067466*c + Cα.
    N, CA, C, O = input_array
    b = CA - N
    c = C - CA
    a = np.cross(b, c)
    CB = np.around(-0.58273431 * a + 0.56802827 * b - 0.54067466 * c + CA, 3)
    return CB

def calculate_cb_distance(structure):
    target_chain = None
    for chain in structure.get_chains():
        if chain.get_id() == 'L':
            target_chain = chain
            break

    if target_chain is None:
        print(f"Warning: Chain 'L' not found in the structure. Skipping this calculation.")
        return None

    residues = list(target_chain.get_residues())
    if not residues:
        print(f"Warning: Chain 'L' has no residues. Skipping this calculation.")
        return None

    first_residue = residues[0]
    last_residue = residues[-1]

    if first_residue.get_resname() == 'GLY':
        tmp_coord = np.array([
            first_residue['N'].get_coord(),
            first_residue['CA'].get_coord(),
            first_residue['C'].get_coord(),
            first_residue['O'].get_coord()
        ])
        first_cb = add_cb(tmp_coord)
    else:
        first_cb = first_residue['CB'].get_coord()

    if last_residue.get_resname() == 'GLY':
        tmp_coord = np.array([
            last_residue['N'].get_coord(),
            last_residue['CA'].get_coord(),
            last_residue['C'].get_coord(),
            last_residue['O'].get_coord()
        ])
        last_cb = add_cb(tmp_coord)
    else:
        last_cb = last_residue['CB'].get_coord()

    distance = np.linalg.norm(first_cb - last_cb)
    return distance


def main():
    parser = argparse.ArgumentParser(description='Calculate CB distances for PDB files in a directory.')
    parser.add_argument('input_dir', type=str, help='Directory containing PDB files.')
    parser.add_argument('output_csv', type=str, help='Output CSV file path.')
    args = parser.parse_args()

    pdb_parser = PDBParser()
    results = []

    pdb_files = []
    for root, dirs, files in os.walk(args.input_dir):
        for file in files:
            if file.endswith('.pdb'):
                pdb_files.append(os.path.join(root, file))

    for file_path in tqdm(pdb_files, desc="Processing PDB files"):
        try:
            structure = pdb_parser.get_structure('pdb', file_path)
            distance = calculate_cb_distance(structure)
            if distance is not None:
                results.append((os.path.basename(file_path), distance))
        except Exception as e:
            print(f"Error processing {os.path.basename(file_path)}: {e}")

    with open(args.output_csv, 'w') as f:
        f.write('Filename,CB_Distance\n')
        for filename, distance in results:
            f.write(f'{filename},{distance}\n')


if __name__ == "__main__":
    main()
