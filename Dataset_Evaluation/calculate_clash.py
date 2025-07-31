import os
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from tqdm import tqdm
import argparse
import gzip
import shutil
import multiprocessing
import csv
from functools import partial

vdw_radii = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80,
    'F': 1.47, 'CL': 1.75, 'BR': 1.85, 'I': 1.98
}

buffer = 0.4 

from scipy.spatial.distance import cdist

def compute_clash(structure, receptor_chain_id='R', ligand_chain_id='L', tolerance=0.4):
    receptor_atoms = []
    ligand_atoms = []
    receptor_elements = []
    ligand_elements = []

    for model in structure:
        for chain in model:
            for atom in chain.get_atoms():
                if atom.element == "H":
                    continue
                parent_res = atom.get_parent()
                if atom.get_id() == "CB" and parent_res.resname == "GLY":
                    continue  # skip virtual CBs on glycine
                if chain.id == receptor_chain_id:
                    receptor_atoms.append(atom)
                    receptor_elements.append(atom.element)
                elif chain.id == ligand_chain_id:
                    ligand_atoms.append(atom)
                    ligand_elements.append(atom.element)

    if not receptor_atoms or not ligand_atoms:
        return 0, 0, 0  # nothing to compare

    ligand_coords = np.array([atom.coord for atom in ligand_atoms])
    receptor_coords = np.array([atom.coord for atom in receptor_atoms])
    distances = cdist(ligand_coords, receptor_coords)

    try:
        src_vdw = np.array([vdw_radii[el] for el in receptor_elements])
        dst_vdw = np.array([vdw_radii[el] for el in ligand_elements])
    except KeyError as e:
        print(f"Missing VDW radius for element: {e}")
        return 0, 0, 0

    vdw_sums = dst_vdw[:, np.newaxis] + src_vdw - tolerance
    clash_mask = distances < vdw_sums
    clash_count = np.sum(clash_mask)
    total_pairs = len(ligand_atoms) * len(receptor_atoms)

    # Compute clash_ratio as number of unique ligand atoms involved in clashes divided by total ligand atoms
    clash_ratio = len(np.unique(np.where(clash_mask)[0])) / len(ligand_atoms) if ligand_atoms else 0.0

    return clash_count, total_pairs, clash_ratio

def process_pdb(pdb_path, pdb_dir):
    if pdb_dir:
        pdb_path = os.path.join(pdb_dir, pdb_path)
    try:
        parser = PDBParser(QUIET=True)
        temp_dir = 'tmp_clash'
        os.makedirs(temp_dir, exist_ok=True)
        if pdb_path.endswith('.pdb.gz'):
            pdb_name = os.path.basename(pdb_path).replace('.gz', '')
            temp_pdb_path = os.path.join(temp_dir, pdb_name)
            with gzip.open(pdb_path, 'rb') as f_in:
                with open(temp_pdb_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            structure = parser.get_structure("complex", temp_pdb_path)
        else:
            structure = parser.get_structure("complex", pdb_path)

        clash_count, total_pairs, clash_ratio = compute_clash(structure)

        result = {
            "pdb_path": pdb_path,
            "clash_count": clash_count,
            "total_pairs": total_pairs,
            "clash_ratio": round(clash_ratio, 4)
        }

        if pdb_path.endswith('.pdb.gz'):
            os.remove(temp_pdb_path)

        return result

    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")
        return None

def write_header(output_csv):
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['pdb_path', 'clash_count', 'total_pairs', 'clash_ratio']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

def append_result(output_csv, result):
    if result:
        with open(output_csv, 'a', newline='') as csvfile:
            fieldnames = ['pdb_path', 'clash_count', 'total_pairs', 'clash_ratio']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writerow(result)

def callback(result, tbar, output_csv):
    append_result(output_csv, result)
    tbar.update(1)

def main():
    parser = argparse.ArgumentParser(description='Calculate clash information from PDB files.')
    parser.add_argument('-i', '--input_csv', required=True, help='Input CSV file containing PDB file names.')
    parser.add_argument('-o', '--output_csv', required=True, help='Output CSV file to save clash results.')
    parser.add_argument('-d', '--pdb_dir', default='', help='Directory of PDB files. Leave blank if PDB paths are complete in the input CSV.')
    parser.add_argument('-c', '--cores', type=int, default=multiprocessing.cpu_count(), help='Number of processes to use.')

    args = parser.parse_args()

    input_csv = args.input_csv
    output_csv = args.output_csv
    pdb_dir = args.pdb_dir
    num_processes = args.cores

    write_header(output_csv)
    df = pd.read_csv(input_csv, header=None)
    pdb_paths = df.iloc[:, 0].tolist()

    tbar = tqdm(total=len(pdb_paths))
    partial_callback = partial(callback, tbar=tbar, output_csv=output_csv)

    with multiprocessing.Pool(processes=num_processes) as pool:
        results = []
        for pdb_path in pdb_paths:
            result = pool.apply_async(
                func=process_pdb, 
                args=(pdb_path, pdb_dir), 
                callback=partial_callback
            )
            results.append(result)
        
        for result in results:
            result.wait()
        
        tbar.close()

    print("Done!")

if __name__ == "__main__":
    main()