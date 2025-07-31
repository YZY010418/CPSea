import os
import numpy as np
from Bio.PDB import PDBParser
import json
import argparse

def get_ca_coords(structure):
    ca_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    ca_coords.append(residue['CA'].coord)
    return ca_coords

def get_resid_id(residue):
    chain_id = residue.get_parent().get_id()
    res_id = residue.get_id()[1]
    return f"{chain_id}{res_id}"

def add_cb(input_array):
    N, CA, C, O = input_array
    b = CA - N
    c = C - CA
    a = np.cross(b, c)
    CB = np.around(-0.58273431*a + 0.56802827*b - 0.54067466*c + CA, 3)
    return CB

def get_cb_coords(residue):
    if residue.get_resname() == 'GLY':
        if 'N' in residue and 'CA' in residue and 'C' in residue and 'O' in residue:
            N = residue['N'].coord
            CA = residue['CA'].coord
            C = residue['C'].coord
            O = residue['O'].coord
            return add_cb([N, CA, C, O])
        else:
            return None
    elif 'CB' in residue:
        return residue['CB'].coord
    else:
        return None

def find_motifs(ligand_coords_cb, rec_residues, pocket_cutoff=10):
    motif_residues = []
    for i in ligand_coords_cb:
        if i is not None:
            for j in rec_residues:
                cb_coord = get_cb_coords(j)
                if cb_coord is not None:
                    dist = np.linalg.norm(cb_coord - i)
                    if dist <= pocket_cutoff and j not in motif_residues:
                        motif_residues.append(j)

    motif_ids = "-".join(sorted([get_resid_id(res) for res in motif_residues]))
    return motif_ids

def process_folder(folder_path):
    result = {}
    subfolders = [f for f in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, f))]
    parser = PDBParser()

    for subfolder in subfolders:
        subfolder_path = os.path.join(folder_path, subfolder)
        peptide_path = os.path.join(subfolder_path, 'peptide.pdb')
        pocket_path = os.path.join(subfolder_path, 'pocket.pdb')

        if os.path.exists(peptide_path) and os.path.exists(pocket_path):
            peptide_structure = parser.get_structure('peptide', peptide_path)
            pocket_structure = parser.get_structure('pocket', pocket_path)

            ligand_coords_cb = [get_cb_coords(res) for model in peptide_structure for chain in model for res in chain]
            rec_residues = [res for model in pocket_structure for chain in model for res in chain]

            motif = find_motifs(ligand_coords_cb, rec_residues)
            result[f"{subfolder}"] = {
                "motif": motif,
            }

    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--output', required=True)
    args = parser.parse_args()

    result = process_folder(args.input)
    with open(args.output, 'w') as f:
        json.dump(result, f, indent=4)