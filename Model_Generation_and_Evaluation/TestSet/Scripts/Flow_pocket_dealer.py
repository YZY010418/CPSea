import os
import shutil
import numpy as np
import argparse
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain

DISTANCE_THRESHOLD = 10

def get_non_hydrogen_atoms(structure):
    non_hydrogen_atoms = []
    for atom in structure.get_atoms():
        if atom.element != 'H':
            non_hydrogen_atoms.append(atom)
    return non_hydrogen_atoms

def add_cb(input_array):
    N, CA, C, O = input_array
    b = CA - N
    c = C - CA
    a = np.cross(b, c)
    CB = np.around(-0.58273431*a + 0.56802827*b - 0.54067466*c + CA, 3)
    return CB

def get_cb_atom(residue):
    if residue.get_resname() == 'GLY':
        try:
            N = np.array(residue['N'].get_coord())
            CA = np.array(residue['CA'].get_coord())
            C = np.array(residue['C'].get_coord())
            O = np.array(residue['O'].get_coord())
            CB_coord = add_cb([N, CA, C, O])
            return CB_coord
        except KeyError:
            return None
    else:
        try:
            return np.array(residue['CB'].get_coord())
        except KeyError:
            return None

def filter_pocket(pocket_structure, peptide_atoms):
    residues_to_keep = []
    peptide_cb_atoms = []
    for residue in peptide_atoms:
        cb_atom = get_cb_atom(residue)
        if cb_atom is not None:
            peptide_cb_atoms.append(cb_atom)

    for residue in pocket_structure.get_residues():
        cb_atom = get_cb_atom(residue)
        if cb_atom is not None:
            for peptide_cb in peptide_cb_atoms:
                distance = np.linalg.norm(cb_atom - peptide_cb)
                if distance < DISTANCE_THRESHOLD:
                    residues_to_keep.append(residue)
                    break

    new_pocket_structure = Structure.Structure('new_pocket')
    model = Model.Model(0)
    new_pocket_structure.add(model)
    chain = Chain.Chain('A')
    model.add(chain)
    for residue in residues_to_keep:
        chain.add(residue.copy())

    return new_pocket_structure

def process_folder(input_folder, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for subfolder in os.listdir(input_folder):
        subfolder_path = os.path.join(input_folder, subfolder)
        if os.path.isdir(subfolder_path):
            pocket_path = os.path.join(subfolder_path, 'pocket.pdb')
            peptide_path = os.path.join(subfolder_path, 'peptide.pdb')

            if os.path.exists(pocket_path) and os.path.exists(peptide_path):
                parser = PDBParser()
                pocket_structure = parser.get_structure('pocket', pocket_path)
                peptide_structure = parser.get_structure('peptide', peptide_path)

                peptide_atoms = [residue for residue in peptide_structure.get_residues()]

                new_pocket_structure = filter_pocket(pocket_structure, peptide_atoms)

                output_subfolder = os.path.join(output_folder, subfolder)
                if not os.path.exists(output_subfolder):
                    os.makedirs(output_subfolder)

                io = PDBIO()
                io.set_structure(new_pocket_structure)
                output_pocket_path = os.path.join(output_subfolder, 'pocket.pdb')
                io.save(output_pocket_path)

                output_peptide_path = os.path.join(output_subfolder, 'peptide.pdb')
                shutil.copyfile(peptide_path, output_peptide_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True, help='Input folder path')
    parser.add_argument('-o', required=True, help='Output folder path')
    args = parser.parse_args()
    process_folder(args.i, args.o)