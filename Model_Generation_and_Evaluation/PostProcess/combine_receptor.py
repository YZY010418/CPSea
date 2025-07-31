import os
import math
from Bio import PDB
import argparse
from tqdm.auto import tqdm

def translate_structure(structure, translation_vector):
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.coord += translation_vector
    return structure

def process_relaxed_file(parser, file_path):
    file_name = os.path.basename(file_path)
    if len(file_name) < 9:
        print(f"Warning: {file_path} has a too short filename, cannot extract identifier, skipped")
        return None, None
    identifier = file_name[5:9]
    structure = parser.get_structure(identifier, file_path)
    
    r_chain = None
    for model in structure:
        for chain in model:
            if chain.id == 'R':
                r_chain = chain
                break
        if r_chain:
            break
    if not r_chain:
        print(f"Warning: No R chain found in {file_path}, skipped")
        return None, None
    
    ca_info_list = []
    for residue in r_chain:
        if 'CA' in residue:
            ca_atom = residue['CA']
            ca_info = {
                'coord': ca_atom.coord.copy(),
                'resname': residue.get_resname(),
                'resid': residue.get_id()[1]
            }
            ca_info_list.append(ca_info)
            if len(ca_info_list) == 4:
                break
    if len(ca_info_list) < 4:
        print(f"Warning: R chain in {file_path} has fewer than 4 residues with CA atoms, skipped")
        return None, None
    
    l_chain = None
    for model in structure:
        for chain in model:
            if chain.id == 'L':
                l_chain = chain
                break
        if l_chain:
            break
    if not l_chain:
        print(f"Warning: No L chain found in {file_path}, skipped")
        return None, None
    
    return (ca_info_list, l_chain), identifier

def process_receptors_file(parser, file_path):
    identifier = os.path.basename(file_path)[:4]
    structure = parser.get_structure(identifier, file_path)
    
    chains = []
    for model in structure:
        chains.extend(model.get_list())
    if len(chains) != 1:
        print(f"Warning: {file_path} contains {len(chains)} chains (expected 1), skipped")
        return None, None
    
    return structure, identifier

def calculate_translation(relaxed_ca, receptor_structure):
    receptor_chain = next(next(receptor_structure.get_models()).get_chains())
    
    receptor_ca_coords = []
    for target in relaxed_ca:
        resid = target['resid']
        if resid not in receptor_chain:
            print(f"Warning: Residue {resid} not found in receptor, skipped")
            return None
        residue = receptor_chain[resid]
        if 'CA' not in residue:
            print(f"Warning: Residue {resid} in receptor lacks CA atom, skipped")
            return None
        receptor_ca_coords.append(residue['CA'].coord)
    
    sum_x, sum_y, sum_z = 0.0, 0.0, 0.0
    count = len(relaxed_ca)
    
    for target, receptor_ca in zip(relaxed_ca, receptor_ca_coords):
        sum_x += target['coord'][0] - receptor_ca[0]
        sum_y += target['coord'][1] - receptor_ca[1]
        sum_z += target['coord'][2] - receptor_ca[2]
    
    return [sum_x/count, sum_y/count, sum_z/count]

def merge_chains(l_chain, receptor_chain):
    new_structure = PDB.Structure.Structure('merged')
    new_model = PDB.Model.Model(0)
    new_structure.add(new_model)
    
    new_model.add(l_chain)
    
    receptor_chain.id = 'R'
    new_model.add(receptor_chain)
    
    return new_structure

def main():
    parser = argparse.ArgumentParser(description='Align and merge PDB structures')
    parser.add_argument('--relaxed', required=True, help='Relaxed folder (contains PDB files with R and L chains)')
    parser.add_argument('--receptors', required=True, help='Receptors folder (contains PDB files with single chain)')
    parser.add_argument('--output', required=True, help='Output folder (save merged structures)')
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    
    pdb_parser = PDB.PDBParser(QUIET=True)
    pdb_writer = PDB.PDBIO()
    
    receptor_files = [f for f in os.listdir(args.receptors) if f.endswith('.pdb')]
    receptor_data = {}
    print("Processing PDB files in receptors folder...")
    for file in tqdm(receptor_files):
        file_path = os.path.join(args.receptors, file)
        structure, identifier = process_receptors_file(pdb_parser, file_path)
        if structure:
            receptor_data[identifier] = (structure, file)
    
    relaxed_files = [f for f in os.listdir(args.relaxed) if f.endswith('.pdb')]
    print("Processing PDB files in relaxed folder and matching...")
    for file in tqdm(relaxed_files):
        file_path = os.path.join(args.relaxed, file)
        data, identifier = process_relaxed_file(pdb_parser, file_path)
        
        if not data:
            continue
        
        if identifier not in receptor_data:
            print(f"Warning: No matching receptor found for {file} (identifier: {identifier}), skipped")
            continue
        
        receptor_struct, receptor_file = receptor_data[identifier]
        relaxed_ca, l_chain = data
        
        translation = calculate_translation(relaxed_ca, receptor_struct)
        if translation is None:
            continue
        
        receptor_copy = receptor_struct.copy()
        translated_struct = translate_structure(receptor_copy, translation)
        
        receptor_chain = next(next(translated_struct.get_models()).get_chains())
        
        merged_struct = merge_chains(l_chain, receptor_chain)
        
        output_file = f"{os.path.splitext(file)[0]}_recon.pdb"
        output_path = os.path.join(args.output, output_file)
        pdb_writer.set_structure(merged_struct)
        pdb_writer.save(output_path)

if __name__ == "__main__":
    main()