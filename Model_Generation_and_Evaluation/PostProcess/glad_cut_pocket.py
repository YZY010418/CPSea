import os
import json
import argparse
from pathlib import Path

def parse_pdb_atom_line(line):
    try:
        atom_serial = int(line[6:11].strip())
        atom_name = line[12:16].strip()
        residue_name = line[17:20].strip()
        chain_id = line[21].strip()
        residue_seq = int(line[22:26].strip())
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        return {
            'atom_serial': atom_serial,
            'atom_name': atom_name,
            'residue_name': residue_name,
            'chain_id': chain_id,
            'residue_seq': residue_seq,
            'x': x,
            'y': y,
            'z': z,
            'original_line': line
        }
    except (ValueError, IndexError):
        return None

def read_json_residues(json_file):
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
            r_residues = []
            for item in data:
                if item[0] == 'R':
                    try:
                        residue_num = int(item[1][0])
                        r_residues.append(residue_num)
                    except (ValueError, IndexError):
                        continue
            return r_residues
    except Exception as e:
        print(f"Error reading {json_file}: {e}")
        return []

def extract_chains(renamed_pdb_file, receptor_pdb_file, json_file, output_file):
    if not os.path.exists(renamed_pdb_file):
        print(f"Renamed PDB file not found: {renamed_pdb_file}")
        return
    
    if not os.path.exists(receptor_pdb_file):
        print(f"Receptor PDB file not found: {receptor_pdb_file}")
        return
    
    if not os.path.exists(json_file):
        print(f"JSON file not found: {json_file}")
        return
    
    r_residues = read_json_residues(json_file)
    if not r_residues:
        print(f"No residues found in JSON file: {json_file}")
        return
    
    l_atoms = []
    with open(renamed_pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                atom = parse_pdb_atom_line(line)
                if atom and atom['chain_id'] == 'L':
                    l_atoms.append(atom)
    
    r_atoms = []
    with open(receptor_pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                atom = parse_pdb_atom_line(line)
                if atom and atom['chain_id'] == 'R' and atom['residue_seq'] in r_residues:
                    r_atoms.append(atom)
    
    with open(output_file, 'w') as f:
        f.write(f"HEADER    Extracted from {os.path.basename(renamed_pdb_file)} and {os.path.basename(receptor_pdb_file)}\n")
        f.write(f"REMARK    R chain residues extracted based on {os.path.basename(json_file)}\n")
        
        for atom in r_atoms:
            f.write(atom['original_line'])
        
        for atom in l_atoms:
            f.write(atom['original_line'])
        
        f.write("END\n")
    
    print(f"Extracted {len(r_atoms)} R-chain atoms and {len(l_atoms)} L-chain atoms to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Extract specific residues from PDB files')
    parser.add_argument('--renamed', default='renamed', help='Path to the renamed PDB folder')
    parser.add_argument('--receptor', default='receptor', help='Path to the receptor PDB folder')
    parser.add_argument('--json', default='../receptor_pocket_0523', help='Path to the JSON folder')
    parser.add_argument('--output', default='glad_cutpocket', help='Path to the output folder')
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    
    for root, dirs, files in os.walk(args.renamed):
        for file in files:
            if file.endswith('.pdb'):
                pdb_file = os.path.join(root, file)
                base_name = os.path.splitext(file)[0].split('_')[1]
                json_file = os.path.join(args.json, f"{base_name}_complex_pocket.json")
                relative_path = os.path.relpath(root, args.renamed)
                output_subfolder = os.path.join(args.output, relative_path)
                os.makedirs(output_subfolder, exist_ok=True)
                output_pdb = os.path.join(output_subfolder, file)
                
                receptor_pdb_files = [f for f in os.listdir(args.receptor) if f.endswith('.pdb') and 'processed' not in f and base_name in f]
                if not receptor_pdb_files:
                    print(f"No matching receptor PDB file found for {pdb_file}")
                    continue
                receptor_pdb_file = os.path.join(args.receptor, receptor_pdb_files[0])
                
                extract_chains(pdb_file, receptor_pdb_file, json_file, output_pdb)

if __name__ == "__main__":
    main()