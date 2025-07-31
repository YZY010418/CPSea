import os
import json
import argparse
from pathlib import Path

def extract_epitope_residues(pdb_file, chain_id, residue_ids):
    residues = {}
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                if len(line) >= 22 and line[21] == chain_id:
                    try:
                        res_seq = int(line[22:26].strip())
                        if res_seq in residue_ids:
                            res_name = line[17:20].strip()
                            atom_name = line[12:16].strip()
                            
                            if res_seq not in residues:
                                residues[res_seq] = {
                                    'res_name': res_name,
                                    'atoms': []
                                }
                            
                            residues[res_seq]['atoms'].append(line)
                    except ValueError:
                        continue
    
    return residues

def extract_l_chain(pdb_file):
    l_chain = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                if len(line) >= 22 and line[21] == 'L':
                    l_chain.append(line)
            elif line.startswith('TER'):
                l_chain.append(line)
    
    return l_chain

def write_combined_pdb(l_chain, epitope_residues, output_file):
    with open(output_file, 'w') as f:
        f.writelines(l_chain)
        
        f.write('TER\n')
        
        for res_seq in sorted(epitope_residues.keys()):
            res_data = epitope_residues[res_seq]
            for atom_line in res_data['atoms']:
                if len(atom_line) >= 22:
                    atom_line = atom_line[:21] + 'R' + atom_line[22:]
                f.write(atom_line)
        
        f.write('TER\n')

def get_receptor_chain_ids(pdb_file):
    chain_ids = set()
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')) and len(line) >= 22:
                chain_id = line[21]
                if chain_id != 'L':
                    chain_ids.add(chain_id)
    return chain_ids

def main():
    parser = argparse.ArgumentParser(description='Combine L chain from generated and specific residues from receptors_approved')
    parser.add_argument('--generated', required=True, help='Directory containing generated PDB files')
    parser.add_argument('--receptors', required=True, help='Directory containing receptors_approved PDB files')
    parser.add_argument('--epitopes', required=True, help='JSON file containing epitope residues')
    parser.add_argument('--output', default='reconstructed', help='Output directory (default: reconstructed)')
    args = parser.parse_args()

    generated_dir = Path(args.generated)
    receptors_dir = Path(args.receptors)
    epitopes_file = Path(args.epitopes)
    
    if not generated_dir.is_dir() or not receptors_dir.is_dir() or not epitopes_file.is_file():
        print(f"Error: Input path is not a valid directory or file")
        return

    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True)

    try:
        with open(epitopes_file, 'r') as f:
            epitope_data = json.load(f)
    except Exception as e:
        print(f"Error: Unable to read epitope JSON file: {str(e)}")
        return

    generated_pdb_files = []
    for root, _, files in os.walk(generated_dir):
        for file in files:
            if file.lower().endswith('.pdb'):
                generated_pdb_files.append(Path(root) / file)

    if not generated_pdb_files:
        print(f"Warning: No PDB files found in generated directory or its subdirectories")
        return

    for generated_pdb in generated_pdb_files:
        try:
            pdb_code = generated_pdb.stem[:9].lower().split('_')[1]
            
            if pdb_code not in epitope_data:
                print(f"Warning: No information for PDB code {pdb_code} found in epitope JSON file")
                continue
            
            motif_str = epitope_data[pdb_code]['motif']
            residue_ids = []
            
            for item in motif_str.split('-'):
                if item.startswith('R'):
                    try:
                        residue_ids.append(int(item[1:]))
                    except ValueError:
                        print(f"Warning: Unable to parse residue number: {item}")
                        continue
            
            if not residue_ids:
                print(f"Warning: Parsed residue list is empty for PDB {pdb_code}")
                continue
            
            receptor_pdb = None
            for file in receptors_dir.glob('*.pdb'):
                if file.stem[:4].lower() == pdb_code:
                    receptor_pdb = file
                    break
            
            if receptor_pdb is None:
                print(f"Warning: No receptor file found for PDB code {pdb_code}")
                continue
            
            chain_ids = get_receptor_chain_ids(receptor_pdb)
            target_chain = None
            if 'A' in chain_ids:
                target_chain = 'A'
            elif 'R' in chain_ids:
                target_chain = 'R'
            else:
                print(f"Warning: No chain A or R found (excluding L chain) in receptor PDB {receptor_pdb}, skipping file")
                continue
            
            print(f"Info: Extracting epitope residues from chain {target_chain} of receptor PDB {receptor_pdb}")
            
            relative_path = generated_pdb.parent.relative_to(generated_dir)
            output_subdir = output_dir / relative_path
            output_subdir.mkdir(parents=True, exist_ok=True)
            
            output_pdb_name = f"{generated_pdb.stem}_recon.pdb"
            output_pdb = output_subdir / output_pdb_name
            
            l_chain = extract_l_chain(generated_pdb)
            
            epitope_residues = extract_epitope_residues(receptor_pdb, target_chain, residue_ids)
            
            if not epitope_residues:
                print(f"Warning: No epitope residues found in chain {target_chain}, skipping file")
                continue
            
            write_combined_pdb(l_chain, epitope_residues, output_pdb)
            
            print(f"Generated: {output_pdb} (contains {len(epitope_residues)} epitope residues)")
            
        except Exception as e:
            print(f"Error processing file {generated_pdb}: {str(e)}")

if __name__ == "__main__":
    main()