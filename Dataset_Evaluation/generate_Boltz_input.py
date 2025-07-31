import argparse
import os
import re
from typing import Dict, List, Tuple

def build_aa_mapping() -> Dict[str, str]:
    return {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        'SEC': 'U',
        'PYL': 'O'
    }

def extract_l_chain_sequence(pdb_path: str, aa_mapping: Dict[str, str]) -> Tuple[str, List[str]]:
    residues = {}
    unrecognized = []
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM  '):
                chain_id = line[21]
                if chain_id == 'L':
                    res_name = line[17:20].strip()
                    res_seq = int(line[22:26].strip())
                    
                    if res_seq not in residues:
                        residues[res_seq] = res_name
                        if res_name not in aa_mapping:
                            unrecognized.append(f"{res_name} (position {res_seq})")
    
    sorted_residues = sorted(residues.items(), key=lambda x: x[0])
    sequence = ''.join([aa_mapping[res] for _, res in sorted_residues if res in aa_mapping])
    
    return sequence, unrecognized

def replace_in_template(template_content: str, sequence: str) -> str:
    new_content = template_content.replace(
        'sequence: SKGNDEISAMGRAVD', 
        f'sequence: {sequence}'
    )
    
    seq_length = len(sequence)
    new_content = re.sub(
        r'(atom2: \[L, )\d+(, "CG"\])',
        f'\\g<1>{seq_length}\\g<2>',
        new_content
    )
    
    return new_content

def main():
    parser = argparse.ArgumentParser(description='Batch generate new YAML files from PDB files and template YAML')
    parser.add_argument('-i', required=True, help='Input folder path containing PDB files')
    parser.add_argument('-o', required=True, help='Output folder path for YAML files')
    parser.add_argument('-t', required=True, help='Path to the template YAML file')
    args = parser.parse_args()
    
    aa_mapping = build_aa_mapping()
    
    os.makedirs(args.o, exist_ok=True)
    
    try:
        with open(args.t, 'r') as f:
            template_content = f.read()
    except Exception as e:
        print(f"Failed to read template YAML: {e}")
        return
    
    pdb_files = [f for f in os.listdir(args.i) if f.lower().endswith('.pdb')]
    if not pdb_files:
        print("No PDB files found in input folder")
        return
    
    for pdb_file in pdb_files:
        pdb_path = os.path.join(args.i, pdb_file)
        print(f"Processing file: {pdb_file}")
        
        sequence, unrecognized = extract_l_chain_sequence(pdb_path, aa_mapping)
        
        if not sequence:
            print(f"  Warning: No valid sequence extracted from chain L, skipping file")
            continue
        if unrecognized:
            print(f"  Warning: Unrecognized residues found (skipped): {', '.join(unrecognized)}")
        
        try:
            new_content = replace_in_template(template_content, sequence)
        except Exception as e:
            print(f"  Failed to generate YAML: {e}")
            continue
        
        output_filename = os.path.splitext(pdb_file)[0] + '.yaml'
        output_path = os.path.join(args.o, output_filename)
        try:
            with open(output_path, 'w') as f:
                f.write(new_content)
            print(f"  Saved: {output_filename} (sequence length: {len(sequence)})")
        except Exception as e:
            print(f"  Failed to save YAML: {e}")
            continue
    
    print("Processing completed")

if __name__ == '__main__':
    main()