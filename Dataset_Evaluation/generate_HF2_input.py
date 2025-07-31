import csv
import os
import argparse
import re

# Three-letter to one-letter amino acid mapping
AA_THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def extract_l_chain_sequence(pdb_file_path):
    """
    Extract amino acid sequence of L chain from PDB file
    
    Parameters:
    pdb_file_path: Path to the PDB file
    
    Returns:
    Extracted L chain sequence (string)
    """
    sequence = []
    prev_residue_num = None
    
    with open(pdb_file_path, 'r', encoding='utf-8', errors='ignore') as pdb_file:
        for line in pdb_file:
            # Process only ATOM records with CA atoms (alpha carbon)
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                chain_id = line[21].strip()
                if chain_id == 'L':
                    res_name = line[17:20].strip()
                    res_num = line[22:26].strip()
                    
                    # Record each residue only once
                    if res_num != prev_residue_num:
                        if res_name in AA_THREE_TO_ONE:
                            sequence.append(AA_THREE_TO_ONE[res_name])
                            prev_residue_num = res_num
                        else:
                            print(f"Warning: Unknown amino acid {res_name} found in {pdb_file_path}, skipped")
    
    return ''.join(sequence)

def generate_a3m_file(file_base_name, sequence, output_dir, sort=False):
    """
    Generate a3m file from sequence
    
    Parameters:
    file_base_name: Base filename
    sequence: Amino acid sequence
    output_dir: Output directory
    sort: Whether to save in subfolders by sequence length
    """
    seq_len = len(sequence)
    
    # Create a3m content with tab-separated numbers
    a3m_content = f"#{seq_len}\t1\n"
    a3m_content += f">101\n{sequence}\n"
    a3m_content += f">101\n{sequence}\n"
    
    a3m_file_name = f"{file_base_name}.a3m"
    
    # Create length-specific subfolder if sorting is enabled
    if sort:
        subdir = f"length_{seq_len}"
        output_dir = os.path.join(output_dir, subdir)
        os.makedirs(output_dir, exist_ok=True)
    
    a3m_file_path = os.path.join(output_dir, a3m_file_name)
    
    with open(a3m_file_path, 'w', encoding='utf-8') as a3mfile:
        a3mfile.write(a3m_content)
    
    return a3m_file_path

def process_csv(csv_file_path, output_dir, sort=False):
    """Process CSV file"""
    with open(csv_file_path, 'r', newline='', encoding='utf-8') as csvfile:
        csv_reader = csv.reader(csvfile)
        
        for row_num, row in enumerate(csv_reader, 1):
            if len(row) < 2:
                print(f"Warning: Incomplete data in row {row_num}, skipped")
                continue
            
            file_base_name = row[0].strip()
            sequence = row[2].strip()
            
            try:
                a3m_file_path = generate_a3m_file(file_base_name, sequence, output_dir, sort)
                print(f"Generated: {a3m_file_path}")
            except AssertionError as e:
                print(f"Warning: Error processing {file_base_name} - {str(e)}, skipped")

def process_pdb_folder(folder_path, output_dir, sort=False):
    """Process folder containing PDB files"""
    for file_name in os.listdir(folder_path):
        if file_name.lower().endswith('.pdb'):
            pdb_file_path = os.path.join(folder_path, file_name)
            file_base_name = os.path.splitext(file_name)[0]
            
            sequence = extract_l_chain_sequence(pdb_file_path)
            
            if not sequence:
                print(f"Warning: No L chain sequence found in {file_name}, skipped")
                continue
            
            try:
                a3m_file_path = generate_a3m_file(file_base_name, sequence, output_dir, sort)
                print(f"Generated: {a3m_file_path}")
            except AssertionError as e:
                print(f"Warning: Error processing {file_name} - {str(e)}, skipped")

def main():
    parser = argparse.ArgumentParser(description='Convert CSV file or folder with PDB files to a3m format')
    parser.add_argument('-i','input', help='Path to input CSV file or folder containing PDB files')
    parser.add_argument('-o','output_dir', help='Directory for output a3m files')
    parser.add_argument('--sort', action='store_true', help='Save in subfolders by sequence length (folder name format: length_sequence length)')
    
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    if os.path.isfile(args.input) and args.input.lower().endswith('.csv'):
        print(f"Detected CSV file, starting processing: {args.input}")
        process_csv(args.input, args.output_dir, args.sort)
    elif os.path.isdir(args.input):
        print(f"Detected folder, starting processing PDB files: {args.input}")
        process_pdb_folder(args.input, args.output_dir, args.sort)
    else:
        print("Error: Input must be a CSV file or folder containing PDB files")
        return
    
    print("Processing completed!")

if __name__ == "__main__":
    main()