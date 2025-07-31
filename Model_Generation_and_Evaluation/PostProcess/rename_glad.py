import os
import re
import argparse
from pathlib import Path

def process_pdb_file(input_path, output_path):
    """Process a single PDB file"""
    with open(input_path, 'r') as f:
        lines = f.readlines()
    
    chain_l_residues = {}  # Store residue number mapping for chain L
    current_resid = None
    new_resid = 0
    
    # First pass: collect residue number mapping for chain L
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            chain_id = line[21].strip()
            resid = line[22:27].strip()
            
            if chain_id == 'L':
                # If it's a new residue number
                if resid != current_resid:
                    current_resid = resid
                    new_resid += 1
                    chain_l_residues[resid] = new_resid
    
    # Second pass: modify chain ID and residue numbers
    modified_lines = []
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            # Modify chain ID: S -> L
            chain_id = line[21]
            if chain_id == 'S':
                line = line[:21] + 'L' + line[22:]
            
            # Modify residue numbers for chain L
            current_chain_id = line[21].strip()
            resid = line[22:27].strip()
            
            if current_chain_id == 'L' and resid in chain_l_residues:
                new_resid = chain_l_residues[resid]
                # Ensure the new residue number has the correct format
                new_resid_str = f"{new_resid:>5}"
                line = line[:22] + new_resid_str + line[27:]
        
        modified_lines.append(line)
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Write the modified file
    with open(output_path, 'w') as f:
        f.writelines(modified_lines)

def rename_file_according_to_rule(original_filename):
    """
    Rename file according to rules: prefix + content before first underscore + number after second underscore
    For example: input_abc_123.pdb -> glad_input_123.pdb
    """
    # Separate filename and extension
    base_name, ext = os.path.splitext(original_filename)
    
    # Split by underscore
    parts = base_name.split('_')
    
    # Process split results
    if len(parts) >= 3:
        # Take the first part and the last part (number part), add prefix
        new_base = f"glad_{parts[0]}_{parts[-1]}"
    else:
        # If the format requirements are not met, use the original filename with prefix
        new_base = f"glad_{base_name}"
    
    return f"{new_base}{ext}"

def main():
    parser = argparse.ArgumentParser(description='Process PDB files, rename chain S to L, and adjust L chain numbering starting from 1')
    parser.add_argument('-i', '--input', required=True, help='Input folder path')
    parser.add_argument('-o', '--output', default='glad_renamed', help='Output folder path, default is glad_renamed')
    args = parser.parse_args()
    
    input_dir = args.input
    output_dir = args.output
    
    # Check if input folder exists
    if not os.path.exists(input_dir):
        print(f"Error: Input folder '{input_dir}' does not exist")
        return
    
    # Traverse all PDB files in the input folder
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.lower().endswith('.pdb'):
                input_path = os.path.join(root, file)
                
                # Build output path without preserving subdirectory structure, apply new naming rules
                new_filename = rename_file_according_to_rule(file)
                output_path = os.path.join(output_dir, new_filename)
                
                try:
                    process_pdb_file(input_path, output_path)
                    print(f"Processed: {input_path} -> {output_path}")
                except Exception as e:
                    print(f"Error processing file {input_path}: {str(e)}")

if __name__ == "__main__":
    main()
