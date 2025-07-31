import os
import argparse
from pathlib import Path

def rename_chains(input_pdb_path, output_pdb_path):
    """
    Read a PDB file, rename chain A to L and chain B to R, then write to a new file
    """
    with open(input_pdb_path, 'r') as f:
        lines = f.readlines()

    with open(output_pdb_path, 'w') as f:
        for line in lines:
            # Only process ATOM and HETATM records
            if line.startswith(('ATOM', 'HETATM')):
                # Chain ID is located at column 22 (index 21)
                if len(line) >= 22 and line[21] == 'A':
                    line = line[:21] + 'L' + line[22:]
                elif len(line) >= 22 and line[21] == 'B':
                    line = line[:21] + 'R' + line[22:]
            f.write(line)

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description='Rename chains in PDB files from A→L and B→R')
    parser.add_argument('-i', '--input', required=True, help='Input directory containing PDB files')
    parser.add_argument('-o', '--output', default='renamed', help='Output directory (default: renamed)')
    args = parser.parse_args()

    # Input directory path
    input_dir = Path(args.input)
    if not input_dir.is_dir():
        print(f"Error: Input path '{input_dir}' is not a valid directory")
        return

    # Output directory path
    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True)

    # Traverse all PDB files in the input directory (including subdirectories)
    pdb_files = []
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.lower().endswith('.pdb'):
                pdb_files.append(Path(root) / file)

    if not pdb_files:
        print(f"Warning: No PDB files found in input directory '{input_dir}' and its subdirectories")
        return

    # Process each PDB file
    for pdb_file in pdb_files:
        try:
            # Get subfolder name
            relative_path = pdb_file.relative_to(input_dir)
            subfolder_name = relative_path.parent.name if relative_path.parent.name else "root"
            
            # Extract the number after the last underscore in the original filename
            original_filename = pdb_file.stem
            last_underscore_index = original_filename.rfind('_')
            if last_underscore_index != -1 and original_filename[last_underscore_index+1:].isdigit():
                number = original_filename[last_underscore_index+1:]
            else:
                # If no valid number found, use a default (e.g., 000)
                number = "000"
                print(f"Warning: No valid number found in filename '{original_filename}', using default {number}")
            
            # Create output filename: diff_{subfolder name}_{number}.pdb
            output_filename = f"diff_{subfolder_name}_{number}.pdb"
            
            # Output directly to the main output directory without subfolders
            output_path = output_dir / output_filename

            rename_chains(pdb_file, output_path)
            print(f"Processed: {pdb_file} → {output_path}")
        except Exception as e:
            print(f"Error processing file {pdb_file}: {str(e)}")

if __name__ == "__main__":
    main()

