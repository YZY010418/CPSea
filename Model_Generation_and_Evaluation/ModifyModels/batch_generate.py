import os
import subprocess
import csv
import argparse

def main():
    parser = argparse.ArgumentParser(description='Batch generate peptides using specified parameters.')
    parser.add_argument('--length_csv', required=True, help='Path to the CSV file containing peptide lengths')
    parser.add_argument('--input_dir', required=True, help='Directory containing receptor PDB files')
    parser.add_argument('--pocket_dir', required=True, help='Directory containing pocket JSON files')
    parser.add_argument('--output_dir', required=True, help='Directory to store output results')
    args = parser.parse_args()

    # Read peptide length mapping from CSV
    length_mapping = {}
    with open(args.length_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pdb_name = row['pdb_name']
            length = row['length']
            length_mapping[pdb_name] = length

    # Ensure output directory exists
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Process each PDB file in input directory
    for filename in os.listdir(args.input_dir):
        if filename.endswith('.pdb') and 'processed' not in filename:
            base_name = os.path.splitext(filename)[0].replace('_receptor', '')
            
            if base_name not in length_mapping:
                print(f"Warning: No length found for {base_name}, skipping processing")
                continue
            peptide_length = length_mapping[base_name]
            
            pdb_path = os.path.join(args.input_dir, filename)
            pocket_json = os.path.join(args.pocket_dir, f'{base_name}_complex_pocket.json')
            out_sub_dir = os.path.join(args.output_dir, base_name)
            if not os.path.exists(out_sub_dir):
                os.makedirs(out_sub_dir)
            
            command = [
                'CUDA_VISIBLE_DEVICES=0',
                'python', '-m', 'api.run',
                '--mode', 'codesign',
                '--pdb', pdb_path,
                '--pocket', pocket_json,
                '--out_dir', out_sub_dir,
                '--length_min', peptide_length,
                '--length_max', peptide_length,
                '--n_samples', '100'
            ]
            
            try:
                subprocess.run(' '.join(command), shell=True, check=True)
                print(f"Successfully processed {filename}, output directory: {out_sub_dir}, peptide length: {peptide_length}")
            except subprocess.CalledProcessError as e:
                print(f"Error processing {filename}: {e}")

if __name__ == "__main__":
    main()