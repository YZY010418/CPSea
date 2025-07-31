import argparse
from tqdm import tqdm
import csv
from Bio.PDB import PDBParser, PPBuilder
from multiprocessing import Pool
import os
import glob
import gzip
import shutil

HYDROPHOBIC_RESIDUES = {'V', 'I', 'L', 'M', 'F', 'W', 'C'}

def calculate_gravy_and_hydrophobic_ratio(sequence):
    hydropathy_values = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }
    
    total_hydropathy = sum(hydropathy_values.get(aa, 0) for aa in sequence)
    gravy = total_hydropathy / len(sequence) if sequence else 0
    
    hydrophobic_count = sum(1 for aa in sequence if aa in HYDROPHOBIC_RESIDUES)
    hydrophobic_ratio = hydrophobic_count / len(sequence) if sequence else 0
    
    return gravy, hydrophobic_ratio

def process_pdb(args):

    pdb_file, chain_id, temp_dir = args
    pid = os.getpid()
    file_id = os.path.basename(pdb_file).split('.')[0]
    
    try:

        is_gzipped = pdb_file.lower().endswith('.gz')
        local_pdb_file = pdb_file
        
        if is_gzipped:

            proc_temp_dir = os.path.join(temp_dir, f"temp_{pid}")
            os.makedirs(proc_temp_dir, exist_ok=True)
            
            local_pdb_file = os.path.join(proc_temp_dir, f"{file_id}.pdb")
            with gzip.open(pdb_file, 'rb') as f_in:
                with open(local_pdb_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("pdb", local_pdb_file)
        
        chain_sequence = ""
        for model in structure:
            if chain_id in model:
                chain = model[chain_id]
                ppb = PPBuilder()
                polypeptides = ppb.build_peptides(chain)
                
                for pp in polypeptides:
                    chain_sequence += str(pp.get_sequence())
                break
        
        if is_gzipped and os.path.exists(local_pdb_file):
            os.remove(local_pdb_file)
            try:
                os.rmdir(proc_temp_dir)
            except OSError:
                pass
        
        if not chain_sequence:
            print(f"Warning: {pdb_file} doesn't have chain {chain_id}")
            return False
        
        gravy, hydrophobic_ratio = calculate_gravy_and_hydrophobic_ratio(chain_sequence)
        
        temp_csv = os.path.join(temp_dir, f"temp_{pid}.csv")
        with open(temp_csv, 'a', newline='') as csvfile:
            writer = csv.DictWriter(
                csvfile, 
                fieldnames=['pdb_file', 'chain', 'sequence', 'length', 'GRAVY', 'hydrophobic_ratio']
            )
            if os.path.getsize(temp_csv) == 0:
                writer.writeheader()
            writer.writerow({
                'pdb_file': pdb_file,
                'chain': chain_id,
                'sequence': chain_sequence,
                'length': len(chain_sequence),
                'GRAVY': gravy,
                'hydrophobic_ratio': hydrophobic_ratio
            })
        return True
    except Exception as e:
        print(f"Error occurs when processing {pdb_file}: {e}")
        return False

def merge_temp_files(temp_dir, output_file):
    temp_files = glob.glob(os.path.join(temp_dir, "temp_*.csv"))
    
    if not temp_files:
        print("Warning: no temp files to merge")
        return
    
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.DictWriter(
            outfile, 
            fieldnames=['pdb_file', 'chain', 'sequence', 'length', 'GRAVY', 'hydrophobic_ratio']
        )
        writer.writeheader()
        
        for temp_file in temp_files:
            try:
                with open(temp_file, 'r') as infile:
                    reader = csv.DictReader(infile)
                    for row in reader:
                        writer.writerow(row)
                os.remove(temp_file)
            except Exception as e:
                print(f"Error when merging {temp_file}: {e}")
    
    try:
        os.rmdir(temp_dir)
    except OSError:
        pass  

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate GRAVY and hydrophobic ratio for a specific chain from PDB files')
    parser.add_argument('-i', '--input', required=True, help='List of PDB file paths (one path per line)')
    parser.add_argument('-o', '--output', required=True, help='Output CSV file')
    parser.add_argument('--chain_id', required=True, help='Chain ID to analyze (e.g., "A")')
    parser.add_argument('-c', '--cores', type=int, default=1, help='Number of CPU processes')
    
    args = parser.parse_args()
    
    # Create a temporary directory
    output_dir = os.path.dirname(os.path.abspath(args.output))
    temp_dir = os.path.join(output_dir, "temp_pdb_processing")
    os.makedirs(temp_dir, exist_ok=True)
    
    # Read the list of PDB files
    with open(args.input) as f:
        pdb_files = [line.strip() for line in f if line.strip()]
    
    print(f"Found {len(pdb_files)} PDB files")
    print(f"Analyzing chain: {args.chain_id}")
    print(f"Output file: {args.output}")
    print(f"Using {args.cores} CPU processes")
    
    # Prepare arguments for the process pool
    task_args = [(pdb_file, args.chain_id, temp_dir) for pdb_file in pdb_files]
    
    # Use a process pool to process PDB files
    with Pool(processes=args.cores) as pool:
        results = list(tqdm(
            pool.imap_unordered(process_pdb, task_args),
            total=len(pdb_files),
            desc="Processing PDB files"
        ))
    
    # Merge all temporary files
    merge_temp_files(temp_dir, args.output)
    
    success_count = sum(results)
    print(f"Processing completed: {success_count}/{len(pdb_files)} files succeeded")