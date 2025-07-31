import os
import gzip
import argparse
from rdkit import Chem
from rdkit.Chem.rdmolfiles import MolFromPDBBlock
from rdkit.Chem.rdMolDescriptors import CalcTPSA, CalcCrippenDescriptors
from rdkit.Chem import SanitizeFlags
from io import StringIO
import csv
from tqdm import tqdm
import warnings
from multiprocessing import Pool
from functools import partial

warnings.filterwarnings("ignore")

def filter_chain(pdb_string, chain_id):
    lines = pdb_string.split('\n')
    chain_lines = []
    for line in lines:
        if line.startswith('ATOM') and line[21:22] == chain_id:
            chain_lines.append(line)
    return '\n'.join(chain_lines)

def sanitize_mol(mol):
    try:
        Chem.SanitizeMol(
            mol, 
            sanitizeOps=(
                SanitizeFlags.SANITIZE_ADJUSTHS |
                SanitizeFlags.SANITIZE_CLEANUP |
                SanitizeFlags.SANITIZE_FINDRADICALS |
                SanitizeFlags.SANITIZE_KEKULIZE |
                SanitizeFlags.SANITIZE_SETCONJUGATION |
                SanitizeFlags.SANITIZE_SETHYBRIDIZATION |
                SanitizeFlags.SANITIZE_SYMMRINGS
            ),
            catchErrors=True
        )
        return mol
    except ValueError as e:
        print(f"Molecule sanitization error: {e}")
        return None

def calculate_tpsa_and_logp(pdb_string, chain_id):

    chain_pdb = filter_chain(pdb_string, chain_id)
    
    pepmol = MolFromPDBBlock(chain_pdb, removeHs=False, sanitize=False)
    
    if pepmol is None:
        print(f"Failed to read PDB block for chain {chain_id}.")
        return None, None, None, None
    
    pepmol = sanitize_mol(pepmol)
    
    if pepmol is None:
        print(f"Failed to sanitize molecule for chain {chain_id}.")
        return None, None, None, None

    try:
        tpsa = CalcTPSA(pepmol)
        logP, _ = CalcCrippenDescriptors(pepmol)
        heavy_atom_count = pepmol.GetNumHeavyAtoms()

        if heavy_atom_count == 0:
            tpsa_over_atom = 0
        else:
            tpsa_over_atom = tpsa / heavy_atom_count

        return tpsa, logP, tpsa_over_atom, heavy_atom_count
    except Exception as e:
        print(f"Error calculating descriptors for chain {chain_id}: {e}")
        return None, None, None, None

def read_pdb_file(pdb_path):
    if pdb_path.endswith('.gz'):
        with gzip.open(pdb_path, 'rt') as f:
            return f.read()
    else:
        with open(pdb_path, 'r') as f:
            return f.read()

def process_one(pdb_path, csv_file_path, chain_id):
    try:
        teststring = read_pdb_file(pdb_path)
        tpsa, logP, tpsa_over_atom, heavy_atom_count = calculate_tpsa_and_logp(teststring, chain_id)

        if tpsa is not None:
            g = os.path.basename(pdb_path)
            with open(csv_file_path, 'a', newline='') as csvfile:
                fieldnames = ['input', 'chain', 'TPSA', 'logP', 'tpsa_over_atom', 'heavy_atom_count']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writerow({
                    'input': g,
                    'chain': chain_id,
                    'TPSA': tpsa,
                    'logP': logP,
                    'tpsa_over_atom': tpsa_over_atom,
                    'heavy_atom_count': heavy_atom_count
                })
        return True
    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")
        return False

def callback(result, tbar):
    tbar.update(1)

def main():
    parser = argparse.ArgumentParser(description='Calculate TPSA and logP for specified chain in PDB files')
    parser.add_argument('-i', '--input_list', required=True, help='Path to the list of PDB file paths')
    parser.add_argument('-o', '--output_csv', required=True, help='Path to the output CSV file')
    parser.add_argument('--chain_id', required=True, help='Chain identifier to process (e.g., A, B, L)')
    parser.add_argument('-c', '--cores', type=int, default=1, help='Number of CPU cores to use')
    args = parser.parse_args()

    input_list_path = args.input_list
    csv_file_path = args.output_csv
    chain_id = args.chain_id
    N_CPU = args.cores

    if not os.path.exists(input_list_path):
        print(f"Error: Input list file '{input_list_path}' does not exist.")
        return

    with open(csv_file_path, 'w', newline='') as csvfile:
        fieldnames = ['input', 'chain', 'TPSA', 'logP', 'tpsa_over_atom', 'heavy_atom_count']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

    try:
        with open(input_list_path, 'r') as f:
            pdb_paths = [line.strip() for line in f if line.strip()]
            total_pdbs = len(pdb_paths)
            print(f"Found {total_pdbs} PDB files in the input list.")
            print(f"Processing chain: {chain_id}")
            
            if total_pdbs == 0:
                print("No PDB files found in the input list.")
                return

            tbar = tqdm(total=total_pdbs)
            partial_callback = partial(callback, tbar=tbar)
            
            with Pool(processes=N_CPU) as pool:
                results = []
                for pdb_path in pdb_paths:
                    if not os.path.exists(pdb_path):
                        print(f"Warning: PDB file '{pdb_path}' does not exist. Skipping...")
                        tbar.update(1)
                        continue
                    
                    result = pool.apply_async(
                        func=process_one, 
                        args=(pdb_path, csv_file_path, chain_id), 
                        callback=partial_callback
                    )
                    results.append(result)
                
                for result in results:
                    result.wait()
                
                tbar.close()
                
            print(f"Processing completed. Results saved to {csv_file_path}")
            
    except Exception as e:
        print(f"Error reading input list: {e}")

if __name__ == '__main__':
    main()