import os
import re
import csv
import sys
import time
import gzip
import shutil
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from functools import partial
import warnings
import xml.etree.ElementTree as ET
from vina import Vina
from Bio import PDB
import subprocess
import uuid

warnings.filterwarnings("ignore")

def process_pdb_file(pdb_path):
    """Process a single PDB file (either .pdb or .pdb.gz), run analysis and return results"""
    try:
        # Create a unique temporary directory for each process
        tmp_dir = os.path.abspath(os.path.join("./tmp_vina", str(uuid.uuid4())))
        os.makedirs(tmp_dir, exist_ok=True)
        
        # Get the file name and base name
        file_name = os.path.basename(pdb_path)
        
        # Determine if file is compressed and get base name
        if file_name.endswith('.pdb.gz'):
            is_compressed = True
            base_name = os.path.splitext(os.path.splitext(file_name)[0])[0]
        elif file_name.endswith('.pdb'):
            is_compressed = False
            base_name = os.path.splitext(file_name)[0]
        else:
            raise ValueError(f"Unsupported file format: {file_name}. Use .pdb or .pdb.gz")
        
        # Handle file decompression if needed
        temp_pdb_file = os.path.join(tmp_dir, f"{base_name}.pdb")
        if is_compressed:
            with gzip.open(pdb_path, 'rb') as f_in:
                with open(temp_pdb_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            shutil.copy2(pdb_path, temp_pdb_file)
        
        parser = PDB.PDBParser(QUIET=True)
        io = PDB.PDBIO()
        structure = parser.get_structure(base_name, temp_pdb_file)
        
        L_path = None
        R_path = None
        for model in structure:
            for chain in model:
                if chain.id == 'L':
                    L_structure = PDB.Structure.Structure(base_name)
                    L_model = PDB.Model.Model(0)
                    L_structure.add(L_model)
                    L_model.add(chain.copy())
                    
                    L_path = os.path.join(tmp_dir, f"{base_name}_L.pdb")
                    io.set_structure(L_structure)
                    io.save(str(L_path))
                elif chain.id == 'R':
                    R_structure = PDB.Structure.Structure(base_name)
                    R_model = PDB.Model.Model(0)
                    R_structure.add(R_model)
                    R_model.add(chain.copy())
                    
                    R_path = os.path.join(tmp_dir, f"{base_name}_R.pdb")
                    io.set_structure(R_structure)
                    io.save(str(R_path))

        if not L_path or not R_path:
            raise ValueError(f"Missing 'L' or 'R' chains in {pdb_path}")

        L_pdbqt_path = os.path.join(tmp_dir, f"{base_name}CP.pdbqt")
        R_pdbqt_path = os.path.join(tmp_dir, f"{base_name}protein.pdbqt")
        ligand_cmd = [
            '/data_hdd/home/yangziyi/Tools/mgltools_x86_64Linux2_1.5.6/bin/pythonsh',
            '/data_hdd/home/yangziyi/Tools/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py',
            '-l', L_path,
            '-o', L_pdbqt_path
        ]
        try:
            subprocess.run(ligand_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"Error running prepare_ligand4.py for {pdb_path}: {e.stderr.decode('utf-8')}", file=sys.stderr)
            return os.path.basename(pdb_path), None
        
        # Prepare the receptor file
        receptor_cmd = [
            '/data_hdd/home/yangziyi/Tools/mgltools_x86_64Linux2_1.5.6/bin/pythonsh',
            '/data_hdd/home/yangziyi/Tools/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py',
            '-r', R_path,
            '-o', R_pdbqt_path
        ]
        try:
            subprocess.run(receptor_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"Error running prepare_receptor4.py for {pdb_path}: {e.stderr.decode('utf-8')}", file=sys.stderr)
            return os.path.basename(pdb_path), None

        vina_cmd = [
            '/data_hdd/home/yangziyi/Tools/autodock_vina_1_1_2_linux_x86/bin/vina',
            '--receptor', R_pdbqt_path,
            '--ligand', L_pdbqt_path,
            '--cpu', '1',
            '--score_only'
        ]
        try:
            result = subprocess.run(vina_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running vina for {pdb_path}: {e.stderr.decode('utf-8')}", file=sys.stderr)
            return os.path.basename(pdb_path), None
        
        # Extract the Affinity value
        affinity_match = re.search(r'Affinity:\s+([-\d.]+)', result.stdout)
        affinity = affinity_match.group(1) if affinity_match else None

        # Clean up temporary files
        for file_path in [L_path, R_path, L_pdbqt_path, R_pdbqt_path, temp_pdb_file]:
            if os.path.exists(file_path):
                os.remove(file_path)

        return base_name, affinity
    except Exception as e:
        print(f"Error processing {pdb_path}: {e}", file=sys.stderr)
        return os.path.basename(pdb_path), None
    finally:
        # Clean up the temporary directory
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir, ignore_errors=True)

def callback(return_data, tbar, metadata_writer):
    base_name, affinity = return_data
    
    row = [base_name, affinity]
    metadata_writer.writerow(row)
    
    tbar.update(1)

def process_pdb_list(pdb_list_file, output_file='check_vina_CPSea.csv', n_cores=1, log_file=None):
    """Process the PDB list file, perform analysis in parallel and write results to CSV"""
    # Determine number of processes
    if n_cores == -1:
        n_cores = max(1, cpu_count() - 1)
    
    # Read PDB file list
    with open(pdb_list_file, 'r') as f:
        pdb_files = [line.strip() for line in f if line.strip()]
    
    total_files = len(pdb_files)
    print(f"Processing {total_files} PDB files with {n_cores} processes...")
    
    # Open CSV file for output
    with open(output_file, 'w', newline='') as metadata_file:
        metadata_writer = csv.writer(metadata_file)
        metadata_writer.writerow(["input", "affinity"])
        tbar = tqdm(total=total_files)  # Fixed progress bar total
        partial_callback = partial(callback, tbar=tbar, metadata_writer=metadata_writer)
        pool = Pool(processes=n_cores)
        for pdb_file in pdb_files:
            pool.apply_async(func=process_pdb_file, args=(pdb_file,), callback=partial_callback)
        pool.close()
        pool.join()
        tbar.close()
    
    print(f"Processing completed. Results saved to {output_file}")

def main():
    """Main function, parse command line arguments and execute processing"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Process PDB files (compressed or uncompressed) with AutoDock Vina')
    parser.add_argument('-i', '--input', default='0427Benchmark50000_decompress.list', help='Input PDB list file')
    parser.add_argument('-o', '--output', default='check_vina_CPSea.csv', help='Output CSV file')
    parser.add_argument('-c', '--cores', type=int, default=1, help='Number of processes to use (-1 for all available minus 1)')
    parser.add_argument('--log', help='Log file path')
    
    args = parser.parse_args()
    
    # Record timing
    start_time = time.time()
    
    # Execute processing
    process_pdb_list(
        pdb_list_file=args.input,
        output_file=args.output,
        n_cores=args.cores,
        log_file=args.log
    )
    
    # Calculate and print total time
    end_time = time.time()
    print(f"Total processing time: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main()