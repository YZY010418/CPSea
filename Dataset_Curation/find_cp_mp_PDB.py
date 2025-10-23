#!/usr/bin/python
# -*- coding:utf-8 -*-
import warnings
import os
import multiprocessing as mp
from functools import partial
from tqdm import tqdm
from Bio.PDB.PDBIO import PDBIO
from find_cp_PDB import find_cyclic_peptide
from pathlib import Path
import gzip
from io import StringIO
import gc

warnings.filterwarnings("ignore")

def process_one(pdb_name, good_folder_path, average_folder_path, root):
    try:
        # Find the input names in root directory, and get the pdb string
        input_pdb_path = os.path.join(root, pdb_name)
        with open(input_pdb_path) as f:
            pdb_string_input = f.read()

        # Perform the cyclic peptide searching
        result = find_cyclic_peptide(pdb_string_input)

        # Deal with the outcome
        metadata = []
        for r in result:
            chain_ID, i, j, pep_seq, length, cb_distance, abs_bsa, rel_bsa, rel_nc_bsa, hydrophobic_ratio, peptide_helix_ratio, peptide_sheet_ratio, peptide_ss, output_structure, cyclization_struct, good_boy = r
            if output_structure is not None:
                io = PDBIO()
                string_io = StringIO()
                io.set_structure(output_structure)
                io.save(string_io)
                pdb_string = string_io.getvalue()
                pdb_bytes = pdb_string.encode('utf-8')
                if good_boy:
                    output_filename = f'{pdb_name[:-4]}_{chain_ID}_{i}_{j}.pdb.gz'
                    with gzip.open(os.path.join(good_folder_path, output_filename), 'wb') as f:
                        f.write(pdb_bytes)
                    metadata.append(f'{pdb_name[pdb_name.find("/")+1:-4]}_{chain_ID}_{i}_{j}.pdb\t{length}\t{cb_distance:.4f}\t{abs_bsa:.1f}\t{rel_bsa:.4f}\t{rel_nc_bsa:.4f}\t{hydrophobic_ratio:.4f}\t{peptide_helix_ratio:.4f}\t{peptide_sheet_ratio:.4f}\t{peptide_ss}\t{cyclization_struct}\tgood_boy\n')
                else:
                    output_filename = f'{pdb_name[:-4]}_{chain_ID}_{i}_{j}.pdb.gz'
                    with gzip.open(os.path.join(average_folder_path, output_filename), 'wb') as f:
                        f.write(pdb_bytes)
                    metadata.append(f'{pdb_name[pdb_name.find("/")+1:-4]}_{chain_ID}_{i}_{j}.pdb\t{length}\t{cb_distance:.4f}\t{abs_bsa:.1f}\t{rel_bsa:.4f}\t{rel_nc_bsa:.4f}\t{hydrophobic_ratio:.4f}\t{peptide_helix_ratio:.4f}\t{peptide_sheet_ratio:.4f}\t{peptide_ss}\t{cyclization_struct}\taverage\n')
                string_io.close()
                del io, string_io, pdb_bytes, output_structure
            del chain_ID, i, j, pep_seq, length, cb_distance, abs_bsa, rel_bsa, rel_nc_bsa, hydrophobic_ratio, peptide_helix_ratio, peptide_sheet_ratio, peptide_ss, cyclization_struct, good_boy
        del result, pdb_string_input
        gc.collect()
        return pdb_name, metadata
    except Exception as e:
        print(f"Error processing {pdb_name}: {e}")
        return pdb_name, [] 
    finally:
        gc.collect() 
    
def callback(return_data, metadata, finished, tbar):
    finished_pdb, individual_metadata = return_data
    tbar.update(1)
    if individual_metadata: 
        metadata.writelines(individual_metadata)
    finished.writelines(f"{finished_pdb}\n")
    metadata.flush()
    finished.flush()

def main(saved_dir, N_CPU, input_list, root, good_folder_path, average_folder_path):
    with open(input_list) as f:
        pdb_list = [line.strip() for line in f.readlines()]
    
    finished_path = os.path.join(saved_dir, 'finished_0725_PDB.csv')
    finished_set = set()
    if os.path.exists(finished_path):
        with open(finished_path) as fin:
            finished_set = set(line.strip() for line in fin)

    todo_list = [pdb for pdb in pdb_list if pdb not in finished_set]
    print(f"Total PDBs: {len(pdb_list)}, Already processed: {len(finished_set)}, To process: {len(todo_list)}")
    
    with open(os.path.join(saved_dir, 'metadata_0725_PDB.csv'), 'a') as metadata, \
         open(finished_path, 'a') as finished:
        
        tbar = tqdm(total=len(todo_list))
        partial_callback = partial(callback, metadata=metadata, finished=finished, tbar=tbar)
        pool = mp.Pool(N_CPU)
        
        for pdb in todo_list:
            pool.apply_async(
                func=process_one,
                args=(pdb, good_folder_path, average_folder_path, root),
                callback=partial_callback
            )
        
        pool.close()
        pool.join()
        tbar.close()

if __name__ == '__main__':
    N_CPU = 32

    input_list = '/data_hdd/home/yangziyi/Projects/CPSea2.0/FindCPBase/PDB70_plain.list'  # A list of PDB names
    root = '/data_hdd/home/yangziyi/PDB70/'  # Path to input PDBs
    saved_dir = '/data_hdd/home/yangziyi/Projects/CyclicPep/FindCycPep/0725_PDB/'  # Path to output unrelaxed peptides.
    good_folder_path = os.path.join(saved_dir, 'good_candidates')
    average_folder_path = os.path.join(saved_dir, 'average_candidates')
    Path(saved_dir).mkdir(parents=True, exist_ok=True)  # Create the directory for saving results.
    Path(good_folder_path).mkdir(parents=True, exist_ok=True)
    Path(average_folder_path).mkdir(parents=True, exist_ok=True)

    main(saved_dir, N_CPU, input_list, root, good_folder_path, average_folder_path)