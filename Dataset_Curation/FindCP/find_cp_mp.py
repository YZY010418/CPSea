#!/usr/bin/python
# -*- coding:utf-8 -*-
import warnings
import os
import multiprocessing as mp
from functools import partial
from tqdm import tqdm
from Bio.PDB.PDBIO import PDBIO
from find_cp import find_cyclic_peptide
from pathlib import Path
import gzip
from io import StringIO
import gc

warnings.filterwarnings("ignore")

def process_one(pdb_name,good_folder_path,average_folder_path,root):

    try:
        # Find the input names in root directory, and get the pdb string
        input_pdb_path = os.path.join(root,pdb_name)
        with open(input_pdb_path) as f:
            pdb_string_input = f.read()

        # Perform the cyclic peptide searching
        result = find_cyclic_peptide(pdb_string_input)

        # Deal with the outcome
        metadata = []
        for r in result:
            i,j,pep_seq,length,cb_distance,worst_pLDDT,abs_bsa,rel_bsa,rel_nc_bsa,hydrophobic_ratio,peptide_helix_ratio,peptide_sheet_ratio,peptide_ss,output_structure,cyclization_struct,good_boy = r
            if output_structure is not None:
                io = PDBIO()
                string_io = StringIO()
                io.set_structure(output_structure)
                io.save(string_io)
                pdb_string = string_io.getvalue()
                pdb_bytes = pdb_string.encode('utf-8')
                if good_boy == True:
                    with gzip.open(os.path.join(good_folder_path, f'{pdb_name[:-4]}_{i}_{j}.pdb.gz'), 'wb') as f:
                        f.write(pdb_bytes)
                    metadata.append(f'{pdb_name[pdb_name.find("/")+1:-4]}_{i}_{j}.pdb\t{length}\t{cb_distance:.4f}\t{worst_pLDDT:.1f}\t{abs_bsa:.1f}\t{rel_bsa:.4f}\t{rel_nc_bsa:.4f}\t{hydrophobic_ratio:.4f}\t{peptide_helix_ratio:.4f}\t{peptide_sheet_ratio:.4f}\t{peptide_ss}\t{cyclization_struct}\tgood_boy\n') 
                else:
                    with gzip.open(os.path.join(average_folder_path, f'{pdb_name[:-4]}_{i}_{j}.pdb.gz'), 'wb') as f:
                        f.write(pdb_bytes)
                    metadata.append(f'{pdb_name[pdb_name.find("/")+1:-4]}_{i}_{j}.pdb\t{length}\t{cb_distance:.4f}\t{worst_pLDDT:.1f}\t{abs_bsa:.1f}\t{rel_bsa:.4f}\t{rel_nc_bsa:.4f}\t{hydrophobic_ratio:.4f}\t{peptide_helix_ratio:.4f}\t{peptide_sheet_ratio:.4f}\t{peptide_ss}\t{cyclization_struct}\taverage\n')
                string_io.close()
                del io, string_io, pdb_bytes, output_structure
            del i,j,pep_seq,length,cb_distance,worst_pLDDT,abs_bsa,rel_bsa,rel_nc_bsa,hydrophobic_ratio,peptide_helix_ratio,peptide_sheet_ratio,peptide_ss,cyclization_struct,good_boy
        del result, pdb_string_input
        gc.collect()
        return pdb_name,metadata
    except Exception as e:
        print(f"{e}")
    
def callback(return_data, metadata, finished, tbar):
    finished_pdb, individual_metadata = return_data
    tbar.update(1)
    metadata.writelines(individual_metadata)
    finished.writelines(f"{finished_pdb}\n")
    metadata.flush()
    finished.flush()

def main(saved_dir, N_CPU, input_list, root, good_folder_path, average_folder_path):
    
    f = open(input_list)
    metadata = open(os.path.join(saved_dir,'metadata_0424_001.tsv'),'a')
    finished = open(os.path.join(saved_dir,'finished_0424_001.txt'),'a')

    tbar = tqdm(total=5988928)
    partial_callback = partial(callback, metadata=metadata, finished=finished, tbar=tbar)
    pool = mp.Pool(N_CPU)
    for pdb in f:
        pool.apply_async(func=process_one, args=(pdb.strip(),good_folder_path,average_folder_path,root), callback=partial_callback)
    pool.close()
    pool.join()

    f.close()
    metadata.close()
    finished.close()

if __name__ == '__main__':
    N_CPU = 32

    input_list = '/data_hdd/home/yangziyi/Projects/CyclicPep/FindCycPep/0424/0424_001_input.list' # A list of PDB names
    root = '/data_hdd/protein/AlphaFoldDB/AF2_domains/' # Path to input PDBs
    saved_dir = '/data_hdd/home/yangziyi/Projects/CyclicPep/FindCycPep/0424/0901_test' # Path to output unrelaxed peptides.
    good_folder_path = os.path.join(saved_dir,'good_candidates')
    average_folder_path = os.path.join(saved_dir,'average_candidates')
    Path(saved_dir).mkdir(parents=True, exist_ok=True) # Create the directory for saving results.
    Path(good_folder_path).mkdir(parents=True, exist_ok=True)
    Path(average_folder_path).mkdir(parents=True, exist_ok=True)

    main(saved_dir, N_CPU, input_list, root, good_folder_path, average_folder_path)
