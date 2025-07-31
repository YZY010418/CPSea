#!/usr/bin/python
# -*- coding:utf-8 -*-
import os
import csv
import time
import datetime
import io
import logging
import numpy as np
import pdbfixer
import openmm
import random
from openmm.app import Modeller
from openmm import app as openmm_app
from openmm import unit
from pathlib import Path
from tqdm import tqdm
from functools import partial
import multiprocessing as mp

ENERGY = unit.kilocalories_per_mole
LENGTH = unit.angstroms

from k_to_de_model import ForceFieldMinimizerKtoDE
from cys_to_cys_model import ForceFieldMinimizerCys
from head_tail_model import ForceFieldMinimizerHeadTail

def determine_cyclization_type(distance):
    if 3 <= distance < 4.5:
        return 'DISULFIDE'
    elif 4.5 <= distance < 6:
        return 'HEADTAIL'
    elif 6 <= distance <= 8:
        return 'ISOPEPTIDE'
    else:
        return 'UNSUPPORTED'

def process_structure(input_data, save_path):
    pdb_path, distance_str, length = input_data
    pdb_filename = os.path.basename(pdb_path)
    output_name = pdb_filename[:-4] + '_relaxed.pdb'
    cb_distance = float(distance_str)
    cyclization_type = determine_cyclization_type(cb_distance)
    if cyclization_type == 'UNSUPPORTED':
        print(f"Unsupported distance {cb_distance} for PDB: {pdb_filename}")
        return output_name, None, (pdb_filename, distance_str, length, cyclization_type)

    if cyclization_type == 'DISULFIDE':
        force_field = ForceFieldMinimizerCys()
        energy_info = force_field(pdb_path, os.path.join(save_path, output_name), 
                              cyclic_chains=['L'], 
                              cyclic_opts=[(('L', 0), ('L', int(length)-1))]
        )
    elif cyclization_type == 'ISOPEPTIDE':
        force_field = ForceFieldMinimizerKtoDE()
        energy_info = force_field(pdb_path, os.path.join(save_path, output_name), 
                              cyclic_chains=['L'], 
                              cyclic_opts=[(('L', 0), ('L', int(length)-1))]
        )
    elif cyclization_type == 'HEADTAIL':
        force_field = ForceFieldMinimizerHeadTail()
        energy_info = force_field(pdb_path, os.path.join(save_path, output_name), 
                              cyclic_chains=['L'], 
                              cyclic_opts=[(('L', 0), ('L', int(length)-1))]
        )

    return output_name, energy_info, (pdb_filename, distance_str, length, cyclization_type)

def write_result(return_data, progress_bar, result_file):
    output_name, energy_info, input_data = return_data
    pdb_name, distance, length, cyclization_type = input_data
    
    if energy_info is not None:
        result_file.write(f'{pdb_name}\t{output_name}\t{cyclization_type}\t{length}\t')
        for value in energy_info:
            result_file.write(f'{value}\t')
        result_file.write('\n')
    else:
        result_file.write(f'{pdb_name}\t{output_name}\t{cyclization_type}\t{length}\tnone\tnone\n')
    
    result_file.flush()
    progress_bar.update(1)

def main():
    config = {
        'num_processes': 64,
        'save_path': '/data_hdd/home/yangziyi/Projects/CyclicPep/evaluate_models/CPSet_Glad/glad_relaxed',
        'metadata_file': '/data_hdd/home/yangziyi/Projects/CyclicPep/evaluate_models/CPSet_Glad/CPSet_Glad_CB_length.csv',
        'output_file': 'energy_data.csv',
        'total_tasks': 1800
    }

    Path(config['save_path']).mkdir(parents=True, exist_ok=True)
    
    with open(config['metadata_file'], 'r') as metadata, \
         open(os.path.join(config['save_path'], config['output_file']), 'a') as result_file:
        
        progress_bar = tqdm(total=config['total_tasks'])
        callback = partial(write_result, progress_bar=progress_bar, result_file=result_file)
        
        with mp.Pool(config['num_processes']) as pool:
            reader = csv.reader(metadata)
            for row in reader:
                input_data = (row[0], row[1], row[2])
                pool.apply_async(func=process_structure, 
                                args=(input_data, config['save_path']), 
                                callback=callback)
            
            pool.close()
            pool.join()

if __name__ == "__main__":
    main()