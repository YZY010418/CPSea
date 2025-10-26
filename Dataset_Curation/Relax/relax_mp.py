#!/usr/bin/python
# -*- coding:utf-8 -*-
import os
import csv
from openmm.app import Modeller
from openmm import app as openmm_app
from openmm import unit
from pathlib import Path
from tqdm import tqdm
from functools import partial
import multiprocessing as mp
ENERGY = unit.kilocalories_per_mole
LENGTH = unit.angstroms

from k_to_de_0424 import ForceFieldMinimizerKtoDE
from cys_to_cys_0424 import ForceFieldMinimizerCys
from head_tail_0424 import ForceFieldMinimizerHeadTail

def cyclization_and_minimize(relax_input, folder_paths, root): # input: (pdb_name, cyclization_type, length)
	
	PDB_name = relax_input[0]
	output_name = PDB_name[:-4] + '_relaxed.pdb.gz'
	if relax_input[1] == 'DISULFIDE':
		force_field = ForceFieldMinimizerCys()
		energy_info = force_field(os.path.join(root,PDB_name), os.path.join(folder_paths[0],output_name), cyclic_chains=['L'], cyclic_opts=[(('L', 1), ('L', int(relax_input[2])))]) 
	elif relax_input[1] == 'ISOPEPTIDE': 
		force_field = ForceFieldMinimizerKtoDE()
		energy_info =force_field(os.path.join(root,PDB_name), os.path.join(folder_paths[1],output_name), cyclic_chains=['L'], cyclic_opts=[(('L', 1), ('L', int(relax_input[2])))]) 
	elif relax_input[1] == 'HEADTAIL': 
		force_field = ForceFieldMinimizerHeadTail()
		energy_info =force_field(os.path.join(root,PDB_name), os.path.join(folder_paths[2],output_name), cyclic_chains=['L'], cyclic_opts=[(('L', 0), ('L', int(relax_input[2])-1))])

	return output_name, energy_info, relax_input

def callback(return_data, tbar, energy_data, abnormal_data):
	output_name, energy_info, relax_input = return_data
	if energy_info is not None:
		if energy_info[1] <= 5000:
			energy_data.write(f'{relax_input[0]}\t{output_name}\t{relax_input[1]}\t{relax_input[2]}\t')
			for i in range(len(energy_info)):
				energy_data.write(f'{energy_info[i]}\t')
			energy_data.write(f'\n')
		else:
			abnormal_data.write(f'{relax_input[0]}\t{output_name}\t{relax_input[1]}\t{relax_input[2]}\t')
			for i in range(len(energy_info)):
				abnormal_data.write(f'{energy_info[i]}\t')
			abnormal_data.write(f'\n') 
	else:
		abnormal_data.write(f'{relax_input[0]}\t{output_name}\t{relax_input[1]}\t{relax_input[2]}\tnone\tnone\n') # Input name, output name, cyc_type, length
	energy_data.flush()
	abnormal_data.flush()
	tbar.update(1)

if __name__ == '__main__':

        N_CPU = 48
        
        # Path to directory where relaxed structs are saved, and metadata for relax input
        save_path = '/data_hdd/home/yangziyi/Projects/CyclicPep/Relax/1015_PDB_Remake'
        root = '/data_hdd/home/yangziyi/Projects/CyclicPep/FindCycPep/0725_PDB/good_candidates/'
        input_metadata = '/data_hdd/home/yangziyi/Projects/CyclicPep/FindCycPep/0725_PDB/metadata_0725_PDB_good.tsv'

        # Make directories for relaxed structs to be saved
        cys_folder_path = os.path.join(save_path,'CysCysRelaxed')
        headtail_folder_path = os.path.join(save_path,'HeadTailRelaxed')
        isopep_folder_path = os.path.join(save_path,'IsoPepRelaxed')
        folder_paths = (cys_folder_path, isopep_folder_path, headtail_folder_path)

        Path(cys_folder_path).mkdir(parents=True, exist_ok=True)
        Path(headtail_folder_path).mkdir(parents=True, exist_ok=True)
        Path(isopep_folder_path).mkdir(parents=True, exist_ok=True)

        metadata = open(input_metadata, 'r')
        energy_data = open(os.path.join(save_path,'energy_data.csv'),'a')
        extremely_high_energy = open(os.path.join(save_path,'extremely_high_energy.csv'),'a')

        tbar = tqdm(total=16000)
        partial_callback = partial(callback, tbar=tbar, energy_data=energy_data, abnormal_data=extremely_high_energy)

        pool = mp.Pool(N_CPU)
        reader = csv.reader(metadata, delimiter = '\t')
        for row in reader:
            one_relax_input = (row[0], row[10], row[1])
            #print(one_relax_input)
            pool.apply_async(func=cyclization_and_minimize, args=(one_relax_input, folder_paths, root), callback=partial_callback)
        pool.close()
        pool.join()

        metadata.close()
        energy_data.close()
        extremely_high_energy.close()


