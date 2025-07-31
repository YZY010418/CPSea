from __future__ import division, print_function
import math
import os
import sys
import argparse
import numpy as np
from Bio import PDB
from multiprocessing import Pool
from tqdm import tqdm
import gzip
import shutil
from functools import partial
import csv

RAMA_PREF_VALUES = None

RAMA_PREFERENCES = {
    "General": {
        "file": os.path.join('data', 'pref_general.data'),
        "bounds": [0, 0.0005, 0.02, 1],
    },
    "GLY": {
        "file": os.path.join('data', 'pref_glycine.data'),
        "bounds": [0, 0.002, 0.02, 1],
    },
    "PRO": {
        "file": os.path.join('data', 'pref_proline.data'),
        "bounds": [0, 0.002, 0.02, 1],
    },
    "PRE-PRO": {
        "file": os.path.join('data', 'pref_preproline.data'),
        "bounds": [0, 0.002, 0.02, 1],
    }
}


def _cache_RAMA_PREF_VALUES():
    global RAMA_PREF_VALUES 
    f_path = "/data_hdd/home/yangziyi/Tools/PyRAMA/pyrama"
    RAMA_PREF_VALUES = {}
    for key, val in RAMA_PREFERENCES.items():
        RAMA_PREF_VALUES[key] = np.full((360, 360), 0, dtype=np.float64)
        with open(os.path.join(f_path, val["file"])) as fn:
            for line in fn:
                if line.startswith("#"):
                    continue
                else:
                    x = int(float(line.split()[1]))
                    y = int(float(line.split()[0]))
                    RAMA_PREF_VALUES[key][x + 180][y + 180] \
                        = RAMA_PREF_VALUES[key][x + 179][y + 179] \
                        = RAMA_PREF_VALUES[key][x + 179][y + 180] \
                        = RAMA_PREF_VALUES[key][x + 180][y + 179] \
                        = float(line.split()[2])
    return RAMA_PREF_VALUES


def process_one(file_name):

    global RAMA_PREF_VALUES

    if RAMA_PREF_VALUES is None:
        RAMA_PREF_VALUES = _cache_RAMA_PREF_VALUES()
        
    tmp_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tmp_rama')
    if not os.path.exists(tmp_folder):
        os.makedirs(tmp_folder)


    if file_name.endswith('.pdb.gz'):
        tmp_pdb_file = os.path.join(tmp_folder, os.path.basename(file_name).replace('.pdb.gz', '.pdb'))
        with gzip.open(file_name, 'rb') as f_in:
            with open(tmp_pdb_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        file_to_process = tmp_pdb_file
    else:
        file_to_process = file_name

    accept_count = 0
    favoured_count = 0
    total_count = 0
    accept_count_nc = 0
    favoured_count_nc = 0
    total_count_nc = 0

    structure = PDB.PDBParser().get_structure('input_structure', file_to_process)
    for model in structure:
        for chain in model:
            if chain.id == 'L':
                polypeptides = PDB.PPBuilder().build_peptides(chain)
                for poly_index, poly in enumerate(polypeptides):
                    phi_psi = poly.get_phi_psi_list()
                    poly_len = len(poly)
                    nc_indices_to_check = [1, 2, poly_len - 1, poly_len - 2]
                    indices_to_check = range(len(poly))

                    for res_index in indices_to_check:
                        residue = poly[res_index]
                        res_name = "{}".format(residue.resname)
                        phi, psi = phi_psi[res_index]
                        if phi and psi:
                            if str(poly[res_index + 1].resname) == "PRO" if res_index < poly_len - 1 else False:
                                aa_type = "PRE-PRO"
                            elif res_name == "PRO":
                                aa_type = "PRO"
                            elif res_name == "GLY":
                                aa_type = "GLY"
                            else:
                                aa_type = "General"
                            phi_deg = math.degrees(phi)
                            psi_deg = math.degrees(psi)
                            acceptable = RAMA_PREF_VALUES[aa_type][int(psi_deg) + 180][int(phi_deg) + 180] >= \
                                         RAMA_PREFERENCES[aa_type]["bounds"][1]
                            favoured = RAMA_PREF_VALUES[aa_type][int(psi_deg) + 180][int(phi_deg) + 180] >= \
                                         RAMA_PREFERENCES[aa_type]["bounds"][2]
                            total_count += 1
                            if acceptable:
                                accept_count += 1
                            if favoured:
                                favoured_count += 1

                    for res_index in indices_to_check:
                        if res_index < 0 or res_index >= poly_len:
                            continue
                        residue = poly[res_index]
                        res_name = "{}".format(residue.resname)
                        phi, psi = phi_psi[res_index]
                        if phi and psi:
                            if str(poly[res_index + 1].resname) == "PRO" if res_index < poly_len - 1 else False:
                                aa_type = "PRE-PRO"
                            elif res_name == "PRO":
                                aa_type = "PRO"
                            elif res_name == "GLY":
                                aa_type = "GLY"
                            else:
                                aa_type = "General"
                            phi_deg = math.degrees(phi)
                            psi_deg = math.degrees(psi)
                            acceptable = RAMA_PREF_VALUES[aa_type][int(psi_deg) + 180][int(phi_deg) + 180] >= \
                                         RAMA_PREFERENCES[aa_type]["bounds"][1]
                            favoured = RAMA_PREF_VALUES[aa_type][int(psi_deg) + 180][int(phi_deg) + 180] >= \
                                         RAMA_PREFERENCES[aa_type]["bounds"][2]
                            total_count_nc += 1
                            if acceptable:
                                accept_count_nc += 1
                            if favoured:
                                favoured_count_nc += 1

    if file_name.endswith('.pdb.gz'):
        os.remove(tmp_pdb_file)

    return file_name, total_count, accept_count, favoured_count, total_count_nc, accept_count_nc, favoured_count_nc


def callback(result, tbar, csvfile):
    file_name, total_count, accept_count, favoured_count, total_count_nc, accept_count_nc, favoured_count_nc = result
    accept_rate = accept_count / total_count if total_count > 0 else 0
    accept_rate_nc = accept_count_nc / total_count_nc if total_count_nc > 0 else 0
    favoured_rate = favoured_count / total_count if total_count > 0 else 0
    favoured_rate_nc = favoured_count_nc / total_count_nc if total_count_nc > 0 else 0

    writer = csv.DictWriter(csvfile, fieldnames=['input', 'total_angle_numbers', 'total_accept', 'accept_rate', 'total_favoured', 'favoured_rate', 'total_terminal_angle_numbers',
                                                'total_accept_nc', 'accept_rate_nc', 'total_favoured_nc', 'favoured_rate_nc'])
    row = {
        'input': file_name,
        'total_angle_numbers': total_count,
        'total_accept': accept_count,
        'accept_rate': accept_rate,
        'total_favoured': favoured_count,
        'favoured_rate': favoured_rate,
        'total_terminal_angle_numbers': total_count_nc,
        'total_accept_nc': accept_count_nc,
        'accept_rate_nc': accept_rate_nc,
        'total_favoured_nc': favoured_count_nc,
        'favoured_rate_nc': favoured_rate_nc
    }
    writer.writerow(row)

    tbar.update(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process PDB files with multi - processing.')
    parser.add_argument('-i', '--input', type=str, help='Path to the file containing PDB file paths.')
    parser.add_argument('-c', '--cores', type=int, default=1, help='Number of CPU cores to use.')
    parser.add_argument('-o', '--output', type=str, help='Path to the output file.')
    parser.add_argument('--terminal', action='store_true', help='Only calculate angles for the terminal amino acids.')
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"Error: The specified list file {args.input} does not exist.")
        sys.exit(1)

    csvfile = open(args.output, 'w', newline='')
    fieldnames = ['input', 'total_angle_numbers', 'total_accept', 'accept_rate', 'total_favoured', 'favoured_rate', 'total_terminal_angle_numbers',
                  'total_accept_nc', 'accept_rate_nc', 'total_favoured_nc', 'favoured_rate_nc']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    f = open(args.input, 'r')
    tbar = tqdm(total=2640000) # For really large datasets
    partial_callback = partial(callback, tbar=tbar, csvfile=csvfile)
    pool = Pool(processes=args.cores)
    for pdb in f:
        pool.apply_async(func=process_one, args=(pdb.strip(),), callback=partial_callback)
    pool.close()
    pool.join()

    csvfile.close()