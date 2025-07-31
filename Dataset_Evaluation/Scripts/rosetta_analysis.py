from multiprocessing import Pool
from functools import partial
import subprocess
from tqdm import tqdm
import sys
import os
import gzip
import shutil
import argparse

def process_one(pdb_file, output_csv):
    """
    Function to run the Rosetta command
    :param pdb_file: Input PDB file
    :param output_csv: Path to the output file of rosetta
    """
    temp_pdb_path = None
    if pdb_file.endswith('.pdb.gz'):
        # Create a temporary directory if it doesn't exist
        temp_dir = 'tmp_rosetta'
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        # Extract the .pdb.gz file to the temporary directory
        base_name = os.path.basename(pdb_file).replace('.gz', '')
        temp_pdb_path = os.path.join(temp_dir, base_name)
        with gzip.open(pdb_file, 'rb') as f_in:
            with open(temp_pdb_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        pdb_file = temp_pdb_path

    command = [
        "/data_hdd/home/yangziyi/Tools/rosetta.binary.ubuntu.release-371/main/source/bin/rosetta_scripts.default.linuxgccrelease",
        "-s", pdb_file,
        "-in:file:native", pdb_file,
        "-parser:protocol", "interface_analyze_ours.xml",
        "-out:file:scorefile", output_csv,
        "-out:file:scorefile_format", "json",
        "-out:file:score_only"
    ]
    try:
        with open('rosetta_stdout.log', 'a') as stdout_file, open('rosetta_stderr.log', 'a') as stderr_file:
            subprocess.run(command, check=True, shell=False, stdout=subprocess.DEVNULL, stderr=stderr_file)
        result = True
    except subprocess.CalledProcessError as e:
        #print(f"Error processing {pdb_file}: {e}")
        result = False

    if temp_pdb_path and os.path.exists(temp_pdb_path):
        os.remove(temp_pdb_path)

    return result

def callback(result, tbar):
    tbar.update(1)

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run Rosetta on a list of PDB files.')
    parser.add_argument('-i', '--input', required=True, help='Path to the input file list')
    parser.add_argument('-o', '--output', required=True, help='Path to the output file of rosetta')
    parser.add_argument('-c', '--cores', type=int, required=True, help='Number of CPU cores to use')
    args = parser.parse_args()

    list_file_path = args.input
    output_csv = args.output
    num_cores = args.cores

    try:
        # Read the list file, specifying the encoding format
        f = open(list_file_path, 'r', encoding='utf-8') 
        tbar = tqdm(total=2640000)
        partial_callback = partial(callback, tbar=tbar)
        pool = Pool(cores=num_cores)
        for pdb in f:
            pool.apply_async(func=process_one, args=(pdb.strip(), output_csv), callback=partial_callback)
        pool.close()
        pool.join()


    except FileNotFoundError:
        print(f"List file not found: {list_file_path}")

if __name__ == "__main__":
    main()