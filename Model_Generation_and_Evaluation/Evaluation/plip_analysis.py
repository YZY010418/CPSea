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
import argparse

warnings.filterwarnings("ignore")

FIXED_HEADERS = [
    'Input',  
    'hydrophobic_interactions', 
    'hydrogen_bonds', 
    'water_bridges', 
    'salt_bridges', 
    'pi_stacks', 
    'pi_cation_interactions', 
    'halogen_bonds', 
    'metal_complexes'
]

def count_interactions(xml_file):
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()

        interaction_counts = {}

        interactions = root.find('bindingsite/interactions')
        if interactions is not None:
            for child in interactions:
                interaction_type = child.tag
                interaction_count = len(child.findall('./*'))
                interaction_counts[interaction_type] = interaction_count

        return interaction_counts
    except Exception as e:
        print(f"Error processing {xml_file}: {e}", file=sys.stderr)
        return {}

def process_pdb_file(pdb_path):
    try:
        tmp_dir = os.path.abspath("./plip_tempfiles")
        os.makedirs(tmp_dir, exist_ok=True)
        
        file_name = os.path.basename(pdb_path)
        base_name = os.path.splitext(os.path.splitext(file_name)[0])[0]  
        
        temp_pdb_file = os.path.join(tmp_dir, f"{base_name}.pdb")
        
        # Check if file is gzipped by extension or magic number
        is_gzipped = pdb_path.endswith('.gz')
        if not is_gzipped:
            with open(pdb_path, 'rb') as f:
                try:
                    header = f.read(2)
                    is_gzipped = (header == b'\x1f\x8b')
                except:
                    pass
        
        # Handle both gzipped and uncompressed files
        if is_gzipped:
            with gzip.open(pdb_path, 'rb') as f_in:
                with open(temp_pdb_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            with open(pdb_path, 'rb') as f_in:
                with open(temp_pdb_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        
        output_xml = os.path.join(tmp_dir, base_name)
        command = f"plip -f {temp_pdb_file} --peptides L -x -o {output_xml} --silent"
        exit_code = os.system(command)
        
        if exit_code != 0:
            print(f"Warning: PLIP command failed for {pdb_path}", file=sys.stderr)
            result = {"file_name": file_name, "interactions": {}}
        else:
            xml_file = f"{output_xml}/report.xml"
            if not os.path.exists(xml_file):
                print(f"Error: XML file {xml_file} not found", file=sys.stderr)
                result = {"file_name": file_name, "interactions": {}}
            else:
                result = {"file_name": file_name, "interactions": count_interactions(xml_file)}
        
        # Cleanup temporary files
        if os.path.exists(temp_pdb_file):
            os.remove(temp_pdb_file)
        if os.path.exists(xml_file):
            os.remove(xml_file)
        
        temp_folder = os.path.join(tmp_dir, base_name)
        if os.path.exists(temp_folder):
            shutil.rmtree(temp_folder)
        
        return result
    except Exception as e:
        print(f"Error processing {pdb_path}: {e}", file=sys.stderr)
        try:
            temp_pdb_file = os.path.join("./tmp", f"{os.path.splitext(os.path.splitext(os.path.basename(pdb_path))[0])[0]}.pdb")
            if os.path.exists(temp_pdb_file):
                os.remove(temp_pdb_file)
            temp_folder = os.path.join("./tmp", os.path.splitext(os.path.splitext(os.path.basename(pdb_path))[0])[0])
            if os.path.exists(temp_folder):
                shutil.rmtree(temp_folder)
        except:
            pass
        return {"file_name": os.path.basename(pdb_path), "interactions": {}}

def callback(return_data, tbar, metadata_writer, log_file):
    file_name = return_data.get("file_name", "")
    interaction_counts = return_data.get("interactions", {})
    
    standardized_counts = standardize_interaction_types(interaction_counts)
    
    row = [file_name]  
    for header in FIXED_HEADERS[1:]:  
        row.append(standardized_counts.get(header, 0))
    
    metadata_writer.writerow(row)
    
    tbar.update(1)
    
    if log_file:
        log_file.write(f"Processed: {file_name} - {sum(standardized_counts.values())} interactions\n")

def standardize_interaction_types(interaction_counts):
    aliases = {
        'hydrophobic': 'hydrophobic_interactions',
        'hydrophobic_interaction': 'hydrophobic_interactions',
        'h_bond': 'hydrogen_bonds',
        'h_bonds': 'hydrogen_bonds',
        'waterbridge': 'water_bridges',
        'water_bridge': 'water_bridges',
        'saltbridge': 'salt_bridges',
        'salt_bridge': 'salt_bridges',
        'pistack': 'pi_stacks',
        'pi_stack': 'pi_stacks',
        'pication': 'pi_cation_interactions',
        'pi_cation': 'pi_cation_interactions',
        'halogenbond': 'halogen_bonds',
        'halogen_bond': 'halogen_bonds',
        'metalcomplex': 'metal_complexes',
        'metal_complex': 'metal_complexes'
    }
    
    standardized = {}
    for key, value in interaction_counts.items():
        normalized_key = key.lower().replace('_', '').rstrip('s')
        standardized_key = aliases.get(normalized_key, key)
        standardized[standardized_key] = value
    
    return standardized

def process_pdb_list(pdb_list_file, output_file='PLIP_metadata_CPSet.csv', n_jobs=-1, log_file=None):
    if n_jobs == -1:
        n_jobs = max(1, cpu_count() - 1)
    
    with open(pdb_list_file, 'r') as f:
        pdb_files = [line.strip() for line in f if line.strip()]
    
    total_files = len(pdb_files)
    print(f"Processing {total_files} PDB files with {n_jobs} processes...")
    
    if log_file:
        log = open(log_file, 'w')
    else:
        log = None
    
    with open(output_file, 'w', newline='') as metadata_file:
        metadata_writer = csv.writer(metadata_file)
        
        metadata_writer.writerow(FIXED_HEADERS)

        with tqdm(total=total_files, desc="Processing PDBs") as tbar:
            partial_callback = partial(
                callback, 
                tbar=tbar, 
                metadata_writer=metadata_writer,
                log_file=log
            )
            
            with Pool(processes=n_jobs) as pool:
                results = []
                for pdb_file in pdb_files:
                    if os.path.exists(pdb_file):
                        result = pool.apply_async(
                            func=process_pdb_file, 
                            args=(pdb_file,), 
                            callback=partial_callback
                        )
                        results.append(result)
                    else:
                        print(f"Warning: File not found: {pdb_file}", file=sys.stderr)
                        tbar.update(1)
                
                for result in results:
                    result.wait()

    if log:
        log.close()
    
    print(f"Processing completed. Results saved to {output_file}")

def main():
    
    parser = argparse.ArgumentParser(description='Process PDB files with PLIP and generate interaction statistics')
    parser.add_argument('-i', '--input', default='PDB_Path.list', help='Input PDB list file')
    parser.add_argument('-o', '--output', default='PLIP_metadata_CPSet.csv', help='Output CSV file')
    parser.add_argument('-c', '--cores', type=int, default=-1, help='Number of processes to use')
    parser.add_argument('--log', help='Log file path')
    
    args = parser.parse_args()
    
    start_time = time.time()
    
    process_pdb_list(
        pdb_list_file=args.input,
        output_file=args.output,
        n_jobs=args.cores,
        log_file=args.log
    )
    
    end_time = time.time()
    print(f"Total processing time: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main()