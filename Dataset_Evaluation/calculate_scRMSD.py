import os
import csv
import shutil
import argparse
from Bio.PDB import PDBParser, Superimposer

def get_identifier(basename):
    marker = "_relaxed_rank_001_alphafold"
    if marker in basename:
        return basename.split(marker)[0]
    else:
        return basename

def get_chain_info(structure):
    model = next(structure.get_models())
    chain_ids = list(model.child_dict.keys())
    return chain_ids, len(chain_ids)

def get_ca_atoms(structure):
    chain_ids, chain_count = get_chain_info(structure)
    ca_atoms = []
    
    model = next(structure.get_models())
    if chain_count > 1:
        if 'L' in chain_ids:
            chain = model['L']
            for residue in chain:
                if 'CA' in residue:
                    ca_atoms.append(residue['CA'])
    else:
        if chain_ids:
            chain = model[chain_ids[0]]
            for residue in chain:
                if 'CA' in residue:
                    ca_atoms.append(residue['CA'])
    
    return ca_atoms, chain_count, chain_ids

def calculate_plddt_average(ca_atoms):
    if not ca_atoms:
        return None
    b_factors = [atom.get_bfactor() for atom in ca_atoms]
    return sum(b_factors) / len(b_factors)

def calculate_superimposed_rmsd(structure1, structure2):
    ca1, _, _ = get_ca_atoms(structure1)
    ca2, _, _ = get_ca_atoms(structure2)
    
    if len(ca1) != len(ca2):
        raise ValueError(f"CA atom count mismatch: {len(ca1)} vs {len(ca2)}")
    if len(ca1) == 0:
        raise ValueError("No valid CA atoms found")
    
    sup = Superimposer()
    sup.set_atoms(ca1, ca2)
    return sup.rms

def find_a3m_files(a3m_folder):
    a3m_files = {}
    for root, _, files in os.walk(a3m_folder):
        for file in files:
            if file.lower().endswith('.a3m'):
                rel_path = os.path.relpath(root, a3m_folder)
                basename = os.path.splitext(file)[0]
                identifier = get_identifier(basename)
                lower_id = identifier.lower()
                
                if lower_id not in a3m_files:
                    a3m_files[lower_id] = []
                a3m_files[lower_id].append((os.path.join(root, file), rel_path))
    
    return a3m_files

def copy_a3m_files(unique_identifiers, a3m_folder, a3m_files):
    unfinished_folder = os.path.join(a3m_folder, 'unfinished')
    os.makedirs(unfinished_folder, exist_ok=True)
    copied_count = 0
    
    for identifier in unique_identifiers:
        lower_id = identifier.lower()
        if lower_id in a3m_files:
            for file_path, rel_path in a3m_files[lower_id]:
                target_dir = os.path.join(unfinished_folder, rel_path)
                os.makedirs(target_dir, exist_ok=True)
                
                target_path = os.path.join(target_dir, os.path.basename(file_path))
                try:
                    shutil.copy2(file_path, target_path)
                    copied_count += 1
                    print(f"Copied: {os.path.basename(file_path)} to {target_dir}")
                except Exception as e:
                    print(f"Failed to copy {os.path.basename(file_path)}: {str(e)}")
    
    print(f"Total {copied_count} a3m files copied to {unfinished_folder}")
    return copied_count

def process_pdb_files(folder1, folder2, a3m_folder, output_csv):
    files1 = [f for f in os.listdir(folder1) if f.endswith('.pdb')]
    files2 = [f for f in os.listdir(folder2) if f.endswith('.pdb')]
    
    pdb_ids1 = {}
    for f in files1:
        basename = os.path.splitext(f)[0]
        id_prefix = get_identifier(basename)
        lower_id = id_prefix.lower()
        pdb_ids1[lower_id] = (f, f)
    
    pdb_ids2 = {}
    for f in files2:
        basename = os.path.splitext(f)[0]
        id_prefix = get_identifier(basename)
        lower_id = id_prefix.lower()
        pdb_ids2[lower_id] = (f, f)
    
    unique_to_folder2 = set(pdb_ids2.keys()) - set(pdb_ids1.keys())
    print(f"Found {len(unique_to_folder2)} identifiers present only in folder2")
    
    if a3m_folder and os.path.isdir(a3m_folder):
        print(f"Searching for files in a3m folder: {a3m_folder}")
        a3m_files = find_a3m_files(a3m_folder)
        copy_a3m_files(unique_to_folder2, a3m_folder, a3m_files)
    elif a3m_folder:
        print(f"Warning: a3m folder does not exist - {a3m_folder}")
    
    common_lower_ids = set(pdb_ids1.keys()) & set(pdb_ids2.keys())
    print(f"Found {len(common_lower_ids)} matching PDB file pairs (new identification rule, case-insensitive)")
    
    parser = PDBParser(QUIET=True)
    results = [
        ["Identifier", "Folder1 Filename", "Folder2 Filename", "RMSD", "CA_Count", 
         "Status", "Folder1 pLDDT Average", "Folder2 pLDDT Average"]
    ]
    
    for lower_id in sorted(common_lower_ids):
        orig_name1, file1 = pdb_ids1[lower_id]
        orig_name2, file2 = pdb_ids2[lower_id]
        
        file1_path = os.path.join(folder1, file1)
        file2_path = os.path.join(folder2, file2)
        
        try:
            structure1 = parser.get_structure(orig_name1, file1_path)
            structure2 = parser.get_structure(orig_name2, file2_path)
            
            ca1, chain_count1, chain_ids1 = get_ca_atoms(structure1)
            ca2, chain_count2, chain_ids2 = get_ca_atoms(structure2)
            
            plddt1 = None
            if orig_name1.startswith('AF') or orig_name1.startswith('af'):
                plddt1 = calculate_plddt_average(ca1)
            
            plddt2 = None
            if orig_name2.startswith('AF') or orig_name2.startswith('af'):
                plddt2 = calculate_plddt_average(ca2)
            
            if len(ca1) == 0 or len(ca2) == 0:
                msg1 = ""
                if len(ca1) == 0:
                    if chain_count1 > 1 and 'L' not in chain_ids1:
                        msg1 = "Folder1: Multi-chain structure without L chain"
                    elif chain_count1 == 1:
                        msg1 = "Folder1: Single-chain structure without CA atoms"
                    else:
                        msg1 = "Folder1: No valid chains or CA atoms"
                
                msg2 = ""
                if len(ca2) == 0:
                    if chain_count2 > 1 and 'L' not in chain_ids2:
                        msg2 = "Folder2: Multi-chain structure without L chain"
                    elif chain_count2 == 1:
                        msg2 = "Folder2: Single-chain structure without CA atoms"
                    else:
                        msg2 = "Folder2: No valid chains or CA atoms"
                
                status = "; ".join(filter(None, [msg1, msg2]))
                plddt1_str = f"{plddt1:.2f}" if plddt1 is not None else ""
                plddt2_str = f"{plddt2:.2f}" if plddt2 is not None else ""
                results.append([
                    lower_id.upper(), orig_name1, orig_name2, "", "", status,
                    plddt1_str, plddt2_str
                ])
                print(f"{lower_id}: {status}")
                continue
            
            rmsd = calculate_superimposed_rmsd(structure1, structure2)
            
            plddt1_str = f"{plddt1:.2f}" if plddt1 is not None else ""
            plddt2_str = f"{plddt2:.2f}" if plddt2 is not None else ""
            results.append([
                lower_id.upper(), orig_name1, orig_name2, f"{rmsd:.4f}", len(ca1),
                "Success", plddt1_str, plddt2_str
            ])
            print(f"{lower_id}: RMSD = {rmsd:.4f} Å, CA atom count = {len(ca1)}")
            
        except Exception as e:
            error_msg = str(e).split('\n')[0]
            results.append([
                lower_id.upper(), orig_name1, orig_name2, "", "", f"Error: {error_msg}",
                "", ""
            ])
            print(f"{lower_id}: Processing error - {error_msg}")
    
    only_in1 = set(pdb_ids1.keys()) - set(pdb_ids2.keys())
    for lower_id in sorted(only_in1):
        orig_name, file1 = pdb_ids1[lower_id]
        file1_path = os.path.join(folder1, file1)
        plddt1_str = ""
        if orig_name.startswith('AF'):
            try:
                structure = parser.get_structure(orig_name, file1_path)
                ca, _, _ = get_ca_atoms(structure)
                plddt = calculate_plddt_average(ca)
                plddt1_str = f"{plddt:.2f}" if plddt is not None else ""
            except:
                plddt1_str = "Calculation failed"
        results.append([
            lower_id.upper(), orig_name, "", "", "", "Present only in Folder1",
            plddt1_str, ""
        ])
        print(f"{lower_id}: Present only in Folder1")
    
    only_in2 = unique_to_folder2
    for lower_id in sorted(only_in2):
        orig_name, file2 = pdb_ids2[lower_id]
        file2_path = os.path.join(folder2, file2)
        plddt2_str = ""
        if orig_name.startswith('AF'):
            try:
                structure = parser.get_structure(orig_name, file2_path)
                ca, _, _ = get_ca_atoms(structure)
                plddt = calculate_plddt_average(ca)
                plddt2_str = f"{plddt:.2f}" if plddt is not None else ""
            except:
                plddt2_str = "Calculation failed"
        results.append([
            lower_id.upper(), "", orig_name, "", "", "Present only in Folder2",
            "", plddt2_str
        ])
        print(f"{lower_id}: Present only in Folder2")
    
    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(results)
    
    print(f"Results saved to {output_csv}")

def main():
    parser = argparse.ArgumentParser(description='Calculate CA atom RMSD values and pLDDT averages for PDB files in two folders, and process related a3m files')
    parser.add_argument('folder1', help='Path to the first folder containing PDB files')
    parser.add_argument('folder2', help='Path to the second folder containing PDB files')
    parser.add_argument('a3m_folder', help='Path to the folder containing a3m files')
    parser.add_argument('-o', '--output', default='ca_rmsd_results.csv', 
                      help='Path for the output CSV file (default: ca_rmsd_results.csv)')
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.folder1):
        print(f"Error: Folder1 does not exist - {args.folder1}")
        return
    if not os.path.isdir(args.folder2):
        print(f"Error: Folder2 does not exist - {args.folder2}")
        return
    if not os.path.isdir(args.a3m_folder):
        print(f"Error: a3m folder does not exist - {args.a3m_folder}")
        return
    
    process_pdb_files(args.folder1, args.folder2, args.a3m_folder, args.output)
    print("Processing completed")

if __name__ == "__main__":
    main()