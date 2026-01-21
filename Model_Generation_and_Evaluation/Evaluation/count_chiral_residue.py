#!/usr/bin/env python3

import os
import argparse
import numpy as np
import biotite.structure as struc
from biotite.structure.io.pdb import PDBFile
from tqdm import tqdm
from multiprocessing import Pool

def analyze_pdb(args):
    pdb_path, target_chain = args
    try:
        pdb = PDBFile.read(pdb_path)
        array = pdb.get_structure()[0]
        
        mask = (array.chain_id == target_chain)
        if not np.any(mask):
            return None
        
        chain = array[mask]
        res_types = []
        
        for residue in struc.residue_iter(chain):
            name = residue.res_name[0]
            if name in ["GLY", "DGLY", "DGL"]:
                res_types.append("G")
                continue
            
            # Extract coordinates for N, CA, C, CB
            coords = {n: None for n in ["N", "CA", "C", "CB"]}
            for i, atom_name in enumerate(residue.atom_name):
                if atom_name in coords:
                    coords[atom_name] = residue.coord[i]
            
            if any(v is None for v in coords.values()):
                res_types.append("None")
                continue
            
            # Chirality calculation using vector cross product
            vN = coords["N"] - coords["CA"]
            vC = coords["C"] - coords["CA"]
            vB = coords["CB"] - coords["CA"]
            normal = np.cross(vN, vC)
            res_types.append("D" if float(np.dot(normal, vB)) < 0.0 else "L")
            
        return {
            "file": os.path.basename(pdb_path),
            "types": res_types,
            "l_count": res_types.count("L")
        }
    except:
        return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-c", "--chain", required=True)
    parser.add_argument("-p", "--processes", type=int, default=4)
    args = parser.parse_args()

    files = []
    for root, _, filenames in os.walk(args.input):
        for f in filenames:
            if f.endswith(".pdb"):
                files.append(os.path.join(root, f))

    print(f"Processing {len(files)} files for chain {args.chain}...")
    
    tasks = [(f, args.chain) for f in files]
    with Pool(args.processes) as pool:
        results = [r for r in tqdm(pool.imap(analyze_pdb, tasks), total=len(tasks)) if r]

    stats = {"L": 0, "D": 0, "G": 0, "None": 0}
    l_heavy_files = []

    for r in results:
        for t in r["types"]:
            if t in stats: stats[t] += 1
        if r["l_count"] > 0:
            l_heavy_files.append(r)

    total = sum(stats.values())
    if total > 0:
        print(f"\nResults for Chain {args.chain}:")
        print(f"Total Residues: {total}")
        for k, v in stats.items():
            print(f"  {k}: {v} ({v/total*100:.2f}%)")
        
        ld = stats["L"] + stats["D"]
        if ld > 0:
            print(f"L/(L+D) Ratio: {stats['L']/ld*100:.2f}%")
        
        if l_heavy_files:
            print(f"\nFiles containing L-amino acids ({len(l_heavy_files)} total):")
            for f in l_heavy_files[:10]:
                print(f"  {f['file']}: {f['l_count']} L-res")
    else:
        print("No valid residues found for the specified chain.")

if __name__ == "__main__":
    main()