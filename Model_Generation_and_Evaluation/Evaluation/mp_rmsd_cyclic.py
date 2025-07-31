import os
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, Select, PDBIO
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
import matplotlib.pyplot as plt
import tempfile
import shutil
import glob
import argparse

from tqdm import tqdm
from multiprocessing import Pool
from functools import partial

def rmsd_pair_analysis(args):
    i, j, path_i, path_j = args
    try:
        u1 = mda.Universe(path_i)
        u2 = mda.Universe(path_j)
        atoms1 = u1.select_atoms("protein and name CA")
        atoms2 = u2.select_atoms("protein and name CA")
        if len(atoms1) != len(atoms2):
            return (i, j, np.nan)
        rmsd_analysis = RMSD(u2, u1, select='name CA', ref_frame=0)
        rmsd_analysis.run()
        rms = rmsd_analysis.results.rmsd[-1, 2]
        return (i, j, rms)
    except Exception as e:
        print(f"Error calculating RMSD {path_i} vs {path_j}: {e}")
        return (i, j, np.nan)

from Bio.SVDSuperimposer import SVDSuperimposer
def rmsd_pair(args):
    i, j, coords_i, coords_j = args
    if coords_i is None or coords_j is None or coords_i.shape[0] != coords_j.shape[0]:
        return (i, j, np.inf)
    try:
        def compute_rmsd(a, b):
            sup = SVDSuperimposer()
            sup.set(a, b)
            sup.run()
            return sup.get_rms()
        def generate_permutations(coords):
            return [np.vstack((coords[k:], coords[:k])) for k in range(len(coords))]
        rmsds = [
            compute_rmsd(pi, pj)
            for pi in generate_permutations(coords_i)
            for pj in generate_permutations(coords_j)
        ]
        return (i, j, min(rmsds))
    except Exception as e:
        print(f"Error calculating RMSD {i} vs {j}: {e}")
        return (i, j, np.inf)

class AChainSelect(Select):
    def accept_chain(self, chain):
        return chain.id == 'L'

def extract_chain_L(input_pdb, output_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb", input_pdb)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, select=AChainSelect())

def extract_one_chain(pdb_path_i_tuple, temp_dir):
    i, pdb_path = pdb_path_i_tuple
    try:
        name = os.path.basename(pdb_path).replace(".pdb", f"_{i}.pdb")
        output_path = os.path.join(temp_dir, name)
        extract_chain_L(pdb_path, output_path)
        return (output_path, os.path.basename(pdb_path))
    except Exception as e:
        print(f"Error processing {pdb_path}: {e}")
        return None

def get_ca_coords(pdb_path):
    u = mda.Universe(pdb_path)
    atoms = u.select_atoms("protein and name CA")
    return atoms.positions.copy()

def main(args):
    rmsd_npy = args.rmsd
    hist_png = args.hist
    cluster_csv = args.cluster
    num_cores = args.num_cores
    precomputed_rmsd = args.precomputed_rmsd

    name_list_path = args.name_list
    with open(name_list_path) as f:
        pdb_paths = [line.strip() for line in f]
    
    valid_pdb_paths = []
    for path in pdb_paths:
        if os.path.exists(path):
            valid_pdb_paths.append(path)
        else:
            print(f"Warning: File does not exist: {path}")

    if not valid_pdb_paths:
        print("Error: No valid PDB files found!")
        return

    print(f"Found {len(valid_pdb_paths)} valid PDB files")

    temp_dir = tempfile.mkdtemp()

    with Pool(processes=num_cores) as pool:
        results = list(tqdm(pool.imap(partial(extract_one_chain, temp_dir=temp_dir), enumerate(valid_pdb_paths)), total=len(valid_pdb_paths), desc="Extracting L chains"))

    extracted_paths = []
    valid_names = []
    for res in results:
        if res is not None:
            extracted_paths.append(res[0])
            valid_names.append(res[1])

    if not extracted_paths:
        print("Error: Failed to extract any chains!")
        shutil.rmtree(temp_dir)
        return

    print(f"Successfully extracted {len(extracted_paths)} chains")

    n = len(extracted_paths)

    if precomputed_rmsd:
        print(f"Loading precomputed RMSD matrix from {precomputed_rmsd}")
        min_rmsd_matrix = np.load(precomputed_rmsd)
    else:
        rmsd_matrix = np.zeros((n, n))
        tasks = [(i, j, extracted_paths[i], extracted_paths[j]) for i in range(n) for j in range(i + 1, n)]
        with Pool(processes=num_cores) as pool:
            for i, j, rms in tqdm(pool.imap(rmsd_pair_analysis, tasks), total=len(tasks), desc="Computing RMSD matrix"):
                rmsd_matrix[i, j] = rmsd_matrix[j, i] = rms

        rmsd_matrix = np.nan_to_num(rmsd_matrix, nan=100.0)

        ca_coords_list = []
        for path in extracted_paths:
            try:
                coords = get_ca_coords(path)
                ca_coords_list.append(coords)
            except Exception as e:
                print(f"Failed to extract CA for {path}: {e}")
                ca_coords_list.append(None)

        tasks = [(i, j, ca_coords_list[i], ca_coords_list[j]) for i in range(n) for j in range(i + 1, n)]
        min_rmsd_matrix = np.full((n, n), np.inf)
        with Pool(processes=num_cores) as pool:
            for i, j, val in tqdm(pool.imap(rmsd_pair, tasks), total=len(tasks), desc="Calculating circular RMSD"):
                min_rmsd_matrix[i, j] = min_rmsd_matrix[j, i] = val

        np.save(rmsd_npy, min_rmsd_matrix)

    upper_vals = min_rmsd_matrix[np.triu_indices(n, k=1)]
    plt.hist(upper_vals[upper_vals < np.inf], bins=50)
    plt.xlabel("Min RMSD (Å)")
    plt.ylabel("Frequency")
    plt.title("Circular-permutation min RMSD distribution")
    plt.savefig(hist_png, dpi=300, bbox_inches='tight')
    plt.close()

    ca_coords_list = []
    for path in extracted_paths:
        try:
            coords = get_ca_coords(path)
            ca_coords_list.append(coords)
        except Exception as e:
            print(f"Failed to extract CA for {path}: {e}")
            ca_coords_list.append(None)
    lengths = [coords.shape[0] if coords is not None else -1 for coords in ca_coords_list]

    threshold = 1.0
    unclustered = list(range(n))
    clusters = []
    for _ in tqdm(range(n), desc="Greedy clustering from RMSD matrix"):
        if not unclustered:
            break
        i = unclustered[0]
        cluster = [i]
        for j in unclustered[1:]:
            if min_rmsd_matrix[i, j] < threshold and lengths[i] == lengths[j]:
                cluster.append(j)
        clusters.append(cluster)
        unclustered = [u for u in unclustered if u not in cluster]

    labels = [0] * len(valid_names)
    for cid, cluster in enumerate(clusters, start=1):
        for i in cluster:
            labels[i] = cid

    result_df = pd.DataFrame({
        "pdb_path": valid_names,
        "cluster": labels
    })
    result_df.to_csv(cluster_csv, index=False)

    print(f"Clustering results saved to: {cluster_csv}")
    print(f"RMSD matrix saved to: {rmsd_npy}")
    print(f"Histogram saved to: {hist_png}")

    shutil.rmtree(temp_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RMSD clustering of cyclic peptides (L chain)")
    parser.add_argument("--rmsd", required=True, help="Output path for RMSD .npy matrix")
    parser.add_argument("--hist", required=True, help="Output path for histogram image")
    parser.add_argument("--cluster", required=True, help="Output path for cluster CSV")
    parser.add_argument("--name_list", required=True, help="Path to name.list file containing full PDB paths")
    parser.add_argument("--num_cores", type=int, default=os.cpu_count(), help="Number of CPU cores to use")
    parser.add_argument("--precomputed_rmsd", default=None, help="Path to precomputed RMSD .npy matrix")
    args = parser.parse_args()
    main(args)