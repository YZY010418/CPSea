import os
from Bio.PDB import PDBParser, Chain, Model, Structure
from utils import add_cb, charged_residues, hydrophobic_residues, aaSMILES, timer, AA3TO1
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import re
from rdkit.Chem.rdMolDescriptors import CalcTPSA
from freesasa import calcBioPDB
from rdkit.Chem import MolFromSmiles
from io import StringIO
from copy import deepcopy
import biotite.structure.io.pdb as bpdb
from biotite.application.dssp import DsspApp
import networkx as nx

class cyclic_filter():
    def __init__(self) -> None:
        self.disconnect_index = 0
    def find_cyclic_peptide(
            self,
            pdb_string,
            ss_info=None,
            min_loop_len = 4,
            max_loop_len = 15, 
            min_BSA = 400, #https://www.sciencedirect.com/science/article/pii/S0959440X1930034X
            max_relncBSA = 0.3,
            max_sheet_percentage = 0.3,
            max_helix_percentage = 0.4,
            connectivity_cutoff = 11,
            ):
        
        #Get PDB from PDB_string
        pdb_handle = StringIO(pdb_string)
        p = PDBParser()
        chain = p.get_structure('0', pdb_handle)[0]['A'] #only deal with A chain
        io = PDBIO()

        # Calculate secondary structure
        if ss_info is None:
            with StringIO(pdb_string) as f:
                pdb_file = bpdb.PDBFile.read(f)
                atom_array = pdb_file.get_structure()

            app = DsspApp(atom_array[0], bin_path='/data_hdd/home/xiehanyuan/miniconda3/envs/my_env/bin/mkdssp')
            app.start()
            app.join()
            ss_array = app.get_sse()
            ss_info = ''.join(ss_array)

        cb_coord = []
        seq = ''
        plddt = []
        res_ids = []  # Store residue IDs for breakpoint detection
        
        for res in chain:
            seq += AA3TO1[res.get_resname()]
            plddt.append(res['CA'].get_bfactor())
            res_ids.append(res.get_id()[1])  # Store residue number
            try:
                cb_coord.append(res['CB'].get_coord())
            except:
                tmp_coord = np.array([
                    res['N'].get_coord(),
                    res['CA'].get_coord(),
                    res['C'].get_coord(),
                    res['O'].get_coord()
                ])
                cb_coord.append(add_cb(tmp_coord))
                
        plddt = np.array(plddt)
        cb_coord = np.array(cb_coord)
        pdb_handle.close()
        
        # Calculate contact matrix and possible cyclic positions
        cb_contact = np.linalg.norm(cb_coord[None,:,:,] - cb_coord[:,None,:], axis=-1)
        possible_pair = (cb_contact >= 3) & (cb_contact <= 8) # based on statistical data
        possible_pair = np.triu(np.tril(possible_pair, max_loop_len), min_loop_len)
        residue_pair = np.where(possible_pair)
        disconnect_index = 0
        accepted = []
        monitor_survival = [0] * 7
        for i, j in zip(residue_pair[0], residue_pair[1]):
            # Initial pool of candidates
            monitor_survival[0] += 1

            # define receptor
            min_dist = np.min(cb_contact[i:j+1], axis=0)
            min_dist[max(i-5,0):min(j+6,len(seq))] = 21
            neighbors_20A = np.where(min_dist < 20)[0]

            # pLDDT filter
            worst_pLDDT = min(plddt[i:j+1])

            if worst_pLDDT < 70:
                monitor_survival[1] += 1
                continue
            
            # Breakpoint filter
            has_breakpoint = False
            for k in range(i, j):
                if res_ids[k+1] - res_ids[k] != 1:
                    has_breakpoint = True
                    break
            if has_breakpoint:
                monitor_survival[2] += 1
                continue
            
            # Secondary structure filter
            peptide_ss = ss_info[i:j+1]
            peptide_helix = peptide_ss.count("H") + peptide_ss.count("G") + peptide_ss.count("I")
            peptide_sheet = peptide_ss.count("E")+peptide_ss.count("B")
            peptide_helix_ratio = peptide_helix / len(peptide_ss)
            peptide_sheet_ratio = peptide_sheet / len(peptide_ss)
            peptide_ss_ratio = (peptide_helix + peptide_sheet) / len(peptide_ss)
            
            
            if peptide_sheet_ratio > max_sheet_percentage or peptide_helix_ratio > max_helix_percentage:
                monitor_survival[3] += 1
                continue


            # Hydrophobic filter
            pep_seq = seq[i:j+1]
            prot_param = ProteinAnalysis(pep_seq)
            aa_percent = prot_param.get_amino_acids_percent()
            hydrophobic_ratio = sum([aa_percent[k] for k in hydrophobic_residues])

            if hydrophobic_ratio > 0.45:
                monitor_survival[4] += 1
                continue

            #build structure and get BSA
            receptor_chain = Chain.Chain('R') 
            ligand_chain = Chain.Chain('L') 

            resids_receptor = [] # For connectivity check, this is the index in contact matrix
            for k,res in enumerate(chain):
                if k>=i and k<=j:
                    ligand_chain.add(res.copy())
                elif k in neighbors_20A:
                    receptor_chain.add(res.copy())
                    resids_receptor.append(k)
            tmp_structure = Structure.Structure('tmp')
            tmp_model = Model.Model(0)
            tmp_structure.add(tmp_model) 

            tmp_model.add(ligand_chain)
            unbounded_SASA = calcBioPDB(tmp_structure)[0].residueAreas()['L']
            unbounded_SASA = [k.total for k in unbounded_SASA.values()]

            tmp_model.add(receptor_chain)
            bounded_SASA = calcBioPDB(tmp_structure)[0].residueAreas()['L']
            bounded_SASA = [k.total for k in bounded_SASA.values()]

            abs_bsa = sum(unbounded_SASA[1:-1]) - sum(bounded_SASA[1:-1])
            rel_bsa = abs_bsa/sum(unbounded_SASA[1:-1])
            rel_nc_bsa = (unbounded_SASA[0]+unbounded_SASA[-1]-bounded_SASA[0]-bounded_SASA[-1])/(unbounded_SASA[0]+unbounded_SASA[-1])

            
            if abs_bsa <= min_BSA or rel_nc_bsa >= max_relncBSA:
                monitor_survival[5] += 1
                continue
            
            #Connectivity Filter
            All_connected = False
            submatrix = cb_contact[np.ix_(resids_receptor,resids_receptor)]
            connectivity = np.where(submatrix <= connectivity_cutoff, 1, 0)
            graph = nx.Graph()
            num_nodes = connectivity.shape[0]
            graph.add_nodes_from(range(num_nodes))
            for i in range(num_nodes):
                for j in range(i + 1, num_nodes):
                    if connectivity[i][j] == 1:
                        graph.add_edge(i, j)
            connected_components = list(nx.connected_components(graph))
            if len(connected_components) == 1:
                All_connected = True
            if All_connected == False:
                io.set_structure(deepcopy(tmp_structure))
                io.save(f"./test_output/version0408_output_test001/disconnected_11A/disconnected_{self.disconnect_index}_11A.pdb")
                self.disconnect_index += 1
                monitor_survival[6] += 1
                continue

            #prepare for output
            cb_distance = cb_contact[i,j]
            length = j-i+1
            output_structure = deepcopy(tmp_structure)
            resids_ligand = [] # For output name, this is the index in pdb structure
            for res in ligand_chain:
                resids_ligand.append(res.id[1])
            i = min(resids_ligand)
            j = max(resids_ligand)
            accepted.append((i,j,pep_seq,length,cb_distance,worst_pLDDT,has_breakpoint,abs_bsa,rel_bsa,rel_nc_bsa,hydrophobic_ratio,peptide_helix_ratio,peptide_sheet_ratio,output_structure,peptide_ss,peptide_ss_ratio))

        return accepted, monitor_survival


if __name__ == '__main__':
    import os
    import csv
    from tqdm import tqdm
    from Bio.PDB.PDBIO import PDBIO
    import warnings

    warnings.filterwarnings("ignore")

    root = '/data_hdd/home/yangziyi/Projects/CyclicPep/check_filter_function/test_input_1000'
    myfilter = cyclic_filter()
    saved = '/data_hdd/home/yangziyi/Projects/CyclicPep/test_output/version0408_output_test001'
    csv_file_path = os.path.join(saved, 'meta_data_0408_001.csv')
    survival_check_path = os.path.join(saved, 'survival_check_0408_001.txt')
    total_survival = [0] * 7

    with open(csv_file_path, 'w', newline='') as csvfile:
        fieldnames = ['input','output','i', 'j', 'pep_seq', 'length', 'cb_distance', 'worst_pLDDT', 'has_breakpoint', 'abs_bsa', 'rel_bsa', 'rel_nc_bsa', 'hydrophobic_ratio', 'peptide_helix_ratio', 'peptide_sheet_ratio',
                      'peptide_ss', 'peptide_ss_ratio']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

    for g in tqdm(os.listdir(root)):
        test = os.path.join(root,g)
        with open(test) as f:
            teststring = f.read()
        result, survival = myfilter.find_cyclic_peptide(teststring)
        for i,j,pep_seq,length,cb_distance,worst_pLDDT,has_breakpoint,abs_bsa,rel_bsa,rel_nc_bsa,hydrophobic_ratio,peptide_helix_ratio,peptide_sheet_ratio,output_structure,peptide_ss,peptide_ss_ratio in result:
            if output_structure is not None:
                io = PDBIO()
                io.set_structure(output_structure)
                io.save(os.path.join(saved,f'{g}_{i}_{j}.pdb'))
                #save metadata       
                with open(csv_file_path, 'a', newline='') as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames = fieldnames)
                    writer.writerow({
                            'input': f'{g}',
                            'output': f'{g}_{i}_{j}.pdb',
                            'i': i,
                            'j': j,
                            'pep_seq': pep_seq,
                            'length': length,
                            'cb_distance': cb_distance,
                            'worst_pLDDT': worst_pLDDT,
                            'has_breakpoint': has_breakpoint,
                            'abs_bsa': abs_bsa,
                            'rel_bsa': rel_bsa,
                            'rel_nc_bsa': rel_nc_bsa,
                            'hydrophobic_ratio': hydrophobic_ratio,
                            'peptide_helix_ratio': peptide_helix_ratio,
                            'peptide_sheet_ratio': peptide_sheet_ratio,
                            'peptide_ss': peptide_ss,
                            'peptide_ss_ratio': peptide_ss_ratio
                        })
        for i in range(len(total_survival)):
            total_survival[i] += survival[i]

    # check survival
    with open(survival_check_path, 'a') as survival_txt:
        survival_txt.write(f"initial candidate pool: {total_survival[0]}\n \
            pLDDT filter kills: {total_survival[1]}\n \
            breakpoint filter kills: {total_survival[2]}\n \
            ss_ratio filter kills: {total_survival[3]} \n \
            hydrophobic filter kills: {total_survival[4]} \n \
            BSA filters kill: {total_survival[5]} \n \
            connectivity filter kills: {total_survival[6]} \n \
            \nParameters: \n \
            CB_distance = 3-8,\n \
            neighbors_20A >= 16,\n \
            worst_pLDDT >= 70, \n \
            hydrophobic_ratio <= 0.45,\n \
            min_loop_len = 4,\n \
            max_loop_len = 15,\n  \
            min_BSA = 400,\n \
            max_relncBSA = 0.3,\n \
            max_sheet_percentage = 0.5,\n \
            max_helix_percentage = 0.4\n ")

                









