#!/usr/bin/python
# -*- coding:utf-8 -*-
from Bio.PDB import PDBParser, Chain, Model, Structure, Residue, Atom
from utils import add_cb, hydrophobic_residues, AA3TO1
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from freesasa import calcBioPDB
from io import StringIO
from copy import deepcopy
import biotite.structure.io.pdb as bpdb
from biotite.application.dssp import DsspApp
import networkx as nx
import gc

def find_cyclic_peptide(
        pdb_string,
        ss_info=None,
        min_loop_len = 4,
        max_loop_len = 19, 
        disulfide_cutoff = 4.5,
        min_BSA = 400, #https://www.sciencedirect.com/science/article/pii/S0959440X1930034X
        max_rBSA = 0.85, 
        max_relncBSA = 0.3,
        good_boy_cutoff = 0.01,
        max_sheet_percentage = 0.34,
        max_helix_percentage = 0.67,
        connectivity_cutoff = 10,
        ):
    
    # STEP 0 -- Get PDB from PDB_string
    pdb_handle = StringIO(pdb_string)
    p = PDBParser()
    chain = p.get_structure('0', pdb_handle)[0]['A'] #only deal with A chain

    # STEP 1 -- Calculate secondary structure of the whole input
    if ss_info is None:
        with StringIO(pdb_string) as f:
            pdb_file = bpdb.PDBFile.read(f)
            atom_array = pdb_file.get_structure()

        app = DsspApp(atom_array[0], bin_path='/data_hdd/home/xiehanyuan/miniconda3/envs/my_env/bin/mkdssp')
        app.start()
        app.join()
        ss_array = app.get_sse()
        ss_info = ''.join(ss_array)

    # STEP 2 -- Get CB coordination (add a virtual CB for GLY)
    cb_coord = []
    seq = ''
    plddt = []
    res_ids = []
    
    for res in chain:
        seq += AA3TO1[res.get_resname()]
        plddt.append(res['CA'].get_bfactor())
        res_ids.append(res.get_id()[1])  # Store residue index in PDB file
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
    
    # STEP 3 -- Calculate CB contact matrix and possible cyclic positions
    cb_contact = np.linalg.norm(cb_coord[None,:,:,] - cb_coord[:,None,:], axis=-1)
    possible_pair = (cb_contact >= 3) & (cb_contact <= 15)
    possible_pair = np.triu(np.tril(possible_pair, max_loop_len), min_loop_len)
    residue_pair = np.where(possible_pair)
    disconnect_index = 0
    accepted = []

    # STEP 4 -- Discuss every possible cyclic pair
    for i, j in zip(residue_pair[0], residue_pair[1]):

        # Initial pool of candidates

        # Define receptor
        min_dist = np.min(cb_contact[i:j+1], axis=0)
        min_dist[max(i-5,0):min(j+6,len(seq))] = 21
        neighbors_20A = np.where(min_dist < 20)[0]

        # pLDDT filter
        worst_pLDDT = min(plddt[i:j+1])

        if worst_pLDDT < 70:
            continue
        
        # Breakpoint filter
        if res_ids[j] - res_ids[i] != j - i:
            continue

        # Existing disulfide filter
        loop_residues = set(range(i+1, j)) # Don't need to check the termini
        neighbor_set = set(neighbors_20A)
        combined_set = loop_residues.union(neighbor_set)
        cys_positions = [x for x in combined_set if seq[x] == 'C']

        has_disulfide = False
        for idx, x in enumerate(cys_positions):
            for y in cys_positions[idx+1:]:
                if cb_contact[x,y] <= disulfide_cutoff:
                    if (x in loop_residues) or (y in loop_residues):
                        has_disulfide = True
                        break
            if has_disulfide:
                break
        if has_disulfide:
            continue
        
        # Secondary structure filter
        peptide_ss = ss_info[i:j+1]
        peptide_helix = peptide_ss.count("H") + peptide_ss.count("G") + peptide_ss.count("I")
        peptide_sheet = peptide_ss.count("E")+peptide_ss.count("B")
        peptide_helix_ratio = peptide_helix / len(peptide_ss)
        peptide_sheet_ratio = peptide_sheet / len(peptide_ss)
        
        if peptide_sheet_ratio > max_sheet_percentage or peptide_helix_ratio > max_helix_percentage:
            continue

        # Hydrophobic filter
        pep_seq = seq[i:j+1]
        prot_param = ProteinAnalysis(pep_seq)
        aa_percent = prot_param.get_amino_acids_percent()
        hydrophobic_ratio = sum([aa_percent[k] for k in hydrophobic_residues])

        if hydrophobic_ratio > 0.45:
            continue

        # Build structure and get BSA
        good_boy = False
        receptor_chain = Chain.Chain('R') 
        ligand_chain = Chain.Chain('L') 

        resids_receptor = [] # For connectivity check, get indexes for edge matrix
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
        
        if abs_bsa <= min_BSA or rel_nc_bsa >= max_relncBSA or rel_bsa >= max_rBSA:
            continue
        if rel_nc_bsa <= good_boy_cutoff:
            good_boy = True
        
        # Connectivity Filter
        submatrix = cb_contact[np.ix_(resids_receptor,resids_receptor)]
        connectivity = np.where(submatrix <= connectivity_cutoff, 1, 0)
        num_nodes = connectivity.shape[0]
        rows, cols = np.nonzero(np.triu(connectivity, k=1))
        edges = zip(rows, cols)
        graph = nx.Graph()
        graph.add_nodes_from(range(num_nodes))
        graph.add_edges_from(edges)
        if not nx.is_connected(graph):
            continue

        # Add terminal structures by cut them off from adjecent residues, if it is not for head-tail cyclization
        cb_distance = cb_contact[i,j]
        if (3 <= cb_distance < 4.5) or (6 <= cb_distance <= 8):
            overall_list = list(chain) 
            ligand_list = list(ligand_chain)
            ligand_index = res_ids[i:j+1]

            ACE_id = (' ', min(ligand_index)-1, ' ')
            NME_id = (' ', max(ligand_index)+1, ' ')

            ACE = Residue.Residue(ACE_id, "ACE", 0) 
            NME = Residue.Residue(NME_id, "NME", 0) # Add ACE and NME as two residues in ligand chain

            # If the cyclization sites are in terminal residues already, throw them away
            if i - 1 >= 0:
                prev_res = overall_list[i - 1]
                for atom_name in ['CA', 'C', 'O']:
                    atom = prev_res[atom_name]
                    new_atom = Atom.Atom(atom_name, atom.get_coord(), atom.get_bfactor(),
                             atom.get_occupancy(), atom.get_altloc(), atom.get_fullname(),
                             atom.get_serial_number(), element=atom.element)
                    ACE.add(new_atom)
            else: 
                continue

            if j + 1 < len(chain):
                next_res = overall_list[j + 1]
                for atom_name in ['N', 'CA']:
                    atom = next_res[atom_name]
                    new_atom = Atom.Atom(atom_name, atom.get_coord(), atom.get_bfactor(),
                             atom.get_occupancy(), atom.get_altloc(), atom.get_fullname(),
                             atom.get_serial_number(), element=atom.element)
                    NME.add(new_atom)
            else: 
                continue

            new_ligand_chain = Chain.Chain('L')
            new_ligand_chain.add(ACE.copy())
            for res in ligand_chain:
                new_ligand_chain.add(res)
            new_ligand_chain.add(NME.copy())

        #prepare for output
        output_structure = Structure.Structure('output')
        output_model = Model.Model(0)
        output_structure.add(tmp_model)
        if (3 <= cb_distance < 4.5) or (6 <= cb_distance <= 8):
            output_model.add(new_ligand_chain)
        else:
            output_model.add(ligand_chain)
        output_model.add(receptor_chain) # Is is necessary to re-generate the structure and model?
        output_structure = deepcopy(output_model)

        length = j-i+1

        cyclization_struct = None
        if cb_distance >= 3 and cb_distance < 4.5: cyclization_struct = 'DISULFIDE'
        elif cb_distance >= 4.5 and cb_distance < 6: cyclization_struct = 'HEADTAIL'
        elif cb_distance >= 6 and cb_distance <= 8: cyclization_struct = 'ISOPEPTIDE'
        
        accepted.append((res_ids[i],res_ids[j],pep_seq,length,cb_distance,worst_pLDDT,abs_bsa,rel_bsa,rel_nc_bsa,hydrophobic_ratio,peptide_helix_ratio,peptide_sheet_ratio,peptide_ss,output_structure,cyclization_struct,good_boy))

        del graph,output_model,output_structure,connectivity,min_dist,tmp_structure,tmp_model
    del cb_contact, possible_pair,app
    gc.collect()

    return accepted
