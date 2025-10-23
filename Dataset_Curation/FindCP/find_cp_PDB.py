#!/usr/bin/python
# -*- coding:utf-8 -*-
from Bio.PDB import PDBParser, Chain, Model, Structure, Residue, Atom
from utils import add_cb, restype_name_to_atom14_names, hydrophobic_residues, AA3TO1
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
        min_loop_len = 4,
        max_loop_len = 15, 
        disulfide_cutoff = 4.5,
        min_BSA = 400, #https://www.sciencedirect.com/science/article/pii/S0959440X1930034X
        max_rBSA = 0.85, 
        max_relncBSA = 0.3,
        good_boy_cutoff = 0.01,
        max_sheet_percentage = 0.34,
        max_helix_percentage = 0.67,
        connectivity_cutoff = 10,
        ):
    
    try:
        pdb_handle = StringIO(pdb_string)
        p = PDBParser()
        struct = p.get_structure('0', pdb_handle)
        model = struct[0] 
        
        all_residues = []
        residue_info = [] 
        global_index = 0
        
        for chain in model:
            chain_id = chain.id
            for res in chain:
                if res.get_id()[0] == ' ':
                    all_residues.append(res)
                    residue_info.append((global_index, chain_id, res.get_id()[1], res))
                    global_index += 1
        total_residues = len(all_residues)
        
        if total_residues < min_loop_len + 3:
            return []
        
        with StringIO(pdb_string) as f:
            pdb_file = bpdb.PDBFile.read(f)
            full_atom_array = pdb_file.get_structure()[0]

        ss_info = [''] * total_residues
        
        for chain_id in set(info[1] for info in residue_info):
            chain_indices = [i for i, info in enumerate(residue_info) if info[1] == chain_id]
            if not chain_indices:
                continue
                
            chain_atom_array = full_atom_array[full_atom_array.chain_id == chain_id]
            if len(chain_atom_array) == 0:
                continue
                
            try:
                app = DsspApp(chain_atom_array, bin_path='/data_hdd/home/xiehanyuan/miniconda3/envs/my_env/bin/mkdssp')
                app.start()
                app.join()
                chain_ss = app.get_sse()
                
                for idx, ss in zip(chain_indices, chain_ss):
                    ss_info[idx] = ss
            except Exception as e:
                print(f"Error calculating DSSP for chain {chain_id}: {e}")
                for idx in chain_indices:
                    ss_info[idx] = ' '
        
        ss_info = ''.join(ss_info)
        
        cb_coord = []
        seq = ''
        
        for res in all_residues:
            res_name = res.get_resname()
            seq += AA3TO1.get(res_name, 'X')
            try:
                cb_coord.append(res['CB'].get_coord())
            except:
                try:
                    tmp_coord = np.array([
                        res['N'].get_coord(),
                        res['CA'].get_coord(),
                        res['C'].get_coord(),
                        res['O'].get_coord()
                    ])
                    cb_coord.append(add_cb(tmp_coord))
                except KeyError:
                    cb_coord.append(res['CA'].get_coord())
        
        cb_coord = np.array(cb_coord)
        
        cb_contact = np.linalg.norm(cb_coord[None, :, :] - cb_coord[:, None, :], axis=-1)
        
        possible_pair = (cb_contact >= 3) & (cb_contact <= 8)
        possible_pair = np.triu(np.tril(possible_pair, max_loop_len), min_loop_len)
        residue_pair = np.where(possible_pair)
        accepted = []
        
        for i, j in zip(residue_pair[0], residue_pair[1]):
            if residue_info[i][1] != residue_info[j][1]:
                continue
            expected_ids = list(range(residue_info[i][2], residue_info[j][2] + 1))
            actual_ids = [info[2] for info in residue_info[i:j+1]]
            
            if expected_ids != actual_ids:
                continue
            
            chain_id = residue_info[i][1]
            
            min_dist = np.min(cb_contact[i:j+1], axis=0)
            min_dist[max(i-5,0):min(j+6,total_residues)] = 21
            neighbors_20A = np.where(min_dist < 20)[0]
            
            valid_atom_composition = True
            check_residues = set(range(i, j+1)) | set(neighbors_20A)
            for k in check_residues:  
                res = all_residues[k]
                res_name = res.get_resname()  
                
                if res_name not in restype_name_to_atom14_names:
                    valid_atom_composition = False
                    break

                expected_atoms = [atom for atom in restype_name_to_atom14_names[res_name] if atom]
                
                actual_atoms = [atom.get_id() for atom in res.get_atoms()]

                if set(actual_atoms) != set(expected_atoms):
                    valid_atom_composition = False
                    break
            if not valid_atom_composition:
                continue
            
            loop_residues = set(range(i+1, j))  
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
            
            peptide_ss = ss_info[i:j+1]
            peptide_helix = peptide_ss.count("H") + peptide_ss.count("G") + peptide_ss.count("I")
            peptide_sheet = peptide_ss.count("E")+peptide_ss.count("B")
            peptide_helix_ratio = peptide_helix / len(peptide_ss) if len(peptide_ss) > 0 else 0
            peptide_sheet_ratio = peptide_sheet / len(peptide_ss) if len(peptide_ss) > 0 else 0
            
            if peptide_sheet_ratio > max_sheet_percentage or peptide_helix_ratio > max_helix_percentage:
                continue

            pep_seq = seq[i:j+1]
            prot_param = ProteinAnalysis(pep_seq)
            aa_percent = prot_param.get_amino_acids_percent()
            hydrophobic_ratio = sum([aa_percent[k] for k in hydrophobic_residues])

            if hydrophobic_ratio > 0.45:
                continue

            good_boy = False
            receptor_chain = Chain.Chain('R') 
            ligand_chain = Chain.Chain('L') 

            resids_receptor = [] 
            for k in range(total_residues):
                if k>=i and k<=j:
                    ligand_chain.add(all_residues[k].copy())
                elif k in neighbors_20A:
                    receptor_chain.add(all_residues[k].copy())
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

            abs_bsa = sum(unbounded_SASA[1:-1]) - sum(bounded_SASA[1:-1]) if len(unbounded_SASA) > 2 and len(bounded_SASA) > 2 else 0
            rel_bsa = abs_bsa/sum(unbounded_SASA[1:-1]) if sum(unbounded_SASA[1:-1]) > 0 else 0
            rel_nc_bsa = (unbounded_SASA[0]+unbounded_SASA[-1]-bounded_SASA[0]-bounded_SASA[-1])/(unbounded_SASA[0]+unbounded_SASA[-1]) if (unbounded_SASA[0]+unbounded_SASA[-1]) > 0 else 0
            
            if abs_bsa <= min_BSA or rel_nc_bsa >= max_relncBSA or rel_bsa >= max_rBSA:
                continue
            if rel_nc_bsa <= good_boy_cutoff:
                good_boy = True
            
            if len(resids_receptor) == 0:
                continue 
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

            cb_distance = cb_contact[i,j]
            new_ligand_chain = None
            if (3 <= cb_distance < 4.5) or (6 <= cb_distance <= 8):
                ligand_index = [residue_info[k][2] for k in range(i, j+1)]
                
                chain_residues = [res for res in model[chain_id]]
                
                ACE_id = (' ', min(ligand_index)-1, ' ')
                NME_id = (' ', max(ligand_index)+1, ' ')

                ACE = Residue.Residue(ACE_id, "ACE", 0) 
                NME = Residue.Residue(NME_id, "NME", 0)

                prev_res = next((res for res in chain_residues if res.get_id()[1] == min(ligand_index)-1), None)
                if prev_res:
                    for atom_name in ['CA', 'C', 'O']:
                        try:
                            atom = prev_res[atom_name]
                            new_atom = Atom.Atom(atom_name, atom.get_coord(), atom.get_bfactor(),
                                    atom.get_occupancy(), atom.get_altloc(), atom.get_fullname(),
                                    atom.get_serial_number(), element=atom.element)
                            ACE.add(new_atom)
                        except KeyError:
                            pass
                else: 
                    continue

                next_res = next((res for res in chain_residues if res.get_id()[1] == max(ligand_index)+1), None)
                if next_res:
                    for atom_name in ['N', 'CA']:
                        try:
                            atom = next_res[atom_name]
                            new_atom = Atom.Atom(atom_name, atom.get_coord(), atom.get_bfactor(),
                                    atom.get_occupancy(), atom.get_altloc(), atom.get_fullname(),
                                    atom.get_serial_number(), element=atom.element)
                            NME.add(new_atom)
                        except KeyError:
                            pass
                else: 
                    continue

                new_ligand_chain = Chain.Chain('L')
                new_ligand_chain.add(ACE.copy())
                for res in ligand_chain:
                    new_ligand_chain.add(res)
                new_ligand_chain.add(NME.copy())

            output_structure = Structure.Structure('output')
            output_model = Model.Model(0)
            output_structure.add(output_model)
            if new_ligand_chain is not None:
                output_model.add(new_ligand_chain)
            else:
                output_model.add(ligand_chain)
            output_model.add(receptor_chain)
            output_structure = deepcopy(output_model)

            length = j-i+1

            cyclization_struct = None
            if 3 <= cb_distance < 4.5:
                cyclization_struct = 'DISULFIDE'
            elif 4.5 <= cb_distance < 6:
                cyclization_struct = 'HEADTAIL'
            elif 6 <= cb_distance <= 8:
                cyclization_struct = 'ISOPEPTIDE'
            
            accepted.append((
                chain_id,
                residue_info[i][2], residue_info[j][2], pep_seq, length, cb_distance,
                abs_bsa, rel_bsa, rel_nc_bsa, hydrophobic_ratio,
                peptide_helix_ratio, peptide_sheet_ratio, peptide_ss,
                output_structure, cyclization_struct, good_boy
            ))
    except Exception as e:
        print(e)
        del graph, output_model, output_structure, connectivity, min_dist, tmp_structure, tmp_model
    pdb_handle.close()
    gc.collect()

    return accepted
