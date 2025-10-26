#!/usr/bin/python
# -*- coding:utf-8 -*-
import os
import time
import io
import logging
import numpy as np
import pdbfixer
import openmm
from openmm.app import Modeller, Topology
from openmm import app as openmm_app
from openmm import unit
ENERGY = unit.kilocalories_per_mole
LENGTH = unit.angstroms

# base_0424.py needs to be in the same directory or accessible in the python path
from base_0424 import ForceFieldMinimizer


def add_cb(input_array):
    #from protein mpnn
    #The virtual Cβ coordinates were calculated using ideal angle and bond length definitions: b = Cα - N, c = C - Cα, a = cross(b, c), Cβ = -0.58273431*a + 0.56802827*b - 0.54067466*c + Cα.
    N, CA, C, O = input_array
    b = CA - N
    c = C - CA
    a = np.cross(b, c)
    CB = np.around(-0.58273431*a + 0.56802827*b - 0.54067466*c + CA, 3)
    return CB

class ForceFieldMinimizerKtoDE(ForceFieldMinimizer):

    def _fix_cyclic(self, fixer, cyclic_chains, cyclic_opts):

        assert cyclic_opts is not None, f'cyclic_opts should not be None, but list of pairs ((chain_id, res_pos), (chain_id, res_pos))'

        main_chain_atoms = ['N', 'CA', 'C', 'O', 'OXT', 'CB'] # important: we need CB to maintain chiralty
        
        cyc_sites = {}
        for resid1, resid2 in cyclic_opts:
            cyc_sites[resid1] = 1
            cyc_sites[resid2] = 1
        
        delete_atoms = {
            'LYS': ['HZ2', 'HZ3'],  # K
            'ASP': ['OD2'],         # D
            'GLU': ['OE2'],         # E
        }
        connect_atoms = {
            'LYS': 'NZ',
            'ASP': 'CG',
            'GLU': 'CD'
        }

        # STEP 0: add CB for Gly
        modeller = Modeller(fixer.topology, fixer.positions)
        if cyclic_opts: 
            gly_to_modify = set()
            end_residues = [cyclic_opts[0][0], cyclic_opts[0][1]]
            for resid_tuple in end_residues:
                chain_id, res_idx_in_chain = resid_tuple
                for chain in modeller.topology.chains():
                    if chain.id == chain_id:
                        for res_idx, res in enumerate(chain.residues()):
                            if res_idx == res_idx_in_chain and res.name == 'GLY':
                                gly_to_modify.add(res.index)
                                break
                        break

            if gly_to_modify:
                new_topology = Topology()
                new_positions = [] * LENGTH
                atom_map = {} 

                for old_chain in modeller.topology.chains():
                    new_chain = new_topology.addChain(old_chain.id)
                    for old_res in old_chain.residues():
                        new_res = new_topology.addResidue(old_res.name, new_chain, old_res.id)

                        for old_atom in old_res.atoms():
                            new_atom = new_topology.addAtom(old_atom.name, old_atom.element, new_res)
                            atom_map[old_atom] = new_atom
                            new_positions.append(modeller.positions[old_atom.index])

                        if old_res.index in gly_to_modify:
                            try:
                                # Get original atoms for coordinate calculation
                                n_atom = next(a for a in old_res.atoms() if a.name == 'N')
                                ca_atom = next(a for a in old_res.atoms() if a.name == 'CA')
                                c_atom = next(a for a in old_res.atoms() if a.name == 'C')
                                o_atom = next(a for a in old_res.atoms() if a.name == 'O')
                                
                                n_pos = modeller.positions[n_atom.index]
                                ca_pos = modeller.positions[ca_atom.index]
                                c_pos = modeller.positions[c_atom.index]
                                o_pos = modeller.positions[o_atom.index]

                                cb_np = add_cb([
                                    n_pos.value_in_unit(LENGTH),
                                    ca_pos.value_in_unit(LENGTH),
                                    c_pos.value_in_unit(LENGTH),
                                    o_pos.value_in_unit(LENGTH)
                                ])
                                cb_vec = openmm.Vec3(cb_np[0], cb_np[1], cb_np[2]) * LENGTH

                                new_cb_atom = new_topology.addAtom('CB', openmm_app.Element.getBySymbol('C'), new_res)
                                atom_map[f'cb_{old_res.index}'] = new_cb_atom
                                new_positions.append(cb_vec)

                            except StopIteration:
                                logging.warning(f"Residue {old_res.id} is GLY but missing main chain atoms, cannot add CB")

                for bond in modeller.topology.bonds():
                    if bond.atom1 in atom_map and bond.atom2 in atom_map:
                        new_topology.addBond(atom_map[bond.atom1], atom_map[bond.atom2])
                
                modeller = Modeller(new_topology, new_positions)

        #STEP1 -- Remove sidechain atoms and Rename residues to Lys and Asp
        for chain in modeller.topology.chains():
            if chain.id not in cyclic_chains: continue
            atoms_to_remove = []
            for i, res in enumerate(chain.residues()):
                if (chain.id, i) not in cyc_sites:
                    continue
                for atom in res.atoms():
                    if atom.name not in main_chain_atoms:
                        atoms_to_remove.append(atom)
                if (chain.id, i) == cyclic_opts[0][0]:
                    res.name = 'LYS'
                elif (chain.id, i) == cyclic_opts[0][1]:
                    res.name = 'ASP'
            modeller.delete(atoms_to_remove)

        #STEP2 -- Fix residues in new names by PDBFixer, and send it back to modeller
        fixer.topology = modeller.topology
        fixer.positions = modeller.positions
        fixer.findMissingResidues()
        fixer.findMissingAtoms() 
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=7.0)
        modeller = Modeller(fixer.topology, fixer.positions)

        #STEP3 -- Delete leaving atoms, and save this to a new fixer
        for chain in modeller.topology.chains():
            if chain.id not in cyclic_chains: continue
            atoms_to_remove = []
            for i, res in enumerate(chain.residues()):
                if (chain.id, i) not in cyc_sites:
                    continue
                atom_names = delete_atoms.get(res.name, [])
                for atom in res.atoms():
                    if atom.name in atom_names:
                        atoms_to_remove.append(atom)
            modeller.delete(atoms_to_remove)
        
        fixer.topology = modeller.topology
        fixer.positions = modeller.positions
        out_handle = io.StringIO()
        openmm_app.PDBFile.writeFile(fixer.topology, fixer.positions, out_handle, keepIds=True)
        pdb_fixed = out_handle.getvalue()
        new_fixer = pdbfixer.PDBFixer(pdbfile=io.StringIO(pdb_fixed))

        #STEP4 -- Find two connecting atoms and add a CONECT in between
        resid2connectatom = {}
        for chain in new_fixer.topology.chains():
            if chain.id not in cyclic_chains: continue
            for i, residue in enumerate(chain.residues()):
                resid = (chain.id, i)
                if resid not in cyc_sites: continue
                atom_name = connect_atoms.get(residue.name)
                if atom_name:
                    for atom in residue.atoms():
                        if atom.name == atom_name: resid2connectatom[resid] = atom
        connects = []
        for res1, res2 in cyclic_opts:
            if res1 in resid2connectatom and res2 in resid2connectatom:
                a1, a2 = resid2connectatom[res1], resid2connectatom[res2]
                connects.append('CONECT' + str(a1.id).rjust(5) + str(a2.id).rjust(5))
                connects.append('CONECT' + str(a2.id).rjust(5) + str(a1.id).rjust(5))
        pdb_fixed = self._add_connects(pdb_fixed, connects)

        return pdb_fixed, connects