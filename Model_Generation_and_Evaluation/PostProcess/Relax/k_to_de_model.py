#!/usr/bin/python
# -*- coding:utf-8 -*-
import os
import time
import io
import logging
import numpy as np
import pdbfixer
import openmm
from openmm.app import Modeller
from openmm import app as openmm_app
from openmm import unit
ENERGY = unit.kilocalories_per_mole
LENGTH = unit.angstroms

from base_model import ForceFieldMinimizer


class ForceFieldMinimizerKtoDE(ForceFieldMinimizer):

    def _fix_cyclic(self, fixer, cyclic_chains, cyclic_opts):

        assert cyclic_opts is not None, f'cyclic_opts should not be None, but list of pairs ((chain_id, res_pos), (chain_id, res_pos))'

        main_chain_atoms = ['N', 'CA', 'C', 'O', 'OXT']
        
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

        #STEP1 -- Remove sidechain atoms and Rename residues to Lys and Asp
        modeller = Modeller(fixer.topology, fixer.positions)
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

        #STEP3 -- Delete leaving atoms, and save this to a new fixer (must be a new one, otherwise cyclic residue cannot be recognized to a template)
        for chain in modeller.topology.chains():
            if chain.id not in cyclic_chains: continue
            atoms_to_remove = []
            for i, res in enumerate(chain.residues()):
                if (chain.id, i) not in cyc_sites:
                    continue
                atom_names = delete_atoms[res.name]
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
                atom_name = connect_atoms[residue.name]
                for atom in residue.atoms():
                    if atom.name == atom_name: resid2connectatom[resid] = atom
        connects = []
        for res1, res2 in cyclic_opts:
            a1, a2 = resid2connectatom[res1], resid2connectatom[res2]
            connects.append('CONECT' + str(a1.id).rjust(5) + str(a2.id).rjust(5))
            connects.append('CONECT' + str(a2.id).rjust(5) + str(a1.id).rjust(5))
        pdb_fixed = self._add_connects(pdb_fixed, connects)

        return pdb_fixed, connects


if __name__ == '__main__':
    import sys
    force_field = ForceFieldMinimizerKtoDE()
    force_field(sys.argv[1], sys.argv[2], cyclic_chains=['L'], cyclic_opts=[(('L', 1), ('L', 10))]) 
    # Careful, the first index 0 and the last index are ACE and NME!!