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

class ForceFieldMinimizerHeadTail(ForceFieldMinimizer):

    def _fix_cyclic(self, fixer, cyclic_chains, cyclic_opts):
        #STEP1 -- Remove terminal atoms
        modeller = Modeller(fixer.topology, fixer.positions)
        for chain in modeller.topology.chains():
            if chain.id not in cyclic_chains: continue
            atoms_to_remove = []
            for i, res in enumerate(chain.residues()):
                if i == 0:
                    for atom in res.atoms():
                        if atom.name == 'H2' or atom.name == 'H3':
                            atoms_to_remove.append(atom)
                elif i == len(chain) - 1:
                    for atom in res.atoms():
                        if atom.name == 'OXT': atoms_to_remove.append(atom)
            modeller.delete(atoms_to_remove)

        #STEP2 -- Fix residues in new names by PDBFixer
        fixer.topology = modeller.topology
        fixer.positions = modeller.positions

        out_handle = io.StringIO()
        openmm_app.PDBFile.writeFile(fixer.topology, fixer.positions, out_handle, keepIds=True)
        pdb_fixed = out_handle.getvalue()
        new_fixer = pdbfixer.PDBFixer(pdbfile=io.StringIO(pdb_fixed))

        #STEP3 -- Find two connecting atoms and add a CONECT in between
        connects = []
        for chain in new_fixer.topology.chains():
            if chain.id not in cyclic_chains: continue
            n_term, c_term = None, None
            for i, res in enumerate(chain.residues()):
                if i == 0:
                    for atom in res.atoms():
                        if atom.name == 'N': n_term = atom.id
                elif i == len(chain) - 1:
                    for atom in res.atoms():
                        if atom.name == 'C': c_term = int(atom.id)
            connects.append('CONECT' + str(n_term).rjust(5) + str(c_term).rjust(5))
            connects.append('CONECT' + str(c_term).rjust(5) + str(n_term).rjust(5))
        
        pdb_fixed = self._add_connects(pdb_fixed, connects)
        return pdb_fixed, connects


if __name__ == '__main__':
    import sys
    force_field = ForceFieldMinimizerHeadTail()
    force_field(sys.argv[1], sys.argv[2], cyclic_chains=['L']) # starts from 0, the i-th residue
    