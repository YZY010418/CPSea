#!/usr/bin/python
# -*- coding:utf-8 -*-
import os
import io
import re
import numpy as np
import pdbfixer
import openmm
from openmm.app import Modeller, Topology
from openmm import app as openmm_app
from openmm import unit
import gzip

from base_0424 import ForceFieldMinimizer

ENERGY = unit.kilocalories_per_mole
LENGTH = unit.angstroms

def _ensure_l_chirality(modeller, residue):

    changed = False

    n_atom, ca_atom, c_atom, o_atom, cb_atom = None, None, None, None, None

    for atom in residue.atoms():
        if atom.name == 'N': n_atom = atom
        elif atom.name == 'CA': ca_atom = atom
        elif atom.name == 'C': c_atom = atom
        elif atom.name == 'O': o_atom = atom
        elif atom.name == 'CB': cb_atom = atom

    if not (n_atom and ca_atom and c_atom and cb_atom):
        return

    pos = modeller.positions
    ca_pos = np.array(pos[ca_atom.index].value_in_unit(LENGTH))
    n_pos = np.array(pos[n_atom.index].value_in_unit(LENGTH))
    c_pos = np.array(pos[c_atom.index].value_in_unit(LENGTH))
    cb_pos = np.array(pos[cb_atom.index].value_in_unit(LENGTH))

    vec_n = n_pos - ca_pos
    vec_c = c_pos - ca_pos
    vec_cb = cb_pos - ca_pos

    signed_volume = np.dot(np.cross(vec_n, vec_c), vec_cb)

    if signed_volume < 0:
        changed = True
        normal = np.cross(vec_n, vec_c)
        normal /= np.linalg.norm(normal)

        proj_length = np.dot(vec_cb, normal)
        
        new_vec_cb = vec_cb - 2 * proj_length * normal
        new_cb_pos = ca_pos + new_vec_cb
        
        modeller.positions[cb_atom.index] = openmm.Vec3(new_cb_pos[0], new_cb_pos[1], new_cb_pos[2]) * LENGTH

        atoms_to_keep_for_rebuild = {n_atom, ca_atom, c_atom, cb_atom, o_atom}

        atoms_in_residue_to_delete = []
        for atom in residue.atoms():
            if atom not in atoms_to_keep_for_rebuild:
                atoms_in_residue_to_delete.append(atom)
        
        modeller.delete(atoms_in_residue_to_delete)

    return changed

class ForceFieldMinimizerHeadTail(ForceFieldMinimizer):

    def _fix_cyclic(self, pdb_str, cyclic_chains=None, cyclic_opts=None):
        fixer = pdbfixer.PDBFixer(pdbfile=io.StringIO(pdb_str))
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=7.0)

        modeller = Modeller(fixer.topology, fixer.positions)
        connects = []

        if cyclic_chains:
            atoms_to_remove = []
            
            for chain in modeller.topology.chains():
                if chain.id not in cyclic_chains:
                    continue
                
                first_res = list(chain.residues())[0]
                for atom in first_res.atoms():
                    if atom.name in ['H2', 'H3']:
                        atoms_to_remove.append(atom)
                
                last_res = list(chain.residues())[-1]
                for atom in last_res.atoms():
                    if atom.name == 'OXT':
                        atoms_to_remove.append(atom)
                        
            modeller.delete(atoms_to_remove)

            fixer.topology = modeller.topology
            fixer.positions = modeller.positions

            out_handle = io.StringIO()
            openmm_app.PDBFile.writeFile(fixer.topology, fixer.positions, out_handle, keepIds=True)
            pdb_fixed = out_handle.getvalue()
            new_fixer = pdbfixer.PDBFixer(pdbfile=io.StringIO(pdb_fixed))

            n_term_atom, c_term_atom = None, None
            for chain in new_fixer.topology.chains():
                if chain.id not in cyclic_chains:
                    continue
                
                for i, res in enumerate(chain.residues()):
                    if i == 0:
                        for atom in res.atoms():
                            if atom.name == 'N': n_term_atom = atom
                    elif i == len(chain) - 1:
                        for atom in res.atoms():
                            if atom.name == 'C': c_term_atom = atom
                
                if n_term_atom and c_term_atom:
                    connects.append('CONECT' + str(c_term_atom.id).rjust(5) + str(n_term_atom.id).rjust(5))
                    connects.append('CONECT' + str(n_term_atom.id).rjust(5) + str(c_term_atom.id).rjust(5))
                else:
                    print("cannot find N or C atom in chain.")

            pdb_fixed = self._add_connects(pdb_fixed, connects)
        
        return pdb_fixed, connects

    def _check_and_correct_chirality(self, pdb_str, cyclic_chains):
        modeller = Modeller(openmm_app.PDBFile(io.StringIO(pdb_str)).topology,
                            openmm_app.PDBFile(io.StringIO(pdb_str)).positions)
        
        any_changed = False
        for chain in modeller.topology.chains():
            if chain.id not in cyclic_chains: continue
            residues_in_chain = list(chain.residues())
            residues_to_check = []
            if residues_in_chain:
                residues_to_check.append(residues_in_chain[0])
                if len(residues_in_chain) > 1:
                    residues_to_check.append(residues_in_chain[-1])

            for res in residues_to_check:
                if res.name not in ['GLY']:
                   changed = _ensure_l_chirality(modeller, res)
                   if changed:
                       any_changed = True
        
        out_handle_corrected = io.StringIO()
        openmm_app.PDBFile.writeFile(modeller.topology, modeller.positions, out_handle_corrected, keepIds=True)
        return out_handle_corrected.getvalue(), any_changed

    def __call__(self, pdb_path, out_path, return_info=True, cyclic_chains=None, cyclic_opts=None):
        try:
            pdb_str = None
            if '\n' not in pdb_path and pdb_path.lower().endswith(".pdb"):
                with open(pdb_path, 'r') as f:
                    pdb_str = f.read()
            else:
                raise ValueError("Input pdb_path is invalid.")
            
            pdb_fixed, connects = self._fix_cyclic(pdb_str, cyclic_chains, cyclic_opts)
            current_pdb_str, _, current_ret = self._minimize(pdb_fixed, out_path)
            
            if current_pdb_str is None or not current_pdb_str.strip():
                raise RuntimeError("Initial minimization failed.")

            for i in range(3): # Loop for three passes of chirality check and minimization
                pdb_corrected_str, changed = self._check_and_correct_chirality(current_pdb_str, cyclic_chains)
                
                if not changed and i > 0: # If no changes and not the very first pass, we can stop early
                    break
                
                # Re-fix the PDB after chirality correction and before next minimization
                re_fixer = pdbfixer.PDBFixer(pdbfile=io.StringIO(pdb_corrected_str))
                re_fixer.findMissingResidues()
                re_fixer.findMissingAtoms()
                re_fixer.addMissingAtoms()
                re_fixer.addMissingHydrogens(pH=7.0)
                
                out_handle_re_fixed = io.StringIO()
                openmm_app.PDBFile.writeFile(re_fixer.topology, re_fixer.positions, out_handle_re_fixed, keepIds=True)
                pdb_re_fixed_str = out_handle_re_fixed.getvalue()
                
                pdb_re_fixed_str_fixed, connects = self._fix_cyclic(pdb_re_fixed_str, cyclic_chains, cyclic_opts)

                # Perform minimization
                minimized_pdb_str, _, min_ret = self._minimize(pdb_re_fixed_str_fixed, out_path)

                if minimized_pdb_str is None or min_ret is None:
                    raise RuntimeError(f"Minimization pass {i+1} failed.")
                
                current_pdb_str = minimized_pdb_str
                current_ret = min_ret
            
            # Final processing after the loop
            pdb_min_final = self._add_connects(current_pdb_str, connects)
            pdb_min_final = self._add_energy_remarks(pdb_min_final, current_ret)
            
            if not os.path.exists(os.path.dirname(out_path)):
                os.makedirs(os.path.dirname(out_path))
            
            pdb_min_bytes = pdb_min_final.encode('utf-8')
            compressed_pdb = gzip.compress(pdb_min_bytes)
            
            with open(out_path, 'wb') as f:
                f.write(compressed_pdb)

            return (current_ret["einit"], current_ret["efinal"],
                current_ret["HarmonicBondForce_init"], current_ret["HarmonicBondForce_final"],
                current_ret["PeriodicTorsionForce_init"], current_ret["PeriodicTorsionForce_final"],
                current_ret["HarmonicAngleForce_init"], current_ret["HarmonicAngleForce_final"],
                current_ret["CMAPTorsionForce_init"], current_ret["CMAPTorsionForce_final"],
                current_ret["CMMotionRemover_init"], current_ret["CMMotionRemover_final"],
                current_ret["CustomTorsionForce_init"], current_ret["CustomTorsionForce_final"],
                current_ret["CustomNonbondedForce_init"], current_ret["CustomNonbondedForce_final"],
                current_ret["CustomBondForce_init"], current_ret["CustomBondForce_final"])
        
        except Exception as e:
            print(f"error: {e}")
            return None
