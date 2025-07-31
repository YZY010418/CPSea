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

os.environ['OPENMM_CPU_THREADS'] = '1'

custom_xml = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    'custom', 'residue.xml'
))

def check_GLY_CB(topology):
    for chain in topology.chains():
        for res in chain.residues():
            if res.name == 'GLY':
                for atom in list(res.atoms()):
                    if atom.name == 'CB':
                        res._atoms.remove(atom)
    return topology

class ForceFieldMinimizer(object):

    def __init__(self, tolerance=2.39*unit.kilocalories_per_mole, platform='CPU'):
        super().__init__()
        converted_tolerance = tolerance.in_units_of(unit.kilojoules_per_mole)
        self.tolerance = converted_tolerance / unit.nanometer
        assert platform in ('CUDA', 'CPU')
        self.platform = platform

    def _fix(self, pdb_str, cyclic_chains, cyclic_opts):

        fixer = pdbfixer.PDBFixer(pdbfile=io.StringIO(pdb_str))
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()

        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens()

        fixer.topology = check_GLY_CB(fixer.topology)

        if cyclic_chains is not None:
            pdb_fixed, connects = self._fix_cyclic(fixer, cyclic_chains, cyclic_opts)
        else:
            out_handle = io.StringIO()
            openmm_app.PDBFile.writeFile(fixer.topology, fixer.positions, out_handle, keepIds=True)
            pdb_fixed = out_handle.getvalue()
            connects = []

        return pdb_fixed, connects
    
    def _fix_cyclic(self, fixer, cyclic_chains, cyclic_opts):

        out_handle = io.StringIO()
        openmm_app.PDBFile.writeFile(fixer.topology, fixer.positions, out_handle, keepIds=True)
        pdb_fixed = out_handle.getvalue()
        connects = []

        return pdb_fixed, connects

    def _get_pdb_string(self, topology, positions):
        with io.StringIO() as f:
            openmm_app.PDBFile.writeFile(topology, positions, f, keepIds=True)
            return f.getvalue()
        
    def _minimize(self, pdb_str, out_path):
        try:
            pdb = openmm_app.PDBFile(io.StringIO(pdb_str))
            force_field = openmm_app.ForceField("charmm36.xml", custom_xml) # referring to http://docs.openmm.org/latest/userguide/application/02_running_sims.html

            # The speed becomes really fast after cancelling H-bonds constraints
            system = force_field.createSystem(pdb.topology, constraints=None)

            # Freeze receptor chain by setting their atom mass to zero, this speed up the calculation a lot
            r_chain_atom_indices = []
            for chain in pdb.topology.chains():
                if chain.id == 'L': continue
                for atom in chain.atoms():
                    system.setParticleMass(atom.index, 0)
                    r_chain_atom_indices.append(atom.index)

            # Disable bond-related forces inside chain R
            for i, force in enumerate(system.getForces()):
                force.setForceGroup(i)
                if isinstance(force, openmm.HarmonicBondForce):
                    for index in range(force.getNumBonds()):
                        p1, p2, length, _ = force.getBondParameters(index)
                        if p1 in r_chain_atom_indices and p2 in r_chain_atom_indices:
                            force.setBondParameters(index, p1, p2, length, 0)
                elif isinstance(force, openmm.HarmonicAngleForce):
                    for index in range(force.getNumAngles()):
                        p1, p2, p3, angle, _ = force.getAngleParameters(index)
                        if p1 in r_chain_atom_indices and p2 in r_chain_atom_indices and p3 in r_chain_atom_indices:
                            force.setAngleParameters(index, p1, p2, p3, angle, 0)
                elif isinstance(force, openmm.PeriodicTorsionForce):
                    for index in range(force.getNumTorsions()):
                        p1, p2, p3, p4, periodicity, phase, _ = force.getTorsionParameters(index)
                        if p1 in r_chain_atom_indices and p2 in r_chain_atom_indices and p3 in r_chain_atom_indices and p4 in r_chain_atom_indices:
                            force.setTorsionParameters(index, p1, p2, p3, p4, periodicity, phase, 0)

            # Set up the simulation
            platform = openmm.Platform.getPlatformByName("CPU")
            simulation = openmm_app.Simulation(pdb.topology, system, openmm.LangevinIntegrator(0, 0.01, 0.0), platform)
            simulation.context.setPositions(pdb.positions)

            # Perform minimization
            state = simulation.context.getState(getEnergy=True, getPositions=True)
            ret = {}
            ret["einit"] = state.getPotentialEnergy().value_in_unit(ENERGY)
            ret["posinit"] = state.getPositions(asNumpy=True).value_in_unit(LENGTH)
            ret["init_pdb"] = self._get_pdb_string(simulation.topology, state.getPositions())


            for i, force in enumerate(system.getForces()):
                force_name = force.__class__.__name__
                state = simulation.context.getState(getEnergy=True, groups=1 << i)
                energy = state.getPotentialEnergy()
                ret[f"{force_name}_init"] = energy.value_in_unit(unit.kilocalories_per_mole)

            #PERTURBATION_SCALE = 0.05
            min_energy_index = 0
            for i in range(3):
                simulation.minimizeEnergy(tolerance=self.tolerance)

                state = simulation.context.getState(getEnergy=True, getPositions=True)
                ret[f"efinal_{i}"] = state.getPotentialEnergy().value_in_unit(ENERGY)
                ret[f"pos_{i}"] = state.getPositions(asNumpy=True).value_in_unit(LENGTH)

                energy = state.getPotentialEnergy().value_in_unit(ENERGY)
                if energy < ret[f"efinal_{min_energy_index}"]:
                    min_energy_index = i

                if ret[f"efinal_{i}"] < 0.0:
                    min_energy_index = i
                    break

                '''perturbation = np.zeros_like(ret[f"pos_{i}"])
                non_r_chain_indices = [j for j in range(len(ret[f"pos_{i}"])) if j not in r_chain_atom_indices]
                perturbation[non_r_chain_indices] = np.random.normal(0, PERTURBATION_SCALE, size=(len(non_r_chain_indices), 3))
                new_positions = ret[f"pos_{i}"] + perturbation
                new_positions_with_units = new_positions * unit.angstroms
                simulation.context.setPositions(new_positions_with_units)'''

            state = simulation.context.getState(getEnergy=True, getPositions=True)
            ret["efinal"] = ret[f"efinal_{min_energy_index}"]
            ret["pos"] = ret[f"pos_{min_energy_index}"]
            ret["min_pdb"] = self._get_pdb_string(simulation.topology, 
                                                 ret[f"pos_{min_energy_index}"] * LENGTH)

            for i, force in enumerate(system.getForces()):
                force_name = force.__class__.__name__
                state = simulation.context.getState(getEnergy=True, groups=1 << i)
                energy = state.getPotentialEnergy()
                ret[f"{force_name}_final"] = energy.value_in_unit(unit.kilocalories_per_mole)

            return ret['min_pdb'], ret['init_pdb'], ret
        except Exception as e:
            print(f"error in minimization: {e}")

    def _add_energy_remarks(self, pdb_str, ret):
        pdb_lines = pdb_str.splitlines()
        pdb_lines.insert(1, "REMARK   1  FINAL ENERGY:   {:.3f} KCAL/MOL".format(ret['efinal']))
        pdb_lines.insert(1, "REMARK   1  INITIAL ENERGY: {:.3f} KCAL/MOL".format(ret['einit']))
        return "\n".join(pdb_lines)

    def _add_connects(self, pdb_str, connects):
        exist_connects = [l for l in pdb_str.split('\n') if 'CONECT' in l]
        connects = [c for c in connects if c not in exist_connects]
        # add CONECT records at the end
        pdb_str = pdb_str.strip().strip('END').strip()
        pdb_str = pdb_str.split('\n')
        pdb_str = pdb_str + connects + ['END\n']
        pdb_str = '\n'.join(pdb_str)
        return pdb_str

    def __call__(self, pdb_path, out_path, return_info=True, cyclic_chains=None, cyclic_opts=None):
        try:
            with open(pdb_path, 'r') as f:
                pdb_str = f.read()
            pdb_fixed, connects = self._fix(pdb_str, cyclic_chains, cyclic_opts)
            pdb_min, pdb_init, ret = self._minimize(pdb_fixed, out_path)
            pdb_min = self._add_connects(pdb_min, connects)
            pdb_min = self._add_energy_remarks(pdb_min, ret)
            if not os.path.exists(os.path.dirname(out_path)):
                os.makedirs(os.path.dirname(out_path))
            with open(out_path, 'w') as f:
                f.write(pdb_min)

            return (ret["einit"], ret["efinal"], \
                ret["HarmonicBondForce_init"], ret["HarmonicBondForce_final"], \
                ret["PeriodicTorsionForce_init"], ret["PeriodicTorsionForce_final"], \
                ret["HarmonicAngleForce_init"], ret["HarmonicAngleForce_final"], \
                ret["CMAPTorsionForce_init"], ret["CMAPTorsionForce_final"], \
                ret["CMMotionRemover_init"], ret["CMMotionRemover_final"], \
                ret["CustomTorsionForce_init"], ret["CustomTorsionForce_final"], \
                ret["CustomNonbondedForce_init"], ret["CustomNonbondedForce_final"], \
                ret["CustomBondForce_init"], ret["CustomBondForce_final"])
        except Exception as e:
            print(f"error: {e}")


if __name__ == '__main__':
    import sys
    force_field = ForceFieldMinimizer()
    force_field(sys.argv[1], sys.argv[2])