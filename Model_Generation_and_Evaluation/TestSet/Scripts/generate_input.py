import argparse
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Chain import Chain
import os
import warnings
from Bio import BiopythonWarning
from Bio.PDB.PDBParser import PDBConstructionWarning


warnings.filterwarnings("ignore", category=BiopythonWarning)
warnings.filterwarnings("ignore", category=PDBConstructionWarning)


def count_amino_acids_in_chain(chain):
    amino_acid_count = 0
    for residue in chain:
        if residue.get_resname() in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                                     'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']:
            amino_acid_count += 1
    return amino_acid_count


def process_pdb_files(input_dir, output_dir):
    parser = PDBParser()
    io = PDBIO()

    receptors_folder = os.path.join(output_dir, 'receptors')
    os.makedirs(receptors_folder, exist_ok=True)

    complexes_folder = os.path.join(output_dir, 'complexes')
    os.makedirs(complexes_folder, exist_ok=True)

    pdb_files = [f for f in os.listdir(input_dir) if f.lower().endswith('.pdb')]
    if not pdb_files:
        print("Warning: No PDB files found in input directory")
        return

    for pdb_file in pdb_files:
        pdb_path = os.path.join(input_dir, pdb_file)
        try:
            structure = parser.get_structure('pdb', pdb_path)
            model = structure[0]
            chains = list(model)
            pdb_name = os.path.basename(pdb_file)
            pdb_id = os.path.splitext(pdb_name)[0]

            l_chain = None
            for chain in chains:
                if chain.id == 'L':
                    l_chain = chain
                    break
            if not l_chain:
                print(f"Error: Chain L not found in {pdb_file}")
                continue

            other_chains = [chain for chain in chains if chain.id != 'L']
            if not other_chains:
                print(f"Error: No chains other than L found in {pdb_file}")
                continue

            r_chain = Chain('R')
            for chain in other_chains:
                for residue in chain:
                    r_chain.add(residue)

            receptor_name = f"{pdb_id}_receptor.pdb"
            receptor_path = os.path.join(receptors_folder, receptor_name)
            io.set_structure(r_chain)
            io.save(receptor_path)

            pdb_id_folder = os.path.join(complexes_folder, pdb_id)
            os.makedirs(pdb_id_folder, exist_ok=True)

            peptide_path = os.path.join(pdb_id_folder, 'peptide.pdb')
            io.set_structure(l_chain)
            io.save(peptide_path)

            protein_path = os.path.join(pdb_id_folder, 'pocket.pdb')
            io.set_structure(r_chain)
            io.save(protein_path)

        except Exception as e:
            print(f"Error processing {pdb_file}: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process PDB files to extract receptors and complexes.')
    parser.add_argument('-i', '--input', required=True, help='Input directory containing PDB files')
    parser.add_argument('-o', '--output', required=True, help='Output directory for results')
    args = parser.parse_args()

    process_pdb_files(args.input, args.output)