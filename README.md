😄 Hi, this is **CPSea**, a large-scale cyclic peptide-protein complex dataset derived from AlphaFold Database. We hope that you can always find the right part out of the CPSea to fulfill your requirements. 🌊

🗺️ This GitHub site mainly focuses on scripts and tutorials dataset generation and evaluation, as well as the evaluation of target-conditioned cyclic peptide design models that are trained on CPSea. Check out our Kaggle and Zenodo sites for clean data and relevant files:

Kaggle: [https://www.kaggle.com/datasets/ziyiyang180104/cpsea](https://www.kaggle.com/datasets/ziyiyang180104/cpsea)

Zenodo: [https://zenodo.org/records/16417466](https://zenodo.org/records/16417466)

---

## Download Current Datasets

We store the main body of datasets in Kaggle and Zenodo. To download clean data of datasets, run:

`curl -C - -o CPSea.zip https://zenodo.org/api/records/16417466/files-archive`

In which are CPSea datasets derived from AFDB and PDB.


## Dataset Generation

🚧 We will release the soure codes for data generation after publishing our paper ~

## Dataset Evaluation and Subset Curation 

### Structure Validity

**For Ramachadran analysis** on main-chain torsion angles, run: 

`python Dataset_Evaluation/rama_analysis.py -i <PDB_path_list> -o <Rama_output> -c <cores>`

The input should be a list file containing paths to PDB in each line, and the same for the following scripts. This script will generate a CSV file that contains the proportion of each structure's Ramachandran angles falling into the Allowed and Favoured regions.

**For PLIP analysis** on interface interaction distribution, run:

`python Dataset_Evaluation/plip_analysis.py -i <PDB_path_list> -o <PlIP_output> -c <cores>`

This will generate a CSV file that stores the counts of various types of identified interactions in the interface of each structure. To obtain statistical information, run:

`python plip_handeler.py <PlIP_output>`

This will generate a CSV file, whose name is the PLIP output file name with an additional "_stats" suffix, which stores statistical information on the aforementioned PLIP analysis results and records the proportions of various types of interactions.

### Binding Affinity

**For Rosetta ddG**, run:

'python rosetta_analysis.py -i <PDB_path_list> -o <Rosetta_sc> -c <cores>'

This will generate a JSON file that records Rosetta's output. We can extract the ddG information from it using the following script:

`python ddG_extractor.py -i <Rosetta_sc> -o <Rosetta_csv>`

This will generate a CSV file that clearly records the Rosetta ddG for each structure.

**For Vina score**, run:

`python vina_analysis.py -i <PDB_path_list> -o <Vina_out> -c <cores>`

This will generate a CSV file that stores the Vina affinity scores for each structure, in the "score only" mode without re-docking and optimization.

### Wet-lab Compatibility

Metrics including GRAVY, logP, rTPSA can be calculated by running:

`python water_or_oil.py -i <PDB_path_list> -o <WoO_out> -c <cores> --chain_id <chain_id>`

`python GRAVY_calculator.py -i <PDB_path_list> -o <GRAVY_out> -c <cores> --chain_id <chain_id>`

All results for these metrics are in the output CSV files.

### Designability, Diversity and Novelty

We employed [**HighFold2**](https://github.com/hongliangduan/HighFold2) and [**Foldseek-multimer**](https://github.com/steineggerlab/foldseek) for these evaluation. 

**Designability**

We refold the head-tail and disulfide cyclic peptides using HighFold2 without MSA input. Herein, we provide a script to generate HighFold2 inputs from PDB files.

`python generate_HF2_input.py -i <pdb_dir> -o <a3m_dir> (--sort)`

The script read a directory of PDB files, and extract their **chain L** to make corresponding MSA files, saving to a new directory. The MSA directory can then be used as the input of HighFold2. Since the refolding of disulfide bonds requires specifying the numbers of the amino acids that form the disulfide bonds, you can use the --sort parameter to save sequences of different lengths in separate subfolders, which makes it more convenient to perform batch refolding.

We have attempted to use AF3 and Boltz for the refolding of cyclic isopeptides, but found it difficult to achieve. Here, we provide scripts for generating input files for AF3 and Boltz, which can produce input folders from PDB folders based on template files. You can modify the atom names in the templates to realize the prediction of different types of cyclic peptides.

`python generate_AF3_input.py -i <pdb_dir> -o <json_dir> -t AF3_template.json`

`python generate_Boltz_input.py -i <pdb_dir> -o <yaml_dir> -t Boltz_template.json`

For scRMSD calculation, we provide a script for automatical comparasion and calculation.
