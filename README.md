😄 Hi, this is **CPSea**, a large-scale cyclic peptide-protein complex dataset derived from AlphaFold Database. We hope that you can always find the right part out of the CPSea to fulfill your requirements. 🌊

🗺️ This GitHub site mainly focuses on scripts and tutorials dataset generation and evaluation, as well as the evaluation of target-conditioned cyclic peptide design models that are trained on CPSea. Check out our [**Kaggle**](https://www.kaggle.com/datasets/ziyiyang180104/cpsea) and [**Zenodo**](https://zenodo.org/records/16417466) sites for clean data, index files and property files.

💡 In this repository, we provide scripts and introductions for dataset generation and evaluation, and also how to use re-trained models to generate cyclic peptides and to evaluate the model outputs.

[**Download Current Datasets**](#download-current-datasets)

[**Dataset Generation**](#dataset-generation)

[**Dataset Evaluation**](#dataset-evaluation)

* [Structure Validity](#structure-validity)
* [Binding Affinity](#binding-affinity)
* [Wet-lab Compatibility](#wet-lab-compatibility)
* [Designability, Diversity and Novelty](#designability-diversity-and-novelty)

[**Model Generation and Evaluation**](#model-generation-and-evaluation)

---

## Download Current Datasets

We store the main body of datasets in Kaggle and Zenodo. To download clean data of datasets, run:

```
curl -C - -o CPSea.zip https://zenodo.org/api/records/16417466/files-archive
```

In which are CPSea datasets derived from AFDB and PDB.


## Dataset Generation

🚧 We will release the soure codes for data generation after publishing our paper ~

## Dataset Evaluation

### Structure Validity

**For Ramachadran analysis**

```
python Dataset_Evaluation/rama_analysis.py -i <PDB_path_list> -o <Rama_output> -c <cores>
```

The input should be a list file containing paths to PDB in each line, and the same for the following scripts. This script will generate a CSV file that contains the proportion of each structure's Ramachandran angles falling into the Allowed and Favoured regions.

**For PLIP analysis**

```
python Dataset_Evaluation/plip_analysis.py -i <PDB_path_list> -o <PlIP_output> -c <cores>
```

This will generate a CSV file that stores the counts of various types of identified interactions in the interface of each structure. To obtain statistical information, run:

```
python Dataset_Evaluation/plip_handeler.py <PlIP_output>
```

This will generate a CSV file, whose name is the PLIP output file name with an additional "_stats" suffix, which stores statistical information on the aforementioned PLIP analysis results and records the proportions of various types of interactions.

### Binding Affinity

**For Rosetta ddG**

```
python Dataset_Evaluation/rosetta_analysis.py -i <PDB_path_list> -o <Rosetta_sc> -c <cores>
```

This will generate a JSON file that records Rosetta's output. We can extract the ddG information from it using the following script:

```
python Dataset_Evaluation/ddG_extractor.py -i <Rosetta_sc> -o <Rosetta_csv>
```

This will generate a CSV file that clearly records the Rosetta ddG for each structure.

**For Vina score**

**Note that** the [Vina](https://vina.scripps.edu) command line in `vina_analysis.py` needs to be changed to fit your installation.

```
python Dataset_Evaluation/vina_analysis.py -i <PDB_path_list> -o <Vina_out> -c <cores>
```

This will generate a CSV file that stores the Vina affinity scores for each structure, in the "score only" mode without re-docking and optimization.

### Wet-lab Compatibility

Metrics including GRAVY, logP, rTPSA can be calculated by running:

```
python Dataset_Evaluation/water_or_oil.py -i <PDB_path_list> -o <WoO_out> -c <cores> --chain_id <chain_id>
python Dataset_Evaluation/GRAVY_calculator.py -i <PDB_path_list> -o <GRAVY_out> -c <cores> --chain_id <chain_id>
```

All results for these metrics are in the output CSV files.

### Designability, Diversity and Novelty

We employed [**HighFold2**](https://github.com/hongliangduan/HighFold2) and [**Foldseek-multimer**](https://github.com/steineggerlab/foldseek) for these evaluation. 

**For Designability**

We refold the head-tail and disulfide cyclic peptides using HighFold2 without MSA input. Herein, we provide a script to generate HighFold2 inputs from PDB files.

```
python Dataset_Evaluation/generate_HF2_input.py -i <pdb_dir> -o <a3m_dir> (--sort)
```

The script read a directory of PDB files, and extract their **chain L** to make corresponding MSA files, saving to a new directory. The MSA directory can then be used as the input of HighFold2. Since the refolding of disulfide bonds requires specifying the numbers of the amino acids that form the disulfide bonds, you can use the --sort parameter to save sequences of different lengths in separate subfolders, which makes it more convenient to perform batch refolding.

We have attempted to use AF3 and Boltz for the refolding of cyclic isopeptides, but found it difficult to achieve. Here, we provide scripts for generating input files for AF3 and Boltz, which can produce input folders from PDB folders based on template files. You can modify the atom names in the templates to realize the prediction of different types of cyclic peptides.

```
python Dataset_Evaluation/generate_AF3_input.py -i <pdb_dir> -o <json_dir> -t AF3_template.json
python Dataset_Evaluation/generate_Boltz_input.py -i <pdb_dir> -o <yaml_dir> -t Boltz_template.yaml
```

For scRMSD calculation, we provide a script for automatical comparasion and calculation.

```
python Dataset_Evaluation/calculate_scRMSD.py <folder1> <folder2> -o <scRMSD_csv>
```

where folder1 and folder2 are PDB directories containing CPSea structures and refolded structures.

**For Diversity and Novelty**

We use foldseek easy-multimercluster and easy-multimersearch to calculate diversity and novelty of our dataset.

```
foldseek easy-multimercluster <dataset_dir> <output_dir>/clu <temp_dir> --alignment-type 2 --cov-mode 0 --min-seq-id 0 --multimer-tm-threshold 0.65 --chain-tm-threshold 0.5 --interface-lddt-threshold 0.65 --threads <threads>
foldseek easy-multimersearch  <dataset_dir> <path_to_pdb> <output_dir>/out --alignment-type 2 --tmscore-threshold 0.0 --max-seqs 1000 --format-output query,target,complexqtmscore,complexttmscore,lddt --threads <threads>
```

Below we provide two codes for calculating the number of clusters and novelty.

```
python Dataset_Evaluation/check_unique.py -i <cluster_tsv> -o <unique_list> 
python Dataset_Evaluation/check_qtm_max_average.py <out_file> (--ignore_R)
```

`check_unique.py` generates a list file containing one file from each cluster, so that the number of lines equals to the number of clusters, and the list file can be used as the index of a non-redundant subset.

`check_qtm_max_average.py` outputs the average qTm max in stdout. You can add the --ignore_R argument to exclude receptors in this analysis.

## Model Generation and Evaluation

Please check the [introduction file](Model_Generation_and_Evaluation/README.md) for details.

Related weights for three re-trained models can be downloaded from [**Kaggle**](https://www.kaggle.com/datasets/ziyiyang180104/cpsea).

## Contact

📞 This repository is still under development. If there are any questions about any protocols, please do not hesitate to contact yangzy23@mails.tsinghua.edu.cn. Thank you for your interest in our CPsea!


