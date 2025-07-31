⏰ This repository is still under building. We are trying our best to set up everything, and will finish this in 24 hours (by August 1st 20:00, China standard time)

---

😊 hi

This is an introduction for replicating the results of model retraining and evaluation in the manuscript of:
CPSea: Large-scale cyclic peptide-protein complex dataset for machine learning in cyclic peptide design.

## Step 0: Model setup

We re-trained three target-conditioned peptide design models in our work, namely [DiffPepBuilder](https://github.com/YuzheWangPKU/DiffPepBuilder), [PepFlow](https://github.com/Maikuraky/PepFlow), and [PepGLAD](https://github.com/THUNLP-MT/PepGLAD). First, please install the three models following their guidance. Then, change these files:

1. move `ModifyModels/DiffPepBuilder/run_inference_new.py` into the DiffPepBuilder directory
2. 

## Step 1: Data preparation

We have pre-processed the input complexes from LNR and CPSet, renaming the binder chain as L chain. After that,
initial complexes are transformed to formats that are compatible with the three models by `generate_input.py`

```
python generate_input.py -i <original_pdb_dir> -o <output_dir>
```

The scripts will generate two folders in <output_dir>. The `complexes` folder stores initial complexes separately in subfolders, with receptors named `pocket.pdb` and binders named `peptide.pdb`. The `receptors` folder stores receptors in initial complexes. This is the format of the PepMerge dataset in the PepFlow paper.

## Step 2: Inference

### 2.1 DiffPepBuilder

**First, define epitopes**

```
python generate_epitopes.py -i <complex_dir> -o epitopes.json
```

This generates a json file containing epitope information for each receptor in the format required by DiffPepBuilder

**Then, process data**

```
python {path_to_diffpepbuilder}/experiments/process_receptor.py --pdb_dir {path_to_data}/receptors --write_dir {path_to_diffpepbuilder}/data/receptor_data --receptor_info_path {path_to_data}/epitopes.json
```

This will generate .pkl files for each target.

**Next**, we can do inference to generate *pseudo* cyclic peptides (Because we have not modified the model to generate cyclic peptides directly, the outputs of these re-trained models are actually linear peptides with their terminals close to each other)

In our manuscript, we let models generate peptides with the same length as in initial complexes. To do this more easily, we modified the run_inference.py by adding a hard-coded path to a csv file that contains desired length for each target, which is in line 288 in `run_inference_new.py`.

After adding the csv file path, you can run:

```
export BASE_PATH={path_to_diffpepbuilder}
torchrun --nproc-per-node=1 experiments/run_inference_new.py data.val_csv_path=data/receptor_data/metadata_test.csv
```

### 2.2 PepFlow

PepFlow seems don't cut out epitopes of receptors during generation, making it easy to occur OOM error. To avoid this, 
first run:

```
python Flow_pocket_dealer.py <old_complex_dir> <new_complex_dir>
```

This will generate a new folder in which the pocket.pdb contains only the epitope.

**To define epitopes**, first generate a list file consists of subfolder names in {new_complex_dir} folder. Then, modify 
`models_con/pep_dataloader.py` as follow: 

1. In line 37, change the list path to the subfolder name list.
2. In line 213, type in the {new_complex_dir} folder path in structure_dir, and the place you want to save the dataset in dataset_dir.
   Also, change the name to whatever you like.
3. The script will not make the directory itself, so we need to make the dataset_dir ourselves.
   
Then, run:

```
python models_con/pep_dataloader.py
```

**Now, we can do inference**. First, modify `configs/learn_angle.yaml`, change the structure_dir, dataset_dir and name of the 
generate dataset. Second, make an output diectory. Then run:

```
python models_con/inference.py --config configs/learn_angle.yaml  --device cuda:0 --ckpt path/to/retrained/ckpt --output output/path \
--num_samples 100
```

**Finally, we can do sampling**.

```
python models_con/sample.py --SAMPLEDIR output/path/of/inference
```

### 2.3 PepGlad

**To define epitopes**, we simply use the json file generated for DiffPepBuilder, but in a different format:

```
python transform_json.py -i {epitope_json_for_DiffPepBuilder} -o {dir_for_epitope_json_of_PepGlad}
```

Now that we can generate by running `batch_generate.py`, modify following content:

1. In line 7, type in the path to the length csv
2. Update input_dir, pocket_dir, and output_dir. Input_dir is the receptor directory generated in Step 1, pocket_dir is the json
directory above.

Then, run:

```
python batch_generate.py
```

## Step 3: Post-processing

As mentioned before, we have not modified these models to make them generate cyclic peptides directly. Therefore, some post-processing
is needed to transform the generated "pseudo" cyclic peptides into real ones.

### 3.1 Rename

To make the rest process smooth, we first rename the outputs from three models, to name the receptor as R and binder as L, and also unify
the file name as {targetID}_{num}.pdb. We suggest to collect all generated pdb files in "XXX_generated" folder first (remain subfolders),
then run:

```
python rename_glad.py -i glad_generated
python rename_flow.py -i flow_generated
python rename_diff.py -i diff_generated
```

These would generate the renamed folders.

### 3.2 Cyclization and Relax

Then, we apply a similar cyclization and relax protocol as that in the CPSea pipeline. First, we need to cut the pocket off for PepGlad 
generated complexes, because PepGlad will fulfill the receptor automatically.

```
python glad_cut_pocket.py --renamed glad_renamed --receptor {receptors} --json CPSet_epitope_for_glad
```

Then, we check the CB and length of each generated peptide, to exclude peptides whose CB distance is not within [3,8] range:

```
python cb_distance_calculater.py glad_cutpocket glad_initial_CB.csv
python filter_CB.py glad_initial_CB.csv glad_cutpocket cyc_failed

python cb_distance_calculater.py flow_renamed flow_initial_CB.csv
python filter_CB.py flow_initial_CB.csv flow_renamed cyc_failed

python cb_distance_calculater.py diff_renamed diff_initial_CB.csv
python filter_CB.py diff_initial_CB.csv diff_renamed cyc_failed
```

Now, before cyclization and relax, it seems better to replace pockets by the same residues in oringinal receptors.

```
python combine_epitope.py --generated diff_renamed --receptors CPSet/clean_receptors --epitopes CPSet/CPSet_epitopes.json \
--output diff_reconstructed
python combine_epitope.py --generated flow_renamed --receptors flow_gt --epitopes CPSet/CPSet_epitopes.json \
--output flow_reconstructed
python combine_epitope.py --generated glad_cutpocket --receptors CPSet/clean_receptors --epitopes CPSet/CPSet_epitopes.json \
--output glad_reconstructed
```

Now we can calculate the success rate for cyclization. Then, we calculate both CB distance and length, as inputs for cyclization and relax:
(needs change path prefix)

```
python cb_distance_and_length.py CPSet_Diff/diff_reconstructed CPSet_Diff/CPSet_Diff_CB_length.csv
python cb_distance_and_length.py CPSet_Flow/flow_reconstructed CPSet_Flow/CPSet_Flow_CB_length.csv
python cb_distance_and_length.py CPSet_Glad/glad_reconstructed CPSet_Glad/CPSet_Glad_CB_length.csv
```

Then run `relax_mp_model.py`, change the file paths and configs in the script.

After running this, we re-combine the cyclic peptides and their full-length receptors by running:

```
python combine_receptor.py --relaxed diff_relaxed/ --receptors ../../evaluate_dataset/CPSet/clean_receptors  --output diff_full
python combine_receptor.py --relaxed flow_relaxed/ --receptors ../../evaluate_dataset/CPSet/clean_receptors  --output flow_full
python combine_receptor.py --relaxed glad_relaxed/ --receptors ../../evaluate_dataset/CPSet/clean_receptors  --output glad_full
```

Now that we finished the target-conditioned cyclic peptide generation!

## Step 4: Evaluation

### 4.1 Energy filter

First, we use rosetta ddG to check if the binding interfaces are designed reasonable

```
find {diff_full} -name "*.pdb" > diff_full.list
python rosetta_analysis.py -i diff_full.list -o rosetta_diff_full.sc -p 64
python ddG_extractor.py -i rosetta_diff_full.sc -o rosetta_ddG_diff_full.csv
```

Now, we can calculate the energy success rate by checking Rosetta ddG. To exclude energy failed designs, run:

```
python filter_ddG.py rosetta_ddG_diff_full.csv {diff_full} energy_failed
```

<diff_full> now contains the final structures of generated cyclic peptides.

### 4.2 Affinity

First, calculate Rosetta ddG and Vina score 

```
find {diff_full} -name "*.pdb" > Diff_Final.list
python rosetta_analysis_mp.py -i Diff_Final.list -o Rosetta_Diff_Final.sc -p 64
python ddG_extractor.py -i Rosetta_Diff_Final.sc -o Rosetta_ddG_Diff_Final.csv
python vina_analysis.py --input Diff_Final.list --output Vina_Diff.csv --jobs 64
```

Then, select the best for each target

```
python select_best.py -i Rosetta_ddG_Diff_Final.csv -o bestddG_Diff.csv
python select_best.py -i Vina_Diff.csv -o bestVina_Diff.csv
```

### 4.3 Structure Validity

We use [**PLIP**](https://plip-tool.biotec.tu-dresden.de/plip-web/plip/index) to analyze interface interactions and summarize the proportions of each types of interaction.

```
python plip_analysis_mp.py --input Diff_Final.list --output PLIP_Diff.csv --jobs 48
python plip_handeler.py PLIP_Diff.csv
```

The statistic review should be print in stdout.

We use Ramachadran plot to evaluate the validity of cyclic peptides themselves.

```
python rama_analysis_mp.py Diff_Final.list -o Rama_Diff.csv -c 48
```

### 4.4 Wet-lab Compatibility

We calculate GRAVY, logP, and rTPSA to assess the synthesis feasiblity and aggregation propensity of cyclic peptides.

```
python water_or_oil_mp.py --input_list Diff_Final.list --output_csv WoO_Diff.csv --chain_id L --n_cpu 4
python GRAVY_calculator.py -i Diff_Final.list -o GRAVY_Diff.csv -c L -p 4
```

### 4.5 Diversity, Novelty and Self-Consistency

For diversity, we cluster the generated structures by cyclic-aware rmsd and Foldseek

```
python mp_rmsd_cyclic.py --rmsd Diff_rmsd.npy --hist Diff_rmsd.png --cluster Diff_rmsd.csv --name_list Diff_Final.list --num_cores 48

foldseek easy-multimercluster {diff_full}  Diff_foldseek_cluster/clu  Diff_foldseek_cluster --multimer-tm-threshold 0.65 \
--chain-tm-threshold 0.5 --interface-lddt-threshold 0.65 --alignment-type 2 --cov-mode 0 --min-seq-id 0 --threads 32
```

For Novelty, we calculate qTm value against PDB database

```
foldseek easy-multimersearch {diff_full} {path/to/PDB/database} Diff_Final_foldseek_novel/out Diff_Final_foldseek_novel --alignment-type 2 \
--tmscore-threshold 0.0 --max-seqs 1000 --format-output query,target,complexqtmscore,complexttmscore,lddt --threads 48
```

Then, analyze the qTm by

```
python check_qtm_max_average.py out --ignore_R --ignore_same_id
```

For self-consistency, we use [**HighFold2**](https://github.com/hongliangduan/HighFold2) to refold head-tail and disulfide cyclic peptides, and calculate the Ca RMSD.


