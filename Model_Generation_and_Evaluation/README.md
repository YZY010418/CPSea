😊 hi

This is an introduction for replicating the results of model retraining and evaluation in the manuscript of:
CPSea: Large-scale cyclic peptide-protein complex dataset for machine learning in cyclic peptide design.

Due to time limit, the protocol described in this introduction has not been fully tested yet, please contact us if there are any problems.

Since the primary focus of CPSea is to construct a cyclic peptid -protein complex dataset to mitigate data scarcity, and developing dedicated cyclic peptide design models is not the core objective of this work, we have not modified existing models to enable direct cyclic peptide generation. We apologize for any inconvenience this may cause. A straightforward yet powerful cyclic peptide design model is currently under development and will be released in the near future.

[**Step 0: Model setup**](#step-0-model-setup)

[**Step 1: Data preparation**](#step-1-data-preparation)

[**Step 2: Inference**](#step-2-inference)

[**Step 3: Post-processing**](#step-3-post-processing)

[**Step 4: Evaluation**](#step-4-evaluation)

## Step 0: Model setup

We re-trained three target-conditioned peptide design models in our work, namely [DiffPepBuilder](https://github.com/YuzheWangPKU/DiffPepBuilder), [PepFlow](https://github.com/Maikuraky/PepFlow), and [PepGLAD](https://github.com/THUNLP-MT/PepGLAD). First, please install the three models following their guidance. Then, change these files:

1. move `ModifyModels/run_inference_new.py` into the {DiffPepBuilder}/experiments directory
2. move `ModifyModels/batch_generate.py` into the {PepGlad} directory
3. move `ModifyModels/run.py` into the {PepGlad}/api directory

Now that we are ready for the following steps.

## Step 1: Data preparation

We have pre-processed the input complexes from LNR and CPSet into `complexes` and `receptors`.

* The `complexes` folder stores complexes separately in subfolders, with receptors named `pocket.pdb` and binders named `peptide.pdb`.
* The `receptors` folder stores receptors in initial complexes.

This is the format of the PepMerge dataset in the PepFlow paper. The pre-processed data are stored in `TestSet/LNR` and `TestSet/CPSet`.

On top of this, we need to define epitope for each target, which can be achieved as follows:

1. For DiffPepBuilder, run

```
python TestSet/Scripts/generate_epitopes.py -i <complex_dir> -o epitopes.json
```

This generates a json file containing epitope information for each receptor in the format required by DiffPepBuilder

2. For PepFlow, because the model could not extract epitopes, we need to make epitope structures mannually, which can be done by:

```
python TestSet/Scripts/Flow_pocket_dealer.py <old_complex_dir> <new_complex_dir>
```

This will generate a new folder in which the pocket.pdb contains only the epitope.

3. For PepGlad, what we need to do is simply change the format of the `epitopes.json` file generated above. 

```
python TestSet/Scripts/transform_json.py -i <epitopes.json> -o <epitopes_glad.json>
```

## Step 2: Inference

During inference, don't forget to change checkpoints into the re-trained versions, as provided in [Kaggle](https://www.kaggle.com/datasets/ziyiyang180104/cpsea).

### 2.1 DiffPepBuilder

First, process the receptor data

```
python {path_to_diffpepbuilder}/experiments/process_receptor.py --pdb_dir {path_to_data}/receptors --write_dir {path_to_diffpepbuilder}/data/receptor_data --receptor_info_path {path_to_data}/epitopes.json
```

This will generate .pkl files for each target.

Next, we can do inference to generate *pseudo* cyclic peptides (Because we have not modified the model to generate cyclic peptides directly, the outputs of these re-trained models are actually linear peptides with their terminals close to each other)

In our manuscript, we let models generate peptides with the same length as in initial complexes. To do this more easily, we modified the run_inference.py by adding a csv file input, named `run_inference_new.py`.

The csv file should have two columns, the first column with the title "pdb_names", and the second column with the title "length".

```
export BASE_PATH={path_to_diffpepbuilder}
torchrun --nproc-per-node=1 experiments/run_inference_new.py data.val_csv_path=data/receptor_data/metadata_test.csv --csv_path <the_length_csv>
```

### 2.2 PepFlow

First, process the receptor data. A list file consists of subfolder names in {new_complex_dir} folder is needed. Then, modify 
`models_con/pep_dataloader.py` as follow: 

1. In line 37, change the list path to the subfolder name list.
2. In line 213, type in the {new_complex_dir} folder path in structure_dir, and the place you want to save the dataset in dataset_dir.
   Also, change the name to whatever you like.
3. NOTE: The script **will not** make the directory itself, so we need to make the dataset_dir ourselves.
   
Then, run:

```
python {path_to_pepflow}/models_con/pep_dataloader.py
```

Now, we can do inference. First, modify `configs/learn_angle.yaml`, change the structure_dir, dataset_dir and name of the 
generate dataset. Second, make an output diectory (since the script does not create the directory itself). Then run:

```
python {path_to_pepflow}/models_con/inference.py --config configs/learn_angle.yaml  --device cuda:0 --ckpt path/to/retrained/ckpt --output output/path \
--num_samples 100
```

Finally, we can do sampling.

```
python {path_to_pepflow}/models_con/sample.py --SAMPLEDIR output/path/of/inference
```

### 2.3 PepGlad

It is very simple to run PepGlad. Similar to DiffPepBuilder, we still need a csv file containing "pdb_names" and "length" columns.

```
python {path_to_pepglad}/batch_generate.py --length_csv <length_csv> --input_dir <receptor_PDB_dir> --pocket_dir <transformed_json> --output_dir <out_dir>
```

## Step 3: Post-processing

As mentioned before, we have not modified these models to make them generate cyclic peptides directly. Therefore, some post-processing
is needed to transform the generated "pseudo" cyclic peptides into real ones.

### 3.1 Rename

To make the rest process smooth, we first rename the outputs from three models, to name the receptor as R and binder as L, and also unify
the file name as {targetID}_{num}.pdb. We suggest to collect all generated pdb files in "XXX_generated" folder first (remain subfolders),
then run:

```
python PostProcess/rename_glad.py -i glad_generated
python PostProcess/rename_flow.py -i flow_generated
python PostProcess/rename_diff.py -i diff_generated
```

These would generate the renamed folders. 

### 3.2 Cyclization and Relax

Then, we apply a similar cyclization and relax protocol as that in the CPSea pipeline. First, we need to cut the pocket off for PepGlad 
generated complexes, because PepGlad will fulfill the receptor automatically.

```
python PostProcess/glad_cut_pocket.py --renamed glad_renamed --receptor {receptors} --json CPSet_epitope_for_glad
```

Then, we check the CB and length of each generated peptide, to exclude peptides whose CB distance is not within [3,8] range:

```
python PostProcess/cb_distance_calculater.py glad_cutpocket glad_initial_CB.csv
python PostProcess/filter_CB.py glad_initial_CB.csv glad_cutpocket cyc_failed

python PostProcess/cb_distance_calculater.py flow_renamed flow_initial_CB.csv
python PostProcess/filter_CB.py flow_initial_CB.csv flow_renamed cyc_failed

python PostProcess/cb_distance_calculater.py diff_renamed diff_initial_CB.csv
python PostProcess/filter_CB.py diff_initial_CB.csv diff_renamed cyc_failed
```

Now, before cyclization and relax, it seems better to replace pockets by the same residues in oringinal receptors.

```
python PostProcess/combine_epitope.py --generated diff_renamed --receptors CPSet/clean_receptors --epitopes CPSet/CPSet_epitopes.json \
--output diff_reconstructed
python PostProcess/combine_epitope.py --generated flow_renamed --receptors flow_gt --epitopes CPSet/CPSet_epitopes.json \
--output flow_reconstructed
python PostProcess/combine_epitope.py --generated glad_cutpocket --receptors CPSet/clean_receptors --epitopes CPSet/CPSet_epitopes.json \
--output glad_reconstructed
```

Now we can calculate the success rate for cyclization. Then, we calculate both CB distance and length, as inputs for cyclization and relax. Please use absolute path here, because the relax script relies on this information to find the file.

```
python PostProcess/cb_distance_and_length.py {abs_path}/CPSet_Diff/diff_reconstructed CPSet_Diff/CPSet_Diff_CB_length.csv
python PostProcess/cb_distance_and_length.py {abs_path}/CPSet_Flow/flow_reconstructed CPSet_Flow/CPSet_Flow_CB_length.csv
python PostProcess/cb_distance_and_length.py {abs_path}/CPSet_Glad/glad_reconstructed CPSet_Glad/CPSet_Glad_CB_length.csv
```

Then run `relax_mp_model.py`, change the file paths and configs in the script as follows:

1. in line 86, set the cores to use (required)
2. in line 87, set the path to save relaxed structures and relevant files (required)
3. in line 88, set the path of <CB_length.csv> to read in (required)
4. in line 89, set the name of the output metadata file (optional)
5. in line 90, set the number of tasks for tqdm bar (optional)

```
python PostProcess/Relax/relax_mp_model.py
```

After running this, we re-combine the cyclic peptides and their full-length receptors by running:

```
python PostProcess/combine_receptor.py --relaxed diff_relaxed/ --receptors ../../evaluate_dataset/CPSet/clean_receptors  --output diff_full
python PostProcess/combine_receptor.py --relaxed flow_relaxed/ --receptors ../../evaluate_dataset/CPSet/clean_receptors  --output flow_full
python PostProcess/combine_receptor.py --relaxed glad_relaxed/ --receptors ../../evaluate_dataset/CPSet/clean_receptors  --output glad_full
```

😉 Now that we finished the target-conditioned cyclic peptide generation!

## Step 4: Evaluation

The evaluation scripts used here are quite similar to those used for dataset evaluation, as shown in [Data_Evaluation](../Dataset_Evaluation/). Nonetheless, we still put the scripts used in model outputs evaluation in [Evaluation](./Evaluation) to make things clear.

### 4.1 Energy filter

First, we use rosetta ddG to check if the binding interfaces are designed reasonable

```
find {diff_full} -name "*.pdb" > diff_full.list
python Evaluation/rosetta_analysis.py -i diff_full.list -o rosetta_diff_full.sc -p 64
python Evaluation/ddG_extractor.py -i rosetta_diff_full.sc -o rosetta_ddG_diff_full.csv
```

Now, we can calculate the energy success rate by checking Rosetta ddG. To exclude energy failed designs, run:

```
python Evaluation/filter_ddG.py rosetta_ddG_diff_full.csv {diff_full} energy_failed
```

<diff_full> now contains the final structures of generated cyclic peptides.

### 4.2 Affinity

First, calculate Rosetta ddG and Vina score. **Note that** the [Vina](https://vina.scripps.edu) command line in `vina_analysis.py` needs to be changed to fit your installation.

```
find {diff_full} -name "*.pdb" > Diff_Final.list
python Evaluation/rosetta_analysis.py -i Diff_Final.list -o Rosetta_Diff_Final.sc -c <cores>
python Evaluation/ddG_extractor.py -i Rosetta_Diff_Final.sc -o Rosetta_ddG_Diff_Final.csv
python Evaluation/vina_analysis.py -i Diff_Final.list -o Vina_Diff.csv -c <cores>
```

Then, select the best for each target

```
python Evaluation/select_best.py -i Rosetta_ddG_Diff_Final.csv -o bestddG_Diff.csv
python Evaluation/select_best.py -i Vina_Diff.csv -o bestVina_Diff.csv
```

### 4.3 Structure Validity

We use [**PLIP**](https://plip-tool.biotec.tu-dresden.de/plip-web/plip/index) to analyze interface interactions and summarize the proportions of each types of interaction.

```
python Evaluation/plip_analysis.py -i Diff_Final.list -o PLIP_Diff.csv -c <cores>
python Evaluation/plip_handeler.py PLIP_Diff.csv
```

The statistic review should be print in stdout.

We use Ramachadran plot to evaluate the validity of cyclic peptides themselves.

```
python Evaluation/rama_analysis.py -i Diff_Final.list -o Rama_Diff.csv -c <cores>
```

### 4.4 Wet-lab Compatibility

We calculate GRAVY, logP, and rTPSA to assess the synthesis feasiblity and aggregation propensity of cyclic peptides.

```
python Evaluation/water_or_oil.py -i Diff_Final.list -o WoO_Diff.csv --chain_id L -c <cores>
python Evaluation/GRAVY_calculator.py -i Diff_Final.list -o GRAVY_Diff.csv --chain_id L -c 4
```

### 4.5 Diversity, Novelty and Self-Consistency

For diversity, we cluster the generated structures by cyclic-aware rmsd and Foldseek

```
python Evaluation/mp_rmsd_cyclic.py --rmsd Diff_rmsd.npy --hist Diff_rmsd.png --cluster Diff_rmsd.csv --name_list Diff_Final.list --num_cores 48

foldseek easy-multimercluster {diff_full}  Diff_foldseek_cluster/clu  Diff_foldseek_cluster --multimer-tm-threshold 0.65 \
--chain-tm-threshold 0.5 --interface-lddt-threshold 0.65 --alignment-type 2 --cov-mode 0 --min-seq-id 0 --threads 32
```

The number of clusters can be checked by `check_unique.py`.

For Novelty, we calculate qTm value against PDB database

```
foldseek easy-multimersearch {diff_full} {path/to/PDB/database} Diff_Final_foldseek_novel/out Diff_Final_foldseek_novel --alignment-type 2 \
--tmscore-threshold 0.0 --max-seqs 1000 --format-output query,target,complexqtmscore,complexttmscore,lddt --threads 48
```

Then, analyze the qTm by

```
python Evaluation/check_qtm_max_average.py out --ignore_R --ignore_same_id
```

For self-consistency, we use [**HighFold2**](https://github.com/hongliangduan/HighFold2) to refold head-tail and disulfide cyclic peptides, and calculate the Ca scRMSD. The input file of HF2 can be batch generated by `generate_HF2_input.py`, and The scRMSD calculation can be finished by `calculate_scRMSD.py`.

`python Evaluation/calculate_scRMSD.py <folder1> <folder2> -o <scRMSD_csv>`

where <folder1> and <folder2> are generated structures and refolded structures.



