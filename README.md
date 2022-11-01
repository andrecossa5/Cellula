# Cellula

A (python) scripts collection for semi-automatic, seamless single-cell analysis (soon a Nextflow pipeline). 

Aim of this project is to provide a toolkit for the exploration of scRNA-seq and (lentiviral-based) single-cell transcriptional lineage tracing data. These tools perform common single-cell analysis tasks (i.e., data pre-processing and integration, cell clustering and annotation (*)) with the following features:

* __Multiple methods__ for each analytical step, to ensure that, for a given task, the best performing method have a fair chance of being selected 
* __Evaluation steps__, to benchmark methods performance at some task (N.B., _without_ any ground truth reference available)
* __Fine control at minimal effort__, thanks to flexible Command Line Interfaces (CLIs) through which users may perform either individual analytical tasks or complete analysis on their data, with minimal (or none at all) coding required.
* __Focus on gene modules and signatures scoring__ 
* __Focus on classification models__, used to prioritize "distinguishing features" (i.e., transcriptional features with high discriminative power in classification tasks defined over categorical cells metadata) and validate clustering results (*)
* __Utils for lentiviral-based single-cell lineage tracing data (sclt) analysis__   
* __Scalability__, to thousands of cells
* __Automatic handling of folders creation/removal__, as individual CLIs are designed to create and populate folders at some user defined location, and to switch among analysis versions (i.e., "steps") without the need for user to handle Input/Output operations manually
* __GUIs__, to visually explore results

## Cellula workflow

For the time being, the main Cellula workflow implements the following tasks:

* (By sample) cell and gene quality Control (QC), followed by expression matrices merging (`0_qc.py`),
data pre-processing (`1_pp.py`) and batch effects assessment (`2_kBET.py`)
* (Optional, if needed) correction of batch effects (`Harmony.py`, `Scanorama.py`, `scVI.py`, and `BBKNN.py` scripts) followed by the assembly of the final pre-processed data (`3_integration_evaluation.py`)
* (Leiden) cell clustering at multiple, tunable resolutions, coupled to cluster markers computation (`4_clustering.py`)
* Clustering solutions evaluation and choice (`5_clustering_diagnostics.py`)
* Signatures (i.e., gene sets, either defined by the user or retrieved by data-driven approaches) scoring (`6_signatures.py`)
* Distinguishing features ranking, through Differential Expression (DE) and classification methods (`7_dist_features.py')
* Interactive:
    * visualization of gene expression programs (`Scorer_app.py`) 
    * distinguishing features interpetation (Gene Set Enrichment and Over-Representation analysis, `Dist_features_app.py`)

`Cellula` has been designed for command-line usage. However, individual functions and classes can be imported individually by users seeking for even more flexibility.
A complete documentation (with nice tutorials and so on) will be provided when the project will reach sufficient stability for being packaged and released (we will get there soon :)).

For now, the following serves as a simple quickstart, including instructions to install (*) and run `Cellula`. To get better understanding on individual CLIs, modules, classes and function features, consider CLIs help messages, source code comments and docstrings, and the (temporary) documentation provided in this repo. 

### Installation (*)

Even if `Cellula` cannot be installed from source, it is already possibile to download its code, and make it work on a local machine or an HPC cluster with few simple commands.

In order to do that, first thing first, clone this repo locally:

```bash
git clone git@github.com:andrecossa5/Cellula.git
```

Then, `cd` to `./Cellula/envs` and create the conda environment for your operating system (N.B.: Cellula has been tested only on Linux and macOS machines. In `Cellula/envs` you can find two _.yml_ files, storing receipts for both environments. `mamba` is used here for performance reasons, but `conda` works fine as well). 
For a Linux machine:

```bash
cd Cellula/envs
mamba env create -f environment_Linux.yml -n cellula_example
```

After that, you have to link the cloned repository path to your newly created environment:

```bash
mamba activate cellula_example
mamba develop Cellula
```

That's it. If your fire up the newly installed python interpreter, you are able to `import Cellula`, you are ready to go. 

### Main folder setup

To begin a new single-cell analysis, `cd` to a location of choice on your machine, and create a new folder that will host the entire project results. Then, `cd` to it. We will refer to this 'main folder' with its absolute path, and we will assign it to an environment variable, `$path_main`.

```bash
cd user_path#
mkdir Cellula_test
cd Cellula_test
path_main=path_to_cellula_test#
```

Once there, we need to initialize the project folder. At the bare minimum, Cellula requires a __two folders__ in $path_main, `matrices` and `data`:

* `matrices` hosts all the sample matrices for the project, as outputted by `CellRanger` or `STARsolo` (along with (optional) lentiviral-clones info, if one has to deal with sclt data). In this repo, the `test_data` folder contains a simple example on how `matrices` needs to be structured for a very simple scRNA-seq analysis, involving two samples _a_ and _b_. Please, follow the same structure for your data. 

* `data` hosts all the main intermediate files for all `Cellula` analysis steps. In the most simple case, it can be initialized as an empty folder that will be automatically filled by each CLI. However, usually a user may want to include optional data that will be used throughout all the following anaysis "steps". For example, one may want to score a list of curated gene sets. In that case, just add this information with every gene set store in .csv format. In this repo, the `test_data` folder contains a simple example on how `data` needs to be structured in this case. Please, follow the same structure for your data.

To setup your `$path_main` folder, first create correctly structured `matrices` and `data` folders. Then,
`cd` to your repo clone, specifically in the `script` folder, and run:

```bash
bash prepare_folder.sh $path_main
```

You should be able to see two new folders created at $path_main: `results_and_plots` and `runs`.
With $path_main correctly configured, we can proceed with the analysis.

### scRNA-seq, no-data-integration example

For this data, we will use data from the `matrices` example in the matrix folder.
Once `$path_main` is ok, we will perform cell (and gene) Quality Control and matrix pre-processing, running

__N.B.__ Cellula works with an __analysis step__ logic. Namely, one Cellula workflow (with its CLIs calls and results) constitutes a step among all others possible single-cell data exploration data available. However, a user is commonly interested in vary these strategies without loosing previous __steps__ results . Therefore, `Cellula`'s CLIs all have the --step argument to activate and write on a specific _step_ folderS.

To perform cell QC, just run

```bash
python 0_qc.py -p $path_main --step 0 --mode seurat --qc_mode filtered_bc_data
```


THEN...


python 1_pp.py -p $path_main --step 0 --norm scanpy --n_HVGs 2000 --score scanpy
python 3_integration_diagnostics.py -p $path_main --step 0 --chosen red_s:original 
python 4_clustering.py -p $path_main --step 0 --range 0.2:1.0 --markers
python 5_clustering_diagnostics.py -p $path_main --step 0  
python 5_clustering_diagnostics.py -p $path_main --step 0 --chosen 1.0
python 6_signatures.py -p $path_main --step 0 --Hotspot --barkley --wu --scoring scanpy
python 7_dist_features.py -p $path_main --step 0 
```

## Repo organization (for whoever wants to contribute :))

This folder is organized as follows:

```bash
.
├── Cellula
├── apps
├── docs
├── envs
├── scripts
└── tests
```

* `envs` contains the .yml file of the conda environment needed for package setup.
* `docs` contais all documentations files.
* `tests` contains all package unit tests.
* `apps` contains the .py scripts that launch `streamlit` GUIs.  
* `scripts` contains all the CLIs which produce Cellula workflow. 
* `Cellula` contains all the modules needed by `scripts`.

### Comments for now
1. This is still a preliminary version of this project, undergoing major and minor refractoring. 
2. The Cellula.drawio.html sketch represents the data flow across Cellula CLIs, along with their dependencies.
3. `tests`, `docs` and `setup.py` needs to be fully implemented yet.
