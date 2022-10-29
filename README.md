# Cellula

A collection of script for semi-automatic, seamless single-cell analysis (soon a Nextflow pipeline). 
Aim of this project is to provide a tool-kit for the exploration of scRNA-seq and (lentiviral-based) single-cell transcriptional lineage tracing data. Cellula's tools performs common tasks as
pre-processing, integration, clustering and cell type annotation (*)), with the following handy features:

* __Multiple methods__ for each analytical step, to guarantee that the best performing method may have the chance of being selected for a certain task
* __Evaluation steps__, to extensively benchmark each method performance at some task, _without_ any ground truth reference available.
* __Fine control at minimal effort__, thanks to flexible Command Line Interfaces (CLIs, soon becoming a Nextflow pipeline) through which users may perform either individual analytical tasks or complete analysis on their data with minimal (or none at all) coding required.
* __Focus on gene modules and signatures scoring__
* __Focus on classification models__, used as a tool to prioritize "distinguishing features" (i.e., transcriptional features with high discriminative power in classification tasks defined over categorical cells metadata) and validate clustering results (*)
* __Utils for lentiviral-based single-cell lineage tracing data (sclt) analysis__   
* __Scalability__, to thousands of cells
* __Automatic handling of folders creation/removal__, as individual CLIs are designed to create and populate folders at a user defined location, and to switch among analytical steps without the need for user to handle this manually
* __GUIs__, to visually explore results

## Cellula workflow

After project folder initialization, the main Cellula (for the time being) implements the following tasks:

* (By sample) cell and gene quality Control (QC), followed by expression matrices merging (`0_qc.py`),
pre-processing (`1_pp.py`) and batch effects assessment (`2_kBET.py`)
* (Optional, if needed) correction of batch effects (`Harmony.py`, `Scanorama.py`, `scVI.py`, and `BBKNN.py` scripts, called if one choose to undergo data integration) followed by the assembly of the final pre-processed data (`3_integration_evaluation.py`)
* (Leiden) cell clustering at multiple, tunable resolutions, coupled to cluster markers computation (`4_clustering.py`)
* Clustering solutions evaluation and choice (`5_integration_diagnostics.py`)
* Signatures (i.e., gene sets, either defined by the user or defined by data-driven approaches) scoring (`6_signatures_scoring.py`)
* Distinguishing features ranking, through Differential Expression (DE) and classification methods (`7_dist_features.py')
* Interactve gene expression programs visualization (`Scorer_app.py`) and distinguishing features interpetation (Gene Set Enrichment and Over-Representation analysis, `Dist_features_app.py`)

A complete documentation, with nice tutorials and stuff, will be provided where this project will reach sufficient stability for being packaged and released and distributed (this is still far from being the case... But we will get there soon :)).

For now, it is provided here a simple usage example, including installation (*) instructions. To get better understanding on individual CLIs functionalities, modules, classes and function, for now consider CLIs help message, along with the temporary documentation provided in this folder. 

### Installation (*)

This project has has not being packaged yet, due to stability reasons. However, it is already possibile to get Cellula code and make it work on a local machine or an HPC cluster with few simple commands.
First thing first, clone this repo locally:

```bash
git clone git@github.com:andrecossa5/Cellula.git
```

Then, `cd` to `envs` and create the conda environment for your operating system (N.B.: Cellula has been tested only on Linux and macOS machines. In `Cellula/envs` you can find two _.yml_ files, storing receipts for both environments. `mamba` is used here for performance reasons, but `conda` works fine as well). 
For a Linux machine:

```bash
cd Cellula/envs
mamba env create -f environment_Linux.yml -n Cellula_example
```

After that, you have to link the cloned repository path to your newly created environment:

```
mamba activate cellula_example
mamba develop Cellula
```

That's it. If your fire up the python interpreter, you are able to `import Cellula`, you are ready to go. 

### Main folder setup

To begin a single-cell analysis, go to your preferred location, create a new folder that will host the entire project results and `cd` to it. We will refer to this folder with its absolute path, and we will assign it to an environment variable, `$path_main`.

```bash
cd ...user_path
mkdir Cellula_test
cd Cellula_test
path_main=pwd
```

Once there, we need to initialize the project folder. At the bare minimum, Cellula requires a single folder in $path_main,`matrices`. `matrices` is the folder hosting sample matrices as outputted by `CellRanger` or `STARsolo` (along with (optional) lentiviral-clones info). In this repo, the `data` folder contains an example on how `matrices` must be structured, for both scRNA-seq and sclt data. 

To setup your `$path_main` folder, then run

```bash
bash prepare_folder.sh $path_main
```

### scRNA-seq, no-data-integration example

For this data, we will use data from the `matrices` example in the matrix folder
Once `$path_main` is ok, we will perform cell (and gene) Quality Control and matrix pre-processing, running

__N.B.__ Cellula works with an __analysis step__ logic. Namely, one Cellula workflow (with its CLIs calls and results) constitutes a step among all others possible single-cell data exploration data available. However, a user is commonly interested in vary these strategies without loosing previous __steps__ results . Therefore, `Cellula`'s CLIs all have the --step argument to activate and write on a specific _step_ folder.

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

## Repo organization (developers/contributors section)

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
1. This is still a preliminary version of this project. 
2. The Cellula.drawio.html sketch represents the data flow across Cellula CLIs, along with their dependencies.
3. `tests`, `docs` and `setup.py` needs to be fully implemented yet.
