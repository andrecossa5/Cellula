# Cellula

A (python) scripts collection for semi-automatic, seamless single-cell analysis (soon a Nextflow pipeline. * denote features which are still pending in the current codebase).

Aim of this project is to provide a toolkit for the exploration of scRNA-seq. These tools perform common single-cell analysis tasks (i.e., data pre-processing and integration, cell clustering and annotation (*)) with the following features:

* __Multiple methods__ for each analytical step, to ensure that, for a given task, the best performing method have a fair chance of being selected 
* __Evaluation steps__, to benchmark methods performance at some task (N.B., _without_ any ground truth reference available)
* __Fine control at minimal effort__, thanks to flexible Command Line Interfaces (CLIs) through which users may perform either individual analytical tasks or complete analysis on their data, with minimal (or none at all) coding required.
* __Focus on gene modules and signatures scoring__ 
* __Focus on classification models__, used to prioritize "distinguishing features" (i.e., transcriptional features with high discriminative power in classification tasks defined over categorical cells metadata) and validate clustering results (*)
* __Utils for lentiviral-based single-cell lineage tracing data (sclt) analysis__   
* __Scalability__, to thousands of cells
* __Automatic handling of folders creation/removal__, as individual CLIs are designed to create and populate folders at some user defined location, and to switch among analysis versions (i.e., "steps") without the need for user to handle Input/Output operations manually
* __Graphical User Interfaces (GUIs)__, to visually explore results

## Cellula workflow

For the time being, the main Cellula workflow implements the following tasks:

* (By sample) cell and gene quality Control (QC), followed by expression matrices merging (`0_qc.py`),
data pre-processing (`1_pp.py`) and batch effects assessment (`2_kBET.py`)
* (Optional, if needed) correction of batch effects (`Harmony.py`, `Scanorama.py`, `scVI.py`, and `BBKNN.py` scripts) followed by the assembly of the final pre-processed data (`3_integration_evaluation.py`)
* (Leiden) cell clustering at multiple, tunable resolutions, coupled to cluster markers computation (`4_clustering.py`)
* Clustering solutions evaluation and choice (`5_clustering_diagnostics.py`)
* Signatures (i.e., gene sets, either defined by the user or retrieved by data-driven approaches) scoring (`6_signatures.py`)
* Distinguishing features ranking, through Differential Expression (DE) and classification methods (`7_dist_features.py`)
* Interactive:
    * visualization of gene expression programs (`Scorer_app.py`) 
    * distinguishing features interpetation (Gene Set Enrichment and Over-Representation analysis, `Dist_features_app.py`)

`Cellula` has been designed for command-line usage. However, individual functions and classes can be imported individually by users that look for even more flexibility.
A complete documentation (with nice tutorials and so on) will be provided when the project will reach sufficient stability for being packaged and released (we will get there soon :)).

For now, the following serves as a simple quickstart, including instructions to install (*) and run `Cellula`. To get better understanding on individual CLIs, modules, classes and function features, consider CLIs help messages, source code comments and docstrings, and the (temporary) documentation provided in this repo. 

### Installation (*)

Even if `Cellula` cannot be installed from source, it is already possibile to download its code, and make it work on a local machine or an HPC cluster with few simple commands.

In order to do that, first thing first, clone this repo locally:

```bash
git clone git@github.com:andrecossa5/Cellula.git
# or git clone https://github.com/andrecossa5/Cellula.git
```

Then, `cd` to `./Cellula/envs` and create the conda environment for your operating system (N.B.: Cellula has been tested only on Linux and macOS machines. In `Cellula/envs` you can find two _.yml_ files, storing receipts for both environments. `mamba` is used here for performance reasons, but `conda` works fine as well). 
For a Linux machine:

```bash
cd ./Cellula/envs
mamba env create -f environment_Linux.yml -n cellula_example
```

After that, you have to link the cloned repository path to your newly created environment:

```bash
mamba activate cellula_example
mamba develop . # you have to be in the cloned repo path
```

That's it. If your fire up the newly installed python interpreter, you are able to `import Cellula`, you are ready to go. 

### Main folder setup

To begin a new single-cell analysis, `cd` to a location of choice on your machine, and create a new folder, This folder will host all data and results of your project. We will refer to this __main folder__ with its absolute path, and we will assign this path to a bash environment variable, `$path_main`.

```bash
export main_folder_name=choose your name
cd #-- your choice here --#
mkdir $main_folder_name
cd $main_folder_name
path_main=`pwd`/
```

Once in `$path_main`, we need to setup this folder in order to begin with the analysis. At the bare minimum, Cellula requires __two folders__ in `$path_main`, `matrices` and `data`:

* `matrices` hosts all sample matrices for the project (i.e., (CellRanger)[https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger] or (STARsolo)[https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md] outputs, including, for each sample, 3 files: barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz). If one has to deal with sclt data, each sample directory need to store lentiviral-clones info (see below). In this repo, the `test_data` folder contains a simple example on how the `matrices` folder needs to be structured for a very simple scRNA-seq analysis, involving two samples _a_ and _b_. This folder will be directly used for the demo below. Please, use the same structure for your data.

* `data` will host all the intermediate files from `Cellula` analysis steps. In the most simple case, one can just initialize this as an empty folder. However, one may want to include other project-specific data (e.g., one a list of curated gene sets to score). In that case, just add this information with every gene set store in `.txt` format. In this repo, the `test_data` folder contains a simple example on how `data` needs to be structured in this case. This folder will be directly used for the demo below. Please, use the same structure for your data.

To setup your `$path_main` folder:

1. __`cd` to `$path_main`__ 
2. __create and fill__ `matrices` and `data` folders (or just copy in `$path_main` `test_data/matrices` and `test_data/data`, for the following demo). 
3. `cd` to your Cellula repository clone, `cd` to the `scripts` folder, and run:

```bash
bash prepare_folder.sh $path_main
```

You should be able to see two new folders created at $path_main: `results_and_plots` and `runs`.
A properly configured `$path_main` folder for a Cellula analysis should look something like this (using `tree`):

```bash
├── data
│   ├── curated_signatures
│   │   ├── Inflammation.txt
│   │   ├── Invasion.txt
│   │   ├── Metastasis.txt
│   │   ├── Proliferation.txt
│   │   ├── Quiescence.txt
│   │   └── Stemness.txt
│   └── removed_cells
├── matrices
│   ├── a
│   │   └── filtered_gene_bc_matrix
│   │       ├── barcodes.tsv.gz
│   │       ├── features.tsv.gz
│   │       └── matrix.mtx.gz
│   └── b
│       └── filtered_gene_bc_matrix
│           ├── barcodes.tsv.gz
│           ├── features.tsv.gz
│           └── matrix.mtx.gz
├── results_and_plots
│   ├── clustering
│   ├── dist_features
│   ├── pp
│   ├── signatures
│   └── vizualization
│       ├── QC
│       ├── clustering
│       ├── dist_features
│       ├── pp
│       └── signatures
└── runs
```

With $path_main correctly configured, we can proceed with the analysis.

### scRNA-seq, no-data-integration example

We will first perform Quality Control and matrix pre-processing.

__N.B. 1__ In a typical Cellula workflow (with all its CLIs calls and results) one choose ... . All these choices constitutes ... . . One is commonly interested in varying these data exploration strategies and compare their results without loosing previously precomputed runs. Therefore, `Cellula` CLIs all have a __--version__ argument to __activate and write on a__ specific _version_ folder. This way, a single place (i.e., the main folder) can store and organize all the results obtained on the same data with different, user-defined strategies/choices. 

# Shrink...
__N.B. 2__ We are still reasoning about the extent to which `Cellula` needs to be automated. Single-cell analysis is explorative in nature, and therefore it may be biased by users subjective choices. Despite our efforts to __alleviate this "subjectivity" problem__, benchmarking methods to guide users choices, a lot of things may need to be adjusted on the go. One might want to inspect every output of a Cellula CLI before running the next one, or might want to run an entire analysis with the least number of CLIs calls possible, inspecting results only at the end. In this quickstart, we propose a recipe for the second scenario, but, now and after, __double checks are by all means suggested and encouraged__.

For now, all CLIs must be called form the `script` directory (i.e., one still has to `cd` to this folder to launch these scripts in a batch job on a HPC cluster).

To perform cell and gene QC followed by data preprocessing, run:

```bash
python 0_qc.py -p $path_main --step 0 --mode filtered --qc_mode seurat
python 1_pp.py -p $path_main --step 0 --norm scanpy --n_HVGs 2000 --score scanpy
```

Here we have specifically activated a "step_0" step. Notice how the related and folders and files have been created and filled within `path_main`. 

After pre-processing, in this case we will skip batch effects evaluation and data integration sections, as _a_ and _b_ samples come from the same experiment, lab and sequencing run (tutorials on how to handle more complicated situations leveraging `Cellula` functionalities at full will be soon available). Here, we will choose to retain the original 'PCA' embedding obtained by reducing (and scaling) the full gene expression matrix to the top 2000 hyper-variable genes (HVGs), a common choice in single_cell analysis (see `1_pp.py`, `2_kBET.py` and integration scripts for further details and alternatives). This data representation will be used for kNN construction, multiple resolution clustering and markers computation. All clustering solutions will be then evaluated for their quality. These three steps (i.e., choice of a cell representation to go with, clustering and initial clustering diagnostics) can be obtained by running:

```bash
python 3_integration_diagnostics.py -p $path_main --step 0 --chosen red_s:original 
python 4_clustering.py -p $path_main --step 0 --range 0.2:1.0 --markers
python 5_clustering_diagnostics.py -p $path_main --step 0  
```

The user can inspects the clustering and clustering visualization folder to visualize properties of the "best" clustering solutions obtained, and then choose one to perform the last steps of Cellula workflow. In this case we will select the 30_NN_30_0.29 solution:

```bash
python 5_clustering_diagnostics.py -p $path_main --step 0 --chosen 30_NN_30_0.29 --kNN 30_NN_30_components --rep original
```

Lastly, we will retrieve and score potentially meaningful gene sets in our data, and we will search for features (i.e., single genes, Principal Components or gene sets) able to distinguishing groups of cells in our data. Specifically, here we will look for distinguishing features discriminating individual __samples__ and __leiden clusters__ (chosen solution) with respect to all the other cells.

```bash
python 6_signatures.py -p $path_main --step 0 --Hotspot --barkley --wu --scoring scanpy
python 7_dist_features.py -p $path_main --step 0 
```

### Setup GUIs to explore Cellula output
...

## Repo organization, for whoever wants to contribute :) 

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
