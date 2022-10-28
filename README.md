# Cellula

A collection of script for semi-automatic, seamless single-cell analysis (soon a Nextflow pipeline). 
Aim of this project is to provide a tool-kit for the exploration of scRNA-seq and (lentiviral-based) single-cell transcriptional lineage tracing data. Cellula's tools performs common tasks as
pre-processing, integration, clustering and cell type annotation (*)), with the following handy features:

* __Multiple methods__ for each analytical step, to guarantee that the best performing method may have the chance of being selected for a certain task
* __Evaluation steps__, to extensively benchmark each method performance at some task, _without_ any ground truth reference available.
* __Fine control at minimal effort__, thanks to flexible Command Line Interfaces (CLIs, soon becoming a Nextflow pipeline) through which users may perform either individual analytical tasks or complete analysis on their data with minimal (or none at all) coding required.
* __Focus on gene modules and signatures scoring__
* __Focus on classification models__, used as a tool to prioritize "distinguishing features" (i.e., transcriptional features with high discriminative power in classification tasks defined over categorical cells metadata) and validate clustering results (*)
* __Utils for lentiviral-based single-cell lineage tracing methods__   
* __Scalability__, to thousands of cells
* __Automatic handling of folders creation/removal__, as individual CLIs are designed to create and populate folders at a user defined location, and to switch among analytical steps without the need for user to handle this manually
* __GUIs__, to visually explore results

## Cellula workflow

After project folder initialization, the Cellula (for the time being) implements the following main tasks:

* (By sample) cell and gene quality Control (QC), followed by expression matrices merging (`0_qc.py`),
pre-processing (`1_pp.py`) and batch effects assessment (`2_kBET.py`)
* (Optional, if needed) correction of batch effects (`Harmony.py`, `Scanorama.py`, `scVI.py`, and `BBKNN.py` scripts, called if one choose to undergo data integration) followed by the assembly of the final pre-processed data (`3_integration_evaluation.py`)
* (Leiden) cell clustering at multiple, tunable resolutions, coupled to cluster markers computation (`4_clustering.py`)
* Clustering solutions evaluation and choice (`5_integration_diagnostics.py`)
* Signatures (i.e., gene sets, either defined by the user or defined by data-driven approaches) scoring (`6_signatures_scoring.py`)
* Distinguishing features ranking, trough Differential Expression (DE) and classification methods (`7_dist_features.py')
* Interactve gene expression programs visualization (`Scorer_app.py`) and distinguishing features interpetation (Gene Set Enrichment and Over-Representation analysis, `Dist_features_app.py`).

A complete documentation, with nice tutorials and stuff, will be provided where this project will reach sufficient stability for being packaged and released and distributed (this is still far from being the case... But we will get there soon :)).

For now, it is provided here a simple usage example, 
### scRNA-seq, no-data-integration example




'''bash
bash prepare_folder.sh $path_main
python 1_pp.py -p $path_main --step 0 --norm scanpy --n_HVGs 2000 --score scanpy
python 3_integration_diagnostics.py -p $path_main --step 0 --chosen red_s:original 
python 4_clustering.py -p $path_main --step 0 --range 0.2:1.0 --markers
python 5_clustering_diagnostics.py -p $path_main --step 0  
python 5_clustering_diagnostics.py -p $path_main --step 0 --chosen ...
python 6_signatures.py -p $path_main --step 0 --Hotspot --barkley --wu --scoring scanpy
python 7_dist_features.py -p $path_main --step 0 
'''

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
2. The Cellula.drawio.png sketch represents the data flow across Cellula CLIs, along with their dependencies.
3. `tests`, `docs` and `setup.py` needs to be implemented yet.
