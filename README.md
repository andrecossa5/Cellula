# Cellula

A collection of script for semi-automatic, seamless single-cell analysis (soon a Nextflow pipeline). Aim of this project is to provide a toolkit for the exploration of scRNA-seq data (including tasks as
pre-processing, integration, clustering and cell type annotation). This toolkit comes with the following, handy features:

* __Multiple methods for__ for each analytical step, to guarantee that the best performing method may have the chance of being selected for a certain task
* __Evaluation steps__, to extensively benchmark each method performance at some task, _without_ any ground truth reference available.
* __Fine control at minimal effort__, thanks to flexible CLIs (soon becoming a Nextflow pipeline) through which users may perform individual analytical tasks on their data with minimal coding required.
* __Focus on gene modules and signatures scoring__
* __Focus on classification models__, used as a tool to prioritize "distinguishing features" (i.e., transcriptional features with the higher discriminative power in a classification task over arbitrary categories derived from cells metadata) and validate clustering
* __Automatic handling of folders creation/removal__, as individual CLIs are designed to create and populate folders at a user defined location, and to switch among analytical steps without the need for user to handle this manually
* __Scalability__, to thousands of cells
* __GUIs__, to visually explore results

## Cellula workflow




## Repo organization

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
