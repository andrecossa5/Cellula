# Cellula

A collection of script for semi-automatic, seamless single-cell analysis (soon a modular Nextflow pipeline).

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
2. The Cellula.drawio.png sketch should represent the logic flow with which Cellula CLIs pipe one into the other, along with their main modules dependencies.
3. `tests`, `docs` and `setup.py` needs to be implemented yet.
