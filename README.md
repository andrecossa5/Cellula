# Cellula

A collection of script for semi-automatic, seamless single-cell analysis (soon a modular Nextflow pipeline).

## Repo organization

This folder is organized as follows:

* ./envs contain the .yml file of the conda environment with all the software needed (macOS)
* ./to_be_added and ./old contain various modules that might be (re-)integrated in the future. 
* ./code contains all the actual code needed.

The ./code folder is organized as follows:

* ./code/Cellula contains the modules with functions and classes called by the CLIs.
* ./code/CLIs contains all the python CLIs that at the time being constitute the 'pipeline'.
* ./code/apps contains the .py script that opens the GUI used to visualize the results of the 7_dist_features.py script, the last CLIs that has to be called for a 'complete' Cellula analysis.
* ./code/tests contains the tests for particularly important functions/scripts.

### Comments for now
1. This is still a preliminary version of this project. 
2. The Cellula.drawio.png sketch should represent the logic flow with which Cellula CLIs should be called. It is a simple draft, for now.
3. Modules import has not being cured so much so far. Basically, all modules import everything.