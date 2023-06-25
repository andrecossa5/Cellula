"""
_dist_features.py: utils for Dist_features.
"""

import os
import yaml
import numpy as np
import pandas as pd
import scanpy as sc
from itertools import product

from ._Contrast import Contrast


##


def prep_one_contrast_jobs(meta, D):
    """
    Prep all the jobs for one contrast. Put em in a list and return it.
    """
    query = D['query']
    c = Contrast(meta, query=query)

    features = []
    models = []
    params = []

    for k, v in D['methods'].items():
        if k == 'DE':
            features.append('genes')
            models.append('wilcoxon')
            params.append(None)
        else:
            ml_pars = v 
            feat_ = ml_pars['features']
            models_ = ml_pars['models']
            mode_ = ml_pars['mode']

            combos = list(product(feat_, models_, [mode_]))

            for combo in combos:
                features.append(combo[0])
                models.append(combo[1])
                params.append(combo[2])

    L = [ 
        { 'features': x, 'model' : y, 'mode' : z} \
        for x, y, z in zip(features, models, params) 
    ]

    return L, c


##


def prep_jobs_contrasts(adata, path, contrasts_name):
    """
    Load contrasts, specified in a .yml file at path_contrasts.
    """
    with open(os.path.join(path, f'{contrasts_name}.yml'), 'r') as f:
        d = yaml.load(f, Loader=yaml.FullLoader)
    
    jobs = {}
    contrasts = {}

    for f in d:
        print(f)
        for k in d[f]:
            D = d[f][k]
            print(k)
            try:
                jobs[k], contrasts[k] = prep_one_contrast_jobs(adata.obs, D)
            except:
                print(f'{f} {k} analysis was impossible to set up, possibly because of wrong contrasts specification')

    return jobs, contrasts


##


def one_hot_from_labels(y):
    """
    My one_hot encoder from a categorical variable.
    """
    if len(y.categories) > 2:
        Y = np.concatenate(
            [ np.where(y == x, 1, 0)[:, np.newaxis] for x in y.categories ],
            axis=1
        )
    else:
        Y = np.where(y == y.categories[0], 1, 0)
    
    return Y
