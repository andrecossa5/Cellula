import pandas as pd 
import numpy as np 
import scanpy as sc 
import pegasus as pg
from sklearn.decomposition import PCA 
from pegasus.tools import predefined_signatures, load_signatures_from_file 
from sklearn.cluster import KMeans  
from Cellula.dist_features._signatures import scanpy_score, wot_zscore, wot_rank
from Cellula._utils import *







    




if not any([ 'cycling' == x for x in adata.obs.columns ]):
    scores = _sig_scores(adata, score_method=score_method, organism=organism)
    adata.obs = adata.obs.join(scores)
    