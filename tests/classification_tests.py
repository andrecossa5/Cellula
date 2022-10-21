# Test classification

########################################################################

# Parsing CLI args 

# Libraries
import sys 
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='test_classification',
    description='''Testing classification script.'''
)

# Add arguments

# Path_main
my_parser.add_argument(
    '--path_main',
    '-p',
    type=str,
    default=None,
    help='Path to main project. Default: None.'
)

my_parser.add_argument(
    '--model', 
    type=str,
    default=None,
    help='Model. Default: None.'
)

my_parser.add_argument(
    '--n_model', 
    type=int,
    default=None,
    help='N of cores to use for a single model instance. Default: None == cpu_count() .'
)

my_parser.add_argument(
    '--n_CV', 
    type=int,
    default=None,
    help='N of cores to use for a single GS CV instance. Default: None == cpu_count() .'
)

my_parser.add_argument(
    '--feat_type', 
    type=str,
    default='genes',
    help='Feature type. Default: Genes.'
)

my_parser.add_argument(
    '--n_features', 
    type=int,
    default=None,
    help='N of features to use. Default: None == all.'
)

my_parser.add_argument(
    '--n_obs', 
    type=int,
    default=None,
    help='N of observation to use. Default: None == all.'
)

my_parser.add_argument(
    '--n_combos', 
    type=int,
    default=50,
    help='N of GS combos to tes to use. Default: 50'
)

my_parser.add_argument(
    '--GS', 
    action='store_true',
    help='Perform CV. Default: False.'
)

args = my_parser.parse_args()

path_main = args.path_main
model = args.model
n_model = args.n_model
n_CV = args.n_CV
feat_type = args.feat_type
n_features = args.n_features
n_obs = args.n_obs
n_combos = args.n_combos

########################################################################

# Libraries

import Cellula
import random
from Cellula._utils import *
from Cellula.dist_features._dist_features import *
from Cellula.ML._ML import *

# Params
models = {

    'logit' : 
    LogisticRegression(n_jobs=n_model, max_iter=10000),

    'xgboost' : 
    LGBMClassifier(n_jobs=n_model, learning_rate=0.01)
}

params = {

    'logit' : 
        
    {
        'logit__penalty' : ['l1', 'l2'],
        'logit__solver' : ['newton-cg', 'lbfgs', 'liblinear'],
        'logit__C' : [100, 10, 1.0, 0.1, 0.01]
    },

    'xgboost' : 

    {
        "xgboost__n_estimators": np.arange(100, 600, 100),
        "xgboost__max_depth": [4, 5, 6, 7, 8, 10],
        "xgboost__importance_type": ["gain", "split"]
    }

}

########################################################################

# Test
def test():

    # Data
    adata = sc.read(path_main + 'data/step_0/lognorm.h5ad')

    # X, y, features_names
    if feat_type == 'genes':

        X = adata.X.toarray()

        if n_features is not None:
            random.seed(1234)
            idx_features = random.sample(range(1, adata.var_names.size), n_features)
            X = X[:, idx_features]
            features_names = adata.var_names[idx_features].to_list()

            if n_obs is not None:
                random.seed(1234)
                idx_obs = random.sample(range(1, adata.obs_names.size), n_obs)
                X = X[idx_obs, :]

    else:

        X = adata.obsm['X_pca']

        if n_obs is not None:
            random.seed(1234)
            idx_obs = random.sample(range(1, adata.obs_names.size), n_obs)
            X = X[idx_obs, :]


    y = np.where(adata.obs['sample'] == 'bulk_d5_tr', 1, 0)

    if n_obs is not None:
        y = y[idx_obs]


    # Here we go
    logs_path = path_main + '/runs/'
    mode = 'w' if not os.path.exists(logs_path + 'logs_classification_test.txt') else 'a'
    logger = set_logger(logs_path, 'logs_classification_test.txt', mode=mode)

    t = Timer()
    t.start()

    classification(X, y, features_names, key=model, GS=args.GS, n_combos=n_combos, 
                    score='f1', cores_GS=n_CV)

    mode = 'full' if args.GS else 'fast'
    logger.info(f'''Run with: model: {model}; {feat_type}; n_features: {n_features}; n_cells: {n_obs}; mode {mode}; {n_model} cores for the model and {n_CV} cores for CV.''')
    logger.info(f'Finished in {round(t.stop() / 60, 2)} min.')

########################################################################

# test
if __name__ == '__main__':
    test()

########################################################################


