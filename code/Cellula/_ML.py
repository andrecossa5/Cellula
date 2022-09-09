from lightgbm import LGBMClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.Linear
from sklearn.pipeline import Pipeline
import numpy as np
from sklearn.model_selection import StratifiedShuffleSplit, RandomizedSearchCV
from joblib import cpu_count
import time


##


def report(results, n_top=3):
    for i in range(1, n_top + 1):
        candidates = np.flatnonzero(results["rank_test_score"] == i)
        for candidate in candidates:
            print("Model with rank: {0}".format(i))
            print(
                "Mean validation score: {0:.3f} (std: {1:.3f})".format(
                    results["mean_test_score"][candidate],
                    results["std_test_score"][candidate],
                )
            )
            print("Parameters: {0}".format(results["params"][candidate]))
            print("")


##


################################################################
n=20
key = 'logit'
################################################################


models = {

    'logit' : 
    LogisticRegression(penalty='l2', n_jobs=cpu_count(), max_iter=10000),

    'xgboost' : 
    LGBMClassifier(n_jobs=cpu_count(), learning_rate=0.01)
}

params = {

    'logit' : 
        
    {
        'logit__solver' : ['newton-cg', 'lbfgs', 'liblinear']
        'logit__c_values' : [100, 10, 1.0, 0.1, 0.01]
    },

    'xgboost' : 

    {
        "xgboost__n_estimators": np.arange(100, 600, 100),
        "xgboost__max_features": ["auto", "sqrt", "log2"],
        "xgboost__max_depth": [4, 5, 6, 7, 8, 10],
        "xgboost__importance_type": ["gain", "split"]
    }

}

pipe = Pipeline( 
    steps=[(
            key, 
            models[key]
        )]
)

random_search = RandomizedSearchCV(
    pipe, 
    param_distributions=params[key], 
    n_iter=n,
    refit=True,
    scoring='f1',
    random_state=1234,
    cv=StratifiedShuffleSplit(n_splits=5)
)

################################################################

# Here we go
import scanpy as sc

path_data = '/Users/IEO5505/Desktop/sc_pipeline_prova/data/'
adata = sc.read(path_data + 'clustered.h5ad')

top_100 = np.argsort(adata.var['var'])[::-1][:100] 
X = adata[:, top_100].X.toarray()
y = np.where(adata.obs['leiden'] == '0', 1, 0)

start = time()
random_search.fit(X, y)
print(
    "RandomizedSearchCV took %.2f seconds for %d candidates parameter settings."
    % ((time() - start), n_iter_search)
)
report(random_search.cv_results_)

################################################################

pipe.fit(X, y)
pipe.score(X, y)

################################################################

# Come lanci tutto con subprocess??

################################################################