"""
ML.py: functions for Machine Learning models training and evaluation. 
"""

import numpy as np
import pandas as pd
from joblib import cpu_count
from scipy.sparse import issparse
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from lightgbm import LGBMClassifier
from sklearn.svm import SVC
from sklearn.pipeline import Pipeline
from sklearn.model_selection import StratifiedShuffleSplit, RandomizedSearchCV
from sklearn.metrics import *

from .._utils import rescale


##

# from sklearn.datasets import make_classification
# X, y = make_classification(1000, 30, n_informative=10, weights=[0.9,0.1])
# feature_names = [ f'f{i}' for i in range(X.shape[1]) ]

##

def classification(X, y, feature_names, key='logit', GS=False, n_combos=5, 
                score='f1', cores_model=8, cores_GS=1):
    """
    Given some input data X y, run a classification analysis in several flavours.
    """

    ########### Params
    models = {

        'logit' : 
        LogisticRegression(penalty='l2', n_jobs=cores_model, max_iter=10000),

        'xgboost' : 
        LGBMClassifier(n_jobs=cores_model, learning_rate=0.01)

    }
    params = {

        'logit' : 

        {
            'logit__solver' : ['newton-cg', 'lbfgs', 'liblinear'],
            'logit__C' : [100, 10, 1.0, 0.1, 0.01]
        },

        'xgboost' : 

        {
            "xgboost__n_estimators" : np.arange(100, 600, 100),
            "xgboost__max_depth" : [4, 5, 6, 7, 8, 10],
            "xgboost__importance_type" : ["gain", "split"]
        }

    }
    ###########

    ##

    # Train-test split 
    rng = np.random.RandomState(1234)
    sss = StratifiedShuffleSplit(n_splits=2, test_size=0.2, random_state=rng)

    if issparse(X):
        X = X.A # Densify if genes as features

    for train_index, test_index in sss.split(X, y):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]

    # Pipe
    pipe = Pipeline( 
        steps=[ 
            
            ('pp', StandardScaler()), # Always scale expression features
            (key, models[key]) 
        ]
    )

    # Train

    # Complex: Randomized grid search CV for hyperparameters tuning
    if GS:
        model = RandomizedSearchCV(
            pipe, 
            param_distributions=params[key], 
            n_iter=n_combos,
            refit=True,
            n_jobs=cores_GS,
            scoring=score,
            random_state=rng,
            cv=StratifiedShuffleSplit(n_splits=5),
            verbose=True
        )
        model.fit(X_train, y_train)
        print(model.best_params_)

    else:
        model = pipe
        model.fit(X_train, y_train)

    # Feature importance
    f = model.best_estimator_[key] if GS else model[key]
    if key == 'logit':
        effect_type = f'{key}_coef'
        effect_size = f.coef_.flatten()

    elif key == 'xgboost':
        effect_type = f'{key}_feature_importance'
        effect_size = f.feature_importances_.flatten()

    # Decision treshold tuning
    precisions, recalls, tresholds = precision_recall_curve(
        y_train, f.predict_proba(X_train)[:,1], 
    )
    f1_scores = 2 * (precisions * recalls) / (precisions + recalls)
    alpha = tresholds[np.argmax(f1_scores)]

    # Test: N.B. X_test should be scaled as X_train
    scaler = StandardScaler()
    X_test = scaler.fit_transform(X_test)

    # Scoring
    y_pred_probabilities = f.predict_proba(X_test)[:,1]
    y_pred = [ 1 if y >= alpha else 0 for y in y_pred_probabilities ]
    scores = {
        'accuracy' : [accuracy_score(y_test, y_pred)] * len(feature_names),
        'balanced_accuracy' : [balanced_accuracy_score(y_test, y_pred)] * len(feature_names),
        'precision' : [precision_score(y_test, y_pred)] * len(feature_names),
        'recall' : [recall_score(y_test, y_pred)] * len(feature_names),
        'evidence' : [f1_score(y_test, y_pred)] * len(feature_names)
    }

    # Format and return 
    df = pd.DataFrame(scores, index=feature_names).assign(
        evidence_type='f1',
        effect_size=effect_size,
        effect_size_rescaled=rescale(effect_size),
        effect_type=effect_type
    )
    df = df.sort_values(by='effect_size_rescaled', ascending=False).assign(
        rank=[ i for i in range(1, df.shape[0]+1) ]
    )

    return df


##


if __name__ == '__main__':
    classification()


