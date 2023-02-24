"""
ML.py: functions for Machine Learning models training and evaluation. 
NB: The params section have to be re-built.
"""

import numpy as np
import pandas as pd
from joblib import cpu_count
from scipy.sparse import issparse
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from lightgbm import LGBMClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.model_selection import StratifiedShuffleSplit, RandomizedSearchCV
from sklearn.metrics import f1_score, balanced_accuracy_score, accuracy_score, precision_score, recall_score
import shap

from .._utils import rescale


##


def classification(X, y, feature_names, key='logit', GS=True, n_combos=5, 
                score='f1', cores_model=8, cores_GS=1):
    """
    Given some input data X y, run a classification analysis in several flavours.
    """

    ########### Params
    models = {

        'logit' : 
        LogisticRegression(penalty='l2', n_jobs=cores_model, max_iter=10000),

        'xgboost' : 
        LGBMClassifier(n_jobs=cores_model, learning_rate=0.01),

        'SVM' : 
        SVC(),

        'kNN' :
        KNeighborsClassifier(n_jobs=cores_model)

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
        },

        'SVM':

        {
            "SVM__kernel" : ["linear", "rbf", "poly"],
            "SVM__gamma" : [0.1, 1, 10, 100],
            "SVM__C" : [0.1, 1, 10, 100, 1000],
            "SVM__degree" : [0, 1, 2, 3, 4, 5, 6]
        },

        'kNN' :

        {
            "kNN__n_neighbors" : np.arange(5, 100, 25)
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

    else:
        model = pipe
        model.fit(X_train, y_train)

    # Calculate (mean) SHAP values
    X_background = shap.utils.sample(X_train, round(X_train.shape[0]*0.2))
    explainer = shap.Explainer(model.predict, X_background, seed=1234)
    shap_values = explainer(X_train, max_evals=10000000)
    mean_values = shap_values.values.mean(axis=0)

    # Test: N.B. X_test should be scaled as X_train
    scaler = StandardScaler()
    X_test = scaler.fit_transform(X_test)

    # Scores
    y_pred = model.predict(X_test)
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
        effect_size=mean_values,
        effect_type='SHAP'
    )
    df = df.sort_values(by='effect_size', ascending=False).assign(
        rank=[ i for i in range(1, df.shape[0]+1) ]
    )

    return df


##


if __name__ == '__main__':
    classification()


