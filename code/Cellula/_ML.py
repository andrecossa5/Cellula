# _ML module

########################################################################

# Libraries
import numpy as np
from joblib import cpu_count
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from lightgbm import LGBMClassifier
from sklearn.pipeline import Pipeline
from sklearn.model_selection import StratifiedShuffleSplit, RandomizedSearchCV
from sklearn.metrics import f1_score, balanced_accuracy_score, accuracy_score
from scipy.sparse import issparse

import sys
sys.path.append('/Users/IEO5505/Desktop/pipeline/code/Cellula/') # Path to pipeline code in docker image
from _utils import *
from _dist_features import *

########################################################################

# Utils


def report_GS(results, n_top=3):
    '''
    Report Grid Search results.
    '''
    for i in range(1, n_top + 1):
        candidates = np.flatnonzero(results["rank_test_score"] == i)
        for candidate in candidates:
            print(f'Model with rank: {i}')
            print(
                f'''
                Mean validation score: 
                {results["mean_test_score"][candidate]:.3f} 
                std: {results["std_test_score"][candidate]:.3f}
                '''
            )
        
            print(f'Parameters: {results["params"][candidate]}')
            print('')


##

 
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

################################################################

# classification()

def classification(X, y, features_names, key='xgboost', GS=True, n_combos=5, 
                score='f1', cores_GS=cpu_count()):
    '''
    Given some input data, run a classification analysis in several flavours.
    '''
    # Train-test split 
    rng = np.random.RandomState(1234)
    sss = StratifiedShuffleSplit(n_splits=2, test_size=0.2, random_state=rng)

    if issparse(X):
        X = X.toarray() # Remove if genes

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
            cv=StratifiedShuffleSplit(n_splits=5)
        )
        model.fit(X_train, y_train)

        if key == 'xgboost':
            effect_size = model.best_estimator_[1].feature_importances_
        elif key == 'logit':
            effect_size = model.best_estimator_[1].coef_.reshape(len(features_names), 1)

    # Simple: once, on all training set, default hyperparameters
    else:
        model = pipe
        model.fit(X_train, y_train)

        if key == 'xgboost':
            effect_size = model.steps[1][1].feature_importances_
        elif key == 'logit':
            effect_size = model.steps[1][1].coef_.reshape(len(features_names), 1)

    # Test: N.B. X_test should be scaled as X_train
    scaler = StandardScaler()
    X_test = scaler.fit_transform(X_test)

    # Scores
    y_pred = model.predict(X_test)
    if score == 'accuracy':
        evidence = accuracy_score(y_test, y_pred)
    elif score == 'balanced_accuracy':
        evidence = balanced_accuracy_score(y_test, y_pred)
    elif score == 'f1':
        evidence = f1_score(y_test, y_pred)
    else:
        raise ValueError('Unknown score for classification.')

    # Format and return 
    df = pd.DataFrame({'evidence': evidence}, index=features_names).assign(
        evidence_type=score,
        effect_size=rescale(effect_size),
        effect_type='importance' if key == 'xgboost' else 'LM_coef'
    )
    df = df.sort_values(by='effect_size', ascending=False).assign(
        rank=[ i for i in range(1, df.shape[0]+1) ]
    )

    return df

################################################################

if __name__ == '__main__':
    classification()

################################################################

