"""
Test functions.
"""

import numpy as np
import pandas as pd
from scipy.special import binom
from scipy.stats import fisher_exact
from statsmodels.sandbox.stats.multicomp import multipletests


#


def test_clone_aggregates(df_, df_gbc_sample):

    df_counts = (
        df_gbc_sample
        .loc[:, pd.IndexSlice['count',:]]
        .droplevel(None, axis=1)
    )
    df_freq = (
        df_gbc_sample
        .loc[:, pd.IndexSlice['prevalence',:]]
        .droplevel(None, axis=1)
    )
    
    assert df_counts.values.sum() == df_.shape[0]
    assert df_freq.values.min() <= 1
    assert df_freq.values.min() >= 0
    
    for s in df_freq.columns:
        clones = df_freq.index[np.random.choice(df_freq.index.size, 10)]
        for clone in clones:
            count_sample = df_.query('sample == @s').shape[0]
            count_clone = df_.query('sample == @s and GBC == @clone').shape[0]
            
            assert count_clone / count_sample == df_freq.loc[clone, s]
            
            
##


def test_clone_common(df_, df_common):

    clones = df_common.index[np.random.choice(df_common.index.size, 10)]
    samples = df_common.columns[:-2]
    
    for c in clones:
        occurrence = df_common.loc[c,'occurrence']
        L = []
        for s in samples:
            L.append(c in df_.query('sample == @s')['GBC'].values)
        
        assert occurrence == np.array(L).sum() 


##


def test_promet(df_, df_promet, dataset):
    
    df_dt = df_.query('dataset == @dataset')
    clones = df_promet.index[:10] if df_promet.shape[0] > 10 else df_promet.index
    
    for c in clones:
        PT,  lung,  met = df_promet.loc[c][-3:]
        
        n_clone_PT = df_dt.query('GBC == @c and origin == "PT"').shape[0]
        n_PT = df_dt.query('origin == "PT"').shape[0]
        n_clone_lung = df_dt.query('GBC == @c and origin == "lung"').shape[0]
        n_lung = df_dt.query('origin == "lung"').shape[0]
        
        assert PT == (n_clone_PT / n_PT)
        assert lung == (n_clone_lung / n_lung)
        assert met == (lung / PT)       


##


def compute_enrichment(df, col1, col2, target):
    """
    Compute -log10(FDR) Fisher's exact test: col1 enrichment towards col2 target category.
    """

    n = df.shape[0]
    groups = np.sort(df[col1].unique())

    target_ratio_array = np.zeros(groups.size)
    oddsratio_array = np.zeros(groups.size)
    pvals = np.zeros(groups.size)

    # Here we go
    for i, g in enumerate(groups):

        test_group = df[col1] == g
        test_target = df[col2] == target

        group_size = test_group.sum()
        group_target_size = (test_group & test_target).sum()
        target_ratio = group_target_size / group_size
        target_ratio_array[i] = target_ratio
        other_groups_state_size = (~test_group & test_target).sum()

        # Fisher
        oddsratio, pvalue = fisher_exact(
            [
                [group_target_size, group_size - group_target_size],
                [other_groups_state_size, n - other_groups_state_size],
            ],
            alternative='greater',
        )
        oddsratio_array[i] = oddsratio
        pvals[i] = pvalue

    # Correct pvals --> FDR
    pvals = multipletests(pvals, alpha=0.05, method="fdr_bh")[1]

    # Results
    results = pd.DataFrame({
        'perc_in_target' : target_ratio_array,
        'odds_ratio' : oddsratio_array,
        'FDR' : pvals,
        'enrichment' : -np.log10(pvals) 
    }).sort_values('enrichment', ascending=False)

    return results