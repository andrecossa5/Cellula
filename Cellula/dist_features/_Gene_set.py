"""
_Gene_set.py: The Gene_set class
"""

import gc
from copy import deepcopy
from joblib import cpu_count
import numpy as np
import pandas as pd
from gseapy import enrichr, prerank


##
 

def rank_top(x, n=None, lowest=False):
    """
    Returns index of top(n) values in some x np.array. Both ascending or descending order 
    can be accessed. If n is not specified, all ranked index are returned.
    """
    if not isinstance(x, np.ndarray):
        x = np.array(x)
    if lowest:
        idx = np.argsort(x)[:n]
    else:
        idx = np.argsort(x)[::-1][:n]

    return idx


##


def check_equal_pars(p_original, p_new, is_tuple=False):
    """
    Check if some default parameters values have been modified or not.
    """
    if is_tuple:
        p_original = np.array([ x[1] for x in p_original.values() ])
        p_new = np.array([ x[1] for x in p_new.values() ])
    else:
        p_original = np.array([ x for x in p_original.values() ])
        p_new = np.array([ x for x in p_new.values() ])

    return (p_original == p_new).all()


##


class Gene_set:
    """
    A class to store and annotate a set of relevant genes.

    Parameters:
    -----------
    results : list-like object or dataframe
        A list-like object or dataframe representing the gene set, according to the method that produced the (ordered)
        or not gene set.
    genes_meta : dataframe
        A dataframe containing information about the genes.
    name : str or None, default=None
        A name for the gene set.
    organism : str, default='human'
        The organism to which the genes belong.

    Attributes:
    -----------
    name : str or None
        The name of the gene set.
    organism : str
        The organism to which the genes belong.
    stats : dataframe
        A dataframe containing statistics about the genes in the set.
    is_ordered : bool
        A boolean value indicating whether the gene set is ordered or not.
    is_filtered : bool
        A boolean value indicating whether the gene set has been filtered or not.
    filtered : dict
        A dictionary containing the filtered gene sets.
    filter_params : dict
        A dictionary containing the parameters for filtering the gene set.
    original_f : dict
        A dictionary containing the original filter parameters.
    rank_sort_params : dict
        A dictionary containing the parameters for ranking and sorting the gene set.
    original_s : dict
        A dictionary containing the original rank and sort parameters.
    ORA : dict
        A dictionary containing the results of the over-representation analysis.
    GSEA : dict
        A dictionary containing the results of the gene set enrichment analysis.

    Methods:
    --------
    __init__(self, results, genes_meta, name=None, organism='human'):
        Set attrs. Results can be a list-like object or a df, according to the method that produced the 
        (ordered) or not gene set.
    filter_rank_genes(self, filtering=False, rank_sort=True, out=True, only_genes=False, filt_kwargs=None, sort_kwargs=None):
        Filters and sorts the gene set based on the specified parameters.
    compute_ORA(self, key='Default_ORA', by='Adjusted P-value', collection='GO_Biological_Process_2021', n_out=50):
        Perform ORA (Over-Representation Analysis)
    compute_GSEA(self, collection='KEGG_2019_Human', n_out=50, permutation_type='phenotype', pval_t=0.05, method='log2_ratio_of_classes')
        Perform GSEA (Gene-Set Enrichment Anlysis).
    
    """

    def __init__(self, results, genes_meta, name=None, organism='human'):
        """
        Set attrs. Results can be a list-like object or a df, according to the method that produced the 
        (ordered) or not gene set.
        """
        self.name = name
        self.organism = organism

        # Checks
        if isinstance(results, list): 
            filtered_results = [ x for x in results if x in genes_meta.index ] # Curated signatures could contain undetected genes
            if len(filtered_results) == 0:
                raise ValueError("All passed genes were not detected in data.")
            self.stats = genes_meta.loc[
                    filtered_results, ['percent_cells', 'highly_variable_features', 'mean', 'var']
            ]
            self.is_ordered = False
            self.filtered = {}
        else:
            self.stats = results.join(
                genes_meta.loc[:, 
                    ['percent_cells', 'highly_variable_features', 'mean', 'var']
                ]
            )
            self.is_ordered = True
            self.is_filtered = False
            self.filtered = {}

        # Params
        self.filter_params = { 
            'effect_size' : [ '>', 0 ], # All of em
            'evidence' : [ '<', 0.1 ], # 10% FDR
            'perc_FC' : [ '>', 1 ] # Whichever difference
        }
        self.original_f = deepcopy(self.filter_params)

        self.rank_sort_params = { 
            'n' : 100,
            'by' : 'effect_size'
        }
        self.original_s = deepcopy(self.rank_sort_params)

        self.ORA = {}
        self.GSEA = {}
    
    ##
    
    def filter_rank_genes(self, filtering=False, rank_sort=True, out=True, only_genes=False,
        filt_kwargs=None, sort_kwargs=None):
        """
        Filters and and sort gs stats. 
        """
        if self.is_ordered:
            effect_type = self.stats['effect_type'].unique()[0]
        else:
            raise ValueError('Filtering and sorting operations can be performed only on ordered gene sets.')

        if effect_type not in ['log2FC', 'loading']:
            raise ValueError('Unknown effect type.')

        # Get stats
        df = self.stats

        # Filter
        if effect_type == 'log2FC' and filtering:
            if filt_kwargs is None:
                pass
            elif isinstance(filt_kwargs, dict): 
                for k, v in filt_kwargs.items():
                    if k in self.filter_params:
                        self.filter_params[k][1] = v 
                    else:
                        print(f'{k} is not available for filter querying.')
            else:
                raise ValueError(
                    '''
                    filtering should be None, or a dictionary with one (or more) of the following key:value pairs:
                    1) effect_size : lower treshold on the es.
                    2) evidence : upper treshold on the evidence. Used only for evidence type: FDR;
                    3) perc_FC : minimum FC between cells percentages expressing an expression feature. 
                    N.B. default filt_kwargs works fine in most cases, and is applied only for DE results. 
                    '''
                )

            query = ' & '.join(f'{k} {op} {tresh}' for k, (op, tresh) in self.filter_params.items())
            df = df.query(query)

        elif effect_type != 'log2FC' and filtering:
            raise ValueError(
                '''
                Filtering implemented only for DE results. Other ordered Gene_sets can 
                only be only sorted and sliced via rank_top.
                '''
            )
        else:
            query = ''

        # Sort 
        if sort_kwargs is None:
            pass
        elif isinstance(sort_kwargs, dict):
            for k, v in sort_kwargs.items():
                if k in self.rank_sort_params:
                    self.rank_sort_params[k] = v 
        else:
            raise ValueError(
                '''
                sort_kwargs should be None, or a dictionary with one (or more) 
                of the following key:value pairs: 
                1) by : stat to rank features by. Default 'effect_size'. 
                2) n : n (top or low) genes to retained. Default: 100. 

                N.B. default sort_kwargs works fine in most cases. 
                '''
            )

        if rank_sort:
            by = self.rank_sort_params['by']
            n = self.rank_sort_params['n']
            lowest = False if by == 'covariate' else True

            idx = rank_top(df[by], n=n, lowest=lowest)
            df = df.iloc[idx, :]

        # Add 
        default = check_equal_pars(self.original_f, self.filter_params, is_tuple=True) 
        default &= check_equal_pars(self.original_s, self.rank_sort_params)

        if default:
            key_to_add = 'Default_ORA'    
        else:
            key_to_add = query
            key_to_add += f'|by: {by}, n: {n}'

        self.filtered[key_to_add] = df

        if only_genes:
            output = df.index.to_list()
        else:
            output = df

        self.is_filtered = True

        if out:
            return output

    ##

    def compute_ORA(self, key='Default_ORA', by='Adjusted P-value', 
        collection='GO_Biological_Process_2021', n_out=50):
        """
        Perform ORA (Over-Representation Analysis)
        """
        if key in self.filtered.keys():
            gene_list = self.filtered[key].index.to_list()
        else:
            gene_list = self.stats.index.to_list()

        results = enrichr(
            gene_list=gene_list,
            gene_sets=[collection],
            organism=self.organism, 
            outdir=None, 
        ).results

        df = results.loc[:, 
            [ 'Term', 'Overlap', 'Adjusted P-value', 'Genes' ]
        ]

        idx = rank_top(df[by], n=n_out, lowest=True)
        filtered_df = df.iloc[idx, :]
        filtered_df = filtered_df.set_index('Term')

        # Add 
        self.ORA[key] = filtered_df

        gc.collect()

    ##

    def compute_GSEA(self, covariate='effect_size', by='Adjusted P-value', 
        collection='GO_Biological_Process_2021', n_out=50):
        """
        Perform GSEA (Gene-Set Enrichment Anlysis).
        """
        if self.is_ordered:
            if self.organism == 'human':
                ranked_gene_list = self.stats[covariate]
            elif self.organism == 'mouse':
                ranked_gene_list = (
                    self.stats
                    .loc[:, [covariate]]
                    .reset_index()
                    .rename(columns={'index':'mouse'})
                )
        else:
            raise ValueError('GSEA can be performed only on ordered gene sets.')

        # Convert if necessary
        if self.organism == 'mouse':

            from gseapy import Biomart
            bm = Biomart()
            m2h = bm.query(
                dataset='mmusculus_gene_ensembl',
                attributes=['external_gene_name', 'hsapiens_homolog_associated_gene_name']
            ).rename(columns={
                'external_gene_name':'mouse', 
                'hsapiens_homolog_associated_gene_name':'human'}
            )

            # Filter and convert
            conversion_df = ranked_gene_list.merge(m2h, on='mouse', how='left').dropna(
                ).drop_duplicates('mouse', keep='last').sort_values(
                    by='effect_size', ascending=False
            )
            ranked_gene_list = conversion_df.set_index('human')['effect_size']

        results = prerank(
            rnk=ranked_gene_list,
            gene_sets=[collection],
            threads=cpu_count(),
            min_size=15,
            max_size=500,
            permutation_num=200, 
            outdir=None, 
            seed=1234,
            verbose=True,
        )

        df = results.res2d.loc[:, 
            [ 'Term', 'ES', 'NES', 'FDR q-val', 'Lead_genes' ]
        ].rename(columns={'FDR q-val' : 'Adjusted P-value'})

        idx = rank_top(df[by], n=n_out, lowest=True)
        filtered_df = df.iloc[idx, :]
        pd.options.mode.chained_assignment = None # Remove warning
        new_term = filtered_df['Term'].map(lambda x: x.split('__')[1])
        filtered_df.loc[:, 'Term'] = new_term
        filtered_df = filtered_df.set_index('Term')

        # Convert back, if necessary
        if self.organism == 'mouse':
            reformat_genes = lambda x: ';'.join([ conversion_df.loc[conversion_df['human'] == y, 'mouse'].values[0] for y in x.split(';') ])
            filtered_df['Lead_genes'] = filtered_df['Lead_genes'].map(reformat_genes)

        # Add 
        self.GSEA['original'] = filtered_df
        gc.collect()
