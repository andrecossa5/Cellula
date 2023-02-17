"""
_Results_app.py: The Results_app class
"""

import numpy as np
import pandas as pd

from ._Gene_set import Gene_set
from ._Contrast import Contrast


##


def report_one(df, gs, comparison_key=None, contrast=None, model=None, n=10, show_genes=False, print_last=True):
    """
    Report ranked features and matched gene sets. If ordered show top10, top5 ORA and GSEA.
    """
    # Filter only desired info
    df_one = df.query('comparison == @comparison_key')

    if not isinstance(list(gs.values())[0], dict):
        assert model == 'wilcoxon'
        print('')
        print(f'For {comparison_key} in {contrast}, show top {n} ranked genes and associated GSEA pathways')
        for k in gs:
            if k == comparison_key:
                g = gs[k]
                print('')
                print('================================================================')
                print('')
                print(f'DE output: top {n} ranked genes')
                print('')
                print(g.stats.loc[:,
                        [ 'effect_size', 'effect_type', 'evidence', 'evidence_type', 'perc_FC' ] 
                    ].head(n)
                )
                print('')
                print('================================================================')
                print('')
                print(f'Ranked genes enrichment (GSEA)')
                print('')
                g.compute_GSEA()
                gsa = g.GSEA['original'].loc[:, ['Adjusted P-value', 'Lead_genes']].head(n)
                if show_genes:
                    for x in gsa.index:
                        print(f'--> {x}')
                        print(gsa.loc[x, 'Lead_genes'])
                        print('')
                else:
                    print(gsa)

    elif isinstance(list(gs.values())[0], dict):
        print('')
        print(f'For {comparison_key} in {contrast}, show classification output and associated ORA/GSEA pathways')
        for k in gs:
            if k == comparison_key:
                g = gs[k]
                print('')
                print('================================================================')
                print('')
                print(f'Classification output, top {n} ranked features:')
                print('')
                print(df_one.loc[:,
                        [ 'effect_size', 'effect_type', 'evidence', 'evidence_type' ] 
                    ].head(n)
                )
                print('')
                print('================================================================')
                print('')
                print(f'Top {n} ranked features annotation:')
                print('')

                last = len(list(gs[k]))-1

                for i, top_feat in enumerate(gs[k]):
                    print(top_feat)
                    print('')
                    g = gs[k][top_feat]
                    if g.is_ordered:
                        g_type = 'ordered' 
                        print(f'{top_feat} gene set type: {g_type}')
                        print(f'{top_feat} associated genes stats:')
                        print('')
                        print(g.stats.head(n))
                        print('')
                        print(f'{top_feat} associated gene set enrichment (GSEA):')
                        print('') 
                        g.compute_GSEA()
                        gsa = g.GSEA['original'].loc[:, ['Adjusted P-value', 'Lead_genes']].head(n)
                        if show_genes:
                            for x in gsa.index:
                                print(f'--> {x}')
                                print(gsa.loc[x, 'Lead_genes'])
                                print('')
                        else:
                            print(gsa)
                    else:
                        g_type = 'unordered' 
                        print(f'{top_feat} gene set type: {g_type}. Associated stats are not displayed.')
                        print('')
                        g.compute_ORA()
                        print(f'{top_feat} associated gene set enrichment (ORA):')
                        print('')
                        gsa = g.ORA['Default_ORA'].loc[:, ['Adjusted P-value', 'Genes']].head(n)
                        if show_genes:
                            for x in gsa.index:
                                print(f'--> {x}')
                                print(gsa.loc[x, 'Genes'])
                                print('')
                        else:
                            print(gsa)
                    print('')
                    
                    do_not = not print_last
            
                    if i == last and do_not:
                        pass
                    else:
                        print('----------------------------------------------------------------')
                        print('')


##


class Results:
    """
    A class to store (and interact) with Dist_features results.
    
    Parameters:
    ----------
    adata : AnnData object
        Annotated data matrix with observations (cells) in rows and features (genes) in columns.
    contrasts : dict
        A dictionary that provides contrast information.
    jobs : dict
        A dictionary that provides information about feature extraction and model specifications.
    
    Attributes:
    ----------
    matrix : AnnData object
        Annotated data matrix with observations (cells) in rows and features (genes) in columns.
    contrasts : dict
        A dictionary that provides contrast information.
    jobs : dict
        A dictionary that provides information about feature extraction and model specifications.
    results : dict
        A dictionary that stores the result data. The keys of the dictionary are in the format 
        `'|'.join([k, x['features'], x['model']])`, where `k` is a job key and `x` is an element of 
        `self.jobs[k]`. The values of the dictionary are also dictionaries, with keys 'df' and 'gs'.
    embeddings : None
    
    Methods:
    --------
    __init__(self, adata, contrasts, jobs):
        Extract features and features metadata from input adata. Prep other attributes.
    add_job_results(self, df, gs, job_key=None):
        Add a result dataframe and a gene set dictionary to the correct job key `self.results` slots.
    get_jobs_keys(self, contrast_key=None, feat_key=None, model_key=None):
        Get jobs from results.
    summary_one_comparison(self, job_key='sample|genes|wilcoxon', comparison_key='bulk_d5_tr_vs_rest', n=10, 
                            show_genes=False, show_contrast=True, print_last=True):
        Print a summary of one comparison.
    summary_one_job(self, job_key='leiden|genes|wilcoxon', n=10, show_genes=False):
        Print a summary of one entire job.
    summary_one_comparison_multiple_jobs(self, contrast_key='leiden', feat_key='PCs', model_key=None, comparison_key=None, 
                                          show_genes=True, n=10):
        Print a summary of all results for a single comparison and a defined set of jobs.
    """

    def __init__(self, adata, contrasts, jobs):
        """
        Extract features and features metadata from input adata. Prep other attributes.
        """
        # Genes
        self.matrix = adata
        self.contrasts = contrasts
        self.jobs = jobs
        self.results = { 
            '|'.join([k, x['features'], x['model']]) : \
            { 'df' : None, 'gs' : None  } \
            for k in self.jobs for x in self.jobs[k] 

        }
        self.embeddings = None
    
    ##

    def add_job_results(self, df, gs, job_key=None):
        """
        Add a result df and a Gene_set dict to the correct job_key self.results slots.
        """
        self.results[job_key]['df'] = df
        self.results[job_key]['gs'] = gs

    ##

    def get_jobs_keys(self, contrast_key=None, feat_key=None, model_key=None):
        """
        Get jobs from results.
        """
        k_list = [ x.split('|') for x in list(self.results.keys()) ]

        expr = ''
        if contrast_key is not None:
            expr += 'x[0] == contrast_key'
        if feat_key is not None:
            expr += ' and x[1] == feat_key'
        if model_key is not None:
            expr += ' and x[2] == model_key'
        if len(expr) == 0:
            expr = 'True'
            print('All jobs selected!')

        local_dict = {'contrast_key':contrast_key, 'feat_key':feat_key, 'model_key':model_key }

        l = list(filter(lambda x: eval(expr, {'x': x}, local_dict), k_list))
        l = [ '|'.join(x) for x in l] 

        return l

    ##

    def summary_one_comparison(self, job_key='sample|genes|wilcoxon', 
        comparison_key='bulk_d5_tr_vs_rest', n=10, show_genes=False, show_contrast=True, print_last=True):
        """
        Print summary one comparison.
        """
        contrast, features, model = job_key.split('|')
        job = self.results[job_key]
        df = job['df']
        gs = job['gs']

        if show_contrast:
            print('')
            print(f'Job {job_key}: analyze contrast {contrast}, using {features} as features and {model} as model.')
            print(f'Comparison: {comparison_key}.')
            print('')
            print('================================================================')

            # Contrast info
            print('')
            c = self.contrasts[contrast]
            print(f'Contrast info: {contrast}')
            print(f'n cells: {c.n_cells}')
            print(f'Query: {c.query}')
            print(f'Type: {c.type}')
            print(f'Groups: {c.category.categories.to_list()}')
            print('Groups frequencies:')
            print('')
            print(c.freqs)
            print('')
            print('================================================================')

        # Show ranked features and associated Gene_sets
        report_one(df, gs, comparison_key=comparison_key, contrast=contrast, 
            model=model, show_genes=show_genes, print_last=print_last, n=n
        )
        print('')

    ##

    def summary_one_job(self, job_key='leiden|genes|wilcoxon', n=10, show_genes=False):
        """
        Print summary one entire job.
        """
        contrast, features, model = job_key.split('|')
        job = self.results[job_key]
        df = job['df']
        gs = job['gs']

        print('')
        print(f'Job {job_key}: analyze contrast {contrast}, using {features} as features and {model} as model.')
        print('')
        print('##')

        # Contrast info
        print('')
        c = self.contrasts[contrast]
        print(f'Contrast info: {contrast}')
        print(f'n cells: {c.n_cells}')
        print(f'Query: {c.query}')
        print(f'Type: {c.type}')
        print(f'Groups: {c.category.categories.to_list()}')
        print('Groups frequencies:')
        print('')
        print(c.freqs)
        print('')
        print('##')
        print('')
        print('')

        # For each comparison, show ranked features and associated Gene_sets
        last = len(list(gs.keys()))-1
        for i, k in enumerate(gs):
            report_one(df, gs, comparison_key=k, contrast=contrast, model=model, 
                show_genes=show_genes, n=n
            )
            print('')
            if i != last:
                print('')
                print('')
                print('')
                print('##')
                print('')
                print('')

    ##

    def summary_one_comparison_multiple_jobs(self, contrast_key='leiden', feat_key='PCs', model_key=None,
        comparison_key=None, show_genes=True, n=10):
        """
        Summary all results for a single comparison and a defined set of jobs.
        """
        keys = self.get_jobs_keys(contrast_key=contrast_key, feat_key=feat_key, model_key=model_key)

        print('')
        print(f'Showing all results for comparison {comparison_key}')
        print('')
        print('##')

        # Contrast info
        contrast, features, model = keys[0].split('|')

        print('')
        c = self.contrasts[contrast]
        print(f'Contrast info: {contrast}')
        print(f'n cells: {c.n_cells}')
        print(f'Query: {c.query}')
        print(f'Type: {c.type}')
        print(f'Groups: {c.category.categories.to_list()}')
        print('Groups frequencies:')
        print('')
        print(c.freqs)
        print('')

        for i, k in enumerate(keys):
            print_last = False

            contrast, features, model = k.split('|')
            print('##')
            print('')
            print('')
            print('')
            print(f'Job {k}: analyze contrast {contrast}, using {features} as features and {model} as model.')
        
            self.summary_one_comparison(
                job_key=k, 
                comparison_key=comparison_key, 
                show_genes=show_genes,
                show_contrast=False, 
                print_last=print_last,
                n=5
            )