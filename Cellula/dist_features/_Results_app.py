"""
_Results_app.py: The Results_app class
"""

import numpy as np
import pandas as pd
import streamlit as st

from ._Gene_set import Gene_set
from ._Contrast import Contrast

##


def restyle(df, mode='GSA'):
    """ 
    Restyle df.
    """
    if mode == 'GSA':
        color = '#7a3121' 
    elif mode == 'dist':
        color = '#155731'
    else:
        color = '#153e57'

    df = df.style.set_properties(**{
                            'background-color': 'white',
                            'color': color,
                            'border-color': 'black'
                            }
                        )
    return df


##


def report_one(df, gs, comparison_key=None, contrast=None, model=None, n=10, 
    collection='GO_Biological_Process_2021', show_genes=False, print_last=True):
    """
    Report ranked features and matched gene sets. If ordered show top10, top5 ORA and GSEA.
    """
    # Filter only desired info
    df_one = df.query('comparison == @comparison_key')

    if not isinstance(list(gs.values())[0], dict):
        assert model == 'wilcoxon'
        st.write('')
        st.markdown(f'For __{comparison_key}__ in __{contrast}__, show top {n} ranked genes and associated GSEA pathways__')
        for k in gs:
            if k == comparison_key:
                g = gs[k]
                st.write('')
                st.write('================================================================')
                st.write('')
                st.markdown(f'__DE output__: top {n} ranked genes')
                st.write('')
                st.dataframe(
                    restyle(
                        g.stats.loc[:,
                            [ 'effect_size', 'effect_type', 'evidence', 'evidence_type', 'perc_FC' ] 
                        ].head(n),
                        mode='dist'
                    )
                )
                st.write('')
                st.write('================================================================')
                st.write('')
                st.markdown(f'Ranked genes enrichment: __GSEA__')
                st.write('')
                g.compute_GSEA(collection=collection)
                gsa = pd.concat(
                    [ 
                        g.GSEA['original'].sort_values('NES', ascending=False).head(n//2),
                        g.GSEA['original'].sort_values('NES', ascending=False).tail(n//2)
                    ],
                    axis=0
                )
                gsa = gsa.loc[:, ['NES', 'Adjusted P-value', 'Lead_genes']]
                if show_genes:
                    for x in gsa.index:
                        st.write(f'--> {x}')
                        st.write(gsa.at[x, 'Lead_genes'])
                        st.write('')
                else:
                    st.write(restyle(gsa, mode='GSA'))

    elif isinstance(list(gs.values())[0], dict):
        st.write('')
        st.markdown(f'For __{comparison_key}__ in __{contrast}__, show classification output and associated ORA/GSEA pathways')
        for k in gs:
            if k == comparison_key:
                g = gs[k]
                st.write('')
                st.write('================================================================')
                st.write('')
                st.markdown(f'__Classification output__, top {n} ranked features:')
                st.write('')
                st.dataframe(
                    restyle(
                        df_one.loc[:,
                            [ 'effect_size', 'effect_type', 'evidence', 'evidence_type' ] 
                        ].head(n),
                        mode='dist'
                    )                    
                )
                st.write('')
                st.write('================================================================')
                st.write('')
                st.write(f'Top {n} ranked features annotation:')
                st.write('')

                last = len(list(gs[k]))-1

                for i, top_feat in enumerate(gs[k]):
                    st.write(top_feat)
                    st.write('')
                    g = gs[k][top_feat]
                    if g.is_ordered:
                        g_type = 'ordered' 
                        st.write(f'{top_feat} gene set type: {g_type}')
                        st.write(f'{top_feat} associated genes stats:')
                        st.write('')
                        st.dataframe(restyle(g.stats.head(n), mode='properties'))
                        st.write('')
                        st.write(f'{top_feat} associated gene set enrichment: __GSEA__')
                        st.write('') 
                        g.compute_GSEA(collection=collection)
                        gsa = pd.concat(
                            [ 
                                g.GSEA['original'].sort_values('NES', ascending=False).head(n//2),
                                g.GSEA['original'].sort_values('NES', ascending=False).tail(n//2)
                            ],
                            axis=0
                        )
                        gsa = gsa.loc[:, ['NES', 'Adjusted P-value', 'Lead_genes']]
                        if show_genes:
                            for x in gsa.index:
                                st.write(f'--> {x}')
                                st.write(gsa.at[x, 'Lead_genes'])
                                st.write('')
                        else:
                            st.dataframe(restyle(gsa, mode='GSA'))
                    else:
                        g_type = 'unordered' 
                        st.write(f'{top_feat} gene set type: {g_type}. Associated stats are not displayed.')
                        st.write('')
                        g.compute_ORA(collection=collection)
                        st.markdown(f'{top_feat} associated gene set enrichment: __ORA__')
                        st.write('')
                        gsa = g.ORA['Default_ORA'].loc[:, ['Overlap', 'Adjusted P-value', 'Genes']].head(n)
                        if show_genes:
                            for x in gsa.index:
                                st.write(f'--> {x}')
                                st.write(gsa.at[x, 'Genes'])
                                st.write('')
                        else:
                            st.dataframe(restyle(gsa, mode='GSA'))
                    st.write('')
                    
                    do_not = not print_last
            
                    if i == last and do_not:
                        pass
                    else:
                        st.write('----------------------------------------------------------------')
                        st.write('')


##


class Results_app:
    """
    A class to store (and interact) with Dist_features results.
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

    def summary_one_comparison(self, job_key='sample|genes|wilcoxon', collection='GO_Biological_Process_2021',
        comparison_key='bulk_d5_tr_vs_rest', n=10, show_genes=False, show_contrast=True, print_last=True):
        """
        Print summary one comparison.
        """
        contrast, features, model = job_key.split('|')
        job = self.results[job_key]
        df = job['df']
        gs = job['gs']

        if show_contrast:
            st.write('')
            st.markdown(f'###### Job {job_key}: analyze contrast {contrast}, using {features} as features and {model} as model.')
            st.markdown(f'###### Comparison: {comparison_key}.')
            st.write('')
            st.write('================================================================')

            # Contrast info
            st.write('')
            c = self.contrasts[contrast]
            st.markdown(f'__Contrast info__: {contrast}')
            st.write(f'n cells: {c.n_cells}')
            st.write(f'Query: {c.query}')
            st.write(f'Type: {c.type}')
            st.write(f'Groups: {c.category.categories.to_list()}')
            st.write('Groups frequencies:')
            st.write('')
            df_freqs = c.freqs
            df_freqs.index = list(df_freqs.index)
            st.dataframe(restyle(df_freqs, mode='other'))
            st.write('')
            st.write('================================================================')

        # Show ranked features and associated Gene_sets
        report_one(df, gs, comparison_key=comparison_key, contrast=contrast, 
            model=model, show_genes=show_genes, print_last=print_last, n=n, collection=collection
        )
        st.write('')

    ##

    # summary_one_comparison(results, job_key='leiden|genes|wilcoxon', comparison_key='1_vs_others', n=10)

    def summary_one_job(self, job_key='leiden|genes|wilcoxon', n=10, show_genes=False,
        collection='GO_Biological_Process_2021'):
        '''
        Print summary one entire job.
        '''
        contrast, features, model = job_key.split('|')
        job = self.results[job_key]
        df = job['df']
        gs = job['gs']

        st.write('')
        st.markdown(f'###### Job {job_key}: analyze contrast {contrast}, using {features} as features and {model} as model.')
        st.write('')
        st.write('##')

        # Contrast info
        st.write('')
        c = self.contrasts[contrast]
        st.markdown(f'__Contrast info__: {contrast}')
        st.write(f'n cells: {c.n_cells}')
        st.write(f'Query: {c.query}')
        st.write(f'Type: {c.type}')
        st.write(f'Groups: {c.category.categories.to_list()}')
        st.write('Groups frequencies:')
        st.write('')
        df_freqs = c.freqs
        df_freqs.index = list(df_freqs.index)
        st.dataframe(restyle(df_freqs, mode='other'))
        st.write('')
        st.write('##')
        st.write('')
        st.write('')

        # For each comparison, show ranked features and associated Gene_sets
        last = len(list(gs.keys()))-1
        for i, k in enumerate(gs):
            report_one(df, gs, comparison_key=k, contrast=contrast, model=model, 
                show_genes=show_genes, n=n, collection=collection
            )
            st.write('')
            if i != last:
                st.write('')
                st.write('')
                st.write('')
                st.write('##')
                st.write('')
                st.write('')

    ##

    def summary_one_comparison_multiple_jobs(self, contrast_key='leiden', feat_key='PCs', model_key=None,
        collection='GO_Biological_Process_2021', comparison_key=None, show_genes=True, n=10):
        """
        Summary all results for a single comparison and a defined set of jobs.
        """
        keys = self.get_jobs_keys(contrast_key=contrast_key, feat_key=feat_key, model_key=model_key)

        st.write('')
        st.markdown(f'##### Showing all results for comparison {comparison_key}')
        st.write('')
        st.write('##')

        # Contrast info
        contrast, features, model = keys[0].split('|')

        st.write('')
        c = self.contrasts[contrast]
        st.markdown(f'__Contrast info__: {contrast}')
        st.write(f'n cells: {c.n_cells}')
        st.write(f'Query: {c.query}')
        st.write(f'Type: {c.type}')
        st.write(f'Groups: {c.category.categories.to_list()}')
        st.write('Groups frequencies:')
        st.write('')
        df_freqs = c.freqs
        df_freqs.index = list(df_freqs.index)
        st.dataframe(restyle(df_freqs, mode='other'))
        st.write('')

        for i, k in enumerate(keys):
            print_last = False

            contrast, features, model = k.split('|')
            st.write('##')
            st.write('')
            st.write('')
            st.write('')
            st.markdown(f'###### Job {k}: analyze contrast {contrast}, using {features} as features and {model} as model.')
        
            self.summary_one_comparison(
                job_key=k, 
                comparison_key=comparison_key, 
                show_genes=show_genes,
                show_contrast=False, 
                print_last=print_last,
                n=5,
                collection=collection
            )
