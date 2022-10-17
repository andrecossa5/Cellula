#/usr/bin/python

import pickle
import pandas as pd
import numpy as np
from Cellula.dist_features._Results_app import Results_app


##


# Title
st.title('Distinguishing features')

# Data
path_results = '/Users/IEO5505/Desktop/sc_pipeline_prova/results_and_plots/dist_features//step_0/' # TO FIX! Must be specified from CL...

with open(path_results + 'dist_features.txt', 'rb') as f:
    results = pickle.load(f)

# Here we go
st.write(f'Choose navigation mode.')

# Query form 
form = st.form(key='Mode')
query = form.radio(
    'Mode',
    ['one comparison and job', 'one job', 'one comparison multiple jobs']
)
submit = form.form_submit_button('Choose')

if submit:
    st.write(f'{query} chosen.')
    st.write('Select options:')

##


# One comparison and job form
if query == 'one comparison and job':

    form_1 = st.form(key='Analysis')
    job_key = form_1.selectbox(
        'Analysis',
        list(results.results.keys()), 
        key='job_key'
    )
    submit_1 = form_1.form_submit_button('Choose')

    form_2 = st.form(key='Comparison + other options')
    comparison = form_2.selectbox(
        'Comparison',
        list(results.results[job_key]['gs'].keys()), 
        key='comparison'
    )
    show_genes = form_2.selectbox(
        'Print enriched pathways gene names in extended mode',
        [False, True],
        key='show_genes'
    )
    show_contrast = form_2.selectbox(
        'Print contrast information', 
        [False, True],
        key='show_contrast'
    )
    n = form_2.selectbox(
        'n features to show',
        list(np.arange(5, 55, 5)),
        key='n'
    )

    submit_2 = form_2.form_submit_button('Run')

    if submit_2:
        results.summary_one_comparison(
            job_key=job_key, 
            comparison_key=comparison, 
            show_genes=show_genes, show_contrast=show_contrast, print_last=False, n=n
        )


##


# One job form
elif query == 'one job':

    form_2 = st.form(key='one job')

    job_key = form_2.selectbox(
        'Analysis',
        list(results.results.keys())
    )

    show_genes = form_2.selectbox(
        'Print enriched pathways gene names in extended mode',
        [False, True]
    )

    n = form_2.selectbox(
        'n features to show',
        list(np.arange(5, 55, 5))
    )

    submit_2 = form_2.form_submit_button('Run')

    if submit_2:
        results.summary_one_job(job_key=job_key, n=n, show_genes=show_genes)


##


# One comparison multiple jobs form
elif query == 'one comparison multiple jobs':
    
    col_1, col_2, col_3, col_4, col_5 = st.columns(5)

    with col_1:
        form_3 = st.form(key='Contrast')
        contrast = form_3.selectbox(
            'Contrast',
            list(np.unique([ x.split('|')[0] for x in list(results.results.keys()) ]))
        )
        submit_3 = form_3.form_submit_button('Choose')
 
    ##

    with col_2:
        form_4 = st.form(key='Features')
        feat_key = form_4.selectbox(
            'Features',
            list(np.unique([ x.split('|')[1] for x in list(results.results.keys()) if x.startswith(contrast) ])) + [None]
        )
        submit_4 = form_4.form_submit_button('Choose')

    ##

    with col_3:
        form_5 = st.form(key='Model')
        model_key = form_5.selectbox(
            'Model',
            list(
                np.unique(
                [ x.split('|')[2] for x in list(results.results.keys()) \
                    if x.split('|')[0] == contrast and x.split('|')[1] == feat_key
                ]
            )) + ['All']
        )
        submit_5 = form_5.form_submit_button('Choose')

    ##

    with col_4:
        job_keys_contrast = [ x for x in results.results.keys() if x.split('|')[0] == contrast ][0]
        alternatives = list(results.results[job_keys_contrast]['gs'].keys())

        form_6 = st.form(key='Comparison')
        comparison = form_6.selectbox(
            'Comparison',
            alternatives
        )
        submit_6 = form_6.form_submit_button('Choose')
    
    ##

    with col_5:
        form_7 = st.form(key='Genes')
        show_genes = st.selectbox(
            'Print genes',
            [False, True]
        )
        submit_7 = form_7.form_submit_button('Choose')

    if submit_7:
        if model_key == 'All':
            model_key = None
        results.summary_one_comparison_multiple_jobs(
        contrast_key=contrast, feat_key=feat_key, model_key=model_key,
        comparison_key=comparison, show_genes=show_genes
    )