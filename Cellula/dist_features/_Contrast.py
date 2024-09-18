"""
_Contrast.py: The Contrast class
"""

import numpy as np
import pandas as pd
from functools import reduce


##


def compo_summary(c):
    """
    Create a df of cell composition, starting from a list-like of labels
    """
    if not isinstance(c, pd.Categorical):
        c = pd.Categorical(c)

    df = pd.DataFrame().assign(
        n_cells=c.value_counts(),
        freq=c.value_counts() / c.value_counts().sum()
    )

    return df


##


class Contrast:
    """
    A class to store all info needed to process a certain contrast with Dist_features.

    Parameters
        ----------
        meta : pandas.DataFrame
            The cell metadata.
        query : str or dict
            The query to be applied to the metadata. If it is a string, it is treated as a column name in the metadata.
            If it is a dictionary, it must be of the form {group_name : query_expression}, where query_expression is a 
            valid Pandas query expression.
        description : str, optional
            A description of the contrast.

    Attributes
    ----------
    meta : pandas.DataFrame
        The cell metadata.
    query : str or dict
        The query to be applied to the metadata. If it is a string, it is treated as a column name in the metadata.
        If it is a dictionary, it must be of the form {group_name : query_expression}, where query_expression is a 
        valid Pandas query expression.
    description : str, optional
        A description of the contrast.

    Methods
    -------
    __init__(self, meta, query, description=None)
        Initializes the contrast instance with the provided metadata and query.
    """

    def __init__(self, meta, query, description=None):
        """
        Set all attrs. Query can be string (i.e., a meta column), or a dict of eval expressions.
        """

        meta = meta.reset_index()
        if isinstance(query, str):

            if query in meta.columns:
                s = meta[query]
                self.n_cells = meta.shape[0]
                self.status = 'original'
                self.query = query
                self.description = description 
                self.type = 'iterative, one vs all'
            else:
                raise KeyError('The query column is not in the provided cell metadata.')

        elif isinstance(query, dict):
            
            try:
                groups_indices = { k : meta.query(query[k]).index for k in query.keys() }
            except:
                raise ValueError('Incorrect queries. Pass correct ones!')

            sizes = np.array([ groups_indices[k].size for k in groups_indices ])
            if np.any(sizes == 0):
                raise ValueError(f'One of the groups has no cells. Change query!')

            int_two = lambda x,y: x&y
            intersection = reduce(int_two, [ set(x.to_list()) for x in groups_indices.values() ])

            if len(intersection) == 0:
                s = np.full(meta.shape[0], fill_value='to_exclude', dtype='O')
                for value_to_add, positions in groups_indices.items():
                    s[positions] = value_to_add
                s = pd.Series(s, index=meta.iloc[:,0]).astype('category')
                self.n_cells = s[s!='to_exclude'].size
                self.status = 'new'
                self.query = query
                self.description = description 
                self.type = f'{len(groups_indices)} groups, one vs each other'
            else:
                raise ValueError('Queries must specify for disjoint sets!')
        
        else:
            raise ValueError('Provide a dict or a string')

        try:
            c = pd.Categorical(s.astype(int)) # For leiden clusters, ensures correct ordering
        except:
            c = pd.Categorical(s)

        self.category = c
        self.codes = c.codes
        self.freqs = compo_summary(c)
