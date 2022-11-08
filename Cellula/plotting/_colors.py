"""
_colors.py stores functions to create Cellula colors.
"""

import pandas as pd
import scanpy as sc
import seaborn as sns


##

def create_one_col(meta, palette, cat, col_list=None):
    col_list =  sns.color_palette(palette, len(meta[cat].cat.categories))
    d = {
        l : col for l, col in zip(meta[cat].cat.categories, col_list)
    }
    return d


##


def create_colors(meta, chosen=None):
    """
    Create base colors. Samples, seq run, and optionally leiden
    """

    # Create a custom dict of colors
    colors = {
        'sample' : create_one_col(meta, 'tab20', 'sample'),
        'seq_run' : {'run_1' : 'b'}
    }
    
    # Add cluster colors, if needed
    if chosen is not None:
        c = sc.pl.palettes.default_20[:len(meta[chosen].cat.categories)]
        colors[chosen] = { cluster : color for cluster, color in zip(meta[chosen].cat.categories, c)}

    return colors

