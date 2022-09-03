#!/usr/bin/python

# CLustering and modules derivation script

########################################################################

# Libraries
import sys
import time
import pickle
import pandas as pd
import numpy as np
import math
import anndata
import scanpy as sc
import wot
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.cm as cm

# Utilities

# Timer(): a timer class
class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""

class Timer:
    def __init__(self):
        self._start_time = None

    def start(self):
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self):
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None

        return round(elapsed_time, 2)


##


# cbar_generator(): funtion to generate a colorbar automathically
def cbar_generator(x, ax, palette='viridis'):
    
    cmap = plt.get_cmap(palette)
    norm = plt.Normalize(x.min(), x.max())
    sm =  cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, location='right', shrink=0.7, pad=0.1)
    cbar.ax.set_title('Trajectory')
    
    return cbar


##


def plot_trajectory(trajectory_ds, d, cell_set, ax, day=None):

    # Scatter cells
    ax.scatter(trajectory_ds.obs['FLE1'], trajectory_ds.obs['FLE2'], c='#f0f0f0',
        s=4, marker=',', edgecolors='none', alpha=0.8)
    # Scatter probability measures
    if day is not None:
        x = d.query(f'day == {day}')['FLE1']
        y = d.query(f'day == {day}')['FLE2']
        z = d.query(f'day == {day}')[f'values_{cell_set}']
    else:
        x = d['FLE1']
        y = d['FLE2']
        z = d[f'values_{cell_set}']
    ax.scatter(x, y, c=z, s=6, marker=',', edgecolors='none', vmax=d[f'values_{cell_set}'].quantile(0.975))
    # Refine
    ax.axis('off')
    if day is not None:
        ax.set_title(f'Cluster: {cell_set}, Day: {int(day)}')
    else:
        ax.set_title(f'Cluster:{cell_set}')

    return ax


##


def plot_expression_trends(trajectory_names, trajectory_trend_datasets, cell_sets, genes, figsize=(7,7)):

    # Separate gene names
    gene_names = genes.split('/')

    # Figure
    fig, ax = plt.subplots(figsize=figsize)
    # Data
    for s in cell_sets:
        idx = trajectory_names.index(s)
        mean = trajectory_trend_datasets[idx]
        mean = mean[:, gene_names]
        timepoints = mean.obs.index.values.astype(float)
        mean.obs.index = mean.obs.index.astype('category')
        # Plot
        if mean.shape[1] > 0:
            for i, g in enumerate(mean.var_names):  # each gene
                mean_i = mean[:, i].X
                # Axes
                ax.plot(timepoints, mean_i, marker='o', label=f'{g}, cluster {s}')
                ax.set(xlabel="Day", ylabel="Expression", title=genes)
                ax.legend(loc="best")

    return fig


##



# fate_viz(): visualize fate matrices
def fate_viz(name1, name2, day, ax):

    # Figure
    fate1 = fate_ds[:,name1][fate_ds.obs['day']==day].X.flatten()
    fate2 = fate_ds[:,name2][fate_ds.obs['day']==day].X.flatten()

    Nrows = len(fate1)
    x = np.zeros(Nrows)
    y = np.zeros(Nrows)
    P = np.array([[1,0],[np.cos(2*math.pi/3),math.sin(2*math.pi/3)],[math.cos(4*math.pi/3),math.sin(4*math.pi/3)]])

    for i in range(0,Nrows):
        ff = np.array([fate1[i],fate2[i],1-(fate1[i]+fate2[i])])
        x[i] = (ff @ P)[0]
        y[i] = (ff @ P)[1]

    vx = P[:,0]
    vy = P[:,1]
    t1 = plt.Polygon(P, color=(0,0,0,0.1))
    ax.add_patch(t1)
    
    ax.scatter(x,y,s=1)
    ax.scatter(vx,vy,s=1)
    ax.text(P[0,0]+.1, P[0,1], name1)
    ax.text(P[1,0]-.1, P[1,1]+.1, name2)
    ax.text(P[2,0]-.1, P[2,1]-.2, 'Other')
    ax.axis('equal')
    ax.axis('off')
    ax.set(title = f'{name1} vs {name2} at {day}')

    return ax


##


# Log-odds viz
def log_odds_viz(fate_ds, set1, set2, ax):
    
    fate1 = fate_ds[:, set1].X
    fate2 = fate_ds[:, set2].X
    p = np.log(1e-9 + np.divide(fate1, fate2, out=np.zeros_like(fate1), where=fate2 != 0)).flatten()
    ax.scatter(fate_ds.obs['day'][p>=0], p[p>=0], s=4, marker=',', color='red')
    ax.scatter(fate_ds.obs['day'][p<0], p[p<0], s=4, marker=',', color='blue')
    ax.axhline(y = 0, color ="black", linestyle ="--")
    ax.set(xlabel='Day', ylabel='Log Odds', title=f'{set1} vs {set2}')
    
    return ax


##


########################################################################

# Set IO paths and options
path_main = sys.argv[1]

path_main = '/Users/IEO5505/Desktop/AML_GFP_retention/'
path_data = path_main + '/data/'
path_results = path_main + '/results_and_plots/trajectories/WOT/'
chosen = 'leiden_0.275'

# Load adata
adata = sc.read(path_data + 'normalized_complete.h5ad')
adata_red = sc.read(path_data + 'clustered.h5ad')
tmaps_path = path_results + 'tmaps/serum'
fle = pd.read_csv(path_data + 'fle_emb.csv', index_col=0)

########################################################################

# Build the TransportMapModel and create populations
tmap_model = wot.tmap.TransportMapModel.from_directory(tmaps_path)
cell_sets = { s : adata.obs.loc[adata.obs[chosen] == s].index.to_list() for s in adata.obs[chosen].cat.categories }
populations = tmap_model.population_from_cell_sets(cell_sets, at_time=5)
# Compute trajectories
trajectory_ds = tmap_model.trajectories(populations)

########################################################################

# Viz trajectories cell sets trajectories

# Interpolate FLE coordinates
nbins = 500
xrange = fle['FLE1'].min(), fle['FLE1'].max()
yrange = fle['FLE2'].min(), fle['FLE2'].max()
fle['FLE1'] = np.floor(np.interp(fle['FLE1'], [xrange[0], xrange[1]], [0, nbins - 1])).astype(int)
fle['FLE2'] = np.floor(np.interp(fle['FLE2'], [yrange[0], yrange[1]], [0, nbins - 1])).astype(int)
# Add to trajectory meta data
trajectory_ds.obs = trajectory_ds.obs.join(fle)

# Visualize all days
with PdfPages(path_results + 'trajectories/united.pdf') as pdf:
    for s in trajectory_ds.var_names:
        # Prep data
        test_cell_set = trajectory_ds.var_names == s
        trajectory_ds.obs[f'values_{s}'] = trajectory_ds[:, test_cell_set].X
        d = trajectory_ds.obs.groupby(['FLE1', 'FLE2'], as_index=False).sum()
        # Figure
        fig, ax = plt.subplots(figsize=(5, 5))
        plot_trajectory(trajectory_ds, d, s, ax)
        cbar_generator(d[f'values_{s}'], ax, palette='viridis')
        # Refine
        fig.tight_layout()
        pdf.savefig()  
        plt.close()

# Visualize separate days
with PdfPages(path_results + 'trajectories/separated.pdf') as pdf:
    for s in trajectory_ds.var_names:
        # Prep data
        test_cell_set = trajectory_ds.var_names == s
        trajectory_ds.obs[f'values_{s}'] = trajectory_ds[:, test_cell_set].X
        d = trajectory_ds.obs.groupby(['FLE1', 'FLE2'], as_index=False).sum()
        # Figure
        #i=0; j=0; nrow=2; ncol=5;
        nrow=1; ncol=4; 
        fig, axs = plt.subplots(nrow, ncol, figsize=(10, 3.5))
        # Axes
        for i, day in enumerate(trajectory_ds.obs['day'].unique()):
            plot_trajectory(trajectory_ds, d, s, axs[i], day=day) #Alternative axs[i, j]
            #j += 1
            #if j >= ncol:
            #    i+=1; j=0;
            #if i >= nrow:
            #    break
        # Refine
        fig.tight_layout()
        pdf.savefig()  
        plt.close()

########################################################################

# Expression trends

#Compute trends for all genes
trajectory_trends = wot.tmap.trajectory_trends_from_trajectory(trajectory_ds, adata)

# Save each trajectory in a separate file
for i in range(len(trajectory_trends)):
    wot.io.write_dataset(trajectory_trends[i], path_results + 'trends/' + trajectory_ds.var_names[i] + '_trends.txt')

# Read
trajectory_trend_datasets = []
trajectory_names = []
for i in range(trajectory_ds.shape[1]):
    trajectory_names.append(trajectory_ds.var.index[i]) 
    trajectory_trend_datasets.append(wot.io.read_dataset(path_results + 'trends/' + trajectory_ds.var_names[i] + '_trends.txt'))

#Visualize
cell_sets = ['0', '2'] # Clusters of choice
genes = 'OAS1/IFI6'   # Genes of choice
plot_expression_trends(trajectory_names, trajectory_trend_datasets, cell_sets, genes, figsize=(7,7))
fig.savefig(path_results + 'trends/example_expression_trends.pdf')

########################################################################

# Ancestor/descendant divergence 

# Compute 
divergence_df = wot.tmap.trajectory_divergence(adata_red, trajectory_ds, distance_metric='total_variation')
divergence_df['name'] = [ str(x) + ' vs ' + str(y) for x, y in zip(divergence_df['name1'], divergence_df['name2']) ]

# Visualize
with PdfPages(path_results + 'trajectories/divergences.pdf') as pdf:
    # Prep data
    for x in sorted(list(set(divergence_df['name1'].unique()).union(set(divergence_df['name2'].unique())))):
        test = ( divergence_df['name1'] == x ) | ( divergence_df['name2'] == x )
        df = divergence_df.loc[test, :]
        # Figure
        fig, ax = plt.subplots(figsize=(7, 6))
        # Axes
        for p, d in df.groupby('name'):
            ax.plot(d['day2'], d['distance'], '-o', label=p)
            ax.set(xlabel='Day', ylabel='Distance', title=x + ' vs others')
            ax.legend(loc='lower left', bbox_to_anchor=(0.05, 0.05), 
                ncol=len(d['name'].unique()), frameon=False)
        fig.tight_layout()
        pdf.savefig()  
        plt.close()

########################################################################

# Compute and visualize fate matrices
target_destinations = tmap_model.population_from_cell_sets(cell_sets, at_time=46)
# Model fates
fate_ds = tmap_model.fates(target_destinations)

# Visualize
day = 5
# Here we go
with PdfPages(path_results + 'fates/fate_maps_d5.pdf') as pdf:
    nrow=2; ncol=3;
    for cell_set in fate_ds.var_names:
        fig, axs = plt.subplots(nrow, ncol, figsize=(10, 6))
        i=0; j=0;
        for x in fate_ds.var_names:
            if x != cell_set:
                fate_viz(cell_set, x, day, axs[i,j])
                if j < ncol-1:
                    j +=1
                else: 
                    j=0; i+=1;
        fig.tight_layout()
        pdf.savefig()  
        plt.close()

########################################################################

# Fates log-odds

# Here we go
with PdfPages(path_results + 'fates/log_odds.pdf') as pdf:
    nrow=2; ncol=3;
    for cell_set in fate_ds.var_names:
        fig, axs = plt.subplots(nrow, ncol, figsize=(10, 6))
        i=0; j=0;
        for x in fate_ds.var_names:
            if x != cell_set:
                log_odds_viz(fate_ds, cell_set, x, axs[i,j])
                if j < ncol-1:
                    j +=1
                else: 
                    j=0; i+=1;
        fig.tight_layout()
        pdf.savefig()  
        plt.close()

########################################################################

# Fates transitions: transition table

# Define start and end populations
with PdfPages(path_results + 'transitions/transitions_tables.pdf') as pdf:
    # For each couple of adjacent timepoints
    for start, end in [ (0, 5), (5, 14), (14, 46) ]:
        # Compute stating and ending pops
        start_populations = tmap_model.population_from_cell_sets(cell_sets, at_time=start)
        end_populations = tmap_model.population_from_cell_sets(cell_sets, at_time=end)
        # Compute transition tables
        transition_table = tmap_model.transition_table(start_populations, end_populations)
        # Figure
        fig, ax = plt.subplots(figsize=(8, 7))
        sns.heatmap(transition_table.X, cmap='viridis', fmt=".2f", robust=True, annot=True,
                    xticklabels=[ x.name for x in end_populations ], 
                    yticklabels=[ x.name for x in start_populations ], ax=ax)
        ax.set(title=f'Day{start}/day{end}')
        ax.tick_params(axis='both', which='both', labelsize=10)
        fig.tight_layout()
        pdf.savefig()  
        plt.close()

########################################################################
