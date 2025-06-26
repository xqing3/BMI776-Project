#load required libraries
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import os
import liana as li
from liana.method import cellphonedb,cellchat
from tqdm import tqdm

# verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

def getData(group_type):
    #load data
    print("Reading data...")
    adata = sc.read_h5ad("Data/Annotated.h5ad")
    
    # Store the counts for later use
    adata.layers["counts"] = adata.X.copy()
    
    # log1p normalize the data
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    
    bdata = adata
    bdata.raw = bdata
    
    unique_mice = bdata.obs[group_type].unique()
    return bdata,unique_mice

def processMice(bdata,mouseName,group_type):
    mouse_data = bdata[bdata.obs[group_type] == mouseName]
    # Run CellPhoneDB on the filtered data
    cellchat(mouse_data,
                 groupby='predicted_labels',
                 resource_name='mouseconsensus',
                 expr_prop=0.1,
                 verbose=True, key_added='cpdb_res')
    
    circle_plot_test = li.pl.circle_plot(adata = mouse_data,
                      uns_key='cpdb_res',
                      groupby='predicted_labels',
                      score_key='lr_means',
                      inverse_score=True,
                      source_labels='Neutrophils',
                      figure_size=(10, 10)
                      )
    fig = circle_plot_test.get_figure()
    # Save the Figure to a file
    img_name = "figures/"+mouseName+"_Neutrophils.png"
    fig.savefig(img_name)
    plt.close(fig)

group_type = 'group'
data,unique_mice = getData(group_type)
for mice in tqdm(unique_mice):
    processMice(data,mice,group_type)

