#load required libraries
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import os
import liana as li
from liana.method import cellphonedb, cellchat
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

def processMiceDB(bdata,mouseName,group_type):
    mouse_data = bdata[bdata.obs[group_type] == mouseName]
    # Run CellPhoneDB on the filtered data
    cellphonedb(mouse_data,
                 groupby='predicted_labels',
                 resource_name='mouseconsensus',
                 expr_prop=0.1,
                 verbose=True, key_added='cpdb_res')
    filename = "IndividualMiceData/"+mouseName+"_CellPhoneDB.h5ad"
    mouse_data.write_h5ad(filename)

def processMiceCC(bdata,mouseName,group_type):
    mouse_data = bdata[bdata.obs[group_type] == mouseName]
    # Run CellPhoneDB on the filtered data
    cellchat(mouse_data,
                 groupby='predicted_labels',
                 resource_name='mouseconsensus',
                 expr_prop=0.1,
                 verbose=True, key_added='cpdb_res')
    filename = "IndividualMiceData/"+mouseName+"_CellChat.h5ad"
    mouse_data.write_h5ad(filename)

group_type = 'orig.ident'
data,unique_mice = getData(group_type)
for mice in tqdm(unique_mice):
    processMiceDB(data,mice,group_type)
    processMiceCC(data,mice,group_type)