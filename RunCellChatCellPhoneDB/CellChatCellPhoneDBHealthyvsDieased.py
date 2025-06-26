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

def getData():
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
    return bdata

def processMiceDB(bdata,name):
    mouse_data = bdata.copy()
    cellphonedb(mouse_data,
                 groupby='predicted_labels',
                 resource_name='mouseconsensus',
                 expr_prop=0.1,
                 verbose=True, key_added='cpdb_res')
    filename = "ProcessedData/"+name+"_CellPhoneDB.h5ad"
    mouse_data.write_h5ad(filename)

def processMiceCC(bdata,name):
    mouse_data = bdata.copy()
    cellchat(mouse_data,
                 groupby='predicted_labels',
                 resource_name='mouseconsensus',
                 expr_prop=0.1,
                 verbose=True, key_added='cpdb_res')
    filename = "ProcessedData/"+name+"_CellChat.h5ad"
    mouse_data.write_h5ad(filename)


data = getData()
# Create adata1 with groups 'GF_HF' and 'SPF_HF'
adata1 = data[data.obs['group'].isin(['GF_HF', 'SPF_HF'])].copy()

# Create adata2 with groups 'GF_HFVHC' and 'SPF_HFVHC'
adata2 = data[data.obs['group'].isin(['GF_HFVHC', 'SPF_HFVHC'])].copy()
for d,n in [[adata1,"Healthy"],[adata2, "Disease"]]:
    print(n)
    processMiceDB(d,n)
    processMiceCC(d,n)