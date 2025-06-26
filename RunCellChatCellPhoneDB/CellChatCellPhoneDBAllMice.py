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

def processMiceDB(bdata):
    mouse_data = bdata.copy()
    cellphonedb(mouse_data,
                 groupby='predicted_labels',
                 resource_name='mouseconsensus',
                 expr_prop=0.1,
                 verbose=True, key_added='cpdb_res')
    filename = "ProcessedData/All_CellPhoneDB.h5ad"
    mouse_data.write_h5ad(filename)

def processMiceCC(bdata):
    mouse_data = bdata.copy()
    cellchat(mouse_data,
                 groupby='predicted_labels',
                 resource_name='mouseconsensus',
                 expr_prop=0.1,
                 verbose=True, key_added='cpdb_res')
    filename = "ProcessedData/All_CellChat.h5ad"
    mouse_data.write_h5ad(filename)


data = getData()
processMiceDB(data)
processMiceCC(data)