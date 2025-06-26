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

def getData(group_name):
    #load data
    print("Reading data...")
    cc_data = sc.read_h5ad("ProcessedData/"+group_name+"_CellChat.h5ad")
    cd_data = sc.read_h5ad("ProcessedData/"+group_name+"_CellPhoneDB.h5ad")
    return cc_data,cd_data

def processMiceDB(mouse_data,group_name):    
    circle_plot_test = li.pl.circle_plot(adata = mouse_data,
                      uns_key='cpdb_res',
                      groupby='predicted_labels',
                      score_key='lr_means',
                      inverse_score=True,
                      target_labels='Neutrophils',
                      figure_size=(10, 10)
                      )
    fig = circle_plot_test.get_figure()
    # Save the Figure to a file
    img_name = "figures/"+group_name+"_target_Neutrophils_CellPhoneDB.png"
    fig.savefig(img_name)
    plt.close(fig)

def processMiceCC(mouse_data,group_name):
    circle_plot_test = li.pl.circle_plot(adata = mouse_data,
                      uns_key='cpdb_res',
                      groupby='predicted_labels',
                      score_key='lr_probs',
                      inverse_score=True,
                      target_labels='Neutrophils',
                      figure_size=(10, 10)
                      )
    fig = circle_plot_test.get_figure()
    # Save the Figure to a file
    img_name = "figures/"+group_name+"_target_Neutrophils_CellChat.png"
    fig.savefig(img_name)
    plt.close(fig)


for group_name in ["GF_HF","GF_HFVHC","SPF_HF","SPF_HFVHC"]:
    cc_data,cd_data = getData(group_name)
    processMiceCC(cc_data,group_name)
    processMiceDB(cd_data,group_name)