import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import os
import liana as li
from liana.method import cellphonedb, cellchat
from tqdm import tqdm
#liana dotplot returns a ggsave object...
from plotnine import ggsave, ggplot

mouse_name = "SPF_HFVHC_136"

cpdb = sc.read_h5ad("IndividualMiceData/"+mouse_name+"_CellPhoneDB.h5ad")
cc = sc.read_h5ad("IndividualMiceData/"+mouse_name+"_CellChat.h5ad")

x = cpdb.uns["cpdb_res"]
both = x[x.ligand == ('H2-Aa' or 'H2-Ab1')]
s = sorted(both.source.unique())
cpdb_plot = li.pl.dotplot(adata = cpdb,
              colour='lr_means',
              size='cellphone_pvals',
              inverse_size=True,
              source_labels=[i for i in s if "DC" in i],
              target_labels=['T cells'],
              figure_size=(25, 7),
              ligand_complex=['H2-Aa', 'H2-Ab1'],
              uns_key='cpdb_res'
             )

img_name = "figures/IndividualMouse_CPDBvsCC/"+mouse_name+"_CellPhoneDB.png"
cpdb_plot.save(img_name)

x = cc.uns["cpdb_res"]
both = x[x.ligand == ('H2-Aa' or 'H2-Ab1')]
s = sorted(both.source.unique())
cc_plot = li.pl.dotplot(adata = cc,
              colour='lr_probs',
              size='cellchat_pvals',
              inverse_size=True,
              source_labels=[i for i in s if "DC" in i],
              target_labels=['T cells'],
              figure_size=(25, 7),
              ligand_complex=['H2-Aa', 'H2-Ab1'],
              uns_key='cpdb_res'
             )

img_name = "figures/IndividualMouse_CPDBvsCC/"+mouse_name+"_CellChat.png"
cc_plot.save(img_name)