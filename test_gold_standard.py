#load required libraries
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

# # verbosity: errors (0), warnings (1), info (2), hints (3)
# sc.settings.verbosity = 3
# sc.settings.set_figure_params(dpi=80, facecolor='white')

#load data
print("Reading data...")
data = sc.read_h5ad("ProcessedData/"+"GF_HF"+"_CellChat.h5ad")

dot_plot = li.pl.dotplot(adata = data,
              colour='lr_probs',
              size='cellchat_pvals',
              inverse_size=True, # we inverse sign since we want small p-values to have large sizes
              source_labels=['Kupffer cells'],
              target_labels=['Monocytes & Monocyte-derived cells'],
              figure_size=(8, 7),
              # ligand_complex=['Dll4'],
              # finally, since cpdbv2 suggests using a filter to FPs
              # we filter the pvals column to <= 0.05
              # filter_fun=lambda x: x['cellphone_pvals'] <= 0.05,
              uns_key='cpdb_res' # uns_key to use, default is 'liana_res'
             )
# fig = dot_plot.get_figure()
#     # Save the Figure to a file
#img_name = "figures/"+"GF_HF"+"_dotplot_goldstandard.png"
#fig.savefig(img_name)
# plt.close(fig)
# print(dot_plot)
# print(type(dot_plot))
print(dot_plot)
dot_plot.save(filename="GF_HF_dotplot_goldstandard.png", dpi=300)
# ggsave(plot=dot_plot, filename="GF_HF_dotplot_goldstandard.png", dpi=300)

#fig = draw(dot_plot)

# Save the figure using matplotlib
#fig.savefig("figures/GF_HF_dotplot_goldstandard.png", dpi=300)
#plt.close(fig)  # Close the figure if you don't want to display it



