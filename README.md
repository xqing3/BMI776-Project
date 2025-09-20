# Project description

This repository contains the analysis code used in the study "Comparing the Mathematical Basis of CellChat and CellPhoneDB for Immune Cell to Cell Communication in Diet-induced Liver Disease" (Qing, Qian, and Gupta, 2025). The code performs cell–cell communication inference using LIANA+, comparing CellChat and CellPhoneDB on single-cell RNA-seq data from liver non-parenchymal cells of diet-induced MASLD/MASH mice. Scripts include scRNA-seq preprocessing (Scanpy), cell type annotation (CellTypist), CCC inference using LIANA+, and downstream comparative analyses focused on scoring functions, expression averages, and significant interaction overlap. Results support methodological insights into CCC tool behavior under varying ligand–receptor expression patterns.

# This folder contains all our Code

- All the folders have sepreate `READE.md` files explaining what code in spresent in the folder.
- `CellChatvsCellPhoneDBMeanScore.ipynb` helps us understand selectiveness of CellChat and CellPhoneBD.
- `CombinedCellPhoneDBCellChatBoxPlot.ipynb` helps us understand shared output of CellChat and CellPhoneBD.
- `ExploreCellChat.ipynb` is when we were exploring different tools in CellChat.
- `test_gold_standard.py` is making a dot plot for a perticular interaction.
- The full paper describing this project can be found in the BMI_776_Final_Paper.pdf

# Contributors

- Crystal (Xin) Qing (qing3@wisc.edu) —  Research Design, Data Curation, Analysis Support
- Siwei Qian (sqian32@wisc.edu) — Research Design, Data Curation, Analysis Support
- Aayush Gupta (gupta385@wisc.edu) - Research Design, Data Curation, Analysis Support
- Dr. Daifeng Wang (daifeng.wang@wisc.edu) — Project Advisor
