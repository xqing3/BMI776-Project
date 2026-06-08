# CellChat vs CellPhoneDB Comparison Project
==============================================
## Table of Contents
---------------
* [Introduction](#introduction)
* [Usage](#usage)
* [Contributors](#contributors)
* [Project Description](#project-description)
## Introduction
------------
This repository contains the analysis code used in the study "Comparing the Mathematical Basis of CellChat and CellPhoneDB for Immune Cell to Cell Communication in Diet-induced Liver Disease" (Qing, Qian, and Gupta, 2025). The code performs cell–cell communication inference using LIANA+, comparing CellChat and CellPhoneDB on single-cell RNA-seq data from liver non-parenchymal cells of diet-induced MASLD/MASH mice.
## Usage
-----
To use this repository, please follow these steps:
1. Clone the repository using `git clone https://github.com/xqing3/CellChat_vs_CellPhoneDB_Using_Experimental_Single-Cell_Data.git`
2. Install the required libraries by running `pip install -r requirements.txt`
3. Run the analysis code using `python RunCellChatCellPhoneDB/CellChatCellPhoneDBAllMice.py`
## Contributors
------------
* Crystal (Xin) Qing (qing3@wisc.edu) — Research Design, Data Curation, Analysis Support
* Siwei Qian (sqian32@wisc.edu) — Research Design, Data Curation, Analysis Support
* Aayush Gupta (gupta385@wisc.edu) - Research Design, Data Curation, Analysis Support
* Dr. Daifeng Wang (daifeng.wang@wisc.edu) — Project Advisor
## Project Description
-------------------
This project compares the performance of CellChat and CellPhoneDB on single-cell RNA-seq data from liver non-parenchymal cells of diet-induced MASLD/MASH mice. The code includes scripts for data preprocessing, cell type annotation, and cell-cell communication inference using LIANA+. The results support methodological insights into CCC tool behavior under varying ligand–receptor expression patterns.
