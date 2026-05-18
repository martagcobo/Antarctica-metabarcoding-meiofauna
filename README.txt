# Antarctic Meiofauna Metabarcoding 

This repository contains the R scripts and datasets used for the analysis of meiofaunal communities in the Antarctic region using 18S rRNA metabarcoding. It includes comparative analyses with previously published surveys from different biogeographical regions.

# Authors: (anonymus until peer-review process is finished).
# Last update: 17/05/2026

## Overview of Analyses

The scripts provided here reproduce the main figures, diversity metrics, and statistical analyses presented in the study, structured as follows:

* **Alpha Diversity:** Assessment of both **taxonomic** and **phylogenetic** alpha diversity across different experimental variables (mesh sizes, habitat types, and depths).
* **Beta Diversity:** Evaluation of community composition shifts and turnover via **taxonomic** and **phylogenetic** beta diversity analyses.
* **Network & Spatial Analysis:** Geographical bipartite networks mapping ASV connectivity (including global and European spatial scales).
* **Cross-Study Comparison:** Statistical benchmarking comparing the proportion of unidentified ASVs (<95% identity in GenBank) in this survey against a compiled dataset of previously published literature.



## Repository Structure

* `scripts/`: Contains all R scripts for data wrangling, mapping, and plotting.
* `data/`: Contains the clean datasets (e.g., `edges_limpios.csv`, `taxa_counts_meiofauna.csv`) required to run the code. *(Note: Raw sequencing data and large files >100MB are hosted separately).*

