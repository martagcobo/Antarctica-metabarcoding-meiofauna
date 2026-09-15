# Antarctic Meiofauna Metabarcoding 

This repository contains the R scripts and datasets used for the analysis of meiofaunal communities in the Antarctic region using 18S rRNA metabarcoding. It includes comparative analyses with previously published surveys from different biogeographical regions.

* **Authors:** Anonymous for peer-review
* **License:** [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
* **Last update:** 15/09/2026 (Updated repository structure and documentation following revision)

---

## Repository Structure

```text
.
├── data/
│   ├── Antarctica marine MBC samples.csv
│   ├── Antarctica_community.csv
│   ├── ASVs_metadata100_gaps.csv
│   ├── Data18S_*.csv (10 cross-study comparative files)
│   ├── edges_clean.csv
│   ├── stations network.csv
│   └── taxa_counts_meiofauna.csv
└── script/
    └── Antarctica_script.R


* **`data/`**: Contains all raw and processed data matrices required to run the workflow.
* **`script/`**: Contains the main R script executing ecological, statistical, and spatial analyses.
*(Bioinformatic processing prior to ecological analysis followed the NGSmeioR workflow: https://github.com/amartinezgarcia/NGSmeioR.git)*

---

## Data Dictionary (`data/`)

1. Primary Antarctic Survey Datasets

Antarctica_community.csv: Matrix of Amplicon Sequence Variant (ASV) read counts per sample for the present study.
Antarctica marine MBC samples.csv: Environmental and experimental sample metadata (e.g., depth, substrate type, mesh size, station IDs).
ASVs_metadata100_gaps.csv: Taxonomic assignment metadata for ASVs, including sequence identity percentages, BLAST hits, and alignment gap scores.

2. Cross-Study Comparative Datasets (18S rRNA Surveys)

Data18S_*.csv (Antarctica_v2, Asinara_v2, Cordier2021_v2, Degenhardt2021_v2, Fonseca2017, Haenel2017_v2, JondeliusAtherton2020, Kapshyna2024_v2, Mazurkiewicz2024, Polinski2019_v2): Compiled 18S read abundance tables from published literature, standardized for comparative richness (Table S7) and proportion of unidentified ASVs.

3. Network & Taxonomic Composition Data

edges_clean.csv: Pairwise edge-list indicating shared ASVs and connectivity metrics between sampling sites.
stations network.csv: Spatial node metadata containing geographical coordinates and region classifications for bipartite network plots.
taxa_counts_meiofauna.csv: Aggregated read and ASV abundance counts categorized by major meiofaunal taxonomic groups (used for composition heatmaps).

4. Code & Documentation

Antarctica_script.R: Main R script executing alpha/beta diversity metrics, spatial networks, cross-study benchmarking, and data visualization.


---

## Software & Environment

All statistical analyses were executed in R. Package versions and environment parameters used for this work are recorded via `sessionInfo()` at the end of `scripts/Antarctica_script.R`.

### Package Dependencies & Installation

This project requires packages from both CRAN and Bioconductor. To set up the environment and install all dependencies, run the following commands in R:

```R
# 1. Standard CRAN packages
cran_packages <- c(
  "tidyverse", "patchwork", "sjPlot", "scales", "ggrepel", 
  "vegan", "BAT", "glmmTMB", "performance", "emmeans", "ggeffects", "iNEXT", 
  "ape", "phytools", "phangorn", 
  "sf", "rnaturalearth", "marmap", "ggmap", 
  "flextable", "officer", "tidygraph", "ggraph", "pheatmap", "writexl", 
  "here"
)
install.packages(setdiff(cran_packages, rownames(installed.packages())))

# 2. Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("Biostrings", "DECIPHER", "msa"))

## License

Data and code are released under the Creative Commons Attribution 4.0 International License (CC BY 4.0).
