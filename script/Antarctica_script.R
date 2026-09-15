#################################
# Antarctic meiofauna metabarcoding paper 
#
# Script by: anynonymus for peer-review.
# Last update: 15/09/2026


##### Set working directory -------------------------------------------------
# The 'here' package locates the project root using the hidden '.here' file

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)



##### Load packages ---------------------------------------------------------
# Data manipulation and visualization
library(tidyverse)
library(patchwork)
library(sjPlot)
library(scales)
library(ggrepel)   

# Community ecology and statistics
library(vegan)
library(BAT)
library(glmmTMB)
library(performance)
library(emmeans)
library(ggeffects)
library(iNEXT)

# Phylogeny and genetics
library(ape)
library(phytools)
library(Biostrings) 
library(DECIPHER)   
library(phangorn)   
library(msa)

# Spatial and mapping
library(sf)
library(rnaturalearth)
library(marmap)
library(ggmap)

# Tables and networks
library(flextable)
library(officer)
library(tidygraph) 
library(ggraph)    
library(pheatmap)
library(writexl)  
library(parallel) 


## Data Preparation -------------------------------------

# 1. LOAD RAW DATA 
ecol_base    <- read.csv2(here::here("data","Antarctica marine MBC samples.csv"))
species_base <- read.csv2(here::here("data","Data18S_Antarctica_v2.csv"))

# Load community matrix with row.names = 1
comm_raw_base <- read.csv2(here::here("data", "Antarctica_community.csv"), row.names = 1)

# Remove internal row index column 'Column1' if present (exported artifact)
if ("Column1" %in% colnames(comm_raw_base)) {
  comm_raw_base$Column1 <- NULL
}

# Explicitly define non-target samples and experimental controls to exclude:
# - Daphnia_: Daphnia obtusa internal positive control (used for tag-jumping noise filtering)
# - Blank_: Extraction and PCR negative control blanks
# - XXXIV...: Preliminary test/pilot sequencing samples outside the benthic meiofauna study scope
non_target_samples <- c(
  "Daphnia_", 
  "Blank_", 
  "XXXIV445_", 
  "XXXIV446_", 
  "XXXIV456_", 
  "XXXIV4457_", 
  "XXXIV461_"
)

# Identify and filter out controls and non-target samples by sample ID
rows_to_remove <- rownames(comm_raw_base) %in% non_target_samples
comm_raw <- comm_raw_base[!rows_to_remove, ]
comm_raw <- comm_raw[, sapply(comm_raw, is.numeric)]

# Print confirmation to console
cat("Removed", sum(rows_to_remove), "non-target control/blank samples.\n")


# 2. SYNCHRONIZE METADATA AND COMMUNITY MATRIX
common_ids <- intersect(rownames(comm_raw), ecol_base$sample_ID)

comm_raw <- comm_raw[common_ids, ]
ecol     <- ecol_base[match(common_ids, ecol_base$sample_ID), ]

stopifnot(all(rownames(comm_raw) == ecol$sample_ID))

# Filter NAs in ecological variables
ecol <- ecol[!is.na(ecol$depth) &
               !is.na(ecol$habitat_norm) &
               !is.na(ecol$Mesh), ]

comm_raw <- comm_raw[match(ecol$sample_ID, rownames(comm_raw)), ]


# 3. PREPARE SPECIES FILTERS & APPLY NOISE THRESHOLD ON RAW DATA
species_clean <- species_base[
  !is.na(species_base$ASVid) &
    !is.na(species_base$sequence) &
    nzchar(species_base$sequence), 
]

species_tot  <- species_clean[which(species_clean$Eval.result == "correct"), ]
species_tot  <- species_tot[which(species_tot$rare.asvs == "FALSE"), ]
species_tot  <- species_tot[!duplicated(species_tot$ASVid), ]

species_meio <- species_tot[which(species_tot$Tax.eco.meio %in% c("TRUE", TRUE, 1, "permanent")), ]

# Subset matrix to Meiofauna ASVs in RAW counts
comm_meio_raw <- comm_raw[, colnames(comm_raw) %in% species_meio$ASVid]

# Apply noise threshold (< 27 reads = 0) on raw reads (before rarefaction)
comm_meio_raw[comm_meio_raw < 27] <- 0

# Remove ASVs that became empty across all samples after noise filtering
comm_meio_raw <- comm_meio_raw[, colSums(comm_meio_raw) > 0]


# 4. FILTER LOW-COVERAGE SAMPLES & SYNCHRONIZE
reads_per_sample <- rowSums(comm_meio_raw)

# Realistic threshold for 18S meiofauna subset (>= 100 reads)
min_meio_reads <- 100
valid_samples  <- names(reads_per_sample[reads_per_sample >= min_meio_reads])

cat("Retained samples (>=", min_meio_reads, "meiofauna reads):", length(valid_samples), "\n")
cat("Excluded samples (low coverage):", length(reads_per_sample) - length(valid_samples), "\n")

# Subset matrices and metadata to valid samples
comm_meio_raw <- comm_meio_raw[valid_samples, ]
ecol          <- ecol[ecol$sample_ID %in% valid_samples, ]

# Synchronize row order
comm_meio_raw <- comm_meio_raw[match(ecol$sample_ID, rownames(comm_meio_raw)), ]


# 5. RAREFY MEIOFAUNA COMMUNITY
set.seed(123) # For reproducibility

min_depth_meio <- min(rowSums(comm_meio_raw))
cat("Rarefying meiofauna community to:", min_depth_meio, "reads per sample\n")

# Global rarefaction on clean meiofauna counts
comm_meio <- vegan::rrarefy(comm_meio_raw, sample = min_depth_meio)


# 6. CALCULATE FINAL METRICS & CLEAN FACTORS
# Align total raw reads with the filtered sample subset
ecol$total_reads_raw   <- rowSums(comm_raw[ecol$sample_ID, ])
ecol$meio_reads_final  <- rowSums(comm_meio)
ecol$richness_meio_rar <- rowSums(comm_meio > 0)

# Presence/Absence Matrix
comm_meio_pa <- comm_meio
comm_meio_pa[comm_meio_pa > 0] <- 1

# Factor formatting
habitat_order     <- c("epilithic", "organic", "spicule", "gravel", "sand", "silt")
ecol$habitat_norm <- factor(ecol$habitat_norm, levels = habitat_order)
ecol$Mesh         <- as.factor(ecol$Mesh)

# Drop unused factor levels to prevent downstream model convergence issues
ecol$habitat_norm <- droplevels(ecol$habitat_norm)
ecol$Mesh         <- droplevels(ecol$Mesh)



## Sensitivity analysis------------------------------------------

# 1. Perform two independent rarefactions on the clean meiofauna matrix
set.seed(123)
rar_1 <- vegan::rrarefy(comm_meio_raw, sample = min_depth_meio)

set.seed(456)
rar_2 <- vegan::rrarefy(comm_meio_raw, sample = min_depth_meio)

# 2. Convert to presence/absence and calculate Jaccard dissimilarity
d1 <- vegan::vegdist((rar_1 > 0) * 1, method = "jaccard")
d2 <- vegan::vegdist((rar_2 > 0) * 1, method = "jaccard")

# 3. Mantel test to confirm beta-diversity matrix correlation
mantel_res <- vegan::mantel(d1, d2, method = "spearman", permutations = 999)
print(mantel_res)



## ALPHA DIVERSITY ----------------------------------------------------------------------------
## A. Taxonomic Alpha Diversity--------------------------------------------------------------

## Meiofauna (Alpha Tax)-----------------------------------------------------------------------
model <- glmmTMB(richness_meio_rar ~ scale(depth) + Mesh + habitat_norm + (1 | ID), data = ecol, family = nbinom2)

# Check model assumptions
performance::check_overdispersion(model)
check_collinearity(model)

# Save model checks to PDF (Relative path used)
pdf("check_model_meio_alpha_tax.pdf", width = 10, height = 7)
performance::check_model(model)
dev.off()

# Summary and ANOVA
summary(model)
car::Anova(model)

# Post-hoc pairwise comparisons
emmeans(model, pairwise ~ Mesh, type="response")
emmeans(model, pairwise ~ habitat_norm, type="response")

## Plots (Predicted values)

# Mesh Size
prediction_mesh <- ggpredict(model, terms = "Mesh", bias_correction = TRUE)

ggplot(prediction_mesh, aes(x = x, y = predicted, group = 1)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#999999") + 
  labs(
    x = "Mesh size (µm)",
    y = "Predicted ASV Richness",
    title = "Meiofauna: Model-Adjusted Effects"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.line = element_line(color = "black")
  )

# Habitat
prediction_habitat <- ggpredict(model, terms = "habitat_norm", bias_correction = TRUE)

ggplot(prediction_habitat, aes(x = x, y = predicted, group = 1)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#999999") + 
  labs(
    x = "Type of habitat",
    y = "Predicted ASV Richness",
    title = "Meiofauna: Model-Adjusted Effects"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.line = element_line(color = "black")
  )


## Copepoda (Alpha Tax) ----------------------------------------------------------

# Identify ASVs
cop_asvs <- species_meio %>%
  filter(Best.group == "Copepoda") %>%
  pull(ASVid)

# Subset meiofauna presence/absence matrix (drop = FALSE keeps it as a matrix)
comm_cop <- comm_meio_pa[, colnames(comm_meio_pa) %in% cop_asvs, drop = FALSE]

# Calculate richness per sample
ecol$richness_cop <- rowSums(comm_cop > 0)

# Model
model_cop <- glmmTMB(richness_cop ~ scale(depth) + Mesh + habitat_norm + (1 | ID), data = ecol, family = poisson) 

# Check assumptions
performance::check_overdispersion(model_cop)
pdf("check_model_cop_output.pdf", width = 10, height = 7)
performance::check_model(model_cop)
dev.off()

# Summary and post-hoc
summary(model_cop)
car::Anova(model_cop)
emmeans(model_cop, pairwise ~ Mesh, type="response")
emmeans(model_cop, pairwise ~ habitat_norm, type="response")

# Predicted Plot: Mesh
prediction_mesh <- ggpredict(model_cop, terms = "Mesh")
ggplot(prediction_mesh, aes(x = x, y = predicted, group = 1)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#012334", alpha=0.8) + 
  labs(x = "Mesh size (µm)", y = "Predicted ASV Richness", title = "Copepoda: Model-Adjusted Effects") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.line = element_line(color = "black"))

# Predicted Plot: Habitat
prediction_habitat <- ggpredict(model_cop, terms = "habitat_norm")
ggplot(prediction_habitat, aes(x = x, y = predicted, group = 1)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#012334") + 
  labs(x = "Type of habitat", y = "Predicted ASV Richness", title = "Copepoda: Model-Adjusted Effects") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.line = element_line(color = "black"))



## Nematoda (Alpha Tax) ------------------------------------------------------------

# Identify ASVs
nem_asvs <- species_meio %>%
  filter(Best.group == "Nematoda") %>%
  pull(ASVid)

# Subset meiofauna presence/absence matrix (drop = FALSE keeps it as a matrix)
comm_nem <- comm_meio_pa[, colnames(comm_meio_pa) %in% nem_asvs, drop = FALSE]

# Calculate richness per sample
ecol$richness_nem <- rowSums(comm_nem > 0)

# model
model_nem <- glmmTMB(richness_nem ~ scale(depth) + Mesh + habitat_norm +  (1 | ID), data = ecol, family = poisson) 

# check the model
performance::check_overdispersion(model_nem)
pdf("check_model_nem_output.pdf", width = 10, height = 7)
performance::check_model(model_nem)
dev.off()

summary(model_nem)
car::Anova(model_nem)
emmeans(model_nem, pairwise ~ Mesh, type="response")
emmeans(model_nem, pairwise ~ habitat_norm, type="response")

prediction_mesh <- ggpredict(model_nem, terms = "Mesh")
ggplot(prediction_mesh, aes(x = x, y = predicted, group = 1)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#15adf8", alpha=0.8) + 
  labs(x = "Mesh size (µm)", y = "Predicted ASV Richness", title = "Nematoda: Model-Adjusted Effects") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.line = element_line(color = "black"))

prediction_habitat <- ggpredict(model_nem, terms = "habitat_norm")
ggplot(prediction_habitat, aes(x = x, y = predicted, group = 1)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#15adf8") + 
  labs(x = "Type of habitat", y = "Predicted ASV Richness", title = "Nematoda: Model-Adjusted Effects") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.line = element_line(color = "black"))


## Platyhelminthes (Alpha Tax) ----------------------------------------------------

# Identify ASVs
plat_asvs <- species_meio %>%
  filter(Best.group == "Platyhelminthes") %>%
  pull(ASVid)

# Subset meiofauna presence/absence matrix (drop = FALSE keeps it as a matrix)
comm_plat <- comm_meio_pa[, colnames(comm_meio_pa) %in% plat_asvs, drop = FALSE]

# Calculate richness per sample
ecol$richness_plat <- rowSums(comm_plat > 0)

# Model
model_plat <- glmmTMB(richness_plat ~ scale(depth) + Mesh + habitat_norm +  (1 | ID), data = ecol, family = poisson) 

# Check model
performance::check_overdispersion(model_plat)
pdf("check_model_plat_output.pdf", width = 10, height = 7)
performance::check_model(model_plat)
dev.off()

summary(model_plat)
car::Anova(model_plat)
emmeans(model_plat, pairwise ~ Mesh, type="response")
emmeans(model_plat, pairwise ~ habitat_norm, type="response")

prediction_mesh <- ggpredict(model_plat, terms = "Mesh")
ggplot(prediction_mesh, aes(x = x, y = predicted, group = 1)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#d47000", alpha=0.8) + 
  labs(x = "Mesh size (µm)", y = "Predicted ASV Richness", title = "Platyhelminthes: Model-Adjusted Effects") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.line = element_line(color = "black"))

prediction_habitat <- ggpredict(model_plat, terms = "habitat_norm")
ggplot(prediction_habitat, aes(x = x, y = predicted, group = 1)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#d47000") + 
  labs(x = "Type of habitat", y = "Predicted ASV Richness", title = "Platyhelminthes: Model-Adjusted Effects") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.line = element_line(color = "black"))



## Table with results alpha tax by group ---------------------------------------


# 1. Put all models in a list with their names
# The name on the left (in quotes) will be used to name the Excel file
model_list <- list(
  "meio" = model,
  "cop"  = model_cop,
  "nem"  = model_nem,
  "plat" = model_plat
)

# 2. Create a LOOP that will iterate through each model one by one
for (taxon_name in names(model_list)) {
  
  # Extract the current model for this loop iteration
  current_model <- model_list[[taxon_name]]
  
  cat("Processing model for:", taxon_name, "...\n")
  
  # 0. Extract basic info
  df_info <- data.frame(
    Parameter = c("Total Observations (N)", "Number of Groups (ID)"),
    Value = c(nobs(current_model), summary(current_model)$ngrps$ID)
  )
  
  # 1. Extract coefficients
  df_summary <- as.data.frame(summary(current_model)$coefficients$cond)
  df_summary$Term <- rownames(df_summary) 
  rownames(df_summary) <- NULL
  df_summary <- df_summary[, c("Term", "Estimate", "Std. Error", "z value", "Pr(>|z|)")]
  
  # 2. Extract performance metrics (with safety fallback for singular models)
  
  # A) R2 Calculation
  r2_calc <- suppressWarnings(tryCatch(performance::r2(current_model), error = function(e) NULL))
  
  # Check safely if R2 returned NA, NULL, or an atomic vector
  if (is.null(r2_calc) || !is.list(r2_calc) || is.null(r2_calc$R2_marginal) || is.na(r2_calc$R2_marginal[1])) {
    
    # Fallback: Fit fixed-effects GLM to extract Deviance Explained (D2)
    form_fixed <- update(formula(current_model), . ~ . - (1 | ID))
    glm_fit    <- tryCatch(
      glm(form_fixed, data = ecol, family = family(current_model)$family),
      error = function(e) NULL
    )
    
    if (!is.null(glm_fit)) {
      d2_val <- 1 - (glm_fit$deviance / glm_fit$null.deviance)
      df_r2  <- data.frame(R2_conditional = d2_val, R2_marginal = d2_val)
    } else {
      df_r2  <- data.frame(R2_conditional = NA, R2_marginal = NA)
    }
    
  } else {
    df_r2 <- as.data.frame(r2_calc)
  }
  
  # B) ICC Calculation
  icc_calc <- suppressWarnings(tryCatch(performance::icc(current_model), error = function(e) NULL))
  
  # Check safely for ICC as well
  if (is.null(icc_calc) || !is.list(icc_calc) || is.null(icc_calc$ICC_adjusted) || is.na(icc_calc$ICC_adjusted[1])) {
    df_icc <- data.frame(ICC = 0) # If random effect variance is zero, ICC is 0
  } else {
    df_icc <- as.data.frame(icc_calc)
  }
  
  # C) AIC/BIC
  df_aic <- data.frame(AIC = AIC(current_model), BIC = BIC(current_model))
  
  # Combine metrics
  df_performance <- cbind(df_aic, df_r2, df_icc)
  
  # 3. Extract Type II ANOVA results
  anova_res <- car::Anova(current_model)
  df_anova <- as.data.frame(anova_res)
  df_anova$Term <- rownames(df_anova)
  rownames(df_anova) <- NULL
  df_anova <- df_anova[, c("Term", "Chisq", "Df", "Pr(>Chisq)")]
  
  # 4. Extract post-hoc comparisons
  em_mesh <- emmeans(current_model, pairwise ~ Mesh, type="response")
  em_hab <- emmeans(current_model, pairwise ~ habitat_norm, type="response")
  df_pair_mesh <- as.data.frame(em_mesh$contrasts)
  df_pair_hab <- as.data.frame(em_hab$contrasts)
  
  # 5. Join everything into a list of tabs
  excel_results_list <- list(
    "0_Model_Info" = df_info,
    "1_Model_Summary" = df_summary,
    "2_Performance_Metrics" = df_performance,
    "3_ANOVA_TypeII" = df_anova,
    "4_Pairwise_Mesh" = df_pair_mesh,
    "5_Pairwise_Habitat" = df_pair_hab
  )
  
  # 6. Create file name dynamically and export
  file_name <- paste0("Meiofauna_alpha_tax_", taxon_name, ".xlsx")
  write_xlsx(excel_results_list, path = file_name)
}



## B. Phylogenetic Alpha Diversity ----------------------------------------------
## Meiofauna (Alpha Phyl)---------------------------------------------------------
# 1. FILTER SPECIES LIST TO ASVs PRESENT IN THE RAREFIED MATRIX
species_meio_phylo <- species_meio[species_meio$ASVid %in% colnames(comm_meio_pa), ]

# Clean sequence strings
species_meio_phylo$sequence <- trimws(as.character(species_meio_phylo$sequence))

species_meio_phylo <- species_meio_phylo[
  !is.na(species_meio_phylo$sequence) &
    nzchar(species_meio_phylo$sequence) &
    species_meio_phylo$sequence != "NA",
]


# 2. FASTA EXPORT & ALIGNMENT
seqs <- Biostrings::DNAStringSet(species_meio_phylo$sequence)
names(seqs) <- species_meio_phylo$ASVid

# Build phylogenetic tree (Neighbor-Joining)
alignment  <- DECIPHER::AlignSeqs(seqs)
phy_data   <- phangorn::phyDat(as.matrix(alignment), type = "DNA")
dist_matrix <- phangorn::dist.ml(phy_data)
tree        <- phangorn::NJ(dist_matrix)


# 3. CALCULATE PHYLOGENETIC DIVERSITY
# Ensure community matrix matching tree tips
comm_meio_pa_phylo <- comm_meio_pa[, colnames(comm_meio_pa) %in% tree$tip.label, drop = FALSE]

# Calculate Alpha Phylogenetic Diversity (BAT package)
pd_sample <- BAT::alpha(comm_meio_pa_phylo, tree = tree)
pd_sample <- as.data.frame(pd_sample)
colnames(pd_sample) <- "phylo.diver"


# 4. ASSIGN TO MAIN METADATA (PRESERVING EXACT ROW ORDER)
# Do NOT use merge() as it reorders rows! Use match() instead:
ecol$phylo.diver <- pd_sample[match(ecol$sample_ID, rownames(pd_sample)), "phylo.diver"]


# 5. MODEL FIT (Gamma family for strictly positive continuous PD)
ecol$depth_sc <- as.numeric(scale(ecol$depth))

mod.phylo <- glmmTMB::glmmTMB(phylo.diver ~ depth_sc + Mesh + habitat_norm + (1 | ID), data   = ecol, family = Gamma(link = "log"))

## Check the model
performance::check_overdispersion(mod.phylo)
performance::check_collinearity(mod.phylo)
performance::check_model(mod.phylo)
# dev.off()

## Summary and post-hoc
car::Anova(mod.phylo)
summary(mod.phylo)

## Post-hoc
emmeans(mod.phylo, pairwise ~ Mesh, type="response")
emmeans(mod.phylo, pairwise ~ habitat_norm, type="response")


## Box-plot (predicted)
# Mesh size

prediction_mesh <- ggpredict(mod.phylo, terms = "Mesh", bias_correction = TRUE)

ggplot(prediction_mesh, aes(x = x, y = predicted, group = 1)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#999999") + 
  labs(
    x = "Mesh size (µm)",
    y = "Predicted Phylogenetic Richness",
    title = "Meiofauna: Model-Adjusted Effects",
    subtitle = ""
  ) +
  
  # Visual cleanup
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.line = element_line(color = "black")
  )


# Habitat
prediction_habitat <- ggpredict(mod.phylo, terms = "habitat_norm", bias_correction = TRUE)

ggplot(prediction_habitat, aes(x = x, y = predicted, group = 1)) +
  
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#999999") + 
  labs(
    x = "Type of habitat",
    y = "Predicted Phylogenetic Richness",
    title = "Meiofauna: Model-Adjusted Effects",
    subtitle = ""
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.line = element_line(color = "black")
  )


## Copepoda (Alpha Phyl) --------------------------------------------------------

# 1. FILTER COPEPODA ASVs PRESENT IN THE RAREFIED MATRIX & CLEAN SEQUENCES
species_cop <- species_meio[
  species_meio$Best.group == "Copepoda" & 
    species_meio$ASVid %in% colnames(comm_meio_pa), 
]

species_cop$sequence <- trimws(as.character(species_cop$sequence))
species_cop <- species_cop[
  !is.na(species_cop$sequence) & 
    nzchar(species_cop$sequence) & 
    species_cop$sequence != "NA", 
]

# 2. FASTA EXPORT & TREE BUILDING
seqs_cop <- Biostrings::DNAStringSet(species_cop$sequence)
names(seqs_cop) <- species_cop$ASVid

alignment_cop   <- DECIPHER::AlignSeqs(seqs_cop)
phy_data_cop    <- phangorn::phyDat(as.matrix(alignment_cop), type = "DNA")
dist_matrix_cop <- phangorn::dist.ml(phy_data_cop)
tree_cop        <- phangorn::NJ(dist_matrix_cop)

# 3. COMMUNITY MATRIX SUBSET & PD CALCULATION
comm_cop_binary <- comm_meio_pa[, colnames(comm_meio_pa) %in% tree_cop$tip.label, drop = FALSE]

pd_cop_res <- BAT::alpha(comm_cop_binary, tree = tree_cop)
pd_cop     <- as.data.frame(pd_cop_res)
colnames(pd_cop) <- "phylo_diver_cop"

# 4. ASSIGN TO METADATA (USING match TO PRESERVE ROW ORDER)
ecol$phylo_diver_cop <- pd_cop[match(ecol$sample_ID, rownames(pd_cop)), "phylo_diver_cop"]

# Replace NAs and negative values with 0 (samples without Copepoda have PD = 0)
ecol$phylo_diver_cop[is.na(ecol$phylo_diver_cop) | ecol$phylo_diver_cop < 0] <- 0

# 5. MODEL FIT (Tweedie handles positive continuous data with exact zeros)
mod.phylo_cop <- glmmTMB::glmmTMB(
  phylo_diver_cop ~ depth_sc + Mesh + habitat_norm + (1 | ID), 
  data   = ecol, 
  family = tweedie(link = "log")
)

## Check model assumptions
performance::check_overdispersion(mod.phylo_cop)
performance::check_collinearity(mod.phylo_cop)
performance::check_model(mod.phylo_cop)

## Summary and post-hoc
car::Anova(mod.phylo_cop)
summary(mod.phylo_cop)

emmeans::emmeans(mod.phylo_cop, pairwise ~ Mesh, type = "response")
emmeans::emmeans(mod.phylo_cop, pairwise ~ habitat_norm, type = "response")

# Predicted Plot: Mesh
prediction_mesh <- ggpredict(mod.phylo_cop, terms = "Mesh")
ggplot(prediction_mesh, aes(x = x, y = predicted, group = 1)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#012334", alpha=0.8) + 
  labs(x = "Mesh size (µm)", y = "Predicted Phylogenetic Richness", title = "Copepoda: Model-Adjusted Effects") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.line = element_line(color = "black"))

# Predicted Plot: Habitat
prediction_habitat <- ggpredict(mod.phylo_cop, terms = "habitat_norm")
ggplot(prediction_habitat, aes(x = x, y = predicted, group = 1)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#012334") + 
  labs(x = "Type of habitat", y = "Predicted Phylogenetic Richness",title = "Copepoda: Model-Adjusted Effects",subtitle = "" ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(),plot.title = element_text(face = "bold"),axis.line = element_line(color = "black"))


## Nematoda (Alpha Phyl) --------------------------------------------------------

# 1. FILTER Nematoda ASVs PRESENT IN THE RAREFIED MATRIX & CLEAN SEQUENCES
species_nem <- species_meio[
  species_meio$Best.group == "Nematoda" & 
    species_meio$ASVid %in% colnames(comm_meio_pa), 
]

species_nem$sequence <- trimws(as.character(species_nem$sequence))
species_nem <- species_nem[
  !is.na(species_nem$sequence) & 
    nzchar(species_nem$sequence) & 
    species_nem$sequence != "NA", 
]

# 2. FASTA EXPORT & TREE BUILDING
seqs_nem <- Biostrings::DNAStringSet(species_nem$sequence)
names(seqs_nem) <- species_nem$ASVid

alignment_nem   <- DECIPHER::AlignSeqs(seqs_nem)
phy_data_nem    <- phangorn::phyDat(as.matrix(alignment_nem), type = "DNA")
dist_matrix_nem <- phangorn::dist.ml(phy_data_nem)
tree_nem        <- phangorn::NJ(dist_matrix_nem)

# 3. COMMUNITY MATRIX SUBSET & PD CALCULATION
comm_nem_binary <- comm_meio_pa[, colnames(comm_meio_pa) %in% tree_nem$tip.label, drop = FALSE]

pd_nem_res <- BAT::alpha(comm_nem_binary, tree = tree_nem)
pd_nem     <- as.data.frame(pd_nem_res)
colnames(pd_nem) <- "phylo_diver_nem"

# 4. ASSIGN TO METADATA (USING match TO PRESERVE ROW ORDER)
ecol$phylo_diver_nem <- pd_nem[match(ecol$sample_ID, rownames(pd_nem)), "phylo_diver_nem"]

# Replace NAs and negative values with 0 (samples without nemepoda have PD = 0)
ecol$phylo_diver_nem[is.na(ecol$phylo_diver_nem) | ecol$phylo_diver_nem < 0] <- 0

# 5. MODEL FIT (Tweedie handles positive continuous data with exact zeros)
mod.phylo_nem <- glmmTMB::glmmTMB(
  phylo_diver_nem ~ depth_sc + Mesh + habitat_norm + (1 | ID), 
  data   = ecol, 
  family = tweedie(link = "log")
)

## Check model assumptions
performance::check_overdispersion(mod.phylo_nem)
performance::check_collinearity(mod.phylo_nem)
performance::check_model(mod.phylo_nem)

## Summary and post-hoc
car::Anova(mod.phylo_nem)
summary(mod.phylo_nem)

emmeans::emmeans(mod.phylo_nem, pairwise ~ Mesh, type = "response")
emmeans::emmeans(mod.phylo_nem, pairwise ~ habitat_norm, type = "response")


## Box-plot (predicted)

# Mesh
# 1. Extract predictions
prediction_mesh <- ggpredict(mod.phylo_nem, terms = "Mesh")
ggplot(prediction_mesh, aes(x = x, y = predicted, group = 1)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#15adf8", alpha=0.8) + 
  labs(
    x = "Mesh size (µm)",
    y = "Predicted Phylogenetic Richness",
    title = "Nematoda: Model-Adjusted Effects",
    subtitle = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.line = element_line(color = "black")
  )


## Platyhelminthes (Alpha Phyl) ------------------------------------------------

# 1. FILTER  ASVs PRESENT IN THE RAREFIED MATRIX & CLEAN SEQUENCES
species_plat <- species_meio[
  species_meio$Best.group == "Platyhelminthes" & 
    species_meio$ASVid %in% colnames(comm_meio_pa), 
]

species_plat$sequence <- trimws(as.character(species_plat$sequence))
species_plat <- species_plat[
  !is.na(species_plat$sequence) & 
    nzchar(species_plat$sequence) & 
    species_plat$sequence != "NA", 
]

# 2. FASTA EXPORT & TREE BUILDING
seqs_plat <- Biostrings::DNAStringSet(species_plat$sequence)
names(seqs_plat) <- species_plat$ASVid

alignment_plat   <- DECIPHER::AlignSeqs(seqs_plat)
phy_data_plat    <- phangorn::phyDat(as.matrix(alignment_plat), type = "DNA")
dist_matrix_plat <- phangorn::dist.ml(phy_data_plat)
tree_plat        <- phangorn::NJ(dist_matrix_plat)

# 3. COMMUNITY MATRIX SUBSET & PD CALCULATION
comm_plat_binary <- comm_meio_pa[, colnames(comm_meio_pa) %in% tree_plat$tip.label, drop = FALSE]

pd_plat_res <- BAT::alpha(comm_plat_binary, tree = tree_plat)
pd_plat     <- as.data.frame(pd_plat_res)
colnames(pd_plat) <- "phylo_diver_plat"

# 4. ASSIGN TO METADATA (USING match TO PRESERVE ROW ORDER)
ecol$phylo_diver_plat <- pd_plat[match(ecol$sample_ID, rownames(pd_plat)), "phylo_diver_plat"]

# Replace NAs and negative values with 0 (samples without platepoda have PD = 0)
ecol$phylo_diver_plat[is.na(ecol$phylo_diver_plat) | ecol$phylo_diver_plat < 0] <- 0

# 5. MODEL FIT (Tweedie handles positive continuous data with exact zeros)
mod.phylo_plat <- glmmTMB::glmmTMB(
  phylo_diver_plat ~ depth_sc + Mesh + habitat_norm + (1 | ID), 
  data   = ecol, 
  family = tweedie(link = "log")
)

## Check model assumptions
performance::check_overdispersion(mod.phylo_plat)
performance::check_collinearity(mod.phylo_plat)
performance::check_model(mod.phylo_plat)

## Summary and post-hoc
car::Anova(mod.phylo_plat)
summary(mod.phylo_plat)

emmeans::emmeans(mod.phylo_plat, pairwise ~ Mesh, type = "response")
emmeans::emmeans(mod.phylo_plat, pairwise ~ habitat_norm, type = "response")


# Habitat
# 1. Extract predictions
prediction_habitat <- ggpredict(mod.phylo_plat, terms = "habitat_norm")

# 2. "Predicted Values" plot, publication style
ggplot(prediction_habitat, aes(x = x, y = predicted, group = 1)) +
  
  # Confidence intervals ("error bars")
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width = 0.1, linewidth = 0.8, color = "#444444") +
  
  # Predicted points (Orange for Platyhelminthes)
  geom_point(size = 4, color = "#d47000") + 
  
  # Labels and formatting
  labs(
    x = "Type of habitat",
    y = "Predicted Phylogenetic Richness",
    title = "Platyhelminthes: Model-Adjusted Effects",
    subtitle = ""
  ) +
  
  # Visual cleanup
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.line = element_line(color = "black")
  )


## Table with results alpha phyl by group --------------------------------------------

# 1. Put all models in a list with their names
# The name on the left (in quotes) will be used to name the Excel file
model_list <- list(
  "meio" = mod.phylo,
  "cop"  = mod.phylo_cop,
  "nem"  = mod.phylo_nem,
  "plat" = mod.phylo_plat
)

# 2. Create a LOOP that will iterate through each model one by one
for (taxon_name in names(model_list)) {
  
  # Extract the current model for this loop iteration
  current_model <- model_list[[taxon_name]]
  
  cat("Processing model for:", taxon_name, "...\n")
  
  # 0. Extract basic info
  df_info <- data.frame(
    Parameter = c("Total Observations (N)", "Number of Groups (ID)"),
    Value = c(nobs(current_model), summary(current_model)$ngrps$ID)
  )
  
  # 1. Extract coefficients
  df_summary <- as.data.frame(summary(current_model)$coefficients$cond)
  df_summary$Term <- rownames(df_summary) 
  rownames(df_summary) <- NULL
  df_summary <- df_summary[, c("Term", "Estimate", "Std. Error", "z value", "Pr(>|z|)")]
  
  # 2. Extract performance metrics (with safety fallback for singular models)
  
  # A) R2 Calculation
  r2_calc <- suppressWarnings(tryCatch(performance::r2(current_model), error = function(e) NULL))
  
  # Check safely if R2 returned NA, NULL, or an atomic vector
  if (is.null(r2_calc) || !is.list(r2_calc) || is.null(r2_calc$R2_marginal) || is.na(r2_calc$R2_marginal[1])) {
    
    # Fallback: Fit fixed-effects GLM to extract Deviance Explained (D2)
    form_fixed <- update(formula(current_model), . ~ . - (1 | ID))
    glm_fit    <- tryCatch(
      glm(form_fixed, data = ecol, family = family(current_model)$family),
      error = function(e) NULL
    )
    
    if (!is.null(glm_fit)) {
      d2_val <- 1 - (glm_fit$deviance / glm_fit$null.deviance)
      df_r2  <- data.frame(R2_conditional = d2_val, R2_marginal = d2_val)
    } else {
      df_r2  <- data.frame(R2_conditional = NA, R2_marginal = NA)
    }
    
  } else {
    df_r2 <- as.data.frame(r2_calc)
  }
  
  # B) ICC Calculation
  icc_calc <- suppressWarnings(tryCatch(performance::icc(current_model), error = function(e) NULL))
  
  # Check safely for ICC as well
  if (is.null(icc_calc) || !is.list(icc_calc) || is.null(icc_calc$ICC_adjusted) || is.na(icc_calc$ICC_adjusted[1])) {
    df_icc <- data.frame(ICC = 0) # If random effect variance is zero, ICC is 0
  } else {
    df_icc <- as.data.frame(icc_calc)
  }
  
  # C) AIC/BIC
  df_aic <- data.frame(AIC = AIC(current_model), BIC = BIC(current_model))
  
  # Combine metrics
  df_performance <- cbind(df_aic, df_r2, df_icc)
  
  # 3. Extract Type II ANOVA results
  anova_res <- car::Anova(current_model)
  df_anova <- as.data.frame(anova_res)
  df_anova$Term <- rownames(df_anova)
  rownames(df_anova) <- NULL
  df_anova <- df_anova[, c("Term", "Chisq", "Df", "Pr(>Chisq)")]
  
  # 4. Extract post-hoc comparisons
  em_mesh <- emmeans(current_model, pairwise ~ Mesh, type="response")
  em_hab <- emmeans(current_model, pairwise ~ habitat_norm, type="response")
  df_pair_mesh <- as.data.frame(em_mesh$contrasts)
  df_pair_hab <- as.data.frame(em_hab$contrasts)
  
  # 5. Join everything into a list of tabs
  excel_results_list <- list(
    "0_Model_Info" = df_info,
    "1_Model_Summary" = df_summary,
    "2_Performance_Metrics" = df_performance,
    "3_ANOVA_TypeII" = df_anova,
    "4_Pairwise_Mesh" = df_pair_mesh,
    "5_Pairwise_Habitat" = df_pair_hab
  )
  
  # 6. Create file name dynamically and export
  file_name <- paste0("Meiofauna_alpha_phyl_", taxon_name, ".xlsx")
  write_xlsx(excel_results_list, path = file_name)
}



## BETA DIVERSITY-------------------------------------------------------------------------

### Filtration: eliminate ASVs that appear only in one sample
asv_occupancy <- colSums(comm_meio_pa)
comm_meio_beta <- comm_meio_pa[, asv_occupancy > 1, drop = FALSE]

# Delete empty samples
valid_samples_meio <- rowSums(comm_meio_beta) > 0
comm_meio_beta <- comm_meio_beta[valid_samples_meio, , drop = FALSE]

# Sincronize metadata
ecol_meio_beta <- ecol[valid_samples_meio, ]
stopifnot(all(rownames(comm_meio_beta) == ecol_meio_beta$sample_ID))


## A. Taxonomic beta diversity ---------------------------------------------------

## Meiofauna (Beta Tax) ----------------------------------------------------------

beta_diversity <- BAT::beta(comm_meio_beta, func = "jaccard")  

ecol_meio_beta$Mesh <- as.factor(ecol_meio_beta$Mesh)

## PERMANOVA
adonis_meio <- adonis2(beta_diversity$Btotal ~ scale(depth) + Mesh + habitat_norm, 
                       data = ecol_meio_beta, strata = ecol_meio_beta$ID, by = "margin")
print(adonis_meio)

# Turnover and nestedness proportions
# Calculate the global mean of each distance matrix
mean_total <- mean(beta_diversity$Btotal)
mean_turnover <- mean(beta_diversity$Brepl)
mean_nestedness <- mean(beta_diversity$Brich)

# Calculate percentages
per_turnover <- (mean_turnover / mean_total) * 100
per_nestedness <- (mean_nestedness / mean_total) * 100

print(paste("Turnover represents:", round(per_turnover, 2), "%"))
print(paste("Nestedness represents:", round(per_nestedness, 2), "%"))


## Copepoda (Beta Tax) ---------------------------------------------------------

cop_asvs <- species_meio %>%
  filter(Best.group == "Copepoda") %>%
  pull(ASVid)

comm_cop <- comm_meio_beta[, colnames(comm_meio_beta) %in% cop_asvs, drop = FALSE]

# Sincronize samples
valid_samples_cop <- rowSums(comm_cop) > 0
comm_cop <- comm_cop[valid_samples_cop, , drop = FALSE]
ecol_cop <- ecol_meio_beta[valid_samples_cop, ]

stopifnot(all(rownames(comm_cop) == ecol_cop$sample_ID))

# Calculate diversity

beta_diversity_cop <- BAT::beta(comm_cop, func = "jaccard")  

## PERMANOVA
ecol_cop$Mesh <- as.factor(ecol_cop$Mesh)

adonis_cop <- adonis2(beta_diversity_cop$Btotal ~ scale(depth) + Mesh + habitat_norm, 
                      data = ecol_cop, strata = ecol_cop$ID, by = "margin")
print(adonis_cop)

# Turnover and nestedness proportions
mean_total_cop <- mean(beta_diversity_cop$Btotal, na.rm = TRUE)
mean_turnover_cop <- mean(beta_diversity_cop$Brepl, na.rm = TRUE)
mean_nestedness_cop <- mean(beta_diversity_cop$Brich, na.rm = TRUE)

per_turnover_cop <- (mean_turnover_cop / mean_total_cop) * 100
per_nestedness_cop <- (mean_nestedness_cop / mean_total_cop) * 100

print(paste("Copepoda - Turnover represents:", round(per_turnover_cop, 2), "%"))
print(paste("Copepoda - Nestedness represents:", round(per_nestedness_cop, 2), "%"))


## Nematoda (Beta Tax) ---------------------------------------------------------

nem_asvs <- species_meio %>%
  filter(Best.group == "Nematoda") %>%
  pull(ASVid)

comm_nem <- comm_meio_beta[, colnames(comm_meio_beta) %in% nem_asvs, drop = FALSE]

# Sincronize samples
valid_samples_nem <- rowSums(comm_nem) > 0
comm_nem <- comm_nem[valid_samples_nem, , drop = FALSE]
ecol_nem <- ecol_meio_beta[valid_samples_nem, ]

stopifnot(all(rownames(comm_nem) == ecol_nem$sample_ID))

# Calculate diversity

beta_diversity_nem <- BAT::beta(comm_nem, func = "jaccard")  

## PERMANOVA
ecol_nem$Mesh <- as.factor(ecol_nem$Mesh)

adonis_nem <- adonis2(beta_diversity_nem$Btotal ~ scale(depth) + Mesh + habitat_norm, 
                      data = ecol_nem, strata = ecol_nem$ID, by = "margin")
print(adonis_nem)

# Turnover and nestedness proportions
mean_total_nem <- mean(beta_diversity_nem$Btotal, na.rm = TRUE)
mean_turnover_nem <- mean(beta_diversity_nem$Brepl, na.rm = TRUE)
mean_nestedness_nem <- mean(beta_diversity_nem$Brich, na.rm = TRUE)

per_turnover_nem <- (mean_turnover_nem / mean_total_nem) * 100
per_nestedness_nem <- (mean_nestedness_nem / mean_total_nem) * 100

print(paste("Nematoda - Turnover represents:", round(per_turnover_nem, 2), "%"))
print(paste("Nematoda - Nestedness represents:", round(per_nestedness_nem, 2), "%"))


## Tables with results beta tax----------------------------------------------

# function
clean_adonis_df <- function(adonis_res) {
  df <- as.data.frame(adonis_res)
  df$Term <- rownames(df)
  rownames(df) <- NULL
  
  # Select and reorder
  cols_present <- colnames(df)
  target_cols <- c("Term", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")
  cols_to_keep <- intersect(target_cols, cols_present)
  df <- df[, cols_to_keep]
  
  # Rename columns
  colnames(df) <- gsub("SumOfSqs", "Sum_of_Sqs", colnames(df))
  colnames(df) <- gsub("Pr\\(>F\\)", "p_value", colnames(df))
  
  # 4 decimals
  num_cols <- sapply(df, is.numeric)
  df[num_cols] <- lapply(df[num_cols], function(x) round(x, 4))
  
  return(df)
}

# 1. C
excel_beta_list <- list(
  "Total_Meiofauna" = clean_adonis_df(adonis_meio),
  "Copepoda"        = clean_adonis_df(adonis_cop),
  "Nematoda"        = clean_adonis_df(adonis_nem)
)

# 2. Export to an excel file
write_xlsx(excel_beta_list, path = "Results_Adonis_Beta_Tax_ALL.xlsx")


## B. Phylogenetic beta diversity---------------------------------------------------------

### HELPER FUNCTION FOR CLEANING ADONIS OUTPUT FOR EXCEL

clean_adonis_df <- function(adonis_res) {
  df <- as.data.frame(adonis_res)
  df$Term <- rownames(df)
  rownames(df) <- NULL
  
  cols_present <- colnames(df)
  target_cols  <- c("Term", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")
  cols_to_keep <- intersect(target_cols, cols_present)
  df           <- df[, cols_to_keep]
  
  colnames(df) <- gsub("SumOfSqs", "Sum_of_Sqs", colnames(df))
  colnames(df) <- gsub("Pr\\(>F\\)", "p_value", colnames(df))
  
  num_cols <- sapply(df, is.numeric)
  df[num_cols] <- lapply(df[num_cols], function(x) round(x, 4))
  
  return(df)
}


### MASTER FILTER & TREE PRUNING 

# A. Filter base community (ASVs in > 1 sample)
asv_occupancy  <- colSums(comm_meio_pa)
comm_meio_beta <- comm_meio_pa[, asv_occupancy > 1, drop = FALSE]

# B. Sincronice samples
valid_samples_meio <- rowSums(comm_meio_beta) > 0
comm_meio_beta     <- comm_meio_beta[valid_samples_meio, , drop = FALSE]
ecol_meio_beta     <- ecol[valid_samples_meio, ]
ecol_meio_beta$Mesh <- as.factor(ecol_meio_beta$Mesh)

stopifnot(all(rownames(comm_meio_beta) == ecol_meio_beta$sample_ID))

# C. Build alignment and phylogenetic master tree
seqs_meio <- DNAStringSet(species_meio$sequence)
names(seqs_meio) <- species_meio$ASVid
writeXStringSet(seqs_meio, filepath = "Sequences_meioAntarctica18S.fasta", format = "fasta")

alignment_meio   <- AlignSeqs(seqs_meio)
phy_data_meio    <- phyDat(as.matrix(alignment_meio), type = "DNA")
dist_matrix_meio <- dist.ml(phy_data_meio) 
tree_master      <- NJ(dist_matrix_meio)

# D. Keep only ASVs from the filtered matrix
tree_meio_beta <- keep.tip(tree_master, intersect(tree_master$tip.label, colnames(comm_meio_beta)))


## Meiofauna (Beta Phyl) -----------------------------------------------------------

# Calculate beta with filtered tree
beta_phylo_meio <- BAT::beta(comm_meio_beta, tree = tree_meio_beta, func = "jaccard", abund = FALSE)

# PERMANOVA
adonis_meio_phyl <- adonis2(
  beta_phylo_meio$Btotal ~ scale(depth) + Mesh + habitat_norm, 
  data   = ecol_meio_beta, 
  strata = ecol_meio_beta$ID, 
  by     = "margin"
)
print(adonis_meio_phyl)

# Proportions Turnover vs. Nestedness
mean_total_meio      <- mean(as.dist(beta_phylo_meio$Btotal), na.rm = TRUE)
mean_turnover_meio   <- mean(as.dist(beta_phylo_meio$Brepl), na.rm = TRUE)
mean_nestedness_meio <- mean(as.dist(beta_phylo_meio$Brich), na.rm = TRUE)

per_turnover_meio   <- (mean_turnover_meio / mean_total_meio) * 100
per_nestedness_meio <- (mean_nestedness_meio / mean_total_meio) * 100

cat("Total Meiofauna Phylogenetic Turnover:", round(per_turnover_meio, 2), "%\n")
cat("Total Meiofauna Phylogenetic Nestedness:", round(per_nestedness_meio, 2), "%\n")


## Copepoda (Beta Phyl) ---------------------------------------------------

cop_asvs <- species_meio %>%
  filter(Best.group == "Copepoda") %>%
  pull(ASVid)

comm_cop_pa <- comm_meio_beta[, colnames(comm_meio_beta) %in% cop_asvs, drop = FALSE]

# Sincronize valid samples 
valid_samples_cop <- rowSums(comm_cop_pa) > 0
comm_cop_pa       <- comm_cop_pa[valid_samples_cop, , drop = FALSE]
ecol_cop          <- ecol_meio_beta[valid_samples_cop, ]
ecol_cop$Mesh     <- as.factor(ecol_cop$Mesh)

# Filter tree
tree_cop <- keep.tip(tree_meio_beta, intersect(tree_meio_beta$tip.label, colnames(comm_cop_pa)))

# Beta phylogenetic and PERMANOVA
beta_phylo_cop <- BAT::beta(comm_cop_pa, tree = tree_cop, func = "jaccard", abund = FALSE)

adonis_cop_phyl <- adonis2(
  beta_phylo_cop$Btotal ~ scale(depth) + Mesh + habitat_norm, 
  data   = ecol_cop, 
  strata = ecol_cop$ID, 
  by     = "margin"
)
print(adonis_cop_phyl)

# Proporitions
mean_total_cop      <- mean(as.dist(beta_phylo_cop$Btotal), na.rm = TRUE)
mean_turnover_cop   <- mean(as.dist(beta_phylo_cop$Brepl), na.rm = TRUE)
mean_nestedness_cop <- mean(as.dist(beta_phylo_cop$Brich), na.rm = TRUE)

per_turnover_cop   <- (mean_turnover_cop / mean_total_cop) * 100
per_nestedness_cop <- (mean_nestedness_cop / mean_total_cop) * 100

cat("Copepoda Phylogenetic Turnover:", round(per_turnover_cop, 2), "%\n")
cat("Copepoda Phylogenetic Nestedness:", round(per_nestedness_cop, 2), "%\n")


## Nematoda (Beta Phyl) ---------------------------------------------------

nem_asvs <- species_meio %>%
  filter(Best.group == "Nematoda") %>%
  pull(ASVid)

comm_nem_pa <- comm_meio_beta[, colnames(comm_meio_beta) %in% nem_asvs, drop = FALSE]

# Sincronize samples
valid_samples_nem <- rowSums(comm_nem_pa) > 0
comm_nem_pa       <- comm_nem_pa[valid_samples_nem, , drop = FALSE]
ecol_nem          <- ecol_meio_beta[valid_samples_nem, ]
ecol_nem$Mesh     <- as.factor(ecol_nem$Mesh)

# Filter tree
tree_nem <- keep.tip(tree_meio_beta, intersect(tree_meio_beta$tip.label, colnames(comm_nem_pa)))

# Beta phylogenetic and PERMANOVA
beta_phylo_nem <- BAT::beta(comm_nem_pa, tree = tree_nem, func = "jaccard", abund = FALSE)

adonis_nem_phyl <- adonis2(
  beta_phylo_nem$Btotal ~ scale(depth) + Mesh + habitat_norm, 
  data   = ecol_nem, 
  strata = ecol_nem$ID, 
  by     = "margin"
)
print(adonis_nem_phyl)

# Proportions
mean_total_nem      <- mean(as.dist(beta_phylo_nem$Btotal), na.rm = TRUE)
mean_turnover_nem   <- mean(as.dist(beta_phylo_nem$Brepl), na.rm = TRUE)
mean_nestedness_nem <- mean(as.dist(beta_phylo_nem$Brich), na.rm = TRUE)

per_turnover_nem   <- (mean_turnover_nem / mean_total_nem) * 100
per_nestedness_nem <- (mean_nestedness_nem / mean_total_nem) * 100

cat("Nematoda Phylogenetic Turnover:", round(per_turnover_nem, 2), "%\n")
cat("Nematoda Phylogenetic Nestedness:", round(per_nestedness_nem, 2), "%\n")


## Tables with results beta phyl----------------------------------------------

partition_summary_df <- data.frame(
  Taxonomic_Group = c("Total Meiofauna", "Copepoda", "Nematoda"),
  Mean_Phylo_Btotal = round(c(mean_total_meio, mean_total_cop, mean_total_nem), 4),
  Mean_Phylo_Brepl  = round(c(mean_turnover_meio, mean_turnover_cop, mean_turnover_nem), 4),
  Mean_Phylo_Brich  = round(c(mean_nestedness_meio, mean_nestedness_cop, mean_nestedness_nem), 4),
  Pct_Turnover      = round(c(per_turnover_meio, per_turnover_cop, per_turnover_nem), 2),
  Pct_Nestedness    = round(c(per_nestedness_meio, per_nestedness_cop, per_nestedness_nem), 2)
)

excel_phylo_list <- list(
  "Partition_Summary" = partition_summary_df,
  "Total_Meiofauna"   = clean_adonis_df(adonis_meio_phyl),
  "Copepoda"          = clean_adonis_df(adonis_cop_phyl),
  "Nematoda"          = clean_adonis_df(adonis_nem_phyl)
)

write_xlsx(excel_phylo_list, path = "Results_Adonis_Beta_Phylo_ALL_Filtered.xlsx")

cat("File 'Results_Adonis_Beta_Phylo_ALL_Filtered.xlsx' created successfully!\n")


## Figure 1 - SAMPLING MAPS-----------------------------------------------------

# 1. Ross Sea Satellite Zoom
ecol_unique <- ecol %>% distinct(longitude, latitude, .keep_all = TRUE)

register_stadiamaps(key = "INSERT_YOUR_STADIA_MAPS_API_KEY_HERE")
satellite_map <- get_stadiamap(
  bbox = c(left = 163.80, bottom = -74.80, right = 164.30, top = -74.65), 
  zoom = 10, 
  maptype = "stamen_terrain"
)

fig1_zoom <- ggmap(satellite_map) +
  geom_point(data = ecol_unique, aes(x = longitude, y = latitude), color = "red", size = 1.85, alpha = 0.8) +
  geom_text(data = ecol_unique, aes(x = longitude, y = latitude, label = ID), color = "black", size = 2.5, vjust = -1) +
  labs(title = "Ross Sea sampling points", x = "Longitude", y = "Latitude") +
  theme_minimal()

# 2. Antarctic Overview Map
antarctica <- ne_countries(scale = "medium", returnclass = "sf") %>% filter(sovereignt == "Antarctica")

fig1_overview <- ggplot() +
  geom_sf(data = antarctica, fill = "white", color = "black") +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, -60)) +
  theme_minimal() +
  labs(title = "Antarctica (Geographic projection)", x = "Longitude", y = "Latitude")



## Figure 4: Phylogenetic tree ---------------------------------------------------

dataframe2fas <- function(x, file) {
  if (!is.data.frame(x)) {
    x <- as.data.frame(x)
  }
  if (ncol(x) != 2) {
    stop("Input dataframe must contain exactly two columns.")
  }
  fasta_lines <- unlist(
    lapply(seq_len(nrow(x)), function(i) {
      c(
        paste0(">", x[i, 1]),
        as.character(x[i, 2])
      )
    })
  )
  writeLines(fasta_lines, file)
  invisible(fasta_lines)
}

# Export sequences to FASTA using species_meio (FIXED: replaced 'species' with 'species_meio')
dataframe2fas(
  species_meio[, c("ASVid", "sequence")],
  file = "sequences.fasta"
)

# Sequence alignment
seqs <- readDNAStringSet("sequences.fasta")

alignment <- msa(
  seqs,
  method = "Muscle",
  type = "dna"
)

alignment_phydat <- msaConvert(
  alignment,
  type = "phangorn::phyDat"
)

# Initial Neighbor-Joining tree
dist_matrix <- dist.ml(alignment_phydat)
tree_start <- NJ(dist_matrix)
tree_start$tip.label <- names(alignment_phydat)

# Constraint tree
constraint_tree <- read.tree(
  text = "(Xenacoelomorpha,
          ((Nematoda,
          (Tardigrada,
          (Acari,
          (Ostracoda,Copepoda)))),
          ((Gnathostomulida,Rotifera),
          (Gastrotricha,Platyhelminthes,Annelida))));"
)

# Maximum likelihood optimization under topology constraint
fit_start <- pml(tree_start, alignment_phydat)

tree_ml <- optim.pml(
  fit_start,
  optNni = TRUE,
  rearrangement = "stochastic",
  control = pml.control(epsilon = 1e-08),
  constraint = constraint_tree
)

# Plot final tree
plot(
  tree_ml$tree,
  type = "fan",
  show.tip.label = FALSE
)

# Ultrametric tree
tree_ultra <- chronos(tree_ml$tree)

plot(
  tree_ultra,
  type = "fan",
  show.tip.label = FALSE
)

# Save results
save(tree_ml, file = "tree_ml.RData")
save(tree_ultra, file = "tree_ultra.RData")


## Figure 5: Histograms ---------------------------------------------------

comm_target <- comm_meio 

comm_long <- comm_target %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample_ID") %>%
  pivot_longer(-sample_ID, names_to = "ASVid", values_to = "count") %>%
  filter(count > 0)

# Levels and colours
taxon_levels <- c(
  "Xenacoelomorpha", "Priapulida", "Nematoda", "Tardigrada", "Acari", "Ostracoda", "Copepoda",
  "Gnathostomulida", "Rotifera", "Gastrotricha", "Platyhelminthes", "Annelida"
)

taxon_colors <- c(
  "Xenacoelomorpha" = "white",
  "Priapulida"      = "#d4f0fe",
  "Nematoda"        = "#15adf8",
  "Tardigrada"      = "#0697e0",
  "Acari"           = "#057eb9",
  "Ostracoda"       = "#046493",
  "Copepoda"        = "#012334",
  "Gnathostomulida" = "#ffe103",
  "Rotifera"        = "#fbc400",
  "Gastrotricha"    = "#fb8500",
  "Platyhelminthes" = "#d47000",
  "Annelida"        = "#723c00"
)

## A. Mesh Size (FIXED: joins directly with species_meio)
df_mesh <- comm_long %>%
  inner_join(species_meio, by = "ASVid") %>%
  inner_join(ecol, by = "sample_ID") %>%
  filter(!is.na(Mesh) & Mesh %in% c(20, 50, 100, 200)) %>%
  group_by(Mesh, Best.group) %>%
  summarise(ASV_count = n_distinct(ASVid), .groups = "drop") %>%
  mutate(Best.group = factor(Best.group, levels = taxon_levels))

plot_mesh <- ggplot(df_mesh, aes(x = factor(Mesh), y = ASV_count, fill = Best.group)) +
  geom_bar(stat = "identity", position = "fill", color = "black", width = 0.7) +
  scale_fill_manual(values = taxon_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  labs(x = "Mesh size (µm)", y = "Proportion of ASVs", fill = "Taxon", title = "A) Mesh Size") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), legend.position = "right")

## B. Habitat Type (FIXED: joins directly with species_meio)
df_habitat <- comm_long %>%
  inner_join(species_meio, by = "ASVid") %>%
  inner_join(ecol, by = "sample_ID") %>%
  filter(!is.na(habitat_norm)) %>%
  group_by(habitat_norm, Best.group) %>%
  summarise(ASV_count = n_distinct(ASVid), .groups = "drop") %>%
  mutate(
    Best.group   = factor(Best.group, levels = taxon_levels),
    habitat_norm = factor(habitat_norm, levels = c("epilithic", "organic", "spicule", "gravel", "sand", "silt"))
  )

plot_habitat <- ggplot(df_habitat, aes(x = habitat_norm, y = ASV_count, fill = Best.group)) +
  geom_bar(stat = "identity", position = "fill", color = "black", width = 0.7) +
  scale_fill_manual(values = taxon_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  labs(x = "Habitat type", y = "Proportion of ASVs", fill = "Taxon", title = "B) Habitat Type") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

## C. Depth (FIXED: joins directly with species_meio)
df_depth <- comm_long %>%
  inner_join(species_meio, by = "ASVid") %>%
  inner_join(ecol, by = "sample_ID") %>%
  filter(!is.na(depth)) %>%
  mutate(depth_bin = cut(depth, breaks = seq(0, 80, by = 10), include.lowest = TRUE)) %>%
  group_by(depth_bin, Best.group) %>%
  summarise(ASV_count = n_distinct(ASVid), .groups = "drop") %>%
  mutate(Best.group = factor(Best.group, levels = taxon_levels))

plot_depth <- ggplot(df_depth, aes(x = depth_bin, y = ASV_count, fill = Best.group)) +
  geom_bar(stat = "identity", position = "fill", color = "black", width = 0.8) +
  scale_fill_manual(values = taxon_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  labs(x = "Depth range (m)", y = "Proportion of ASVs", fill = "Taxon", title = "C) Depth") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")


plot_mesh
plot_habitat
plot_depth


## C. Metabarcoding papers comparison--------------------------

# Unified Master Color Palette
paper_colors <- c(
  "Ref1"   = "#fb8500", "Ref2"  = "#fbc400", "Ref3"   = "#cb4c07",
  "Ref4"   = "#15adf8", "Ref5"  = "#723c00", "Ref6"   = "#0697e0",
  "Ref7"   = "#057eb9", "Ref6+7"= "#068fdc", "Ref8"   = "#046493",
  "Ref9"   = "#012334", "Ref10" = "#ffe103"
)

## Section 1: MASTER DATA PREPARATION & NOVELTY TABLES (1, 2, & 3)-------

# 1. External Dataset Mapping
datasets_ext <- tibble(
  ref_id = paste0("Ref", 2:10),
  filename = c(
    "Data18S_Asinara_v2.csv", "Data18S_Cordier2021_v2.csv", "Data18S_Degenhardt2021_v2.csv", 
    "Data18S_Fonseca2017.csv", "Data18S_Haenel2017_v2.csv", "Data18S_JondeliusAtherton2020.csv", 
    "Data18S_Kapshyna2024_v2.csv", "Data18S_Mazurkiewicz2024.csv", "Data18S_Polinski2019_v2.csv"
  ),
  dataset_name = c(
    "Asinara", "Cordier2021", "Degenhardt2021", "Fonseca2017", 
    "Haenel2017", "Jondelius2020", "Kapshyna2024", "Mazurkiewicz2024", "Polinski2019"
  )
)

all_global_list  <- list()
all_phylum_list  <- list()
meio_global_list <- list()
meio_phylum_list <- list()
ext_meio_master  <- list()

# 2. Process External Datasets
for (i in seq_len(nrow(datasets_ext))) {
  
  data <- read.csv(here::here("data", datasets_ext$filename[i]), sep = ";", stringsAsFactors = FALSE)
  
  # Extraer vector de similitud (buscando en columna directa o en 'label')
  raw_sim <- if ("similitud" %in% colnames(data)) {
    data$similitud
  } else if ("pident.max" %in% colnames(data)) {
    data$pident.max
  } else if ("pident" %in% colnames(data)) {
    data$pident
  } else if ("label" %in% colnames(data)) {
    str_split_i(data$label, "\\|", 4)
  } else {
    NA
  }
  
  # ASIGNAR Y LIMPIAR COMAS ANTES DEL FILTRO para no romper dimensiones
  data$similitud <- as.numeric(gsub(",", ".", as.character(raw_sim)))
  
  clean_data <- data %>%
    filter(Eval.result == "correct") %>%
    distinct(ASVid, .keep_all = TRUE) %>%
    mutate(
      phylum = if("label" %in% colnames(data)) str_split_i(label, "\\|", 1) else Best.group,
      nreads = as.numeric(nreads)
    )
  
  # A. Summary Lists: All ASVs
  all_global_list[[i]] <- data.frame(
    Dataset     = datasets_ext$dataset_name[i],
    Total_Reads = sum(clean_data$nreads, na.rm = TRUE),
    Total_ASVs  = nrow(clean_data)
  )
  
  all_phylum_list[[i]] <- clean_data %>%
    filter(!is.na(phylum) & phylum != "") %>%
    group_by(phylum) %>%
    summarise(Total_Unique_ASVs = n(), .groups = "drop") %>%
    mutate(Dataset = datasets_ext$dataset_name[i])
  
  # B. Summary Lists: Meiofauna (Incluyendo porcentaje total global <95%)
  meio_data <- clean_data %>%
    filter(Tax.eco.meio == "TRUE" | Tax.eco.meio == TRUE)
  
  meio_global_list[[i]] <- data.frame(
    Dataset      = datasets_ext$dataset_name[i],
    Total_Reads  = sum(meio_data$nreads, na.rm = TRUE),
    Total_ASVs   = nrow(meio_data),
    ASVs_under95 = sum(meio_data$similitud < 95.0, na.rm = TRUE),
    Perc_under95 = mean(meio_data$similitud < 95.0, na.rm = TRUE) * 100
  )
  
  meio_phylum_list[[i]] <- meio_data %>%
    filter(!is.na(phylum) & phylum != "") %>%
    group_by(phylum) %>%
    summarise(
      Total_Unique_ASVs = n(),
      Perc_under_95     = mean(similitud < 95.0, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    mutate(Dataset = datasets_ext$dataset_name[i])
  
  # C. Store raw meiofauna for Master Ref Dataframe
  ext_meio_master[[i]] <- meio_data %>%
    select(ASVid, Best.group = phylum, nreads) %>%
    mutate(dataset = datasets_ext$ref_id[i])
}


# 3. Process Antarctic Data (Ref1: Quality-Filtered & Post-27 Noise Threshold, Unrarefied)

# Apply noise threshold (< 27 reads per sample = 0) on unrarefied comm_raw
asvs_validos_ant <- intersect(colnames(comm_raw), species_tot$ASVid)
comm_ant_27 <- comm_raw[, asvs_validos_ant, drop = FALSE]
comm_ant_27[comm_ant_27 < 27] <- 0

# Extract active ASVs (> 0 reads post-threshold)
ant_active_reads <- colSums(comm_ant_27)
ant_active_ids   <- names(ant_active_reads[ant_active_reads > 0])

# Subset species metadata to active post-27 ASVs
ant_asvs_clean <- species_tot %>%
  filter(ASVid %in% ant_active_ids, !is.na(Best.group))

ant_meio_clean <- ant_asvs_clean %>%
  filter(Tax.eco.meio %in% c("TRUE", TRUE, 1, "permanent"))

# Table Summaries for Antarctica (All ASVs)
ant_global_all <- tibble(
  Dataset     = "Antarctica", 
  Total_Reads = sum(comm_ant_27[, ant_asvs_clean$ASVid]), 
  Total_ASVs  = n_distinct(ant_asvs_clean$ASVid)
)

ant_phylum_all <- ant_asvs_clean %>% 
  group_by(phylum = Best.group) %>% 
  summarise(Total_Unique_ASVs = n_distinct(ASVid), .groups = "drop") %>% 
  mutate(Dataset = "Antarctica")

# Table Summaries for Antarctica (Meiofauna - con novedad total)
comm_ant_meio <- comm_ant_27[, ant_meio_clean$ASVid, drop = FALSE]

ant_global_meio <- tibble(
  Dataset      = "Antarctica", 
  Total_Reads  = sum(comm_ant_meio), 
  Total_ASVs   = n_distinct(ant_meio_clean$ASVid),
  ASVs_under95 = sum(ant_meio_clean$pident.max < 95, na.rm = TRUE),
  Perc_under95 = mean(ant_meio_clean$pident.max < 95, na.rm = TRUE) * 100
)

ant_phylum_meio <- ant_meio_clean %>%
  group_by(phylum = Best.group) %>%
  summarise(
    Total_Unique_ASVs = n_distinct(ASVid),
    Perc_under_95     = mean(pident.max < 95, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(Dataset = "Antarctica")

# Master Antarctic Meiofauna Dataframe (Ref1) for Ref_maestra & iNEXT
ant_meio_reads <- colSums(comm_ant_meio)
ant_meio_df <- ant_meio_clean %>%
  mutate(
    dataset = "Ref1",
    nreads = ant_meio_reads[ASVid]
  ) %>%
  select(dataset, ASVid, Best.group, nreads)


# 4. Construct Single Master Dataframe (Ref_maestra)
Ref_maestra <- bind_rows(ant_meio_df, bind_rows(ext_meio_master))


# 5. Format and Export Tables 1, 2, and 3
all_global_clean  <- bind_rows(c(all_global_list, list(ant_global_all)))
all_phylum_clean  <- bind_rows(c(all_phylum_list, list(ant_phylum_all)))
meio_global_clean <- bind_rows(c(meio_global_list, list(ant_global_meio)))
meio_phylum_clean <- bind_rows(c(meio_phylum_list, list(ant_phylum_meio)))

format_summary_table <- function(global_df, phylum_df) {
  wide_phylum <- phylum_df %>%
    select(Dataset, phylum, Total_Unique_ASVs) %>%
    pivot_wider(names_from = Dataset, values_from = Total_Unique_ASVs, values_fill = 0) %>%
    arrange(phylum)
  
  reads_row <- global_df %>% select(Dataset, Total_Reads) %>% pivot_wider(names_from = Dataset, values_from = Total_Reads) %>% mutate(phylum = "total number of reads")
  asvs_row  <- global_df %>% select(Dataset, Total_ASVs) %>% pivot_wider(names_from = Dataset, values_from = Total_ASVs) %>% mutate(phylum = "total number of ASVs")
  
  bind_rows(reads_row, asvs_row, wide_phylum)
}

table1_all       <- format_summary_table(all_global_clean, all_phylum_clean)
table2_meiofauna <- format_summary_table(meio_global_clean, meio_phylum_clean)

# Table 3: Taxonomic novelty  
overall_novelty_row <- meio_global_clean %>%
  select(Dataset, Perc_under95) %>%
  pivot_wider(names_from = Dataset, values_from = Perc_under95) %>%
  mutate(phylum = "Overall ASVs <95% identity (%)")

phylum_novelty_matrix <- meio_phylum_clean %>%
  select(Dataset, phylum, Perc_under_95) %>%
  pivot_wider(names_from = Dataset, values_from = Perc_under_95, values_fill = 0) %>%
  arrange(phylum)

table3_under95 <- bind_rows(overall_novelty_row, phylum_novelty_matrix)

# 
write_csv(table1_all, "Table_All_ASVs_and_Reads.csv")
write_csv(table2_meiofauna, "TableS2_Meiofauna_ASVs_and_Reads.csv")
write_csv(table3_under95, "TableS7_Phylum_Matrix_Unique_Under95.csv")


# 6. Export Filtered Species Metadata 

species_ant_filtered_all <- species_tot %>%
  filter(ASVid %in% ant_asvs_clean$ASVid) %>%
  mutate(nreads_post27 = ant_active_reads[ASVid])

species_ant_filtered_meio <- species_tot %>%
  filter(ASVid %in% ant_meio_clean$ASVid) %>%
  mutate(nreads_post27 = ant_meio_reads[ASVid])

write_csv(species_ant_filtered_all, "Species_Antarctica_Filtered_TableS9_All.csv")
write_csv(species_ant_filtered_meio, "Species_Antarctica_Filtered_Table_Meiofauna.csv")


## Section 2: NETWORK COMPARISON & SHARED ASVs----------------------------------------


# 1. Geographic Node Preparation
world_map <- map_data("world")

ecol_network <- read.csv(here::here("data","stations network.csv"), sep = ";") %>%
  mutate(
    lat = as.numeric(gsub(",", ".", latitude)),
    lon = as.numeric(gsub(",", ".", longitude)),
    lat = ifelse(regional_location == "Antarctica", -abs(lat), lat),
    studyid = as.character(studyid)
  ) %>% 
  filter(!is.na(lat))

coord_combined <- ecol_network %>%
  filter(studyid %in% c("Ref6", "Ref7")) %>%
  summarise(studyid = "Ref6+7", lon = mean(lon), lat = mean(lat))

ecol_comb <- ecol_network %>%
  filter(!studyid %in% c("Ref6", "Ref7")) %>%
  bind_rows(coord_combined) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

nodos_comb <- ecol_comb %>%
  st_drop_geometry() %>%
  select(id = studyid, lon, lat) %>%
  distinct(id, .keep_all = TRUE) %>%
  mutate(lon_jitter = jitter(lon, amount = 2), lat_jitter = jitter(lat, amount = 2))

edges_base_limpia <- read.csv(here::here("data", "edges_clean.csv"), sep = ";") %>%
  mutate(from = as.character(from), to = as.character(to))

# 2. Global Network Map
graph_geo_comb <- tbl_graph(nodes = nodos_comb, edges = edges_base_limpia %>% uncount(weight), directed = FALSE)

set.seed(42)
plot_geo_map_comb <- ggraph(graph_geo_comb, layout = "manual", x = lon_jitter, y = lat_jitter) +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "#e8e8e8", color = NA) +
  geom_edge_fan(color = "#2c3e50", width = 0.1, alpha = 0.5, spread = 1, show.legend = FALSE) +
  geom_node_point(aes(color = id), size = 4, show.legend = FALSE) +
  geom_node_label(aes(label = id, fill = id), repel = TRUE, size = 3.5, fontface = "bold", color = "black", alpha = 0.8, max.overlaps = Inf, show.legend = FALSE) +
  scale_fill_manual(values = paper_colors) +
  scale_color_manual(values = paper_colors) +
  coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-90, 90)) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=10)))

# 3. Europe Zoom Map
refs_interes <- c("Ref2", "Ref3", "Ref4", "Ref6+7", "Ref8", "Ref9")

edges_eu_comb <- edges_base_limpia %>%
  filter(from %in% refs_interes & to %in% refs_interes) %>%
  uncount(weight)

nodos_eu_comb <- ecol_comb %>%
  st_drop_geometry() %>%
  select(id = studyid, lon, lat) %>%
  distinct(id, .keep_all = TRUE) %>%
  filter(id %in% refs_interes) %>%
  mutate(lon_jitter = jitter(lon, amount = 0.5), lat_jitter = jitter(lat, amount = 0.5))

graph_geo_eu_comb <- tbl_graph(nodes = nodos_eu_comb, edges = edges_eu_comb, directed = FALSE)

set.seed(42)
plot_eu_comb <- ggraph(graph_geo_eu_comb, layout = "manual", x = lon_jitter, y = lat_jitter) +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "#e8e8e8", color = NA) +
  geom_edge_fan(color = "#2c3e50", width = 0.1, alpha = 0.5, spread = 1, show.legend = FALSE) +
  geom_node_point(aes(color = id), size = 4, show.legend = FALSE) +
  geom_node_label(aes(label = id, fill = id), repel = TRUE, size = 3.5, fontface = "bold", color = "black", alpha = 0.8, max.overlaps = Inf, show.legend = FALSE) +
  scale_fill_manual(values = paper_colors) +
  scale_color_manual(values = paper_colors) +
  coord_fixed(ratio = 1.3, xlim = c(-2, 20), ylim = c(36, 70)) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=10)),
    labs(title = "ASV Connectivity: North Sea and Mediterranean", subtitle = "Ref6 and Ref7 combined (Individual ASV threads)")
  )

# 4. Clean Shared ASVs Metadata
df_meta <- read.csv(here::here("data","ASVs_metadata100_gaps.csv"), sep = ";")
df_meta_clean <- df_meta %>%
  group_by(study1, study2, ASVid2_1) %>%
  summarise(
    ASVs_destino_combinados = paste(unique(ASVid2_2), collapse = " | "),
    Numero_de_Matches       = n_distinct(ASVid2_2),
    Grupo_Principal         = paste(unique(ASVid2_1_bestgroup), collapse = " | "),
    Taxonomia_Detallada     = paste(unique(ASVid2_1_sistergroup), collapse = " | "),
    .groups = "drop"
  )

write.csv(df_meta_clean, "Metadata_ASVs_clean.csv", row.names = FALSE)


## Section 3: COVERAGE-BASED RAREFACTION & EXTRAPOLATION (iNEXT)------------------

nombres_legibles <- c(
  "Ref1"  = "Antarctica (Present study)",
  "Ref2"  = "Asinara (Martínez et al. 2020)",
  "Ref3"  = "Cordier et al. (2021)",
  "Ref4"  = "Degenhardt et al. (2021)",
  "Ref5"  = "Fonseca et al. (2017)",
  "Ref6"  = "Haenel et al. (2017)",
  "Ref7"  = "Atherton & Jondelius (2020)",
  "Ref8"  = "Kapshyna et al. (2024)",
  "Ref9"  = "Mazurkiewicz et al. (2024)",
  "Ref10" = "Polinski et al. (2019)"
)

list_studies <- Ref_maestra %>%
  filter(nreads > 0) %>%
  mutate(dataset_name = ifelse(dataset %in% names(nombres_legibles), nombres_legibles[dataset], dataset)) %>%
  split(.$dataset_name) %>%
  lapply(function(df) df$nreads)

# Parallel Processing Setup
num_cores <- max(1, min(length(list_studies), parallel::detectCores() - 2))

# 1. Exact Estimation at C = 0.95 using estimateD()
study_names <- names(list_studies)

cl <- parallel::makeCluster(num_cores)
parallel::clusterEvalQ(cl, library(iNEXT))
parallel::clusterExport(cl, c("list_studies", "study_names"))

estimate_list <- parallel::parLapply(cl, study_names, function(nm) {
  res <- iNEXT::estimateD(
    list_studies[[nm]], 
    q = 0, 
    datatype = "abundance", 
    base = "coverage", 
    level = 0.95, 
    nboot = 100, 
    conf = 0.95
  )
  res$Assemblage <- nm  # Asigna explícitamente el nombre del dataset
  return(res)
})

parallel::stopCluster(cl)

# Consolidate Table S5 directly from estimateD
table_s5_richness <- do.call(rbind, estimate_list) %>%
  as.data.frame() %>%
  select(
    `Study Dataset`              = Assemblage,
    `Coverage Level`             = SC,
    `Sample Size (m)`            = m,
    `Method`                     = Method,
    `Standardized Richness (qD)` = qD,
    `95% CI Lower (LCL)`         = qD.LCL,
    `95% CI Upper (UCL)`         = qD.UCL
  ) %>%
  arrange(desc(`Standardized Richness (qD)`))

write.csv(table_s5_richness, "Table_S5_Standardized_Richness_q0.csv", row.names = FALSE)


## Section 4: FIGURE 2 - TAXONOMIC HEATMAP--------------------------------------

# Prepare Heatmap Data Directly from Ref_maestra (No intermediate CSV re-reading)
df_taxa_counts <- Ref_maestra %>%
  filter(nreads > 0, !is.na(Best.group) & Best.group != "") %>%
  group_by(dataset, taxa = Best.group) %>%
  summarise(ASV_count = n_distinct(ASVid), .groups = "drop")

totals_df <- df_taxa_counts %>%
  group_by(Ref = dataset) %>%
  summarise(Total_ASVs = sum(ASV_count, na.rm = TRUE), .groups = "drop")

df_final <- df_taxa_counts %>%
  rename(Ref = dataset) %>%
  left_join(totals_df, by = "Ref") %>%
  mutate(
    percentage = (ASV_count / Total_ASVs) * 100,
    alpha_val  = na_if(percentage, 0),
    label_text = ifelse(is.na(alpha_val), "", paste0(round(percentage, 1)))
  )

totals_row <- totals_df %>%
  mutate(
    taxa = "Total ASVs",
    ASV_count = Total_ASVs,
    percentage = NA,
    alpha_val = NA,
    label_text = as.character(Total_ASVs)
  )

df_plot <- bind_rows(df_final, totals_row)

orden_refs <- paste0("Ref", 1:10)
df_plot$Ref <- factor(df_plot$Ref, levels = orden_refs)

taxa_levels <- c(sort(unique(df_taxa_counts$taxa)), "Total ASVs")
df_plot$taxa <- factor(df_plot$taxa, levels = taxa_levels)

plot_heatmap <- ggplot(df_plot, aes(x = Ref, y = taxa)) +
  geom_tile(aes(fill = Ref, alpha = alpha_val), color = "white", linewidth = 0.5) +
  geom_text(aes(label = label_text), color = "black", size = 3) +
  scale_fill_manual(values = paper_colors, guide = "none") + 
  scale_alpha_continuous(range = c(0.2, 1), trans = "sqrt", name = "Relative Richness\n(% of ASVs)", na.value = 0) +
  scale_y_discrete(limits = rev) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = paper_colors[levels(df_plot$Ref)]),
    axis.text.y = element_text(size = 10, face = ifelse(rev(levels(df_plot$taxa)) == "Total ASVs", "bold", "plain")),
    panel.grid = element_blank(), 
    panel.background = element_rect(fill = "#fdfdfd", color = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b=15))
  ) +
  labs(title = "ASV Composition by Taxonomic Group", x = "Reference", y = "Taxonomic Group")

print(plot_heatmap)



## Session info---------------------------------------------------------------

# R version 4.3.3 (2024-02-29 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 11 x64 (build 26200)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=Spanish_Spain.utf8  LC_CTYPE=Spanish_Spain.utf8    LC_MONETARY=Spanish_Spain.utf8
# [4] LC_NUMERIC=C                   LC_TIME=Spanish_Spain.utf8    
# 
# time zone: Europe/Madrid
# tzcode source: internal
# 
# attached base packages:
# [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] msa_1.34.0          writexl_1.5.3       pheatmap_1.0.12     ggraph_2.2.1        tidygraph_1.3.1    
#  [6] officer_0.7.3       flextable_0.9.10    ggmap_4.0.2         marmap_1.0.12       rnaturalearth_1.0.1
# [11] sf_1.0-19           phangorn_2.12.1     DECIPHER_2.30.0     RSQLite_2.4.3       Biostrings_2.70.3  
# [16] GenomeInfoDb_1.38.8 XVector_0.42.0      IRanges_2.36.0      S4Vectors_0.40.2    BiocGenerics_0.48.1
# [21] phytools_2.4-4      maps_3.4.2.1        ape_5.8-1           iNEXT_3.0.2         ggeffects_2.3.2    
# [26] emmeans_1.10.7      performance_0.15.1  glmmTMB_1.1.10      BAT_2.11.0          vegan_2.6-10       
# [31] lattice_0.22-5      permute_0.9-7       ggrepel_0.9.6       scales_1.4.0        sjPlot_2.9.0       
# [36] patchwork_1.3.2     lubridate_1.9.4     forcats_1.0.0       stringr_1.5.1       dplyr_1.1.4        
# [41] purrr_1.0.4         readr_2.1.5         tidyr_1.3.1         tibble_3.2.1        ggplot2_4.0.0      
# [46] tidyverse_2.0.0     here_1.0.2         
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-9            httr_1.4.7              RColorBrewer_1.1-3      insight_1.4.2          
#   [5] doParallel_1.0.17       numDeriv_2016.8-1.1     tools_4.3.3             sjlabelled_1.2.0       
#   [9] R6_2.6.1                mgcv_1.9-1              withr_3.0.2             sp_2.2-0               
#  [13] gridExtra_2.3           prettyunits_1.2.0       textshaping_1.0.0       cli_3.6.4              
#  [17] labeling_0.4.3          mvtnorm_1.3-3           S7_0.2.0                proxy_0.4-27           
#  [21] askpass_1.2.1           systemfonts_1.3.1       R.utils_2.13.0          DHARMa_0.4.7           
#  [25] DEoptim_2.2-8           parallelly_1.45.1       pdist_1.2.1             rstudioapi_0.17.1      
#  [29] optimParallel_1.0-2     generics_0.1.3          shape_1.4.6.1           combinat_0.0-8         
#  [33] vroom_1.6.6             zip_2.3.2               car_3.1-3               Matrix_1.6-5           
#  [37] abind_1.4-8             R.methodsS3_1.8.2       terra_1.8-21            lifecycle_1.0.4        
#  [41] scatterplot3d_0.3-44    carData_3.0-5           clusterGeneration_1.3.8 recipes_1.3.1          
#  [45] grid_4.3.3              blob_1.2.4              PlotTools_0.3.1         crayon_1.5.3           
#  [49] haven_2.5.5             pillar_1.10.1           knitr_1.49              boot_1.3-29            
#  [53] estimability_1.5.1      future.apply_1.20.0     codetools_0.2-19        fastmatch_1.1-6        
#  [57] glue_1.8.0              fontLiberation_0.1.0    data.table_1.17.0       vctrs_0.6.5            
#  [61] png_0.1-8               Rdpack_2.6.2            gtable_0.3.6            datawizard_1.2.0       
#  [65] cachem_1.1.0            ks_1.14.3               gower_1.0.2             xfun_0.51              
#  [69] rbibutils_2.3           prodlim_2024.06.25      pracma_2.4.4            coda_0.19-4.1          
#  [73] reformulas_0.4.0        survival_3.5-8          ncdf4_1.24              geometry_0.5.2         
#  [77] timeDate_4041.110       iterators_1.0.14        hardhat_1.4.2           units_0.8-5            
#  [81] lava_1.8.1              ipred_0.9-15            nlme_3.1-164            bit64_4.6.0-1          
#  [85] fontquiver_0.2.1        progress_1.2.3          rprojroot_2.1.1         R.cache_0.17.0         
#  [89] TMB_1.9.16              TreeTools_1.13.1        KernSmooth_2.23-22      rpart_4.1.23           
#  [93] colorspace_2.1-1        DBI_1.2.3               raster_3.6-32           nnet_7.3-19            
#  [97] mnormt_2.1.1            tidyselect_1.2.1        curl_6.2.1              bit_4.5.0.1            
# [101] compiler_4.3.3          see_0.12.0              expm_1.0-0              xml2_1.4.0             
# [105] fontBitstreamVera_0.1.1 bayestestR_0.17.0       classInt_0.4-11         quadprog_1.5-8         
# [109] palmerpenguins_0.1.1    digest_0.6.37           minqa_1.2.8             rmarkdown_2.29         
# [113] htmltools_0.5.8.1       pkgconfig_2.0.3         jpeg_0.1-10             lme4_1.1-36            
# [117] fastmap_1.2.0           rlang_1.1.5             farver_2.1.2            adehabitatMA_0.3.17    
# [121] jsonlite_1.9.0          mclust_6.1.1            ModelMetrics_1.2.2.2    R.oo_1.27.1            
# [125] nls2_0.3-4              RCurl_1.98-1.17         magrittr_2.0.3          Formula_1.2-5          
# [129] GenomeInfoDbData_1.2.11 gdistance_1.6.5         Rcpp_1.0.14             viridis_0.6.5          
# [133] proto_1.0.0             gdtools_0.4.4           stringi_1.8.4           pROC_1.18.5            
# [137] zlibbioc_1.48.2         MASS_7.3-60.0.1         plyr_1.8.9              hypervolume_3.1.5      
# [141] listenv_0.9.1           graphlayouts_1.2.2      splines_4.3.3           hms_1.1.3              
# [145] uuid_1.2-1              igraph_2.1.4            fastcluster_1.2.6       reshape2_1.4.4         
# [149] magic_1.6-1             evaluate_1.0.3          tweenr_2.0.3            nloptr_2.1.1           
# [153] tzdb_0.5.0              foreach_1.5.2           polyclip_1.10-7         openssl_2.3.2          
# [157] future_1.67.0           ggforce_0.5.0           xtable_1.8-4            e1071_1.7-16           
# [161] viridisLite_0.4.2       ragg_1.3.3              class_7.3-22            memoise_2.0.1          
# [165] cluster_2.1.6           timechange_0.3.0        globals_0.18.0          caret_7.0-1  
# ==============================================================================
