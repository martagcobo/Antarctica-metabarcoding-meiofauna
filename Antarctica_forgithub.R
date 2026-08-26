#################################
# Antarctic meiofauna metabarcoding paper
#
# Script by: anonymus until peer-review is finished
# Last update: 26/08/2026

###########################


##### 1. Set working directory -------------------------------------------------
# Set to your local folder (Update this path as needed)
setwd("C:")

##### 2. Load packages ---------------------------------------------------------
# Data manipulation and visualization
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(sjPlot)
library(tibble)
library(scales)
library(VennDiagram)
library(grid)

# Community ecology and statistics
library(vegan)
library(BAT)
library(glmmTMB)
library(performance)
library(lme4)
library(car)
library(emmeans)
library(pairwiseAdonis)
library(ggeffects)

# Phylogeny and genetics
library(ape)
library(phytools)
library(Biostrings) # For the function readDNAStringSet
library(DECIPHER)   # For the function AlignSeqs
library(phangorn)   # For dist.ml
library(msa)

# Spatial and mapping
library(sf)
library(rnaturalearth)
library(marmap)
library(ggmap)

# Tables and networks
library(flextable)
library(officer)
library(rlang)
library(ggrepel)   # For paper labels
library(tidygraph) # To handle networks
library(ggraph)    # To draw networks
library(pheatmap)
library(writexl)   # ADDED: Required for write_xlsx

library(iNEXT) # rarefaction
library(parallel) # for computational effort


#### Data Preparation -------------------------------------

# 1. LOAD RAW DATA 
ecol_base    <- read.csv2("Antarctica marine MBC samples.csv")
species_base <- read.csv2("Data18S_Antarctica_v2.csv")

# Load community matrix with row.names = 1 to keep ALL raw counts untouched
comm_raw_base <- read.csv2("Antarctica_community.csv", row.names = 1)

# Remove non-target rows if applicable
comm_raw <- comm_raw_base[-(211:218), ]
comm_raw <- comm_raw[-(81:88), ]

# Ensure all columns are numeric ASV counts
comm_raw <- comm_raw[, sapply(comm_raw, is.numeric)]


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


# 3. RAREFY THE COMPLETE RAW COMMUNITY
set.seed(123) # For reproducibility

min_depth_tot <- min(rowSums(comm_raw))
cat("Rarefying total raw community to:", min_depth_tot, "reads\n")

# Global rarefaction
comm_tot_rar <- rrarefy(comm_raw, sample = min_depth_tot)


# 4. PREPARE SPECIES FILTERS AND APPLY POST-RAREFACTION NOISE FILTER
species_clean <- species_base[
  !is.na(species_base$ASVid) &
    !is.na(species_base$sequence) &
    nzchar(species_base$sequence), 
]

species_tot <- species_clean[which(species_clean$Eval.result == "correct"), ]
species_tot <- species_tot[which(species_tot$rare.asvs == "FALSE"), ]
species_tot <- species_tot[!duplicated(species_tot$ASVid), ]

species_meio <- species_tot[which(species_tot$Tax.eco.meio %in% c("TRUE", TRUE, 1, "permanent")), ]

# Filter correct & non-rare ASVs
comm_tot_valid <- comm_tot_rar[, colnames(comm_tot_rar) %in% species_tot$ASVid]

# Apply noise threshold (< 27 reads per sample = 0)
comm_tot_valid[comm_tot_valid < 27] <- 0

# Subset to Meiofauna ASVs
comm_meio <- comm_tot_valid[, colnames(comm_tot_valid) %in% species_meio$ASVid]


# 5. CALCULATE FINAL METRICS
ecol$total_reads_raw   <- rowSums(comm_raw)
ecol$meio_reads_final  <- rowSums(comm_meio)
ecol$richness_meio_rar <- rowSums(comm_meio > 0)

# Presence/Absence Matrix
comm_meio_pa <- comm_meio
comm_meio_pa[comm_meio_pa > 0] <- 1

# Factor formatting
habitat_order     <- c("epilithic", "organic", "spicule", "gravel", "sand", "silt")
ecol$habitat_norm <- factor(ecol$habitat_norm, levels = habitat_order)
ecol$Mesh         <- as.factor(ecol$Mesh)



## Sensitivity analysis------------------------------------------

# 1. Do two independent rarefactions
set.seed(123)
rar_1 <- rrarefy(comm_raw, sample = min(rowSums(comm_raw)))
set.seed(456)
rar_2 <- rrarefy(comm_raw, sample = min(rowSums(comm_raw)))

# 2. Convert to presence/ absence and calculate Jaccard distance
d1 <- vegdist((rar_1 > 0)*1, method = "jaccard")
d2 <- vegdist((rar_2 > 0)*1, method = "jaccard")

# 3. Mantel test
mantel(d1, d2)



################## ALPHA DIVERSITY #############################---------------------
############# A. Taxonomic Alpha Diversity #####################

#### Meiofauna (Alpha Tax)-----------------------------------------
model <- glmmTMB(richness_meio_rar ~ scale(depth) + Mesh + habitat_norm + (1 | ID), data = ecol, family = poisson)

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


#### Copepoda (Alpha Tax) ------------------------------------------------

# Identify ASVs
cop_asvs <- species_meio %>%
  filter(Best.group == "Copepoda") %>%
  pull(ASVid)

# Subset meiofauna presence/absence matrix (drop = FALSE keeps it as a matrix)
comm_cop <- comm_meio_pa[, colnames(comm_meio_pa) %in% cop_asvs, drop = FALSE]

# Calculate richness per sample
ecol$richness_cop <- rowSums(comm_cop > 0)

# Model
model_cop <- glmmTMB(richness_cop ~ scale(depth) + Mesh + habitat_norm + (1 | ID), data = ecol, family = nbinom2) 

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



#### Nematoda (Alpha Tax) ------------------------------------------------------

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


#### Platyhelminthes (Alpha Tax) -----------------------------------------------

# Identify ASVs
plat_asvs <- species_meio %>%
  filter(Best.group == "Platyhelminthes") %>%
  pull(ASVid)

# Subset meiofauna presence/absence matrix (drop = FALSE keeps it as a matrix)
comm_plat <- comm_meio_pa[, colnames(comm_meio_pa) %in% plat_asvs, drop = FALSE]

# Calculate richness per sample
ecol$richness_plat <- rowSums(comm_plat > 0)

# Model
model_plat <- glmmTMB(richness_plat ~ scale(depth) + Mesh + habitat_norm +  (1 | ID), data = ecol, family = nbinom2) 

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



############# B. Phylogenetic Alpha Diversity ####################################---------------------------------------

## Meiofauna (Alpha Phyl)----------------------------------------
# Clean and test sequences (using the Meiofauna subset)
species_meio$sequence <- trimws(as.character(species_meio$sequence))

species_meio <- species_meio[
  !is.na(species_meio$sequence) &
    nzchar(species_meio$sequence) &
    species_meio$sequence != "NA",
]

# Verify sequence cleaning (should all return FALSE)
any(is.na(species_meio$sequence))
any(species_meio$sequence == "")
any(grepl("^\\s+$", species_meio$sequence))
any(species_meio$sequence == "NA")

# Export and read FASTA file
seqs <- DNAStringSet(species_meio$sequence)
names(seqs) <- species_meio$ASVid
writeXStringSet(seqs, filepath = "Sequences_meioAntarctica18S.fasta", format = "fasta")

seqs <- readDNAStringSet("Sequences_meioAntarctica18S.fasta")

# Tree building
alignment <- AlignSeqs(seqs)

# Convert alignment to a phyDat object
phy_data <- phyDat(as.matrix(alignment), type = "DNA")

# Make phylogenetic tree with Neighbor-Joining
dist_matrix <- dist.ml(phy_data)  # Distance matrix
tree <- NJ(dist_matrix)           # Tree with Neighbor-Joining

# Calculate Phylogenetic Diversity (using rarefied meiofauna matrix)
# Filter matrix to only include ASVs present in the tree
comm_meio_pa_phylo <- comm_meio_pa[, colnames(comm_meio_pa) %in% tree$tip.label, drop = FALSE]

# Calculate alpha phylogenetic diversity (BAT package)
pd_sample <- BAT::alpha(comm_meio_pa_phylo, tree = tree)
pd_sample <- as.data.frame(pd_sample)
colnames(pd_sample) <- "phylo.diver"
pd_sample$sample_ID <- rownames(pd_sample)

# Merge PD results directly into the main 'ecol' dataframe
ecol <- merge(ecol, pd_sample, by = "sample_ID", all.x = TRUE)

# Replace NAs with 0 (samples with 0 Meiofauna ASVs have 0 Phylogenetic Diversity)
ecol$phylo.diver[is.na(ecol$phylo.diver)] <- 0

# Transform negative values into 0
ecol$phylo.diver[ecol$phylo.diver < 0] <- 0


# Model
mod.phylo <- glmmTMB(phylo.diver ~ scale(depth) + Mesh + habitat_norm + (1 | ID), ziformula = ~ 1,
  data   = ecol, family = tweedie(link = "log"))

## Check the model
performance::check_overdispersion(mod.phylo)
performance::check_collinearity(mod.phylo)
# pdf("C:/Users/Lab 22/OneDrive - Universidad Complutense de Madrid (UCM)/Mi unidad/Proyectos/Emilia-Romagna/check_model_2.pdf", width = 10, height = 7) 
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


## Copepoda (Alpha Phyl) ------------------------------------------------

# Filter sequences for Copepoda
species_cop <- species_meio %>% 
  dplyr::filter(Best.group == "Copepoda")

seqs_cop <- DNAStringSet(species_cop$sequence)
names(seqs_cop) <- species_cop$ASVid
writeXStringSet(seqs_cop, filepath = "Sequences_copAntarctica18S.fasta", format = "fasta")

# Alignment and NJ tree
alignment_cop <- AlignSeqs(seqs_cop)
phy_data_cop  <- phyDat(as.matrix(alignment_cop), type = "DNA")
dist_matrix_cop <- dist.ml(phy_data_cop)
tree_cop      <- NJ(dist_matrix_cop)

# Community Matrix (Presence/Absence)
asvs_cop <- species_cop$ASVid
# Use drop = FALSE to ensure it stays a matrix even if there is only 1 ASV
comm_cop_binary <- comm_meio_pa[, colnames(comm_meio_pa) %in% asvs_cop, drop = FALSE]

# Phylogenetic richness (PD)
pd_cop_res <- BAT::alpha(comm_cop_binary, tree = tree_cop)
pd_cop <- as.data.frame(pd_cop_res)

# Name it specifically for Copepoda to avoid overwriting total meiofauna PD
colnames(pd_cop) <- "phylo_diver_cop"
pd_cop$sample_ID <- rownames(pd_cop)

# Merge with the main 'ecol' dataframe
ecol <- merge(ecol, pd_cop, by = "sample_ID", all.x = TRUE)

# Substitute NAs with 0 (samples without Copepoda have PD = 0)
ecol$phylo_diver_cop[is.na(ecol$phylo_diver_cop)] <- 0

# Transform negative values into 0
ecol$phylo_diver_cop[ecol$phylo_diver_cop < 0] <- 0

# Phylogenetic richness model
# Note: Using ziformula = ~ 1 handles excess structural zeros perfectly
mod.phylo_cop <- glmmTMB(phylo_diver_cop ~ scale(depth) + Mesh + habitat_norm + (1 | ID), 
  data      = ecol, family    = tweedie(link = "log"))

## Check the model
performance::check_overdispersion(mod.phylo_cop)
performance::check_collinearity(mod.phylo_cop)
# pdf(...)
performance::check_model(mod.phylo_cop)
# dev.off()

## Summary and post-hoc
car::Anova(mod.phylo_cop)
summary(mod.phylo_cop)

## Post-hoc
emmeans(mod.phylo_cop, pairwise ~ Mesh, type="response")
emmeans(mod.phylo_cop, pairwise ~ habitat_norm, type="response")


# Predicted Plot: Mesh
prediction_mesh <- ggpredict(mod.phylo_cop, terms = "Mesh")
ggplot(prediction_mesh, aes(x = x, y = predicted, group = 1)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#012334", alpha=0.8) + 
  labs(x = "Mesh size (µm)", y = "Predicted Phylogenetic Richness", title = "Copepoda: Model-Adjusted Effects") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.line = element_line(color = "black"))



## Nematoda (Alpha Phyl) -------------------------------------------------------

# Filter sequences for Nematoda
species_nem <- species_meio %>% 
  dplyr::filter(Best.group == "Nematoda")

seqs_nem <- DNAStringSet(species_nem$sequence)
names(seqs_nem) <- species_nem$ASVid
writeXStringSet(seqs_nem, filepath = "Sequences_nemAntarctica18S.fasta", format = "fasta")

# Alignment and NJ tree
alignment_nem <- AlignSeqs(seqs_nem)
phy_data_nem  <- phyDat(as.matrix(alignment_nem), type = "DNA")
dist_matrix_nem <- dist.ml(phy_data_nem)
tree_nem      <- NJ(dist_matrix_nem)

# Community Matrix (Presence/Absence)
asvs_nem <- species_nem$ASVid
# Use drop = FALSE to ensure it stays a matrix even if there is only 1 ASV
comm_nem_binary <- comm_meio_pa[, colnames(comm_meio_pa) %in% asvs_nem, drop = FALSE]

# Phylogenetic richness (PD)
pd_nem_res <- BAT::alpha(comm_nem_binary, tree = tree_nem)
pd_nem <- as.data.frame(pd_nem_res)

# Name it specifically for Nematoda to avoid overwriting total meiofauna PD
colnames(pd_nem) <- "phylo_diver_nem"
pd_nem$sample_ID <- rownames(pd_nem)

# Merge with the main 'ecol' dataframe
ecol <- merge(ecol, pd_nem, by = "sample_ID", all.x = TRUE)

# Substitute NAs with 0 (samples without Nematoda have PD = 0)
ecol$phylo_diver_nem[is.na(ecol$phylo_diver_nem)] <- 0

# Transform negative values into 0
ecol$phylo_diver_nem[ecol$phylo_diver_nem < 0] <- 0

# Phylogenetic richness model
mod.phylo_nem <- glmmTMB(phylo_diver_nem ~ scale(depth) + Mesh + habitat_norm + (1 | ID), 
                         data      = ecol, family    = tweedie(link = "log"))

## Check the model
performance::check_overdispersion(mod.phylo_nem)
performance::check_collinearity(mod.phylo_nem)

# pdf(...)
performance::check_model(mod.phylo_nem)

# dev.off()

## Summary and post-hoc
car::Anova(mod.phylo_nem)
summary(mod.phylo_nem)

## Post-hoc
emmeans(mod.phylo_nem, pairwise ~ Mesh, type="response")
emmeans(mod.phylo_nem, pairwise ~ habitat_norm, type="response")

## Box-plot (predicted)

# Mesh
# 1. Extract predictions
prediction_mesh <- ggpredict(mod.phylo_nem, terms = "Mesh")

# 2. "Predicted Values" plot, publication style
ggplot(prediction_mesh, aes(x = x, y = predicted, group = 1)) +
  
  # Confidence intervals ("error bars")
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width = 0.1, linewidth = 0.8, color = "#444444") +
  
  # Predicted points (Blue for Nematoda)
  geom_point(size = 4, color = "#15adf8", alpha=0.8) + 
  
  # Labels and formatting
  labs(
    x = "Mesh size (µm)",
    y = "Predicted Phylogenetic Richness",
    title = "Nematoda: Model-Adjusted Effects",
    subtitle = ""
  ) +
  
  # Visual cleanup
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.line = element_line(color = "black")
  )


## Platyhelminthes (Alpha Phyl) ------------------------------------------------

# Filter sequences for Platyhelminthes
species_plat <- species_meio %>% 
  dplyr::filter(Best.group == "Platyhelminthes")

seqs_plat <- DNAStringSet(species_plat$sequence)
names(seqs_plat) <- species_plat$ASVid
writeXStringSet(seqs_plat, filepath = "Sequences_platAntarctica18S.fasta", format = "fasta")

# Alignment and NJ tree
alignment_plat <- AlignSeqs(seqs_plat)
phy_data_plat  <- phyDat(as.matrix(alignment_plat), type = "DNA")
dist_matrix_plat <- dist.ml(phy_data_plat)
tree_plat      <- NJ(dist_matrix_plat)

# Community Matrix (Presence/Absence)
asvs_plat <- species_plat$ASVid
# Use drop = FALSE to ensure it stays a matrix even if there is only 1 ASV
comm_plat_binary <- comm_meio_pa[, colnames(comm_meio_pa) %in% asvs_plat, drop = FALSE]

# Phylogenetic richness (PD)
pd_plat_res <- BAT::alpha(comm_plat_binary, tree = tree_plat)
pd_plat <- as.data.frame(pd_plat_res)

# Name it specifically for Platyhelminthes to avoid overwriting total meiofauna PD
colnames(pd_plat) <- "phylo_diver_plat"
pd_plat$sample_ID <- rownames(pd_plat)

# Merge with the main 'ecol' dataframe
ecol <- merge(ecol, pd_plat, by = "sample_ID", all.x = TRUE)

# Substitute NAs with 0 (samples without Platyhelminthes have PD = 0)
ecol$phylo_diver_plat[is.na(ecol$phylo_diver_plat)] <- 0

# Transform negative values into 0
ecol$phylo_diver_plat[ecol$phylo_diver_plat < 0] <- 0

# Phylogenetic richness model
# Note: Using ziformula = ~ 1 handles excess structural zeros perfectly
mod.phylo_plat <- glmmTMB(phylo_diver_plat ~ scale(depth) + Mesh + habitat_norm + (1 | ID), 
                        data      = ecol, family    = tweedie(link = "log"))
## Check the model
performance::check_overdispersion(mod.phylo_plat)
performance::check_collinearity(mod.phylo_plat)

# pdf(...)
performance::check_model(mod.phylo_plat)
# dev.off()

## Summary and post-hoc
car::Anova(mod.phylo_plat)
summary(mod.phylo_plat)

## Post-hoc
emmeans(mod.phylo_plat, pairwise ~ Mesh, type="response")
emmeans(mod.phylo_plat, pairwise ~ habitat_norm, type="response")

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



################### BETA DIVERSITY #############################################

### Filtration: eliminate ASVs that appear only in one sample
asv_occupancy <- colSums(comm_meio_pa)
comm_meio_beta <- comm_meio_pa[, asv_occupancy > 1, drop = FALSE]

# Delete empty samples
valid_samples_meio <- rowSums(comm_meio_beta) > 0
comm_meio_beta <- comm_meio_beta[valid_samples_meio, , drop = FALSE]

# Sincronize metadata
ecol_meio_beta <- ecol[valid_samples_meio, ]
stopifnot(all(rownames(comm_meio_beta) == ecol_meio_beta$sample_ID))


############## A. Taxonomic beta diversity ---------------------------------------

## Meiofauna (Beta Tax) ---------------------------------------------

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


############# B. Phylogenetic beta diversity #####################################

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


## Meiofauna (Beta Phyl) ----------------------------------------

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
  Taxonomic_Group = c("Total Meiofauna", "Copepoda", "Nematoda", "Platyhelminthes"),
  Mean_Phylo_Btotal = round(c(mean_total_meio, mean_total_cop, mean_total_nem, mean_total_plat), 4),
  Mean_Phylo_Brepl  = round(c(mean_turnover_meio, mean_turnover_cop, mean_turnover_nem, mean_turnover_plat), 4),
  Mean_Phylo_Brich  = round(c(mean_nestedness_meio, mean_nestedness_cop, mean_nestedness_nem, mean_nestedness_plat), 4),
  Pct_Turnover      = round(c(per_turnover_meio, per_turnover_cop, per_turnover_nem, per_turnover_plat), 2),
  Pct_Nestedness    = round(c(per_nestedness_plat, per_nestedness_cop, per_nestedness_nem, per_nestedness_plat), 2)
)

excel_phylo_list <- list(
  "Partition_Summary" = partition_summary_df,
  "Total_Meiofauna"   = clean_adonis_df(adonis_meio_phyl),
  "Copepoda"          = clean_adonis_df(adonis_cop_phyl),
  "Nematoda"          = clean_adonis_df(adonis_nem_phyl),
  "Platyhelminthes"   = clean_adonis_df(adonis_plat_phyl)
)

write_xlsx(excel_phylo_list, path = "Results_Adonis_Beta_Phylo_ALL_Filtered.xlsx")

cat("File 'Results_Adonis_Beta_Phylo_ALL_Filtered.xlsx' created successfully!\n")



### Figure 4: Phylogenetic tree ---------------------------------------------------

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


### Figure 5: Histograms ---------------------------------------------------

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


### Figure 6: Venn diagram ----------------------------------------------------------

comm_venn <- comm_meio_pa %>%
  as.data.frame() %>%
  mutate(Mesh = factor(ecol$Mesh, levels = c("20", "50", "100", "200")))

asv_sets <- comm_venn %>%
  group_by(Mesh) %>%
  summarise(across(where(is.numeric), ~ any(. > 0)), .groups = "drop") %>%
  as.data.frame()

set_list <- lapply(split(asv_sets[, -1], asv_sets$Mesh), function(x) {
  colnames(x)[which(as.logical(x))]
})

set_list <- set_list[c("20", "50", "100", "200")]

mesh_colors <- c(
  "20"  = "#cce4f6",
  "50"  = "#66b3e6",
  "100" = "#1f78b4",
  "200" = "#012334"
)

grid.newpage()

venn.plot <- venn.diagram(
  x = set_list,
  filename = NULL,
  category.names = names(set_list),
  fill = mesh_colors[names(set_list)],
  alpha = 0.5,
  cex = 1.3,
  cat.cex = 1.2,
  cat.pos = 0
)

grid.draw(venn.plot)


#### C. Metabarcoding papers comparison--------------------------

# Unified Master Color Palette
paper_colors <- c(
  "Ref1"   = "#fb8500", "Ref2"  = "#fbc400", "Ref3"   = "#cb4c07",
  "Ref4"   = "#15adf8", "Ref5"  = "#723c00", "Ref6"   = "#0697e0",
  "Ref7"   = "#057eb9", "Ref6+7"= "#068fdc", "Ref8"   = "#046493",
  "Ref9"   = "#012334", "Ref10" = "#ffe103"
)

# Section 1: MASTER DATA PREPARATION & NOVELTY TABLES (1, 2, & 3)-------

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
  
  data <- read.csv(datasets_ext$filename[i], sep = ";", stringsAsFactors = FALSE)
  
  clean_data <- data %>%
    filter(Eval.result == "correct") %>%
    distinct(ASVid, .keep_all = TRUE) %>%
    mutate(
      similitud = as.numeric(str_split_i(label, "\\|", 4)),
      phylum    = str_split_i(label, "\\|", 1),
      nreads    = as.numeric(nreads)
    )
  
  # A. Summary Lists: All ASVs
  all_global_list[[i]] <- data.frame(
    Dataset = datasets_ext$dataset_name[i],
    Total_Reads = sum(clean_data$nreads, na.rm = TRUE),
    Total_ASVs = nrow(clean_data)
  )
  
  all_phylum_list[[i]] <- clean_data %>%
    filter(!is.na(phylum) & phylum != "") %>%
    group_by(phylum) %>%
    summarise(Total_Unique_ASVs = n(), .groups = "drop") %>%
    mutate(Dataset = datasets_ext$dataset_name[i])
  
  # B. Summary Lists: Meiofauna
  meio_data <- clean_data %>%
    filter(Tax.eco.meio == "TRUE" | Tax.eco.meio == TRUE)
  
  meio_global_list[[i]] <- data.frame(
    Dataset = datasets_ext$dataset_name[i],
    Total_Reads = sum(meio_data$nreads, na.rm = TRUE),
    Total_ASVs = nrow(meio_data)
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

# Table Summaries for Antarctica (Post-27 noise filter, unrarefied)
ant_global_all <- tibble(
  Dataset = "Antarctica", 
  Total_Reads = sum(comm_ant_27[, ant_asvs_clean$ASVid]), 
  Total_ASVs = n_distinct(ant_asvs_clean$ASVid)
)
ant_phylum_all <- ant_asvs_clean %>% 
  group_by(phylum = Best.group) %>% 
  summarise(Total_Unique_ASVs = n_distinct(ASVid), .groups = "drop") %>% 
  mutate(Dataset = "Antarctica")

comm_ant_meio <- comm_ant_27[, ant_meio_clean$ASVid, drop = FALSE]

ant_global_meio <- tibble(
  Dataset = "Antarctica", 
  Total_Reads = sum(comm_ant_meio), 
  Total_ASVs = n_distinct(ant_meio_clean$ASVid)
)

ant_phylum_meio <- ant_meio_clean %>%
  group_by(phylum = Best.group) %>%
  summarise(
    Total_Unique_ASVs = n_distinct(ASVid),
    Perc_under_95 = mean(pident.max < 95, na.rm = TRUE) * 100,
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
table3_under95   <- meio_phylum_clean %>%
  select(Dataset, phylum, Perc_under_95) %>%
  pivot_wider(names_from = Dataset, values_from = Perc_under_95, values_fill = 0) %>%
  arrange(phylum)

write_csv(table1_all, "Table1_All_ASVs_and_Reads.csv")
write_csv(table2_meiofauna, "Table2_Meiofauna_ASVs_and_Reads.csv")
write_csv(table3_under95, "Table3_Phylum_Matrix_Unique_Under95.csv")

# 6. Export Filtered Species Metadata (All original columns retained) ----------------

# A. Complete filtered dataset for Antarctica (Table 1: All ASVs post-27 threshold)
species_ant_filtered_all <- species_tot %>%
  filter(ASVid %in% ant_asvs_clean$ASVid) %>%
  mutate(nreads_post27 = ant_active_reads[ASVid]) # Añade el total de lecturas pos-filtro

# B. Meiofauna-only filtered dataset for Antarctica (Table 2)
species_ant_filtered_meio <- species_tot %>%
  filter(ASVid %in% ant_meio_clean$ASVid) %>%
  mutate(nreads_post27 = ant_meio_reads[ASVid])

# Export to CSV
write_csv(species_ant_filtered_all, "Species_Antarctica_Filtered_Table1_All.csv")
write_csv(species_ant_filtered_meio, "Species_Antarctica_Filtered_Table2_Meiofauna.csv")

# Section 2: NETWORK COMPARISON & SHARED ASVs----------------------------------------


# 1. Geographic Node Preparation
world_map <- map_data("world")

ecol_raw <- read.csv("stations network.csv", sep = ";") %>%
  mutate(
    lat = as.numeric(gsub(",", ".", latitude)),
    lon = as.numeric(gsub(",", ".", longitude)),
    lat = ifelse(regional_location == "Antarctica", -abs(lat), lat),
    studyid = as.character(studyid)
  ) %>% 
  filter(!is.na(lat))

coord_combined <- ecol_raw %>%
  filter(studyid %in% c("Ref6", "Ref7")) %>%
  summarise(studyid = "Ref6+7", lon = mean(lon), lat = mean(lat))

ecol_comb <- ecol_raw %>%
  filter(!studyid %in% c("Ref6", "Ref7")) %>%
  bind_rows(coord_combined) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

nodos_comb <- ecol_comb %>%
  st_drop_geometry() %>%
  select(id = studyid, lon, lat) %>%
  distinct(id, .keep_all = TRUE) %>%
  mutate(lon_jitter = jitter(lon, amount = 2), lat_jitter = jitter(lat, amount = 2))

edges_base_limpia <- read.csv("edges_clean.csv", sep = ";") %>%
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
df_meta <- read.csv("ASVs_metadata100_gaps.csv", sep = ";")
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


# Section 3: COVERAGE-BASED RAREFACTION & EXTRAPOLATION (iNEXT)-----------------------------------------------------

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
num_cores <- max(1, min(length(list_studies), detectCores() - 2))
cl <- makeCluster(num_cores)
clusterEvalQ(cl, library(iNEXT))
clusterExport(cl, "list_studies")

# 1. iNEXT Curves Parallel Computation
inext_parallel_list <- parLapply(cl, list_studies, function(study_data) {
  iNEXT(study_data, q = 0, datatype = "abundance", knots = 100, se = TRUE, conf = 0.95, nboot = 30)
})

stopCluster(cl)

# 2. Consolidate iNEXT Objects
study_names <- names(list_studies)
names(inext_parallel_list) <- study_names
for(nm in study_names) {
  inext_parallel_list[[nm]]$DataInfo$Assemblage <- nm
  inext_parallel_list[[nm]]$AsyEst$Assemblage <- nm
  if (is.data.frame(inext_parallel_list[[nm]]$iNextEst)) {
    inext_parallel_list[[nm]]$iNextEst$Assemblage <- nm
  } else {
    inext_parallel_list[[nm]]$iNextEst$size_based$Assemblage <- nm
    inext_parallel_list[[nm]]$iNextEst$coverage_based$Assemblage <- nm
  }
}

merge_inext_objects <- function(inext_list) {
  combined_DataInfo <- do.call(rbind, lapply(inext_list, function(x) x$DataInfo))
  combined_AsyEst   <- do.call(rbind, lapply(inext_list, function(x) x$AsyEst))
  
  if ("iNextEst" %in% names(inext_list[[1]])) {
    if (is.data.frame(inext_list[[1]]$iNextEst)) {
      combined_iNextEst <- do.call(rbind, lapply(inext_list, function(x) x$iNextEst))
    } else {
      combined_iNextEst <- list(
        size_based     = do.call(rbind, lapply(inext_list, function(x) x$iNextEst$size_based)),
        coverage_based = do.call(rbind, lapply(inext_list, function(x) x$iNextEst$coverage_based))
      )
    }
  }
  
  res <- list(DataInfo = combined_DataInfo, iNextEst = combined_iNextEst, AsyEst = combined_AsyEst)
  class(res) <- "iNEXT"
  return(res)
}

inext_out <- merge_inext_objects(inext_parallel_list)

# 3. Extract Table S5 directly from inext_out (Guarantees 100% match with Figure S2)
cov_data <- inext_out$iNextEst$coverage_based

table_s5_richness <- cov_data %>%
  filter(Order.q == 0) %>%
  group_by(Assemblage) %>%
  slice(which.min(abs(SC - 0.95))) %>%
  ungroup() %>%
  select(
    `Study Dataset` = Assemblage, 
    `Coverage Level` = SC, 
    `Standardized Richness (qD)` = qD, 
    `95% CI Lower (LCL)` = qD.LCL, 
    `95% CI Upper (UCL)` = qD.UCL
  ) %>%
  arrange(desc(`Standardized Richness (qD)`))

write.csv(table_s5_richness, "Table_S5_Standardized_Richness_q0.csv", row.names = FALSE)


# iNEXT Curves Parallel Computation
inext_parallel_list <- parLapply(cl, list_studies, function(study_data) {
  iNEXT(study_data, q = 0, datatype = "abundance", knots = 100, se = TRUE, conf = 0.95, nboot = 30)
})

stopCluster(cl)

# Consolidate iNEXT Objects
names(inext_parallel_list) <- study_names
for(nm in study_names) {
  inext_parallel_list[[nm]]$DataInfo$Assemblage <- nm
  inext_parallel_list[[nm]]$AsyEst$Assemblage <- nm
  if (is.data.frame(inext_parallel_list[[nm]]$iNextEst)) {
    inext_parallel_list[[nm]]$iNextEst$Assemblage <- nm
  } else {
    inext_parallel_list[[nm]]$iNextEst$size_based$Assemblage <- nm
    inext_parallel_list[[nm]]$iNextEst$coverage_based$Assemblage <- nm
  }
}

merge_inext_objects <- function(inext_list) {
  combined_DataInfo <- do.call(rbind, lapply(inext_list, function(x) x$DataInfo))
  combined_AsyEst   <- do.call(rbind, lapply(inext_list, function(x) x$AsyEst))
  
  if ("iNextEst" %in% names(inext_list[[1]])) {
    if (is.data.frame(inext_list[[1]]$iNextEst)) {
      combined_iNextEst <- do.call(rbind, lapply(inext_list, function(x) x$iNextEst))
    } else {
      combined_iNextEst <- list(
        size_based     = do.call(rbind, lapply(inext_list, function(x) x$iNextEst$size_based)),
        coverage_based = do.call(rbind, lapply(inext_list, function(x) x$iNextEst$coverage_based))
      )
    }
  }
  
  res <- list(DataInfo = combined_DataInfo, iNextEst = combined_iNextEst, AsyEst = combined_AsyEst)
  class(res) <- "iNEXT"
  return(res)
}

inext_out <- merge_inext_objects(inext_parallel_list)

# Figure S2:  Coverage-Based Rarefaction Plot
figure_S2 <- ggiNEXT(inext_out, type = 3) + 
  coord_cartesian(xlim = c(0.80, 1.00)) +
  scale_x_continuous(breaks = seq(0.80, 1.00, by = 0.05)) +
  theme_bw() +
  labs(
    x = "Sample Coverage", 
    y = "Standardized ASV Richness (q = 0)", 
    title = "Coverage-based Rarefaction and Extrapolation"
  ) +
  theme(
    legend.position = "right", 
    plot.title = element_text(face = "bold", size = 12)
  )

ggsave("Figure_S2_Rarefaction.png", plot = figure_S2, width = 8, height = 6, dpi = 300)

# Section 4: FIGURE 1 - SAMPLING MAPS------------------------------

# 1. Ross Sea Satellite Zoom
ecol_unique <- ecol_raw %>% distinct(longitude, latitude, .keep_all = TRUE)

register_stadiamaps(key = "e629ead7-e2f7-4cbb-9529-d204a2aaf84d")
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


# Section 5: FIGURE 2 - TAXONOMIC HEATMAP--------------------------------------

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




### Figure S1: Gaussian Distribution of ASVs < 95% Identity ----------------

# 1. Extract Global Percentages (<95% identity) directly from Section 1 Data

# External published studies
published_vals <- sapply(seq_len(nrow(datasets_ext)), function(i) {
  data <- read.csv(datasets_ext$filename[i], sep = ";", stringsAsFactors = FALSE)
  meio <- data %>% 
    filter(Eval.result == "correct", Tax.eco.meio %in% c("TRUE", TRUE)) %>%
    mutate(similitud = as.numeric(str_split_i(label, "\\|", 4)))
  mean(meio$similitud < 95.0, na.rm = TRUE) * 100
})

# Present study (Antarctica, post-27 noise threshold)
our_study_val <- mean(ant_meio_clean$pident.max < 95.0, na.rm = TRUE) * 100


# 2. Statistical Calculations (Fitted Gaussian Curve)
mean_val <- mean(published_vals, na.rm = TRUE)
sd_val   <- sd(published_vals, na.rm = TRUE)

# Dynamic X-axis range for curve rendering
x_min <- max(0, min(c(published_vals, our_study_val)) - 10)
x_max <- min(100, max(c(published_vals, our_study_val)) + 10)
x_val <- seq(x_min, x_max, length.out = 300)

df_curve     <- data.frame(x = x_val, y = dnorm(x_val, mean = mean_val, sd = sd_val))
df_published <- data.frame(x = published_vals, y = dnorm(published_vals, mean = mean_val, sd = sd_val))
df_our       <- data.frame(x = our_study_val, y = dnorm(our_study_val, mean = mean_val, sd = sd_val))


# 3. Render and Save Figure S1
figure_S1 <- ggplot() +
  # Gaussian Curve
  geom_line(data = df_curve, aes(x = x, y = y), linewidth = 1, color = "black") +
  # Benchmark Studies (Red Dots)
  geom_point(data = df_published, aes(x = x, y = y), color = "red", size = 3.5) +
  # Present Study (Blue Triangle)
  geom_point(data = df_our, aes(x = x, y = y), shape = 17, color = "blue", size = 5) +
  # Theme and Formatting
  labs(
    x = "percentage",
    y = "density"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black")
  )

ggsave("Figure_S1_Gaussian_Novelty.png", plot = figure_S1, width = 7, height = 6, dpi = 300)
