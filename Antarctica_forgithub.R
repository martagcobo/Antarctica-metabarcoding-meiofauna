#################################
# Antarctic meiofauna metabarcoding paper
#
# Script by: anonymus until peer-review is finished
# Last update: 17/05/2026

###########################


##### 1. Set working directory -------------------------------------------------
# Set to your local folder (Update this path as needed)
setwd("C:/")

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



#### 3. Prepare files ----------------------------------------------------------

# Load base data
ecol_base <- read.csv2("Antarctica marine MBC samples.csv")
species_base <- read.csv2("Data18S_Antarctica_v2.csv")
comm_base <- read.csv2("Antarctica_community.csv")

# Clean sequences
species_base <- species_base[
  !is.na(species_base$ASVid) &
    !is.na(species_base$sequence) &
    nzchar(species_base$sequence), 
]

# Prepare species dataframe (using which() to avoid NA issues)
species_tot <- species_base[which(species_base$Eval.result == "correct"), ]
species_tot <- species_tot[which(species_tot$rare.asvs == "FALSE"), ]
species_tot <- species_tot[!duplicated(species_tot$ASVid), ]

# Meiofauna subset
species <- species_tot[which(species_tot$Tax.eco.meio == "TRUE"), ]

# Prepare community matrix
rownames(comm_base) <- comm_base$sample_ID

# Remove specific rows/cols (Update indices if data changes)
comm_base <- comm_base[-(211:218), -c(1,2)] 
comm_base <- comm_base[-(81:88), ]

# Remove misstags
comm_base[comm_base < 27] <- 0

# --- TOTAL COMMUNITY ---
comm_tot <- comm_base[, colnames(comm_base) %in% species_tot$ASVid]

# Keep only common samples
common_ids_tot <- intersect(rownames(comm_tot), ecol_base$sample_ID)

comm_tot <- comm_tot[common_ids_tot, ]
ecol_tot <- ecol_base[match(common_ids_tot, ecol_base$sample_ID), ]

stopifnot(all(rownames(comm_tot) == ecol_tot$sample_ID))

# Total reads
ecol_tot$total_reads <- rowSums(comm_tot)

# Filter 0 reads
ecol_tot <- ecol_tot[ecol_tot$total_reads > 0, ]
comm_tot <- comm_tot[match(ecol_tot$sample_ID, rownames(comm_tot)), ]

# Presence/absence transformation
comm_tot[comm_tot > 0] <- 1

# Calculate Richness
ecol_tot$richness <- rowSums(comm_tot)

# Remove NAs in ecological variables
ecol_tot <- ecol_tot[!is.na(ecol_tot$depth) &
                       !is.na(ecol_tot$habitat_norm) &
                       !is.na(ecol_tot$Mesh), ]

comm_tot <- comm_tot[match(ecol_tot$sample_ID, rownames(comm_tot)), ]

# --- MEIOFAUNA SUBSET ---
comm <- comm_base[, colnames(comm_base) %in% species$ASVid]

# Keep only common samples
common_ids <- intersect(rownames(comm), ecol_base$sample_ID)

comm <- comm[common_ids, ]
ecol <- ecol_base[match(common_ids, ecol_base$sample_ID), ]

stopifnot(all(rownames(comm) == ecol$sample_ID))

# Total reads
ecol$total_reads <- rowSums(comm)

# Filter 0 reads
ecol <- ecol[ecol$total_reads > 0, ]
comm <- comm[match(ecol$sample_ID, rownames(comm)), ]

# Presence/absence transformation
comm[comm > 0] <- 1

# Calculate Richness
ecol$richness <- rowSums(comm)

# Remove NAs in ecological variables
ecol <- ecol[!is.na(ecol$depth) &
               !is.na(ecol$habitat_norm) &
               !is.na(ecol$Mesh), ]

comm <- comm[match(ecol$sample_ID, rownames(comm)), ]

# Set factors
habitat_order <- c("epilithic", "organic", "spicule", "gravel", "sand", "silt")
ecol$habitat_norm <- factor(ecol$habitat_norm, levels = habitat_order)
ecol$Mesh <- as.factor(ecol$Mesh)


#### 4. Descriptive Data -------------------------------------------------------

# 1. Total ASVs (Total Community)
cat("Total ASVs (all metazoans):", ncol(comm_tot), "\n")

cat("\nAll metazoans by group:\n")
asvs_per_group_tot <- table(species_tot$Best.group)
print(asvs_per_group_tot)

# 2. Total ASVs (Meiofauna subset)
cat("\nTotal ASVs (Meiofauna):", ncol(comm), "\n")

# 3. ASVs by group in Meiofauna
cat("\nMeiofauna ASVs by group:\n")
asvs_per_group <- table(species$Best.group)
print(asvs_per_group)

# 4. Meiofauna richness by sample
min_asvs <- min(ecol$richness)
max_asvs <- max(ecol$richness)
mean_asvs <- round(mean(ecol$richness), 1)

cat("\nMinimum ASVs per sample:", min_asvs, "\n")
cat("Maximum ASVs per sample:", max_asvs, "\n")
cat("Mean ASVs per sample:", mean_asvs, "\n")

# Identify specific samples with min/max ASVs
sample_min <- ecol$sample_ID[which.min(ecol$richness)]
sample_max <- ecol$sample_ID[which.max(ecol$richness)]

cat("\nThe sample with the fewest ASVs is:", sample_min, "\n")
cat("The sample with the most ASVs is:", sample_max, "\n")


################## ALPHA DIVERSITY #############################
############# A. Taxonomic Alpha Diversity #####################

#### MEIOFAUNA MODEL
model <- glmmTMB(richness ~ scale(depth) + Mesh + habitat_norm + offset(log(total_reads)) + (1 | ID), data = ecol, family = poisson)

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
  geom_line(color = "grey70", linetype = "dashed", linewidth = 0.5) +
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
  geom_line(color = "grey70", linetype = "dashed", linewidth = 0.5) +
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

## Tables with Alpha Taxonomic Results

# 0. Extract basic model info (N and Groups)
df_info <- data.frame(
  Parameter = c("Total Observations (N)", "Number of Groups (ID)"),
  Value = c(nobs(model), summary(model)$ngrps$ID)
)

# 1. Extract model coefficients
df_summary <- as.data.frame(summary(model)$coefficients$cond)
df_summary$Term <- rownames(df_summary) 
rownames(df_summary) <- NULL
df_summary <- df_summary[, c("Term", "Estimate", "Std. Error", "z value", "Pr(>|z|)")]

# 2. Extract performance metrics
df_r2 <- as.data.frame(performance::r2(model))
df_icc <- as.data.frame(performance::icc(model))
df_aic <- data.frame(AIC = AIC(model), BIC = BIC(model))
df_performance <- cbind(df_aic, df_r2, df_icc)

# 3. Extract Type II ANOVA results
anova_res <- car::Anova(model)
df_anova <- as.data.frame(anova_res)
df_anova$Term <- rownames(df_anova)
rownames(df_anova) <- NULL
df_anova <- df_anova[, c("Term", "Chisq", "Df", "Pr(>Chisq)")]

# 4. Extract post-hoc pairwise comparisons
em_mesh <- emmeans(model, pairwise ~ Mesh, type="response")
em_hab <- emmeans(model, pairwise ~ habitat_norm, type="response")

df_pair_mesh <- as.data.frame(em_mesh$contrasts)
df_pair_hab <- as.data.frame(em_hab$contrasts)

# 5. Compile into a list and export to Excel
excel_results_list <- list(
  "0_Model_Info" = df_info,          
  "1_Model_Summary" = df_summary,
  "2_Performance_Metrics" = df_performance,
  "3_ANOVA_TypeII" = df_anova,
  "4_Pairwise_Mesh" = df_pair_mesh,
  "5_Pairwise_Habitat" = df_pair_hab
)

write_xlsx(excel_results_list, path = "Meiofauna_alpha_tax.xlsx")


#### Copepoda (Alpha Tax) ------------------------------------------------------

# Filter community matrix
cop_asvs <- species %>%
  filter(Best.group == "Copepoda") %>%
  pull(ASVid)

comm_cop <- comm[, colnames(comm) %in% cop_asvs]

# Calculate richness
if(is.null(dim(comm_cop))) {
  ecol$richness_cop <- (comm_cop > 0) * 1
} else {
  ecol$richness_cop <- rowSums(comm_cop > 0)
}

# Model
model_cop <- glmmTMB(richness_cop ~ scale(depth) + Mesh + habitat_norm + offset(log(total_reads)) + (1 | ID), data = ecol, family = poisson) 

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
  geom_line(color = "#012334", linetype = "dashed", linewidth = 0.5) +
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

nem_asvs <- species %>% filter(Best.group == "Nematoda") %>% pull(ASVid)
comm_nem <- comm[, colnames(comm) %in% nem_asvs]

if(is.null(dim(comm_nem))) {
  ecol$richness_nem <- (comm_nem > 0) * 1
} else {
  ecol$richness_nem <- rowSums(comm_nem > 0)
}

model_nem <- glmmTMB(richness_nem ~ scale(depth) + Mesh + habitat_norm + offset(log(total_reads)) + (1 | ID), data = ecol, family = poisson) 

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
  geom_line(color = "#15adf8", linetype = "dashed", linewidth = 0.5) +
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

plat_asvs <- species %>% filter(Best.group == "Platyhelminthes") %>% pull(ASVid)
comm_plat <- comm[, colnames(comm) %in% plat_asvs]

if(is.null(dim(comm_plat))) {
  ecol$richness_plat <- (comm_plat > 0) * 1
} else {
  ecol$richness_plat <- rowSums(comm_plat > 0)
}

model_plat <- glmmTMB(richness_plat ~ scale(depth) + Mesh + habitat_norm + offset(log(total_reads)) + (1 | ID), data = ecol, family = poisson) 

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
  geom_line(color = "#d47000", linetype = "dashed", linewidth = 0.5) +
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


#### Ostracoda (Alpha Tax) -----------------------------------------------------

os_asvs <- species %>% filter(Best.group == "Ostracoda") %>% pull(ASVid)
comm_os <- comm[, colnames(comm) %in% os_asvs]

if(is.null(dim(comm_os))) {
  ecol$richness_os <- (comm_os > 0) * 1
} else {
  ecol$richness_os <- rowSums(comm_os > 0)
}

model_os <- glmmTMB(richness_os ~ scale(depth) + Mesh + habitat_norm + offset(log(total_reads)) + (1 | ID), data = ecol, family = poisson) 

performance::check_overdispersion(model_os)
pdf("check_model_os_output.pdf", width = 10, height = 7)
performance::check_model(model_os)
dev.off()

summary(model_os)
car::Anova(model_os)
emmeans(model_os, pairwise ~ Mesh, type="response")
emmeans(model_os, pairwise ~ habitat_norm, type="response")

prediction_mesh <- ggpredict(model_os, terms = "Mesh")
ggplot(prediction_mesh, aes(x = x, y = predicted, group = 1)) +
  geom_line(color = "#046493", linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#046493", alpha=0.8) + 
  labs(x = "Mesh size (µm)", y = "Predicted ASV Richness", title = "Ostracoda: Model-Adjusted Effects") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.line = element_line(color = "black"))

prediction_habitat <- ggpredict(model_os, terms = "habitat_norm")
ggplot(prediction_habitat, aes(x = x, y = predicted, group = 1)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#046493") + 
  labs(x = "Type of habitat", y = "Predicted ASV Richness", title = "Ostracoda: Model-Adjusted Effects") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.line = element_line(color = "black"))


#### Annelida (Alpha Tax) ------------------------------------------------------

ann_asvs <- species %>% filter(Best.group == "Annelida") %>% pull(ASVid)
comm_ann <- comm[, colnames(comm) %in% ann_asvs]

if(is.null(dim(comm_ann))) {
  ecol$richness_ann <- (comm_ann > 0) * 1
} else {
  ecol$richness_ann <- rowSums(comm_ann > 0)
}

model_ann <- glmmTMB(richness_ann ~ scale(depth) + Mesh + habitat_norm + offset(log(total_reads)) + (1 | ID), data = ecol,  family = poisson) 

performance::check_overdispersion(model_ann)
pdf("check_model_ann_output.pdf", width = 10, height = 7)
performance::check_model(model_ann)
dev.off()

summary(model_ann)
car::Anova(model_ann)
emmeans(model_ann, pairwise ~ Mesh, type="response")
emmeans(model_ann, pairwise ~ habitat_norm, type="response")

prediction_mesh <- ggpredict(model_ann, terms = "Mesh")
ggplot(prediction_mesh, aes(x = x, y = predicted, group = 1)) +
  geom_line(color = "#723c00", linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#723c00", alpha=0.8) + 
  labs(x = "Mesh size (µm)", y = "Predicted ASV Richness", title = "Annelida: Model-Adjusted Effects") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.line = element_line(color = "black"))

prediction_habitat <- ggpredict(model_ann, terms = "habitat_norm")
ggplot(prediction_habitat, aes(x = x, y = predicted, group = 1)) +
  geom_line(color = "#723c00", linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#723c00") + 
  labs(x = "Type of habitat", y = "Predicted ASV Richness", title = "Annelida: Model-Adjusted Effects") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.line = element_line(color = "black"))


#### Gastrotricha (Alpha Tax) --------------------------------------------------

gas_asvs <- species %>% filter(Best.group == "Gastrotricha") %>% pull(ASVid)
comm_gas <- comm[, colnames(comm) %in% gas_asvs]

if(is.null(dim(comm_gas))) {
  ecol$richness_gas <- (comm_gas > 0) * 1
} else {
  ecol$richness_gas <- rowSums(comm_gas > 0)
}

# Note: Uses nbinom2 instead of poisson
model_gas <- glmmTMB(richness_gas ~ scale(depth) + Mesh + habitat_norm + offset(log(total_reads)) + (1 | ID), data = ecol,  family = nbinom2) 

performance::check_overdispersion(model_gas)
pdf("check_model_gas_output.pdf", width = 10, height = 7)
performance::check_model(model_gas)
dev.off()

summary(model_gas)
car::Anova(model_gas)
emmeans(model_gas, pairwise ~ Mesh, type="response")
emmeans(model_gas, pairwise ~ habitat_norm, type="response")

prediction_mesh <- ggpredict(model_gas, terms = "Mesh")
ggplot(prediction_mesh, aes(x = x, y = predicted, group = 1)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#fb8500", alpha=0.8) + 
  labs(x = "Mesh size (µm)", y = "Predicted ASV Richness", title = "Gastrotricha: Model-Adjusted Effects") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.line = element_line(color = "black"))

prediction_habitat <- ggpredict(model_gas, terms = "habitat_norm")
ggplot(prediction_habitat, aes(x = x, y = predicted, group = 1)) +
  geom_line(color = "#fb8500", linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "#fb8500") + 
  labs(x = "Type of habitat", y = "Predicted ASV Richness", title = "Gastrotricha: Model-Adjusted Effects") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.line = element_line(color = "black"))


#### Xenacoelomorpha (Alpha Tax) -----------------------------------------------

xen_asvs <- species %>% filter(Best.group == "Xenacoelomorpha") %>% pull(ASVid)
comm_xen <- comm[, colnames(comm) %in% xen_asvs]

if(is.null(dim(comm_xen))) {
  ecol$richness_xen <- (comm_xen > 0) * 1
} else {
  ecol$richness_xen <- rowSums(comm_xen > 0)
}

model_xen <- glmmTMB(richness_xen ~ scale(depth) + Mesh + habitat_norm + offset(log(total_reads)) + (1 | ID), data = ecol, family = poisson) 

performance::check_overdispersion(model_xen)
pdf("check_model_xen_output.pdf", width = 10, height = 7)
performance::check_model(model_xen)
dev.off()

summary(model_xen)
car::Anova(model_xen)
emmeans(model_xen, pairwise ~ Mesh, type="response")
emmeans(model_xen, pairwise ~ habitat_norm, type="response")


prediction_mesh <- ggpredict(model_xen, terms = "Mesh")
ggplot(prediction_mesh, aes(x = x, y = predicted, group = 1)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "grey", alpha=0.8) + 
  labs(x = "Mesh size (µm)", y = "Predicted ASV Richness", title = "Xenacoelomorpha: Model-Adjusted Effects") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.line = element_line(color = "black"))

prediction_habitat <- ggpredict(model_xen, terms = "habitat_norm")
ggplot(prediction_habitat, aes(x = x, y = predicted, group = 1)) +
  geom_line(color = "grey", linetype = "dashed", linewidth = 0.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1, linewidth = 0.8, color = "#444444") +
  geom_point(size = 4, color = "grey") + 
  labs(x = "Type of habitat", y = "Predicted ASV Richness", title = "Xenacoelomorpha: Model-Adjusted Effects") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"), axis.line = element_line(color = "black"))



## Table with results alpha tax by group ---------------------------------------


# 1. Put all models in a list with their names
# The name on the left (in quotes) will be used to name the Excel file
model_list <- list(
  "meio" = model,
  "cop"  = model_cop,
  "nem"  = model_nem,
  "plat" = model_plat,
  "os"   = model_os,
  "ann"  = model_ann,
  "gas"  = model_gas,
  "xen"  = model_xen
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
  
  # 2. Extract performance metrics
  df_r2 <- as.data.frame(performance::r2(current_model))
  df_icc <- as.data.frame(performance::icc(current_model))
  df_aic <- data.frame(AIC = AIC(current_model), BIC = BIC(current_model))
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


############# B. Phylogenetic alpha diversity ####################################

# Clean and test sequences
species$sequence <- trimws(as.character(species$sequence))

species <- species[
  !is.na(species$sequence) &
    nzchar(species$sequence) &
    species$sequence != "NA",
]

any(is.na(species$sequence))
any(species$sequence == "")
any(grepl("^\\s+$", species$sequence))
any(species$sequence == "NA")

# Assuming 'species' has ASVid and sequence columns
seqs <- DNAStringSet(species$sequence)
names(seqs) <- species$ASVid
writeXStringSet(seqs, filepath = "Sequences_meioAntarctica18S.fasta", format = "fasta")

# Read FASTA file from the specified path
seqs <- readDNAStringSet("Sequences_meioAntarctica18S.fasta")

# Check the first lines
head(seqs)

## Tree building
alignment <- AlignSeqs(seqs)

# Convert alignment to a phyDat object
phy_data <- phyDat(as.matrix(alignment), type = "DNA")

# Make phylogenetic tree with Neighbor-Joining
dist_matrix <- dist.ml(phy_data)  # Distance matrix
tree <- NJ(dist_matrix)           # Tree with Neighbor-Joining

comm <- (comm > 0) * 1

pd_sample <- BAT::alpha(comm, tree = tree)
pd_sample <- as.data.frame(pd_sample)
colnames(pd_sample) <- "phylo.diver"
pd_sample$sample_ID <- rownames(pd_sample)

stations <- merge(ecol, pd_sample, by = "sample_ID", all.x = TRUE)
stations$phylo.diver[is.na(stations$phylo.diver)] <- 0


## Model
stations$Mesh <- as.factor(stations$Mesh)

mod.phylo <- glmmTMB(phylo.diver ~ scale(depth) + Mesh + habitat_norm + offset(log(total_reads)) + (1 | ID), 
                     data = stations, family = tweedie(link="log"))

## Check the model
performance::check_overdispersion(mod.phylo)
performance::check_collinearity(mod.phylo)
plot(mod.phylo)

# pdf("C:/Users/DEEPGAMING/OneDrive - Universidad Complutense de Madrid (UCM)/Mi unidad/Proyectos/Emilia-Romagna/check_model_2.pdf", width = 10, height = 7) 
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
# 1. Extract predictions
prediction_mesh <- ggpredict(mod.phylo, terms = "Mesh", bias_correction = TRUE)

# 2. "Predicted Values" plot, publication style
ggplot(prediction_mesh, aes(x = x, y = predicted, group = 1)) +
  
  # Confidence intervals ("error bars")
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width = 0.1, linewidth = 0.8, color = "#444444") +
  
  # Predicted points (Grey for Meiofauna)
  geom_point(size = 4, color = "#999999") + 
  
  # Labels and formatting
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
# 1. Extract predictions
prediction_habitat <- ggpredict(mod.phylo, terms = "habitat_norm", bias_correction = TRUE)

# 2. "Predicted Values" plot, publication style
ggplot(prediction_habitat, aes(x = x, y = predicted, group = 1)) +
  # Add a soft line to connect the trends (optional)
  geom_line(color = "grey70", linetype = "dashed", linewidth = 0.5) +
  
  # Confidence intervals ("error bars")
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width = 0.1, linewidth = 0.8, color = "#444444") +
  
  # Predicted points (Grey for Meiofauna)
  geom_point(size = 4, color = "#999999") + 
  
  # Labels and formatting
  labs(
    x = "Type of habitat",
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


## Tables with results alpha phyl 

# 0. Extract basic info from mod.phylo (N and Groups)
# nobs() gets total observations and summary()$ngrps gets random factor levels
df_info <- data.frame(
  Parameter = c("Total Observations (N)", "Number of Groups (ID)"),
  Value = c(nobs(mod.phylo), summary(mod.phylo)$ngrps$ID)
)

# 1. Extract model coefficients (Estimates, SE, z-values, p-values)
df_summary <- as.data.frame(summary(mod.phylo)$coefficients$cond)
df_summary$Term <- rownames(df_summary) # Move row names to a column
rownames(df_summary) <- NULL
df_summary <- df_summary[, c("Term", "Estimate", "Std. Error", "z value", "Pr(>|z|)")]

# 2. Extract performance metrics (Explicitly forcing R2, ICC, and AIC)
df_r2 <- as.data.frame(performance::r2(mod.phylo))
df_icc <- as.data.frame(performance::icc(mod.phylo))
df_aic <- data.frame(AIC = AIC(mod.phylo), BIC = BIC(mod.phylo))

# Combine into a single row (ensuring no duplicate columns)
df_performance <- cbind(df_aic, df_r2, df_icc)

# 3. Extract Type II ANOVA results
anova_res <- car::Anova(mod.phylo)
df_anova <- as.data.frame(anova_res)
df_anova$Term <- rownames(df_anova)
rownames(df_anova) <- NULL
df_anova <- df_anova[, c("Term", "Chisq", "Df", "Pr(>Chisq)")]

# 4. Extract post-hoc (Pairwise comparisons)
em_mesh <- emmeans(mod.phylo, pairwise ~ Mesh, type="response")
em_hab <- emmeans(mod.phylo, pairwise ~ habitat_norm, type="response")

df_pair_mesh <- as.data.frame(em_mesh$contrasts)
df_pair_hab <- as.data.frame(em_hab$contrasts)

# 5. Join everything into a list and export to a single Excel with multiple tabs
excel_results_list <- list(
  "0_mod.phylo_Info" = df_info,         
  "1_mod.phylo_Summary" = df_summary,
  "2_Performance_Metrics" = df_performance,
  "3_ANOVA_TypeII" = df_anova,
  "4_Pairwise_Mesh" = df_pair_mesh,
  "5_Pairwise_Habitat" = df_pair_hab
)

# Export file to working directory
write_xlsx(excel_results_list, path = "Meiofauna_alpha_phyl.xlsx")


## Phylogenetic richness by taxa

## Copepoda (alpha phyl) -------------------------------------------------------

species_cop <- species %>% dplyr::filter(Best.group == "Copepoda")

seqs <- DNAStringSet(species_cop$sequence)
names(seqs) <- species_cop$ASVid
writeXStringSet(seqs, filepath = "Sequences_copAntarctica18S.fasta", format = "fasta")

# Alignment and NJ tree
alignment <- AlignSeqs(seqs)
phy_data <- phyDat(as.matrix(alignment), type = "DNA")
dist_matrix <- dist.ml(phy_data)
tree_cop <- NJ(dist_matrix)

## Community Matrix
asvs_cop <- species_cop$ASVid
comm_cop_subset <- comm[, colnames(comm) %in% asvs_cop]

# Presence/Absence
comm_cop_binary <- as.matrix((comm_cop_subset > 0) * 1)

## Phylogenetic richness (PD)
pd_cop_res <- BAT::alpha(comm_cop_binary, tree = tree_cop)
pd_cop <- data.frame(
  sample_ID = rownames(pd_cop_res),
  phylo.diver = as.numeric(pd_cop_res[,1])
)

## Merge with ecol 
stations_cop <- merge(ecol, pd_cop, by = "sample_ID", all.x = TRUE)

# Substitute NAs with 0
stations_cop$phylo.diver[is.na(stations_cop$phylo.diver)] <- 0

## Phylogenetic richness model
stations_cop$Mesh <- as.factor(stations_cop$Mesh)

mod.phylo_cop <- glmmTMB(phylo.diver ~ scale(depth) + Mesh + habitat_norm + offset(log(total_reads)) + (1 | ID), 
                         data = stations_cop, family = tweedie(link="log"))

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


## Nematoda (alpha phyl) -------------------------------------------------------

species_nem <- species %>% dplyr::filter(Best.group == "Nematoda")

seqs <- DNAStringSet(species_nem$sequence)
names(seqs) <- species_nem$ASVid
writeXStringSet(seqs, filepath = "Sequences_nemAntarctica18S.fasta", format = "fasta")

# Alignment and NJ tree
alignment <- AlignSeqs(seqs)
phy_data <- phyDat(as.matrix(alignment), type = "DNA")
dist_matrix <- dist.ml(phy_data)
tree_nem <- NJ(dist_matrix)

## Community Matrix
asvs_nem <- species_nem$ASVid
comm_nem_subset <- comm[, colnames(comm) %in% asvs_nem]

# Presence/Absence
comm_nem_binary <- as.matrix((comm_nem_subset > 0) * 1)

## Phylogenetic richness (PD)
pd_nem_res <- BAT::alpha(comm_nem_binary, tree = tree_nem)
pd_nem <- data.frame(
  sample_ID = rownames(pd_nem_res),
  phylo.diver = as.numeric(pd_nem_res[,1])
)

## Merge with ecol 
stations_nem <- merge(ecol, pd_nem, by = "sample_ID", all.x = TRUE)

# Substitute NAs with 0
stations_nem$phylo.diver[is.na(stations_nem$phylo.diver)] <- 0

## Phylogenetic richness model
stations_nem$Mesh <- as.factor(stations_nem$Mesh)

mod.phylo_nem <- glmmTMB(phylo.diver ~ scale(depth) + Mesh + habitat_norm + offset(log(total_reads)) + (1 | ID), 
                         data = stations_nem, family = tweedie(link="log"))

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


## Platyhelminthes (alpha phyl) ------------------------------------------------

species_plat <- species %>% dplyr::filter(Best.group == "Platyhelminthes")

seqs <- DNAStringSet(species_plat$sequence)
names(seqs) <- species_plat$ASVid
writeXStringSet(seqs, filepath = "Sequences_platAntarctica18S.fasta", format = "fasta")

# Alignment and NJ tree
alignment <- AlignSeqs(seqs)
phy_data <- phyDat(as.matrix(alignment), type = "DNA")
dist_matrix <- dist.ml(phy_data)
tree_plat <- NJ(dist_matrix)

## Community Matrix
asvs_plat <- species_plat$ASVid
comm_plat_subset <- comm[, colnames(comm) %in% asvs_plat]

# Presence/Absence
comm_plat_binary <- as.matrix((comm_plat_subset > 0) * 1)

## Phylogenetic richness (PD)
pd_plat_res <- BAT::alpha(comm_plat_binary, tree = tree_plat)
pd_plat <- data.frame(
  sample_ID = rownames(pd_plat_res),
  phylo.diver = as.numeric(pd_plat_res[,1])
)

## Merge with ecol 
stations_plat <- merge(ecol, pd_plat, by = "sample_ID", all.x = TRUE)

# Substitute NAs with 0
stations_plat$phylo.diver[is.na(stations_plat$phylo.diver)] <- 0

## Phylogenetic richness model
stations_plat$Mesh <- as.factor(stations_plat$Mesh)

mod.phylo_plat <- glmmTMB(phylo.diver ~ scale(depth) + Mesh + habitat_norm + offset(log(total_reads)) + (1 | ID), 
                          data = stations_plat, family = tweedie(link="log"))

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


## Ostracoda (alpha phyl) ------------------------------------------------------

species_os <- species %>% dplyr::filter(Best.group == "Ostracoda")

seqs <- DNAStringSet(species_os$sequence)
names(seqs) <- species_os$ASVid
writeXStringSet(seqs, filepath = "Sequences_osAntarctica18S.fasta", format = "fasta")

# Alignment and NJ tree
alignment <- AlignSeqs(seqs)
phy_data <- phyDat(as.matrix(alignment), type = "DNA")
dist_matrix <- dist.ml(phy_data)
tree_os <- NJ(dist_matrix)

## Community Matrix
asvs_os <- species_os$ASVid
comm_os_subset <- comm[, colnames(comm) %in% asvs_os]

# Presence/Absence
comm_os_binary <- as.matrix((comm_os_subset > 0) * 1)

## Phylogenetic richness (PD)
pd_os_res <- BAT::alpha(comm_os_binary, tree = tree_os)
pd_os <- data.frame(
  sample_ID = rownames(pd_os_res),
  phylo.diver = as.numeric(pd_os_res[,1])
)

## Merge with ecol 
stations_os <- merge(ecol, pd_os, by = "sample_ID", all.x = TRUE)

# Substitute NAs with 0
stations_os$phylo.diver[is.na(stations_os$phylo.diver)] <- 0

## Phylogenetic richness model
stations_os$Mesh <- as.factor(stations_os$Mesh)

mod.phylo_os <- glmmTMB(phylo.diver ~ scale(depth) + Mesh + habitat_norm + offset(log(total_reads)) + (1 | ID), 
                        data = stations_os, family = tweedie(link="log"))

## Check the model
performance::check_overdispersion(mod.phylo_os)
performance::check_collinearity(mod.phylo_os)

# pdf(...)
performance::check_model(mod.phylo_os)
# dev.off()

## Summary and post-hoc
car::Anova(mod.phylo_os)
summary(mod.phylo_os)

## Post-hoc
emmeans(mod.phylo_os, pairwise ~ Mesh, type="response")
emmeans(mod.phylo_os, pairwise ~ habitat_norm, type="response")

# Habitat
# 1. Extract predictions (using mod.phylo_os)
prediction_habitat <- ggpredict(mod.phylo_os, terms = "habitat_norm")

# 2. "Predicted Values" plot, publication style
ggplot(prediction_habitat, aes(x = x, y = predicted, group = 1)) +
  
  # Confidence intervals ("error bars")
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width = 0.1, linewidth = 0.8, color = "#444444") +
  
  # Predicted points (Dark blue for Ostracoda)
  geom_point(size = 4, color = "#046493") + 
  
  # Labels and formatting
  labs(
    x = "Type of habitat",
    y = "Predicted Phylogenetic Richness",
    title = "Ostracoda: Model-Adjusted Effects",
    subtitle = ""
  ) +
  
  # Visual cleanup
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.line = element_line(color = "black")
  )


## Annelida (alpha phyl) -------------------------------------------------------

species_ann <- species %>% dplyr::filter(Best.group == "Annelida")

seqs <- DNAStringSet(species_ann$sequence)
names(seqs) <- species_ann$ASVid
writeXStringSet(seqs, filepath = "Sequences_annAntarctica18S.fasta", format = "fasta")

# Alignment and NJ tree
alignment <- AlignSeqs(seqs)
phy_data <- phyDat(as.matrix(alignment), type = "DNA")
dist_matrix <- dist.ml(phy_data)
tree_ann <- NJ(dist_matrix)

## Community Matrix
asvs_ann <- species_ann$ASVid
comm_ann_subset <- comm[, colnames(comm) %in% asvs_ann]

# Presence/Absence
comm_ann_binary <- as.matrix((comm_ann_subset > 0) * 1)

## Phylogenetic richness (PD)
pd_ann_res <- BAT::alpha(comm_ann_binary, tree = tree_ann)
pd_ann <- data.frame(
  sample_ID = rownames(pd_ann_res),
  phylo.diver = as.numeric(pd_ann_res[,1])
)

## Merge with ecol 
stations_ann <- merge(ecol, pd_ann, by = "sample_ID", all.x = TRUE)

# Substitute NAs with 0
stations_ann$phylo.diver[is.na(stations_ann$phylo.diver)] <- 0

## Phylogenetic richness model
stations_ann$Mesh <- as.factor(stations_ann$Mesh)

mod.phylo_ann <- glmmTMB(phylo.diver ~ scale(depth) + Mesh + habitat_norm + offset(log(total_reads)) + (1 | ID), 
                         data = stations_ann, family = tweedie(link="log"))

## Check the model
performance::check_overdispersion(mod.phylo_ann)
performance::check_collinearity(mod.phylo_ann)

# pdf(...)
performance::check_model(mod.phylo_ann)
# dev.off()

## Summary and post-hoc
car::Anova(mod.phylo_ann)
summary(mod.phylo_ann)

## Post-hoc
emmeans(mod.phylo_ann, pairwise ~ Mesh, type="response")
emmeans(mod.phylo_ann, pairwise ~ habitat_norm, type="response")


## Gastrotricha (alpha phyl) ---------------------------------------------------

species_gas <- species %>% dplyr::filter(Best.group == "Gastrotricha")

seqs <- DNAStringSet(species_gas$sequence)
names(seqs) <- species_gas$ASVid
writeXStringSet(seqs, filepath = "Sequences_gasAntarctica18S.fasta", format = "fasta")

# Alignment and NJ tree
alignment <- AlignSeqs(seqs)
phy_data <- phyDat(as.matrix(alignment), type = "DNA")
dist_matrix <- dist.ml(phy_data)
tree_gas <- NJ(dist_matrix)

## Community Matrix
asvs_gas <- species_gas$ASVid
comm_gas_subset <- comm[, colnames(comm) %in% asvs_gas]

# Presence/Absence
comm_gas_binary <- as.matrix((comm_gas_subset > 0) * 1)

## Phylogenetic richness (PD)
pd_gas_res <- BAT::alpha(comm_gas_binary, tree = tree_gas)
pd_gas <- data.frame(
  sample_ID = rownames(pd_gas_res),
  phylo.diver = as.numeric(pd_gas_res[,1])
)

## Merge with ecol 
stations_gas <- merge(ecol, pd_gas, by = "sample_ID", all.x = TRUE)

# Substitute NAs with 0
stations_gas$phylo.diver[is.na(stations_gas$phylo.diver)] <- 0

## Phylogenetic richness model
stations_gas$Mesh <- as.factor(stations_gas$Mesh)

mod.phylo_gas <- glmmTMB(phylo.diver ~ scale(depth) + Mesh + habitat_norm + offset(log(total_reads)) + (1 | ID), 
                         data = stations_gas, family = tweedie(link="log"))

## Check the model
performance::check_overdispersion(mod.phylo_gas)
performance::check_collinearity(mod.phylo_gas)

# pdf(...)
performance::check_model(mod.phylo_gas)
# dev.off()

## Summary and post-hoc
car::Anova(mod.phylo_gas)
summary(mod.phylo_gas)

## Post-hoc
emmeans(mod.phylo_gas, pairwise ~ Mesh, type="response")
emmeans(mod.phylo_gas, pairwise ~ habitat_norm, type="response")

## Box plot (predicted)
# Mesh
# 1. Extract predictions
prediction_mesh <- ggpredict(mod.phylo_gas, terms = "Mesh")

# 2. "Predicted Values" plot, publication style
ggplot(prediction_mesh, aes(x = x, y = predicted, group = 1)) +
  
  # Confidence intervals ("error bars")
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width = 0.1, linewidth = 0.8, color = "#444444") +
  
  # Predicted points (Light orange for Gastrotricha)
  geom_point(size = 4, color = "#fb8500", alpha=0.8) + 
  
  # Labels and formatting
  labs(
    x = "Mesh size (µm)",
    y = "Predicted Phylogenetic Richness",
    title = "Gastrotricha: Model-Adjusted Effects",
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

# 1. Put all phylogeny models in a list (excluding mites)
phylo_model_list <- list(
  "cop"  = mod.phylo_cop,
  "nem"  = mod.phylo_nem,
  "plat" = mod.phylo_plat,
  "os"   = mod.phylo_os,
  "gas"  = mod.phylo_gas
)

# 2. Create a fail-safe LOOP for the phylogenetic models
for (taxon_name in names(phylo_model_list)) {
  
  current_model <- phylo_model_list[[taxon_name]]
  cat("Processing phylogenetic model for:", taxon_name, "...\n")
  
  # 0. Extract basic info (Protected against models without variance in ID)
  groups_ID <- summary(current_model)$ngrps$ID
  if(is.null(groups_ID)) groups_ID <- 0 
  
  df_info <- data.frame(
    Parameter = c("Total Observations (N)", "Number of Groups (ID)"),
    Value = c(nobs(current_model), groups_ID)
  )
  
  # 1. Extract model coefficients
  df_summary <- as.data.frame(summary(current_model)$coefficients$cond)
  df_summary$Term <- rownames(df_summary) 
  rownames(df_summary) <- NULL
  df_summary <- df_summary[, c("Term", "Estimate", "Std. Error", "z value", "Pr(>|z|)")]
  
  # 2. Extract performance metrics using tryCatch as a failsafe
  df_aic <- data.frame(AIC = AIC(current_model), BIC = BIC(current_model))
  
  # Try to extract R2 (sometimes fails in tweedie if random variance is 0)
  df_r2 <- tryCatch({
    as.data.frame(performance::r2(current_model))
  }, error = function(e) {
    data.frame(R2_conditional = NA, R2_marginal = NA, Note = "Math error in R2 (Tweedie)")
  })
  
  # Try to extract ICC
  df_icc <- tryCatch({
    as.data.frame(performance::icc(current_model))
  }, error = function(e) {
    data.frame(ICC = NA)
  })
  
  # Merge everything into a single table
  df_performance <- cbind(df_aic, df_r2, df_icc)
  
  # 3. Extract Type II ANOVA results
  anova_res <- car::Anova(current_model)
  df_anova <- as.data.frame(anova_res)
  df_anova$Term <- rownames(df_anova)
  rownames(df_anova) <- NULL
  df_anova <- df_anova[, c("Term", "Chisq", "Df", "Pr(>Chisq)")]
  
  # 4. Extract post-hoc (Pairwise)
  # Use tryCatch here too in case any model fails during contrasts
  df_pair_mesh <- tryCatch({
    em_mesh <- emmeans(current_model, pairwise ~ Mesh, type="response")
    as.data.frame(em_mesh$contrasts)
  }, error = function(e) data.frame(Note = "Error in pairwise for Mesh"))
  
  df_pair_hab <- tryCatch({
    em_hab <- emmeans(current_model, pairwise ~ habitat_norm, type="response")
    as.data.frame(em_hab$contrasts)
  }, error = function(e) data.frame(Note = "Error in pairwise for Habitat"))
  
  # 5. Combine everything into a list of tabs for Excel
  excel_results_list <- list(
    "0_Model_Info" = df_info,
    "1_Model_Summary" = df_summary,
    "2_Performance_Metrics" = df_performance,
    "3_ANOVA_TypeII" = df_anova,
    "4_Pairwise_Mesh" = df_pair_mesh,
    "5_Pairwise_Habitat" = df_pair_hab
  )
  
  # 6. Create the filename dynamically (Note: 'phyl' added to the name)
  file_name <- paste0("Meiofauna_alpha_phyl_", taxon_name, ".xlsx")
  
  # Export to Excel
  write_xlsx(excel_results_list, path = file_name)
}


################### BETA DIVERSITY #############################################

## Helper Function for Exporting Results
## Defined once here to avoid repeating it in every taxon block


clean_adonis_table <- function(adonis_model, table_title) {
  
  # 1. Convert to dataframe and extract row names
  df <- as.data.frame(adonis_model)
  df$Variable <- rownames(df)
  
  # 2. Manual renaming (Bulletproof approach)
  # Look for the "Pr(>F)" column and rename it to "P_value"
  names(df)[names(df) == "Pr(>F)"] <- "P_value"
  
  # Look for the "F" column and rename it to "F_Model" 
  # (Avoids R interpreting it as FALSE)
  names(df)[names(df) == "F"] <- "F_Model"
  
  # 3. Selection and rounding
  # Ensure we only take columns that actually exist
  final_cols <- c("Variable", "Df", "R2", "F_Model", "P_value")
  
  # Filter only the columns present in the dataframe
  existing_cols <- intersect(final_cols, names(df))
  df_clean <- df[, existing_cols]
  
  # Round numeric columns
  df_clean <- df_clean %>%
    mutate(across(where(is.numeric), \(x) round(x, 3)))
  
  # 4. Create the flextable
  ft <- flextable(df_clean) %>%
    set_caption(caption = table_title) %>%
    autofit() %>%
    bold(part = "header") 
  
  return(ft)
}


############## A. Taxonomic beta diversity ---------------------------------------

### Meiofauna community (beta tax) ---------------------------------------------

beta_diversity <- BAT::beta(comm, func = "jaccard")  

str(beta_diversity)

ecol$Mesh <- as.factor(ecol$Mesh)

## PERMANOVA
adonis_meio <- adonis2(beta_diversity$Btotal ~ scale(depth) + Mesh + habitat_norm + log(total_reads+1), 
                       data = ecol, strata = ecol$ID, by = "margin")
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

## Tables for results
table1 <- clean_adonis_table(adonis_meio, "Table 1: Adonis Beta Total - Meiofauna")

# Save table
save_as_docx(
  table1, 
  path = "Results_Adonis_Meio_Tax.docx"
)

# Estimates for "Mesh" and "Habitat"
pairwise.adonis2(comm ~ Mesh, data = ecol, method = "jaccard")
pairwise.adonis2(comm ~ habitat_norm, data = ecol, method = "jaccard")

# Check collinearity (requires adonis model, make sure 'adonis' was defined previously if using car::vif)
# car::vif(adonis_meio)


### Community composition by taxa 

## Copepoda (beta tax) ---------------------------------------------------------

cop_asvs <- species %>%
  filter(Best.group == "Copepoda") %>%
  pull(ASVid)

comm_cop <- comm[, colnames(comm) %in% cop_asvs]

# Copepoda richness
if(is.null(dim(comm_cop))) {
  ecol$richness_cop <- (comm_cop > 0) * 1
} else {
  ecol$richness_cop <- rowSums(comm_cop > 0)
}

# Calculate beta diversity (Jaccard)
beta_diversity_cop <- beta(comm_cop, func = "jaccard")  
str(beta_diversity_cop)

## PERMANOVA
adonis_cop <- adonis2(beta_diversity_cop$Btotal ~ scale(depth) + Mesh + habitat_norm + log(total_reads+1), 
                      data = ecol, strata = ecol$ID, by = "margin")
print(adonis_cop)

# Turnover and nestedness proportions
mean_total <- mean(beta_diversity_cop$Btotal)
mean_turnover <- mean(beta_diversity_cop$Brepl)
mean_nestedness <- mean(beta_diversity_cop$Brich)

per_turnover <- (mean_turnover / mean_total) * 100
per_nestedness <- (mean_nestedness / mean_total) * 100

print(paste("Turnover represents:", round(per_turnover, 2), "%"))
print(paste("Nestedness represents:", round(per_nestedness, 2), "%"))

## Tables for results
table1 <- clean_adonis_table(adonis_cop, "Table 1: Adonis Beta Total - Copepoda")

# Save table
save_as_docx(
  table1, 
  path = "Results_Adonis_Cop_Tax.docx"
)

pairwise.adonis2(comm_cop ~ Mesh, data = ecol, method = "jaccard")
pairwise.adonis2(comm_cop ~ habitat_norm, data = ecol, method = "jaccard")


## Nematoda (beta tax) ---------------------------------------------------------

nem_asvs <- species %>%
  filter(Best.group == "Nematoda") %>%
  pull(ASVid)

comm_nem <- comm[, colnames(comm) %in% nem_asvs]

# Nematoda richness
if(is.null(dim(comm_nem))) {
  ecol$richness_nem <- (comm_nem > 0) * 1
} else {
  ecol$richness_nem <- rowSums(comm_nem > 0)
}

# Calculate beta diversity (Jaccard)
beta_diversity_nem <- beta(comm_nem, func = "jaccard")  
str(beta_diversity_nem)

## PERMANOVA
adonis_nem <- adonis2(beta_diversity_nem$Btotal ~ scale(depth) + Mesh + habitat_norm + log(total_reads+1), 
                      data = ecol, strata = ecol$ID, by = "margin")
print(adonis_nem)

# Turnover and nestedness proportions
mean_total <- mean(beta_diversity_nem$Btotal)
mean_turnover <- mean(beta_diversity_nem$Brepl)
mean_nestedness <- mean(beta_diversity_nem$Brich)

per_turnover <- (mean_turnover / mean_total) * 100
per_nestedness <- (mean_nestedness / mean_total) * 100

print(paste("Turnover represents:", round(per_turnover, 2), "%"))
print(paste("Nestedness represents:", round(per_nestedness, 2), "%"))

## Tables for results
table1 <- clean_adonis_table(adonis_nem, "Table 1: Adonis Beta Total - Nematoda")

# Save table
save_as_docx(
  table1, 
  path = "Results_Adonis_Nem_Tax.docx"
)

pairwise.adonis2(comm_nem ~ Mesh, data = ecol, method = "jaccard")
pairwise.adonis2(comm_nem ~ habitat_norm, data = ecol, method = "jaccard")


## Platyhelminthes (beta tax) --------------------------------------------------

plat_asvs <- species %>%
  filter(Best.group == "Platyhelminthes") %>%
  pull(ASVid)

comm_plat <- comm[, colnames(comm) %in% plat_asvs]

# Platyhelminthes richness
if(is.null(dim(comm_plat))) {
  ecol$richness_plat <- (comm_plat > 0) * 1
} else {
  ecol$richness_plat <- rowSums(comm_plat > 0)
}

# Calculate beta diversity (Jaccard)
beta_diversity_plat <- beta(comm_plat, func = "jaccard")  
str(beta_diversity_plat)

## PERMANOVA
adonis_plat <- adonis2(beta_diversity_plat$Btotal ~ scale(depth) + Mesh + habitat_norm + log(total_reads+1), 
                       data = ecol, strata = ecol$ID, by = "margin")
print(adonis_plat)

# Turnover and nestedness proportions
mean_total <- mean(beta_diversity_plat$Btotal)
mean_turnover <- mean(beta_diversity_plat$Brepl)
mean_nestedness <- mean(beta_diversity_plat$Brich)

per_turnover <- (mean_turnover / mean_total) * 100
per_nestedness <- (mean_nestedness / mean_total) * 100

print(paste("Turnover represents:", round(per_turnover, 2), "%"))
print(paste("Nestedness represents:", round(per_nestedness, 2), "%"))

## Tables for results
table1 <- clean_adonis_table(adonis_plat, "Table 1: Adonis Beta Total - Platyhelminthes")

# Save table
save_as_docx(
  table1, 
  path = "Results_Adonis_Plat_Tax.docx"
)

pairwise.adonis2(comm_plat ~ Mesh, data = ecol, method = "jaccard")
pairwise.adonis2(comm_plat ~ habitat_norm, data = ecol, method = "jaccard")


## Ostracoda (beta tax) --------------------------------------------------------

os_asvs <- species %>%
  filter(Best.group == "Ostracoda") %>%
  pull(ASVid)

comm_os <- comm[, colnames(comm) %in% os_asvs]

# Ostracoda richness
if(is.null(dim(comm_os))) {
  ecol$richness_os <- (comm_os > 0) * 1
} else {
  ecol$richness_os <- rowSums(comm_os > 0)
}

# Calculate beta diversity (Jaccard)
beta_diversity_os <- beta(comm_os, func = "jaccard")  
str(beta_diversity_os)

## PERMANOVA
adonis_os <- adonis2(beta_diversity_os$Btotal ~ scale(depth) + Mesh + habitat_norm + log(total_reads+1), 
                     data = ecol, strata = ecol$ID, by = "margin")
print(adonis_os)

# Turnover and nestedness proportions
mean_total <- mean(beta_diversity_os$Btotal)
mean_turnover <- mean(beta_diversity_os$Brepl)
mean_nestedness <- mean(beta_diversity_os$Brich)

per_turnover <- (mean_turnover / mean_total) * 100
per_nestedness <- (mean_nestedness / mean_total) * 100

print(paste("Turnover represents:", round(per_turnover, 2), "%"))
print(paste("Nestedness represents:", round(per_nestedness, 2), "%"))

## Tables for results
table1 <- clean_adonis_table(adonis_os, "Table 1: Adonis Beta Total - Ostracoda")

# Save table
save_as_docx(
  table1, 
  path = "Results_Adonis_Os_Tax.docx"
)

pairwise.adonis2(comm_os ~ Mesh, data = ecol, method = "jaccard")
pairwise.adonis2(comm_os ~ habitat_norm, data = ecol, method = "jaccard")


## Annelida (beta tax) ---------------------------------------------------------

ann_asvs <- species %>%
  filter(Best.group == "Annelida") %>%
  pull(ASVid)

comm_ann <- comm[, colnames(comm) %in% ann_asvs]

# Annelida richness
if(is.null(dim(comm_ann))) {
  ecol$richness_ann <- (comm_ann > 0) * 1
} else {
  ecol$richness_ann <- rowSums(comm_ann > 0)
}

# Calculate beta diversity (Jaccard)
beta_diversity_ann <- beta(comm_ann, func = "jaccard")  
str(beta_diversity_ann)

## PERMANOVA
adonis_ann <- adonis2(beta_diversity_ann$Btotal ~ scale(depth) + Mesh + habitat_norm + log(total_reads+1), 
                      data = ecol, strata = ecol$ID, by = "margin")
print(adonis_ann)

# Turnover and nestedness proportions
mean_total <- mean(beta_diversity_ann$Btotal)
mean_turnover <- mean(beta_diversity_ann$Brepl)
mean_nestedness <- mean(beta_diversity_ann$Brich)

per_turnover <- (mean_turnover / mean_total) * 100
per_nestedness <- (mean_nestedness / mean_total) * 100

print(paste("Turnover represents:", round(per_turnover, 2), "%"))
print(paste("Nestedness represents:", round(per_nestedness, 2), "%"))

## Tables for results
table1 <- clean_adonis_table(adonis_ann, "Table 1: Adonis Beta Total - Annelida")

# Save table
save_as_docx(
  table1, 
  path = "Results_Adonis_Ann_Tax.docx"
)


## Gastrotricha (beta tax) -----------------------------------------------------

gas_asvs <- species %>%
  filter(Best.group == "Gastrotricha") %>%
  pull(ASVid)

comm_gas <- comm[, colnames(comm) %in% gas_asvs]

# Gastrotricha richness
if(is.null(dim(comm_gas))) {
  ecol$richness_gas <- (comm_gas > 0) * 1
} else {
  ecol$richness_gas <- rowSums(comm_gas > 0)
}

# Calculate beta diversity (Jaccard)
beta_diversity_gas <- beta(comm_gas, func = "jaccard")  
str(beta_diversity_gas)

## PERMANOVA
adonis_gas <- adonis2(beta_diversity_gas$Btotal ~ scale(depth) + Mesh + habitat_norm + log(total_reads+1), 
                      data = ecol, strata = ecol$ID, by = "margin")
print(adonis_gas)

# Turnover and nestedness proportions
mean_total <- mean(beta_diversity_gas$Btotal)
mean_turnover <- mean(beta_diversity_gas$Brepl)
mean_nestedness <- mean(beta_diversity_gas$Brich)

per_turnover <- (mean_turnover / mean_total) * 100
per_nestedness <- (mean_nestedness / mean_total) * 100

print(paste("Turnover represents:", round(per_turnover, 2), "%"))
print(paste("Nestedness represents:", round(per_nestedness, 2), "%"))

## Tables for results
table1 <- clean_adonis_table(adonis_gas, "Table 1: Adonis Beta Total - Gastrotricha")

# Save table
save_as_docx(
  table1, 
  path = "Results_Adonis_Gas_Tax.docx"
)

pairwise.adonis2(comm_gas ~ Mesh, data = ecol, method = "jaccard")
pairwise.adonis2(comm_gas ~ habitat_norm, data = ecol, method = "jaccard")


## Xenacoelomorpha (beta tax) --------------------------------------------------

xen_asvs <- species %>%
  filter(Best.group == "Xenacoelomorpha") %>%
  pull(ASVid)

comm_xen <- comm[, colnames(comm) %in% xen_asvs]

# Xenacoelomorpha richness
if(is.null(dim(comm_xen))) {
  ecol$richness_xen <- (comm_xen > 0) * 1
} else {
  ecol$richness_xen <- rowSums(comm_xen > 0)
}

# Calculate beta diversity (Jaccard)
beta_diversity_xen <- beta(comm_xen, func = "jaccard")  
str(beta_diversity_xen)

## PERMANOVA
adonis_xen <- adonis2(beta_diversity_xen$Btotal ~ scale(depth) + Mesh + habitat_norm + log(total_reads+1), 
                      data = ecol, strata = ecol$ID, by = "margin")
print(adonis_xen)

# Turnover and nestedness proportions
mean_total <- mean(beta_diversity_xen$Btotal)
mean_turnover <- mean(beta_diversity_xen$Brepl)
mean_nestedness <- mean(beta_diversity_xen$Brich)

per_turnover <- (mean_turnover / mean_total) * 100
per_nestedness <- (mean_nestedness / mean_total) * 100

print(paste("Turnover represents:", round(per_turnover, 2), "%"))
print(paste("Nestedness represents:", round(per_nestedness, 2), "%"))

## Tables for results
table1 <- clean_adonis_table(adonis_xen, "Table 1: Adonis Beta Total - Xenacoelomorpha")
# NOTE: The turnover and nestedness models were not calculated in the code above. 
# Uncomment the following lines only if you calculate 'adonis_xen_turnover' and 'adonis_xen_nestedness'
# table2 <- clean_adonis_table(adonis_xen_turnover, "Table 2: Adonis Turnover")
# table3 <- clean_adonis_table(adonis_xen_nestedness, "Table 3: Adonis Nestedness")

# Save table
save_as_docx(
  table1, 
  # table2, 
  # table3, 
  path = "Results_Adonis_Xen_Tax.docx"
)

pairwise.adonis2(comm_xen ~ Mesh, data = ecol, method = "jaccard")
pairwise.adonis2(comm_xen ~ habitat_norm, data = ecol, method = "jaccard")


## Acari (beta tax) ------------------------------------------------------------

aca_asvs <- species %>%
  filter(Best.group == "Acari") %>%
  pull(ASVid)

comm_aca <- comm[, colnames(comm) %in% aca_asvs]

# Acari richness
if(is.null(dim(comm_aca))) {
  ecol$richness_aca <- (comm_aca > 0) * 1
} else {
  ecol$richness_aca <- rowSums(comm_aca > 0)
}

# Calculate beta diversity (Jaccard)
beta_diversity_aca <- beta(comm_aca, func = "jaccard")  
str(beta_diversity_aca)

## PERMANOVA
adonis_aca <- adonis2(beta_diversity_aca$Btotal ~ scale(depth) + Mesh + habitat_norm + log(total_reads+1), 
                      data = ecol, strata = ecol$ID, by = "margin")
print(adonis_aca)

# Turnover and nestedness proportions
mean_total <- mean(beta_diversity_aca$Btotal)
mean_turnover <- mean(beta_diversity_aca$Brepl)
mean_nestedness <- mean(beta_diversity_aca$Brich)

per_turnover <- (mean_turnover / mean_total) * 100
per_nestedness <- (mean_nestedness / mean_total) * 100

print(paste("Turnover represents:", round(per_turnover, 2), "%"))
print(paste("Nestedness represents:", round(per_nestedness, 2), "%"))

## Tables for results
table1 <- clean_adonis_table(adonis_aca, "Table 1: Adonis Beta Total - Acari")
# NOTE: Same as above, uncomment if turnover/nestedness models are calculated.
# table2 <- clean_adonis_table(adonis_aca_turnover, "Table 2: Adonis Turnover")
# table3 <- clean_adonis_table(adonis_aca_nestedness, "Table 3: Adonis Nestedness")

# Save table
save_as_docx(
  table1, 
  # table2, 
  # table3, 
  path = "Results_Adonis_Aca_Tax.docx"
)

pairwise.adonis2(comm_aca ~ Mesh, data = ecol, method = "jaccard")
pairwise.adonis2(comm_aca ~ habitat_norm, data = ecol, method = "jaccard")


############# B. Phylogenetic beta diversity #####################################

## Meiofauna (beta phyl) -------------------------------------------------------

library(Biostrings)
library(DECIPHER)
library(phangorn)

# Assuming species has 'ASVid' and 'sequence'
seqs <- DNAStringSet(species$sequence)
names(seqs) <- species$ASVid
writeXStringSet(seqs, filepath = "Sequences_meioAntarctica18S.fasta", format = "fasta")

# Read FASTA file from the specified path
seqs <- readDNAStringSet("Sequences_meioAntarctica18S.fasta")

# Inspect the first lines
head(seqs)

## Tree building
alignment <- AlignSeqs(seqs)

# Convert alignment to a phyDat object
phy_data <- phyDat(as.matrix(alignment), type = "DNA")

# Make phylogenetic tree with Neighbor-Joining
dist_matrix <- dist.ml(phy_data)  # Distance matrix
tree <- NJ(dist_matrix)           # Tree with Neighbor-Joining

comm <- (comm > 0) * 1

# Phylogenetic Richness (Alpha) included in this section
pd_sample <- BAT::alpha(comm, tree = tree)
pd_sample <- as.data.frame(pd_sample)
colnames(pd_sample) <- "phylo.diver"
pd_sample$sample_ID <- rownames(pd_sample)

stations <- merge(ecol, pd_sample, by = "sample_ID", all.x = TRUE)
stations$phylo.diver[is.na(stations$phylo.diver)] <- 0

stations$Mesh <- as.factor(stations$Mesh)

## Beta diversity (Phylogenetic)
beta.phylo <- BAT::beta(comm, tree = tree)

## PERMANOVA
adonis_meio_phyl <- adonis2(beta.phylo$Btotal ~ scale(depth) + Mesh + habitat_norm + log(total_reads+1), 
                            data = ecol, strata = ecol$ID, by = "margin")
print(adonis_meio_phyl)

# Turnover and nestedness proportions
mean_total <- mean(beta.phylo$Btotal)
mean_turnover <- mean(beta.phylo$Brepl)
mean_nestedness <- mean(beta.phylo$Brich)

per_turnover <- (mean_turnover / mean_total) * 100
per_nestedness <- (mean_nestedness / mean_total) * 100

print(paste("Turnover represents:", round(per_turnover, 2), "%"))
print(paste("Nestedness represents:", round(per_nestedness, 2), "%"))

## Tables for results
table1 <- clean_adonis_table(adonis_meio_phyl, "Table 1: Adonis Beta Total (Phylogenetic) - Meiofauna")

# Save table
save_as_docx(
  table1, 
  path = "Results_Adonis_Meio_Phyl.docx"
)


## Copepoda (beta phyl) --------------------------------------------------------

species_cop <- species %>% dplyr::filter(Best.group == "Copepoda")

seqs <- DNAStringSet(species_cop$sequence)
names(seqs) <- species_cop$ASVid
writeXStringSet(seqs, filepath = "Sequences_copAntarctica18S.fasta", format = "fasta")

# Alignment and NJ tree
alignment <- AlignSeqs(seqs)
phy_data <- phyDat(as.matrix(alignment), type = "DNA")
dist_matrix <- dist.ml(phy_data)
tree_cop <- NJ(dist_matrix)

## Community Matrix
asvs_cop <- species_cop$ASVid
comm_cop <- comm[, colnames(comm) %in% asvs_cop]

# Presence/Absence
comm_cop <- as.matrix((comm_cop > 0) * 1)

## Phylogenetic richness (PD)
pd_cop_res <- BAT::alpha(comm_cop, tree = tree_cop)
pd_cop <- data.frame(
  sample_ID = rownames(pd_cop_res),
  phylo.diver = as.numeric(pd_cop_res[,1])
)

## Merge with ecol 
stations_cop <- merge(ecol, pd_cop, by = "sample_ID", all.x = TRUE)

# Substitute NAs with 0
stations_cop$phylo.diver[is.na(stations_cop$phylo.diver)] <- 0

## Beta diversity (Phylogenetic)
beta.phylo_cop <- BAT::beta(comm_cop, tree = tree_cop)

## PERMANOVA
adonis_cop_phyl <- adonis2(beta.phylo_cop$Btotal ~ scale(depth) + Mesh + habitat_norm + log(total_reads+1), 
                           data = ecol, strata = ecol$ID, by = "margin")
print(adonis_cop_phyl)

# Turnover and nestedness proportions
mean_total <- mean(beta.phylo_cop$Btotal)
mean_turnover <- mean(beta.phylo_cop$Brepl)
mean_nestedness <- mean(beta.phylo_cop$Brich)

per_turnover <- (mean_turnover / mean_total) * 100
per_nestedness <- (mean_nestedness / mean_total) * 100

print(paste("Turnover represents:", round(per_turnover, 2), "%"))
print(paste("Nestedness represents:", round(per_nestedness, 2), "%"))

## Tables for results
table1 <- clean_adonis_table(adonis_cop_phyl, "Table 1: Adonis Beta Total (Phylogenetic) - Copepoda")

# Save table
save_as_docx(
  table1, 
  path = "Results_Adonis_Cop_Phyl.docx"
)



## Nematoda (beta phyl) --------------------------------------------------------

species_nem <- species %>% dplyr::filter(Best.group == "Nematoda")

seqs <- DNAStringSet(species_nem$sequence)
names(seqs) <- species_nem$ASVid
writeXStringSet(seqs, filepath = "Sequences_nemAntarctica18S.fasta", format = "fasta")

# Alignment and NJ tree
alignment <- AlignSeqs(seqs)
phy_data <- phyDat(as.matrix(alignment), type = "DNA")
dist_matrix <- dist.ml(phy_data)
tree_nem <- NJ(dist_matrix)

## Community Matrix
asvs_nem <- species_nem$ASVid
comm_nem <- comm[, colnames(comm) %in% asvs_nem]

# Presence/Absence
comm_nem <- as.matrix((comm_nem > 0) * 1)

## Phylogenetic richness (PD)
pd_nem_res <- BAT::alpha(comm_nem, tree = tree_nem)
pd_nem <- data.frame(
  sample_ID = rownames(pd_nem_res),
  phylo.diver = as.numeric(pd_nem_res[,1])
)

## Merge with ecol 
stations_nem <- merge(ecol, pd_nem, by = "sample_ID", all.x = TRUE)

# Substitute NAs with 0
stations_nem$phylo.diver[is.na(stations_nem$phylo.diver)] <- 0

## Beta diversity (Phylogenetic)
beta.phylo_nem <- BAT::beta(comm_nem, tree = tree_nem)

## PERMANOVA
adonis_nem_phyl <- adonis2(beta.phylo_nem$Btotal ~ scale(depth) + Mesh + habitat_norm + log(total_reads+1), 
                           data = ecol, strata = ecol$ID, by = "margin")
print(adonis_nem_phyl)

# Turnover and nestedness proportions
mean_total <- mean(beta.phylo_nem$Btotal)
mean_turnover <- mean(beta.phylo_nem$Brepl)
mean_nestedness <- mean(beta.phylo_nem$Brich)

per_turnover <- (mean_turnover / mean_total) * 100
per_nestedness <- (mean_nestedness / mean_total) * 100

print(paste("Turnover represents:", round(per_turnover, 2), "%"))
print(paste("Nestedness represents:", round(per_nestedness, 2), "%"))

## Tables for results
table1 <- clean_adonis_table(adonis_nem_phyl, "Table 1: Adonis Beta Total (Phylogenetic) - nemepoda")

# Save table
save_as_docx(
  table1, 
  path = "Results_Adonis_nem_Phyl.docx"
)



## Platyhelminthes (beta phyl) --------------------------------------------------------

species_plat <- species %>% dplyr::filter(Best.group == "Platyhelminthes")

seqs <- DNAStringSet(species_plat$sequence)
names(seqs) <- species_plat$ASVid
writeXStringSet(seqs, filepath = "Sequences_platAntarctica18S.fasta", format = "fasta")

# Alignment and NJ tree
alignment <- AlignSeqs(seqs)
phy_data <- phyDat(as.matrix(alignment), type = "DNA")
dist_matrix <- dist.ml(phy_data)
tree_plat <- NJ(dist_matrix)

## Community Matrix
asvs_plat <- species_plat$ASVid
comm_plat <- comm[, colnames(comm) %in% asvs_plat]

# Presence/Absence
comm_plat <- as.matrix((comm_plat > 0) * 1)

## Phylogenetic richness (PD)
pd_plat_res <- BAT::alpha(comm_plat, tree = tree_plat)
pd_plat <- data.frame(
  sample_ID = rownames(pd_plat_res),
  phylo.diver = as.numeric(pd_plat_res[,1])
)

## Merge with ecol 
stations_plat <- merge(ecol, pd_plat, by = "sample_ID", all.x = TRUE)

# Substitute NAs with 0
stations_plat$phylo.diver[is.na(stations_plat$phylo.diver)] <- 0

## Beta diversity (Phylogenetic)
beta.phylo_plat <- BAT::beta(comm_plat, tree = tree_plat)

## PERMANOVA
adonis_plat_phyl <- adonis2(beta.phylo_plat$Btotal ~ scale(depth) + Mesh + habitat_norm + log(total_reads+1), 
                           data = ecol, strata = ecol$ID, by = "margin")
print(adonis_plat_phyl)

# Turnover and nestedness proportions
mean_total <- mean(beta.phylo_plat$Btotal)
mean_turnover <- mean(beta.phylo_plat$Brepl)
mean_nestedness <- mean(beta.phylo_plat$Brich)

per_turnover <- (mean_turnover / mean_total) * 100
per_nestedness <- (mean_nestedness / mean_total) * 100

print(paste("Turnover represents:", round(per_turnover, 2), "%"))
print(paste("Nestedness represents:", round(per_nestedness, 2), "%"))

## Tables for results
table1 <- clean_adonis_table(adonis_plat_phyl, "Table 1: Adonis Beta Total (Phylogenetic) - platepoda")

# Save table
save_as_docx(
  table1, 
  path = "Results_Adonis_plat_Phyl.docx"
)



## Ostracoda (beta phyl) --------------------------------------------------------

species_os <- species %>% dplyr::filter(Best.group == "Ostracoda")

seqs <- DNAStringSet(species_os$sequence)
names(seqs) <- species_os$ASVid
writeXStringSet(seqs, filepath = "Sequences_osAntarctica18S.fasta", format = "fasta")

# Alignment and NJ tree
alignment <- AlignSeqs(seqs)
phy_data <- phyDat(as.matrix(alignment), type = "DNA")
dist_matrix <- dist.ml(phy_data)
tree_os <- NJ(dist_matrix)

## Community Matrix
asvs_os <- species_os$ASVid
comm_os <- comm[, colnames(comm) %in% asvs_os]

# Presence/Absence
comm_os <- as.matrix((comm_os > 0) * 1)

## Phylogenetic richness (PD)
pd_os_res <- BAT::alpha(comm_os, tree = tree_os)
pd_os <- data.frame(
  sample_ID = rownames(pd_os_res),
  phylo.diver = as.numeric(pd_os_res[,1])
)

## Merge with ecol 
stations_os <- merge(ecol, pd_os, by = "sample_ID", all.x = TRUE)

# Substitute NAs with 0
stations_os$phylo.diver[is.na(stations_os$phylo.diver)] <- 0

## Beta diversity (Phylogenetic)
beta.phylo_os <- BAT::beta(comm_os, tree = tree_os)

## PERMANOVA
adonis_os_phyl <- adonis2(beta.phylo_os$Btotal ~ scale(depth) + Mesh + habitat_norm + log(total_reads+1), 
                            data = ecol, strata = ecol$ID, by = "margin")
print(adonis_os_phyl)

# Turnover and nestedness proportions
mean_total <- mean(beta.phylo_os$Btotal)
mean_turnover <- mean(beta.phylo_os$Brepl)
mean_nestedness <- mean(beta.phylo_os$Brich)

per_turnover <- (mean_turnover / mean_total) * 100
per_nestedness <- (mean_nestedness / mean_total) * 100

print(paste("Turnover represents:", round(per_turnover, 2), "%"))
print(paste("Nestedness represents:", round(per_nestedness, 2), "%"))

## Tables for results
table1 <- clean_adonis_table(adonis_os_phyl, "Table 1: Adonis Beta Total (Phylogenetic) - osepoda")

# Save table
save_as_docx(
  table1, 
  path = "Results_Adonis_os_Phyl.docx"
)


## Annelida (beta phyl) --------------------------------------------------------

species_ann <- species %>% dplyr::filter(Best.group == "Annelida")

seqs <- DNAStringSet(species_ann$sequence)
names(seqs) <- species_ann$ASVid
writeXStringSet(seqs, filepath = "Sequences_annAntarctica18S.fasta", format = "fasta")

# Alignment and NJ tree
alignment <- AlignSeqs(seqs)
phy_data <- phyDat(as.matrix(alignment), type = "DNA")
dist_matrix <- dist.ml(phy_data)
tree_ann <- NJ(dist_matrix)

## Community Matrix
asvs_ann <- species_ann$ASVid
comm_ann <- comm[, colnames(comm) %in% asvs_ann]

# Presence/Absence
comm_ann <- as.matrix((comm_ann > 0) * 1)

## Phylogenetic richness (PD)
pd_ann_res <- BAT::alpha(comm_ann, tree = tree_ann)
pd_ann <- data.frame(
  sample_ID = rownames(pd_ann_res),
  phylo.diver = as.numeric(pd_ann_res[,1])
)

## Merge with ecol 
stations_ann <- merge(ecol, pd_ann, by = "sample_ID", all.x = TRUE)

# Substitute NAs with 0
stations_ann$phylo.diver[is.na(stations_ann$phylo.diver)] <- 0

## Beta diversity (Phylogenetic)
beta.phylo_ann <- BAT::beta(comm_ann, tree = tree_ann)

## PERMANOVA
adonis_ann_phyl <- adonis2(beta.phylo_ann$Btotal ~ scale(depth) + Mesh + habitat_norm + log(total_reads+1), 
                            data = ecol, strata = ecol$ID, by = "margin")
print(adonis_ann_phyl)

# Turnover and nestedness proportions
mean_total <- mean(beta.phylo_ann$Btotal)
mean_turnover <- mean(beta.phylo_ann$Brepl)
mean_nestedness <- mean(beta.phylo_ann$Brich)

per_turnover <- (mean_turnover / mean_total) * 100
per_nestedness <- (mean_nestedness / mean_total) * 100

print(paste("Turnover represents:", round(per_turnover, 2), "%"))
print(paste("Nestedness represents:", round(per_nestedness, 2), "%"))

## Tables for results
table1 <- clean_adonis_table(adonis_ann_phyl, "Table 1: Adonis Beta Total (Phylogenetic) - annepoda")

# Save table
save_as_docx(
  table1, 
  path = "Results_Adonis_ann_Phyl.docx"
)



## Gastrotricha (beta phyl) --------------------------------------------------------

species_gas <- species %>% dplyr::filter(Best.group == "Gastrotricha")

seqs <- DNAStringSet(species_gas$sequence)
names(seqs) <- species_gas$ASVid
writeXStringSet(seqs, filepath = "Sequences_gasAntarctica18S.fasta", format = "fasta")

# Alignment and NJ tree
alignment <- AlignSeqs(seqs)
phy_data <- phyDat(as.matrix(alignment), type = "DNA")
dist_matrix <- dist.ml(phy_data)
tree_gas <- NJ(dist_matrix)

## Community Matrix
asvs_gas <- species_gas$ASVid
comm_gas <- comm[, colnames(comm) %in% asvs_gas]

# Presence/Absence
comm_gas <- as.matrix((comm_gas > 0) * 1)

## Phylogenetic richness (PD)
pd_gas_res <- BAT::alpha(comm_gas, tree = tree_gas)
pd_gas <- data.frame(
  sample_ID = rownames(pd_gas_res),
  phylo.diver = as.numeric(pd_gas_res[,1])
)

## Merge with ecol 
stations_gas <- merge(ecol, pd_gas, by = "sample_ID", all.x = TRUE)

# Substitute NAs with 0
stations_gas$phylo.diver[is.na(stations_gas$phylo.diver)] <- 0

## Beta diversity (Phylogenetic)
beta.phylo_gas <- BAT::beta(comm_gas, tree = tree_gas)

## PERMANOVA
adonis_gas_phyl <- adonis2(beta.phylo_gas$Btotal ~ scale(depth) + Mesh + habitat_norm + log(total_reads+1), 
                            data = ecol, strata = ecol$ID, by = "margin")
print(adonis_gas_phyl)

# Turnover and nestedness proportions
mean_total <- mean(beta.phylo_gas$Btotal)
mean_turnover <- mean(beta.phylo_gas$Brepl)
mean_nestedness <- mean(beta.phylo_gas$Brich)

per_turnover <- (mean_turnover / mean_total) * 100
per_nestedness <- (mean_nestedness / mean_total) * 100

print(paste("Turnover represents:", round(per_turnover, 2), "%"))
print(paste("Nestedness represents:", round(per_nestedness, 2), "%"))

## Tables for results
table1 <- clean_adonis_table(adonis_gas_phyl, "Table 1: Adonis Beta Total (Phylogenetic) - gasepoda")

# Save table
save_as_docx(
  table1, 
  path = "Results_Adonis_gas_Phyl.docx"
)


## Acari (beta phyl) --------------------------------------------------------

species_aca <- species %>% dplyr::filter(Best.group == "Acari")

seqs <- DNAStringSet(species_aca$sequence)
names(seqs) <- species_aca$ASVid
writeXStringSet(seqs, filepath = "Sequences_acaAntarctica18S.fasta", format = "fasta")

# Alignment and NJ tree
alignment <- AlignSeqs(seqs)
phy_data <- phyDat(as.matrix(alignment), type = "DNA")
dist_matrix <- dist.ml(phy_data)
tree_aca <- NJ(dist_matrix)

## Community Matrix
asvs_aca <- species_aca$ASVid
comm_aca <- comm[, colnames(comm) %in% asvs_aca]

# Presence/Absence
comm_aca <- as.matrix((comm_aca > 0) * 1)

## Phylogenetic richness (PD)
pd_aca_res <- BAT::alpha(comm_aca, tree = tree_aca)
pd_aca <- data.frame(
  sample_ID = rownames(pd_aca_res),
  phylo.diver = as.numeric(pd_aca_res[,1])
)

## Merge with ecol 
stations_aca <- merge(ecol, pd_aca, by = "sample_ID", all.x = TRUE)

# Substitute NAs with 0
stations_aca$phylo.diver[is.na(stations_aca$phylo.diver)] <- 0

## Beta diversity (Phylogenetic)
beta.phylo_aca <- BAT::beta(comm_aca, tree = tree_aca)

## PERMANOVA
adonis_aca_phyl <- adonis2(beta.phylo_aca$Btotal ~ scale(depth) + Mesh + habitat_norm + log(total_reads+1), 
                           data = ecol, strata = ecol$ID, by = "margin")
print(adonis_aca_phyl)

# Turnover and nestedness proportions
mean_total <- mean(beta.phylo_aca$Btotal)
mean_turnover <- mean(beta.phylo_aca$Brepl)
mean_nestedness <- mean(beta.phylo_aca$Brich)

per_turnover <- (mean_turnover / mean_total) * 100
per_nestedness <- (mean_nestedness / mean_total) * 100

print(paste("Turnover represents:", round(per_turnover, 2), "%"))
print(paste("Nestedness represents:", round(per_nestedness, 2), "%"))

## Tables for results
table1 <- clean_adonis_table(adonis_aca_phyl, "Table 1: Adonis Beta Total (Phylogenetic) - acaepoda")

# Save table
save_as_docx(
  table1, 
  path = "Results_Adonis_aca_Phyl.docx"
)



################### FIGURES #############################
############## Figure 1: Sampling map---------------------------------

## 1. Zoomed-in map (satellite)
# Erase duplicated coordinates

ecol_unique <- ecol %>%
  distinct(longitude, latitude, .keep_all = TRUE)

print(ecol_unique)

register_stadiamaps(key = "e629ead7-e2f7-4cbb-9529-d204a2aaf84d")

satellite_map <- get_stadiamap(
  bbox = c(left = 163.80, bottom = -74.80, right = 164.30, top = -74.65), 
  zoom = 10, 
  maptype = "stamen_terrain" # Options: "terrain", "toner", "watercolor"
)

ggmap(satellite_map)

ggmap(satellite_map) +
  geom_point(data = ecol_unique, aes(x = longitude, y = latitude),
             color = "red", size = 1.85, alpha = 0.8) +  # Sampling points
  geom_text(data = ecol_unique, aes(x = longitude, y = latitude, label = ID),
            color = "black", size = 2.5, vjust = -1) +
  labs(title = "Ross Sea sampling points",
       x = "Longitude", y = "Latitude") +
  theme_minimal()



## 2. Map of all Antarctica
### Geographic projection

# Mapa de la Antártida sin proyección especial
antarctica <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(sovereignt == "Antarctica")

ggplot() +
  geom_sf(data = antarctica, fill = "white", color = "black") +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, -60)) +  # Mismo sistema que ggmap
  theme_minimal() +
  labs(title = "Antarctica (Geographic projection)", x = "Longitude", y = "Latitude")




############## Figure 2: Metabarcoding papers comparison--------------------------

# 1. COLORS
# Define custom colors for the papers
paper_colors <- c(
  "Ref3"  = "#d4f0fe", "Ref4"  = "#15adf8", "Ref6"   = "#0697e0",
  "Ref7"  = "#057eb9", "Ref6+7"= "#068fdc", "Ref8"   = "#046493",
  "Ref9"  = "#012334", "Ref10" = "#ffe103", "Ref2"   = "#fbc400",
  "Ref1"  = "#fb8500", "Ref5"  = "#723c00"
)

world_map <- map_data("world")

# 2. NODE PREPARATION (ECOL) 
ecol_raw <- read.csv("stations network.csv", sep = ";") %>%
  mutate(
    lat = as.numeric(gsub(",", ".", latitude)),
    lon = as.numeric(gsub(",", ".", longitude)),
    lat = ifelse(regional_location == "Antarctica", -abs(lat), lat),
    studyid = as.character(studyid)
  ) %>% 
  filter(!is.na(lat))

# Combine coordinates for Ref6 and Ref7
coord_combined <- ecol_raw %>%
  filter(studyid %in% c("Ref6", "Ref7")) %>%
  summarise(studyid = "Ref6+7", lon = mean(lon), lat = mean(lat))

# Create the combined ecological dataset with geometries
ecol_comb <- ecol_raw %>%
  filter(!studyid %in% c("Ref6", "Ref7")) %>%
  bind_rows(coord_combined) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

# Generate nodes with a bit of jitter to prevent overlapping
nodos_comb <- ecol_comb %>%
  st_drop_geometry() %>%
  select(id = studyid, lon, lat) %>%
  distinct(id, .keep_all = TRUE) %>%
  mutate(
    lon_jitter = jitter(lon, amount = 2),
    lat_jitter = jitter(lat, amount = 2)
  )

# 3. READ CLEAN EDGES TABLE 
# Load the manually cleaned CSV and ensure 'from'/'to' are characters
edges_base_limpia <- read.csv("edges_limpios.csv", sep = ";") %>%
  mutate(from = as.character(from), to = as.character(to))


# 4. MAP 1: GLOBAL 

# Multiply edges directly from the clean table based on their weight
edges_global <- edges_base_limpia %>% uncount(weight)

# Create the graph object
graph_geo_comb <- tbl_graph(nodes = nodos_comb, edges = edges_global, directed = FALSE)

set.seed(42)
plot_geo_map_comb <- ggraph(graph_geo_comb, layout = "manual", x = lon_jitter, y = lat_jitter) +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "#e8e8e8", color = NA) +
  geom_edge_fan(color = "#2c3e50", width = 0.1, alpha = 0.5, spread = 1, show.legend = FALSE) +
  geom_node_point(aes(color = id), size = 4, show.legend = FALSE) +
  geom_node_label(aes(label = id, fill = id),
                  repel = TRUE, size = 3.5, fontface = "bold",
                  color = "black", alpha = 0.8, max.overlaps = Inf,
                  show.legend = FALSE) +
  scale_fill_manual(values = paper_colors) +
  scale_color_manual(values = paper_colors) +
  coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-90, 90)) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=10)))

plot_geo_map_comb


# 5. MAP 2: EUROPE ZOOM 

refs_interes <- c("Ref2", "Ref3", "Ref4", "Ref6+7", "Ref8", "Ref9")

# Filter the clean edges table for Europe and THEN multiply by weight
edges_eu_comb <- edges_base_limpia %>%
  filter(from %in% refs_interes & to %in% refs_interes) %>%
  uncount(weight)

# Filter nodes for Europe with adjusted jitter
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
  geom_node_label(aes(label = id, fill = id),
                  repel = TRUE, size = 3.5, fontface = "bold",
                  color = "black", alpha = 0.8, max.overlaps = Inf,
                  show.legend = FALSE) +
  scale_fill_manual(values = paper_colors) +
  scale_color_manual(values = paper_colors) +
  coord_fixed(ratio = 1.3, xlim = c(-2, 20), ylim = c(36, 70)) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b=10))) +
  labs(title = "ASV Connectivity: North Sea and Mediterranean",
       subtitle = "Ref6 and Ref7 combined (Individual ASV threads)")

plot_eu_comb



#### Shared ASVs

# 1. READ METADATA FILE
df_meta <- read.csv("ASVs_metadata100_gaps.csv", sep = ",")

# 2. CLEAN AND COLLAPSE REDUNDANT ASVs
df_meta_clean <- df_meta %>%
  # Group by the two connected studies and the main ASV (ASVid2_1)
  group_by(study1, study2, ASVid2_1) %>%
  summarise(
    # Collapse all matching secondary ASVs into a single string
    ASVs_destino_combinados = paste(unique(ASVid2_2), collapse = " | "),
    
    # Count distinct ASV matches to assess redundancy
    Numero_de_Matches = n_distinct(ASVid2_2),
    
    # Keep taxonomy (collapsing to handle minor variations if any)
    Grupo_Principal = paste(unique(ASVid2_1_bestgroup), collapse = " | "),
    Taxonomia_Detallada = paste(unique(ASVid2_1_sistergroup), collapse = " | "),
    
    .groups = "drop"
  )

# View the first few rows in the console
head(df_meta_clean)

# 3. SAVE THE CLEAN DATASET
write.csv(df_meta_clean, "Metadatos_ASVs_Limpios.csv", row.names = FALSE)


## Novelty of our paper

library(dplyr)
library(stringr)
library(tidyr)

# 1. Define filenames and dataset names
archivos <- c(
  "Data18S_Antarctica_v2.csv", "Data18S_Asinara_v2.csv", "Data18S_Cordier2021_v2.csv", 
  "Data18S_Degenhardt2021_v2.csv", "Data18S_Fonseca2017.csv", "Data18S_Haenel2017_v2.csv", 
  "Data18S_JondeliusAtherton2020.csv", "Data18S_Kapshyna2024_v2.csv", 
  "Data18S_Mazurkiewicz2024.csv", "Data18S_Polinski2019_v2.csv"
)

nombres_datasets <- c(
  "Antarctica", "Asinara", "Cordier2021", "Degenhardt2021", "Fonseca2017", 
  "Haenel2017", "Jondelius2020", "Kapshyna2024", "Mazurkiewicz2024", "Polinski2019"
)

# 2. Initialize lists to store results
lista_global <- list()
lista_phylum <- list()

# 3. UPDATED LOOP (WITH UNIQUE ASV FILTER)
for (i in seq_along(archivos)) {
  
  datos <- read.csv(archivos[i], sep = ";", stringsAsFactors = FALSE)
  
  datos_procesados <- datos %>%
    # Filter for correct assignments and meiofauna only
    filter(Eval.result == "correct" & Tax.eco.meio == "TRUE") %>%
    # KEY STEP: Keep only unique ASVs
    distinct(ASVid, .keep_all = TRUE) %>%
    mutate(
      similitud = as.numeric(str_split_i(label, "\\|", 4)),
      phylum = str_split_i(label, "\\|", 1)
    )
  
  # Global table calculating <97% and <95% thresholds
  res_global <- datos_procesados %>%
    summarise(
      Dataset = nombres_datasets[i],
      Total_ASVs_Unicos = n(), 
      Porc_menor_97 = mean(similitud < 97.0, na.rm = TRUE) * 100,
      Porc_menor_95 = mean(similitud < 95.0, na.rm = TRUE) * 100
    )
  
  lista_global[[i]] <- res_global
  
  # Table by phylum calculating <97% and <95% thresholds
  res_phylum <- datos_procesados %>%
    filter(!is.na(phylum) & phylum != "") %>%
    group_by(phylum) %>%
    summarise(
      Total_ASVs_Unicos = n(),
      Porc_menor_97 = mean(similitud < 97.0, na.rm = TRUE) * 100,
      Porc_menor_95 = mean(similitud < 95.0, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    mutate(Dataset = nombres_datasets[i]) %>%
    select(Dataset, phylum, Total_ASVs_Unicos, Porc_menor_97, Porc_menor_95)
  
  lista_phylum[[i]] <- res_phylum
}

# 4. Combine all tables
tabla_global_final <- bind_rows(lista_global)
tabla_phylum_final <- bind_rows(lista_phylum)

# 5. CREATE WIDE MATRICES
tabla_ancha_97 <- tabla_phylum_final %>%
  select(Dataset, phylum, Porc_menor_97) %>%
  pivot_wider(names_from = Dataset, values_from = Porc_menor_97, values_fill = 0)

tabla_ancha_95 <- tabla_phylum_final %>%
  select(Dataset, phylum, Porc_menor_95) %>%
  pivot_wider(names_from = Dataset, values_from = Porc_menor_95, values_fill = 0)

# --- VISUALIZATION ---
print("=== GLOBAL SUMMARY (MEIOFAUNA - UNIQUE ASVs) ===")
print(tabla_global_final)

# --- EXPORT TO CSV ---
write.csv(tabla_global_final, "Resumen_Global_Meiofauna_Unicos_95_97.csv", row.names = FALSE)
write.csv(tabla_ancha_97, "Matriz_Phylum_Papers_Unicos_Menor97.csv", row.names = FALSE)
write.csv(tabla_ancha_95, "Matriz_Phylum_Papers_Unicos_Menor95.csv", row.names = FALSE)



## Taxa counts

# 1. Map filenames to their respective Ref IDs 
# (This avoids having to manually read and create Ref1, Ref2, etc., one by one)
dataset_files <- c(
  "Ref1"  = "Data18S_Antarctica_v2.csv",
  "Ref2"  = "Data18S_Asinara_v2.csv",
  "Ref3"  = "Data18S_Cordier2021_v2.csv",
  "Ref4"  = "Data18S_Degenhardt2021_v2.csv",
  "Ref5"  = "Data18S_Fonseca2017.csv",
  "Ref6"  = "Data18S_Haenel2017_v2.csv",
  "Ref7"  = "Data18S_JondeliusAtherton2020.csv",
  "Ref8"  = "Data18S_Kapshyna2024_v2.csv",
  "Ref9"  = "Data18S_Mazurkiewicz2024.csv",
  "Ref10" = "Data18S_Polinski2019_v2.csv"
)

# Load all datasets directly into a list
lista_Refs <- lapply(dataset_files, read.csv, sep = ";")

# 2. Filter datasets and bind them into a single master dataframe
Ref_maestra <- bind_rows(
  lapply(lista_Refs, function(Ref) {
    
    # Check if 'rare.asvs' column exists and filter if it does
    if("rare.asvs" %in% colnames(Ref)) {
      Ref <- Ref %>% filter(rare.asvs %in% c(FALSE, "FALSE", "false", "False", 0))
    }
    
    Ref %>% 
      filter(
        # Filter only meiofauna
        Tax.eco.meio %in% c(TRUE, "TRUE", "true", "True"),
        # Only correct assignments
        Eval.result == "correct"
      ) %>%
      # Crucial step! Remove duplicate ASVs, keeping the first occurrence
      distinct(ASVid, .keep_all = TRUE) %>%
      # Select relevant columns
      select(ASVid, Best.group, nreads) 
  }),
  .id = "dataset"
)

# 3. Calculate total reads per dataset
total_reads_Ref <- Ref_maestra %>%
  group_by(dataset) %>%
  summarise(count = sum(nreads, na.rm = TRUE), .groups = "drop") %>%
  mutate(taxa = "total_reads") %>%
  pivot_wider(names_from = dataset, values_from = count, values_fill = 0)

# 4. Calculate ASVs by taxonomic group
asvs_por_grupo_Ref <- Ref_maestra %>%
  filter(nreads > 0) %>% 
  group_by(dataset, Best.group) %>%
  summarise(count = n_distinct(ASVid), .groups = "drop") %>%
  rename(taxa = Best.group) %>%
  pivot_wider(names_from = dataset, values_from = count, values_fill = 0) %>%
  arrange(taxa)

# 5. Unify and export
final_table <- bind_rows(total_reads_Ref, asvs_por_grupo_Ref)

write_csv2(final_table, "taxa_counts_meiofauna.csv")

# View the result
print(final_table)



## HEATMAP WITH METABARCODING PAPERS 

# 1. CUSTOM COLORS (Re-declared here for self-contained execution)
paper_colors_heatmap <- c(
  "Ref3"  = "#d4f0fe", "Ref4"  = "#15adf8", "Ref6"   = "#0697e0",
  "Ref7"  = "#057eb9", "Ref8"  = "#046493", "Ref9"   = "#012334",
  "Ref10" = "#ffe103", "Ref2"  = "#fbc400", "Ref1"   = "#fb8500",
  "Ref5"  = "#723c00"
)

# 2. LOAD AND PREPARE DATA
df_raw <- read.csv("taxa_counts_meiofauna.csv", sep = ";")

# Remove the total_reads row (not used for this calculation)
df_taxa <- df_raw %>% filter(taxa != "total_reads")

# Convert to long format
df_long <- df_taxa %>%
  pivot_longer(cols = -taxa, names_to = "Ref", values_to = "ASV_count") %>%
  mutate(ASV_count = as.numeric(ASV_count))

# Calculate the ACTUAL TOTAL ASVs for each reference
totals_df <- df_long %>%
  group_by(Ref) %>%
  summarise(Total_ASVs = sum(ASV_count, na.rm = TRUE), .groups = "drop")

# 3. NORMALIZE DATA (Percentages)
df_final <- df_long %>%
  left_join(totals_df, by = "Ref") %>%
  mutate(
    # Calculate the percentage of total ASVs
    percentage = (ASV_count / Total_ASVs) * 100,
    # Set 0 to NA for transparency (white in heatmap)
    alpha_val = na_if(percentage, 0),
    # Create text label with rounded percentage
    label_text = ifelse(is.na(alpha_val), "", paste0(round(percentage, 1)))
  )

# 4. CREATE ABSOLUTE TOTALS ROW
totals_row <- totals_df %>%
  mutate(
    taxa = "Total ASVs",
    ASV_count = Total_ASVs,
    percentage = NA,
    alpha_val = NA, # Since it's NA, R won't color it (it will remain white)
    label_text = as.character(Total_ASVs) # Here we only display the absolute number
  )

# Bind the normal data with the totals row
df_plot <- bind_rows(df_final, totals_row)

# 5. ORDER FACTORS (For X and Y axes)
orden_refs <- c("Ref1", "Ref5", "Ref2", "Ref10", "Ref3", "Ref4", "Ref6", "Ref7", "Ref8", "Ref9")
df_plot$Ref <- factor(df_plot$Ref, levels = orden_refs)

# Order the Y axis alphabetically, but force "Total ASVs" to the end
taxa_levels <- c(sort(unique(df_taxa$taxa)), "Total ASVs")
df_plot$taxa <- factor(df_plot$taxa, levels = taxa_levels)

# 6. PLOT HEATMAP
plot_heatmap <- ggplot(df_plot, aes(x = Ref, y = taxa)) +
  
  # Color layer (now based on percentage)
  geom_tile(aes(fill = Ref, alpha = alpha_val), color = "white", linewidth = 0.5) +
  
  # Text layer (Percentages on top, absolute numbers on the last row)
  geom_text(aes(label = label_text), color = "black", size = 3) +
  
  scale_fill_manual(values = paper_colors_heatmap, guide = "none") + 
  
  # Transparency scale (Linear, no log needed since it's 0-100%)
  scale_alpha_continuous(
    range = c(0.2, 1), 
    trans = "sqrt",
    name = "Relative Richness\n(% of ASVs)",
    na.value = 0 # This makes zeros and the Totals row white
  ) +
  
  # Reverse Y axis so it reads top-to-bottom
  scale_y_discrete(limits = rev) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1, 
      face = "bold",
      color = paper_colors_heatmap[levels(df_plot$Ref)]
    ),
    # Make the "Total ASVs" label bold on the Y axis for emphasis
    axis.text.y = element_text(
      size = 10, 
      face = ifelse(rev(levels(df_plot$taxa)) == "Total ASVs", "bold", "plain")
    ),
    panel.grid = element_blank(), 
    panel.background = element_rect(fill = "#fdfdfd", color = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b=15))
  ) +
  labs(
    title = "ASV Composition by Taxonomic Group",
    x = "Reference",
    y = "Taxonomic Group"
  )

plot_heatmap

### Figure 4: Phylogenetic tree ---------------------------------------------------

# Function: dataframe to FASTA

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


# Export sequences to FASTA
dataframe2fas(
  species[, c("ASVid", "sequence")],
  file = "sequences.fasta"
)

# Sequence alignment
seqs <- readDNAStringSet("sequences.fasta")

alignment <- msa(
  seqs,
  method = "Muscle",
  type = "dna"
)

alignment_phydat <- msaConvert(       # this step took several hours
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


# Optional: ultrametric tree
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

# 1. Convert the community matrix to long format
# We do this globally just once to make the code cleaner and faster
comm_long <- comm %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample_ID") %>%
  pivot_longer(-sample_ID, names_to = "ASVid", values_to = "count") %>%
  filter(count > 0)  # Filter to keep only present ASVs

# Define taxon levels to maintain a consistent order across all plots
taxon_levels <- c(
  "Xenacoelomorpha", "Priapulida", "Nematoda", "Tardigrada", "Acari", "Ostracoda", "Copepoda",
  "Gnathostomulida", "Rotifera", "Gastrotricha", "Platyhelminthes", "Annelida"
)

# Assign colors to each taxon according to the desired groups
taxon_colors <- c(
  "Xenacoelomorpha" = "white",      # Group 1 (Xenacoelomorpha)
  "Priapulida"      = "#d4f0fe",    # Group 2 (Ecdysozoa)
  "Nematoda"        = "#15adf8",    # Group 2 (Ecdysozoa)
  "Tardigrada"      = "#0697e0",    # Group 2 (Ecdysozoa)
  "Acari"           = "#057eb9",    # Group 2 (Ecdysozoa)
  "Ostracoda"       = "#046493",    # Group 2 (Ecdysozoa)
  "Copepoda"        = "#012334",    # Group 2 (Ecdysozoa)
  "Gnathostomulida" = "#ffe103",    # Group 3 (Lophotrochozoa)
  "Rotifera"        = "#fbc400",    # Group 3 (Lophotrochozoa)
  "Gastrotricha"    = "#fb8500",    # Group 3 (Lophotrochozoa)
  "Platyhelminthes" = "#d47000",    # Group 3 (Lophotrochozoa)
  "Annelida"        = "#723c00"     # Group 3 (Lophotrochozoa)
)


## A. Mesh Size 

# Merge ecological and species info
df_mesh <- comm_long %>%
  inner_join(species, by = "ASVid") %>%
  inner_join(ecol, by = "sample_ID") %>%
  filter(!is.na(Mesh)) %>%                  # Remove NAs in Mesh
  filter(Mesh %in% c(20, 50, 100, 200)) %>% # Keep only target mesh sizes
  group_by(Mesh, Best.group) %>%
  summarise(ASV_count = n(), .groups = "drop")

# Order the Best.group factor levels
df_mesh$Best.group <- factor(df_mesh$Best.group, levels = taxon_levels)
ecol$Mesh <- as.factor(ecol$Mesh)

# Create the stacked horizontal bar chart
ggplot(df_mesh, aes(x = factor(Mesh), y = ASV_count, fill = Best.group)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_fill_manual(values = taxon_colors) +
  labs(
    x = "Mesh size",
    y = "Proportion of ASVs",
    fill = "Taxon",
    title = "Proportion of ASVs by Taxon and Mesh Size"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


## B. Habitat Type 

# Merge ecological and species info
df_habitat <- comm_long %>%
  inner_join(species, by = "ASVid") %>%
  inner_join(ecol, by = "sample_ID") %>%
  filter(!is.na(habitat_norm)) %>%          # Remove NAs in habitat
  group_by(habitat_norm, Best.group) %>%
  summarise(ASV_count = n(), .groups = "drop")

# Order factor levels for Best.group and habitat_norm
df_habitat$Best.group <- factor(df_habitat$Best.group, levels = taxon_levels)

df_habitat$habitat_norm <- factor(df_habitat$habitat_norm, levels = c(
  "epilithic", "organic", "spicule", "gravel", "sand", "silt"
))

# Create the stacked horizontal bar chart
ggplot(df_habitat, aes(x = habitat_norm, y = ASV_count, fill = Best.group)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_fill_manual(values = taxon_colors) +
  labs(
    x = "Habitat type",
    y = "Proportion of ASVs",
    fill = "Taxon",
    title = "Proportion of ASVs by Taxon and Habitat"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


## C. Depth 

# Merge ecological and species info
df_depth <- comm_long %>%
  inner_join(species, by = "ASVid") %>%
  inner_join(ecol, by = "sample_ID") %>%
  filter(!is.na(depth)) %>%                 # Remove NAs in depth
  group_by(depth, Best.group) %>%
  summarise(ASV_count = n(), .groups = "drop")

# Order the Best.group factor levels
df_depth$Best.group <- factor(df_depth$Best.group, levels = taxon_levels)

# Plot 1: Stacked bar chart (Proportions)
ggplot(df_depth, aes(x = depth, y = ASV_count, fill = Best.group)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_fill_manual(values = taxon_colors) +
  labs(
    x = "Depth (m)",
    y = "Proportion of ASVs",
    fill = "Taxon",
    title = "Proportion of ASVs by Taxon and Depth"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# Plot 2: Stacked area chart (Absolute counts / Richness)
ggplot(df_depth, aes(x = depth, y = ASV_count, fill = Best.group)) +
  geom_area(position = "stack", color = "black") +
  scale_fill_manual(values = taxon_colors) +
  labs(
    x = "Depth (m)",
    y = "Richness of ASVs",
    fill = "Taxon",
    title = "Richness of ASVs across Depth"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# Plot 3: Stacked area chart (Proportions)

# 1. Prepare data by calculating proportions for each depth level
df_prop <- df_depth %>%
  group_by(depth) %>%
  mutate(proportion = ASV_count / sum(ASV_count)) %>%
  ungroup()

# 2. Stacked Area Plot (Proportional)
ggplot(df_prop, aes(x = depth, y = proportion, fill = Best.group)) +
  # Use geom_area with position = "fill" to ensure it reaches 100%
  geom_area(position = "fill", alpha = 0.9) + 
  # Optional: add subtle lines between areas
  geom_line(position = "fill", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = taxon_colors) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(
    x = "Depth (m)",
    y = "Proportion of ASVs (%)",
    fill = "Taxon",
    title = "Relative Taxonomic Composition across Depth"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )



### Figure 6: Venn diagram ----------------------------------------------------------

# Create a temporary dataframe to avoid modifying the original 'comm' matrix.
# We append the 'Mesh' variable from the 'ecol' dataset.
comm_venn <- comm %>%
  as.data.frame() %>%
  mutate(Mesh = factor(ecol$Mesh, levels = c("20", "50", "100", "200")))

# Group by Mesh and detect presence of ASVs 
# (Using > 0 is safer than == 1 in case of abundance matrices)
asv_sets <- comm_venn %>%
  group_by(Mesh) %>%
  summarise(across(where(is.numeric), ~ any(. > 0)), .groups = "drop") %>%
  as.data.frame()

# Generate a list of ASVs present for each mesh size
set_list <- lapply(split(asv_sets[, -1], asv_sets$Mesh), function(x) {
  colnames(x)[which(as.logical(x))]
})

# Ensure the names in set_list are explicitly ordered
set_list <- set_list[c("20", "50", "100", "200")]

# Define colors by mesh size (ordered from lightest to darkest)
mesh_colors <- c(
  "20"  = "#cce4f6",
  "50"  = "#66b3e6",
  "100" = "#1f78b4",
  "200" = "#012334"
)

# Clear any previous plots to prevent overlapping in the viewing window
grid.newpage()

# Draw the Venn diagram with the desired colors and order
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

# Render the plot
grid.draw(venn.plot)



### Figure S1 ----------------------------------------------------------

# 1. READ THE FILE
# We use check.names = FALSE so R doesn't replace spaces/parentheses with dots
df_s6 <- read.csv("Table_S6.csv", sep = ";", stringsAsFactors = FALSE, check.names = FALSE)

# Rename the first column generically to "Variable" for easier filtering
colnames(df_s6)[1] <- "Variable"

# 2. EXTRACT AND CLEAN THE DATA
# Find the exact row containing the global percentages and transform the table
df_percentages <- df_s6 %>%
  # Filter the row that contains the percentage text (ignoring the "by group" rows)
  filter(
    str_detect(Variable, "Percentage of ASVs with <95% identity") & 
      !str_detect(Variable, "by group")
  ) %>%
  # Pivot the study columns into rows (long format)
  pivot_longer(cols = -Variable, names_to = "Study", values_to = "Percentage") %>%
  # Ensure the percentage is numeric
  mutate(Percentage = as.numeric(Percentage))

# 3. SEPARATE "THIS STUDY" FROM THE REST
# Extract the value for your study (looks for any column containing "This study")
our_study_val <- df_percentages %>%
  filter(str_detect(Study, "(?i)This study")) %>%
  pull(Percentage)

# Extract the values for the other published studies
published_vals <- df_percentages %>%
  filter(!str_detect(Study, "(?i)This study")) %>%
  pull(Percentage)

# 4. STATISTICAL CALCULATIONS (Gaussian Curve)
mean_val <- mean(published_vals, na.rm = TRUE)
sd_val <- sd(published_vals, na.rm = TRUE)

# Generate the coordinates for the curve line
x_val <- seq(5, 85, length.out = 200)
y_val <- dnorm(x_val, mean = mean_val, sd = sd_val)

# Prepare the exact dataframes for ggplot
df_curve <- data.frame(x = x_val, y = y_val)
df_published <- data.frame(x = published_vals, y = dnorm(published_vals, mean = mean_val, sd = sd_val))
df_our <- data.frame(x = our_study_val, y = dnorm(our_study_val, mean = mean_val, sd = sd_val))

# 5. DRAW THE FIGURE
ggplot() +
  # Main curve
  geom_line(data = df_curve, aes(x = x, y = y), linewidth = 1) +
  # Red dots (previous studies)
  geom_point(data = df_published, aes(x = x, y = y), color = "red", size = 4) +
  # Blue triangle (your study)
  geom_point(data = df_our, aes(x = x, y = y), shape = 17, color = "blue", size = 6) +
  
  # Aesthetics and labels
  theme_bw(base_size = 14) +
  labs(
    x = "percentage",
    y = "density"
  ) +
  theme(
    panel.grid.minor = element_blank()
  )