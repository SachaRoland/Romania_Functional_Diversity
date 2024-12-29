
# PCAs
library(vegan)

#

# FSC variables and species USING RODENT FUNCTIONAL UNITS ! 

# Convert factors to numeric as needed by PCA. 
traits_matrix_updated$ActivityPattern <- as.numeric(traits_matrix_updated$ActivityPattern)
traits_matrix_updated$TrophicLevel <- as.numeric(traits_matrix_updated$TrophicLevel)

trait_pca <- rda(traits_matrix_updated, scale = TRUE)  # Vegan package
str(traits_matrix_updated)
biplot(trait_pca, scaling = 2)

# Extract PCA loadings
loadings <- scores(trait_pca, display = "species")
print(loadings)

# View loadings for PC1
loadings_PC1 <- loadings[, 1]  # The first column corresponds to PC1
print(loadings_PC1)


# Import abundances
setwd("/Users/sacharoland/Desktop/WUR/THESIS/FSC and vertebrates- Romania/Data Analysis/Data")
RAI <- read.csv("All_densities.csv")
colnames(RAI)
RAI <- RAI %>% select(-"X", ,
                          -"Glis_glis" ) # no real abundance for glis.

# Adapt plotname columns and save names in seprate object.  
rownames(RAI) <- RAI$Plotname
RAI <- RAI[, -which(names(RAI) == "Plotname")]
print(RAI)
colnames(RAI)

# aggregates Apedomus species abudance. 
# Add a new column "Apodemus_combined" by summing "Apodemus_sp_" and "Apodemus_agrarius"
RAI$Apodemus_sp. <- RAI$Apodemus_sp_ + RAI$Apodemus_agrarius

# Optionally, you can remove the original columns if you no longer need them
RAI <- RAI %>% select(-Apodemus_sp_, -Apodemus_agrarius)

# View the updated data frame
head(RAI)


#Ensure the RAI_noR data frame is fully numeric
str(RAI)

RAI <- RAI %>%
  mutate(
    Martes = rowSums(across(c(Martes, Martes_martes, Martes_foina)), na.rm = TRUE)  # Combine values
  ) %>%
  dplyr::select(-Martes_foina, -Martes_martes)  # Drop the unused columns

# Standardize column names in species_only for matching
colnames(RAI) <- colnames(RAI) %>%
  gsub("_", " ", .)  # Replace underscores with spaces


# calculate the weighted trait values for each plot
# # Find common species
# common_species <- intersect(colnames(RAI_noR), rownames(traits_matrix_noR))
#
# # Subset both datasets to include only common species
# RAI_noR <- RAI_noR[, common_species, drop = FALSE]
# traits_matrix_noR <- traits_matrix_noR[common_species, , drop = FALSE]
#
# #Multiply RAI_noR (plots × species) with traits_matrix_noR (species × traits)
#
# # Convert RAI_noR and traits_matrix_noR to matrices
# RAI_noR_matrix <- as.matrix(RAI_noR)
# traits_matrix_noR_matrix <- as.matrix(traits_matrix_noR)
#
# # Calculate weighted mean trait values for each plot
# # which represents the weighted mean trait values for each plot based on species densities
# plot_trait_values <- RAI_noR_matrix %*% traits_matrix_noR_matrix / rowSums(RAI_noR_matrix)
# # Check the structure of the output
# str(plot_trait_values)
# print(plot_trait_values)

# # View results
# print(plot_trait_values)
#
# FD_noR_summary$PC1
# plot_trait_values

# #Correlate Weighted Trait Values with PC1 ####
# # Ensure dimensions match before computing correlation
# if (nrow(plot_trait_values) == length(FD_noR_summary$PC1)) {
#   # Calculate correlation for each trait
#   cor_PC1_traits_noR <- sapply(1:ncol(plot_trait_values), function(i) {
#     cor(plot_trait_values[, i], FD_noR_summary$PC1, method = "pearson")
#   })
#   
#   # Assign trait names to the correlation results
#   names(cor_PC1_traits_noR) <- colnames(plot_trait_values)
#   
#   # View results
#   print(cor_PC1_traits_noR)
# } else {
#   stop("Row dimensions of plot_trait_values and FD_noR_summary$PC1 do not match.")
# }
# 
# # values close to 1 or -1 suggest strong relationship. 
# barplot(loadings_PC1, main = "Trait Contributions to PC1", las = 2, col = "blue")
# barplot(cor_PC1_traits_noR, main = "Correlation of Traits with PC1", las = 2, col = "green")


# RUN FSC-PCA.R FIRST ####
# Include both traits and sites (plots) in the PCA biplot, combining forest structural 
# complexity (FSC) variables and traits.

# Load required libraries
library(vegan)
library(ggplot2)
library(dplyr)
library(grid)

# STEP 1: Prepare Data -------------------------------------------------------

# Combine FSC variables (pca_data_scaled) and traits (plot_trait_values)
fsc_species_data <- cbind(pca_data_scaled, RAI)
#view(fsc_species_data)

# Ensure all data is numeric and scaled
fsc_species_data_scaled <- scale(fsc_species_data)

# STEP 2: Perform PCA --------------------------------------------------------

# Perform PCA on the combined dataset
pca_combined <- rda(fsc_species_data_scaled)

# Summary of PCA to get variance explained
pca_summary <- summary(pca_combined)

# STEP 3: Extract Scores -----------------------------------------------------

# Extract plot scores (site positions)
site_scores <- scores(pca_combined, display = "sites", scaling = 2)
site_scores <- as.data.frame(site_scores)
site_scores$Plotname <- plot_names  # Add plot names for labeling

# Extract variable scores (FSC variables and traits)
variable_scores <- scores(pca_combined, display = "species", scaling = 2)
arrow_data <- as.data.frame(variable_scores[, 1:2])
colnames(arrow_data) <- c("PC1", "PC2")
arrow_data$Variable <- rownames(variable_scores)  # Add variable names

# STEP 4: Select Top 10 FSC Variables ----------------------------------------

# Calculate contributions of variables to PC1 and PC2
arrow_data$Contribution_PC1 <- arrow_data$PC1^2 / sum(arrow_data$PC1^2)
arrow_data$Contribution_PC2 <- arrow_data$PC2^2 / sum(arrow_data$PC2^2)
arrow_data$Total_Contribution <- arrow_data$Contribution_PC1 + arrow_data$Contribution_PC2

# Select the top 10 FSC variables
top_fsc_variables <- arrow_data %>%
  filter(!Variable %in% colnames(RAI)) %>%  # Exclude traits
  arrange(desc(Total_Contribution)) %>%
  slice_head(n = 10)

# Filter FSC variables and include all traits
arrow_data_filtered <- arrow_data %>%
  filter(Variable %in% top_fsc_variables$Variable | Variable %in% colnames(RAI))

# STEP 5: Add Plot Groups (e.g., Old-Growth) ---------------------------------

# Define old-growth plot names
old_growth_plots <- c("UpperMoria", "Rivendell", "Isengard", 
                      "WolfTracks", "SalamanderValley", 
                      "UpperArcer", "LowerMoria")

# Assign categories based on plot names
site_scores <- site_scores %>%
  mutate(Group = ifelse(Plotname %in% old_growth_plots, 
                        "old_growth", 
                        "comparison"))

# STEP 6: Differentiate FSC Variables and Species -----------------------------

# Label variable types
arrow_data_filtered$Type <- ifelse(arrow_data_filtered$Variable %in% colnames(RAI), 
                                   "Species", 
                                   "FSC Variable")

# STEP 7: Create Biplot with ggplot2 -----------------------------------------
# install.packages("ggforce")
# install.packages("ggrepel")
library(ggforce)
library(ggrepel)


plot_FSC_Species <- ggplot() +
  # Add convex hulls for groups
  #geom_mark_hull(data = site_scores, aes(x = PC1, y = PC2, fill = Group, label = Group),
  #concavity = 3, alpha = 0.2, label.size = 0, show.legend = FALSE) +
  # Plot site points with transparency
  geom_point(data = site_scores, aes(x = PC1, y = PC2, color = Group), size = 2, alpha = 0.7) +
  # Add plot names with repelled labels
  geom_text_repel(data = site_scores, aes(x = PC1, y = PC2, label = Plotname), 
                  size = 2.5, color = "black", max.overlaps = 15) +
  # Add FSC variable and trait arrows, scaled
  geom_segment(data = arrow_data_filtered,
               aes(x = 0, y = 0, xend = PC1 * 0.8, yend = PC2 * 0.8, color = Type),
               arrow = arrow(length = unit(0.3, "cm")), linewidth = 1) +
  # Add variable names with repelled labels
  geom_text_repel(data = arrow_data_filtered,
                  aes(x = PC1, y = PC2, label = Variable, color = Type),
                  size = 3, max.overlaps = 15) +
  # Titles and axis labels
  labs(title = "PCA Biplot with Improved Readability",
       x = paste0("PC1 (", round(pca_summary$cont$importance[2, 1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(pca_summary$cont$importance[2, 2] * 100, 2), "%)"),
       color = "Variable Type / Group") +
  # Customize colors
  scale_color_manual(values = c("FSC Variable" = "blue", "Species" = "red", 
                                "old_growth" = "green", "comparison" = "purple")) +
  scale_fill_manual(values = c("old_growth" = "green", "comparison" = "purple")) +
  # Minimal theme with adjustments
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "right") +
  # Optional: Zoom into main cluster
  coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5))

print(plot_FSC_Species)

#ggsave("PCA FSC-SPECIES_plot_A4.png", plot_FSC_Species, width = 8.27, height = 11.69, units = "in", dpi = 300)



#-------------------------------------------
# Perform PCA on species abundance data to visualize species  clustering along PC1:###

# Variance of traits
plot_RAI_values_scaled <- scale(RAI)
RAI_variance <- apply(RAI, 2, var)
print(plot_RAI_values_scaled)
view(plot_RAI_values_scaled)

# Add PC1 and PC2 to FD summary 
FD_summary$PC1_FSC <- PC1_scores # From PCA for FSC --> so the axis expaling mos tof the variance 
FD_summary$PC2_FSC <- PC2_scores

cor_PC1_traits <- cor(FD_summary$PC1_FSC, RAI, method = "pearson")
print(cor_PC1_traits)

cor_PC2_traits <- cor(FD_summary$PC2_FSC, RAI, method = "pearson")
print(cor_PC2_traits)

combined_data <- cbind( RAI, FD_summary[, c("PC1_FSC", "PC2_FSC")])
combined_data_scaled <- scale(combined_data)

# Split into response and predictor variables
response_scaled_RAI <- combined_data_scaled[, colnames(RAI)]


# Threshold for filtering RAI
#threshold <- 0.2
# Identify informative traits
# informative_RAI <- colnames(RAI_noR)[apply(abs(cor_PC1_traits), 2, max) > threshold |
#                                                     apply(abs(cor_PC2_traits), 2, max) > threshold]

informative_RAI <- colnames(RAI)[apply(abs(cor_PC1_traits), 2, max) |
                                       apply(abs(cor_PC2_traits), 2, max)]
# Subset the trait data
#filtered_response_scaled_RAI <- response_scaled_RAI [, informative_RAI]

predictors_scaled <- combined_data_scaled[, c("PC1_FSC", "PC2_FSC")]
# Split into response and predictor variables
response_scaled_RAI <- combined_data_scaled[, colnames(RAI)]

# Rerun RDA
RAI_rda <- rda(response_scaled_RAI  ~ PC1_FSC + PC2_FSC, data = as.data.frame(predictors_scaled))

# Extract the summary of the RDA object
rda_summary <- summary(RAI_rda)


# Result Table ####
library(flextable)
# STEP 1: Prepare Data for the Table -----------------------------------------

# Partitioning of Variance
variance_partitioning <- data.frame(
  Component = c("Total", "Constrained", "Unconstrained"),
  Inertia = c(rda_summary$tot.chi, rda_summary$constr.chi, rda_summary$unconst.chi),
  Proportion = c(1, rda_summary$constr.chi / rda_summary$tot.chi, rda_summary$unconst.chi / rda_summary$tot.chi)
)

# Eigenvalues and Variance Explained
eigenvalues <- as.data.frame(t(rda_summary$cont$importance))
eigenvalues <- tibble::rownames_to_column(eigenvalues, "Component")

# STEP 2: Combine Data for APA-Style Table -----------------------------------
# Combine Partitioning and Eigenvalues for a single summary table
summary_table <- list(
  "Partitioning of Variance" = variance_partitioning,
  "Eigenvalues and Variance Explained" = eigenvalues
)

# STEP 3: Create an APA-Styled Table -----------------------------------------

# Create flextable for Partitioning of Variance
variance_table <- flextable(variance_partitioning) %>%
  set_header_labels(
    Component = "Component",
    Inertia = "Inertia",
    Proportion = "Proportion"
  ) %>%
  colformat_double(digits = 4) %>%
  add_footer_lines("Partitioning of variance in the RDA model.") %>%
  autofit()

# Create flextable for Eigenvalues
eigenvalues_table <- flextable(eigenvalues) %>%
  set_header_labels(
    Component = "Component",
    RDA1 = "RDA1",
    RDA2 = "RDA2",
    PC1 = "PC1",
    PC2 = "PC2"
  ) %>%
  colformat_double(digits = 4) %>%
  add_footer_lines("Eigenvalues and proportion of variance explained by the components.") %>%
  autofit()

print(variance_table)
print(eigenvalues_table)


# # STEP 4: Export the Tables to Word  ----------------------------------
# # Save Partitioning Table
# save_as_docx(variance_table, path = "RDA_Species_Partitioning_Table.docx")
# # Save Eigenvalues Table
# save_as_docx(eigenvalues_table, path = "RDA_Species_Eigenvalues_Table.docx")


# Calculate the percentage of variance explained by each axis
variance_explained <- rda_summary$cont$importance[2, ] * 100  # Proportion explained * 100

# Plot with custom axis labels showing % variance explained
plot(RAI_rda, scaling = 2, 
     main = "Species Clustering with FSC",
     display = c("bp", "sites"),
     xlab = paste0("RDA1 (", round(variance_explained[1], 2), "%)"),
     ylab = paste0("RDA2 (", round(variance_explained[2], 2), "%)"))

# Add custom arrows for species
arrows(0, 0, 
       scores(RAI_rda, display = "species")[, 1] * 3, 
       scores(RAI_rda, display = "species")[, 2] * 3, 
       col = "red", length = 0.1)

# Add custom labels for species
text(scores(RAI_rda, display = "species")[, 1] * 3, 
     scores(RAI_rda, display = "species")[, 2] * 3, 
     labels = rownames(scores(RAI_rda, display = "species")), 
     col = "red", pos = 4, cex = 0.8)


#Explore correlations between PC1 (FSC) and species:
# Extract scores
rda_env_scores <- as.data.frame(scores(RAI_rda, display = "bp", scaling = 2))  # FSC predictors
rda_RAI_scores <- as.data.frame(scores(RAI_rda, display = "species", scaling = 2))  # Species



# Exploring unconstrained Variance ####
# **** RUN PCA FSC and TRAITS + RFU.R before ***** 
combined_species_traits_FSC_data <- cbind(FD_summary[, c("PC1_FSC", "PC2_FSC")], plot_RAI_values_scaled, 
                                          response_scaled[, c("AdultBodyMass_g", "HomeRange_km2", 
                                                              "TrophicLevel", "LRO", "Social", "DietBreadth",
                                                              "ActivityPattern")])
)
            
                                                                                                              
# response_scaled loaded via other script. 
# not the same as response_scaled_RAI. 

str(combined_species_traits_FSC_data)
print(combined_species_traits_FSC_data)


cor_matrix_species <- cor(combined_species_traits_FSC_data, method = "pearson")

#install.packages("pheatmap")
library(pheatmap)
library(ggcorrplot)
ggcorrplot(cor_matrix_species, 
           method = "circle", 
           type = "lower", 
           lab = TRUE, 
           title = "Correlation Heatmap: FSC and Species and Traits")


# test for significance
# Initialize a matrix for p-values with the same dimensions as the correlation matrix
p_value_matrix <- matrix(NA, nrow = nrow(cor_matrix_species), ncol = ncol(cor_matrix_species))
rownames(p_value_matrix) <- rownames(cor_matrix_species)
colnames(p_value_matrix) <- colnames(cor_matrix_species)

# Iterate over each pair of variables
for (i in 1:nrow(cor_matrix_species)) {
  for (j in 1:ncol(cor_matrix_species)) {
    if (i != j) { # Avoid self-correlation
      test <- cor.test(cor_matrix_species[, i], cor_matrix_species[, j], method = "pearson")
      p_value_matrix[i, j] <- test$p.value
    }
  }
}

# Apply multiple testing correction (e.g., Holm method)
p_value_matrix_corrected <- p.adjust(p_value_matrix, method = "holm")

# Reshape to a matrix
p_value_matrix_corrected <- matrix(p_value_matrix_corrected, nrow = nrow(p_value_matrix))
rownames(p_value_matrix_corrected) <- rownames(cor_matrix_species)
colnames(p_value_matrix_corrected) <- colnames(cor_matrix_species)

# Visualize significant correlations
library(pheatmap)
pheatmap(p_value_matrix_corrected, cluster_rows = FALSE, cluster_cols = FALSE, 
         main = "Heatmap of Adjusted P-values", display_numbers = TRUE)



# Get into the significant results
# Create a data frame to store significant results
significant_results <- data.frame(
  Variable1 = character(),
  Variable2 = character(),
  Correlation = numeric(),
  P_Value = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  stringsAsFactors = FALSE
)

# Iterate over each pair of variables
for (i in 1:(ncol(cor_matrix_species) - 1)) {
  for (j in (i + 1):ncol(cor_matrix_species)) {  # Avoid duplicates
    # Compute correlation and test
    test <- cor.test(cor_matrix_species[, i], cor_matrix_species[, j], method = "pearson")
    
    # If significant, add to the results
    if (test$p.value < 0.05) {  # Threshold for significance
      significant_results <- rbind(significant_results, data.frame(
        Variable1 = colnames(cor_matrix_species)[i],
        Variable2 = colnames(cor_matrix_species)[j],
        Correlation = test$estimate,
        P_Value = test$p.value,
        CI_Lower = test$conf.int[1],
        CI_Upper = test$conf.int[2]
      ))
    }
  }
}

# Adjust p-values for multiple testing
significant_results$Adjusted_P_Value <- p.adjust(significant_results$P_Value, method = "holm")

# Filter for relationships still significant after correction
significant_results <- significant_results[significant_results$Adjusted_P_Value < 0.05, ]

# View the significant relationships
library(knitr)
kable(significant_results, caption = "Significant Correlations with Statistics")

      

