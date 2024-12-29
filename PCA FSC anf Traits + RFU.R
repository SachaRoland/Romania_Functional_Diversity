# PCAs
library(vegan)

# FSC variables and species
# with Rodent Functional Units (RFU)

# FSC variables and traits
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

species_only


# # Import abundances
# setwd("/Users/sacharoland/Desktop/WUR/THESIS/FSC and vertebrates- Romania/Data Analysis/Data")
# RAI <- read.csv("All_densities.csv")
# colnames(RAI)
# RAI_noR <- RAI %>% select(-"X", -"Apodemus_sp_",-"Muscardinus_avellanarius", -"Myodes_glareolus", -"Shrew",
#                           -"Glis_glis", -"Apodemus_agrarius")
# 
# # Adapt plotname columns and save names in seprate object.  
# rownames(RAI_noR) <- RAI_noR$Plotname
# RAI_noR <- RAI_noR[, -which(names(RAI_noR) == "Plotname")]
# 
# print(RAI_noR)
# 
# #Ensure the RAI_noR data frame is fully numeric
# str(RAI_noR)
# 
# RAI_noR <- RAI_noR %>%
#   mutate(
#     Martes = rowSums(across(c(Martes, Martes_martes, Martes_foina)), na.rm = TRUE)  # Combine values
#   ) %>%
#   dplyr::select(-Martes_foina, -Martes_martes)  # Drop the unused columns
# 
# # Standardize column names in species_only for matching
# colnames(RAI_noR) <- colnames(RAI_noR) %>%
#   gsub("_", " ", .)  # Replace underscores with spaces
# 
# 
# # calculate the weighted trait values for each plot
# # Find common species
# common_species <- intersect(colnames(RAI_noR), rownames(traits_matrix_noR))
# 
# # Subset both datasets to include only common species
# RAI_noR <- RAI_noR[, common_species, drop = FALSE]
# traits_matrix_noR <- traits_matrix_noR[common_species, , drop = FALSE]

#Multiply RAI_noR (plots × species) with traits_matrix_noR (species × traits)

# Convert RAI_noR and traits_matrix_noR to matrices
species_only_matrix <- as.matrix(species_only)
traits_matrix_updated_matrix <- as.matrix(traits_matrix_updated)

# Calculate weighted mean trait values for each plot
# which represents the weighted mean trait values for each plot based on species densities
plot_trait_values <- species_only_matrix %*% traits_matrix_updated_matrix  / rowSums(species_only_matrix)
# Check the structure of the output
str(plot_trait_values)

# View results
print(plot_trait_values)

FD_summary$PC1_FSC <- PC1_scores # From PCA for FSC --> so the axis explaining most of the variance of the variance 
FD_summary$PC2_FSC <- PC2_scores 


# values close to 1 or -1 suggest strong relationship. 
barplot(loadings_PC1, main = "Trait Contributions to PC1", las = 2, col = "blue")


# RUN FSC-PCA.R FIRST ####
# Include both traits and sites (plots) in the PCA biplot, combining forest structural 
# complexity (FSC) variables and traits.

# Load required libraries
library(vegan)
library(ggplot2)
library(dplyr)
library(grid)

# Load required libraries
library(vegan)
library(ggplot2)
library(dplyr)
library(grid)

# STEP 1: Prepare Data -------------------------------------------------------

# Combine FSC variables (pca_data_scaled) and traits (plot_trait_values)
fsc_traits_data <- cbind(pca_data_scaled, plot_trait_values)

# Ensure all data is numeric and scaled
fsc_traits_scaled <- scale(fsc_traits_data)

# STEP 2: Perform PCA --------------------------------------------------------

# Perform PCA on the combined dataset
pca_combined <- rda(fsc_traits_scaled)

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
  filter(!Variable %in% colnames(plot_trait_values)) %>%  # Exclude traits
  arrange(desc(Total_Contribution)) %>%
  slice_head(n = 10)

# Filter FSC variables and include all traits
arrow_data_filtered <- arrow_data %>%
  filter(Variable %in% top_fsc_variables$Variable | Variable %in% colnames(plot_trait_values))

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

# STEP 6: Differentiate FSC Variables and Traits -----------------------------

# Label variable types
arrow_data_filtered$Type <- ifelse(arrow_data_filtered$Variable %in% colnames(plot_trait_values), 
                                   "Trait", 
                                   "FSC Variable")

# STEP 7: Create Biplot with ggplot2 -----------------------------------------
# install.packages("ggforce")
# install.packages("ggrepel")
library(ggforce)
library(ggrepel)


plot_FSC_Traits <-ggplot() +
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
  scale_color_manual(values = c("FSC Variable" = "blue", "Trait" = "red", 
                                "old_growth" = "green", "comparison" = "purple")) +
  scale_fill_manual(values = c("old_growth" = "green", "comparison" = "purple")) +
  # Minimal theme with adjustments
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "right") +
  # Optional: Zoom into main cluster
  coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5))
print(plot_FSC_Traits)

#ggsave("PCA FSC-TRAITS_plot_A4.png", plot_FSC_Traits, width = 8.27, height = 11.69, units = "in", dpi = 300)



#-------------------------------------------
# Perform PCA on trait data to visualize trait clustering along PC1:###

# Variance of traits
plot_trait_values_scaled <- scale(plot_trait_values)
trait_variance <- apply(plot_trait_values_scaled, 2, var)
print(trait_variance)

cor_PC1_traits <- cor(FD_summary$PC1_FSC, plot_trait_values, method = "pearson")
print(cor_PC1_traits)

cor_PC2_traits <- cor(FD_summary$PC2_FSC, plot_trait_values, method = "pearson")
print(cor_PC2_traits)


combined_data <- cbind(plot_trait_values, FD_summary[, c("PC1_FSC", "PC2_FSC")])
combined_data_scaled <- scale(combined_data)
view(combined_data_scaled)

# Split into response and predictor variables
response_scaled <- combined_data_scaled[, colnames(plot_trait_values)]
predictors_scaled <- combined_data_scaled[, c("PC1_FSC", "PC2_FSC")]


# Threshold for filtering traits
# threshold <- 0.2
# # Identify informative traits
# informative_traits <- colnames(plot_trait_values)[apply(abs(cor_PC1_traits), 2, max) > threshold |
#                                                     apply(abs(cor_PC2_traits), 2, max) > threshold]
# 
# # Subset the trait data
# filtered_trait_response_scaled <- response_scaled [, informative_traits]


# Rerun RDA
trait_rda <- rda(response_scaled ~ PC1_FSC + PC2_FSC, data = as.data.frame(predictors_scaled))
colnames(combined_data_scaled)

# Extract the summary of the RDA object
rda_summary <- summary(trait_rda)


# Result Tbale ####
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
# STEP 4: Export the Tables to Word  ----------------------------------
# Save Partitioning Table
# save_as_docx(variance_table, path = "RDA_Partitioning_Table.docx")
# # Save Eigenvalues Table
# save_as_docx(eigenvalues_table, path = "RDA_Eigenvalues_Table.docx")


# Calculate the percentage of variance explained by each axis
# Extract constrained eigenvalues and proportions
# Extract constrained eigenvalues and proportions
str(rda_summary)
constrained_importance <- rda_summary$concont$importance

# Extract the percentage of variance explained by RDA1 and RDA2
variance_explained <- constrained_importance["Proportion Explained", ] * 100
names(variance_explained) <- c("RDA1", "RDA2")  # Name the axes for clarity

# Print results
print(variance_explained)


  # Extract species (traits) and site scores
      trait_scores <- scores(trait_rda, display = "species", scaling = 2)  # Traits (response)
      site_scores <- scores(trait_rda, display = "sites", scaling = 2)    # Sites (plots)
      predictor_scores <- scores(trait_rda, display = "bp", scaling = 2)  # Predictors (FSC)
      
      # Create base plot
      plot(trait_rda, scaling = 2, 
           main = "RDA Biplot: Traits and FSC Predictors", 
           display = c("sites", "species"),
           xlab = paste0("RDA1 (", round(variance_explained[1], 2), "%)"),
           ylab = paste0("RDA2 (", round(variance_explained[2], 2), "%)"),
           type = "n")  # Suppress default plotting
      
      # Add site points
      points(site_scores[, 1], site_scores[, 2], pch = 21, bg = "white", cex = 1.2)
      # Add plot names
      text(site_scores[, 1], site_scores[, 2], 
           labels = rownames(site_scores), col = "black", cex = 0.8, pos = 4)
      
      # Add trait arrows (response variables)
      arrows(0, 0, 
             trait_scores[, 1] * 1.5,  # Scale for visibility
             trait_scores[, 2] * 1.5, 
             col = "red", length = 0.1, lwd = 1.5)
      text(trait_scores[, 1] * 1.7, 
           trait_scores[, 2] * 1.7, 
           labels = rownames(trait_scores), col = "red", cex = 0.9)
      
      # Add predictor arrows (explanatory variables)
      arrows(0, 0, 
             predictor_scores[, 1] * 1.5, 
             predictor_scores[, 2] * 1.5, 
             col = "blue", length = 0.1, lwd = 1.5)
      text(predictor_scores[, 1] * 1.7, 
           predictor_scores[, 2] * 1.7, 
           labels = rownames(predictor_scores), col = "blue", cex = 0.9)
      
      # Add legend
      legend("topright", legend = c("Traits", "Predictors"), 
             col = c("red", "blue"), lty = 1, cex = 0.9)
      
      
  
# explore loading to RDA1 (70%) ####
# Extract loadings for traits and predictors
trait_loadings <- rda_summary$species[, 1:2]  # RDA1 and RDA2 for traits
predictor_loadings <- rda_summary$biplot[, 1:2]  # RDA1 and RDA2 for predictors

# Print the loadings
print("Trait Loadings:")
print(trait_loadings)

print("Predictor Loadings:")
print(predictor_loadings)

# PC1_FSC dominates RDA1, meaning that traits aligned with RDA1 are primarily driven 
# by overall forest structural complexity as represented by PC1.
# PC2_FSC dominates RDA2, suggesting it explains secondary patterns in trait variation.


#  Explore Unconstrained Variance (89.44%) ####

# Correlation of traits with unconstrained axes (PC1, PC2)
unconstrained_loadings <- rda_summary$species[, 3:4]  # PC1 and PC2 (unconstrained)
print("Unconstrained Trait Loadings:")
print(unconstrained_loadings)

# Correlation matrix of unconstrained axes with traits
unconstrained_correlations <- cor(plot_trait_values, rda_summary$sites[, 3:4])
print("Correlations with Unconstrained Axes:")
print(unconstrained_correlations)



# Enhanced RDA plot
library(ggplot2)
library(ggrepel)

# Convert loadings to data frames for ggplot
traits_df <- as.data.frame(trait_loadings)
traits_df$Trait <- rownames(trait_loadings)

predictors_df <- as.data.frame(predictor_loadings)
predictors_df$Predictor <- rownames(predictor_loadings)

sites_df <- as.data.frame(rda_summary$sites[, 1:2])
sites_df$Site <- rownames(rda_summary$sites)

# Plot better 
rda_plot <- ggplot() +
  # Add sites
  geom_point(data = sites_df, aes(x = RDA1, y = RDA2), color = "gray", size = 2, alpha = 0.7) +
  geom_text_repel(data = sites_df, aes(x = RDA1, y = RDA2, label = Site), size = 3) +
  # Add trait arrows
  geom_segment(data = traits_df, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.3, "cm")), color = "red", size = 1) +
  geom_text_repel(data = traits_df, aes(x = RDA1, y = RDA2, label = Trait), color = "red", size = 3) +
  # Add predictor arrows
  geom_segment(data = predictors_df, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.3, "cm")), color = "blue", size = 1) +
  geom_text_repel(data = predictors_df, aes(x = RDA1, y = RDA2, label = Predictor), color = "blue", size = 3) +
  # Customize plot
  labs(title = "Enhanced RDA Biplot", x = "RDA1", y = "RDA2") +
  theme_minimal()

print(rda_plot)


# Examine Trait Contributions to FD Indices ####

# Weighted trait contributions for each FD index
FD_results

# Prepare FD indices and CWM data
fd_indices <- as.data.frame(FD_results[c("FRic", "FEve", "FDiv")])
cwm_traits <- as.data.frame(FD_results$CWM)
# Check structure of the data frames
str(cwm_traits)
str(fd_indices)
# Convert non-numeric columns to numeric (if applicable)
cwm_traits <- cwm_traits %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate(across(where(is.factor), as.numeric))

fd_indices <- fd_indices %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate(across(where(is.factor), as.numeric))


# Ensure alignment of plots in both datasets
fd_indices <- fd_indices[rownames(cwm_traits), ]

# Compute correlations between CWM traits and FD indices
cor_results <- sapply(colnames(cwm_traits), function(trait) {
  sapply(colnames(fd_indices), function(fd_index) {
    cor(cwm_traits[[trait]], fd_indices[[fd_index]], method = "pearson")
  })
})


# Identify columns with zero variance
zero_variance_cols <- sapply(cwm_traits, function(col) sd(col, na.rm = TRUE) == 0)
print(zero_variance_cols)
# Remove columns with zero variance
cwm_traits <- cwm_traits[, !zero_variance_cols]


# Compute correlations again after removing zero-variance columns
cor_results <- sapply(colnames(cwm_traits), function(trait) {
  sapply(colnames(fd_indices), function(fd_index) {
    cor(cwm_traits[[trait]], fd_indices[[fd_index]], method = "pearson")
  })
})

# Convert the results to a readable data frame
cor_results_df <- as.data.frame(cor_results)
rownames(cor_results_df) <- colnames(fd_indices)

print(cor_results_df)


# Scatterplots of Traits vs. FD Indices ####
library(ggplot2)
# Example: Scatterplot for a specific trait and FD index
ggplot(data.frame(FRic = fd_indices$FRic, Trait = cwm_traits$AdultBodyMass_g), aes(x = Trait, y = FRic)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(x = "Adult Body Mass (g)", y = "FRic", title = "Adult Body Mass vs. FRic")

library(reshape2)

# Reshape for heatmap
## Add row names as a column for identification
cor_results_df$FD_Index <- rownames(cor_results_df)

# Melt the data frame into long format
cor_results_melted <- reshape2::melt(cor_results_df, id.vars = "FD_Index")

# Rename columns for clarity
colnames(cor_results_melted) <- c("FD_Index", "Trait", "Correlation")

# View the melted data
head(cor_results_melted)

cor_results_melted <- melt(cor_results_df)

# Heatmap ####
ggplot(cor_results_melted, aes(x = variable, y = FD_Index, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    breaks = seq(-1, 1, 0.2), limits = c(-1, 1)
  ) +
  labs(
    title = "Correlation Heatmap",
    x = "Traits",
    y = "Functional Diversity Indices",
    fill = "Correlation"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Test statistical significance ####
strong_correlations <- cor_results_melted %>%
  filter(abs(value) > 0.3)


# Initialize a data frame to store the results
cor_test_results <- data.frame()

# Loop through each combination of FD indices and traits
for (fd_index in unique(cor_results_melted$FD_Index)) {
  for (trait in unique(cor_results_melted$variable)) {
    # Extract the relevant columns for testing
    x <- cwm_traits[[trait]]
    y <- fd_indices[[fd_index]]
    
    # Perform the correlation test
    test <- cor.test(x, y, method = "pearson")
    
    # Store the results
    cor_test_results <- rbind(
      cor_test_results,
      data.frame(
        FD_Index = fd_index,
        Trait = trait,
        Correlation = test$estimate,
        P_Value = test$p.value,
        Conf_Low = test$conf.int[1],
        Conf_High = test$conf.int[2]
      )
    )
  }
}


# Filter significant correlations
significant_results <- cor_test_results %>% filter(P_Value < 0.05)

# View significant results
print(significant_results)
# No trait significantly drive FD indices along PC1_FSC



# Perhaps individual traits are not strongly related to FD indices, but combinations
# of traits are. Perform PCA or create composite indices for traits.

# install.packages("FactoMineR")
# install.packages("factoextra")
library(FactoMineR)
library(factoextra)

# Perform PCA on traits
pca_results <- PCA(cwm_traits, scale.unit = TRUE)

# Visualize PCA
fviz_pca_biplot(pca_results, repel = TRUE)

# Use PCA scores for correlation
pca_scores <- as.data.frame(pca_results$ind$coord)
cor_test_results_pca <- sapply(colnames(pca_scores), function(pc) {
  sapply(colnames(fd_indices), function(fd_index) {
    cor.test(pca_scores[[pc]], fd_indices[[fd_index]], method = "pearson")
  })
})

      
# Extract correlations, p-values, and confidence intervals
cor_values <- as.numeric(cor_test_results_pca[1, ])
p_values <- as.numeric(cor_test_results_pca[3, ])
conf_low <- as.numeric(cor_test_results_pca[10, ])
conf_high <- as.numeric(cor_test_results_pca[13, ])

# Create a tidy data frame
cor_summary <- data.frame(
  Dimension = rep(colnames(fd_indices), each = ncol(pca_scores)),
  PCA_Component = rep(paste0("Dim.", 1:ncol(pca_scores)), times = ncol(fd_indices)),
  Correlation = cor_values,
  P_Value = p_values,
  Conf_Low = conf_low,
  Conf_High = conf_high
)

# View results
print(cor_summary)

# Filter significant correlations
significant_results <- cor_summary %>%
  filter(P_Value < 0.05)

# View significant results
print(significant_results)
# none... 
      











      
      #Explore correlations between PC1 (FSC) and traits:
# Extract scores
rda_env_scores <- as.data.frame(scores(trait_rda, display = "bp", scaling = 2))  # FSC predictors
rda_trait_scores <- as.data.frame(scores(trait_rda, display = "species", scaling = 2))  # Traits

# Plot with ggplot2
# library(ggplot2)
# ggplot() +
#   geom_segment(data = rda_env_scores, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
#                arrow = arrow(length = unit(0.2, "cm")), color = "blue", size = 1) +
#   geom_text(data = rda_env_scores, aes(x = RDA1, y = RDA2, label = rownames(rda_env_scores)), 
#             color = "blue", vjust = -0.5) +
#   geom_segment(data = rda_trait_scores, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
#                arrow = arrow(length = unit(0.2, "cm")), color = "red", size = 1) +
#   geom_text(data = rda_trait_scores, aes(x = RDA1, y = RDA2, label = rownames(rda_trait_scores)), 
#             color = "red", vjust = -0.5) +
#   labs(title = "Trait Clustering with FSC", x = "RDA1", y = "RDA2") +
#   theme_minimal()





# Visualize Traits as Points #### 
# That results is already showing in the RDA plots. 
# trophic Level, Terrestriality and LRO are strongly correlated with FSC. 

# # Extract scores
# trait_scores <- as.data.frame(scores(trait_rda, display = "species", scaling = 2))  # Traits
# explanatory_scores <- as.data.frame(scores(trait_rda, display = "bp", scaling = 2))  # PC1_FSC and PC2_FSC
# 
# # Plot with ggplot2
# library(ggplot2)
# ggplot() +
#   # Trait points
#   geom_point(data = trait_scores, aes(x = RDA1, y = RDA2), color = "blue", size = 3) +
#   geom_text(data = trait_scores, aes(x = RDA1, y = RDA2, label = rownames(trait_scores)),
#             color = "blue", size = 3, vjust = -0.5) +
#   # FSC arrows
#   geom_segment(data = explanatory_scores, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
#                arrow = arrow(length = unit(0.3, "cm")), color = "red", size = 1) +
#   geom_text(data = explanatory_scores, aes(x = RDA1, y = RDA2, label = rownames(explanatory_scores)),
#             color = "red", size = 3, vjust = -0.5) +
#   # Labels and theme
#   labs(title = "Trait Clustering Along FSC",
#        x = "RDA1 (FSC Gradient 1)", y = "RDA2 (FSC Gradient 2)") +
#   theme_minimal()

# 
# # Correlation heatmap #
# # Combine PC1_FSC, PC2_FSC with informative traits
# correlation_data <- cbind(FD_summary[, c("PC1_FSC", "PC2_FSC")], combined_data_scaled)
# 
# # Calculate correlation matrix
# cor_matrix <- cor(correlation_data))
# 
# # Visualize as heatmap
# library(ggcorrplot)
# ggcorrplot(cor_matrix, 
#            method = "circle", 
#            type = "lower", 
#            lab = TRUE, 
#            title = "Correlation Heatmap: FSC and Traits")
# 
# # Test relationship between each trait and PC1_FSC ####
# # To confirm how specific traits relate to FSC, run linear models or correlations for each trait:
# # Use linear models for traits with normally distributed residuals.
# # Use Spearman’s rank correlation for monotonic relationships or non-normal residuals.
# 
# 
# # Connect observed clustering to functional diversity metrics like FDiv or FEve:
# #  Are traits related to FDiv (e.g., resource partitioning) clustering closer to FSC gradients?
# #  Do traits related to FEve (e.g., even distribution) spread uniformly?
# 
# 
# # Test whether traits associated with PC1_FSC (e.g., TrophicLevel, LRO) are 
# # significant predictors of FDiv.
# 
# # check distirbution of FDiv
# hist(FD_noR_summary$FDiv, breaks = 10, main = "Distribution of FDiv", xlab = "FDiv")
# shapiro.test(FD_noR_summary$FDiv)  # normal --> go for linear model
# 
# # Prep one df for the analysis. 
# traits_df <- as.data.frame(filtered_response_scaled)
# traits_df$Plot <- rownames(traits_df)  # Add plot names as a column
# all(rownames(FD_noR_summary) %in% traits_df$Plot)  # Should return TRUE
# 
# # Join df
# # Join the datasets
# traits_FDiv_combined_data <- FD_noR_summary %>%
#   left_join(traits_df, by = "Plot")
# 
# str(traits_FDiv_combined_data)
# print(traits_FDiv_combined_data)
# 
# model_FDiv <- lm(FDiv ~ PC1_FSC + TrophicLevel + LRO , data = traits_FDiv_combined_data)
# summary(model_FDiv)  # View model results
# 
# # Residual diagnostics
# plot(model_FDiv)
# # Test for heteroscedasticity
# bptest(model_FDiv)  # all good; no evidnec eof heteroscadicity. 
# 
# # Suggests that when accounting for traits (TrophicLevel, LRO, Terrestriality), 
# # PC1_FSC still has a marginally negative relationship with FDiv, but its 
# # significance is reduced because the additional predictors dilute its effect
# 
# # check for multicollinearity.
# # no multicollinearity 
# car::vif(model_FDiv)  # Variance inflation factor
# 
# 
# 
# # FDiv and traits 
# # FDiv and species
# 
