
# PCA with species Relative Abundance Index (RAI) and forest plots is a good idea 
# to explore community composition and identify clusters of plots with similar species
# Does NOT include rodents.

library(vegan)
library(dplyr)

species_RAI <- read.csv("All_densities.csv")
species_RAI <- species_RAI %>%
  dplyr::select(-"Apodemus_sp_", -"Muscardinus_avellanarius", "Myodes_glareolus",
                -"Shrew",  -"Glis_glis", -"Apodemus_agrarius", -"X")
colnames(species_RAI)

#Deal with Martes 
species_RAI <- species_RAI %>%
  mutate(
    Martes = rowSums(across(c(Martes, Martes_martes, Martes_foina)), na.rm = TRUE)  # Combine values
  ) %>%
  dplyr::select(-Martes_foina, -Martes_martes)  # Drop the unused columns

# Ensure only numeric columns (species RAI)
species_RAI_numeric <- species_RAI %>%
  dplyr::select(where(is.numeric))

# Check for missing values and handle them
species_RAI_numeric <- na.omit(species_RAI_numeric)

#If the species RAI data contains many zeros, use Hellinger transformation to reduce their impact:
species_RAI_hellinger <- decostand(species_RAI_numeric, method = "hellinger")

# Run PCA on transformed data
species_pca <- rda(species_RAI_hellinger, scale = TRUE)

# Summarize PCA
summary(species_pca)

# Visualize PCA with forest plot labels
ggplot(plot_scores, aes(x = PC1, y = PC2, label = PlotName)) +
  geom_point(size = 3, aes(color = forest_plot_attribute)) +
  geom_text(vjust = -0.5, size = 3) +
  labs(title = "PCA of Species Composition and Forest Plots",
       x = paste0("PC1 (", round(species_pca$sdev[1]^2 / sum(species_pca$sdev^2) * 100, 1), "% Variance)"),
       y = paste0("PC2 (", round(species_pca$sdev[2]^2 / sum(species_pca$sdev^2) * 100, 1), "% Variance)")) +
  theme_minimal()

# Check proportion of variance explained by each PC
eigenvalues <- eigenvals(species_pca)
variance_explained <- eigenvalues / sum(eigenvalues) * 100
print(variance_explained)

# SCREE PLOT ####
scree_plot <- screeplot(species_pca , bstick = TRUE, main = "Scree Plot with Broken Stick Model")
# Broken stick model
broken_stick <- bstick(length(eigenvalues))
# Consider the 2 first axes as there are PCs that have eigen values larger than the 
# broken stick model




# Exploring Connection to Community Turnover ####
# Calculate Bray-Curtis dissimilarity
community_dist <- vegdist(species_RAI, method = "bray")

# Test correlation with PC1 or FEve
mantel_test <- mantel(community_dist, vegdist(FEve, method = "euclidean"))
print(mantel_test)










# BIPLOT ####
# Extract site scores (plot positions)
site_scores <- scores(pca_result, display = "sites", scaling = 2)
site_scores <- as.data.frame(site_scores)
# Add correct Plotname column
site_scores$Plotname <- plot_names  # Use the preserved plot names

## Extract PC1 scores for plots #####
PC1_scores <-  site_scores$PC1
names(PC1_scores) <- site_scores$Plotname

# Extract species scores (arrows for variables)
variable_scores <- scores(pca_result, display = "species", scaling = 2)
arrow_data <- as.data.frame(variable_scores[, 1:2])  # Extract PC1 and PC2
colnames(arrow_data) <- c("PC1", "PC2")  # Rename columns
arrow_data$Variable <- rownames(variable_scores)  # Add variable names as a column

# Calculate contributions for PC1 (optional, not needed for plot)
arrow_data$Contribution_PC1 <- arrow_data$PC1^2 / sum(arrow_data$PC1^2)

# Keep only top 10 contributors
# top_contributors <- arrow_data %>%
#   arrange(desc(Contribution_PC1)) %>%
#   slice_head(n = 10)

# Combined PCA Biplot
ggplot() +
  # Plot site points (plots)
  geom_point(data = site_scores, aes(x = PC1, y = PC2), color = "black", size = 1) +
  geom_text(data = site_scores, aes(x = PC1, y = PC2, label = Plotname), 
            vjust = -1, hjust = 0.5, size = 3, color = "black") +
  # Draw arrows for variables
  geom_segment(data = arrow_data,
               aes(x = 0, y = 0, xend = PC1, yend = PC2, color = Contribution_PC1),
               arrow = arrow(length = unit(0.3, "cm")), size = 1) +
  # Add variable names at arrow tips (same color as arrows)
  geom_text(data = arrow_data,
            aes(x = PC1, y = PC2, label = Variable, color = Contribution_PC1),
            hjust = 0.5, vjust = -0.5, size = 4) +
  # Customize color scale for arrows and labels
  scale_color_gradient(low = "blue", high = "red") +
  # Add titles and axis labels
  labs(title = "PCA Biplot with Plot and Variable Names",
       x = "PC1 (35.27%)",
       y = "PC2 (13.31%)",
       color = "Contribution to PC1") +
  # Apply minimal theme
  theme_minimal() +
  # Adjust theme for readability
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "right")

# Variable (species) loadings
loadings <- scores(pca_result, display = "species")
print(loadings)

#variables could use some renaming... 

# Site scores (sample positions)
site_scores <- scores(pca_result, display = "sites")
print(site_scores)
  
  
  
