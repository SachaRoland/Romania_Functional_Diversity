
# PCA 

#install.packages("mvabund")
#install.packages(AER)
#install.packages(glmmTMB)
#install.packages("flextable")
library(flextable)
library(officer)  # For fp_border
#install.packages("semPlot")
library(semPlot)
# Fit a structural equation model for visualization
library(lavaan)
# Visualize the correlation matrix as a heatmap
library(ggplot2)
library(reshape2)
library(car)
library(gridExtra)

pacman::p_load(
  tidyverse, 
  DHARMa,     #residual diagnostics for hierarchical (multilevel/mixed) regression models
  emmeans,    #created estimated marginal means
  DescTools,  #collection of miscellaneous basic statistic functions
  lme4,       #Fit linear and generalized linear mixed-effects models.
  lubridate,  #handles dates and times
  nlme,       #Linear and Nonlinear Mixed Effects Models
  stringi,
  ggcorrplot,
  readxl,
  arsenal,    #for comparing dataframes
  multcomp, #Multiple Comparisons Using R
  vegan,  #for running PCA
  AER, # for the test of dispersion (for GLM only!)
glmmTMB) # the run a gebative binomial GLMM in case of overdispersion. 

select <- dplyr::select #always use select() in dplyr

setwd("/Users/sacharoland/Desktop/WUR/THESIS/FSC and vertebrates- Romania/Data Analysis/Data")
combined_data <- read_csv("FSC_combined_data_2.csv")
combined_data <- combined_data %>% separate(PlotID, c("Plotname", "Subplot"), "-")
colnames(combined_data)

# Calculate mean for all variables grouped by Plot
data <- combined_data %>%
  group_by(Plotname) %>%
  summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))
# View the result
print(data)


# Use first level variables 
# Select only numeric columns for PCA
pca_data <- data[, c("stem_count_standing_trees_small", "n_standing_trees_big", 
                              "stem_count_lying_trees_medium", "stem_count_lying_trees_big", 
                              "total_volume_m3", "mean_decay_stage", "vol_late_decay",
                              "stem_count_living_trees", "total_basal_area_m2", "diameter_span",                  
                              "min_diameter", "max_diameter", "diameter_mean_basal_area_cm",
                              "trees_70cm", "n_forest_stages", "freq_CF",                     
                              "freq_ML", "freq_TC", "freq_TFC",                      
                              "n_tree_species", "num_trees_0_5_to_1_3m", "num_trees_over_1_3",
                     "total_regeneration", "tree_cover", "shrub_cover", "brushwood_cover",
                     "herb_cover"
                        )]

new_names <- c("standing_trees (7-20cm)", "standing_trees (>20cm)", "lying_trees (20-50cm)", 
               "lying_trees (>30cm)", "total_deadwood_volume", "mean_decay_stage", "volume_late_decay_stage",
               "living_tree_count", "basal_area", "diameter_span",                  
               "min_diameter", "max_diameter", "mean_diameter_ba", "trees (>70cm)",
               "forest_stages_count", "CF_frequency", "ML_frequency", 
               "TC_frequency", "TFC_frequency", "tree_species_count", 
               "regen (0.5-1.3cm)", "regen (>1.3cm)", "total_regeneration_count", 
               "tree_cover", "shrub_cover", "brushwood_cover", "herb_cover")

colnames(pca_data) <- new_names

# Standardize the data to ensure all variables contribute equally
# Scale the data: mean = 0, sd = 1
pca_data_scaled <- scale(pca_data)
colnames(pca_data_scaled)

# Check PCA Assumptions ####

## Linearity
# Compute a correlation matrix
# Correlations near ±1 indicate strong linear relationships (positive or negative).
# Correlations near 0 suggest weak or no linear relationships.
cor_matrix <- cor(pca_data_scaled, use = "pairwise.complete.obs")

# Convert correlation matrix to long format for ggplot2
cor_long <- melt(cor_matrix)

# Create heatmap; Maybe not the best way to see what's correlated. 
ggplot(cor_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  labs(title = "Correlation Heatmap", x = "Variables", y = "Variables") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(caret)
# Detect linear dependencies
linear_combos <- findLinearCombos(as.matrix(pca_data_scaled))

# Check which variables are redundant
if (!is.null(linear_combos$remove)) {
  print("Variables causing linear dependencies:")
  print(colnames(pca_data_scaled)[linear_combos$remove])
} else {
  print("No linear dependencies found.")
}


# max_diameter causes linear dependencies --> remove 
# Use select() to remove the column 
pca_data_scaled <- as.data.frame(pca_data_scaled) %>% 
  select(-"max_diameter")


# Preserve Plotname column before scaling
plot_names <- data$Plotname  # Extract plot names

# Perform PCA
pca_result <- rda(pca_data_scaled)

# Summary of PCA
summary(pca_result)

# SCREE PLOT ####
# retrieving the eigenvalues of the principal components
eigenvalues <- eigenvals(pca_result)

# computing a summary of the eigenvalues
summary(eigenvals(pca_result))
# We would need to retain the 10 first PC to keep 90% of the variance explained. 

#to better link the Inertia/eigenvalues to the properties of the FSC
# data, compute the sum of value of all eigenvalues
sum(eigenvals(pca_result))
sum(diag(var(pca_data_scaled)))

scree_plot <- screeplot(pca_result, bstick = TRUE, main = "Scree Plot with Broken Stick Model")
# Broken stick model
broken_stick <- bstick(length(eigenvalues))
# Consider the 2 first axes as there are PCs that have eigen values larger than the 
# broken stick model


# BIPLOT ####
# Extract site scores (plot positions)
site_scores <- scores(pca_result, display = "sites", scaling = 2)
site_scores <- as.data.frame(site_scores)
# Add correct Plotname column
site_scores$Plotname <- plot_names  # Use the preserved plot names

## Extract PC1 scores for plots #####
PC1_scores <-  site_scores$PC1
PC2_scores <- site_scores$PC2

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

library(ggforce)
library(ggrepel)

library(ggplot2)
library(ggforce)
library(ggrepel)


# # Combined PCA Biplot
# ggplot() +
#   # Plot site points (plots)
#   geom_point(data = site_scores, aes(x = PC1, y = PC2), color = "black", size = 1) +
#   geom_text(data = site_scores, aes(x = PC1, y = PC2, label = Plotname), 
#             vjust = -1, hjust = 0.5, size = 3, color = "black") +
#   # Draw arrows for variables
#   geom_segment(data = arrow_data,
#                aes(x = 0, y = 0, xend = PC1, yend = PC2, color = Contribution_PC1),
#                arrow = arrow(length = unit(0.3, "cm")), size = 1) +
#   # Add variable names at arrow tips (same color as arrows)
#  geom_text(data = arrow_data,
#  aes(x = PC1, y = PC2, label = Variable, color = Contribution_PC1),
#  hjust = 0.5, vjust = -0.5, size = 4) +
#   # Customize color scale for arrows and labels
#   scale_color_gradient(low = "blue", high = "red") +
#   # Add titles and axis labels
#   labs(title = "PCA Biplot with Plot and Variable Names",
#        x = "PC1 (35.27%)",
#        y = "PC2 (13.31%)",
#        color = "Contribution to PC1") +
#   # Apply minimal theme
#   theme_minimal() +
#   # Adjust theme for readability
#   theme(plot.title = element_text(hjust = 0.5, size = 14),
#         legend.position = "right") +
#   # Optional: Zoom into main cluster
#   coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2))

# Variable (species) loadings
loadings <- scores(pca_result, display = "species")
print(loadings)

#variables could use some renaming... 

# Site scores (sample positions)
site_scores <- scores(pca_result, display = "sites")
print(site_scores)



################################################################################
# Combine PCA Biplot

# Extract loadings for species (variables)
loadings <- as.data.frame(scores(pca_result, display = "species", choices = 1:2))
loadings$Variable <- rownames(loadings)

# Rename columns for clarity
colnames(loadings) <- c("PC1", "PC2", "Variable")

# Calculate contributions for PC1 and PC2
loadings$Contribution_PC1 <- (loadings$PC1^2) / sum(loadings$PC1^2) * 100
loadings$Contribution_PC2 <- (loadings$PC2^2) / sum(loadings$PC2^2) * 100
# Print the loadings with contributions
print(loadings)

# Sort the table by absolute PC1 values
loadings_table <- loadings[order(-abs(loadings$PC1)), ]

# Round to two decimal places for better presentation
loadings_table$PC1 <- round(loadings_table$PC1, 2)
loadings_table$PC2 <- round(loadings_table$PC2, 2)
# Rearrange columns: Variable first, followed by PC1 and PC2
loadings_table <- loadings_table[, c("Variable", "PC1", "PC2")]

# Make a flextable object
ft <- flextable(loadings_table)

# Format the table in APA style
ft <- ft %>%
  set_header_labels(Variable = "Variable", PC1 = "Loading on PC1", PC2 = "Loading on PC2") %>%
  autofit() %>%
  border_remove() %>%
  border_outer(part = "all", border = fp_border(color = "black", width = 1)) %>%
  hline_top(border = fp_border(color = "black", width = 1)) %>%
  hline_bottom(border = fp_border(color = "black", width = 1)) %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  fontsize(size = 11, part = "all")

# Save the table to Word
#save_as_docx(ft, path = "PCA_Loadings_Table_APA.docx")




# Extract site scores
site_scores <- as.data.frame(scores(pca_result, display = "sites", choices = 1:2))
# Add site labels (if needed)
site_scores$Plotname <- rownames(site_scores)

# Rename columns for clarity
colnames(site_scores) <- c("PC1", "PC2", "Plotname")
plot_names <- data$Plotname
site_scores$Plotname <- plot_names
print(site_scores)

# Arrow data for ggplot
arrow_data <- loadings[, c("Variable", "PC1", "PC2", "Contribution_PC1", "Contribution_PC2")]
# View arrow data
print(arrow_data)

# Combine contributions into a total contribution score
arrow_data$Total_Contribution <- arrow_data$Contribution_PC1 + arrow_data$Contribution_PC2

# Select the top 10 variables based on total contribution
top_variables <- arrow_data[order(-arrow_data$Total_Contribution), ][1:10, ]
# View the filtered arrow data
print(top_variables)

# Create a long format for the full bar graph (PC1 and PC2 contributions)
arrow_data_long <- data.frame(
  Variable = rep(arrow_data$Variable, 2),
  Contribution = c(arrow_data$Contribution_PC1, arrow_data$Contribution_PC2),
  PC = rep(c("PC1", "PC2"), each = nrow(arrow_data))
)


# Biplot with all Variables; place by top variables if want tope 10 
pca_biplot <- ggplot() +
  # Plot site points
  geom_point(data = site_scores, aes(x = PC1, y = PC2), color = "gray30", size = 2) +
  geom_text(data = site_scores, aes(x = PC1, y = PC2, label = Plotname), 
            vjust = -1, size = 3, color = "black") +
  # Draw scaled arrows for top variables
  geom_segment(data = arrow_data, 
               aes(x = 0, y = 0, xend = PC1 * 1.5, yend = PC2 * 1.5, color = Contribution_PC1),
               arrow = arrow(length = unit(0.3, "cm")), size = 1) +
  # Add variable names at arrow tips
  geom_text(data = arrow_data, 
            aes(x = PC1 * 1.8, y = PC2 * 1.8, label = Variable, color = Contribution_PC1),
            size = 4, vjust = -0.5) +
  scale_color_gradient(low = "blue", high = "red", name = "Contribution to PC1") +
  # Add dynamic axis labels with variance explained
  # Add titles and axis labels
  labs(title = "PCA Biplot with Plot and Variable Names",
       x = "PC1 (35.27%)",
       y = "PC2 (13.31%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14))
print(pca_biplot)

# Contribution Bar Plot for All Variables
contribution_plot <- ggplot(arrow_data_long, aes(x = reorder(Variable, -Contribution), y = Contribution, fill = PC)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = c("PC1" = "blue", "PC2" = "red"), name = "Principal Component") +
  labs(title = "Variable Contributions to PC1 and PC2",
       x = "Variables", y = "Contribution (%)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "right")
print(contribution_plot)

# Save the PCA Biplot (Top 10 Variables)
# ggsave("pca_biplot_A4.png", pca_biplot, width = 8.27, height = 5.83, dpi = 300)  # A4 landscape: 8.27 x 5.83 inches

# Save the Bar Plot (All Variables)
# ggsave("contribution_barplot_A4.png", contribution_plot, width = 8.27, height = 5.83, dpi = 300)


# CLUSTERING IN THE BIPLOT --> NOT NEEDED...; 

#  data <- data %>%
#    mutate(Stand_type = ifelse(str_detect(Plotname, "UpperMoria|Rivendell|Isengard|WolfTracks|SalamanderValley|UpperArcer|LowerMoria"),
#                              "old_growth", 
#                               "comparison"))
# # 
# # # Create a grouping variable (e.g., based on forest type)
#  groups <- data$Stand_type # Replace with your grouping column
# # # Add site points
# # points(site_scores[, 1], site_scores[, 2], col = as.factor(groups), pch = 1)
# 
# # Draw hulls around groups
# # Define a color palette
# ordihull(site_scores, groups = groups, draw = "lines", col = c("black", "red"), lwd = 1)
# # Add legend for groups
# legend("topright", legend = unique(groups), col = c("black", "red"), pch = 19)


# Conclusion: Deadwood total volume explains much of the variance among plots. 
# Key Finding: Deadwood volume and late decay wood are structuring the plots, 
# suggesting they are critical variables differentiating forest conditions.
# Ecological Insight: These variables are proxies for forest maturity and structural 
# complexity, which have profound implications for biodiversity, habitat quality, 
# and ecosystem function.


# Include environmental gradients --> SLOPE for example
# # Fit environmental vectors
# envfit_result <- envfit(pca_result, your_environmental_data, scaling = 2)
# 
# # Add vectors to the biplot
# plot(pca_result, scaling = 2, type = "n", main = "PCA with Environmental Gradients")
# points(scores(pca_result, display = "sites", scaling = 2), 
#        pch = 19, col = as.factor(groups))
# arrows(0, 0, envfit_result$vectors[, 1], envfit_result$vectors[, 2], 
#        length = 0.1, col = "red")
# 
# # Add labels for environmental variables
# text(envfit_result$vectors[, 1], envfit_result$vectors[, 2], 
#      labels = rownames(envfit_result$vectors), col = "red", pos = 4)

# Next step add environmental gradient such as slope!!!


# PCA scores as Explanatory variables####
# Extract site scores for PC1 and PC2
pca_scores <- scores(pca_result, display = "sites", scaling = 2)
pca_scores <- as.data.frame(pca_scores) # Convert to data frame
pca_scores$Plotname <- data$Plotname # Add plot identifiers (if applicable)

# Retain only PC1 and PC2
pca_scores <- pca_scores[, c("PC1", "PC2", "Plotname")]






# Import Ticks data. 
tick_data <- read_csv("Ticks_data.csv")

# Assuming functional diversity and tick abundance data has a 'Plot' column
merged_data <- merge(pca_scores, tick_data, by = "Plotname")


# DIRECT Effect - FSC --> Ticks #####
## GLMM - larvae #### 
# GLMM is used due to tcik data being count. --> poison distribution applies. 
larvae_glmm <- glmer(Larvae_total ~ PC1 +  (1|Plotname), 
                    family = poisson(link = "log"), data = merged_data)

# # Run a test of overdispersion manually. 
# # Calculate residual deviance and degrees of freedom
# residual_deviance <- sum(resid(larvae_glmm, type = "pearson")^2)
# df_residual <- nrow(merged_data) - length(fixef(larvae_glmm)) # Observations minus fixed effects
# 
# # Overdispersion ratio
# overdispersion <- residual_deviance / df_residual
# print(overdispersion)
# # Interpretation
# if (overdispersion > 1.5) {
#   print("Evidence of overdispersion. Consider a negative binomial model.")
# } else {
#   print("No significant overdispersion detected.")
# }

## GLMM - Nymphs ####
nymphs_glmm <- glmer(Nymphs_total ~ PC1 +  (1|Plotname), 
                     family = poisson(link = "log"), data = merged_data)
## GLMM - Adults ####
adults_glmm <- glmer(Adults_total ~ PC1 +  (1|Plotname), 
                     family = poisson(link = "log"), data = merged_data)

## GLMM - Total ticks  ####
ticks_total_glmm <- glmer(Ticks_total ~ PC1 +  (1|Plotname), 
                     family = poisson(link = "log"), data = merged_data)

 
summary(adults_glmm) # no significance
summary(nymphs_glmm) # no significance
summary(larvae_glmm) # no significance 
summary(ticks_total_glmm) # no significance 


# DIRECT Effect - Saturation Deficit --> Ticks #####
setwd("/Users/sacharoland/Desktop/WUR/THESIS/FSC and vertebrates- Romania/Data Analysis/Data")
SD <- read_csv("diurnal_variance_saturation_deficit.csv")

# Need to deal with "per day- per plot" for saturation deficit. 
# average day-night variance of saturation deficit.
# Calculate the mean saturation deficit per plot
mean_SD <- SD %>%
  group_by(Plot) %>%
  summarize(Mean_Variance_Saturation_Deficit = mean(Diurnal_Variance, na.rm = TRUE),
            SD_Variance_Saturation_Deficit = sd(Diurnal_Variance, na.rm = TRUE)
  )%>%
  rename(Plotname = Plot)

# Line plot with error bars (showing mean ± SD)
ggplot(mean_SD, aes(reorder(x = Plotname, -Mean_Variance_Saturation_Deficit), y = Mean_Variance_Saturation_Deficit)) +
  geom_point(size = 4, color = "blue") +  # Add points
  geom_errorbar(aes(
    ymin = Mean_Variance_Saturation_Deficit - SD_Variance_Saturation_Deficit,
    ymax = Mean_Variance_Saturation_Deficit + SD_Variance_Saturation_Deficit
  ), width = 0.2, color = "lightblue") +  # Add error bars
  labs(
    title = "Diurnal Variance of Saturation Deficit by Plot",
    x = "Plot",
    y = "Mean Diurnal Variance (kPa)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

SEM_merged_data <- merge(merged_data, mean_SD, by = "Plotname")


## GLMM - Nymphs ####
larvae_glmm_SD <- glmer(Larvae_total ~ Mean_Variance_Saturation_Deficit +  (1|Plotname), 
                        family = poisson(link = "log"), data = SEM_merged_data)

## GLMM - Nymphs ####
nymphs_glmm_SD <- glmer(Nymphs_total ~ Mean_Variance_Saturation_Deficit +  (1|Plotname), 
                     family = poisson(link = "log"), data = SEM_merged_data)

## GLMM - Adults ####
adults_glmm_SD <- glmer(Adults_total ~ Mean_Variance_Saturation_Deficit +  (1|Plotname), 
                        family = poisson(link = "log"), data = SEM_merged_data)

## GLMM - Total ticks  ####
ticks_total_glmm_SD <- glmer(Ticks_total ~ Mean_Variance_Saturation_Deficit +  (1|Plotname), 
                             family = poisson(link = "log"), data = SEM_merged_data)

summary(adults_glmm_SD) # significant
summary(nymphs_glmm_SD) # no significance
summary(larvae_glmm_SD) # no significance 
summary(ticks_total_glmm_SD) # no significance 


# DIRECT Effect - FSC --> Saturation Deficit #####
## GLMM ####

shapiro.test(SEM_merged_data$Mean_Variance_Saturation_Deficit)
# data is approximatly normal

SEM_merged_data$Plotname <- as.factor(SEM_merged_data$Plotname)

FSC_SD_lm <- lm(Mean_Variance_Saturation_Deficit ~ PC1, data = SEM_merged_data)

plot(residuals(FSC_SD_lm) ~ fitted(FSC_SD_lm))
abline(h = 0, col = "red")
# Perform the Breusch-Pagan test
bptest(FSC_SD_lm) # p-value > 0.05: Fail to reject H₀, indicating homoscedasticity.

summary(FSC_SD_lm)


# INDIRECT effect via Saturation Deficit ####
#install.packages("mediation")
library(mediation)

# Mediator Model: Saturation Deficit ~ FSC (PC1)
# Keep for all life stages models
mediator_model <- lm(Mean_Variance_Saturation_Deficit ~ PC1, data = SEM_merged_data)

# SEM ####
SEM_merged_data <- merge(merged_data, mean_SD, by = "Plotname")
## Larvae ####
# Outcome Model: Tick Abundance ~ Saturation Deficit + FSC (PC1)
larvae_outcome_model <- glm(Larvae_total ~ Mean_Variance_Saturation_Deficit + PC1, 
                     family = poisson(link = "log"), data = SEM_merged_data)

larvae_mediation_result <- mediate(mediator_model, larvae_outcome_model, 
                            treat = "PC1", mediator = "Mean_Variance_Saturation_Deficit",
                            boot = TRUE, sims = 5000) # Use bootstrapping for confidence intervals

## Nymphs ####
nymphs_outcome_model <- glm(Nymphs_total ~ Mean_Variance_Saturation_Deficit + PC1, 
                            family = poisson(link = "log"), data = SEM_merged_data)
# Run mediation analysis
nymphs_mediation_result <- mediate(mediator_model, nymphs_outcome_model, 
                                   treat = "PC1", mediator = "Mean_Variance_Saturation_Deficit",
                                   boot = TRUE, sims = 5000) # Use bootstrapping for confidence intervals

## Adults ####
adults_outcome_model <- glm(Adults_total ~ Mean_Variance_Saturation_Deficit + PC1, 
                            family = poisson(link = "log"), data = SEM_merged_data)
# Run mediation analysis
adults_mediation_result <- mediate(mediator_model, adults_outcome_model, 
                                   treat = "PC1", mediator = "Mean_Variance_Saturation_Deficit",
                                   boot = TRUE, sims = 5000) # Use bootstrapping for confidence intervals

## Ticks total ####
tick_total_outcome_model <- glm(Ticks_total ~ Mean_Variance_Saturation_Deficit + PC1, 
                            family = poisson(link = "log"), data = SEM_merged_data)
# Run mediation analysis
tick_total_mediation_result <- mediate(mediator_model, tick_total_outcome_model, 
                                   treat = "PC1", mediator = "Mean_Variance_Saturation_Deficit",
                                   boot = TRUE, sims = 5000) # Use bootstrapping for confidence intervals

# Summary of mediation results
summary(larvae_mediation_result)
summary(nymphs_mediation_result)
summary(adults_mediation_result)
summary(tick_total_mediation_result)

# Interpret the Results
# The output will include:
# ACME (Average Causal Mediation Effect): The indirect effect of FSC on tick abundance via saturation deficit.
# ADE (Average Direct Effect): The direct effect of FSC on tick abundance, independent of saturation deficit.
# Total Effect: The combined effect of FSC on tick abundance.
# Proportion Mediated: The proportion of the total effect explained by the indirect pathway.

# Nothing measured effect is statistically significant. 

# Visualize the Paths ####
# You can use the semPlot package for a graphical representation of the mediation model.

## Larvae ####
# Define the model
larvae_sem_model <- '
  Mean_Variance_Saturation_Deficit ~ PC1
 Larvae_total ~ Mean_Variance_Saturation_Deficit + PC1
'
# Fit the model
fit <- sem(larvae_sem_model, data = SEM_merged_data)
# Plot the model
semPaths(fit, what = "std", layout = "circle", edge.label.cex = 1.2)

## Nymphs ####
# Define the model
nymphs_sem_model <- '
  Mean_Variance_Saturation_Deficit ~ PC1
 Nymphs_total ~ Mean_Variance_Saturation_Deficit + PC1
'
# Fit the model
fit <- sem(nymphs_sem_model, data = SEM_merged_data)
# Plot the model
semPaths(fit, what = "std", layout = "circle", edge.label.cex = 1.2)

## Adults ####
# Define the model
adults_sem_model <- '
  Mean_Variance_Saturation_Deficit ~ PC1
 Adults_total ~ Mean_Variance_Saturation_Deficit + PC1
'
# Fit the model
fit <- sem(adults_sem_model, data = SEM_merged_data)
# Plot the model
semPaths(fit, what = "std", layout = "circle", edge.label.cex = 1.2)

## Adults ####
# Define the model
ticks_total_sem_model <- '
  Mean_Variance_Saturation_Deficit ~ PC1
 Ticks_total ~ Mean_Variance_Saturation_Deficit + PC1
'
# Fit the model
fit <- sem(ticks_total_sem_model, data = SEM_merged_data)
# Plot the model
semPaths(fit, what = "std", layout = "circle", edge.label.cex = 1.2)


# c) Multivariate Analysis (e.g., RDA):
# If you want to explore how PC1 and PC2 influence multiple response variables 
# (e.g., vertebrate traits and tick stages):







