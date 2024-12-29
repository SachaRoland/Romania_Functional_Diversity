

# Ticks BINOMIAL and Saturation Deficit ####
# RATIONALE : detection for larvae is definitely limited. Finding 1 automatcially means many 
# Nymphs can stay as counts 
# Adults are rare --> binomial.


# Setting the working directory
setwd("/Users/sacharoland/Desktop/WUR/THESIS/FSC and vertebrates- Romania/Data Analysis/Data")

# Reading a semicolon-delimited file
data <- read.table("tick_transects.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE)
colnames(data)
str(data)

data$Temperature <- as.numeric(data$Temperature)
data$LL_thickness <- as.numeric(data$LL_thickness)
data$Humidity <- as.numeric(data$Humidity)
data$Larvae <- as.numeric(data$Larvae)
data$Nymphs <- as.numeric(data$Nymphs)
data$Males <- as.numeric(data$Males)
data$Females <- as.numeric(data$Females)

# Replace "N/A" with NA in the entire data frame
data[data == "N/A"] <- NA

# View the data frame with SD calculated
print(head(data))
data <- data[!is.na(data$Temperature), ]

# Calculate Saturation Vapor Pressure (SVP) in kPa
data$SVP <- 0.61078 * exp((17.27 * data$Temperature) / (data$Temperature + 237.3))

# Calculate Saturation Deficit (SD) ####
data$Saturation_Deficit <- (1 - (data$Humidity / 100)) * data$SVP

# View the data frame with SD calculated
print(head(data))
data <- data[!is.na(data$Temperature), ]

# Create Adults column
data <- data %>%
  mutate(Adults = Males + Females)
view(data)

# Transform larvae and adults to binomial (1/0)
# Create a binomial larv
data$Larvae <- ifelse(data$Larvae > 0, 1, 0)
data$Adults <- ifelse(data$Adults > 0, 1, 0)
head(data)

write.csv(data, "Ticks_data.csv")
head(data)


# think of models adapted to binomial data. 
# If your data is grouped (e.g., multiple transects per plot), a mixed-effects model 
# accounts for random effects such as plot variability.
library(lme4)


## Larvae  ####
# -> no effect 
larvae_model_mixed <- glmer(Larvae ~ Saturation_Deficit + LL_thickness + (1 | Plotname), 
                     family = binomial, 
                     data = data)
summary(larvae_model_mixed)

## Adults ####
# --> no effect
adults_model_mixed <- glmer(Adults ~ Saturation_Deficit + LL_thickness + (1 | Plotname), 
                            family = binomial, 
                            data = data)
summary(adults_model_mixed)


## Nymphs ####
## model_nymphs can stay as count. 
# Does not include many zeros so maybe apply another model. 

# Check for overdispersion 
model_nymphs <- glm(Nymphs ~ Saturation_Deficit + LL_thickness, family = poisson, data = data)
summary(model_nymphs) #VARIANCE > DF 

#If the dispersion statistic > 1, the data are overdispersed, and a negative binomial model may be better....

dispersion <- sum(residuals(model_nymphs, type = "pearson")^2) / df.residual(model_nymphs)
cat("Dispersion statistic:", dispersion)  ### > 1 --> overdispersion 

# Compare models for Nymphs
#Binomial
library(MASS)
model_nbr_nymphs <- glm.nb(Nymphs ~ Saturation_Deficit + LL_thickness, data = data)
summary(model_nbr_nymphs)

# zero inflated 
library(pscl)
model_zinb_nymphs <- zeroinfl(Nymphs~ Saturation_Deficit + LL_thickness | Saturation_Deficit + LL_thickness, 
                              data = data, dist = "negbin")
summary(model_zinb_nymphs)

# Mixed-Effects Negative Binomial Regression
library(glmmTMB)
model_nymph <- glmmTMB(Nymphs ~ Saturation_Deficit + LL_thickness + (1 | Plotname), 
                       family = nbinom2, 
                       data = data)
summary(model_nymph)

AIC(model_nbr_nymphs, model_nymph)
# Mixed-Effects Negative Binomial Regression has best results
# no significant effect. 

#install.packages("effects")
library(effects)
plot(allEffects(model_nymph))

################################################################################
# Load necessary libraries
#install.packages("broom.mixed") 
library(broom.mixed)
library(dplyr)  # For data wrangling
library(glmmTMB)

# Extract summaries of each model ####
# Larvae model
larvae_summary <- tidy(larvae_model_mixed) %>%
  filter(term != "(Intercept)") %>%  # Exclude intercept for clarity
  mutate(Response = "Larvae", Model = "Binomial GLMM")

# Adults model
adults_summary <- tidy(adults_model_mixed) %>%
  filter(term != "(Intercept)") %>%
  mutate(Response = "Adults", Model = "Binomial GLMM")

# Nymphs: Negative Binomial Mixed Model
nymphs_summary <- tidy(model_nymph) %>%
  filter(term != "(Intercept)") %>%
  mutate(Response = "Nymphs", Model = "Negative Binomial GLMM")

# Combine all summaries into one table
model_summary <- bind_rows(larvae_summary, adults_summary, nymphs_summary) %>%
  select(Response, Model, term, estimate, std.error, p.value) %>%
  rename(
    Predictor = term,
    Estimate = estimate,
    StdError = std.error,
    PValue = p.value
  )

# Add Dispersion statistic for Poisson model (nymphs Poisson)
dispersion_stat <- data.frame(
  Response = "Nymphs",
  Model = "Poisson GLM",
  Predictor = "Dispersion",
  Estimate = dispersion,
  StdError = NA,
  PValue = NA
)

# Final combined summary with dispersion
final_summary <- bind_rows(model_summary, dispersion_stat)

# View the table
print(final_summary)


# Refined Results Table
# Refine and clean the summary table
library(dplyr)
library(ggplot2)

# Refine and clean the summary table
final_summary_clean <- final_summary %>%
  filter(!Predictor %in% c("sd__(Intercept)", "Dispersion")) %>%  # Remove random effects and dispersion rows
  mutate(Significance = case_when(
    PValue < 0.001 ~ "***",
    PValue < 0.01 ~ "**",
    PValue < 0.05 ~ "*",
    TRUE ~ ""
  )) %>%  # Add significance stars
  rename(
    `Response Variable` = Response,
    `Model Type` = Model,
    `Predictor Variable` = Predictor,
    `Coefficient Estimate` = Estimate,
    `Standard Error` = StdError,
    `P-Value` = PValue
  ) %>%
  select(`Response Variable`, `Model Type`, `Predictor Variable`,
         `Coefficient Estimate`, `Standard Error`, `P-Value`, `Significance`)

# View the refined table
print(final_summary_clean)

library(flextable)
library(officer)
# Create a flextable
ft <- flextable(final_summary_clean)

# Add formatting (optional)
ft <- ft %>%
  autofit() %>%  # Adjust column widths
  theme_vanilla()  # Apply a clean theme

# Save as a Word document
doc <- read_docx() %>%
  body_add_flextable(ft) %>%
  body_add_par("Table: Effect of Predictors on Tick Life Stages", style = "centered") %>%
  print(target = "Microclimate on ticks.docx")



# Create a coefficient plot
library(ggplot2)

# Adjusted coefficient plot
ggplot(final_summary_clean, aes(x = `Predictor Variable`, y = `Coefficient Estimate`, fill = `Response Variable`)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = `Coefficient Estimate` - 1.96 * `Standard Error`,
                    ymax = `Coefficient Estimate` + 1.96 * `Standard Error`),
                width = 0.15, position = position_dodge(width = 0.7)) +
  scale_fill_brewer(palette = "Set2") +  # Professional color palette
  theme_classic() +  # Cleaner theme
  labs(
    title = "Effect of Predictors on Tick Life Stages",
    y = "Coefficient Estimate (±95% CI)",
    x = "Predictor Variable",
    fill = "Life Stage"
  ) +
  coord_flip() +  # Flip coordinates for better readability
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

# No statistically significant effect of LL_thickness of Saturation deficit on the abundance 
# or absence probability of ticks at all life stages....



# MODELS ####
# Since FSC is measured at the plot level and I have multiple transects per plot, 
# I should aggregate transect-level data to the plot level and merge it with FSC scores.

# Aggregate data to plot level
plot_data <- data %>%
  group_by(Plotname) %>%
  summarise(
    Mean_LL = mean(LL_thickness, na.rm = TRUE),
    Mean_SD = mean(Saturation_Deficit, na.rm = TRUE),
    Larvae_Presence = mean(Larvae > 0),  # Proportion of transects with larvae
    Nymphs_Count = sum(Nymphs),          # Total nymph count per plot
    Adults_Presence = mean(Adults > 0),  # Proportion of transects with adults
    .groups = 'drop'
  )

# Add FSC scores to the aggregated data
plot_data <- plot_data %>%
  mutate(
    PC1_scores = PC1_scores,  # Assuming PC1_scores are in the same order as plots
    PC2_scores = PC2_scores
  )

colnames(plot_data)

# FSC -> Saturation deficit ###################################################
model_SD <- lm(Mean_SD ~ PC1_scores + PC2_scores, data = plot_data)
summary(model_SD)

# FSC -> LL_thickness ###################################################
model_LL <- lm(Mean_LL ~ PC1_scores + PC2_scores, data = plot_data)
summary(model_LL)


# Residual diagnostics
par(mfrow = c(2, 2))
plot(model_SD)  # For Saturation Deficit
plot(model_LL)  # For Leaf Litter 
# --> looking good. 

library(dplyr)
library(flextable)
library(officer)
# Microclimate models data
# Microclimate models data
library(dplyr)

# Microclimate models data
library(dplyr)

# Create the data frame without backticks
microclimate_data <- data.frame(
  Response_Variable = c("Saturation Deficit", "Saturation Deficit", 
                        "Leaf Litter Thickness", "Leaf Litter Thickness"),
  Model_Type = c("Linear Regression", "Linear Regression", 
                 "Linear Regression", "Linear Regression"),
  Predictor_Variable = c("PC1_scores", "PC2_scores", 
                         "PC1_scores", "PC2_scores"),
  Coefficient_Estimate = c(-0.093, -0.145, 0.103, -0.169),
  Standard_Error = c(0.074, 0.074, 0.246, 0.235),
  P_Value = c(0.221, 0.062, 0.681, 0.480)
)

# Add Significance column
microclimate_data <- microclimate_data %>%
  mutate(Significance = case_when(
    P_Value < 0.001 ~ "***",
    P_Value < 0.01 ~ "**",
    P_Value < 0.05 ~ "*",
    P_Value < 0.1 ~ ".",
    TRUE ~ ""
  ))

# View the refined table
print(microclimate_data)

# Create flextables
ft_microclimate <- flextable(microclimate_data) %>%
  autofit() %>%
  theme_vanilla() %>%
  set_caption(caption = "Effect of FSC on Microclimate Variables")




# FSC -> Tick abundance #######################################################

# for larvae and adults --> problem with that fact when now have proportion of transects where 
# presence was detected, which is not binomial anymore. 
# Could try to gte them back to binominal befining a thresholds.
# will need to ensure that there is enough varibility to detect potential effect. 

# Ensure response variables are binary (0 or 1)
plot_data <- plot_data %>%
  mutate(
    Larvae_Presence = ifelse(Larvae_Presence > 0, 1, 0),  # Convert proportions to binary
    Adults_Presence = ifelse(Adults_Presence > 0, 1, 0)   # Convert proportions to binary
  )

view(plot_data) # looks ok... 

## larvae ####
model_larvae <- glm(Larvae_Presence ~ PC1_scores + PC2_scores, family = binomial, data = plot_data)
summary(model_larvae)

#Residual Deviance
deviance(model_larvae) / df.residual(model_larvae) #close to 1 is good fit 

#Evaluate the model's ability to discriminate between presence (1) and absence (0)
# using a Receiver Operating Characteristic (ROC) curve:
library(pROC)
roc_curve <- roc(plot_data$Larvae_Presence, predict(model_larvae, type = "response"))
plot(roc_curve)
auc(roc_curve)

vif(model_larvae)


## adults ####
model_adults <- glm(Adults_Presence ~ PC1_scores + PC2_scores, family = binomial, data = plot_data)
summary(model_adults)
#Residual Deviance
deviance(model_adults) / df.residual(model_adults) #close to 1 is good fit 

#Evaluate the model's ability to discriminate between presence (1) and absence (0)
# using a Receiver Operating Characteristic (ROC) curve:
roc_curve <- roc(plot_data$Adults_Presence, predict(model_adults, type = "response"))
plot(roc_curve)
auc(roc_curve)

vif(model_adults)

## Assumptions met !



## nymphs ####
# FSC → Tick Abundance (Nymphs as count)
library(MASS)
model_nymphs <- glm.nb(Nymphs_Count ~ PC1_scores + PC2_scores, data = plot_data)

#Check for overdispersion before deciding on negative binomial vs. Poisson:
dispersion <- sum(residuals(glm(Nymphs_Count ~ PC1_scores + PC2_scores, family = poisson, data = plot_data),
                            type = "pearson")^2) / df.residual(glm(Nymphs_Count ~ PC1_scores + PC2_scores, 
                                                                   family = poisson, data = plot_data))
print(dispersion)
# If dispersion > 1, use the negative binomial model.


# Summarize results
summary(model_larvae)
summary(model_nymphs)
summary(model_adults)


# check models assumptions ####
#For each model, check residuals for normality and heteroscedasticity:
par(mfrow = c(2, 2))
plot(model_nymphs)  # Replace with your model of interest

library(lmtest)
bptest(model_nymphs) # not heteroskedasticity.

# Use AIC or adjusted R² to assess the model's explanatory power:
AIC(model_SD, model_LL, model_larvae, model_nymphs, model_adults)


# Tick models data
tick_data <- data.frame(
  Response_Variable = c("Larvae Presence", "Larvae Presence", 
                          "Adults Presence", "Adults Presence", 
                          "Nymphs Count", "Nymphs Count"),
  Model_Type = c("Binomial GLM", "Binomial GLM", 
                   "Binomial GLM", "Binomial GLM", 
                   "Negative Binomial GLMM", "Negative Binomial GLMM"),
  Predictor_Variable = c("PC1_scores", "PC2_scores", 
                           "PC1_scores", "PC2_scores", 
                           "PC1_scores", "PC2_scores"),
  Coefficient_Estimate = c(0.136, -0.567, 0.273, -0.373, -0.080, -0.131),
  Standard_Error = c(0.412, 0.487, 0.412, 0.477, 0.148, 0.148),
  P_Value = c(0.742, 0.244, 0.508, 0.434, 0.590, 0.376)
) %>%
  mutate(Significance = case_when(
    P_Value < 0.001 ~ "***",
    P_Value < 0.01 ~ "**",
    P_Value < 0.05 ~ "*",
    P_Value < 0.1 ~ ".",
    TRUE ~ ""
  ))

# Create flextables ####
ft_ticks <- flextable(tick_data) %>%
  autofit() %>%
  theme_vanilla() %>%
  set_caption(caption = "Effect of FSC on Ticks Variables")


# Save tables as Word document
# doc <- read_docx() %>%
#   body_add_flextable(ft_microclimate) %>%
#   body_add_par("", style = "Normal") %>%
#   body_add_flextable(ft_ticks) %>%
#   print(target = "FSC_Results_Tables.docx")


# Summarize results
summary(model_larvae)
summary(model_nymphs)
summary(model_adults)
summary(model_SD)
summary(model_LL)


# SEM #########################################################################
#install.packages("piecewiseSEM")
library(piecewiseSEM)

sem_model <- psem(
  lm(Mean_SD ~ PC1_scores + PC2_scores, data = plot_data),
  lm(Mean_LL ~ PC1_scores + PC2_scores, data = plot_data),
  glm(Larvae_Presence ~ PC1_scores + PC2_scores + Mean_SD, family = binomial, data = plot_data)
)

summary(sem_model)













