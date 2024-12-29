

# Functional traits matrix with rodent functional units (RFU)

# Funtional Diversity Indices by Villéger et al. 2010
# https://github.com/conservationscience/functionaltraits
# https://github.com/funecology/fundiversity

# https://github.com/ropensci/traits !!!!!
library(devtools)
library(remotes)
# remotes::install_github("ropensci/bold")
# remotes::install_github("ropensci/taxize")
# remotes::install_github("ropensci/traits")
library(traits)
library(dplyr)
library(ggplot2)


# Download PanTHERIA ####
# https://esapubs.org/archive/ecol/E090/184/metadata.htm
# read the tab-delimited text file
# Specify the path to the downloaded file
file_path <- "/Users/sacharoland/Desktop/WUR/THESIS/FSC and vertebrates- Romania/Data Analysis/Data/AGOUTI/PanTHERIA_1-0_WR93_Aug2008.txt" 

# Load the dataset
pantheria <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# View the first few rows
head(pantheria)

# View unique species names in the dataset
unique(pantheria$MSW93_Binomial)


# Target species 
species_list <- c(
  "Capreolus_capreolus", "Cervus_elaphus", "Sus_scrofa", 
  "Mustela_putorius", "Martes_martes", "Felis_silvestris", 
  "Vulpes_vulpes", "Canis_lupus", "Martes_foina", 
  "Lynx_lynx", "Ursus_arctos", "Meles_meles", 
  "Erinaceus_europaeus", "Lepus_europaeus", "Sciurus_vulgaris", 
  "Apodemus_flavicollis",  # yellow-nekced mouse as representative (45% of community)
  "Muscardinus_avellanarius", "Myodes_glareolus", 
  "Sorex_araneus",  # Common shrew as representative
  "Glis_glis", "Apodemus_agrarius"
)


# Replace underscores with spaces in species_list
species_list_cleaned <- gsub("_", " ", species_list)

# Convert to lowercase to ensure consistency
#species_list_cleaned <- tolower(species_list_cleaned)

# Check unmatched species again
unmatched_species <- species_list_cleaned[!species_list_cleaned %in% pantheria$MSW93_Binomial]
print(unmatched_species)

# Search for "myodes glareolus"
grep("Myodes", pantheria$MSW93_Binomial, value = TRUE)
# Search for "glis glis"
grep("glis", pantheria$MSW93_Binomial, value = TRUE)
species_list_cleaned[species_list_cleaned == "Glis glis"] <- "Myoxus glis"

# Check the dataset structure
str(pantheria)
# Filter for your species
mammal_traits <- pantheria[pantheria$MSW93_Binomial %in% species_list_cleaned, ]

# Check the filtered data
print(mammal_traits)
colnames(mammal_traits)

# Select relevant columns
selected_columns <- c(
  "MSW93_Order",                 # Order name
  "MSW93_Family",                # Family name
  "MSW93_Genus",                 # Genus name
  "MSW93_Species",               # Species name
  "MSW93_Binomial",              # Species binomial name
  "X5.1_AdultBodyMass_g",        # Body mass
  "X12.1_HabitatBreadth",        # Habitat breadth
  "X22.1_HomeRange_km2",         # Home range size
  "X1.1_ActivityCycle",          # Activity pattern
  "X6.1_DietBreadth",            # Diet breadth
  "X6.2_TrophicLevel",           # Trophic level
  "X10.2_SocialGrpSize",         # Social group size
  "X16.1_LittersPerYear",        # Litters per year
  "X15.1_LitterSize",            # Litter size
  "X17.1_MaxLongevity_m",        # Maximum longevity (in months)
  "X12.2_Terrestriality"
)


mammal_traits_selected <- mammal_traits[, selected_columns]
# View the selected traits
print(mammal_traits_selected)

# Replace -999 (missing value code in PanTHERIA) with NA
mammal_traits_selected[mammal_traits_selected == -999] <- NA


# # Map trophic level codes to descriptions
# trophic_level_map <- c("1" = "Herbivore", "2" = "Omnivore", "3" = "Carnivore")
# mammal_traits_selected$X12.1_TrophicLevel <- trophic_level_map[as.character(mammal_traits_selected$X12.1_TrophicLevel)]
# 
# # Map activity cycle codes to descriptions
# activity_cycle_map <- c("1" = "Nocturnal", "2" = "Diurnal", "3" = "Crepuscular") #"4" = "Cathemeral")
# mammal_traits_selected$X1.1_ActivityCycle <- activity_cycle_map[as.character(mammal_traits_selected$X1.1_ActivityCycle)]

# Rename columns
colnames(mammal_traits_selected) <- c(
  "Order", "Family", "Genus", "Species", "Binomial",
  "AdultBodyMass_g", "HabitatBreadth", "HomeRange_km2", "ActivityPattern", "DietBreadth",
  "TrophicLevel", "SocialGroupSize", "LittersPerYear", "LitterSize", "MaxLongevity_m", "Terrestriality")

#Fill in missing values for species + get rid of weird traits; 

#Remove HabitatBreadth... too arbitrary. 
mammal_traits_selected <- subset(mammal_traits_selected, select = - HabitatBreadth)


# HomeRange ####
mammal_traits_selected[, c("Binomial", "HomeRange_km2")]
# Fill missing HomeRange values based on suggested estimates
mammal_traits_selected$HomeRange_km2[mammal_traits_selected$Binomial == "Apodemus agrarius"] <- 0.3


# ActivityCycles ####

# Assign suggested activity pattern values
mammal_traits_selected$ActivityPattern <- c(
  "Apodemus agrarius" = 2,         # Nocturnal (Rodentia: Apodemus)
  "Apodemus flavicollis" = 1,      # Nocturnal (Rodentia: Apodemus)
  "Canis lupus" = 2,               # Cathemeral (Carnivora: Canidae)
  "Capreolus capreolus" = 2,       # Crepuscular (Artiodactyla: Cervidae)
  "Cervus elaphus" = 2,            # Crepuscular (Artiodactyla: Cervidae)
  "Erinaceus europaeus" = 1,       # Nocturnal (Eulipotyphla: Erinaceidae)
  "Felis silvestris" = 2,          # Nocturnal/crepuscular (Carnivora: Felidae)
  "Lepus europaeus" = 2,           # Crepuscular (Lagomorpha: Leporidae)
  "Lynx lynx" = 2,                 # Nocturnal/crepuscular (Carnivora: Felidae)
  "Martes foina" = 2,              # Nocturnal/crepuscular (Carnivora: Mustelidae)
  "Martes martes" = 1,             # Nocturnal (Carnivora: Mustelidae)
  "Meles meles" = 2,               # Cathemeral (Carnivora: Mustelidae)
  "Myoxus glis" = 1,               # Nocturnal (Rodentia: Gliridae) - Synonym for Glis glis
  "Muscardinus avellanarius" = 1,  # Nocturnal (Rodentia: Gliridae)
  "Mustela putorius" = 1,          # Nocturnal (Carnivora: Mustelidae)
  "Sorex araneus" = 1,             # Nocturnal (Eulipotyphla: Soricidae)
  "Sciurus vulgaris" = 2,          # Diurnal/crepuscular (Rodentia: Sciuridae)
  "Sus scrofa" = 2,                # Cathemeral (Artiodactyla: Suidae)
  "Ursus arctos" = 2,              # Cathemeral (Carnivora: Ursidae)
  "Vulpes vulpes" = 2              # Cathemeral (Carnivora: Canidae)
)


# DietBreadth ####
mammal_traits_selected$Binomial
# Add Diet Breadth values to mammal_traits_selected
diet_breadth_values <- c(
  "Apodemus agrarius" = 4,      #Feeds on seeds, fruits, green plant material, and insects
  "Apodemus flavicollis" = 4,   #Consumes seeds, fruits, invertebrates, and some plant material.
  "Canis lupus" = 5,            #Opportunistic omnivore: includes vertebrates, invertebrates, fruits, grasses, and roots.
  "Capreolus capreolus" = 3,    #Browses primarily on leaves, branches, and some seeds.
  "Cervus elaphus" = 3,         #Grazes on grasses, leaves, and seeds.
  "Erinaceus europaeus" = 5,    #Consumes invertebrates, fruits, roots, seeds, and some plant material.
  "Felis silvestris" = 3,       #Predominantly preys on vertebrates (small mammals and birds), occasionally consumes invertebrates.
  "Lepus europaeus" = 3,        #Feeds on grasses, leaves, and seeds
  "Lynx lynx" = 3,              #Primarily carnivorous: vertebrates (e.g., ungulates, small mammals), occasionally fruits or grasses.
  "Martes foina" = 5,           #Omnivorous: consumes vertebrates, invertebrates, fruits, seeds, and occasionally plant material.
  "Martes martes" = 5,          #Omnivorous: feeds on small mammals (vertebrates), invertebrates, fruits, seeds, and leaves.
  "Meles meles" = 6,            #Highly omnivorous: consumes vertebrates, invertebrates, fruits, seeds, roots, and grasses.
  "Muscardinus avellanarius" = 2,  #Primarily feeds on seeds and fruits.
  "Mustela putorius" = 3,          #Diet includes small vertebrates (mammals and amphibians), invertebrates, and occasional plant material.
  "Myoxus glis" = 4,               #Omnivorous: feeds on seeds, fruits, invertebrates, and leaves.
  "Sciurus vulgaris" = 3,          #Primarily seeds and nuts, occasionally fruits.
  "Sorex araneus" = 2,             #Insectivorous: consumes invertebrates and occasionally vertebrates.
  "Sus scrofa" = 7,                #Highly omnivorous: consumes vertebrates, invertebrates, fruits, seeds, roots/tubers, grasses, and plant material.
  "Ursus arctos" = 7,              #Highly omnivorous: consumes vertebrates, invertebrates, fruits, seeds, grasses, roots, and plant material.
  "Vulpes vulpes" = 6              #Generalist: feeds on vertebrates, invertebrates, fruits, seeds, grasses, and plant material.
)

# Map values to the Binomial column
mammal_traits_selected$DietBreadth <- diet_breadth_values[mammal_traits_selected$Binomial]
# View the updated data frame
print(mammal_traits_selected)


# TrophicLevel ####
mammal_traits_selected[, c("Binomial", "TrophicLevel")]
# 1 = Herbivore: Consumes only plant material
# 2 = Omnivore: Includes both animal (vertebrates/invertebrates) and plant material
# 3 = Carnivore: Feeds exclusively on vertebrates and/or invertebrates

# Create a named vector for trophic levels
trophic_level_values <- c(
  "Apodemus agrarius" = 2,        # Omnivore: Eats seeds, fruits, and invertebrates
  "Apodemus flavicollis" = 2,     # Omnivore: Similar to A. agrarius
  "Canis lupus" = 3,              # Carnivore: Primarily preys on vertebrates
  "Capreolus capreolus" = 1,      # Herbivore: Feeds on leaves, shoots, and shrubs
  "Cervus elaphus" = 1,           # Herbivore: Grazes on grasses and shrubs
  "Erinaceus europaeus" = 2,      # Omnivore: Eats invertebrates, fruits, and plant material
  "Felis silvestris" = 3,         # Carnivore: Hunts small vertebrates and invertebrates
  "Lepus europaeus" = 1,          # Herbivore: Feeds on grasses and herbs
  "Lynx lynx" = 3,                # Carnivore: Specializes in vertebrates
  "Martes foina" = 2,             # Omnivore: Mix of vertebrates, invertebrates, fruits
  "Martes martes" = 2,            # Omnivore: Similar diet to M. foina
  "Meles meles" = 2,              # Omnivore: Highly varied diet
  "Muscardinus avellanarius" = 1, # Herbivore: Feeds on nuts and fruits
  "Mustela putorius" = 3,         # Carnivore: Preys on small vertebrates and invertebrates
  "Myoxus glis" = 2,              # Omnivore: Combines plant materials and invertebrates
  "Sciurus vulgaris" = 1,         # Herbivore: Feeds on seeds, nuts, and plant material
  "Sorex araneus" = 3,            # Carnivore: Insectivorous; consumes invertebrates
  "Sus scrofa" = 2,               # Omnivore: Eats plants, invertebrates, and vertebrates
  "Ursus arctos" = 2,             # Omnivore: Highly varied diet
  "Vulpes vulpes" = 2             # Omnivore: Mix of vertebrates, invertebrates, and plants
)

mammal_traits_selected$TrophicLevel <- trophic_level_values[mammal_traits_selected$Binomial]


# SocialGroupSize ####
# Define binary social behavior: 1 = Social, 0 = Solitary
# Including notes for each species to justify the classification
social_behavior <- c(
  "Apodemus agrarius" = 0,        # Solitary: Overlapping nests, but no true social structure
  "Apodemus flavicollis" = 0,     # Solitary: Limited social interaction
  "Canis lupus" = 1,              # Social: Lives in packs with complex social structures
  "Capreolus capreolus" = 0,      # Solitary: Loosely gregarious only in winter
  "Cervus elaphus" = 1,           # Social: Forms harems during rutting and herds otherwise
  "Erinaceus europaeus" = 0,      # Solitary: Highly territorial
  "Felis silvestris" = 0,         # Solitary: Territorial, except for maternal care
  "Lepus europaeus" = 0,          # Solitary: Interaction limited to breeding season
  "Lynx lynx" = 0,                # Solitary: Territorial predator
  "Martes foina" = 0,             # Solitary: Minimal social interaction
  "Martes martes" = 0,            # Solitary: Similar to M. foina
  "Meles meles" = 1,              # Social: Lives in clans with shared burrow systems
  "Muscardinus avellanarius" = 0, # Solitary: Minimal social interaction, except during breeding
  "Mustela putorius" = 0,         # Solitary: Territorial predator
  "Myoxus glis" = 0,              # Solitary: Occasionally shares winter dens
  "Sciurus vulgaris" = 0,         # Solitary: Territorial; nests may overlap, but no social structure
  "Sorex araneus" = 0,            # Solitary: Highly territorial insectivore
  "Sus scrofa" = 1,               # Social: Lives in sounders (female-led groups)
  "Ursus arctos" = 0,             # Solitary: Interaction limited to mating and maternal care
  "Vulpes vulpes" = 0             # Solitary: Occasionally forms loose family groups during breeding
)

# Map binary social behavior to the Binomial column in mammal_traits_selected
# Assuming "Binomial" column contains the species names
mammal_traits_selected$Social <- social_behavior[mammal_traits_selected$Binomial]

# Remove the SocialGroupSize column (if it is no longer relevant)
mammal_traits_selected <- subset(mammal_traits_selected, select = -SocialGroupSize)

# View the updated data frame to confirm the changes
print(mammal_traits_selected)

# Check for any species that were not assigned a Social value
unassigned <- mammal_traits_selected$Binomial[is.na(mammal_traits_selected$Social)]
if (length(unassigned) > 0) {
  message("The following species were not assigned Social values:")
  print(unassigned)
} else {
  message("All species successfully assigned Social values!")
}



# LittersPerYear ####
mammal_traits_selected[, c("Binomial", "LittersPerYear")]
# Fill missing LittersPerYear values based on suggested estimates
mammal_traits_selected$LittersPerYear[mammal_traits_selected$Binomial == "Capreolus capreolus"] <- 1.00
mammal_traits_selected$LittersPerYear[mammal_traits_selected$Binomial == "Cervus elaphus"] <- 1.00
mammal_traits_selected$LittersPerYear[mammal_traits_selected$Binomial == "Martes foina"] <- 1.00
mammal_traits_selected$LittersPerYear[mammal_traits_selected$Binomial == "Meles meles"] <- 0.80
mammal_traits_selected$LittersPerYear[mammal_traits_selected$Binomial == "Ursus arctos"] <- 0.33
mammal_traits_selected$LittersPerYear[mammal_traits_selected$Binomial == "Vulpes vulpes"] <- 1.00

# LitterSize ####
# All good 


# Lifetime Reproductive output (LRO) ####
# Find average Longevity(year) instead of Max longevity (month)
# Create a data frame with species and average longevity
# https://www.demogr.mpg.de/longevityrecords/0203.htm

AverageLongevity <- data.frame(
  Binomial = c(
    "Apodemus agrarius", "Apodemus flavicollis", "Canis lupus", 
    "Capreolus capreolus", "Cervus elaphus", "Erinaceus europaeus", 
    "Felis silvestris", "Lepus europaeus", "Lynx lynx", 
    "Martes foina", "Martes martes", "Meles meles", 
    "Muscardinus avellanarius", "Mustela putorius", "Myoxus glis", 
    "Sciurus vulgaris", "Sorex araneus", "Sus scrofa", 
    "Ursus arctos", "Vulpes vulpes"
  ),
  AverageLongevity_y = c(
    2, 2, 8,    # Apodemus species and Canis lupus
    10, 15, 3,  # Capreolus capreolus, Cervus elaphus, Erinaceus europaeus
    7, 4, 10,   # Felis silvestris, Lepus europaeus, Lynx lynx
    8, 8, 14,   # Martes foina, Martes martes, Meles meles
    3, 5, 6,    # Muscardinus avellanarius, Mustela putorius, Myoxus glis
    6, 1.5, 10, # Sciurus vulgaris, Sorex araneus, Sus scrofa
    20, 5       # Ursus arctos, Vulpes vulpes
  ) # in years
)

# View the updated data frame
print(AverageLongevity)

# Merge average longevity into the existing dataset
mammal_traits_selected <- merge(mammal_traits_selected, AverageLongevity , by = "Binomial", all.x = TRUE)


# Calculate LRO
mammal_traits_selected$LRO <- mammal_traits_selected$LittersPerYear *
  mammal_traits_selected$LitterSize *
  mammal_traits_selected$AverageLongevity_y

mammal_traits_selected[, c("Binomial", "LRO")]


# Remove unnecessary columns
colnames(mammal_traits_selected)
# Remove the columns using dplyr::select()
mammal_traits_selected <- mammal_traits_selected |> 
  dplyr::select(-any_of(c("MaxLongevity_m")))

# View the cleaned data
print(mammal_traits_selected)

# Terrestriality ####
mammal_traits_selected$Terrestriality[mammal_traits_selected$Binomial == "Apodemus agrarius"] <- 1
mammal_traits_selected$Terrestriality[mammal_traits_selected$Binomial == "Apodemus flavicollis"] <- 1
mammal_traits_selected$Terrestriality[mammal_traits_selected$Binomial == "Lepus europaeus"] <- 1
mammal_traits_selected$Terrestriality[mammal_traits_selected$Binomial == "Sus scrofa"] <- 1


# Add Myodes glareolus ####
# Myodes glareolus is not in the PanTHERA data set.
# Create a new row for Myodes glareolus
new_species <- data.frame(
  Order = "Rodentia",
  Family = "Cricetidae",
  Genus = "Myodes",
  Species = "glareolus",
  Binomial = "Myodes glareolus",
  AdultBodyMass_g = 18.7,        # Grams
  HomeRange_km2 = 0.0015,        # Square kilometers
  ActivityPattern = 1,           # Nocturnal
  DietBreadth = 4,               # Generalist feeder
  TrophicLevel = 2,              # Omnivorous
  Social = 1,                    # Solitary
  LittersPerYear = 4,            # 3-5 litters/year
  LitterSize = 5.6,              # Average litter size
  AverageLongevity_y = 1.5,      # Average longevity in years
  Terrestriality = 1,            # Fully terrestrial
  LRO = 33.6
)

# Append the new row to the existing data frame
mammal_traits_selected <- rbind(mammal_traits_selected, new_species)
colnames(mammal_traits_selected)
colnames(new_species)
# Sort the data frame alphabetically by the column 'Binomial'
mammal_traits_selected <- mammal_traits_selected[order(mammal_traits_selected$Binomial), ]

# Calculate LRO
# 4 * 5.6 * 1.5



# RELATIVE ABUNDANCE INDEX ####
# Example RAI data (replace with your actual data)
RAI <- read.csv("All_densities.csv")
colnames(RAI)

mammal_traits_selected$Binomial

# NEED TO MATCH NAMES 
library(dplyr)

# Extract species columns from RAI, excluding non-species columns like "X" and "Plotname"
species_columns <- colnames(RAI)[!colnames(RAI) %in% c("X", "Plotname")]

# Standardize the names for matching
species_columns_clean <- gsub("_", " ", species_columns)    # Replace underscores with spaces
#species_columns_clean <- tolower(species_columns_clean)     # Convert to lowercase
#binomial_names_clean <- tolower(mammal_traits_selected$Binomial)

# Create a mapping between RAI species columns and Binomial names
species_mapping <- setNames(mammal_traits_selected$Binomial[match(species_columns_clean, mammal_traits_selected$Binomial)], species_columns)

# Check for unmatched species
unmatched_species <- names(species_mapping)[is.na(species_mapping)]
if (length(unmatched_species) > 0) {
  message("The following species in RAI did not match any Binomial names:")
  print(unmatched_species)
}

# The following species in RAI did not match any Binomial names:
#   [1] "Martes"       "Apodemus_sp_" "Shrew"        "Glis_glis" 


# Martes problem ####
# 2 species, one family. 
# Combine RAI values for Martes foina and Martes in the RAI data frame
# Combine the RAI values for Martes, Martes_martes, and Martes_foina
RAI <- RAI %>%
  mutate(
    Martes = rowSums(across(c(Martes, Martes_martes, Martes_foina)), na.rm = TRUE)  # Combine values
  ) %>%
  dplyr::select(-Martes_foina, -Martes_martes)  # Drop the unused columns

# Retain only one row for Martes in mammal_traits_selected
mammal_traits_selected <- mammal_traits_selected %>%
  filter(Binomial != "Martes martes") %>%  # Remove Martes martes
  mutate(
    Binomial = ifelse(Binomial == "Martes foina", "Martes", Binomial)  # Rename Martes foina to Martes
  )

# Only Martes in both df. 


# Deal with "Apodemus_sp_", "Shrew", "Glis_glis"
# Rename unmatched species in RAI
colnames(RAI) <- colnames(RAI) %>%
  gsub("Apodemus_sp_", "Apodemus flavicollis", .) %>%  # Map to a specific Apodemus species or leave as-is
  gsub("Shrew", "Sorex araneus", .) %>%            # Map to a representative shrew species
  gsub("Glis_glis", "Myoxus glis", .)

# # Standardize names for matching
species_columns_clean <- gsub("_", " ", colnames(RAI))
# # species_columns_clean <- tolower(species_columns_clean)
# # binomial_names_clean <- tolower(mammal_traits_selected$Binomial)

# Check for unmatched species
unmatched_species <- species_columns_clean[!species_columns_clean %in% mammal_traits_selected$Binomial]
if (length(unmatched_species) > 0) {
  message("The following species in RAI still do not match any Binomial names:")
  print(unmatched_species)
} else {
  message("All species in RAI now match Binomial names!")
}

print(unmatched_species) ## Need to deal ith non-species columns here... 


# Exclude non-species columns from RAI data
species_only <- RAI %>%
  dplyr::select(-c(X, Plotname))  # Drop the non-species columns


# Put species names same order across df
# Filter and reorder species columns in species_only to match Binomial names in mammal_traits_selected
species_only <- as.data.frame(species_only)
# Standardize column names in species_only for matching
colnames(species_only) <- colnames(species_only) %>%
  gsub("_", " ", .)  # Replace underscores with spaces

mismatched_species <- setdiff(mammal_traits_selected$Binomial, colnames(species_only))
# Print mismatched species
if (length(mismatched_species) > 0) {
  message("The following species are in mammal_traits_selected$Binomial but not in species_only:")
  print(mismatched_species)
} else {
  message("All species match between species_only and mammal_traits_selected$Binomial!")
}



# Align Column Names ####

# 1) standardize Binomial names for comparison ####
binomial_clean <- mammal_traits_selected$Binomial
view(mammal_traits_selected)
# Check again for mismatches
mismatched_species <- setdiff(binomial_clean, colnames(species_only))

if (length(mismatched_species) > 0) {
  message("The following species are still mismatched:")
  print(mismatched_species)
} else {
  message("Species are now aligned!")
}


# 2) Match rows and columns ####
# Set Plotname as row names
rownames(species_only) <- RAI$Plotname

# Drop the Plotname column if it exists
species_only <- species_only[, colnames(species_only) != "Plotname"]

# Reorder species_only columns to match Binomial order in mammal_traits_selected
species_only <- species_only[, binomial_clean]
view(species_only)

# Prepare Trait MATRIX 
# Set Binomial as row names
traits_matrix <- mammal_traits_selected %>%
  column_to_rownames(var = "Binomial") %>%
  dplyr::select(-Order, -Family, -Genus, -Species)  # Keep only numeric trait columns
view(traits_matrix)
print(traits_matrix)



# Functional Unit for Rodents #####
# Rodent species to group
group_small_terrestrial <- c("Apodemus agrarius", "Apodemus flavicollis", "Myodes glareolus")
group_arboreal <- c("Muscardinus avellanarius", "Myoxus glis")


# Aggregate RAIs for functional units in species_only
# Aggregate RAIs for functional units by summing across columns (plots)
small_terrestrial_RAI <- rowSums(species_only[, group_small_terrestrial], na.rm = TRUE)
arboreal_RAI <- rowSums(species_only[, group_arboreal], na.rm = TRUE)

# Remove original grouped rodent species
species_only <- species_only[, !(colnames(species_only) %in% c(
  group_small_terrestrial, group_arboreal
))]

# Add new columns for functional units
species_only$`Small Terrestrial Rodents` <- small_terrestrial_RAI
species_only$`Arboreal Rodents` <- arboreal_RAI


# View the updated abundance matrix # hell yeah
print(species_only) 
view(species_only)


# Calculate means for each group in traits matrix 
group_means <- function(species, traits_matrix) {
  traits_matrix %>%
    filter(rownames(traits_matrix) %in% species) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
}

# Functional units
small_terrestrial_traits <- group_means(group_small_terrestrial, traits_matrix)
arboreal_traits <- group_means(group_arboreal, traits_matrix)


# Assign meaningful row names for functional units
rownames(small_terrestrial_traits) <- "Small Terrestrial Rodents"
rownames(arboreal_traits) <- "Arboreal Rodents"


# Remove original grouped species from the trait matrix
traits_matrix_updated <- traits_matrix %>%
  filter(!rownames(traits_matrix) %in% c(group_small_terrestrial, group_arboreal))

# Combine back into the matrix with functional units as rows
traits_matrix_updated <- rbind(
  traits_matrix_updated,
  small_terrestrial_traits,
  arboreal_traits
)

# View the updated matrix
print(traits_matrix_updated)

# review some of the categorical traits forr the RFU
traits_matrix_updated["Small Terrestrial Rodents", "ActivityPattern"] <- 1
traits_matrix_updated["Small Terrestrial Rodents", "Social"] <- 0
traits_matrix_updated["Arboreal Rodents", "TrophicLevel"] <- 2



# Load necessary library
library(FD)

# Updated Trait Matrix: Remove unwanted columns
traits_matrix_updated <- traits_matrix_updated %>%
  dplyr::select(-LittersPerYear, -LitterSize, -AverageLongevity_y)

# View the updated matrix
#View(traits_matrix_updated)

# Identify continuous variables
continuous_vars <- c("AdultBodyMass_g", "HomeRange_km2", "DietBreadth", "LRO")

# Standardize continuous variables (mean = 0, sd = 1)
traits_matrix_updated[continuous_vars] <- scale(traits_matrix_updated[continuous_vars])

# Ensure categorical variables are treated correctly
traits_matrix_updated$ActivityPattern <- as.factor(traits_matrix_updated$ActivityPattern)
traits_matrix_updated$TrophicLevel <- as.factor(traits_matrix_updated$TrophicLevel)
traits_matrix_updated$Social <- as.numeric(as.character(traits_matrix_updated$Social))
traits_matrix_updated$Terrestriality <- as.numeric(as.character(traits_matrix_updated$Terrestriality))

# View structure of the matrix
str(traits_matrix_updated)
str(species_only)

# Check for mismatches
setdiff(rownames(traits_matrix_updated), colnames(species_only))  # Traits not in abundance
setdiff(colnames(species_only), rownames(traits_matrix_updated))  # Abundance not in traits

# Align names if needed
species_only <- species_only[, colnames(species_only) %in% rownames(traits_matrix_updated)]
traits_matrix_updated <- traits_matrix_updated[rownames(traits_matrix_updated) %in% colnames(species_only), ]



# Functional Diversity Analysis ####
FD_results <- dbFD(
  x = traits_matrix_updated, 
  a = species_only)

# View results
print(FD_results)


# Extract Functional Diversity Metrics
FD_summary <- data.frame(
  Plot = rownames(species_only),
  FRic = FD_results$FRic,
  FEve = FD_results$FEve,
  FDiv = FD_results$FDiv,
  Species_Richness = FD_results$sing.sp
)

# Visualize the results
library(ggplot2)
library(tidyr)

# Convert FD_summary to long format for visualization
FD_long <- FD_summary %>%
  pivot_longer(cols = FRic:Species_Richness, names_to = "Metric", values_to = "Value")

# Create a bar plot
ggplot(FD_long, aes(x = Plot, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Functional Diversity Metrics per Plot",
    x = "Plot",
    y = "Value",
    fill = "Metric"
  )


# Modeling Relationships ####
# Run  Script FSC-PCA.R
# Add PC1 data to FD_summary
FD_summary$PC1 <- PC1_scores 

# Convert FD_summary to long format for ggplot2
FD_long <- FD_summary %>%
  pivot_longer(cols = FRic:Species_Richness, names_to = "Metric", values_to = "Value")
View(FD_long)

# Species Richness - GLM ####
# Count data
ggplot(FD_summary, aes(x = PC1, y = Species_Richness)) +
  geom_point() +
  geom_smooth(method = "lm", linetype = "dashed", color = "blue") +  # Linear trend
  geom_smooth(method = "loess", linetype = "solid", color = "red") +  # LOESS (non-linear trend)
  theme_minimal() +
  labs(title = "FRic vs. PC1: Checking Linearity", x = "PC1", y = "FRic")
# Not linear ?

# Fit a linear model
lm_Species_Richness <- lm(Species_Richness ~ PC1, data = FD_summary)
# Fit a quadratic model (includes squared PC1 term)
lm_Species_Richness_quad <- lm(Species_Richness ~ poly(PC1, 2), data = FD_summary)
# Perform ANOVA to test for lack of fit
# If the p-value is significant, the quadratic model fits the data better, indicating a non-linear relationship
anova(lm_Species_Richness, lm_Species_Richness_quad)
# no evidence for non-linear relationships.

# # Check for unimodal relationhsip 
# # Add a quadratic term for PC1
# FD_summary$PC1_squared <- FD_summary$PC1^2
# # Run the quadratic model
# lm_quad_Species_Richness <- lm(Species_Richness ~ PC1 + PC1_squared, data = FD_summary)
# # Summarize the model
# summary(lm_quad_Species_Richness)

# Test for distribution of residuals.
shapiro.test(lm_Species_Richness$residuals) # normal residuals 

plot(lm_Species_Richness, which = 1)  # Residuals vs Fitted plot

#heteroscadicity
bptest(lm_Species_Richness)  # Breusch-Pagan test for constant variance
# no significant heteroscedasticity!


#glm
glm_Species_Richness<- glm(Species_Richness ~ PC1, family = Gamma(link = "log"), data = FD_summary)
summary(glm_Species_Richness)

AIC(lm_Species_Richness, glm_Species_Richness)
# stick to lm
summary(lm_Species_Richness)



# FRIC - Linear ####
ggplot(FD_summary, aes(x = PC1, y = FRic)) +
  geom_point() +
  geom_smooth(method = "lm", linetype = "dashed", color = "blue") +  # Linear trend
  geom_smooth(method = "loess", linetype = "solid", color = "red") +  # LOESS (non-linear trend)
  theme_minimal() +
  labs(title = "FRic vs. PC1: Checking Linearity", x = "PC1", y = "FRic")
# Not linear

# Fit a linear model
lm_FRic <- lm(FRic ~ PC1, data = FD_summary)
# Fit a quadratic model (includes squared PC1 term)
lm_FRic_quad <- lm(FRic ~ poly(PC1, 2), data = FD_summary)
# Perform ANOVA to test for lack of fit
# If the p-value is significant, the quadratic model fits the data better, indicating a non-linear relationship
anova(lm_FRic, lm_FRic_quad)
# no evidence for non-linear relationships. Improvement via quadratic is minimal. 

# Test for distribution of residuals.
shapiro.test(lm_FRic$residuals) # non-normal residuals 

# Log-transform FRic
FD_summary$log_FRic <- log(FD_summary$FRic + 1)  # Add 1 to avoid log(0)
# Re-run the linear model with transformed FRic
lm_log_FRic <- lm(log_FRic ~ PC1, data = FD_summary)
summary(lm_log_FRic)

# Check residuals again
shapiro.test(residuals(lm_log_FRic)) # transformation did not improve the model!
# residuals are still not normal 

# Try GLM 
# plot distribution of residuals
plot(lm_log_FRic)
#glm
glm_FRic <- glm(FRic ~ PC1, family = Gamma(link = "log"), data = FD_summary)
summary(glm_FRic)

AIC(lm_log_FRic, glm_FRic)
# lm is better but residuals are still not normal... can'ty use it. 

# Now run a non parametric 
# Non-parametric Alternative: If transformations fail to normalize the residuals, 
# use a Spearman's rank correlation to test the relationship between PC1 and FRic
spear_FRic_noR <- cor.test( FD_noR_summary$FRic,  FD_noR_summary$PC1, method = "spearman")
# no significant relationship.

# Check for unimodal relationship 
# Add a quadratic term for PC1
FD_noR_summary$PC1_squared <- FD_noR_summary$PC1^2

# Run the quadratic model
lm_quad_FRic <- lm(log_FRic ~ PC1 + PC1_squared, data = FD_noR_summary)
# Summarize the model # R2 decrease don't keep. 
summary(lm_quad_FRic)

# Plot the observed data and the quadratic curve
ggplot(FD_noR_summary, aes(x = PC1, y = log_FRic)) +
  geom_point() +  # Scatter plot of data points
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), color = "blue", se = TRUE) +
  theme_minimal() +
  labs(
    title = "Unimodal Relationship Between FSC (PC1) and Functional Diversity (FRic)",
    x = "Forest Structural Complexity (PC1)",
    y = "Log(Functional Richness)"
  )

# Not significant. 
# testing for non linear relationship 
library(mgcv)
gam_FRic <- gam(FRic ~ s(PC1), data = FD_summary)
summary(gam_FRic)
plot(gam_FRic, pages = 1, rug = TRUE, shade = TRUE)

AIC(gam_FRic, lm_FRic) # since the smoothing term is not significant, let stick linearity 
gam.check(gam_FRic)




# FEve - GAM ####
ggplot(FD_summary, aes(x = PC1, y = FEve)) +
  geom_point() +
  geom_smooth(method = "lm", linetype = "dashed", color = "blue") +  # Linear trend
  geom_smooth(method = "loess", linetype = "solid", color = "red") +  # LOESS (non-linear trend)
  theme_minimal() +
  labs(title = "FEve vs. PC1: Checking Linearity", x = "PC1", y = "FRic")
# Not linear

# Fit a linear model
lm_FEve <- lm(FEve ~ PC1, data = FD_summary)

# Fit a quadratic model (includes squared PC1 term)
lm_FEve_quad <- lm(FEve ~ poly(PC1, 2), data = FD_summary)

# Perform ANOVA to test for lack of fit
# If the p-value is significant, the quadratic model fits the data better, indicating a non-linear relationship
anova(lm_FEve, lm_FEve_quad)
# Quadratic model is not better 

shapiro.test(lm_FEve$residuals) # normal residuals 
summary(lm_FEve) # no significant effect 

# Check for unimodal relationship 
# Run the quadratic model
FD_summary$PC1_squared <-   FD_summary$PC1 ^ 2
lm_quad_FEve <- lm(FEve ~ PC1 + PC1_squared, data = FD_summary)
# Summarize the model # R2 decrease don't keep. 
summary(lm_quad_FEve)


# testing for non linear relationship 
library(mgcv)
gam_FEve <- gam(FEve ~ s(PC1), data = FD_summary)
summary(gam_FEve)
plot(gam_FEve, pages = 1, rug = TRUE, shade = TRUE)

AIC(gam_FEve , lm_FEve, lm_quad_FEve)
gam.check(gam_FEve)


# # FDis ####
# ggplot(FD_noR_summary, aes(x = PC1, y = FDis)) +
#   geom_point() +
#   geom_smooth(method = "lm", linetype = "dashed", color = "blue") +  # Linear trend
#   geom_smooth(method = "loess", linetype = "solid", color = "red") +  # LOESS (non-linear trend)
#   theme_minimal() +
#   labs(title = "FDis vs. PC1: Checking Linearity", x = "PC1", y = "FRic")
# # Sorta linear??? 
# # No evidence for non-linear relationships.?
# 
# # Fit a linear model
# lm_FDis <- lm(FDis ~ PC1, data = FD_noR_summary)
# # Fit a quadratic model (includes squared PC1 term)
# lm_FDis_quad <- lm(FDis ~ poly(PC1, 2), data = FD_noR_summary)
# # Perform ANOVA to test for lack of fit
# anova(lm_FDis, lm_FDis_quad)
# # polynomial did not improve model. 
# 
# shapiro.test(lm_FDis$residuals) # normal residuals 
# summary(lm_FDis) #no significnat relationship



# FDiv - LM ####
ggplot(FD_summary, aes(x = PC1, y = FDiv)) +
  geom_point() +
  geom_smooth(method = "lm", linetype = "dashed", color = "blue") +  # Linear trend
  geom_smooth(method = "loess", linetype = "solid", color = "red") +  # LOESS (non-linear trend)
  theme_minimal() +
  labs(title = "FDiv vs. PC1: Checking Linearity", x = "PC1", y = "FRic")
# Not linear

# Fit a linear model
lm_FDiv <- lm(FDiv ~ PC1, data = FD_summary)
# Fit a quadratic model (includes squared PC1 term)
lm_FDiv_quad <- lm(FDiv ~ poly(PC1, 2), data = FD_summary)
# Fit a cubic model
lm_FDiv_cubic <- lm(FDiv ~ poly(PC1, 3), data = FD_summary)

# Compare models
anova(lm_FDiv,lm_FDiv_cubic)
# Perform ANOVA to test for lack of fit
anova(lm_FDiv, lm_FDiv_quad) 
# No evidence for non-linear relationships. stick to linear. 
summary(lm_FDiv)

# testing for non linear relationship 
library(mgcv)
gam_FDiv<- gam(FDiv ~ s(PC1), data = FD_summary)
summary(gam_FDiv)
plot(gam_FDiv, pages = 1, rug = TRUE, shade = TRUE) # gam fits a lienar model. 

AIC(gam_FDiv , lm_FDiv)
gam.check(gam_FEve) 




# Model used for FD indices

# Species richness
summary(glm_Species_Richness)
# Calculate pseudo R-squared
null_model <- glm(Species_Richness ~ 1, family = Gamma(link = "log"), data = FD_summary)
pseudo_r2 <- 1 - (logLik(glm_Species_Richness) / logLik(null_model))
pseudo_r2_value <- as.numeric(pseudo_r2)

# Proportion of deviance explained
deviance_explained <- 1 - (glm_Species_Richness$deviance / glm_Species_Richness$null.deviance)

# FRic
summary(lm_FRic) 
#FEve
summary(lm_FEve)
#FDiv
summary(lm_FDiv)


# Results table ####
results_summary <- data.frame(
  Index = c("Species_Richness", "FRic", "FEve", "FDiv"),
  Model = c("GLM", "LM", "LM", "LM"),
  p_value = c(
    summary(glm_Species_Richness)$coefficients[2, 4],
    summary(lm_FRic)$coefficients[2, 4],
    summary(lm_FEve)$coefficients[2, 4],
    summary(lm_FDiv)$coefficients[2, 4]
  ),
  R_squared = c(
    deviance_explained,  # Or use pseudo_r2_value here
    summary(lm_FRic)$r.squared,
    summary(lm_FEve)$r.squared,
    summary(lm_FDiv)$r.squared
  )
)


print(results_summary)

# Make a flextable object
ft <- flextable(results_summary)

# Format the table in APA style
ft <- ft %>%
  set_header_labels(Index = "Index", Model = "Model", 
                    p_value = "P-value", R_squared = "R²") %>%
  autofit() %>%
  border_remove() %>%
  border_outer(part = "all", border = fp_border(color = "black", width = 1)) %>%
  hline_top(border = fp_border(color = "black", width = 1)) %>%
  hline_bottom(border = fp_border(color = "black", width = 1)) %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  fontsize(size = 11, part = "all")
save_as_docx(ft, path = "LM_PC1_FD_RFU_Table_APA.docx")



# PLOT ####
# Add PC1 data to FD_summary
FD_summary$PC1 <- PC1_scores
# Convert FD_summary to long format for ggplot2
FD_long <- FD_summary %>%
  pivot_longer(cols = FRic:Species_Richness, names_to = "Metric", values_to = "Value")
View(FD_long)
unique(FD_long$Metric)

# Update Metric Labels
FD_long$Metric <- factor(FD_long$Metric,
                             levels = c("Species_Richness", "FRic", "FEve", "FDiv" #, "FDis"
                             ),
                             labels = c("Species Richness",
                                        "Functional Richness (FRic)",
                                        "Functional Evenness (FEve)",
                                        "Functional Divergence (FDiv)"
                                        #"Functional Dispersion (FDis)"
                             ))

# Add a column to specify model type based on the results table
FD_long <- FD_long %>%
  mutate(ModelType = case_when(
    Metric == "Species Richness" ~ "Generalized Linear",
    Metric == "Functional Richness (FRic)" ~ "Linear",
    Metric == "Functional Evenness (FEve)" ~ "Linear",
    Metric == "Functional Divergence (FDiv)" ~ "Linear"
  ))

# Add a text label with P-values (example with hypothetical results)
FD_long <- FD_long  %>%
  mutate(P_value = case_when(
    Metric == "Species Richness" ~ 0.240,
    Metric == "Functional Richness (FRic)" ~ 0.351,
    Metric == "Functional Evenness (FEve)" ~ 0.627,
    Metric == "Functional Divergence (FDiv)" ~ 0.099
    #,Metric == "FDis" ~ 0.48
  ))

FD_long  <- FD_long  %>%
  mutate(Signif = case_when(
    P_value < 0.05 ~ "*",       # Significant
    P_value < 0.1 ~ ".",        # Marginally significant
    TRUE ~ "ns"                 # Not significant
  ))



# Check and clean data
# Ensure Metric is a factor
FD_long$Metric <- as.factor(FD_long$Metric)

# Remove rows where Metric is "FDis"
FD_long <- FD_long %>%
  filter(Metric != "FDis")

# Drop unused levels
FD_long$Metric <- droplevels(FD_long$Metric)

# Check the updated table
table(FD_long$Metric, useNA = "ifany")

# remove all rows for Function Distance. 
# Using dplyr
# FD_noR_long <- FD_noR_long %>% drop_na(ModelType)
# colnames(FD_noR_long)
# View(FD_noR_long)

library(ggplot2)
library(viridis)

# Separate the data for linear and spearman models
linear_data <- FD_long %>% filter(ModelType == "Linear")
glm_data <- FD_long %>% filter(ModelType == "Generalized Linear")



combined_plot <- ggplot(FD_long, aes(x = PC1, y = Value, color = Metric)) +
  # Add points for all data
  geom_point(size = 3, alpha = 0.8) +
  
  # Add linear trend lines
  geom_smooth(
    aes(linetype = ModelType),
    method = "lm", formula = y ~ x, se = TRUE,
    data = FD_long %>% filter(ModelType == "Linear"),
    color = "black", linewidth = 1.2
  ) +
  
  # Add GLM trend lines
  geom_smooth(
    aes(linetype = ModelType),
    method = "lm", formula = y ~ x, se = TRUE,
    data = FD_long %>% filter(ModelType == "Generalized Linear"),
    color = "black", linewidth = 1.2
  ) +
  
  # Add significance labels for linear models
  geom_text(
    data = linear_data, 
    aes(x = Inf, y = Inf, label = Signif),
    hjust = 1.2, vjust = 1.2, size = 5, color = "black", fontface = "bold", inherit.aes = FALSE
  ) +
  
  # Add significance labels for Spearman models
  geom_text(
    data = glm_data, 
    aes(x = Inf, y = Inf, label = Signif),
    hjust = 1.2, vjust = 1.2, size = 5, color = "black", fontface = "bold", inherit.aes = FALSE
  ) +
  
  # Facet the plot by Metric
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  
  # Add titles and axis labels
  labs(
    title = "Functional Diversity Metrics and Forest Structural Complexity",
    x = "Forest Structural Complexity (PC1)",
    y = "Functional Diversity Metrics",
    linetype = "Model Type",
    color = "Metric"
  ) +
  
  # Use modern color palette
  scale_color_viridis_d(option = "plasma") +
  
  # Refine theme for better readability
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 12, color = "white"),
    strip.background = element_rect(fill = "gray30"),
    legend.position = "right",
    legend.box = "vertical",
    panel.grid.major = element_line(color = "gray85", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  )

# Print the plot
print(combined_plot)

# Save the plot in A4 size
ggsave("combined_plot_A4.png", combined_plot, width = 8.27, height = 11.69, units = "in", dpi = 300)



# FD - Ticks relationships ####

# Nymphs
# Regression models for count data, such as generalized linear models (GLMs) or 
# generalized linear mixed models (GLMMs).

# Larvae and Adults
# Presence/Absence Data:
# If your data is binary (0 = absence, 1 = presence), use logistic regression (binary GLMs/GLMMs).

library(lme4)
library(MASS)
library(car)
library(broom.mixed)
library(dplyr)
library(knitr)

# Setting the working directory
setwd("/Users/sacharoland/Desktop/WUR/THESIS/FSC and vertebrates- Romania/Data Analysis/Data")
ticks <- read.csv("Ticks_data.csv")
view(ticks)


# Aggregate data to plot level
plot_data <- ticks %>%
  group_by(Plotname) %>%
  summarise(
    Mean_LL = mean(LL_thickness, na.rm = TRUE),
    Mean_SD = mean(Saturation_Deficit, na.rm = TRUE),
    Larvae_Presence = mean(Larvae > 0),  # Proportion of transects with larvae
    Nymphs_Count = sum(Nymphs),          # Total nymph count per plot
    Adults_Presence = mean(Adults > 0),  # Proportion of transects with adults
    .groups = 'drop'
  )

#LowerArcer does not have data for LL_thickness. 

# Ensure response variables are binary (0 or 1)
plot_data <- plot_data %>%
  mutate(
    Larvae_Presence = ifelse(Larvae_Presence > 0, 1, 0),  # Convert proportions to binary
    Adults_Presence = ifelse(Adults_Presence > 0, 1, 0)   # Convert proportions to binary
  )

view(plot_data) # looks ok... 

# Add FSC scores to the aggregated data
plot_data <- plot_data %>%
  mutate(
    PC1_scores = PC1_scores,  # Assuming PC1_scores are in the same order as plots
    PC2_scores = PC2_scores
  )

# Check column names; both need col Plotname.
colnames(plot_data)


# Get columns for FD_summary and ticks in order. 
# FD_summary <- FD_summary %>% rename(Plotname = Plot)
colnames(FD_summary)

# Check for matching plot names
common_plots <- intersect(plot_data$Plotname, FD_summary$Plotname)
print(common_plots) #yes

# Merge DF 
FD_summary_ticks <- plot_data %>%
  left_join(FD_summary, by = "Plotname")
colnames(FD_summary_ticks)

# MODELS WITH SP. RICHNESS AS COUNT ###########################################
# Standardize the functional diversity indices (mean = 0; sd = 1)
# ensures that all predictors have equal weight numerically, helping the model converge more easily.
FD_summary_ticks <- FD_summary_ticks %>%
  mutate(
    # Species_Richness = scale(Species_Richness), #eventhough it is count data? 
    # do with and without scaling species_richness; 
    FRic_scaled = scale(FRic),
    FEve_scaled = scale(FEve),
    FDiv_scaled = scale(FDiv)
  )

view(FD_summary_ticks)

#Consider fitting the model both with scaled and unscaled versions of Species_Richness 
# to check how the scaling affects your results. If the model results differ substantially, 
# interpret why this might be the case.


colnames(FD_summary_ticks)

# check for multicollinearity
vif_check <- lm(FRic_scaled ~ FEve_scaled + FDiv_scaled, data = FD_summary_ticks)
vif(vif_check)  # VIF > 5 indicates high multicollinearity

# Example: Mixed model for larvae abundance
larvae_model_traits <- glmer(Larvae_Presence ~ FRic_scaled + FEve_scaled + FDiv_scaled + Species_Richness + (1 | Plotname), 
                             family = binomial, data = FD_summary_ticks)
summary(larvae_model_traits)
# change if Plotname 
simpler_model <- glm(Larvae_Presence ~ FRic_scaled + FEve_scaled + FDiv_scaled + Species_Richness, 
                     family = binomial, data = FD_summary_ticks)
AIC(larvae_model_traits, simpler_model) # Plotname does improve the model significantly.  




# Example: Mixed model for nymph abundance
nymph_model_traits <- glmmTMB(Nymphs_Count ~ FRic_scaled + FEve_scaled +  FDiv_scaled + (1 | Plotname), 
                              family = nbinom2, data = FD_summary_ticks)

#Warning message:
#In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#Model convergence problem; non-positive-definite Hessian matrix.

vif_check <- lm(Nymphs_Count ~ FRic_scaled + FEve_scaled + FDiv_scaled + Species_Richness, data = FD_summary_ticks)
vif(vif_check) #no multicollinearity. 

# juts try to make the model more simple. 
nymph_model_traits <- glmmTMB(Nymphs_Count ~ FRic_scaled + FEve_scaled +  (1 | Plotname), 
                              family = nbinom2, data = FD_summary_ticks)
summary(nymph_model_traits)
plot(residuals(nymph_model_traits, type = "pearson")) #looks ok. 


nymph_model_traits <- glmmTMB(Nymphs_Count ~ FRic_scaled + FDiv_scaled + (1 | Plotname), 
                              family = nbinom2, data = FD_summary_ticks)
summary(nymph_model_traits)
plot(residuals(nymph_model_traits, type = "pearson")) #looks ok. 


# incloduing species richness --> run into issues. 
nymph_model_traits <- glmmTMB(Nymphs_Count ~ FRic_scaled + Species_Richness + (1 | Plotname), 
                              family = nbinom2, data = FD_summary_ticks)
summary(nymph_model_traits)
plot(residuals(nymph_model_traits, type = "pearson")) #looks ok. 


# FRic_scaled has an effect on nymphs count
# not FEve or FDiv


# Example: Mixed model for adults
adults_model_traits <- glmer(Adults_Presence ~ FRic_scaled + FEve_scaled + FDiv_scaled + (1 | Plotname), 
                             family = binomial, data = FD_summary_ticks)
summary(adults_model_traits)
VarCorr(adults_model_traits) # Plotname is very small (or exactly zero), it suggests 
# that Plotname may not contribute much variability to the model.

#simpler model 
adults_model_traits <- glm(Adults_Presence ~ FRic_scaled + FEve_scaled + FDiv_scaled + Species_Richness, 
                              family = binomial, data = FD_summary_ticks)

summary(adults_model_traits)


# CHECKING ASSUMPTIONS
# For mixed models, assess residuals and random effects
# nymphs
plot(residuals(nymph_model_traits, type = "pearson"))
qqnorm(residuals(nymph_model_traits))
qqline(residuals(nymph_model_traits))


# FOR fixed effects, evaluate multicollinearity:
vif_check <- lm(Nymphs_Count ~ FRic_scaled + FEve_scaled + FDiv_scaled + Species_Richness, data = FD_summary_ticks)
vif(vif_check) #OOKK





# Summarize the results
summary(larvae_model_traits) #No significnat results
summary(nymph_model_traits) # FRich is negatively associated with nymphs counts! 
summary(adults_model_traits) # FEve is positvely associated with Adults presence 




# Make nice results tables
# Tidy summaries for each model
larvae_summary <- tidy(larvae_model_traits, effects = "fixed") %>%
  mutate(Response = "Larvae", Model = "Binomial GLMM")

nymphs_summary <- tidy(nymph_model_traits, effects = "fixed") %>%
  mutate(Response = "Nymphs", Model = "Negative Binomial GLMM")

adults_summary <- tidy(adults_model_traits, effects = "fixed") %>%
  mutate(Response = "Adults", Model = "Binomial GLMM")

# Combine all results into one table
combined_summary <- bind_rows(larvae_summary, nymphs_summary, adults_summary) %>%
  filter(term != "(Intercept)") %>%  # Exclude intercept for clarity
  select(Response, Model, term, estimate, std.error, p.value) %>%
  rename(
    Predictor = term,
    Estimate = estimate,
    StdError = std.error,
    PValue = p.value
  )

# Add significance stars
combined_summary <- combined_summary %>%
  mutate(Significance = case_when(
    PValue < 0.001 ~ "***",
    PValue < 0.01 ~ "**",
    PValue < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Print the table in a clean format
kable(combined_summary, caption = "Summary of Model Results for Tick Life Stages")



# MODELS WITH SP. RICHNESS AS CONTINOUS  ###########################################
# Standardize the functional diversity indices (mean = 0; sd = 1)
# ensures that all predictors have equal weight numerically, helping the model converge more easily.
FD_summary_ticks <- FD_summary_ticks %>%
  mutate(
    Species_Richness_scaled = scale(Species_Richness),
    FRic_scaled = scale(FRic),
    FEve_scaled = scale(FEve),
    FDiv_scaled = scale(FDiv)
  )

view(FD_summary_ticks)

#Consider fitting the model both with scaled and unscaled versions of Species_Richness 
# to check how the scaling affects your results. If the model results differ substantially, 
# interpret why this might be the case.


colnames(FD_summary_ticks)

# check for multicollinearity
vif_check <- lm(FRic_scaled ~ FEve_scaled + FDiv_scaled, data = FD_summary_ticks)
vif(vif_check)  # VIF > 5 indicates high multicollinearity

# Example: Mixed model for larvae abundance
larvae_model_traits <- glmer(Larvae_Presence ~ FRic_scaled + FEve_scaled + FDiv_scaled + Species_Richness + (1 | Plotname), 
                             family = binomial, data = FD_summary_ticks)
summary(larvae_model_traits)
# change if Plotname 
simpler_model <- glm(Larvae_Presence ~ FRic_scaled + FEve_scaled + FDiv_scaled + Species_Richness, 
                     family = binomial, data = FD_summary_ticks)
AIC(larvae_model_traits, simpler_model) # Plotname does improve the model significantly.  




# Example: Mixed model for nymph abundance
nymph_model_traits <- glmmTMB(Nymphs_Count ~ FRic_scaled + FEve_scaled +  FDiv_scaled + (1 | Plotname), 
                              family = nbinom2, data = FD_summary_ticks)

#Warning message:
#In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#Model convergence problem; non-positive-definite Hessian matrix.

vif_check <- lm(Nymphs_Count ~ FRic_scaled + FEve_scaled + FDiv_scaled + Species_Richness, data = FD_summary_ticks)
vif(vif_check) #no multicollinearity. 

# juts try to make the model more simple. 
nymph_model_traits <- glmmTMB(Nymphs_Count ~ FRic_scaled + FEve_scaled +  (1 | Plotname), 
                              family = nbinom2, data = FD_summary_ticks)
summary(nymph_model_traits)
plot(residuals(nymph_model_traits, type = "pearson")) #looks ok. 


nymph_model_traits <- glmmTMB(Nymphs_Count ~ FRic_scaled + FDiv_scaled + (1 | Plotname), 
                              family = nbinom2, data = FD_summary_ticks)
summary(nymph_model_traits)
plot(residuals(nymph_model_traits, type = "pearson")) #looks ok. 


# incloduing species richness --> run into issues. 
nymph_model_traits <- glmmTMB(Nymphs_Count ~ FRic_scaled + Species_Richness + (1 | Plotname), 
                              family = nbinom2, data = FD_summary_ticks)
summary(nymph_model_traits)
plot(residuals(nymph_model_traits, type = "pearson")) #looks ok. 


# FRic_scaled has an effect on nymphs count
# not FEve or FDiv


# Example: Mixed model for adults
adults_model_traits <- glmer(Adults_Presence ~ FRic_scaled + FEve_scaled + FDiv_scaled + (1 | Plotname), 
                             family = binomial, data = FD_summary_ticks)
summary(adults_model_traits)
VarCorr(adults_model_traits) # Plotname is very small (or exactly zero), it suggests 
# that Plotname may not contribute much variability to the model.

#simpler model 
adults_model_traits <- glm(Adults_Presence ~ FRic_scaled + FEve_scaled + FDiv_scaled + Species_Richness, 
                           family = binomial, data = FD_summary_ticks)

summary(adults_model_traits)


# CHECKING ASSUMPTIONS
# For mixed models, assess residuals and random effects
# nymphs
plot(residuals(nymph_model_traits, type = "pearson"))
qqnorm(residuals(nymph_model_traits))
qqline(residuals(nymph_model_traits))


# FOR fixed effects, evaluate multicollinearity:
vif_check <- lm(Nymphs_Count ~ FRic_scaled + FEve_scaled + FDiv_scaled + Species_Richness, data = FD_summary_ticks)
vif(vif_check) #OOKK


# Summarize the results
summary(larvae_model_traits) #No significnat results
summary(nymph_model_traits) # FRich is negatively associated with nymphs counts! 
summary(adults_model_traits) # FEve is positvely associated with Adults presence 


# Make nice results tables
# Tidy summaries for each model
larvae_summary <- tidy(larvae_model_traits, effects = "fixed") %>%
  mutate(Response = "Larvae", Model = "Binomial GLMM")

nymphs_summary <- tidy(nymph_model_traits, effects = "fixed") %>%
  mutate(Response = "Nymphs", Model = "Negative Binomial GLMM")

adults_summary <- tidy(adults_model_traits, effects = "fixed") %>%
  mutate(Response = "Adults", Model = "Binomial GLMM")

# Combine all results into one table
combined_summary <- bind_rows(larvae_summary, nymphs_summary, adults_summary) %>%
  filter(term != "(Intercept)") %>%  # Exclude intercept for clarity
  select(Response, Model, term, estimate, std.error, p.value) %>%
  rename(
    Predictor = term,
    Estimate = estimate,
    StdError = std.error,
    PValue = p.value
  )

# Add significance stars
combined_summary <- combined_summary %>%
  mutate(Significance = case_when(
    PValue < 0.001 ~ "***",
    PValue < 0.01 ~ "**",
    PValue < 0.05 ~ "*",
    TRUE ~ ""
  ))

# Print the table in a clean format
kable(combined_summary, caption = "Summary of Model Results for Tick Life Stages")



