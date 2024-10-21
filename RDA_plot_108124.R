library(vegan)
library(phyloseq)
library(ggplot2)
library(ggvegan)
# Install ggvegan if not already installed
#devtools::install_github("gavinsimpson/ggvegan")
setwd("C:/Users/Lindsay Alma/Dropbox/Eelgrass/West Coast TagSeq primers/FullDataset")
# Load the phyloseq object
ps <- readRDS("ps_18S_lulu.RDS")

# Filter out specific persons
ps <- subset_samples(ps, Person != "Olivia" & Person != "CB")

# Subset out unwanted domains and kingdoms
ps <- subset_taxa(ps, is.na(domain) | domain != "Bacteria")
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Chloroplastida")
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Animalia")
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Alveolata")
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Amoebozoa")
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Holozoa")
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Rhodophyceae")
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Diatomea")
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Ochrophyta")

# Subset samples by Lesion Type
ps <- subset_samples(ps, Lesion_Type == "Crescent" | Lesion_Type == "Not crescent")

# Read the new ASV matrix from CSV
new_asv_data <- read.csv("otu_normalized_base_72624.csv", row.names = 1)

# Convert to matrix
new_asv_matrix <- as.matrix(new_asv_data)

# Update the otu_table with the new matrix
otu_table(ps) <- otu_table(new_asv_matrix, taxa_are_rows = TRUE)

# Extract and transpose the OTU table
otu_table_matrix <- as(otu_table(ps), "matrix")
otu_table_transposed <- t(otu_table_matrix)

#################################################################################33

sum(is.na(sample_data_ps$State))  # Check how many NAs in State

# Filter sample_data_ps to remove samples with missing State
sample_data_ps_clean <- sample_data_ps[!is.na(sample_data_ps$State), ]

# Ensure row names match between OTU table and sample data
otu_samples <- rownames(otu_table_clean)
meta_samples <- rownames(sample_data_ps_clean)

# Find the intersection of row names (samples that are present in both)
common_samples <- intersect(otu_samples, meta_samples)

# Subset the OTU table and sample data to include only the common samples
otu_table_clean <- otu_table_clean[common_samples, ]
sample_data_clean <- sample_data_ps_clean[common_samples, ]

# Extract the cleaned 'State' variable
state_clean <- as.factor(sample_data_clean$State)

# Perform RDA with cleaned OTU data and State variable
rda_result <- rda(otu_table_clean ~ state_clean)

# Summarize and plot the RDA result
summary(rda_result)

#########################################################################################

# List of specified ASVs to include
selected_asvs <- c("ASV_18S_1", "ASV_18S_2", "ASV_18S_4", "ASV_18S_5", "ASV_18S_6")

# Extract sample and species scores from the RDA result
rda_scores <- scores(rda_result, display = c("sites", "species"))

# Create a data frame for the samples (sites)
samples_df <- as.data.frame(rda_scores$sites)
samples_df$State <- state_clean  # Add the State variable to the sample data

# Create a data frame for the species (ASVs)
species_df <- as.data.frame(rda_scores$species)
species_df$ASV <- rownames(species_df)  # Add ASV names

# Subset the species data to include only the selected ASVs
species_df_filtered <- species_df[species_df$ASV %in% selected_asvs, ]

# Remove "_18S" from ASV names for labeling
species_df_filtered$ASV <- gsub("_18S", "", species_df_filtered$ASV)

# Calculate the percentage of variance explained by each axis
explained_variance <- summary(rda_result)$cont$importance[2, 1:2] * 100
axis_labels <- paste0("RDA", 1:2, " (", round(explained_variance, 1), "%)")

# Create the ggplot
p <- ggplot() +
  # Plot samples (sites) colored by State
  geom_point(data = samples_df, aes(x = RDA1, y = RDA2, color = State), size = 3) +
  # Add arrows for ASVs (species)
  geom_segment(data = species_df_filtered, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linetype = "dashed") +
  # Add ASV labels without "_18S"
  geom_text(data = species_df_filtered, aes(x = RDA1, y = RDA2, label = ASV), color = "black", hjust = -0.2, vjust = 0.5) +
   # Customize plot
  labs(title = "RDA of Selected ASVs by State", x = axis_labels[1], y = axis_labels[2], color = "State") +
  scale_color_manual(values = rainbow(length(unique(samples_df$State)))) +  # Custom colors for States
  theme_minimal()  # Clean theme

# Display the plot
print(p)


