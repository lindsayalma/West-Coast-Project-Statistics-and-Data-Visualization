library(readr)
library(tidymodels)
library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)


ps <- readRDS("ps_18S_lulu.RDS")

ps <- subset_samples(ps, Person != "Olivia" & Person != "CB")


# Subset out unwanted domains and kingdoms
ps <- subset_taxa(ps, is.na(domain) | domain != "Bacteria")         # Remove Bacteria
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Chloroplastida")  # Remove Chloroplastida
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Animalia")        # Remove Animalia
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Alveolata")       # Remove Alveolata
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Amoebozoa")       # Remove Amoebozoa
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Holozoa")         # Remove Holozoa
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Rhodophyceae")    # Remove Rhodophyceae
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Diatomea")        # Remove Diatomea
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Ochrophyta")      # Remove Ochrophyta

ps <- subset_samples(ps, Lesion_Type=="Crescent" | Lesion_Type=="Not crescent")

# Read the new ASV matrix from CSV
new_asv_data <- read.csv("otu_normalized_base_72624.csv", row.names = 1)
print(new_asv_data)

# Convert to matrix
new_asv_matrix <- as.matrix(new_asv_data)

# Assuming your phyloseq object is named `ps`
# Update the otu_table with the new matrix
ps_new <- ps
otu_table(ps_new) <- otu_table(new_asv_matrix, taxa_are_rows = TRUE)

meta<- sample_data(ps)
ASV<-otu_table(ps_new)
####################################################################################


# Assuming your phyloseq object is named `ps_new`
# Get the ASV counts and sample metadata
asv_counts <- as.data.frame(otu_table(ps_new))
meta <- as.data.frame(sample_data(ps_new))

# Transpose the ASV counts
asv_counts_transposed <- t(asv_counts)

# Combine metadata and ASV counts
combined_data <- cbind(meta, asv_counts_transposed)

# Gather the data to long format
long_data <- combined_data %>%
  pivot_longer(cols = starts_with("ASV"), names_to = "ASV", values_to = "Count")

# Filter for State == "AK"
long_data_ak <- long_data %>%
  filter(State == "AK")#########################################################################################################

# Group by Site and ASV, then summarize to get counts for each ASV at each sample
dominant_asvs_ak <- long_data_ak %>%
  group_by(Site, ASV) %>%
  summarize(Dominance = sum(Count), .groups = 'drop') %>%  # Sum the counts for dominance
  group_by(Site) %>%
  slice_max(Dominance, n = 1, with_ties = FALSE) %>%  # Get the single dominant ASV for each site
  ungroup() %>%
  arrange(Site)  # Sort by Site

# Print the dominant ASVs at each site
print(dominant_asvs_ak)

# Ensure Site is a factor for proper ordering
dominant_asvs_ak$Site <- factor(dominant_asvs_ak$Site)

# Create the plot
ggplot(dominant_asvs_ak, aes(x = Site, y = Dominance, fill = ASV)) +  # Use Site directly
  geom_bar(stat = "identity", position = "dodge") +
  geom_vline(xintercept = seq(1.5, length(unique(dominant_asvs_ak$Site)) - 0.5), color = "black", linetype = "solid") +  # Vertical lines
  labs(title = "Dominant ASV by Sample (Alaska)", x = "Site", y = "Dominance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

####################################################################################

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(phyloseq)  # Only if you haven't already loaded phyloseq
library(tibble)    # Needed for rownames_to_column

# Convert metadata and ASV data
meta_df <- data.frame(sample_data(meta))  # Convert phyloseq metadata to data frame
ASV_df <- as.data.frame(otu_table(ASV))  # Convert phyloseq ASV data to data frame
ASV_df <- rownames_to_column(ASV_df, var = "ASV")  # Reset row names to a column

# Transpose the ASV data for easier manipulation
ASV_long <- ASV_df %>%
  pivot_longer(cols = -ASV, names_to = "Samples", values_to = "Percentage")


# Filter metadata for Alaska samples and add row names as a column
meta_ak <- meta_df %>%
  filter(State == "SD") %>%################################################################################################
  rownames_to_column(var = "Samples")  # Add row names as a new column called 'Samples'

# Join ASV data with metadata for Alaska
ASV_meta_ak <- ASV_long %>%
  inner_join(meta_ak, by = "Samples")  # Ensure you join on the correct column

# Identify the dominant ASV for each sample
dominant_asvs_ak <- ASV_meta_ak %>%
  group_by(Sample) %>%
  slice_max(order_by = Percentage, n = 1, with_ties = FALSE) %>%
  ungroup()  # Ungroup for further processing


# Assuming ASV_meta_ak has columns: Sample, ASV, Percentage

# Step 1: Group by Sample and find the ASV with the highest percentage for each sample
dominant_asvs_ak <- ASV_meta_ak %>%
  group_by(Sample) %>%
  slice_max(order_by = Percentage, n = 1, with_ties = FALSE) %>%  # Keep only the ASV with the highest percentage
  ungroup()  # Ungroup the data for further use

# Step 2: View the resulting dominant ASV dataframe
print(dominant_asvs_ak)

# Remove rows with NA in the Percentage column
dominant_asvs_ak <- dominant_asvs_ak %>%
  filter(!is.na(Percentage))  # Keep only rows where Percentage is not NA

# Step 2: View the resulting dominant ASV dataframe
print(dominant_asvs_ak)

# Remove rows with NA or 0 in the Percentage column
dominant_asvs_ak <- dominant_asvs_ak %>%
  filter(!is.na(Percentage) & Percentage != 0)  # Keep rows where Percentage is not NA and not 0

# Step 2: View the resulting dominant ASV dataframe
print(dominant_asvs_ak)

# Remove the string "_18S_" from the ASV column in dominant_asvs_ak
dominant_asvs_ak$ASV <- gsub("_18S_", "", dominant_asvs_ak$ASV)

# Check the updated dataframe
print(dominant_asvs_ak)

# Step 1: Reorder the Sample column based on the Site column
dominant_asvs_ak <- dominant_asvs_ak %>%
  arrange(Site)

# Create a new column that combines Sample and Site
dominant_asvs_ak$Sample_Site <- paste(dominant_asvs_ak$Sample, "(", dominant_asvs_ak$Site, ")", sep = "")

# Define a named vector of specific colors for ASVs
asv_colors <- c(
  "ASV1" = "#1f77b4",  # Example: ASV1 is blue
  "ASV2" = "#ff7f0e",  # Example: ASV2 is orange
  "ASV4" = "lightskyblue",  # Example: ASV4 is red
  "ASV5" = "#9467bd",   # Example: ASV5 is purple
  "ASV6" = "lightgreen",   # Example: ASV5 is purple
  "ASV8" = "gold4",   # Example: ASV5 is purple
  "ASV9" = "hotpink2",   # Example: ASV5 is purple
  "ASV10" = "aquamarine",   # Example: ASV5 is purple
  "ASV11" = "aquamarine4",   # Example: ASV5 is purple
  "ASV13" = "azure4",   # Example: ASV5 is purple
  "ASV14" = "bisque3",   # Example: ASV5 is purple
  "ASV16" = "olivedrab",   # Example: ASV5 is purple
  "ASV23" = "gold1",   # Example: ASV5 is purple
  "ASV26" = "gray",   # Example: ASV5 is purple
  "ASV27" = "darkorange3",   # Example: ASV5 is purple
  "ASV28" = "blueviolet",   # Example: ASV5 is purple
  "ASV31" = "goldenrod2",   # Example: ASV5 is purple
  "ASV41" = "lightblue4",   # Example: ASV5 is purple
  "ASV24" = "yellowgreen",   # Example: ASV5 is purple
  "ASV46" = "darkmagenta",   # Example: ASV5 is purple
  "ASV68" = "dodgerblue",   # Example: ASV5 is purple
  "ASV74" = "lightslateblue",   # Example: ASV5 is purple
  "ASV106" = "yellow",   # Example: ASV5 is purple
  "ASV175" = "blue3",   # Example: ASV5 is purple
  "ASV1574" = "darkblue",   # Example: ASV5 is purple
  "ASV2494" = "khaki"   # Example: ASV5 is purple
)

# Create the bar plot with specific ASV colors
ggplot(dominant_asvs_ak, aes(x = reorder(Sample, Sample), y = Percentage, fill = ASV)) + 
  geom_bar(stat = "identity") + 
  labs(title = "Percentage of Dominant ASV by Sample San Diego", x = "", y = "Percentage") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = asv_colors)  # Use specific ASV colors
