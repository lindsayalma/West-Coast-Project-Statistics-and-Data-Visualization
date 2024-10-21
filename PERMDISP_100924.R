library(vegan)
library(phyloseq)
librar(dplyr)
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



##################################################################33

# Extract ASV data and metadata from phyloseq object
asv_data <- otu_table(ps)  # Extract OTU table (ASV matrix)
metadata <- sample_data(ps)  # Extract metadata

# Convert ASV data to matrix if it is not already
asv_matrix <- as.matrix(asv_data)

# Transpose the ASV matrix
asv_matrix <- t(asv_matrix)

# Check and align sample names
asv_samples <- rownames(asv_matrix)
metadata_samples <- rownames(metadata)

# Identify samples in ASV matrix not in metadata
missing_in_metadata <- setdiff(asv_samples, metadata_samples)

# Identify samples in metadata but not in ASV matrix
missing_in_asv <- setdiff(metadata_samples, asv_samples)

if (length(missing_in_metadata) > 0 | length(missing_in_asv) > 0) {
  stop("Sample names in ASV matrix do not match metadata.")
}

# Subset metadata to only include samples present in ASV matrix
metadata <- metadata[asv_samples, ]

#################################################################3
asv_matrix_clean <- na.omit(asv_matrix)
distance_matrix <- vegdist(asv_matrix_clean, method = "bray")
asv_matrix[is.na(asv_matrix)] <- 0

# Prepare grouping variable
#put the variable here to change from satte to GvsL to crescent
group <- metadata$Lesion_Type  # Replace `group` with your actual grouping variable


# Ensure the grouping variable is a factor and matches the number of samples
if (length(group) != nrow(asv_matrix)) {
  stop("The length of the grouping variable does not match the number of samples in the ASV matrix.")
}

group <- factor(group)  # Convert to factor if not already
any(is.na(group))
length(group) == nrow(asv_matrix)
any(is.na(distance_matrix))

# Identify valid samples (no missing data)
valid_samples <- !is.na(group) & complete.cases(asv_matrix)

# Subset the distance matrix and group
distance_matrix_subset <- as.dist(as.matrix(distance_matrix)[valid_samples, valid_samples])
group_subset <- group[valid_samples]


# Identify rows (samples) that are all zeros
empty_rows <- rowSums(asv_matrix_clean) == 0
sum(empty_rows)  # Number of empty samples removed
# Remove empty rows
asv_matrix_clean <- asv_matrix_clean[!empty_rows, ]
distance_matrix <- vegdist(asv_matrix_clean, method = "bray")


# Run betadisper with the subset distance matrix and group
permdisp_result <- betadisper(distance_matrix_subset, group_subset)

# View results
summary(permdisp_result)  # Summary of PERMDISP results

# Set margins and plot size (adjust as necessary)
par(mar = c(5, 5, 2, 2))  # Set margins (bottom, left, top, right)
plot(permdisp_result)


#####################################################################3

# Plot distances to centroids
plot(permdisp_result$distances, col = group, pch = 19, 
     xlab = "Sample", ylab = "Distance to Centroid", 
     main = "Distances to Group Centroids",
     xlim = c(1, length(permdisp_result$distances)))  # Adjust x-axis limits if needed

# Adjust plot margins to make space for legend
par(mar = c(5, 7, 7, 9) + 0.1)

# Add legend outside the plot (right side)
# The `x` and `y` arguments are relative to the plotting region and will need to be adjusted manually
legend("topright", legend = levels(group), col = 1:length(levels(group)), pch = 19, 
       inset = c(-0.5, 0),  # Adjust `inset` to position the legend further right
       xpd = TRUE)  # Allow plotting outside the plot region


# Perform PERMDISP test (ANOVA-like)
anova_result <- anova(permdisp_result)
print(anova_result)


#RESULTS
#G vs L p=0.01308
# Lesion_Type p=0.1756
#crescent vs not p=0.5879


###################################Lesion_Types alone#####################################33



# Extract ASV data and metadata from phyloseq object
asv_data <- otu_table(ps)  # Extract OTU table (ASV matrix)
metadata <- sample_data(ps)  # Extract metadata

# Convert ASV data to matrix if it is not already
asv_matrix <- as.matrix(asv_data)

# Transpose the ASV matrix
asv_matrix <- t(asv_matrix)

# Filter metadata to include only rows where Lesion_Type is "AK"
metadata_ak <- metadata[metadata$Lesion_Type == "WA", ]

# Filter ASV matrix to include only samples from the metadata with Lesion_Type "AK"
asv_matrix <- asv_matrix[rownames(asv_matrix) %in% rownames(metadata_ak), ]

# Re-align metadata to match the filtered ASV matrix
metadata_ak <- metadata_ak[rownames(metadata_ak) %in% rownames(asv_matrix), ]

# Check and align sample names after filtering
asv_samples <- rownames(asv_matrix)
metadata_samples <- rownames(metadata_ak)

# Verify alignment of sample names
if (!all(asv_samples %in% metadata_samples)) {
  stop("Sample names in ASV matrix do not match metadata after filtering.")
}

# Subset metadata to only include samples present in ASV matrix
metadata_ak <- metadata_ak[asv_samples, ]

# Compute distance matrix
distance_matrix <- vegdist(asv_matrix, method = "bray")

# Prepare grouping variable (assuming 'Site' column is to be used for grouping)
group <- metadata_ak$Site  # Replace `Site` with your actual grouping variable

# Ensure the grouping variable is a factor and matches the number of samples
if (length(group) != nrow(asv_matrix)) {
  stop("The length of the grouping variable does not match the number of samples in the ASV matrix.")
}

group <- factor(group)  # Convert to factor if not already

# Run PERMDISP
permdisp_result <- betadisper(distance_matrix, group)

# View results
summary(permdisp_result)  # Summary of PERMDISP results

# Plot distances to centroids
plot(permdisp_result$distances, col = group, pch = 19, 
     xlab = "Sample", ylab = "Distance to Centroid", 
     main = "Distances to Group Centroids-Washington",
     xlim = c(1, length(permdisp_result$distances)))  # Adjust x-axis limits if needed

# Adjust plot margins to make space for legend
par(mar = c(5, 4, 4, 9) + 0.1)

# Add legend outside the plot (right side)
legend("topright", legend = levels(group), col = 1:length(levels(group)), pch = 19, 
       inset = c(-0.2, 0),  # Adjust `inset` to position the legend further right
       xpd = TRUE)  # Allow plotting outside the plot region

# Perform PERMDISP test (ANOVA-like)
anova_result <- anova(permdisp_result)
print(anova_result)



#RESULTS
#AK$Site p= 1.71e-9
#BB$Site p = 0.1752
#BC$Site p = 0.8846
#OR$Site p = 0.04039
#SD$Site p = 0.2007
#WA$Site p = 0.333



################################ASVs (this one works)######################################
# Load necessary libraries
library(phyloseq)
library(vegan)

# Check the structure of the phyloseq object
print(ps_transposed)  # Check the overall structure
otu_table_transposed <- otu_table(ps_transposed)
print(otu_table_transposed)  # Check the OTU table structure

#  Check if the OTU table is numeric and has data
if (!is.numeric(otu_table_transposed)) {
  stop("OTU table is not numeric; please check your data.")
}

# Extract the grouping variable from sample data
grouping_variable <- sample_data(ps_transposed)$State


#  Create a distance matrix
distance_matrix <- vegdist(otu_table_transposed, method = "bray")

# Run betadisper
dispersion_result <- betadisper(distance_matrix, grouping_variable)

# Print the results
print(dispersion_result)

# To perform permutation tests, you can use the anova function
anova_result <- anova(dispersion_result)

# Print the anova results
print(anova_result)

# Plot the results (optional)
plot(dispersion_result)

