library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(vegan)


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

# Convert to matrix
new_asv_matrix <- as.matrix(new_asv_data)

# Assuming your phyloseq object is named `ps`
# Update the otu_table with the new matrix
ps_new <- ps
otu_table(ps_new) <- otu_table(new_asv_matrix, taxa_are_rows = TRUE)

meta<- sample_data(ps)

# Create the new columns by splitting the 'Sample' column
meta$Sample_match <- sub("^(.*?)-(.*?)-.*$", "\\1-\\2", meta$Sample)
meta$match <- sub("^(.*?)-(.*?)-(.*)$", "\\3", meta$Sample)

# Reorder the columns to place 'Sample_match' and 'match' next to each other
meta <- meta[, c(setdiff(names(meta), c("Sample_match", "match")), "Sample_match", "match")]

# Extract the numbers after the first character in the 'match' column
meta$new_column <- sub("^.", "", meta$match)

# Reorder the columns to place the new column next to 'Sample_match' and 'match'
meta <- meta[, c(setdiff(names(meta), "new_column"), "Sample_match", "match", "new_column")]

# Create a binary match column
meta$binary_match <- ave(seq_along(meta$Sample_match), meta$Sample_match, meta$new_column, FUN = function(x) ifelse(length(x) > 1, 1, 0))

# Keep only rows where binary_match is 1
meta <- meta[meta$binary_match == 1, ]

# Sort by new_column and then by Sample_match
meta <- meta[order(meta$new_column, meta$Sample_match), ]

# Save the meta dataframe as a CSV file
write.csv(meta, file = "meta.csv", row.names = FALSE)
metamatch <- read.csv("meta.csv", header = TRUE)

# Ensure the Sample columns exist
if ("Sample" %in% colnames(meta) && "Sample" %in% colnames(metamatch)) {
  
  # Match the Sample column in meta with the Sample column in metamatch
  # Create a new column in metamatch to store the corresponding row names from meta
  metamatch$Sequence <- NA  # Initialize the new column with the name 'Sequence'
  
  # Use match to find the indices of the matched Samples
  match_indices <- match(metamatch$Sample, meta$Sample)
  
  # Transfer over the row names from meta
  metamatch$Sequence[!is.na(match_indices)] <- rownames(meta)[match_indices[!is.na(match_indices)]]
} else {
  message("Sample column not found in one or both dataframes.")
}

# Ensure that the Sequence column exists in metamatch
if ("Sequence" %in% colnames(metamatch)) {
  
  # Get the unique values from the Sequence column
  matching_sequences <- unique(metamatch$Sequence)
  
  # Subset new_asv_data to keep only the columns that match the Sequence values
  new_asv_data <- new_asv_data[, colnames(new_asv_data) %in% matching_sequences, drop = FALSE]
} else {
  message("Sequence column not found in metamatch.")
}

new_asv_matrix <- as.matrix(new_asv_matrix)
otu_table <- otu_table(new_asv_matrix, taxa_are_rows = TRUE)
ps <- phyloseq(otu_table, metamatch)

# Check and handle missing values in the OTU table
otu_table(ps)[is.na(otu_table(ps))] <- 0  # Replace NAs with 0

# Remove samples with all zeros
ps <- prune_samples(colSums(otu_table(ps)) > 0, ps)

# Optionally, verify the number of remaining samples
cat("Number of remaining samples:", nsamples(ps), "\n")

# Calculate Distance Matrix
distance_matrix <- distance(ps, method = "bray")

# Run PCoA
pcoa_result <- cmdscale(distance_matrix, k = 2, eig = TRUE)

# Extract Coordinates
pcoa_df <- as.data.frame(pcoa_result$points)
pcoa_df$SampleID <- rownames(pcoa_df)


# Change the column name from 'Sequence' to 'SampleID'
colnames(metamatch)[colnames(metamatch) == "Sequence"] <- "SampleID"

# Add Sample Metadata
sample_data_df <- data.frame(metamatch)
pcoa_df <- merge(pcoa_df, metamatch, by = "SampleID")

print(colnames(pcoa_df))
sample_data_df <- data.frame(metamatch)
print(head(sample_data_df))

# Plot the PCoA Results
ggplot(pcoa_df, aes(x = V1, y = V2, color = State)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCoA of Phyloseq Data", x = "PCoA 1", y = "PCoA 2") +
  theme(legend.position = "right")


pcoa_df <- pcoa_df %>%
  left_join(sample_data_df[, c("SampleID", "match_number")], by = "SampleID")
print(head(pcoa_df))

## Keep only samples that have at least one other match_number in common
matched_samples <- pcoa_df %>%
  group_by(match_number.y) %>%
  filter(n() > 1) %>%  # Keep groups with more than 1 sample
  ungroup()

# Create the ggplot with points and lines connecting matches
ggplot(matched_samples, aes(x = V1, y = V2, color = State)) + 
  geom_point(size = 3) +  # Add points
  geom_line(aes(group = match_number.y), size = 1, alpha = 0.5) +  # Draw lines connecting matches
  theme_minimal() + 
  labs(title = "PCoA of Phyloseq Data with Match Lines", x = "PCoA 1", y = "PCoA 2") + 
  theme(legend.position = "right")

ggplot(matched_samples, aes(x = V1, y = V2, color = Sample_Type)) + 
  geom_point(size = 3) +  # Add points
  geom_line(aes(group = match_number.y), size = 1, color = "black", alpha = 0.5) +  # Draw black lines connecting matches
  theme_minimal() + 
  labs(title = "PCoA of Phyloseq Data with Match Lines", x = "PCoA 1", y = "PCoA 2") + 
  theme(legend.position = "right")





