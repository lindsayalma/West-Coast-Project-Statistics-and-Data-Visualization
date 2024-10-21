library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
packageVersion("phyloseq")
packageVersion("vegan")
citation("vegan")

setwd("C:/Users/Lindsay Alma/Dropbox/Eelgrass/West Coast TagSeq primers/FullDataset")

# Load the phyloseq object
ps <- readRDS("ps_18S_lulu.RDS")

# Filter out specific persons
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

# Subset samples by Lesion Type
ps <- subset_samples(ps, Lesion_Type == "Crescent" | Lesion_Type == "Not crescent")

# Read the new ASV matrix from CSV
new_asv_data <- read.csv("otu_normalized_base_72624.csv", row.names = 1)

# Convert to matrix
new_asv_matrix <- as.matrix(new_asv_data)

# Update the otu_table with the new matrix
otu_table(ps) <- otu_table(new_asv_matrix, taxa_are_rows = TRUE)

# Extract the metadata from the phyloseq object
meta <- sample_data(ps)

# Extract ASVs
asv_matrix <- otu_table(ps)

# Extract the normalized percentage values for ASVs 1, 2, 4, 5, and 6
asv_18S_1 <- as.vector(asv_matrix["ASV_18S_1", ])  # Ensure asv_18S_1 is a vector
asv_18S_2 <- as.vector(asv_matrix["ASV_18S_2", ])  # Ensure asv_18S_2 is a vector
asv_18S_4 <- as.vector(asv_matrix["ASV_18S_4", ])  # Ensure asv_18S_4 is a vector
asv_18S_5 <- as.vector(asv_matrix["ASV_18S_5", ])  # Ensure asv_18S_5 is a vector
asv_18S_6 <- as.vector(asv_matrix["ASV_18S_6", ])  # Ensure asv_18S_6 is a vector

# Check if lengths match for all ASVs and metadata
if (length(asv_18S_1) == nrow(meta) && length(asv_18S_2) == nrow(meta) && 
    length(asv_18S_4) == nrow(meta) && length(asv_18S_5) == nrow(meta) && 
    length(asv_18S_6) == nrow(meta)) {
  
  # Create a new dominant column that considers ASVs 1, 2, 4, 5, and 6
  meta$dominant <- ifelse(asv_18S_1 > 0.10, "ASV1", 
                          ifelse(asv_18S_2 > 0.10, "ASV2",
                                 ifelse(asv_18S_4 > 0.10, "ASV4", 
                                        ifelse(asv_18S_5 > 0.10, "ASV5", 
                                               ifelse(asv_18S_6 > 0.10, "ASV6", NA)))))
  
  # Replace NA values with "coinfection"
  meta$dominant[is.na(meta$dominant)] <- "coinfection"
  
} else {
  stop("Length of ASVs does not match number of samples in metadata.")
}

# View the updated metadata to check the 'dominant' column
head(meta)

##########################################################################################

# Update the sample_data with the new metadata
sample_data(ps) <- sample_data(meta)

# Check for zero-sum samples
sample_sums <- sample_sums(ps)
zero_sum_samples <- names(sample_sums[sample_sums == 0])

# Remove zero-sum samples
ps_filtered <- prune_samples(sample_sums != 0, ps)

# Re-run PCoA
erie_pcoa <- ordinate(
  physeq = ps_filtered, 
  method = "PCoA", 
  distance = "bray"
)

# Modify the 'dominant' variable to ensure specific ASVs are on top
sample_data(ps_filtered)$dominant <- factor(sample_data(ps_filtered)$dominant, 
                                            levels = c("ASV6", "ASV5", "ASV4", "ASV2", "ASV1", "coinfection"))

# Now plot the ordination with the added ASVs
plot_ordination(
  physeq = ps_filtered,
  ordination = erie_pcoa,
  title = "PCoA"
) + 
  geom_point(aes(color = dominant), alpha = 0.9, size = 3, stroke = 0) +
  scale_color_manual(values = c(
    "coinfection" = "darkorchid4",
    "ASV1" = "cyan3",
    "ASV2" = "darkorange2",
    "ASV4" = "green3",
    "ASV5" = "blue",
    "ASV6" = "hotpink"
  ))




##############ASV Barplot horizontal###############################################################################3


# Load the phyloseq object
ps <- readRDS("ps_18S_lulu.RDS")

# Filter out specific persons
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

# Subset samples by Lesion Type
ps <- subset_samples(ps, Lesion_Type == "Crescent" | Lesion_Type == "Not crescent")

# Read the new ASV matrix from CSV
new_asv_data <- read.csv("otu_normalized_base_72624.csv", row.names = 1)

# Convert to matrix
new_asv_matrix <- as.matrix(new_asv_data)

# Update the otu_table with the new matrix
otu_table(ps) <- otu_table(new_asv_matrix, taxa_are_rows = TRUE)

# Extract the metadata from the phyloseq object
meta <- sample_data(ps)

# Extract ASVs
asv_matrix <- otu_table(ps)

# Extract the normalized percentage values for ASVs 1, 2, 4, 5, and 6
asv_18S_1 <- as.vector(asv_matrix["ASV_18S_1", ])  # Ensure asv_18S_1 is a vector
asv_18S_2 <- as.vector(asv_matrix["ASV_18S_2", ])  # Ensure asv_18S_2 is a vector
asv_18S_4 <- as.vector(asv_matrix["ASV_18S_4", ])  # Ensure asv_18S_4 is a vector
asv_18S_5 <- as.vector(asv_matrix["ASV_18S_5", ])  # Ensure asv_18S_5 is a vector
asv_18S_6 <- as.vector(asv_matrix["ASV_18S_6", ])  # Ensure asv_18S_6 is a vector

# Check if lengths match for all ASVs and metadata
if (length(asv_18S_1) == nrow(meta) && length(asv_18S_2) == nrow(meta) && 
    length(asv_18S_4) == nrow(meta) && length(asv_18S_5) == nrow(meta) && 
    length(asv_18S_6) == nrow(meta)) {
  
  # Create a new dominant column that considers ASVs 1, 2, 4, 5, and 6
  meta$dominant <- ifelse(asv_18S_1 > 0.10, "ASV1", 
                          ifelse(asv_18S_2 > 0.10, "ASV2",
                                 ifelse(asv_18S_4 > 0.10, "ASV4", 
                                        ifelse(asv_18S_5 > 0.10, "ASV5", 
                                               ifelse(asv_18S_6 > 0.10, "ASV6", NA)))))
  
  # Replace NA values with "coinfection"
  meta$dominant[is.na(meta$dominant)] <- "coinfection"
  
} else {
  stop("Length of ASVs does not match number of samples in metadata.")
}

# View the updated metadata to check the 'dominant' column
head(meta)

##########################################################################################

# Update the sample_data with the new metadata
sample_data(ps) <- sample_data(meta)

# Check for zero-sum samples
sample_sums <- sample_sums(ps)
zero_sum_samples <- names(sample_sums[sample_sums == 0])

# Remove zero-sum samples
ps_filtered <- prune_samples(sample_sums != 0, ps)

# Re-run PCoA
erie_pcoa <- ordinate(
  physeq = ps_filtered, 
  method = "PCoA", 
  distance = "bray"
)

# Modify the 'dominant' variable to ensure specific ASVs are on top
sample_data(ps_filtered)$dominant <- factor(sample_data(ps_filtered)$dominant, 
                                            levels = c("ASV6", "ASV5", "ASV4", "ASV2", "ASV1", "coinfection"))

# Now plot the ordination with the added ASVs
plot_ordination(
  physeq = ps_filtered,
  ordination = erie_pcoa,
  title = "PCoA"
) + 
  geom_point(aes(color = dominant), alpha = 0.9, size = 3, stroke = 0) +
  scale_color_manual(values = c(
    "coinfection" = "darkorchid4",
    "ASV1" = "cyan3",
    "ASV2" = "darkorange2",
    "ASV4" = "green3",
    "ASV5" = "blue",
    "ASV6" = "hotpink"
  ))



##############ASV Barplot horizontal###############################################################################3


# Load the phyloseq object
ps <- readRDS("ps_18S_lulu.RDS")

# Filter out specific persons
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

# Subset samples by Lesion Type
ps <- subset_samples(ps, Lesion_Type == "Crescent" | Lesion_Type == "Not crescent")

# Read the new ASV matrix from CSV
new_asv_data <- read.csv("otu_normalized_base_72624.csv", row.names = 1)

# Convert to matrix
new_asv_matrix <- as.matrix(new_asv_data)

# Update the otu_table with the new matrix
otu_table(ps) <- otu_table(new_asv_matrix, taxa_are_rows = TRUE)

# Extract the metadata from the phyloseq object
meta <- sample_data(ps)

# Extract ASVs
asv_matrix <- otu_table(ps)

# Extract the normalized percentage values for ASVs 1, 2, 4, 5, and 6
asv_18S_1 <- as.vector(asv_matrix["ASV_18S_1", ])  # Ensure asv_18S_1 is a vector
asv_18S_2 <- as.vector(asv_matrix["ASV_18S_2", ])  # Ensure asv_18S_2 is a vector
asv_18S_4 <- as.vector(asv_matrix["ASV_18S_4", ])  # Ensure asv_18S_4 is a vector
asv_18S_5 <- as.vector(asv_matrix["ASV_18S_5", ])  # Ensure asv_18S_5 is a vector
asv_18S_6 <- as.vector(asv_matrix["ASV_18S_6", ])  # Ensure asv_18S_6 is a vector

# Check if lengths match for all ASVs and metadata
if (length(asv_18S_1) == nrow(meta) && length(asv_18S_2) == nrow(meta) && 
    length(asv_18S_4) == nrow(meta) && length(asv_18S_5) == nrow(meta) && 
    length(asv_18S_6) == nrow(meta)) {
  
  # Create a new dominant column that considers ASVs 1, 2, 4, 5, and 6
  meta$dominant <- ifelse(asv_18S_1 > 0.10, "ASV1", 
                          ifelse(asv_18S_2 > 0.10, "ASV2",
                                 ifelse(asv_18S_4 > 0.10, "ASV4", 
                                        ifelse(asv_18S_5 > 0.10, "ASV5", 
                                               ifelse(asv_18S_6 > 0.10, "ASV6", NA)))))
  
  # Replace NA values with "coinfection"
  meta$dominant[is.na(meta$dominant)] <- "coinfection"
  
} else {
  stop("Length of ASVs does not match number of samples in metadata.")
}

# View the updated metadata to check the 'dominant' column
head(meta)

##########################################################################################

# Update the sample_data with the new metadata
sample_data(ps) <- sample_data(meta)

# Check for zero-sum samples
sample_sums <- sample_sums(ps)
zero_sum_samples <- names(sample_sums[sample_sums == 0])

# Remove zero-sum samples
ps_filtered <- prune_samples(sample_sums != 0, ps)

# Re-run PCoA
erie_pcoa <- ordinate(
  physeq = ps_filtered, 
  method = "PCoA", 
  distance = "bray"
)

# Modify the 'dominant' variable to ensure specific ASVs are on top
sample_data(ps_filtered)$dominant <- factor(sample_data(ps_filtered)$dominant, 
                                            levels = c("ASV6", "ASV5", "ASV4", "ASV2", "ASV1", "coinfection"))

# Now plot the ordination with the added ASVs
plot_ordination(
  physeq = ps_filtered,
  ordination = erie_pcoa,
  title = "PCoA"
) + 
  geom_point(aes(color = dominant), alpha = 0.9, size = 3, stroke = 0) +
  scale_color_manual(values = c(
    "coinfection" = "darkorchid4",
    "ASV1" = "cyan3",
    "ASV2" = "darkorange2",
    "ASV4" = "green3",
    "ASV5" = "blue",
    "ASV6" = "hotpink"
  ))



##############ASV Barplot horizontal###############################################################################3

# Ensure 'dominant' contains ASV names as factors
meta$dominant <- factor(meta$dominant, levels = c("ASV1", "ASV2", "ASV4", "ASV5", "ASV6", "coinfection"))

# Plot ASVs 1, 2, 4, 5, 6 with specified colors for each state
ggplot(meta, aes(y = dominant, fill = State)) +
  geom_bar() +
  coord_flip() +
  theme_minimal() +
  labs(x = "Count", y = "ASV", fill = "State", title = "Distribution of Dominant ASVs by State") +
  scale_fill_manual(values = c(
    "AK" = "orange",
    "BB" = "gold",
    "BC" = "darkgreen",
    "OR" = "blue",
    "SD" = "lightblue",
    "WA" = "purple"
  ))


# Plot ASVs 1, 2, 4, 5, 6 with specified colors for each state
ggplot(meta, aes(y = dominant, fill = Lesion_Type)) +
  geom_bar() +
  coord_flip() +
  theme_minimal() +
  labs(x = "Count", y = "ASV", fill = "Lesion_Type", title = "Distribution of Dominant ASVs by State") +
  scale_fill_manual(values = c(
    "Crescent" = "darkorange3",
    "Not crescent" = "forestgreen"
      ))




############### STATE PCOA##################################################################################

# Update the sample_data with the new metadata
sample_data(ps) <- sample_data(meta)

# Check for zero-sum samples
sample_sums <- sample_sums(ps)
zero_sum_samples <- names(sample_sums[sample_sums == 0])

# Remove zero-sum samples
ps_filtered <- prune_samples(sample_sums != 0, ps)

# Re-run PCoA
erie_pcoa <- ordinate(
  physeq = ps_filtered, 
  method = "PCoA", 
  distance = "bray"
)

# Modify the 'State' variable to ensure specific states are on top
sample_data(ps_filtered)$State <- factor(sample_data(ps_filtered)$State, 
                                         levels = c("AK", "BB", "BC", "OR", "SD", "WA"))

# Now plot the ordination with the 'State' variable
plot_ordination(
  physeq = ps_filtered,
  ordination = erie_pcoa,
  title = "PCoA by State"
) + 
  geom_point(aes(color = State), alpha = 0.9, size = 3, stroke = 0) +
  scale_color_manual(values = c(
    "AK" = "orange",
    "BB" = "gold",
    "BC" = "darkgreen",
    "OR" = "blue",
    "SD" = "lightblue",
    "WA" = "purple"
  ))




######################################### PCOA L vs G################################3
# Update the sample_data with the new metadata
sample_data(ps) <- sample_data(meta)

# Check for zero-sum samples
sample_sums <- sample_sums(ps)
zero_sum_samples <- names(sample_sums[sample_sums == 0])

# Remove zero-sum samples
ps_filtered <- prune_samples(sample_sums != 0, ps)

# Re-run PCoA
erie_pcoa <- ordinate(
  physeq = ps_filtered, 
  method = "PCoA", 
  distance = "bray"
)

# Modify the 'Sample_Type' variable to ensure specific types are on top
sample_data(ps_filtered)$Sample_Type <- factor(sample_data(ps_filtered)$Sample_Type, 
                                               levels = c("L", "G"))

# Now plot the ordination with the 'Sample_Type' variable
plot_ordination(
  physeq = ps_filtered,
  ordination = erie_pcoa,
  title = "PCoA by Sample Type"
) + 
  geom_point(aes(color = Sample_Type), alpha = 0.9, size = 3, stroke = 0) +
  scale_color_manual(values = c(
    "L" = "orange",
    "G" = "green"
  ))


####################### PCOA Crescent vs not ###################################

# Update the sample_data with the new metadata
sample_data(ps) <- sample_data(meta)

# Check for zero-sum samples
sample_sums <- sample_sums(ps)
zero_sum_samples <- names(sample_sums[sample_sums == 0])

# Remove zero-sum samples
ps_filtered <- prune_samples(sample_sums != 0, ps)

# Re-run PCoA
erie_pcoa <- ordinate(
  physeq = ps_filtered, 
  method = "PCoA", 
  distance = "bray"
)

# Modify the 'Lesion_Type' variable to ensure specific types are on top
sample_data(ps_filtered)$Lesion_Type <- factor(sample_data(ps_filtered)$Lesion_Type, 
                                               levels = c("Crescent", "Not crescent"))

# Now plot the ordination with the 'Lesion_Type' variable
plot_ordination(
  physeq = ps_filtered,
  ordination = erie_pcoa,
  title = "PCoA by Lesion Type"
) + 
  geom_point(aes(color = Lesion_Type), alpha = 0.9, size = 3, stroke = 0) +
  scale_color_manual(values = c(
    "Crescent" = "darkorange3",
    "Not crescent" = "darkgreen"
  ))


##################################### dominant ASV1 + crescent#########################
dominant_samples <- meta[meta$dominant == "ASV1" & meta$Lesion_Type == "Crescent", ]
dominant_samples
