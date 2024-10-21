library(vegan)
library(phyloseq)
library(ggplot2)
library(reshape2)
library(dplyr)

ps <- readRDS("ps_18S_lulu.RDS")

ps <- subset_samples(ps, Person != "Olivia" & Person != "CB")

# Remove non-target taxa based on 'kingdom' and 'domain'
ps <- subset_taxa(ps, is.na(domain) | domain != "Bacteria")  # Remove Bacteria
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Chloroplastida")  # Remove Chloroplastida
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Animalia")  # Remove Animalia
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Alveolata")  # Remove Alveolata
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Amoebozoa")  # Remove Amoebozoa
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Holozoa")  # Remove Holozoa
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Rhodophyceae")  # Remove Rhodophyceae
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Diatomea")  # Remove Diatomea
ps <- subset_taxa(ps, is.na(kingdom) | kingdom != "Ochrophyta")  # Remove Ochrophyta

#ps <- subset_samples(ps, State == "Crescent" | State == "Not crescent")
# Check the resulting sample and taxa counts
cat("Final sample count:", nsamples(ps), "\n")
cat("Final taxa count:", ntaxa(ps), "\n")


# Extract OTU table and sample data
otu_matrix <- as(otu_table(ps), "matrix")

# Transpose the OTU matrix
otu_matrix_transposed <- t(otu_matrix)

#  Extract the sample data and convert it to a data frame
sample_data_df <- as(sample_data(ps), "data.frame")

#  Convert the data frame back to a sample_data object
new_sample_data <- sample_data(sample_data_df)

#  Update the sample data in the phyloseq object
sample_data(ps) <- new_sample_data

View(sample_data(ps))
View(sample_data(ps))


# Extract sample data
sample_data <- as(sample_data(ps), "data.frame")


# Ensure row names are properly set
if (is.null(rownames(otu_matrix_transposed))) rownames(otu_matrix_transposed) <- seq_len(nrow(otu_matrix_transposed))
if (is.null(rownames(sample_data))) rownames(sample_data) <- seq_len(nrow(sample_data))

# Ensure the OTU matrix and sample data have matching samples
common_samples <- intersect(rownames(otu_matrix_transposed), rownames(sample_data))
print(length(common_samples)) # Number of common samples

# Filter the transposed OTU matrix and sample data
otu_matrix_filtered <- otu_matrix_transposed[common_samples, , drop = FALSE]
sample_data_filtered <- sample_data[common_samples, , drop = FALSE]

# Check dimensions after filtering
print("Dimensions after filtering:")
print(dim(otu_matrix_filtered))
print(dim(sample_data_filtered))


# Calculate alpha diversity
shannon_diversity <- diversity(otu_matrix_filtered, index = "shannon")
simpson_diversity <- diversity(otu_matrix_filtered, index = "simpson")
observed_richness <- rowSums(otu_matrix_filtered > 0)

# Combine results into a data frame
alpha_diversity <- data.frame(
  Shannon = shannon_diversity,
  Simpson = simpson_diversity,
  Observed = observed_richness,
  row.names = common_samples
)

# Combine with sample data
alpha_diversity_combined <- cbind(sample_data_filtered[rownames(alpha_diversity), ], alpha_diversity)


# Convert the alpha diversity data to long format for ggplot2
alpha_diversity_long <- melt(alpha_diversity_combined, id.vars = colnames(sample_data_filtered))

# Define colors for the boxplots
box_colors <- c("Not crescent" = "darkgreen", "Crescent" = "darkorange3")  # Replace with your actual lesion type names


# Plot Shannon Diversity
ggplot(alpha_diversity_long[alpha_diversity_long$variable == "Shannon", ], aes(x = State, y = value, fill = State)) +
  geom_boxplot() +
  labs(title = "Shannon Diversity Index by State",
       x = "State",
       y = "Shannon Diversity Index") +
  scale_fill_manual(values = box_colors) +
  theme_minimal()

# Plot Simpson Diversity
ggplot(alpha_diversity_long[alpha_diversity_long$variable == "Simpson", ], aes(x = State, y = value, fill = State)) +
  geom_boxplot() +
  labs(title = "Simpson Diversity Index by State",
       x = "State",
       y = "Simpson Diversity Index") +
  scale_fill_manual(values = box_colors) +
  theme_minimal()

# Plot Observed Richness
ggplot(alpha_diversity_long[alpha_diversity_long$variable == "Observed", ], aes(x = State, y = value, fill = State)) +
  geom_boxplot() +
  labs(title = "Observed Richness by State",
       x = "State",
       y = "Observed Richness") +
  scale_fill_manual(values = box_colors) +
  theme_minimal()

# Define colors for the boxplots with specific colors for each State
box_colors <- c(    "AK" = "orange",
                    "BB" = "gold",
                    "BC" = "darkgreen",
                    "OR" = "blue",
                    "SD" = "lightblue",
                    "WA" = "purple")


# Create combined plot with specified colors
ggplot(alpha_diversity_long, aes(x = State, y = value, fill = State)) +
  geom_boxplot() +
  labs(x = "State", y = "Diversity Index Value") +
  facet_wrap(~variable, scales = "free_y", ncol = 1) +  # Facet by diversity index type
  scale_fill_manual(values = box_colors) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))  # Bold facet labels




# List of diversity indices to test
diversity_indices <- c("Shannon", "Simpson", "Observed")

# Perform Kruskal-Wallis test for all indices in one line
kruskal_test_results <- lapply(diversity_indices, function(index) {
  kruskal.test(as.formula(paste(index, "~ State")), data = alpha_diversity_combined)
})


# Print the results
names(kruskal_test_results) <- diversity_indices  # Label the results by index name
kruskal_test_results

# List of diversity metrics
diversity_metrics <- c("Shannon", "Simpson", "Observed")

# Perform pairwise Wilcoxon test for each diversity metric and store the results in a list
pairwise_wilcox_results <- lapply(diversity_metrics, function(metric) {
  pairwise.wilcox.test(alpha_diversity_combined[[metric]], alpha_diversity_combined$State, p.adjust.method = "BH")
})

# Print the results for each diversity metric
names(pairwise_wilcox_results) <- diversity_metrics
pairwise_wilcox_results

