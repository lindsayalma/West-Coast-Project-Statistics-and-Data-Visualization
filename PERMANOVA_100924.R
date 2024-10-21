library(vegan)
library(phyloseq)
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

# Remove only completely empty rows and columns (i.e., those that are zero across all samples)
otu_table_transposed_cleaned <- otu_table_transposed[rowSums(otu_table_transposed) != 0, ]


# Check for remaining missing values and handle them as necessary (optional, depending on analysis)
otu_table_transposed_cleaned[is.na(otu_table_transposed_cleaned)] <- 0


# Ensure that the metadata aligns with the OTU table samples
metadata_cleaned <- as(sample_data(ps), "data.frame")

# Remove samples in metadata that are not in the OTU table
metadata_cleaned <- metadata_cleaned[rownames(metadata_cleaned) %in% rownames(otu_table_transposed_cleaned), ]

# Also subset OTU table to match metadata
otu_table_transposed_cleaned <- otu_table_transposed_cleaned[rownames(metadata_cleaned), ]

# Check if row names match between OTU table and metadata
if (!identical(rownames(metadata_cleaned), rownames(otu_table_transposed_cleaned))) {
  stop("Mismatch between OTU table and metadata row names")
}

# Ensure 'Lesion_Type' is a factor
metadata_cleaned$Lesion_Type <- as.factor(metadata_cleaned$Lesion_Type)


# Recalculate the distance matrix using Bray-Curtis, after ensuring no empty rows
dist_matrix_cleaned <- vegdist(otu_table_transposed_cleaned, method = "bray")

# Check for missing values in the distance matrix
if (any(is.na(as.matrix(dist_matrix_cleaned)))) {
  stop("Missing values found in the distance matrix.")
}


# Run PERMANOVA on Lesion_Type variable using the distance matrix
permanova_result <- adonis2(
  dist_matrix_cleaned ~ Lesion_Type,         # Formula
  data = metadata_cleaned,             # Data frame containing the variables
  permutations = 10000,                # Number of permutations
  method = "bray"                      # Ensure this matches your distance method
)

# Display the results
print(permanova_result)




######################################################################
#########INDIVIDUAL STATES##################################3

ps_state <- subset_samples(ps, State == "AK")
#ps_state <- subset_samples(ps, State == "BC")
#ps_state <- subset_samples(ps, State == "WA")
#ps_state <- subset_samples(ps, State == "OR")
#ps_state <- subset_samples(ps, State == "BB")
#ps_state <- subset_samples(ps, State == "SD")
sample_data(ps_state)
table(sample_data(ps_state)$State)
nsamples(ps_state)
ntaxa(ps_state)

# Extract OTU table and transpose
otu_table_matrix <- as(otu_table(ps_state), "matrix")
otu_table_transposed <- t(otu_table_matrix)

# Remove empty rows and columns
otu_table_transposed_cleaned <- otu_table_transposed[rowSums(otu_table_transposed) > 0, ]
otu_table_transposed_cleaned <- otu_table_transposed_cleaned[, colSums(otu_table_transposed_cleaned) > 0]

# Handle missing values by replacing them with zeros
otu_table_transposed_cleaned[is.na(otu_table_transposed_cleaned)] <- 0

# Recalculate the distance matrix
dist_matrix_cleaned <- vegdist(otu_table_transposed_cleaned, method = "bray")
as.matrix(dist_matrix_cleaned)[1:10, 1:10]  # Display first 5 rows and columns
# Check class and dimensions of the distance matrix
cat("Distance Matrix Class: ", class(dist_matrix_cleaned), "\n")
cat("Distance Matrix Dimensions: ", attr(dist_matrix_cleaned, "Size"), "\n")

# Ensure metadata is aligned
metadata_cleaned <- sample_data(ps_state)[row.names(sample_data(ps_state)) %in% row.names(otu_table_transposed_cleaned), ]
metadata_cleaned <- as.data.frame(sample_data(ps_state))
# Subset metadata to match the OTU table
metadata_cleaned <- metadata_cleaned[rownames(otu_table_transposed_cleaned), ]
# Check if the row names match the samples in the cleaned OTU table
all(rownames(metadata_cleaned) %in% rownames(otu_table_transposed_cleaned))  # Should be TRUE

# Check the alignment
head(metadata_cleaned)
# Verify the alignment of metadata with OTU table
identical(rownames(metadata_cleaned), rownames(otu_table_transposed_cleaned))  # Should be TRUE

# Check dimensions to confirm alignment
cat("OTU Table Samples: ", nrow(otu_table_transposed_cleaned), "\n")
cat("Metadata Samples: ", nrow(metadata_cleaned), "\n")

# Check column names of the metadata
colnames(metadata_cleaned)

# Check levels of the 'State' variable
table(metadata_cleaned$State)

# Specify the formula explicitly


# Convert 'State' to a factor if it's not already
metadata_cleaned$State <- as.factor(metadata_cleaned$State)

# Check the levels of the 'State' factor
levels(metadata_cleaned$State)
# Extract sample data and ensure it's a data frame
metadata_cleaned <- as.data.frame(sample_data(ps_state))
# Check first few rows of the metadata
head(metadata_cleaned)
metadata_cleaned <- metadata_cleaned[rownames(otu_table_transposed_cleaned), ]
metadata_cleaned$State <- as.factor(metadata_cleaned$State)
# Check dimensions and contents
cat("Distance Matrix Dimensions:", attr(dist_matrix_cleaned, "Size"), "\n")
cat("Metadata Dimensions:", nrow(metadata_cleaned), "\n")
head(metadata_cleaned)
# Convert sample data to a data frame
metadata_cleaned <- as(sample_data(ps_state), "data.frame")

# Align metadata with OTU table samples
metadata_cleaned <- metadata_cleaned[rownames(otu_table_transposed_cleaned), ]

# Verify that metadata and OTU table samples match
identical(rownames(metadata_cleaned), rownames(otu_table_transposed_cleaned))  # Should be TRUE

# Check the structure of the metadata
str(metadata_cleaned)
head(metadata_cleaned)



# Run PERMANOVA

# Example code with adonis2
permanova_result <- adonis2(
  dist_matrix_cleaned ~ Site,         # Formula
  data = metadata_cleaned,             # Data frame containing the variables
  permutations = 10000,                  # Number of permutations
  method = "bray"                      # Ensure this matches your distance method
)

# Display the results
print(permanova_result)

#
#AKp=9.999e-5
#BC=0.0021
#WA= 0.0477
#OR= 9.999e-5
#BB=0.004
#SD= 2e-4


########################################
############L vs G#################
# Load or create a phyloseq object
ps <- readRDS("ps_18S_lulu.RDS")

# Extract OTU table and transpose
otu_table_matrix <- as(otu_table(ps), "matrix")
otu_table_transposed <- t(otu_table_matrix)

# Remove empty rows and columns
otu_table_transposed_cleaned <- otu_table_transposed[rowSums(otu_table_transposed) > 0, ]
otu_table_transposed_cleaned <- otu_table_transposed_cleaned[, colSums(otu_table_transposed_cleaned) > 0]

# Handle missing values by replacing them with zeros
otu_table_transposed_cleaned[is.na(otu_table_transposed_cleaned)] <- 0

# Recalculate the distance matrix
dist_matrix_cleaned <- vegdist(otu_table_transposed_cleaned, method = "bray")
as.matrix(dist_matrix_cleaned)[1:10, 1:10]  # Display first 5 rows and columns
# Check class and dimensions of the distance matrix
cat("Distance Matrix Class: ", class(dist_matrix_cleaned), "\n")
cat("Distance Matrix Dimensions: ", attr(dist_matrix_cleaned, "Size"), "\n")

# Ensure metadata is aligned
metadata_cleaned <- sample_data(ps)[row.names(sample_data(ps)) %in% row.names(otu_table_transposed_cleaned), ]
metadata_cleaned <- as.data.frame(sample_data(ps))
# Subset metadata to match the OTU table
metadata_cleaned <- metadata_cleaned[rownames(otu_table_transposed_cleaned), ]

# Check if the row names match the samples in the cleaned OTU table
all(rownames(metadata_cleaned) %in% rownames(otu_table_transposed_cleaned))  # Should be TRUE

# Check the alignment
head(metadata_cleaned)
# Verify the alignment of metadata with OTU table
identical(rownames(metadata_cleaned), rownames(otu_table_transposed_cleaned))  # Should be TRUE

# Check dimensions to confirm alignment
cat("OTU Table Samples: ", nrow(otu_table_transposed_cleaned), "\n")
cat("Metadata Samples: ", nrow(metadata_cleaned), "\n")

# Check column names of the metadata
colnames(metadata_cleaned)

# Check levels of the 'State' variable
table(metadata_cleaned$Sample_Type)

# Specify the formula explicitly
formula <- as.formula("dist_matrix_cleaned ~ Sample_Type")

# Convert 'State' to a factor if it's not already
metadata_cleaned$Sample_type <- as.factor(metadata_cleaned$Sample_Type)

# Check the levels of the 'State' factor
levels(metadata_cleaned$Sample_Type)
# Extract sample data and ensure it's a data frame
metadata_cleaned <- as.data.frame(sample_data(ps))
# Check first few rows of the metadata
head(metadata_cleaned)
metadata_cleaned <- metadata_cleaned[rownames(otu_table_transposed_cleaned), ]
metadata_cleaned$Sample_Type <- as.factor(metadata_cleaned$Sample_Type)
# Check dimensions and contents
cat("Distance Matrix Dimensions:", attr(dist_matrix_cleaned, "Size"), "\n")
cat("Metadata Dimensions:", nrow(metadata_cleaned), "\n")
head(metadata_cleaned)
# Convert sample data to a data frame
metadata_cleaned <- as(sample_data(ps), "data.frame")

# Align metadata with OTU table samples
metadata_cleaned <- metadata_cleaned[rownames(otu_table_transposed_cleaned), ]

# Verify that metadata and OTU table samples match
identical(rownames(metadata_cleaned), rownames(otu_table_transposed_cleaned))  # Should be TRUE

# Check the structure of the metadata
str(metadata_cleaned)
head(metadata_cleaned)



# Run PERMANOVA
permanova_result <- adonis2(formula, data = metadata_cleaned)
# Example code with adonis2
permanova_result <- adonis2(
  dist_matrix_cleaned ~ Sample_Type,         # Formula
  data = metadata_cleaned,             # Data frame containing the variables
  permutations = 10000,                  # Number of permutations
  method = "bray"                      # Ensure this matches your distance method
)

# Display the results
print(permanova_result)

# Optional: Visualize with MDS
mds_result <- cmdscale(dist_matrix_cleaned, k = 2)
plot(mds_result, main = "MDS of Bray-Curtis Dissimilarity", xlab = "MDS1", ylab = "MDS2")


#L vs G = 9.999e-5

######################## ASVs#######################################



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
new_asv_matrix <- as.matrix(meta)

# Update the otu_table with the new matrix
otu_table(ps) <- otu_table(new_asv_matrix, taxa_are_rows = TRUE)

# Extract and transpose the OTU table
otu_table_matrix <- as(otu_table(ps), "matrix")
otu_table_transposed <- t(otu_table_matrix)

# Remove only completely empty rows and columns (i.e., those that are zero across all samples)
otu_table_transposed_cleaned <- otu_table_transposed[rowSums(otu_table_transposed) != 0, ]


# Check for remaining missing values and handle them as necessary (optional, depending on analysis)
otu_table_transposed_cleaned[is.na(otu_table_transposed_cleaned)] <- 0


# Ensure that the metadata aligns with the OTU table samples
metadata_cleaned <- as(sample_data(ps), "data.frame")

# Remove samples in metadata that are not in the OTU table
metadata_cleaned <- metadata_cleaned[rownames(metadata_cleaned) %in% rownames(otu_table_transposed_cleaned), ]

# Also subset OTU table to match metadata
otu_table_transposed_cleaned <- otu_table_transposed_cleaned[rownames(metadata_cleaned), ]

# Check if row names match between OTU table and metadata
if (!identical(rownames(metadata_cleaned), rownames(otu_table_transposed_cleaned))) {
  stop("Mismatch between OTU table and metadata row names")
}

# Ensure 'Lesion_Type' is a factor
metadata_cleaned$Lesion_Type <- as.factor(metadata_cleaned$Lesion_Type)


# Recalculate the distance matrix using Bray-Curtis, after ensuring no empty rows
dist_matrix_cleaned <- vegdist(otu_table_transposed_cleaned, method = "bray")

# Check for missing values in the distance matrix
if (any(is.na(as.matrix(dist_matrix_cleaned)))) {
  stop("Missing values found in the distance matrix.")
}


# Run PERMANOVA on Lesion_Type variable using the distance matrix
permanova_result <- adonis2(
  dist_matrix_cleaned ~ Lesion_Type,         # Formula
  data = metadata_cleaned,             # Data frame containing the variables
  permutations = 10000,                # Number of permutations
  method = "bray"                      # Ensure this matches your distance method
)

# Display the results
print(permanova_result)

#p=9.99e-5