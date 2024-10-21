# Load or create a phyloseq object
ps <- readRDS("ps_18S_lulu.RDS")
# Assume `ps` is your phyloseq object
# Extract OTU table as a matrix
otu_matrix <- as(otu_table(ps), "matrix")


# Convert to a data frame if necessary
otu_df <- as.data.frame(otu_matrix)

# Check the first few rows
head(otu_df)

# Sum each column of the OTU table
column_sums <- colSums(otu_df, na.rm = TRUE)

# Convert to a data frame
sums_df_base <- data.frame(OTU = names(column_sums), Sum = column_sums)

# Print the resulting data frame
print(sums_df_base)
View(sums_df_base)

# Extract the OTU table as a matrix
otu_matrix <- as(otu_table(ps), "matrix")

# Convert to a data frame
otu_df <- as.data.frame(otu_matrix)

# Display the original OTU table
print("Original OTU Table:")
View(otu_df)

# Calculate column sums
column_sums <- colSums(otu_df, na.rm = TRUE)
# Normalize each element by its column sum
otu_normalized_base <- sweep(otu_df, 2, column_sums, FUN = "/")
# Print the normalized OTU table
View(otu_normalized_base)
# Write the OTU_normalized_base to a CSV file
write.csv(otu_normalized_base, "otu_normalized_base.csv", row.names = TRUE)

