library(phyloseq)
library(ggplot2)
library(dplyr)

ps <- readRDS("ps_18S_lulu.RDS")

# Extract the sample data
sample_data_df <- sample_data(ps)

# Convert to a data frame for editing
sample_data_df <- as.data.frame(sample_data_df)

sample_data_df$Copies.mg[is.na(sample_data_df$Copies.mg)] <- 0

# Create the "qPCR" column based on the "Copies.mg" values
sample_data_df$qPCR <- ifelse(sample_data_df$Copies.mg > 0, "positive", "negative")

# Create a new column combining "qPCR" and "Sample_Type"
sample_data_df$qPCR_Sample_Type <- paste(sample_data_df$qPCR, sample_data_df$Sample_Type, sep = "_")

# Create a new column combining "qPCR" and "Lesion_Type"
sample_data_df$qPCR_Lesion_Type <- paste(sample_data_df$qPCR, sample_data_df$Lesion_Type, sep = "_")

View(sample_data_df)

# Convert back to a phyloseq sample data object
sample_data(ps) <- sample_data(sample_data_df)

# Export sample_data_df to a CSV file in a specific directory
write.csv(sample_data_df, file = "qPCR_sample_data_df.csv", row.names = TRUE)


################################################################################
asv_data <- otu_table(ps)  # Extract OTU table (ASV matrix)
metadata <- sample_data(ps)  # Extract metadata     

# Convert the OTU table to a data frame (if necessary)
asv_data_df <- as.data.frame(asv_data)

# Write the data frame to a CSV file
write.csv(asv_data_df, file = "asv_data_raw_082624.csv", row.names = TRUE)

asv_df<- read.csv("C:/Users/Lindsay Alma/Dropbox/Eelgrass/West Coast TagSeq primers/FullDataset/asv_data_raw_top20_082624.csv")

# Transpose the data frame
transposed_asv_df <- t(asv_df)

# Convert the transposed matrix back to a data frame
transposed_asv_df <- as.data.frame(transposed_asv_df)
new_header <- transposed_asv_df[1, ]
transposed_asv_df <- transposed_asv_df[-1, ]
colnames(transposed_asv_df) <- new_header

#write.csv(transposed_asv_df, file = "asv_data_raw_transposed_082624.csv", row.names = TRUE)
asv<- read.csv("asv_data_raw_transposed_082824.csv")
View(asv)

##########################################################################################

# Create a bar plot of the qPCR column
ggplot(asv, aes(x = qPCR)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(title = "Distribution of qPCR Results",
       x = "qPCR Result",
       y = "Count") +
  theme_minimal()


# Create a bar plot of the Sample_Type_qPCR column
ggplot(asv, aes(x = qPCR_Sample_Type)) +
  geom_bar(fill = "lightgreen", color = "black") +
  labs(title = "Distribution of Sample_Type_qPCR Results",
       x = "Sample_Type_qPCR",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Plot the data
ggplot(asv_summary, aes(x = qPCR_Lesion_Type, y = Count)) +
  geom_bar(stat = "identity", fill = "lightgreen", color = "black") +
  labs(title = "Distribution of qPCR and Lesion Type",
       x = "qPCR_Lesion_Type",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Create a summary table with counts
asv_summary <- as.data.frame(table(asv$qPCR_Lesion_Type))
colnames(asv_summary) <- c("qPCR_Lesion_Type", "Count")

# Add the zero count for "positive_crescent"
if (!"positive_crescent" %in% asv_summary$qPCR_Lesion_Type) {
  asv_summary <- rbind(asv_summary, data.frame(qPCR_Lesion_Type = "positive_crescent", Count = 0))
}

# Define the desired order of levels
desired_levels <- c(
  setdiff(asv_summary$qPCR_Lesion_Type, "positive_crescent"),  # All levels except "positive_crescent"
  "positive_crescent"  # Add "positive_crescent" last
)

# Ensure "positive_crescent" is the third column
# Modify this order if there are more specific levels to place before it
if (length(desired_levels) > 2) {
  desired_levels <- c(
    desired_levels[1:2],  # The first two levels
    "positive_crescent",  # The third level
    setdiff(desired_levels, c(desired_levels[1:2], "positive_crescent"))  # The rest
  )
}

# Reorder the factor levels in the summary table
asv_summary$qPCR_Lesion_Type <- factor(asv_summary$qPCR_Lesion_Type, levels = desired_levels)

# Plot the data
ggplot(asv_summary, aes(x = qPCR_Lesion_Type, y = Count)) +
  geom_bar(stat = "identity", fill = "pink", color = "black") +
  labs(title = "Distribution of qPCR and Lesion Type",
       x = "qPCR_Lesion_Type",
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Count the number of samples for each State
state_counts_df <- asv %>%
  group_by(State) %>%
  summarize(Count = n())  # Count occurrences of samples for each State

# Plot the counts
ggplot(state_counts_df, aes(x = State, y = Count, fill = State)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Counts of Samples by State",
       x = "State",
       y = "Count") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "none"  # Remove the legend
  )


# Calculate the number of samples for each State
sample_counts_df <- asv %>%
  group_by(State) %>%
  summarize(Number_of_Samples = n())  # Count the number of samples for each State

# Print the result
print(sample_counts_df)

# Calculate the number of positive samples for each State
positive_samples_df <- asv %>%
  filter(qPCR == "positive") %>%  # Filter rows where aPCR is positive
  group_by(State) %>%              # Group by State
  summarize(Number_of_Positive_Samples = n())  # Count the number of positive samples

# Print the result
print(positive_samples_df)

# Merge the data frames by State
merged_df <- left_join(sample_counts_df, positive_samples_df, by = "State")

# Calculate the ratio of total samples to positive samples
merged_df <- merged_df %>%
  mutate(Ratio = Number_of_Positive_Samples/Number_of_Samples)

# Print the result
print(merged_df)


# Plot the ratio of total samples to positive samples by State
ggplot(merged_df, aes(x = State, y = Ratio, fill = State)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Ratio of Total Samples to Positive Samples by State",
       x = "State",
       y = "Ratio") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "none"  # Remove the legend
  )


###########HEAT MAPS########################################
#filter by state
ak <- asv %>%
  filter(qPCR == "negative")############################################3


# Select only numeric columns
# Identify numeric columns in the 'ak' data frame
numeric_columns <- ak[sapply(ak, is.numeric)]

# Remove the column "Copies.mg" from numeric columns
numeric_columns <- numeric_columns[, !names(numeric_columns) %in% "Copies.mg"]
numeric_columns
# Sum all numeric values in the data frame
total_sum <- sum(numeric_columns, na.rm = TRUE)

# Print the total sum
print(total_sum)

# Calculate the sums for numeric columns
sums <- sapply(ak, function(x) if(is.numeric(x)) sum(x, na.rm = TRUE) else NA)

# Remove the sum for the column "Copies.mg"
sums <- sums[names(sums) != "Copies.mg"]

# Convert the sums to a data frame and add as the first row
sums_df <- as.data.frame(t(sums), stringsAsFactors = FALSE)

# Print the resulting data frame
print(sums_df)

ak <- ak[, !names(ak) %in% "Copies.mg"]
# Combine the sums row with the original data frame
ak_with_sums <- rbind(sums_df, ak)



# Step 2: Divide each sum by the total_sum
normalized_sums <- sums / total_sum

# Step 3: Convert the normalized sums to a data frame and add as the first row
normalized_sums_df <- as.data.frame(t(normalized_sums), stringsAsFactors = FALSE)


# Add a row name for the new row (optional)
rownames(normalized_sums_df) <- "Normalized"


# Step 4: Combine the new row with the original data frame by adding it to the top
ak_with_normalized <- rbind(normalized_sums_df, ak)

# View the updated data frame
View(ak_with_normalized)

normalized_row <- as.numeric(ak_with_normalized["Normalized", ])

# Calculate the sum of the 'Normalized' values
total_normalized <- sum(normalized_row, na.rm = TRUE)

# Print the total sum-should be 1
print(total_normalized)

# Extract the first row from the original data frame
first_row <- ak_with_normalized[1, , drop = FALSE]

# Remove columns with any NA values
cleaned_df_negative <- first_row[, colSums(is.na(first_row)) == 0]

# Verify the cleaned data frame
print(cleaned_df_negative)


######################################################################

combined_df <- rbind(cleaned_df_negative_Not_crescent,  cleaned_df_positive_Not_crescent,cleaned_df_negative_Crescent,cleaned_df_positive_Crescent)


# Define the new row names
new_row_names <- c("Negative_Not_crescent", "Positve_Not_crescent", "Negative_Crescent", "Positive_Crescent"
)

# Assign the new row names to the data frame
rownames(combined_df) <- new_row_names

# Check the updated data frame
print(combined_df)

# Convert to matrix format
heatmap_matrix <- as.matrix(combined_df)

# Convert matrix to long format for ggplot2
heatmap_long <- melt(heatmap_matrix, varnames = c("Row", "Column"))

heatmap_plot <- ggplot(heatmap_long, aes(x = Row, y = Column, fill = value)) +
  geom_tile(color = "black") +  # Add black outlines around each tile
  scale_fill_gradient(low = "white", high = "blue") +  # Adjust colors as needed
  scale_x_discrete(labels = function(x) x) +  # Use row names as x-axis labels
  labs(title = "qPCR Data", x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # Rotate x-axis labels vertically
    legend.title = element_blank()  # Remove legend title
  ) +
  coord_fixed()  # Ensure each tile is a square

print(heatmap_plot)

###################################################################################
