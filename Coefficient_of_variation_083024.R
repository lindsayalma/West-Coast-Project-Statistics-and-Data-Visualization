asv_table <- read.csv("asv_data_raw_top20_082624.csv", row.names=1)
metadata <- read.csv("sapmple_data_ps_083024.csv")

# Function to calculate CV
calculate_cv <- function(x) {
  mean_x <- mean(x)
  sd_x <- sd(x)
  cv <- (sd_x / mean_x) * 100
  return(cv)
}

# Apply function to each column (sample) in the ASV table
cv_results <- apply(asv_table, 2, calculate_cv)

# Print the CV for each sample
print(cv_results)

# Convert cv_results to a matrix
cv_results_matrix <- matrix(cv_results, nrow = 1, dimnames = list(NULL, names(cv_results)))

cv_results_matrix_transposed <- t(cv_results_matrix)

# Convert the transposed matrix to a data frame
cv_results_df <- as.data.frame(cv_results_matrix_transposed)

# Rename the column to 'CV'
colnames(cv_results_df) <- 'CV'

View(cv_results_df)

# Ensure 'Sample' column in cv_results_df
cv_results_df$Sample <- rownames(cv_results_df)

# Reorder columns for merging
cv_results_df <- cv_results_df[, c("Sample", "CV")]

# Merge with metadata
CV <- merge(cv_results_df, metadata, by = "Sample")
CV <- CV[is.finite(CV$CV), ]

View(CV)

# Create the boxplot
ggplot(CV, aes(x = Lesion_Type, y = CV)) +
  geom_boxplot() +
  labs(title = "Coefficient of Variation by Lesion_Type",
       x = "Lesion_Type",
       y = "Coefficient of Variation (%)") +
  theme_minimal()



# 3. Non-parametric alternative: Kruskal-Wallis Test
kruskal_test_result <- kruskal.test(CV ~ Lesion_Type, data = CV)
print(kruskal_test_result)


# 1. Calculate summary statistics: mean and standard error for each state
cv_summary <- CV %>%
  group_by(State) %>%
  summarise(
    mean_CV = mean(CV, na.rm = TRUE),
    se_CV = sd(CV, na.rm = TRUE) / sqrt(n())
  )

# 2. Create the bar plot with error bars
ggplot(cv_summary, aes(x = State, y = mean_CV)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  geom_errorbar(aes(ymin = mean_CV - se_CV, ymax = mean_CV + se_CV), width = 0.2) +
  labs(title = "Coefficient of Variation by State",
       x = "State",
       y = "Mean Coefficient of Variation (%)") +
  theme_minimal()


