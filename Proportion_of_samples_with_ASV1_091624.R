setwd("C:/Users/Lindsay Alma/Dropbox/Eelgrass/West Coast TagSeq primers/FullDataset")

library(readr)
library(ggplot2)
library(dplyr)


meta <- read_csv("sample_data_ps_091324.csv", show_col_types = FALSE)
View(meta)


# Calculate the number of samples in each State
state_counts <- meta %>%
  group_by(State) %>%
  summarize(total_samples = n())
state_counts

# Determine the number of samples in each State that have a value of 1 for ASV1
asv1_value_counts <- meta %>%
  filter(ASV1 == 1) %>%
  group_by(State) %>%
  summarize(value_one_samples = n())

# Join the two dataframes to calculate the proportion
sample_proportions <- state_counts %>%
  left_join(asv1_value_counts, by = "State") %>%
  mutate(value_one_samples = ifelse(is.na(value_one_samples), 0, value_one_samples),
         proportion = value_one_samples / total_samples)

# Print to check the dataframe
print(sample_proportions)

# Plot with no grid lines and black axis lines
ggplot(sample_proportions, aes(x = State, y = proportion, fill = color)) +
  geom_col() +
  scale_fill_manual(values = c("green", "gray")) +
  labs(x = "Sample Type", y = "Proportion", title = "Proportion of Each Sample Type") +
  theme_minimal() +
  theme(
    legend.position = "none",                # Remove the legend
    panel.grid.major = element_blank(),      # Remove major grid lines
    panel.grid.minor = element_blank(),      # Remove minor grid lines
    axis.line = element_line(color = "black") # Add black lines for each axis
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 100))  # Format y-axis as percentage

################################################################################

# Calculate the number of samples in each Sample_Type
sample_type_counts <- meta %>%
  group_by(Sample_Type) %>%
  summarize(total_samples = n())
sample_type_counts

# Determine the number of samples in each Sample_Type that have a value of 1 for ASV1
asv1_value_counts <- meta %>%
  filter(ASV1 == 1) %>%
  group_by(Sample_Type) %>%
  summarize(value_one_samples = n())

# Join the two dataframes to calculate the proportion
sample_proportions <- sample_type_counts %>%
  left_join(asv1_value_counts, by = "Sample_Type") %>%
  mutate(value_one_samples = ifelse(is.na(value_one_samples), 0, value_one_samples),
         proportion = value_one_samples / total_samples)

# Print to check the dataframe
print(sample_proportions)

#############################################################################
# Calculate the number of samples in each Lesion_Type
sample_type_counts <- meta %>%
  group_by(Lesion_Type) %>%
  summarize(total_samples = n())
sample_type_counts

# Determine the number of samples in each Sample_Type that have a value of 1 for ASV1
asv1_value_counts <- meta %>%
  filter(ASV1 == 1) %>%
  group_by(Lesion_Type) %>%
  summarize(value_one_samples = n())

# Join the two dataframes to calculate the proportion
sample_proportions <- sample_type_counts %>%
  left_join(asv1_value_counts, by = "Lesion_Type") %>%
  mutate(value_one_samples = ifelse(is.na(value_one_samples), 0, value_one_samples),
         proportion = value_one_samples / total_samples)

# Print to check the dataframe
print(sample_proportions)
################################################################################

meta <- meta %>%
  mutate(qPCRASV = case_when(
    P2 == 0 & U2 == 0 ~ "negqPCR negASV1",
    P2 == 0 & U2 == 1 ~ "negqPCR posASV1",
    P2 == 1 & U2 == 0 ~ "posqPCR negASV1",
    P2 == 1 & U2 == 1 ~ "posqPCR posASV1",
    TRUE ~ "Unknown"  # Handle unexpected cases
  ))

# Calculate the proportions of each qPCRASV within each State
sample_proportions <- meta %>%
  group_by(State, qPCRASV1) %>%
  tally() %>%
  group_by(State) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()
sample_proportions
print(sample_proportions, n = Inf)

