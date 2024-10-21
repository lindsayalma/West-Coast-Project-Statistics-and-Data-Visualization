# Load necessary libraries
library(phyloseq)


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
################################################################################################3


# Step 1: Extract OTU (ASV) data and convert it to a data frame
otu_data <- as.data.frame(otu_table(ps))

# Ensure OTUs are columns (transpose if necessary)
if(taxa_are_rows(ps)) {
  otu_data <- t(otu_data)
}

# Step 2: Extract sample data and convert to a data frame
sample_data_df <- as.data.frame(sample_data(ps))
View(sample_data_df)

# Step 3: Combine OTU data and sample data into one data frame
# The rownames of both data frames should match for correct merging
combined_data <- cbind(sample_data_df, otu_data)


###Sample_Type  ###################################################################################3

# Step 4: Ensure that Sample_Type is a factor
combined_data$Sample_Type <- factor(combined_data$Sample_Type, levels = c("L", "G"))

# Step 5: Fit the logistic regression model using the ASVs as predictors
logit_model <- glm(Sample_Type ~ ASV_18S_1 + ASV_18S_2+ ASV_18S_4+ASV_18S_5 + ASV_18S_6, 
                   data = combined_data, 
                   family = binomial(link = "logit"))

# Step 6: View the model summary
summary(logit_model)

# Step 7: Get odds ratios (exponentiated coefficients)
exp(coef(logit_model))

# Step 8: Predict probabilities for each sample being "L"
predicted_probs <- predict(logit_model, type = "response")

# Step 9: Classify samples based on a threshold (e.g., 0.5)
predicted_classes <- ifelse(predicted_probs > 0.5, "L", "G")

# View the predicted classes
predicted_classes


######## Lesion_Type ###################################################################################333


# Step 6: Fit the logistic regression model using the ASVs as predictors
logit_model <- glm(Lesion_Type ~ ASV_18S_6 , 
                   data = combined_data_balanced, 
                   family = binomial(link = "logit"))

# Step 7: View the model summary
summary(logit_model)

# Step 8: Get odds ratios (exponentiated coefficients)
exp(coef(logit_model))

# Step 9: Predict probabilities for each sample being "Crescent"
predicted_probs <- predict(logit_model, type = "response")

# Step 10: Classify samples based on a threshold (e.g., 0.5)
predicted_classes <- ifelse(predicted_probs > 0.5, "Crescent", "Not crescent")

# View the predicted classes
predicted_classes


######## qPCR ###################################################################################333

# Step 4: Ensure that Sample_Type is a factor
combined_data$laby_Positive_Negative <- factor(combined_data$laby_Positive_Negative, levels = c("Positive", "Negative"))

# Step 5: Fit the logistic regression model using the ASVs as predictors
logit_model <- glm(laby_Positive_Negative ~ ASV_18S_1 + ASV_18S_2+ ASV_18S_4+ASV_18S_5 + ASV_18S_6 , 
                   data = combined_data, 
                   family = binomial(link = "logit"))

# Step 6: View the model summary
summary(logit_model)

# Step 7: Get odds ratios (exponentiated coefficients)
exp(coef(logit_model))

# Step 8: Predict probabilities for each sample being "L"
predicted_probs <- predict(logit_model, type = "response")

# Step 9: Classify samples based on a threshold (e.g., 0.5)
predicted_classes <- ifelse(predicted_probs > 0.5, "L", "G")

# View the predicted classes
predicted_classes

########## STATES ###########################################3

# Step 1: Ensure 'State' is a factor
combined_data$State <- factor(combined_data$State, 
                              levels = c("AK", "BB", "BC", "OR", "SD", "WA"))

# Step 2: Fit the logistic regression model with multiple ASV variables and 'State'
logit_model <- glm(laby_Positive_Negative ~ ASV_18S_1 + ASV_18S_2 + ASV_18S_4 + 
                     ASV_18S_5 + ASV_18S_6 + State, 
                   data = combined_data, 
                   family = binomial(link = "logit"))

# Step 3: View the summary of the model
summary(logit_model)


############ PLOTS###############################

# Step 1: Ensure Sample_Type is a factor with the appropriate levels
combined_data$Sample_Type <- factor(combined_data$Sample_Type, levels = c("L", "G"))

# Step 2: Fit the logistic regression model using the ASVs as predictors
logit_model <- glm(Sample_Type ~ ASV_18S_1 + ASV_18S_2 + ASV_18S_4 + ASV_18S_5 + ASV_18S_6, 
                   data = combined_data, 
                   family = binomial(link = "logit"),
                   na.action = na.exclude)  # Exclude NA cases during fitting

# Step 3: Predict probabilities using the full data set
predicted_probs <- predict(logit_model, type = "response", newdata = combined_data)

# Step 4: Create the plot_data data frame
plot_data <- data.frame(
  Sample_Type = combined_data$Sample_Type,
  Predicted_Probabilities = predicted_probs
)

# Step 5: Load ggplot2 and create the plot
library(ggplot2)

ggplot(plot_data, aes(x = Sample_Type, y = Predicted_Probabilities, color = Sample_Type)) +
  geom_point(size = 3) +  # Use points to visualize predicted probabilities
  geom_jitter(width = 0.1, height = 0) +  # Jitter for better visibility
  scale_color_manual(values = c("L" = "purple", "G" = "orange")) +  # Customize colors
  theme_minimal() +
  labs(title = "Predicted Probabilities by Sample Type",
       x = "Sample Type",
       y = "Predicted Probability of Being 'L'") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +  # Add threshold line
  theme(legend.position = "none")  # Remove legend if not needed


# Load necessary libraries
library(ggplot2)

# Create a boxplot of predicted probabilities by Sample Type
ggplot(plot_data, aes(x = Sample_Type, y = Predicted_Probabilities, fill = Sample_Type)) +
  geom_boxplot() +  # Create boxplots
  scale_fill_manual(values = c("L" = "blue", "G" = "orange")) +  # Set custom colors
  theme_minimal() +
  labs(title = "Predicted Probabilities by Sample Type",
       x = "Sample Type",
       y = "Predicted Probability of Being 'L'") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red")  # Add threshold line


# Create a violin plot of predicted probabilities by Sample Type
ggplot(plot_data, aes(x = Sample_Type, y = Predicted_Probabilities, fill = Sample_Type)) +
  geom_violin(trim = FALSE) +  # Create violin plots
  scale_fill_manual(values = c("L" = "blue", "G" = "orange")) +  # Set custom colors
  theme_minimal() +
  labs(title = "Predicted Probabilities by Sample Type",
       x = "Sample Type",
       y = "Predicted Probability of Being 'L'") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red")  # Add threshold line

# Create a scatter plot with jitter
ggplot(plot_data, aes(x = Sample_Type, y = Predicted_Probabilities, color = Sample_Type)) +
  geom_jitter(width = 0.2, alpha = 0.6) +  # Add jitter to points
  scale_color_manual(values = c("L" = "blue", "G" = "orange")) +  # Set custom colors
  theme_minimal() +
  labs(title = "Predicted Probabilities by Sample Type",
       x = "Sample Type",
       y = "Predicted Probability of Being 'L'") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red")  # Add threshold line


# Create a combined plot
ggplot(plot_data, aes(x = Sample_Type, y = Predicted_Probabilities)) +
  geom_boxplot(aes(fill = Sample_Type), outlier.shape = NA) +  # Boxplot without outliers
  geom_jitter(aes(color = Sample_Type), width = 0.2, alpha = 0.6) +  # Add jittered points
  scale_fill_manual(values = c("L" = "blue", "G" = "orange")) +  # Set custom colors
  scale_color_manual(values = c("L" = "blue", "G" = "orange")) +  # Set custom colors
  theme_minimal() +
  labs(title = "Predicted Probabilities by Sample Type",
       x = "Sample Type",
       y = "Predicted Probability of Being 'L'") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red")  # Add threshold line


########### S curve plots #########################################################################3


print(table(combined_data$Sample_Type))
View(combined_data)

combined_data <- subset(combined_data, !is.na(ASV_18S_1))
combined_data <- subset(combined_data, !is.na(Sample_Type))

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Step 1: Remove rows with NA in ASV_18S_1 and Sample_Type
combined_data <- subset(combined_data, !is.na(ASV_18S_1) & !is.na(Sample_Type))

# Step 2: Ensure Sample_Type is a factor
combined_data$Sample_Type <- factor(combined_data$Sample_Type, levels = c("G", "L"))

# Step 3: Create a binary response variable for logistic regression
combined_data$Sample_Type_binary <- as.numeric(combined_data$Sample_Type == "L")  # 1 for "L", 0 for "G"

# Check the distribution of the new binary variable
print(table(combined_data$Sample_Type_binary))

# Step 4: Fit the logistic regression model using the binary response variable
logit_model <- glm(Sample_Type_binary ~ ASV_18S_1 + ASV_18S_2 + ASV_18S_4 + ASV_18S_5 + ASV_18S_6, 
                   data = combined_data, 
                   family = binomial(link = "logit"))

# Step 5: Predict probabilities based on the fitted model
combined_data$predicted_prob <- predict(logit_model, type = "response")



#####################S curve for ASV vs Sample_Type ###################################################################3


# Assume the logistic model is already fitted as logit_model
# Create a sequence of ASV_18S_1 values for predictions
asv_values <- seq(min(combined_data$ASV_18S_1, na.rm = TRUE), 
                  max(combined_data$ASV_18S_1, na.rm = TRUE), 
                  length.out = 100)

# Prepare a new data frame for predictions
new_data <- data.frame(ASV_18S_1 = asv_values,
                       ASV_18S_2 = mean(combined_data$ASV_18S_2, na.rm = TRUE), 
                       ASV_18S_4 = mean(combined_data$ASV_18S_4, na.rm = TRUE),
                       ASV_18S_5 = mean(combined_data$ASV_18S_5, na.rm = TRUE),
                       ASV_18S_6 = mean(combined_data$ASV_18S_6, na.rm = TRUE))

# Get predicted probabilities and standard errors
pred_probs <- predict(logit_model, newdata = new_data, type = "response", se.fit = TRUE)

# Extract the fitted probabilities and standard errors
pm <- pred_probs$fit
pu <- pm + (1.96 * pred_probs$se.fit)  # Upper 95% CI
pl <- pm - (1.96 * pred_probs$se.fit)  # Lower 95% CI

# Create the plot
plot(combined_data$ASV_18S_1, combined_data$Sample_Type_binary,
     pch = 16,
     cex = 1,
     ylab = "Probability of Sample Type (L)",
     xlab = "ASV_18S_1",
     ylim = c(0, 1),  # Set y limits from 0 to 1
     main = "Predicted Probability of Sample Type by ASV_18S_1"
)

# Add grid
grid()

# Add confidence interval as a polygon
polygon(c(rev(asv_values), asv_values), 
        c(rev(pl), pu), 
        col = "grey90", border = NA)

# Add fitted line for predicted probabilities
lines(asv_values, pm, lwd = 2)
lines(asv_values, pu, lwd = 2, col = "red")
lines(asv_values, pl, lwd = 2, col = "red")

# Add horizontal lines for probability thresholds
abline(h = 0.1, lty = 2)
abline(h = 0.5, lty = 2)
abline(h = 0.9, lty = 2)

########################################################################################


# Load necessary libraries
library(ggplot2)
library(dplyr)

# Check the distribution of the original Sample_Type
print(table(combined_data$Sample_Type))
View(combined_data)

# Step 1: Remove rows with NA in ASV_18S_1 and Lesion_Type
combined_data <- subset(combined_data, !is.na(ASV_18S_1) & !is.na(Sample_Type))

# Step 2: Replace Sample_Type with Lesion_Type
combined_data$Lesion_Type <- factor(combined_data$Sample_Type, levels = c("Crescent", "Not crescent"))

# Step 3: Create a binary response variable for logistic regression
combined_data$Lesion_Type_binary <- as.numeric(combined_data$Lesion_Type == "Crescent")  # 1 for "Crescent", 0 for "Not crescent"

# Check the distribution of the new binary variable
print(table(combined_data$Lesion_Type_binary))

# Step 4: Fit the logistic regression model using the binary response variable
logit_model <- glm(Lesion_Type_binary ~ ASV_18S_1 + ASV_18S_1 + ASV_18S_4 + ASV_18S_5 + ASV_18S_6, 
                   data = combined_data, 
                   family = binomial(link = "logit"))

# Step 5: Predict probabilities based on the fitted model
combined_data$predicted_prob <- predict(logit_model, type = "response")

##################### S curve for ASV vs Lesion_Type #############################

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Check the distribution of the original Sample_Type
print(table(combined_data$Sample_Type))


# Step 1: Remove rows with NA in ASV_18S_1 and Sample_Type
combined_data <- subset(combined_data, !is.na(ASV_18S_1) & !is.na(Sample_Type))

# Step 2: Ensure Sample_Type is a factor
combined_data$Sample_Type <- factor(combined_data$Sample_Type, levels = c("G", "L"))

# Step 3: Create a binary response variable for logistic regression
combined_data$Sample_Type_binary <- as.numeric(combined_data$Sample_Type == "L")  # 1 for "L", 0 for "G"

# Check the distribution of the new binary variable
print(table(combined_data$Sample_Type_binary))

# Step 4: Fit the logistic regression model using the binary response variable
logit_model <- glm(Sample_Type_binary ~ ASV_18S_1, 
                   data = combined_data, 
                   family = binomial(link = "logit"))

# Step 5: Predict probabilities based on the fitted model
combined_data$predicted_prob <- predict(logit_model, type = "response")

##################### S curve for ASV vs Sample_Type #############################

# Create a sequence of ASV_18S_1 values for predictions
asv_values <- seq(min(combined_data$ASV_18S_1, na.rm = TRUE), 
                  max(combined_data$ASV_18S_1, na.rm = TRUE), 
                  length.out = 100)

# Prepare a new data frame for predictions
new_data <- data.frame(ASV_18S_1 = asv_values)

# Get predicted probabilities and standard errors
pred_probs <- predict(logit_model, newdata = new_data, type = "response", se.fit = TRUE)

# Extract the fitted probabilities and standard errors
pm <- pred_probs$fit
pu <- pm + (1.96 * pred_probs$se.fit)  # Upper 95% CI
pl <- pm - (1.96 * pred_probs$se.fit)  # Lower 95% CI

# Create the plot
plot(combined_data$ASV_18S_1, combined_data$Sample_Type_binary,
     pch = 16,
     cex = 1,
     ylab = "Probability of Sample Type (L)",
     xlab = "ASV_18S_1",
     ylim = c(0, 1),  # Set y limits from 0 to 1
     main = "Predicted Probability of Sample Type by ASV_18S_1"
)

# Add grid
grid()

# Add confidence interval as a polygon
polygon(c(rev(asv_values), asv_values), 
        c(rev(pl), pu), 
        col = "grey90", border = NA)

# Add fitted line for predicted probabilities
lines(asv_values, pm, lwd = 2)
lines(asv_values, pu, lwd = 2, col = "red")
lines(asv_values, pl, lwd = 2, col = "red")

# Add horizontal lines for probability thresholds
abline(h = 0.1, lty = 2)
abline(h = 0.5, lty = 2)
abline(h = 0.9, lty = 2)



######### multuple lines for multuple asvs##################################################333

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Step 1: Remove rows with NA in ASV columns and Sample_Type
combined_data <- subset(combined_data, 
                        !is.na(ASV_18S_1) & 
                          !is.na(ASV_18S_2) & 
                          !is.na(ASV_18S_4) & 
                          !is.na(ASV_18S_5) & 
                          !is.na(ASV_18S_6) & 
                          !is.na(Sample_Type))

# Step 2: Ensure Sample_Type is a factor
combined_data$Sample_Type <- factor(combined_data$Sample_Type, levels = c("G", "L"))

# Step 3: Create a binary response variable for logistic regression
combined_data$Sample_Type_binary <- as.numeric(combined_data$Sample_Type == "L")  # 1 for "L", 0 for "G"

# Step 4: Fit the logistic regression model for each ASV variable and store results
asv_vars <- c("ASV_18S_1", "ASV_18S_2", "ASV_18S_4", "ASV_18S_5", "ASV_18S_6")
colors <- c("purple", "orange", "cyan", "forestgreen", "blue")  # Use specified colors

# Adjust margins: bottom, left, top, right
par(mar = c(5, 5, 4, 5))  # Increase bottom margin for legend

# Create an empty plot
plot(NULL, xlim = c(min(combined_data$ASV_18S_1, na.rm = TRUE), 
                    max(combined_data$ASV_18S_1, na.rm = TRUE)), 
     ylim = c(0, 1), 
     xlab = "ASV Values", 
     ylab = "Probability of Sample Type (L)", 
     main = "Predicted Probability of Sample Type by ASV")

# Add grid
grid()

# Loop through each ASV variable to fit a model and add points to the plot
for (i in seq_along(asv_vars)) {
  asv_var <- asv_vars[i]
  
  # Plot points for "L" samples
  points(combined_data[[asv_var]][combined_data$Sample_Type == "L"], 
         rep(1, sum(combined_data$Sample_Type == "L")),  # Set y-axis to 1 for "L"
         pch = 16, 
         col = adjustcolor(colors[i], alpha.f = 0.7))  # Color by ASV
  
  # Plot points for "G" samples
  points(combined_data[[asv_var]][combined_data$Sample_Type == "G"], 
         rep(0, sum(combined_data$Sample_Type == "G")),  # Set y-axis to 0 for "G"
         pch = 16, 
         col = adjustcolor(colors[i], alpha.f = 0.7))  # Color by ASV
  
  # Fit the logistic regression model
  logit_model <- glm(Sample_Type_binary ~ get(asv_var), 
                     data = combined_data, 
                     family = binomial(link = "logit"))
  
  # Create a sequence of ASV values for predictions
  asv_values <- seq(min(combined_data[[asv_var]], na.rm = TRUE), 
                    max(combined_data[[asv_var]], na.rm = TRUE), 
                    length.out = 100)
  
  # Prepare a new data frame for predictions
  new_data <- data.frame(asv_value = asv_values)
  
  # Create a new column in new_data for the current ASV variable
  new_data[[asv_var]] <- new_data$asv_value
  
  # Get predicted probabilities
  pred_probs <- predict(logit_model, newdata = new_data, type = "response", se.fit = TRUE)
  
  # Extract the fitted probabilities and standard errors
  pm <- pred_probs$fit
  pu <- pm + (1.96 * pred_probs$se.fit)  # Upper 95% CI
  pl <- pm - (1.96 * pred_probs$se.fit)  # Lower 95% CI
  
  # Add confidence interval as a polygon
  polygon(c(rev(asv_values), asv_values), 
          c(rev(pl), pu), 
          col = adjustcolor(colors[i], alpha.f = 0.2), border = NA)
  
  # Add fitted line for predicted probabilities
  lines(asv_values, pm, lwd = 2, col = colors[i])
}

# Add horizontal lines for probability thresholds
abline(h = 0.1, lty = 2)
abline(h = 0.5, lty = 2)
abline(h = 0.9, lty = 2)

# Add a legend positioned to the left of the plot
legend("topright", 
       legend = asv_vars, 
       col = colors, 
       lty = 1, 
       lwd = 2, 
       bty = "n", 
       inset = c(-0.15, 0),  # Move the legend to the left
       xpd = TRUE,           # Allow drawing outside plot region
       cex = 0.5)            # Set the legend text size smaller


