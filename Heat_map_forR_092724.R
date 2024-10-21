library(pheatmap)
library(phyloseq)


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

# Verify the changes in the tax_table
#View(tax_table(ps))

ps <- subset_samples(ps, Lesion_Type=="Crescent" | Lesion_Type=="Not crescent")

# Step 2: Read the new OTU table
# Replace with your actual file path and make sure it's read in correctly
new_otu_table <- read.csv("asv_data_raw_top20_percent_083024.csv", row.names = 1)

# Step 3: Convert to matrix (if necessary)
otu_matrix <- as.matrix(new_otu_table)

# Step 4: Create a new otu_table object
# Set the taxa as rows (features) and samples as columns
otu_table_new <- otu_table(otu_matrix, taxa_are_rows = TRUE)

# Step 5: Extract relevant sample names from the new OTU table
# Assuming the sample names are in the column names of the new OTU table
# Use string manipulation to match the sample names after the last underscore
new_sample_names <- sub(".*_(.*)", "\\1", colnames(otu_matrix))

# Update the sample names in the new OTU table
colnames(otu_table_new) <- new_sample_names

ps <- merge_phyloseq(otu_table_new, tax_table(ps), sample_data(ps))

# Step 7: Check if the new OTU table has been imported correctly
print(otu_table(ps))

################################BB####################################################################33
# Subset by State to keep only samples where State is ""
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

# Verify the changes in the tax_table
#View(tax_table(ps))

ps <- subset_samples(ps, Lesion_Type=="Crescent" | Lesion_Type=="Not crescent")

# Step 2: Read the new OTU table
# Replace with your actual file path and make sure it's read in correctly
new_otu_table <- read.csv("asv_data_raw_top20_percent_083024.csv", row.names = 1)

# Step 3: Convert to matrix (if necessary)
otu_matrix <- as.matrix(new_otu_table)

# Step 4: Create a new otu_table object
# Set the taxa as rows (features) and samples as columns
otu_table_new <- otu_table(otu_matrix, taxa_are_rows = TRUE)

# Step 5: Extract relevant sample names from the new OTU table
# Assuming the sample names are in the column names of the new OTU table
# Use string manipulation to match the sample names after the last underscore
new_sample_names <- sub(".*_(.*)", "\\1", colnames(otu_matrix))

# Update the sample names in the new OTU table
colnames(otu_table_new) <- new_sample_names

ps <- merge_phyloseq(otu_table_new, tax_table(ps), sample_data(ps))

# Step 7: Check if the new OTU table has been imported correctly
print(otu_table(ps))

ps <- subset_samples(ps, State == "BB")

#Extract OTU (abundance) matrix
otu_matrix <- as.matrix(otu_table(ps))

# Extract taxonomic information at the ASV level (or desired level)
asv_names <- taxa_names(ps)  # Get ASV names

# Remove "18S_" from ASV names
asv_names <- gsub("18S_", "", asv_names)
rownames(otu_matrix) <- make.unique(as.character(asv_names))  # Set row names as cleaned ASVs

# Extract metadata for annotations
metadata <- data.frame(sample_data(ps))

# Prepare column annotations
annotation_col <- data.frame(
  Site = metadata$Site,
  Sample_Type = metadata$Sample_Type,
  Lesion_Type = metadata$Lesion_Type
)
rownames(annotation_col) <- rownames(metadata)
unique(metadata$Site)

# Define colors for annotations
ann_colors <- list(
  Site = c("A" = "yellow", "D" = "plum1", "F" = "cyan"),
  Sample_Type = c("L" = "orange4", "G" = "green"),  
  Lesion_Type = c("Crescent" = "cornsilk4", "Not crescent" = "forestgreen")
)

# Define the custom color palette for the heatmap
color_palette <- colorRampPalette(c("white", "purple",  "darkorange"))(100)


pheatmap(
  otu_matrix, 
  annotation_col = annotation_col, 
  annotation_colors = ann_colors,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,  # Show ASV names on the y-axis
  show_colnames = FALSE,  # Remove x-axis labels
  fontsize_row = 10,  
  fontsize_col = 10,
  main = "Bodega Bay",  # Add title here
  color = color_palette  # Apply the custom color palette
)


###############################BC####################################################################33
# Subset by State to keep only samples where State is ""
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

# Verify the changes in the tax_table
#View(tax_table(ps))

ps <- subset_samples(ps, Lesion_Type=="Crescent" | Lesion_Type=="Not crescent")

# Step 2: Read the new OTU table
# Replace with your actual file path and make sure it's read in correctly
new_otu_table <- read.csv("asv_data_raw_top20_percent_083024.csv", row.names = 1)

# Step 3: Convert to matrix (if necessary)
otu_matrix <- as.matrix(new_otu_table)

# Step 4: Create a new otu_table object
# Set the taxa as rows (features) and samples as columns
otu_table_new <- otu_table(otu_matrix, taxa_are_rows = TRUE)

# Step 5: Extract relevant sample names from the new OTU table
# Assuming the sample names are in the column names of the new OTU table
# Use string manipulation to match the sample names after the last underscore
new_sample_names <- sub(".*_(.*)", "\\1", colnames(otu_matrix))

# Update the sample names in the new OTU table
colnames(otu_table_new) <- new_sample_names

ps <- merge_phyloseq(otu_table_new, tax_table(ps), sample_data(ps))

# Step 7: Check if the new OTU table has been imported correctly
print(otu_table(ps))

ps <- subset_samples(ps, State == "BC")

#Extract OTU (abundance) matrix
otu_matrix <- as.matrix(otu_table(ps))

# Extract taxonomic information at the ASV level (or desired level)
asv_names <- taxa_names(ps)  # Get ASV names

# Remove "18S_" from ASV names
asv_names <- gsub("18S_", "", asv_names)
rownames(otu_matrix) <- make.unique(as.character(asv_names))  # Set row names as cleaned ASVs

# Extract metadata for annotations
metadata <- data.frame(sample_data(ps))

# Prepare column annotations
annotation_col <- data.frame(
  Site = metadata$Site,
  Sample_Type = metadata$Sample_Type,
  Lesion_Type = metadata$Lesion_Type
)
rownames(annotation_col) <- rownames(metadata)
unique(metadata$Site)

# Define colors for annotations
ann_colors <- list(
  Site = c("A" = "yellow", "D" = "plum1", "E" = "cyan", "H" = "blue", "Q" = "darkcyan"),
  Sample_Type = c("L" = "orange4", "G" = "green"),  
  Lesion_Type = c("Crescent" = "cornsilk4", "Not crescent" = "forestgreen")
)

# Define the custom color palette for the heatmap
color_palette <- colorRampPalette(c("white", "purple",  "darkorange"))(100)


pheatmap(
  otu_matrix, 
  annotation_col = annotation_col, 
  annotation_colors = ann_colors,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,  # Show ASV names on the y-axis
  show_colnames = FALSE,  # Remove x-axis labels
  fontsize_row = 10,  
  fontsize_col = 10,
  main = "British Columbia",  # Add title here
  color = color_palette  # Apply the custom color palette
)



########################OR#################################


# Subset by State to keep only samples where State is ""
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

# Verify the changes in the tax_table
#View(tax_table(ps))

ps <- subset_samples(ps, Lesion_Type=="Crescent" | Lesion_Type=="Not crescent")

# Step 2: Read the new OTU table
# Replace with your actual file path and make sure it's read in correctly
new_otu_table <- read.csv("asv_data_raw_top20_percent_083024.csv", row.names = 1)

# Step 3: Convert to matrix (if necessary)
otu_matrix <- as.matrix(new_otu_table)

# Step 4: Create a new otu_table object
# Set the taxa as rows (features) and samples as columns
otu_table_new <- otu_table(otu_matrix, taxa_are_rows = TRUE)

# Step 5: Extract relevant sample names from the new OTU table
# Assuming the sample names are in the column names of the new OTU table
# Use string manipulation to match the sample names after the last underscore
new_sample_names <- sub(".*_(.*)", "\\1", colnames(otu_matrix))

# Update the sample names in the new OTU table
colnames(otu_table_new) <- new_sample_names

ps <- merge_phyloseq(otu_table_new, tax_table(ps), sample_data(ps))

# Step 7: Check if the new OTU table has been imported correctly
print(otu_table(ps))

ps <- subset_samples(ps, State == "OR")

#Extract OTU (abundance) matrix
otu_matrix <- as.matrix(otu_table(ps))

# Extract taxonomic information at the ASV level (or desired level)
asv_names <- taxa_names(ps)  # Get ASV names

# Remove "18S_" from ASV names
asv_names <- gsub("18S_", "", asv_names)
rownames(otu_matrix) <- make.unique(as.character(asv_names))  # Set row names as cleaned ASVs

# Extract metadata for annotations
metadata <- data.frame(sample_data(ps))

# Prepare column annotations
annotation_col <- data.frame(
  Site = metadata$Site,
  Sample_Type = metadata$Sample_Type,
  Lesion_Type = metadata$Lesion_Type
)
rownames(annotation_col) <- rownames(metadata)
unique(metadata$Site)

# Define colors for annotations
ann_colors <- list(
  Site = c("B" = "yellow", "C" = "plum1", "D" = "cyan"),
  Sample_Type = c("L" = "orange4", "G" = "green"),  
  Lesion_Type = c("Crescent" = "cornsilk4", "Not crescent" = "forestgreen")
)

# Define the custom color palette for the heatmap
color_palette <- colorRampPalette(c("white", "purple",  "darkorange"))(100)


pheatmap(
  otu_matrix, 
  annotation_col = annotation_col, 
  annotation_colors = ann_colors,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,  # Show ASV names on the y-axis
  show_colnames = FALSE,  # Remove x-axis labels
  fontsize_row = 10,  
  fontsize_col = 10,
  main = "Oregon",  # Add title here
  color = color_palette  # Apply the custom color palette
)




########################SD#################################


# Subset by State to keep only samples where State is ""
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

# Verify the changes in the tax_table
#View(tax_table(ps))

ps <- subset_samples(ps, Lesion_Type=="Crescent" | Lesion_Type=="Not crescent")

# Step 2: Read the new OTU table
# Replace with your actual file path and make sure it's read in correctly
new_otu_table <- read.csv("asv_data_raw_top20_percent_083024.csv", row.names = 1)

# Step 3: Convert to matrix (if necessary)
otu_matrix <- as.matrix(new_otu_table)

# Step 4: Create a new otu_table object
# Set the taxa as rows (features) and samples as columns
otu_table_new <- otu_table(otu_matrix, taxa_are_rows = TRUE)

# Step 5: Extract relevant sample names from the new OTU table
# Assuming the sample names are in the column names of the new OTU table
# Use string manipulation to match the sample names after the last underscore
new_sample_names <- sub(".*_(.*)", "\\1", colnames(otu_matrix))

# Update the sample names in the new OTU table
colnames(otu_table_new) <- new_sample_names

ps <- merge_phyloseq(otu_table_new, tax_table(ps), sample_data(ps))

# Step 7: Check if the new OTU table has been imported correctly
print(otu_table(ps))

ps <- subset_samples(ps, State == "SD")

#Extract OTU (abundance) matrix
otu_matrix <- as.matrix(otu_table(ps))

# Extract taxonomic information at the ASV level (or desired level)
asv_names <- taxa_names(ps)  # Get ASV names

# Remove "18S_" from ASV names
asv_names <- gsub("18S_", "", asv_names)
rownames(otu_matrix) <- make.unique(as.character(asv_names))  # Set row names as cleaned ASVs

# Extract metadata for annotations
metadata <- data.frame(sample_data(ps))

# Prepare column annotations
annotation_col <- data.frame(
  Site = metadata$Site,
  Sample_Type = metadata$Sample_Type,
  Lesion_Type = metadata$Lesion_Type
)
rownames(annotation_col) <- rownames(metadata)
unique(metadata$Site)

# Define colors for annotations
ann_colors <- list(
  Site = c("A" = "yellow", "B" = "plum1", "C" = "cyan"),
  Sample_Type = c("L" = "orange4", "G" = "green"),  
  Lesion_Type = c("Crescent" = "cornsilk4", "Not crescent" = "forestgreen")
)

# Define the custom color palette for the heatmap
color_palette <- colorRampPalette(c("white", "purple",  "darkorange"))(100)


pheatmap(
  otu_matrix, 
  annotation_col = annotation_col, 
  annotation_colors = ann_colors,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,  # Show ASV names on the y-axis
  show_colnames = FALSE,  # Remove x-axis labels
  fontsize_row = 10,  
  fontsize_col = 10,
  main = "San Diego",  # Add title here
  color = color_palette  # Apply the custom color palette
)




########################WA#################################


# Subset by State to keep only samples where State is ""
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

# Verify the changes in the tax_table
#View(tax_table(ps))

ps <- subset_samples(ps, Lesion_Type=="Crescent" | Lesion_Type=="Not crescent")

# Step 2: Read the new OTU table
# Replace with your actual file path and make sure it's read in correctly
new_otu_table <- read.csv("asv_data_raw_top20_percent_083024.csv", row.names = 1)

# Step 3: Convert to matrix (if necessary)
otu_matrix <- as.matrix(new_otu_table)

# Step 4: Create a new otu_table object
# Set the taxa as rows (features) and samples as columns
otu_table_new <- otu_table(otu_matrix, taxa_are_rows = TRUE)

# Step 5: Extract relevant sample names from the new OTU table
# Assuming the sample names are in the column names of the new OTU table
# Use string manipulation to match the sample names after the last underscore
new_sample_names <- sub(".*_(.*)", "\\1", colnames(otu_matrix))

# Update the sample names in the new OTU table
colnames(otu_table_new) <- new_sample_names

ps <- merge_phyloseq(otu_table_new, tax_table(ps), sample_data(ps))

# Step 7: Check if the new OTU table has been imported correctly
print(otu_table(ps))

ps <- subset_samples(ps, State == "WA")

#Extract OTU (abundance) matrix
otu_matrix <- as.matrix(otu_table(ps))

# Extract taxonomic information at the ASV level (or desired level)
asv_names <- taxa_names(ps)  # Get ASV names

# Remove "18S_" from ASV names
asv_names <- gsub("18S_", "", asv_names)
rownames(otu_matrix) <- make.unique(as.character(asv_names))  # Set row names as cleaned ASVs

# Extract metadata for annotations
metadata <- data.frame(sample_data(ps))

# Prepare column annotations
annotation_col <- data.frame(
  Site = metadata$Site,
  Sample_Type = metadata$Sample_Type,
  Lesion_Type = metadata$Lesion_Type
)
rownames(annotation_col) <- rownames(metadata)
unique(metadata$Site)

# Define colors for annotations
ann_colors <- list(
  Site = c("A" = "yellow", "B" = "plum1", "E" = "cyan", "F"="blue"),
  Sample_Type = c("L" = "orange4", "G" = "green"),  
  Lesion_Type = c("Crescent" = "cornsilk4", "Not crescent" = "forestgreen")
)

# Define the custom color palette for the heatmap
color_palette <- colorRampPalette(c("white", "purple",  "darkorange"))(100)


pheatmap(
  otu_matrix, 
  annotation_col = annotation_col, 
  annotation_colors = ann_colors,
  cluster_rows = TRUE, 
  cluster_cols = TRUE,
  show_rownames = TRUE,  # Show ASV names on the y-axis
  show_colnames = FALSE,  # Remove x-axis labels
  fontsize_row = 10,  
  fontsize_col = 10,
  main = "Washington",  # Add title here
  color = color_palette  # Apply the custom color palette
)
