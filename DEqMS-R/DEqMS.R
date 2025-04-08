
# Load required libraries
library(DEqMS)
library(limma)
library(tidyverse)

# Load the metadata file
metadata <- read.csv("meta_data2.csv", row.names = 1)

### refrence WTSham
# Convert Genotype and Treatment into factors with appropriate reference levels
metadata$Genotype <- factor(metadata$Genotype, levels = c("WT", "Btg", "BK"))  # WT as reference
metadata$Treatment <- factor(metadata$Treatment, levels = c("Sham", "ShamTB", "LPSTB", "LPS"))  # Sham as reference

# Create the design matrix based on Genotype and Treatment
design <- model.matrix(~ Genotype + Treatment, data = metadata)

# Load your normalized expression data (rows are genes/proteins, columns are samples)
expression_data <- read.csv("Mito_norm_total.csv", row.names = 1)

# The row names of expression_data are your Accession numbers (protein IDs)
accession_numbers <- rownames(expression_data)

# Fit the linear model using the expression data and design matrix
fit <- lmFit(expression_data, design)

# Create the design matrix for Genotype and Treatment interactions
design_group <- model.matrix(~ Genotype:Treatment , data = metadata)

# Check and clean column names of the design matrix
colnames(design_group) <- make.names(colnames(design_group))

# Define the contrast matrix using only estimable coefficients
contrast.matrix <- makeContrasts(
  BKvsWTInSham = GenotypeBK.TreatmentSham - GenotypeWT.TreatmentSham,
  BtgvsWTInSham = GenotypeBtg.TreatmentSham - GenotypeWT.TreatmentSham,
  ShamTBvsShamInWT = GenotypeWT.TreatmentShamTB - GenotypeWT.TreatmentSham,
  LPSvsShamInWT = GenotypeWT.TreatmentLPS - GenotypeWT.TreatmentSham,
  levels = design_group
)

# Fit the model using the expression data and design matrix
fit <- lmFit(expression_data, design_group)

# Apply the contrast matrix
fit2 <- contrasts.fit(fit, contrast.matrix)

# Apply empirical Bayes moderation
fit2 <- eBayes(fit2)

# Extract results for each contrast and add the Accession number (protein IDs) to the results
results_BKvsWTinSham <- topTable(fit2, coef="BKvsWTInSham", adjust="BH",number = Inf) %>%
  rownames_to_column(var = "Accession")

results_BtgvsWTinSham <- topTable(fit2, coef="BtgvsWTInSham", adjust="BH",number = Inf) %>%
  rownames_to_column(var = "Accession")

results_ShamTBvsShamInWT <- topTable(fit2, coef="ShamTBvsShamInWT", adjust="BH",number = Inf) %>%
  rownames_to_column(var = "Accession")

results_LPSvsShamInWT <- topTable(fit2, coef="LPSvsShamInWT", adjust="BH",number = Inf) %>%
  rownames_to_column(var = "Accession")

# Save results as CSV files with Accession number included
write.csv(results_BKvsWTinSham, file="results_BKvsWTInSham.csv", row.names=FALSE)
write.csv(results_BtgvsWTinSham, file="results_BtgvsWTInSham.csv", row.names=FALSE)
write.csv(results_ShamTBvsShamInWT, file="results_ShamTBvsShamInWT.csv", row.names=FALSE)
write.csv(results_LPSvsShamInWT, file="results_LPSvsShamInWT.csv", row.names=FALSE)


######

### refrence BtgSham


# Load required libraries
library(DEqMS)
library(limma)
library(tidyverse)

# Load the metadata file
metadata <- read.csv("meta_data2.csv", row.names = 1)

# Convert Genotype and Treatment into factors with appropriate reference levels
metadata$Genotype <- factor(metadata$Genotype, levels = c("Btg", "WT", "BK"))  # Btg as reference
metadata$Treatment <- factor(metadata$Treatment, levels = c("Sham", "ShamTB", "LPSTB", "LPS"))  # Sham as reference

# Create the design matrix based on Genotype and Treatment
design <- model.matrix(~ Genotype + Treatment, data = metadata)

# Load your normalized expression data (rows are genes/proteins, columns are samples)
expression_data <- read.csv("Mito_norm_total.csv", row.names = 1)

# The row names of expression_data are your Accession numbers (protein IDs)
accession_numbers <- rownames(expression_data)

# Fit the linear model using the expression data and design matrix
fit <- lmFit(expression_data, design)

# Create the design matrix for Genotype and Treatment interactions
design_group <- model.matrix(~ Genotype:Treatment , data = metadata)

# Check and clean column names of the design matrix
colnames(design_group) <- make.names(colnames(design_group))

# Define the contrast matrix using only estimable coefficients
contrast.matrix <- makeContrasts(
  LPSvsShamInBtg = GenotypeBtg.TreatmentLPS - GenotypeBtg.TreatmentSham,
  levels = design_group
)

# Fit the model using the expression data and design matrix
fit <- lmFit(expression_data, design_group)

# Apply the contrast matrix
fit2 <- contrasts.fit(fit, contrast.matrix)

# Apply empirical Bayes moderation
fit2 <- eBayes(fit2)

# Extract results for each contrast and add the Accession number (protein IDs) to the results
results_LPSvsShamInBtg <- topTable(fit2, coef="LPSvsShamInBtg", adjust="BH",number = Inf) %>%
  rownames_to_column(var = "Accession")

# Save results as CSV files with Accession number included
write.csv(results_LPSvsShamInBtg, file="results_LPSvsShamInBtg.csv", row.names=FALSE)
#######

### refrence BKSham

# Load required libraries
library(DEqMS)
library(limma)
library(tidyverse)

# Load the metadata file
metadata <- read.csv("meta_data2.csv", row.names = 1)

# Convert Genotype and Treatment into factors with appropriate reference levels
metadata$Genotype <- factor(metadata$Genotype, levels = c("BK", "WT", "Btg"))  # BK as reference
metadata$Treatment <- factor(metadata$Treatment, levels = c("Sham", "ShamTB", "LPSTB", "LPS"))  # Sham as reference

# Create the design matrix based on Genotype and Treatment
design <- model.matrix(~ Genotype + Treatment, data = metadata)

# Load your normalized expression data (rows are genes/proteins, columns are samples)
expression_data <- read.csv("Mito_norm_total.csv", row.names = 1)

# The row names of expression_data are your Accession numbers (protein IDs)
accession_numbers <- rownames(expression_data)

# Fit the linear model using the expression data and design matrix
fit <- lmFit(expression_data, design)

# Create the design matrix for Genotype and Treatment interactions
design_group <- model.matrix(~ Genotype:Treatment , data = metadata)

# Check and clean column names of the design matrix
colnames(design_group) <- make.names(colnames(design_group))

# Define the contrast matrix using only estimable coefficients
contrast.matrix <- makeContrasts(
  LPSvsShamInBK = GenotypeBK.TreatmentLPS - GenotypeBK.TreatmentSham,
  levels = design_group
)

# Fit the model using the expression data and design matrix
fit <- lmFit(expression_data, design_group)

# Apply the contrast matrix
fit2 <- contrasts.fit(fit, contrast.matrix)

# Apply empirical Bayes moderation
fit2 <- eBayes(fit2)

# Extract results for each contrast and add the Accession number (protein IDs) to the results
results_LPSvsShamInBK <- topTable(fit2, coef="LPSvsShamInBK", adjust="BH",number = Inf) %>%
  rownames_to_column(var = "Accession")

# Save results as CSV files with Accession number included
write.csv(results_LPSvsShamInBK, file="results_LPSvsShamInBK.csv", row.names=FALSE)
#######


###WTTBSHAM refrence
# Load required libraries
library(DEqMS)
library(limma)
library(tidyverse)

# Load the metadata file
metadata <- read.csv("meta_data2.csv", row.names = 1)

# Convert Genotype and Treatment into factors with appropriate reference levels
#metadata$Genotype <- factor(metadata$Genotype, levels = c("WT", "BK", "Btg"))  # WTTB as reference
#metadata$Treatment <- factor(metadata$Treatment, levels = c("TBSham", "Sham", "LPSTB", "LPS"))  # Sham as reference

metadata$group <- factor(metadata$group, levels = c("WTShamTB","WTLPSTB","WTSham","WTLPS",
                                                    "BtgSham","BtgLPS","BKSham","BKLPS"))
 
# Create the design matrix based on Genotype and Treatment
design <- model.matrix(~ group, data = metadata)

# Load your normalized expression data (rows are genes/proteins, columns are samples)
expression_data <- read.csv("Mito_norm_total.csv", row.names = 1)

# The row names of expression_data are your Accession numbers (protein IDs)
accession_numbers <- rownames(expression_data)

# Fit the linear model using the expression data and design matrix
fit <- lmFit(expression_data, design)

# Create the design matrix for Genotype and Treatment interactions
design_group <- model.matrix(~ group , data = metadata)

# Check and clean column names of the design matrix
colnames(design_group) <- make.names(colnames(design_group))

# Define the contrast matrix using only estimable coefficients
contrast.matrix <- makeContrasts(
  TBLPSvsTBShamInWT = groupWTLPSTB,
  TBLPSvsTWSham=groupWTLPSTB,
  levels = design_group
)

# Fit the model using the expression data and design matrix
fit <- lmFit(expression_data, design_group)

# Apply the contrast matrix
fit2 <- contrasts.fit(fit, contrast.matrix)

# Apply empirical Bayes moderation
fit2 <- eBayes(fit2)

# Extract results for each contrast and add the Accession number (protein IDs) to the results
results_TBLPSvsTBShamInWT <- topTable(fit2, coef="TBLPSvsTBShamInWT", adjust="BH",number = Inf) %>%
  rownames_to_column(var = "Accession")

# Save results as CSV files with Accession number included
write.csv(results_TBLPSvsTBShamInWT, file="results_TBLPSvsTBShamInWT.csv", row.names=FALSE)
#######
#########


###WTSHAM refrence
# Load required libraries
library(DEqMS)
library(limma)
library(tidyverse)

# Load the metadata file
metadata <- read.csv("meta_data2.csv", row.names = 1)

# Convert Genotype and Treatment into factors with appropriate reference levels
#metadata$Genotype <- factor(metadata$Genotype, levels = c("WT", "BK", "Btg"))  # WTTB as reference
#metadata$Treatment <- factor(metadata$Treatment, levels = c("TBSham", "Sham", "LPSTB", "LPS"))  # Sham as reference

metadata$group <- factor(metadata$group, levels = c("WTShamTB","WTLPSTB","WTSham","WTLPS",
                                                    "BtgSham","BtgLPS","BKSham","BKLPS"))

# Create the design matrix based on Genotype and Treatment
design <- model.matrix(~ group, data = metadata)

# Load your normalized expression data (rows are genes/proteins, columns are samples)
expression_data <- read.csv("Mito_norm_total.csv", row.names = 1)

# The row names of expression_data are your Accession numbers (protein IDs)
accession_numbers <- rownames(expression_data)

# Fit the linear model using the expression data and design matrix
fit <- lmFit(expression_data, design)

# Create the design matrix for Genotype and Treatment interactions
design_group <- model.matrix(~ group , data = metadata)

# Check and clean column names of the design matrix
colnames(design_group) <- make.names(colnames(design_group))

# Define the contrast matrix using only estimable coefficients
contrast.matrix <- makeContrasts(
  TBLPSvsTWSham=groupWTLPSTB,
  levels = design_group
)

# Fit the model using the expression data and design matrix
fit <- lmFit(expression_data, design_group)

# Apply the contrast matrix
fit2 <- contrasts.fit(fit, contrast.matrix)

# Apply empirical Bayes moderation
fit2 <- eBayes(fit2)

# Extract results for each contrast and add the Accession number (protein IDs) to the results
results_TBLPSvsWTShamInWT <- topTable(fit2, coef="WTTBLPSvsWTShamInWT", adjust="BH",number = Inf) %>%
  rownames_to_column(var = "Accession")

# Save results as CSV files with Accession number included
write.csv(results_TBLPSvsShamInWT, file="results_TBLPSvsShamInWT.csv", row.names=FALSE)
#######



####new ref WTSHAM refrence

# Load required libraries
library(DEqMS)
library(limma)
library(tidyverse)

# Step 1: Load Metadata
metadata <- read.csv("meta_data2.csv", row.names = 1)

# Step 2: Convert columns into factors with appropriate reference levels
metadata$group <- factor(metadata$group, levels = c("WTShamTB", "WTLPSTB", "WTSham", "WTLPS",
                                                    "BtgSham", "BtgLPS", "BKSham", "BKLPS"))

# Step 3: Create the design matrix based on the group variable
design <- model.matrix(~ 0 + group, data = metadata)  # 0 removes the intercept to avoid reference level confusion
colnames(design) <- make.names(colnames(design))  # Clean column names to be syntactically valid

# Step 4: Load normalized expression data (rows = proteins, columns = samples)
expression_data <- read.csv("Mito_norm_total.csv", row.names = 1)

# Step 5: Fit the linear model using the expression data and the design matrix
fit <- lmFit(expression_data, design)

# Step 6: Define the contrast matrix for the comparison of interest
contrast.matrix <- makeContrasts(
  TBLPSvsWTSham = groupWTLPSTB - groupWTSham,
  levels = design
)

# Step 7: Apply the contrast matrix and empirical Bayes moderation
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Step 8: Extract the results for the contrast and add the Accession number (protein IDs)
results_TBLPSvsWTSham <- topTable(fit2, coef = "TBLPSvsWTSham", adjust = "BH", number = Inf) %>%
  rownames_to_column(var = "Accession")

# Step 9: Save the results to a CSV file
write.csv(results_TBLPSvsWTSham, file = "results_TBLPSvsWTSham.csv", row.names = FALSE)

# Print a summary of the results
head(results_TBLPSvsWTSham)


#### refrence WTLPS

# Load required libraries
library(DEqMS)
library(limma)
library(tidyverse)

# Step 1: Load Metadata
metadata <- read.csv("meta_data2.csv", row.names = 1)

# Step 2: Convert columns into factors with WTLPS as the reference level
metadata$group <- factor(metadata$group, levels = c("WTLPS", "BtgLPS", "BKLPS", "WTLPSTB", "WTSham", "WTShamTB", "BtgSham", "BKSham"))

# Step 3: Create the design matrix based on the group variable
design <- model.matrix(~ 0 + group, data = metadata)  # 0 removes the intercept to avoid reference level confusion
colnames(design) <- make.names(colnames(design))  # Clean column names

# Step 4: Load normalized expression data (rows = proteins, columns = samples)
expression_data <- read.csv("Mito_norm_total.csv", row.names = 1)

# Step 5: Fit the linear model using the expression data and the design matrix
fit <- lmFit(expression_data, design)

# Step 6: Define the contrast matrix for the desired comparisons
contrast.matrix <- makeContrasts(
  BtgLPSvsWTLPS = groupBtgLPS - groupWTLPS,
  BKLPSvsWTLPS = groupBKLPS - groupWTLPS,
  WTTBLPSvsWTLPS = groupWTLPSTB - groupWTLPS,
  levels = design
)

# Step 7: Apply the contrast matrix and empirical Bayes moderation
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Step 8: Extract the results for each contrast and add the Accession number (protein IDs)
results_BtgLPSvsWTLPS <- topTable(fit2, coef = "BtgLPSvsWTLPS", adjust = "BH", number = Inf) %>%
  rownames_to_column(var = "Accession")

results_BKLPSvsWTLPS <- topTable(fit2, coef = "BKLPSvsWTLPS", adjust = "BH", number = Inf) %>%
  rownames_to_column(var = "Accession")

results_WTTBLPSvsWTLPS <- topTable(fit2, coef = "WTTBLPSvsWTLPS", adjust = "BH", number = Inf) %>%
  rownames_to_column(var = "Accession")

# Step 9: Save the results to CSV files
write.csv(results_BtgLPSvsWTLPS, file = "results_BtgLPSvsWTLPS.csv", row.names = FALSE)
write.csv(results_BKLPSvsWTLPS, file = "results_BKLPSvsWTLPS.csv", row.names = FALSE)
write.csv(results_WTTBLPSvsWTLPS, file = "results_WTTBLPSvsWTLPS.csv", row.names = FALSE)

# Print a summary of the results
head(results_BtgLPSvsWTLPS)
head(results_BKLPSvsWTLPS)
head(results_WTTBLPSvsWTLPS)

