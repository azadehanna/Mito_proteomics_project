Differential Expression Analysis of Mitochondrial Proteomics Data Using limma and DEqMS
This repository contains R scripts for performing differential protein expression analysis on normalized mitochondrial proteomics data using the limma and DEqMS packages. The analysis is based on various combinations of genotypes and treatments in mice, with different reference levels to explore the impact of treatment and genetic background.

-Files and Inputs
Mito_norm_total.csv: Normalized expression matrix (rows = protein IDs, columns = samples).
meta_data2.csv: Metadata file describing each sample, including genotype, treatment, and a combined group variable.
-Required R Packages -> R version 4.4.2
Make sure the following R packages are installed and loaded:
library(DEqMS) version: 1.24.0
library(limma) version: 3.62.2
library(tidyverse) version: 2.0.0

-Analysis Overview
The analysis involves several reference-level comparisons using linear modeling and contrasts:

1. WT Sham as Reference
2. Btg Sham as Reference
3. BK Sham as Reference
4. WT ShamTB as Reference (using group variable)
5. WT Sham as Reference (using group variable)
6. WT LPS as Reference (using group variable)

-Outputs:
results_XXvsYY.csv
XX and YY is different based on which situation would like to test.
-Workflow Summary
Load metadata and expression data.
Set factor levels to define reference conditions.
Construct the design matrix for modeling.
Fit linear models using limma::lmFit.
Define contrast matrices for the desired comparisons.
Apply contrasts and empirical Bayes smoothing (eBayes).
Export results to CSV files with protein accession numbers.

-Notes
Ensure consistency between sample order in metadata and expression matrix.
You may adjust reference levels in the metadata depending on the contrast of interest.
Design matrices are recreated for each comparison to accommodate changes in reference levels.

Author
Developed by Azadeh Nikouee for mitochondrial proteomics analysis in genetically and pharmacologically treated mice.


