Heatmap Generator for Differential Expression Analysis
This R script generates individual heatmaps from multiple differential expression result files. It highlights the top upregulated and downregulated mitochondrial proteins based on logFC values, annotated by their mitochondrial sub-localization.

-File Overview
Script: Heatmaploop.R
Input: CSV files with differential expression results with add information from mitosublocalization
Output: PDF files with annotated heatmaps

-Input Requirements
Each input CSV file should contain at least the following columns:
Protein_Name — Name of the protein.
logFC — Log fold-change value (numeric).
sub_mito — Mitochondrial sub-compartment classification (MIM, IMS, Matrix, MOM).

-Heatmap Description
For each file:
Selects up to 10 most upregulated and 10 most downregulated proteins.
Generates a heatmap of logFC values.
Adds row annotations for sub_mito categories using custom colors.
Saves output as a PDF file (one per input CSV).
Unknown compartments are labeled as "Unknown".
-Required R Packages  R version 4.4.2 
ComplexHeatmap version  2.22.0
Circlize  version 0.4.16  

dplyr version  1.1.4
readr version  2.1.5
ggplot2 version  3.5.1
tibble version3.2.1
-How to Run

Place all required .csv files in the working directory.
Run the script
Rscript Heatmaploop.R
PDF files named {contrast_name}_heatmap.pdf will be generated for each input.
-Input Files Defined in Script
The script processes the following files:
All filtered_results_XXvsYY.csv
-Notes
The script automatically handles non-numeric or missing values in logFC.
Each heatmap is not clustered to preserve comparison-specific order.
Modify the file_list in the script to add/remove comparisons.

