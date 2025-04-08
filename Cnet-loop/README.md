GO Enrichment and Cnetplot Visualization in Mouse Gene Data
This R script (Cnetloop.R) automates Gene Ontology (GO) enrichment analysis for the Biological Process (BP) category and generates cnetplots for significant gene sets using clusterProfiler. It supports batch processing of multiple differential expression result files.

-Required R Packages
Make sure to install and load the following packages before running the script:
library(DOSE)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(igraph)
-Input File Format
Each input CSV file must contain the following columns:
Gene_Name – Gene symbols (mouse)
logFC – Log fold-change
P.Value – P-value for differential expression
Files are provided in the format:
filtered_results_<Comparison>.csv
filtered_results_XXvsYY.csv

- Workflow Overview
For each file in the file_list, the script performs:
1. Filtering Genes
Criteria: |logFC| ≥ 0.58  , P.Value ≤ 0.05

2. Gene ID Conversion
Converts SYMBOL → ENTREZID using org.Mm.eg.db.
3. GO Enrichment Analysis
Conducted for Biological Process (BP) ontology.
Result saved as:
*_go_enrich_BP.csv
4. Cnetplot Generation
Converts enriched gene lists back to gene symbols.
Generates a cnetplot with:
Blue GO terms
Red genes
Layout: igraph::layout_with_kk
Plot saved as PDF:
*_cnetplot_BP.pdf

-Output Files
For each comparison, the script generates:

Output Type	File Example
GO BP Enrichment	LPSvsShamInWT_go_enrich_BP.csv
Cnetplot PDF	LPSvsShamInWT_cnetplot_BP.pdf
- Batch Processing
The file_list object maps all CSV input files to output prefixes. The script loops through each file and automatically applies the full pipeline:

- Author
Developed by Azadeh Nikouee to visualize enriched biological processes from differential expression analyses in mouse proteomics studies.


