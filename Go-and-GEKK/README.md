Functional Enrichment Analysis Using clusterProfiler
This project automates Gene Ontology (GO) and KEGG enrichment analysis for multiple sets of differentially expressed genes using the clusterProfiler suite in R. The pipeline includes statistical filtering, ID conversion, enrichment testing, and visualization.

-Required R Packages
The following R packages must be installed before running the script:
library(DOSE)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggrepel)
library(gridExtra)
-Input Files
Each input file should be a CSV with at least the following columns:
Gene_Name: gene symbols
logFC: log fold change
P.Value: raw p-value
The script reads a list of filtered differential expression results from multiple comparisons:
"filtered_results_XXvsYY.csv"
These files are listed in the file_list variable in the script.

-Analysis Workflow
1. Filter Significant Genes
Criteria: |logFC| ≥ 0.58 and P.Value ≤ 0.05

2. Convert Gene Symbols to Entrez IDs
Uses bitr() from clusterProfiler with org.Mm.eg.db as the annotation database.

3. GO Enrichment
Enrichment is performed for:
Biological Process (BP)
Molecular Function (MF)
Results are saved as:
*_go_enrich_BP.csv
*_go_enrich_MF.csv

4. KEGG Pathway Enrichment
Enrichment is performed for Mus musculus (organism = 'mmu')
If significant KEGG pathways are found:
Top 5 pathways are plotted
Saved as: *_KEGG_top5_barplot.pdf

5. Plotting
GO term enrichment plots saved as:
*_go_enrich_BP_barplot.pdf
*_go_enrich_MF_barplot.pdf

Bar plots use custom theming and display the top 10 terms by significance.

- Output Files
For each input file, the script outputs:

Type	Output File Example
GO BP Results	BKLPSvsWTLPS_go_enrich_BP.csv
GO MF Results	BKLPSvsWTLPS_go_enrich_MF.csv
GO BP Barplot	BKLPSvsWTLPS_go_enrich_BP_barplot.pdf
GO MF Barplot	BKLPSvsWTLPS_go_enrich_MF_barplot.pdf
KEGG Pathway Plot	BKLPSvsWTLPS_KEGG_top5_barplot.pdf
-Automation
The script loops over the predefined file_list:
Each file is processed using the process_gene_data() function, and results are saved automatically.

-Author
Developed by Azadeh Nikouee for pathway enrichment and gene ontology analysis of proteomics-based differential gene expression results in mice.

