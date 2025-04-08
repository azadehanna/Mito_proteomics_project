# Load required libraries
library(DOSE)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(igraph)

# Define function to process each dataset
process_gene_data <- function(csv_file, output_prefix) {
  # Load gene data
  gene_data <- read.csv(csv_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Filter significant genes based on logFC and P-value
  gene_data_filtered <- subset(gene_data, abs(logFC) >= 0.58 & P.Value <= 0.05)
  
  # Extract significant gene names
  significant_genes <- gene_data_filtered$Gene_Name
  
  # Convert gene names to Entrez IDs
  gene_ids <- bitr(significant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  # Perform GO enrichment analysis for BP
  go_enrich_BP <- enrichGO(
    gene = gene_ids$ENTREZID, OrgDb = org.Mm.eg.db, keyType = "ENTREZID",
    ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05
  )
  
  # Save GO enrichment results to CSV file
  write.csv(go_enrich_BP@result, paste0(output_prefix, "_go_enrich_BP.csv"), row.names = FALSE)
  
  # Perform and visualize Cnetplot for BP ontology
  gene_symbol_mapping <- bitr(gene_ids$ENTREZID, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
  gene_names_per_term <- lapply(go_enrich_BP@result$geneID, function(ids) {
    genes <- unlist(strsplit(ids, "/"))
    mapped_genes <- gene_symbol_mapping$SYMBOL[match(genes, gene_symbol_mapping$ENTREZID)]
    paste(mapped_genes, collapse = "/")
  })
  
  go_enrich_BP@result$geneID <- gene_names_per_term
  go_enrich_BP@result <- na.omit(go_enrich_BP@result)
  
  pc <- cnetplot(
    go_enrich_BP,
    layout = igraph::layout_with_kk,
    showCategory = 20,
    color_category = "blue",
    size_category = 1,
    color_item = "red",
    size_item = 0.9,
    color_edge = "darkgray",
    size_edge = 0.5,
    node_label = "all"
  )
  
  # Save the Cnetplot
  ggsave(paste0(output_prefix, "_cnetplot_BP.pdf"), plot = pc, width = 10, height = 8)
}

# Define file list
file_list <- list(
  "filtered_results_BKLPSvsWTLPS.csv" = "BKLPSvsWTLPS",
  "filtered_results_BKvsWTInSham.csv" = "BKvsWTInSham",
  "filtered_results_BtgLPSvsWTLPS.csv" = "BtgLPSvsWTLPS",
  "filtered_results_BtgvsWTInSham.csv" = "BtgvsWTInSham",
  "filtered_results_LPSvsShamInBK.csv" = "LPSvsShamInBK",
  "filtered_results_LPSvsShamInBtg.csv" = "LPSvsShamInBtg",
  "filtered_results_LPSvsShamInWT.csv" = "LPSvsShamInWT",
  "filtered_results_ShamTBvsShamInWT.csv" = "ShamTBvsShamInWT",
  "filtered_results_TBLPSvsTBShamInWT.csv" = "TBLPSvsTBShamInWT",
  "filtered_results_TBLPSvsWTSham.csv" = "TBLPSvsWTSham",
  "filtered_results_WTTBLPSvsWTLPS.csv" = "WTTBLPSvsWTLPS"
)

# Loop through all files
for (csv_file in names(file_list)) {
  output_prefix <- file_list[[csv_file]]
  process_gene_data(csv_file, output_prefix)
}




