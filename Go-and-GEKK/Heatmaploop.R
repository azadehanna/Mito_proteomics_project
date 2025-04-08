# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readr)
library(ggplot2)
library(tibble)

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

# Define colors for `sub_mito` annotation
sub_mito_colors <- c("MIM" = "yellow", "IMS" = "lightblue", "Matrix" = "lightgreen", "MOM" = "pink")

# Loop through all files to generate separate heatmaps
for (csv_file in names(file_list)) {
  output_prefix <- file_list[[csv_file]]
  
  # Load the dataset
  df <- read_csv(csv_file)
  
  # Ensure 'logFC' is numeric and remove NaN/Inf values
  df <- df %>% mutate(logFC = as.numeric(logFC)) %>% filter(!is.na(sub_mito), !is.na(logFC))
  
 
  # Select up to 10 most upregulated and 10 most downregulated proteins (or fewer if unavailable)
  num_up <- min(10, sum(df$logFC > 0, na.rm = TRUE))  # Count available upregulated proteins
  num_down <- min(10, sum(df$logFC < 0, na.rm = TRUE))  # Count available downregulated proteins
  
  df_up <- df %>% filter(logFC > 0) %>% arrange(desc(logFC)) %>% slice_head(n = num_up)
  df_down <- df %>% filter(logFC < 0) %>% arrange(logFC) %>% slice_head(n = num_down)
  
  # Combine selected proteins
  df_selected <- bind_rows(df_up, df_down)
  
  
  # Ensure unique row names for heatmap
  df_selected <- df_selected %>% mutate(Protein_Name = make.unique(as.character(Protein_Name)))
  
  # Prepare heatmap data (logFC values)
  heatmap_data <- df_selected %>% select(Protein_Name, logFC) %>% column_to_rownames(var = "Protein_Name")
  
  # Convert `sub_mito` into a correctly formatted annotation
  df_selected$sub_mito[!(df_selected$sub_mito %in% names(sub_mito_colors))] <- "Unknown"
  sub_mito_annot <- df_selected %>% select(Protein_Name, sub_mito) %>% column_to_rownames(var = "Protein_Name")
  
  # Create annotation
  row_ha <- rowAnnotation(
    sub_mito = factor(sub_mito_annot$sub_mito, levels = names(sub_mito_colors)),
    col = list(sub_mito = sub_mito_colors),
    annotation_legend_param = list(title = "sub_mito"),
    show_annotation_name = TRUE
  )
  
  # Generate and save individual heatmap as PDF
  pdf(file = paste0(output_prefix, "_heatmap.pdf"), width = 6, height = 8)
  draw(
    Heatmap(
      as.matrix(heatmap_data),
      name = "logFC",
      col = colorRamp2(c(min(heatmap_data, na.rm = TRUE), 0, max(heatmap_data, na.rm = TRUE)), c("blue", "white", "red")),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = TRUE,
      show_column_names = FALSE,
      left_annotation = row_ha,
      row_split = factor(sub_mito_annot$sub_mito, levels = names(sub_mito_colors))
    )
  )
  dev.off()
}


