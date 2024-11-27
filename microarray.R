# Setup and package loading

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(pheatmap)
  library(ggplot2)
  library(ggrepel)
})

# function to read and process data
process_microarray_data <- function(file_path) {
  raw_data <- read.table(file_path, header = T, row.names = 1, stringsAsFactors = FALSE)
  log2_data <- log2(raw_data)
  data_norm <- as.matrix(normalizeQuantiles(log2_data))
  return(data_norm)
}

# function for differential expression analysis
run_differential_analysis <- function(norm_data, design_matrix) {
  fit <- limFit(norm_data, design_matrix)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = 2, adjust.method = 'BH', sort.by = 'p', number=Inf)
  return(results)
}

# Enhanced volcano plot function
create_volcano_plot <- function(de_results, fc_threshold = 1, p_threshold = 0.05) {
  de_results <- de_results %>%
    mutate(
      significance = case_when(
        abs(logFC) >= fc_threshold & adj.P.Val < p_threshold ~ "Significant",
        TRUE ~ "Not Significant"
      ),
      direction = case_when(
        logFC >= fc_threshold & adj.P.Val < p_threshold ~ "Up",
        logFC <= -fc_threshold & adj.P.Val < p_threshold ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  ggplot(de_results, aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
    geom_point(size = 2, alpha = 0.6) +
    scale_color_manual(values = c("Down" = "blue", "Up" = "red", "NS" = "grey")) +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed") +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed") +
    theme_minimal() +
    labs(
      title = "Volcano Plot of Differential Expression",
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    )
}

# Function to analyze GO terms
analyze_go_terms <- function(normalized_data, gene_list, term_name) {
  genes_in_term <- unique(unlist(gene_list))
  term_data <- normalized_data[match(genes_in_term, rownames(normalized_data)), ]
  
  pheatmap(term_data,
           scale = "row",
           cluster_cols = FALSE,
           border_color = NA,
           main = paste("Heatmap:", term_name),
           show_rownames = TRUE,
           fontsize_row = 8)
}

# Main analysis pipeline
main_analysis <- function(raw_data_path, go_terms_path) {
  # Process data
  norm_data <- process_microarray_data(raw_data_path)
  
  # Create design matrix (modify as needed)
  design_mat <- cbind(1, c(0, 0, 1, 1))
  
  # Run differential expression
  de_results <- run_differential_analysis(norm_data, design_mat)
  
  # Create and save volcano plot
  volcano_plot <- create_volcano_plot(de_results)
  ggsave("volcano_plot.png", volcano_plot, dpi = 300, width = 8, height = 6)
  
  # Save results
  write.csv(de_results, "differential_expression_results.csv")
  write.csv(norm_data, "normalized_data.csv")
  
  return(list(
    normalized_data = norm_data,
    de_results = de_results
  ))
}
