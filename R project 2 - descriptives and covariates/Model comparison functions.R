

evaluate_deseq_models <- function(count_data, col_data, design_formulas, alpha = 0.05, logfc_cutoff = 2) {
  results_list <- list()
  
  for (i in seq_along(design_formulas)) {
    formula <- design_formulas[[i]]
    model_name <- paste0("model_", i, "_", gsub("[~+ ]", "_", deparse(formula)))
    
    dds <- DESeqDataSetFromMatrix(
      countData = count_data,
      colData = col_data,
      design = formula
    )
    
    dds <- DESeq(dds)
    res <- results(dds, alpha = alpha)
    
    # Create dataframe from results
    res_df <- as.data.frame(res) %>%
      filter(!is.na(padj), !is.na(log2FoldChange)) %>%
      mutate(
        significant = padj < alpha,
        high_logfc = abs(log2FoldChange) > logfc_cutoff,
        direction = ifelse(log2FoldChange > 0, "Up", "Down"),
        category = paste0(
          ifelse(significant, "Sig", "NotSig"), "_",
          ifelse(high_logfc, "HighFC", "LowFC"), "_",
          direction
        )
      )
    
    # Count occurrences in each category
    category_counts <- table(res_df$category) %>% as.data.frame()
    colnames(category_counts) <- c("category", "count")
    category_matrix <- tidyr::pivot_wider(category_counts, names_from = category, values_from = count, values_fill = 0)
    category_matrix$model <- model_name
    
    # Main metrics
    sig_genes <- sum(res_df$significant)
    sig_highfc_genes <- sum(res_df$significant & res_df$high_logfc)
    
    results_list[[model_name]] <- list(
      dds = dds,
      results = res,
      significant_gene_count = sig_genes,
      significant_highFC_count = sig_highfc_genes,
      breakdown = category_matrix
    )
  }
  
  return(results_list)
}


get_formula_combinations <- function(base, covariates) {
  n <- length(covariates)
  formula_list <- list()
  
  # Loop over all subset sizes (1 to n)
  for (k in 0:n) {
    combos <- combn(covariates, k, simplify = FALSE)
    for (combo in combos) {
      rhs <- paste(c(base, combo), collapse = "+")
      formula_list[[length(formula_list) + 1]] <- as.formula(paste("~", rhs))
    }
  }
  
  return(formula_list)
}



make_volcano_plot <- function(res, model_name, xlim_v, ylim_v) {
  res_df <- as.data.frame(res)
  res_df <- res_df %>%
    mutate(
      log10padj = -log10(padj),
      sig = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "Not Significant")
    )

  ggplot(res_df, aes(x = log2FoldChange, y = log10padj)) +
    geom_point(aes(color = sig), alpha = 0.6) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
    labs(title = model_name, x = "Log2 Fold Change", y = "-log10(padj)") +
    theme_minimal() +
    theme(legend.position = "none") +
    xlim(xlim_v) +
    ylim(ylim_v)
}


plot_TE_bar_horizontal <- function(df, group_var, title_text) {
  df %>%
    group_by(.data[[group_var]], direction) %>%
    summarise(count = n(), .groups = "drop") %>%
    ggplot(aes(x = count, y = .data[[group_var]], fill = direction)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = title_text, y = group_var, x = "Count", fill = "Direction") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8)) +
    scale_fill_manual(values = c("Upregulated" = "yellow", "Downregulated" = "#b51963"))
}





