# Compare full and reduced models using DESeq2 LRT
compare_deseq_models_lrt <- function(count_data, col_data, full_formula, reduced_formula, 
                                     alpha = 0.05, logfc_cutoff = 2, test_coef = NULL) {
  
  # Create DESeqDataSet with full model
  dds_full <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = col_data,
    design = full_formula
  )
  
  # Run DESeq with LRT comparing full vs reduced model
  dds_lrt <- DESeq(dds_full, test = "LRT", reduced = reduced_formula)
  
  # Get LRT results (tests if full model is significantly better than reduced)
  res_lrt <- results(dds_lrt, alpha = alpha)
  
  # Get Wald test results for effect sizes (if specific coefficient specified)
  res_wald <- NULL
  if (!is.null(test_coef)) {
    dds_wald <- DESeq(dds_full, test = "Wald")
    res_wald <- results(dds_wald, name = test_coef, alpha = alpha)
  } else {
    # Use last coefficient from full model
    dds_wald <- DESeq(dds_full, test = "Wald")
    res_wald <- results(dds_wald, alpha = alpha)
  }
  
  # Create comprehensive results dataframe
  res_df <- as.data.frame(res_lrt) %>%
    rownames_to_column("gene_id") %>%
    filter(!is.na(padj)) %>%
    mutate(
      # LRT significance (model comparison)
      lrt_significant = padj < alpha,
      lrt_log10padj = -log10(padj)
    )
  
  # Add Wald test results for fold changes
  if (!is.null(res_wald)) {
    wald_df <- as.data.frame(res_wald) %>%
      rownames_to_column("gene_id") %>%
      select(gene_id, log2FoldChange, wald_padj = padj) %>%
      filter(!is.na(log2FoldChange))
    
    res_df <- res_df %>%
      left_join(wald_df, by = "gene_id") %>%
      filter(!is.na(log2FoldChange)) %>%
      mutate(
        # Wald test significance (effect size)
        wald_significant = wald_padj < alpha,
        high_logfc = abs(log2FoldChange) > logfc_cutoff,
        direction = ifelse(log2FoldChange > 0, "Up", "Down"),
        
        # Combined classification
        category = case_when(
          lrt_significant & wald_significant & high_logfc ~ paste0("LRT_Wald_HighFC_", direction),
          lrt_significant & wald_significant & !high_logfc ~ paste0("LRT_Wald_LowFC_", direction),
          lrt_significant & !wald_significant ~ "LRT_only",
          !lrt_significant & wald_significant & high_logfc ~ paste0("Wald_only_HighFC_", direction),
          !lrt_significant & wald_significant & !high_logfc ~ paste0("Wald_only_LowFC_", direction),
          TRUE ~ "Not_significant"
        )
      )
  }
  
  # Summary statistics
  lrt_sig_genes <- sum(res_df$lrt_significant, na.rm = TRUE)
  wald_sig_genes <- sum(res_df$wald_significant, na.rm = TRUE)
  both_sig_genes <- sum(res_df$lrt_significant & res_df$wald_significant, na.rm = TRUE)
  both_sig_highfc <- sum(res_df$lrt_significant & res_df$wald_significant & res_df$high_logfc, na.rm = TRUE)
  
  # Category breakdown
  category_counts <- table(res_df$category) %>% 
    as.data.frame() %>%
    setNames(c("category", "count"))
  
  category_matrix <- tidyr::pivot_wider(category_counts, 
                                        names_from = category, 
                                        values_from = count, 
                                        values_fill = 0)
  
  # Model information
  full_terms <- all.vars(full_formula)
  reduced_terms <- all.vars(reduced_formula)
  tested_terms <- setdiff(full_terms, reduced_terms)
  
  results <- list(
    dds_lrt = dds_lrt,
    dds_wald = if(exists("dds_wald")) dds_wald else NULL,
    lrt_results = res_lrt,
    wald_results = res_wald,
    combined_results = res_df,
    
    # Model info
    full_formula = full_formula,
    reduced_formula = reduced_formula,
    tested_terms = tested_terms,
    test_coefficient = test_coef,
    
    # Summary stats
    lrt_significant_count = lrt_sig_genes,
    wald_significant_count = wald_sig_genes,
    both_significant_count = both_sig_genes,
    both_significant_highFC_count = both_sig_highfc,
    
    # Breakdown
    category_breakdown = category_matrix,
    category_counts = category_counts
  )
  
  return(results)
}

# Batch compare multiple model pairs
batch_compare_lrt <- function(count_data, col_data, model_pairs, alpha = 0.05, logfc_cutoff = 2) {
  results_list <- list()
  
  for (i in seq_along(model_pairs)) {
    pair_name <- names(model_pairs)[i]
    if (is.null(pair_name)) pair_name <- paste0("comparison_", i)
    
    full_formula <- model_pairs[[i]]$full
    reduced_formula <- model_pairs[[i]]$reduced
    test_coef <- model_pairs[[i]]$test_coef
    
    cat("Running LRT comparison:", pair_name, "\n")
    cat("Full model:", deparse(full_formula), "\n")
    cat("Reduced model:", deparse(reduced_formula), "\n")
    
    results_list[[pair_name]] <- compare_deseq_models_lrt(
      count_data = count_data,
      col_data = col_data,
      full_formula = full_formula,
      reduced_formula = reduced_formula,
      alpha = alpha,
      logfc_cutoff = logfc_cutoff,
      test_coef = test_coef
    )
    
    cat("LRT significant genes:", results_list[[pair_name]]$lrt_significant_count, "\n")
    cat("Both LRT & Wald significant:", results_list[[pair_name]]$both_significant_count, "\n\n")
  }
  
  return(results_list)
}

# Create model pairs for common comparisons
create_model_pairs <- function(base_vars, test_vars, covariates = NULL) {
  model_pairs <- list()
  
  for (test_var in test_vars) {
    # Full model includes test variable
    full_terms <- c(base_vars, test_var)
    if (!is.null(covariates)) full_terms <- c(full_terms, covariates)
    full_formula <- as.formula(paste("~", paste(full_terms, collapse = " + ")))
    
    # Reduced model excludes test variable
    reduced_terms <- base_vars
    if (!is.null(covariates)) reduced_terms <- c(reduced_terms, covariates)
    reduced_formula <- as.formula(paste("~", paste(reduced_terms, collapse = " + ")))
    
    pair_name <- paste0("test_", test_var)
    model_pairs[[pair_name]] <- list(
      full = full_formula,
      reduced = reduced_formula,
      test_coef = test_var  # coefficient to test in Wald test
    )
  }
  
  return(model_pairs)
}

# Enhanced volcano plot for LRT results
make_lrt_volcano_plot <- function(lrt_comparison, plot_type = "lrt", model_name = "", 
                                  xlim_v = c(-10, 10), ylim_v = c(0, 50)) {
  
  res_df <- lrt_comparison$combined_results
  
  if (plot_type == "lrt") {
    # Plot LRT p-values vs fold changes
    p <- ggplot(res_df, aes(x = log2FoldChange, y = lrt_log10padj)) +
      geom_point(aes(color = lrt_significant), alpha = 0.6) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
      labs(title = paste("LRT Test:", model_name), 
           x = "Log2 Fold Change (Wald)", 
           y = "-log10(LRT padj)",
           subtitle = paste("Tested terms:", paste(lrt_comparison$tested_terms, collapse = ", "))) +
      theme_minimal() +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
      geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "blue")
    
  } else if (plot_type == "wald") {
    # Standard volcano plot with Wald test results
    res_df <- res_df %>%
      mutate(wald_log10padj = -log10(wald_padj))
    
    p <- ggplot(res_df, aes(x = log2FoldChange, y = wald_log10padj)) +
      geom_point(aes(color = wald_significant), alpha = 0.6) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
      labs(title = paste("Wald Test:", model_name), 
           x = "Log2 Fold Change", 
           y = "-log10(Wald padj)") +
      theme_minimal() +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
      geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "blue")
    
  } else {
    # Combined plot showing both tests
    p <- ggplot(res_df, aes(x = log2FoldChange, y = lrt_log10padj)) +
      geom_point(aes(color = category), alpha = 0.7) +
      labs(title = paste("Combined LRT & Wald:", model_name), 
           x = "Log2 Fold Change (Wald)", 
           y = "-log10(LRT padj)",
           color = "Significance") +
      theme_minimal() +
      theme(legend.position = "bottom") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
      geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "blue")
  }
  
  return(p + xlim(xlim_v) + ylim(ylim_v))
}

# Summary plot for LRT comparisons
plot_lrt_summary <- function(lrt_results_list) {
  summary_data <- map_dfr(lrt_results_list, function(x) {
    data.frame(
      lrt_significant = x$lrt_significant_count,
      wald_significant = x$wald_significant_count,
      both_significant = x$both_significant_count,
      both_sig_highfc = x$both_significant_highFC_count,
      tested_terms = paste(x$tested_terms, collapse = ", ")
    )
  }, .id = "comparison")
  
  summary_long <- summary_data %>%
    select(-tested_terms) %>%
    pivot_longer(-comparison, names_to = "test_type", values_to = "count")
  
  ggplot(summary_long, aes(x = comparison, y = count, fill = test_type)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "DESeq2 LRT vs Wald Test Comparison",
         x = "Model Comparison", y = "Significant Genes", fill = "Test Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(type = "qual", palette = "Set2")
}