te_gene_linear_modeling <- function(te_gene_pairs, 
                                    gene_expression_matrix, 
                                    te_expression_matrix, 
                                    participant_data,
                                    data_type = "normalized", # "counts", "normalized", or "auto"
                                    count_distribution = "negative_binomial", # "poisson" or "negative_binomial"
                                    covariates = NULL,
                                    p_value_threshold = 0.05,
                                    adjust_method = "BH",
                                    library_size_normalization = TRUE) {
  
  # Load required libraries
  require(dplyr)
  require(broom)
  if (data_type == "counts" || data_type == "auto") {
    require(MASS)  # For negative binomial GLM
    require(edgeR) # For count data utilities
  }
  
  # Input validation
  if (!is.data.frame(te_gene_pairs)) {
    stop("te_gene_pairs must be a data frame")
  }
  if (!all(c("gene_id", "te_id") %in% colnames(te_gene_pairs))) {
    stop("te_gene_pairs must contain 'gene_id' and 'te_id' columns")
  }
  if (!is.matrix(gene_expression_matrix) && !is.data.frame(gene_expression_matrix)) {
    stop("gene_expression_matrix must be a matrix or data frame")
  }
  if (!is.matrix(te_expression_matrix) && !is.data.frame(te_expression_matrix)) {
    stop("te_expression_matrix must be a matrix or data frame")
  }
  if (!is.data.frame(participant_data)) {
    stop("participant_data must be a data frame")
  }
  
  # Ensure matrices have sample names
  if (is.null(colnames(gene_expression_matrix))) {
    stop("gene_expression_matrix must have column names (sample IDs)")
  }
  if (is.null(colnames(te_expression_matrix))) {
    stop("te_expression_matrix must have column names (sample IDs)")
  }
  
  # Auto-detect data type
  if (data_type == "auto") {
    gene_sample <- as.numeric(gene_expression_matrix[1:min(100, nrow(gene_expression_matrix)), 1:min(10, ncol(gene_expression_matrix))])
    te_sample <- as.numeric(te_expression_matrix[1:min(100, nrow(te_expression_matrix)), 1:min(10, ncol(te_expression_matrix))])
    
    all_integers <- all(gene_sample == round(gene_sample), na.rm = TRUE) && 
      all(te_sample == round(te_sample), na.rm = TRUE)
    all_non_negative <- all(gene_sample >= 0, na.rm = TRUE) && 
      all(te_sample >= 0, na.rm = TRUE)
    
    if (all_integers && all_non_negative) {
      data_type <- "counts"
      cat("Auto-detected count data. Using GLM approach.\n")
    } else {
      data_type <- "normalized"
      cat("Auto-detected normalized data. Using linear model approach.\n")
    }
  }
  
  # Calculate library sizes for normalization (counts only)
  if (data_type == "counts" && library_size_normalization) {
    gene_lib_sizes <- colSums(gene_expression_matrix, na.rm = TRUE)
    te_lib_sizes <- colSums(te_expression_matrix, na.rm = TRUE)
    cat("Calculated library sizes for normalization\n")
  }
  
  # Match common samples
  gene_samples <- colnames(gene_expression_matrix)
  te_samples <- colnames(te_expression_matrix)
  participant_samples <- rownames(participant_data)
  
  if (is.null(participant_samples)) {
    participant_samples <- participant_data[[1]]
    rownames(participant_data) <- participant_samples
  }
  
  common_samples <- Reduce(intersect, list(gene_samples, te_samples, participant_samples))
  if (length(common_samples) == 0) {
    stop("No common samples found across gene expression, TE expression, and participant data")
  }
  
  cat("Found", length(common_samples), "common samples for analysis\n")
  
  # Subset
  gene_expr_subset <- gene_expression_matrix[, common_samples, drop = FALSE]
  te_expr_subset <- te_expression_matrix[, common_samples, drop = FALSE]
  participant_subset <- participant_data[common_samples, , drop = FALSE]
  
  # Results table
  results <- data.frame(
    gene_id = character(),
    te_id = character(),
    te_coefficient = numeric(),
    te_std_error = numeric(),
    te_z_or_t_value = numeric(),
    te_p_value = numeric(),
    te_adj_p_value = numeric(),
    model_deviance = numeric(),
    model_aic = numeric(),
    model_df_residual = numeric(),
    overdispersion_param = numeric(),
    n_samples = numeric(),
    model_type = character(),
    stringsAsFactors = FALSE
  )
  
  cat("Testing", nrow(te_gene_pairs), "TE-gene pairs...\n")
  pb <- txtProgressBar(min = 0, max = nrow(te_gene_pairs), style = 3)
  
  # Loop through TE-gene pairs
  for (i in 1:nrow(te_gene_pairs)) {
    setTxtProgressBar(pb, i)
    current_gene <- te_gene_pairs$gene_id[i]
    current_te <- te_gene_pairs$te_id[i]
    
    if (!current_gene %in% rownames(gene_expr_subset)) {
      warning(paste("Gene", current_gene, "not found in gene expression matrix"))
      next
    }
    if (!current_te %in% rownames(te_expr_subset)) {
      warning(paste("TE", current_te, "not found in TE expression matrix"))
      next
    }
    
    # Extract expression
    gene_expr <- as.numeric(gene_expr_subset[current_gene, ])
    te_expr <- as.numeric(te_expr_subset[current_te, ])
    
    model_data <- data.frame(
      sample_id = common_samples,
      gene_expression = gene_expr,
      te_expression = te_expr,
      participant_subset,
      stringsAsFactors = FALSE
    )
    
    # Add lib sizes if counts
    if (data_type == "counts" && library_size_normalization) {
      model_data$gene_lib_size <- gene_lib_sizes[common_samples]
    }
    
    model_data <- model_data[complete.cases(model_data), ]
    if (nrow(model_data) < 3) {
      warning(paste("Insufficient data for gene", current_gene, "and TE", current_te))
      next
    }
    
    # Build formula
    formula_terms <- "gene_expression ~ te_expression"
    if (!is.null(covariates)) {
      available_covariates <- covariates[covariates %in% colnames(model_data)]
      if (length(available_covariates) > 0) {
        formula_terms <- paste(formula_terms, "+", paste(available_covariates, collapse = " + "))
      }
    }
    if (data_type == "counts" && library_size_normalization) {
      formula_terms <- paste0(formula_terms, " + offset(log(gene_lib_size))")
    }
    model_formula <- as.formula(formula_terms)
    
    # Fit models
    tryCatch({
      if (data_type == "counts") {
        if (count_distribution == "poisson") {
          glm_model <- glm(model_formula, data = model_data, family = poisson(link = "log"))
          model_summary <- summary(glm_model)
          overdispersion <- NA
          model_type <- "Poisson GLM"
          
        } else if (count_distribution == "negative_binomial") {
          nb_model <- glm.nb(model_formula, data = model_data)
          glm_model <- nb_model
          model_summary <- summary(nb_model)
          overdispersion <- nb_model$theta
          model_type <- "Negative Binomial GLM"
        }
        
        te_coef_info <- model_summary$coefficients["te_expression", ]
        
        results <- rbind(results, data.frame(
          gene_id = current_gene,
          te_id = current_te,
          te_coefficient = te_coef_info["Estimate"],
          te_std_error = te_coef_info["Std. Error"],
          te_z_or_t_value = ifelse(model_type == "Poisson GLM", te_coef_info["z value"], te_coef_info["z value"]),
          te_p_value = te_coef_info[length(te_coef_info)],  # last element is p-value
          te_adj_p_value = NA,
          model_deviance = glm_model$deviance,
          model_aic = glm_model$aic,
          model_df_residual = glm_model$df.residual,
          overdispersion_param = ifelse(is.na(overdispersion), NA, overdispersion),
          n_samples = nrow(model_data),
          model_type = model_type,
          stringsAsFactors = FALSE
        ))
        
      } else {
        lm_model <- lm(model_formula, data = model_data)
        model_summary <- summary(lm_model)
        te_coef_info <- model_summary$coefficients["te_expression", ]
        
        results <- rbind(results, data.frame(
          gene_id = current_gene,
          te_id = current_te,
          te_coefficient = te_coef_info["Estimate"],
          te_std_error = te_coef_info["Std. Error"],
          te_z_or_t_value = te_coef_info["t value"],
          te_p_value = te_coef_info["Pr(>|t|)"],
          te_adj_p_value = NA,
          model_deviance = deviance(lm_model),
          model_aic = AIC(lm_model),
          model_df_residual = lm_model$df.residual,
          overdispersion_param = NA,
          n_samples = nrow(model_data),
          model_type = "Linear Model",
          stringsAsFactors = FALSE
        ))
      }
    }, error = function(e) {
      warning(paste("Error fitting model for gene", current_gene, "and TE", current_te, ":", e$message))
    })
  }
  
  close(pb)
  
  # Multiple testing correction
  if (nrow(results) > 0) {
    results$te_adj_p_value <- p.adjust(results$te_p_value, method = adjust_method)
    results$te_significant <- results$te_adj_p_value < p_value_threshold
    results$te_effect_direction <- ifelse(results$te_coefficient > 0, "positive", "negative")
    results <- results[order(results$te_adj_p_value), ]
    
    cat("\n--- Analysis Summary ---\n")
    cat("Data type used:", unique(results$model_type), "\n")
    cat("Total TE-gene pairs tested:", nrow(results), "\n")
    cat("Significant associations (adj. p <", p_value_threshold, "):", 
        sum(results$te_significant), "\n")
    cat("Positive effects:", sum(results$te_significant & results$te_effect_direction == "positive"), "\n")
    cat("Negative effects:", sum(results$te_significant & results$te_effect_direction == "negative"), "\n")
    cat("Multiple testing correction method:", adjust_method, "\n")
    if (data_type == "counts") {
      cat("Distribution used:", count_distribution, "\n")
      cat("Library size normalization:", library_size_normalization, "\n")
    }
  } else {
    cat("\nNo successful model fits. Please check your input data.\n")
  }
  
  return(results)
}
