
parse_splice_annotations <- function(df) {
  # Make sure input has rownames
  if (is.null(rownames(df))) {
    stop("The input data frame must have rownames formatted with '|' separators.")
  }
  
  # Split rownames on "|"
  splice_list <- strsplit(rownames(df), "\\|")
  
  # Check length of each split
  n_fields <- sapply(splice_list, length)
  if (any(n_fields < 6)) {
    warning("Some rownames have fewer than 6 fields. Missing values will be filled with NA.")
  }
  
  safe_extract <- function(lst, idx) {
    sapply(lst, function(x) ifelse(length(x) >= idx, x[idx], NA))
  }
  
  # Extract chromosome and coordinates
  df$chromosome <- safe_extract(splice_list, 1)
  df$start_chord <- safe_extract(splice_list, 2)
  df$end_chord <- safe_extract(splice_list, 3)
  
  # TE annotation
  df$TE_annotation <- safe_extract(splice_list, 4)
  te_anno_list <- strsplit(df$TE_annotation, ":")
  
  df$TE_name   <- safe_extract(te_anno_list, 1)
  df$TE_family <- safe_extract(te_anno_list, 2)
  df$TE_class  <- safe_extract(te_anno_list, 3)
  
  # TE length and strand orientation
  df$te_length <- safe_extract(splice_list, 5)
  df$strand_orientation <- safe_extract(splice_list, 6)
  
  return(df)
}
