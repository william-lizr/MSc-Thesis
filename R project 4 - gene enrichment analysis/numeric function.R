# Function to safely convert all numeric-like columns into numeric
convert_to_numeric <- function(df) {
  df[] <- lapply(df, function(col) {
    # If column is factor, convert to character first
    if (is.factor(col)) {
      col <- as.character(col)
    }
    
    # If column is character, check if it can be converted to numeric
    if (is.character(col)) {
      suppressWarnings(num_col <- as.numeric(col))
      
      # If conversion succeeds for at least one non-NA value, replace
      if (any(!is.na(num_col))) {
        return(num_col)
      } else {
        return(col) # leave as character
      }
    }
    
    # If already numeric or integer, leave unchanged
    if (is.numeric(col) || is.integer(col)) {
      return(as.numeric(col))
    }
    
    # Otherwise, return original
    return(col)
  })
  
  return(df)
}
