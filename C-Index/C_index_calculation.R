# 1. Extract C-index directly from the model summary stats
# Dxy is Somers' D rank correlation
dxy <- f$stats["Dxy"]

# 2. Convert Dxy to C-index using the standard formula: C = 0.5 + Dxy/2
c_index_final <- as.numeric(0.5 + (dxy / 2))

# 3. Print the result in a clean format
cat("\n================================================\n")
cat("      EXTERNAL VALIDATION ACCURACY METRICS      \n")
cat("================================================\n")
cat("Harrell's C-index: ", round(c_index_final, 4), "\n")
cat("================================================\n")
