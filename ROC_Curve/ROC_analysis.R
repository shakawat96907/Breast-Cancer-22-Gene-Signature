# 1. Open TIFF device (1000 DPI)
tiff("ROC_External_Validation_Final.tiff", 
     width = 7, height = 7, units = 'in', res = 1000, compression = "lzw")

# 2. Run ROC calculation (using your merged data)
library(timeROC)
ROC_res <- timeROC(T = final_merged$time,
                   delta = final_merged$status,
                   marker = final_merged$Risk_Score,
                   cause = 1,
                   times = c(365, 1095, 1825), 
                   iid = TRUE)

# 3. Plotting the curves
# We use 'xlab' and 'ylab' inside plot() and set 'title = FALSE' to prevent overlap
plot(ROC_res, time = 365, col = "red", title = FALSE, lwd = 3, lty = 1,
     xlab = "False Positive Rate (1-Specificity)", 
     ylab = "True Positive Rate (Sensitivity)")

# Add 3-Year (Blue)
plot(ROC_res, time = 1095, add = TRUE, col = "blue", lwd = 3, lty = 1)

# Add 5-Year (Dark Green)
plot(ROC_res, time = 1825, add = TRUE, col = "darkgreen", lwd = 3, lty = 1)

# 4. Add Reference line and Legend
abline(a = 0, b = 1, lty = 2, col = "darkgray", lwd = 2)

legend("bottomright", 
       c(paste0("1-Year AUC: ", round(ROC_res$AUC[1], 3)),
         paste0("3-Year AUC: ", round(ROC_res$AUC[2], 3)),
         paste0("5-Year AUC: ", round(ROC_res$AUC[3], 3))),
       col = c("red", "blue", "darkgreen"), 
       lty = 1, lwd = 3, bty = "n", cex = 1.1)

# 5. Add Main Title and Grid
mtext("Time-Dependent ROC Curves (External Validation)", side = 3, line = 1.5, cex = 1.2, font = 2)
grid(col = "lightgray", lty = "dotted")

# 6. Close the device
dev.off()
