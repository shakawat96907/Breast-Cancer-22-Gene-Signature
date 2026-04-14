# 1. Load Libraries
if (!require("rms")) install.packages("rms")
library(rms)
library(survival)

# 2. Load Data correctly
lasso_info <- read.csv("LASSO_signature.csv")
# Load expression data without setting row.names yet
raw_expr <- read.csv("TCGA_BRCA_Full_Expression_1094.csv", check.names = FALSE)
clinical_data <- read.csv("TCGA_BRCA_Clinical_Final_1094.csv")


# 3. Transpose Expression Data (So Barcodes become Rows)
# Assuming 1st column is 'Gene' or 'Symbol'
gene_names <- raw_expr[, 1]
expr_data <- as.data.frame(t(raw_expr[, -1])) # Remove 1st col and Flip
colnames(expr_data) <- gene_names # Set genes as columns

# 4. Standardize Barcodes (12 characters, Uppercase)
# Now rownames(expr_data) contains the TCGA barcodes
expr_clean <- toupper(gsub("\\.", "-", rownames(expr_data)))
expr_data$Match_ID <- substr(expr_clean, 1, 12)
expr_data <- expr_data[!duplicated(expr_data$Match_ID), ]

clin_clean <- toupper(gsub("\\.", "-", clinical_data$bcr_patient_barcode))
clinical_data$Match_ID <- substr(clin_clean, 1, 12)

# 5. Risk Score Calculation
target_genes <- trimws(lasso_info$Gene)
available_genes <- intersect(target_genes, colnames(expr_data))

# Check if genes matched
cat("Genes matched:", length(available_genes), "out of", length(target_genes), "\n")

risk_matrix <- as.matrix(expr_data[, available_genes])
match_idx <- match(available_genes, lasso_info$Gene)
expr_data$Final_Risk_Score <- as.numeric(risk_matrix %*% lasso_info$Coefficient[match_idx])

# 6. Final Merge
final_merged <- merge(clinical_data, 
                      expr_data[, c("Final_Risk_Score", available_genes, "Match_ID")], 
                      by = "Match_ID")

cat("Final matching patients:", nrow(final_merged), "\n")



# 7. Automatic Column Name Detection
# Check for Age column
all_cols <- colnames(final_merged)
age_col <- all_cols[grep("age", all_cols, ignore.case = TRUE)][1]
stage_col <- all_cols[grep("stage", all_cols, ignore.case = TRUE)][1]

cat("Detected Age column:", age_col, "\n")
cat("Detected Stage column:", stage_col, "\n")

# 8. Rename columns to standard names for the model
# We create new columns so the original data remains safe
final_merged$Age_Standard <- as.numeric(final_merged[[age_col]])
final_merged$Stage_Standard <- as.factor(final_merged[[stage_col]])

# 9. Setup datadist again with new columns
dd <- datadist(final_merged)
options(datadist='dd')

# 10. Final Cox Model with Standardized Names
# We use Final_Risk_Score + Age + Stage
f <- cph(Surv(time, status) ~ Final_Risk_Score + Age_Standard + Stage_Standard, 
         data = final_merged, x=TRUE, y=TRUE, surv=TRUE)



# 1. Filter out non-positive survival time and non-numeric ages
# Removing any rows where time is 0 or negative
final_merged <- final_merged[final_merged$time > 0, ]

# 2. Ensure Age is numeric and Stage is a factor without missing levels
final_merged$Age_Standard <- as.numeric(final_merged[[age_col]])
final_merged$Stage_Standard <- as.factor(final_merged[[stage_col]])

# Remove any remaining NA rows in key columns
final_merged <- na.omit(final_merged[, c("time", "status", "Final_Risk_Score", "Age_Standard", "Stage_Standard")])

# 3. Setup datadist again
dd <- datadist(final_merged)
options(datadist='dd')

# 4. Final Cox Model
f <- cph(Surv(time, status) ~ Final_Risk_Score + Age_Standard + Stage_Standard, 
         data = final_merged, x=TRUE, y=TRUE, surv=TRUE)






# 1. Scaling Risk Score for better visualization (Optional but recommended)
# This will make the scale look like -18 to +6 instead of -18000 to 6000
final_merged$Final_Risk_Score_Scaled <- final_merged$Final_Risk_Score / 1000

# 2. Clean Stage Names (Shortening for space)
final_merged$Stage_Standard <- gsub("Stage ", "S-", final_merged$Stage_Standard)
final_merged$Stage_Standard <- as.factor(final_merged$Stage_Standard)

# 3. Update datadist with the scaled score
dd <- datadist(final_merged)
options(datadist='dd')

# 4. Re-run Cox Model with the Scaled Score
f <- cph(Surv(time, status) ~ Final_Risk_Score_Scaled + Age_Standard + Stage_Standard, 
         data = final_merged, x=TRUE, y=TRUE, surv=TRUE)

# 5. Save Professional Plots
pdf("Final_Professional_Validation.pdf", width = 11, height = 13)

# --- Nomogram Fix ---
surv_func <- Survival(f)
nom <- nomogram(f, fun=list(function(x) surv_func(365, x), 
                            function(x) surv_func(365*3, x),
                            function(x) surv_func(365*5, x)),
                funlabel=c("1-Year OS Prob", "3-Year OS Prob", "5-Year OS Prob"))

# Use plot() to control font sizes and spacing
plot(nom, 
     xfrac=0.3,           # Space for labels
     cex.axis=0.7,        # Font size for axis numbers
     cex.var=0.9,         # Font size for variable names (Age, Stage)
     lmgp=0.3,            # Margin between text and axis
     label.every=1)       # Show every label

# --- Calibration Plots ---
par(mfrow=c(3,1), mar=c(5,5,3,2))

# Function for clean calibration plots
plot_cal_final <- function(cal, title) {
  plot(cal, lwd=2, lty=1, conf.int=TRUE,
       errbar.col="blue", col="red",
       xlab="Nomogram Predicted Survival", ylab="Actual Observed Survival",
       main=title, xlim=c(0,1), ylim=c(0,1))
  abline(0, 1, lty=3, lwd=2, col="grey")
}

# 1-Year Calibration
cal1 <- calibrate(f, cmethod='KM', method="boot", u=365, m=100, B=1000)
plot_cal_final(cal1, "1-Year OS Calibration Curve")

# 3-Year Calibration
cal3 <- calibrate(f, cmethod='KM', method="boot", u=1095, m=100, B=1000)
plot_cal_final(cal3, "3-Year OS Calibration Curve")

# 5-Year Calibration
cal5 <- calibrate(f, cmethod='KM', method="boot", u=1825, m=100, B=1000)
plot_cal_final(cal5, "5-Year OS Calibration Curve")

dev.off()

cat("Everything finished successfully! Check 'Final_Professional_Validation.pdf'.\n")
