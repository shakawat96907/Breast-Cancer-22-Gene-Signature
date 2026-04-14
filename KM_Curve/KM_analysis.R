# library load
library(survival)
library(survminer)

# file load
clinical <- read.csv("TCGA_BRCA_Clinical_Final_1094.csv")
expression <- read.csv("TCGA_BRCA_Full_Expression_1094.csv")

# convert dot sign to hypen sign
colnames(expression) <- gsub("\\.", "-", colnames(expression))

# 22 signature genes with their coefficient
sig_data <- data.frame(
  Gene = c("AB306139", "ADAMTS15", "AFF3", "CCL8", "CHAD", "CLIC6", "CXCL13", "CYP4F22", 
           "ELOVL2", "FCER1A", "GRPR", "GSTM5", "MARCO", "MIR635", "NAT1", "PTGER3", 
           "SCARA5", "SLC7A2", "TCRVB", "TFF1", "TK1", "TP63"),
  Coeff = c(-0.059702430, -0.033946095, -0.008751713, 0.019218459, -0.006918279, 
            -0.025429251, -0.039243328, -0.003191574, -0.030036518, -0.008773601, 
            -0.013292914, -0.004513125, 0.057587821, -0.003089808, -0.020223151, 
            -0.010548777, -0.029573798, -0.044844060, -0.017301796, -0.015482446, 
            0.069595101, -0.053156012)
)

# check which genes are not available
common_genes <- intersect(sig_data$Gene, expression$gene_name)
missing_genes <- setdiff(sig_data$Gene, expression$gene_name)

cat("availb=able genes:", length(common_genes), "\n")
if(length(missing_genes) > 0) {
  cat(" Unavailable genes:", paste(missing_genes, collapse=", "), "\n")
}

# create data with available genes
selected_exp <- expression[expression$gene_name %in% common_genes, ]
rownames(selected_exp) <- selected_exp$gene_name
exp_matrix <- t(as.matrix(selected_exp[,-1])) 

# filter properly with 22 genes 
sig_subset <- sig_data[sig_data$Gene %in% common_genes, ]
exp_matrix <- exp_matrix[, sig_subset$Gene]

# risk score calculation 
risk_scores <- exp_matrix %*% sig_subset$Coeff

# merge with clinical data
risk_df <- data.frame(bcr_patient_barcode = rownames(risk_scores), Risk_Score = as.numeric(risk_scores))
final_merged <- merge(risk_df, clinical, by = "bcr_patient_barcode")
final_merged$Group <- ifelse(final_merged$Risk_Score > median(final_merged$Risk_Score), "High-Risk", "Low-Risk")

# ৯. KM Plot create 
fit <- survfit(Surv(time, status) ~ Group, data = final_merged)
ggsurvplot(fit, data = final_merged, pval = TRUE, risk.table = TRUE,
           title = "External Validation (TCGA-BRCA)",
           palette = c("#E41A1C", "#377EB8"), legend.labs = c("High-Risk", "Low-Risk"))
