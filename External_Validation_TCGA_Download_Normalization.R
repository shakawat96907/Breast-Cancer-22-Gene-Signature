if (!require("tidyverse")) install.packages("tidyverse")
if (!require("data.table")) install.packages("data.table")

library(tidyverse)
library(data.table)


# select folder directory
parent_folder <- "D:/Research Project/Breast Cancer (Done)/final/breast/GDCdata/TCGA-BRCA/Transcriptome_Profiling"

# check folder
dir.exists(parent_folder)
list.files(parent_folder)

# inside correct folder
parent_folder <- file.path(parent_folder, "Gene_Expression_Quantification")


# file detect
all_files <- list.files(parent_folder, recursive = TRUE, full.names = TRUE)

length(all_files)
head(all_files, 10)

# just .tsv filter
all_files <- all_files[grepl("\\.tsv$", all_files, ignore.case = TRUE)]

length(all_files)
head(all_files)


length(all_files)
head(all_files)


# FIXED: Optimized TCGA-BRCA Merging Script (Handles Duplicates)

library(data.table)
library(dplyr)

# 1. Path settings
parent_folder <- "D:/Research Project/Breast Cancer (Done)/final/breast/GDCdata/TCGA-BRCA/Transcriptome_Profiling/Gene_Expression_Quantification"
output_path <- "D:/Research Project/Breast Cancer (Done)/final/breast/TCGA_BRCA_Merged_Counts_Final.csv"

# 2. Detect Files
all_files <- list.files(path = parent_folder, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)
cat("Files detected:", length(all_files), "\n")

# 3. Reading and Cleaning with Duplicate Removal
datalist <- list()

cat("Reading and Cleaning files...\n")
for (i in 1:length(all_files)) {
  sample_id <- basename(dirname(all_files[i]))
  
  # Read and clean
  df <- fread(all_files[i], skip = 1, select = c("gene_name", "unstranded"))
  df <- df[-(1:5), ]
  
  # CRITICAL: Remove duplicate gene names within the same file before merging
  df <- df[!duplicated(gene_name), ]
  
  setnames(df, "unstranded", sample_id)
  datalist[[i]] <- df
  
  if (i %% 100 == 0) cat("Processed", i, "files...\n")
}

# 4. Merging with allow.cartesian=TRUE
cat("Merging into a single matrix. This might take a moment...\n")

# Using Reduce with a custom merge function that allows cartesian joins
final_matrix <- Reduce(function(x, y) {
  merge(x, y, by = "gene_name", all = TRUE, allow.cartesian = TRUE)
}, datalist)

# 5. Export Final Matrix
cat("Exporting to CSV...\n")
fwrite(final_matrix, output_path)

# 6. Success Check
if(file.exists(output_path)) {
  cat("\nSuccess! File saved at:", output_path, "\n")
} else {
  cat("\nError: File not found. Check permissions.\n")
}



# check the row and column 
dim(final_matrix)




# clinical data 
if (!require("TCGAbiolinks")) BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

# Direct TCGA-BRCA clinical dataset download
query_clin <- GDCquery(project = "TCGA-BRCA", 
                       data.category = "Clinical", 
                       data.type = "Clinical Supplement", 
                       data.format = "BCR Biotab")
GDCdownload(query_clin)
clinical_list <- GDCprepare(query_clin)



# column arrangement
clinical_full_data <- clinical_list$clinical_patient_brca
write.csv(clinical_full_data, "TCGA_BRCA_Full_Clinical_Data.csv", row.names = FALSE)



# more specific dataset for TCGA
clinical_subset <- clinical_list$clinical_patient_brca %>%
  select(
    bcr_patient_barcode, 
    gender, 
    race,
    vital_status, 
    death_days_to,           
    last_contact_days_to,    
    age_at_diagnosis,
    ajcc_pathologic_tumor_stage,
    er_status_by_ihc,        # Estrogen Receptor status
    pr_status_by_ihc,        # Progesterone Receptor status
    her2_status_by_ihc,      # HER2 status
    histological_type
  )

# Save
write.csv(clinical_subset, "TCGA_BRCA_Clinical_Processed.csv", row.names = FALSE)
head(clinical_subset)




meta <- getResults(query_meta)
map_data <- data.frame(file_id = meta$file_id, barcode = meta$cases)

current_names <- colnames(final_matrix)
new_names <- map_data$barcode[match(current_names, map_data$file_id)]

new_names[is.na(new_names)] <- current_names[is.na(new_names)]

colnames(final_matrix) <- new_names

head(colnames(final_matrix), 10)


# save
write.csv(final_matrix, "TCGA_BRCA_Final_Matrix_Mapped.csv", row.names = FALSE)



# final checking for both dataset
# 1. Trim Expression Barcodes to 12 characters (Match Clinical Format)
colnames(final_matrix) <- substr(colnames(final_matrix), 1, 12)

# 2. Re-check the matching
common_samples <- intersect(colnames(final_matrix), clinical_clean$bcr_patient_barcode)

cat("\n--- Updated Matching Statistics ---\n")
cat("Number of patients found in both datasets:", length(common_samples), "\n")

# 3. See the first few matched barcodes
print(head(common_samples, 10))




# 1. Identify common patient IDs
common_ids <- intersect(colnames(final_matrix), clinical_clean$bcr_patient_barcode)

# 2. Filter Expression Data for these 1,094 patients
# Keeping 'gene_name' as the first column
final_expression_export <- final_matrix[, c("gene_name", common_ids)]

# 3. Filter Clinical Data for these 1,094 patients
final_clinical_export <- clinical_clean[clinical_clean$bcr_patient_barcode %in% common_ids, ]

# 4. Save to CSV files for future research
write.csv(final_expression_export, "TCGA_BRCA_Expression_Final_1094.csv", row.names = FALSE)
write.csv(final_clinical_export, "TCGA_BRCA_Clinical_Final_1094.csv", row.names = FALSE)

# 5. Confirmation message
print("Success! Both files have been saved as 'TCGA_BRCA_..._1094.csv'")



# 1. Check if final_matrix exists in your environment


# 1. Ensure common_ids are set
common_ids <- intersect(colnames(final_matrix), clinical_clean$bcr_patient_barcode)

# 2. Force final_matrix to be a data frame before subsetting
# This step fixes the NULL dimension error
final_matrix_df <- as.data.frame(final_matrix)

# 3. Now try to extract the columns
final_expression_final <- final_matrix_df[, c("gene_name", common_ids)]

# 4. Verification - This MUST show two numbers now
print("--- Check Dimension Now ---")
print(dim(final_expression_final))

# 5. Save if it's finally a table
if(!is.null(dim(final_expression_final))) {
  write.csv(final_expression_final, "TCGA_BRCA_Full_Expression_Final.csv", row.names = FALSE)
  print("SUCCESS: Full table saved to your PC!")
} else {
  print("STILL NULL: Please run 'class(final_matrix)' and tell me what it says.")
}








# 1. Clinical Data Check
clinical_final <- read.csv("TCGA_BRCA_Clinical_Final_1094.csv")
cat("--- Clinical Data (Survival Info) ---\n")
cat("Dimensions (Rows x Columns):", dim(clinical_final), "\n")
cat("Column Names:\n")
print(colnames(clinical_final))
cat("First 5 Patient IDs:\n")
print(head(clinical_final$bcr_patient_barcode, 5))




# 2. Expression Data Check
# Note: This file is large, so we only check the structure
expression_final <- read.csv("TCGA_BRCA_Full_Expression_1094.csv")
cat("\n--- Expression Data (Gene Values) ---\n")
cat("Dimensions (Genes x Samples):", dim(expression_final), "\n")
cat("First 5 Column Names (Should be 'gene_name' + Barcodes):\n")
print(colnames(expression_final)[1:5])
cat("First 5 Gene Names:\n")
print(head(expression_final$gene_name, 5))




signature_final <- read.csv("LASSO_signature.csv")
print(signature_final)



