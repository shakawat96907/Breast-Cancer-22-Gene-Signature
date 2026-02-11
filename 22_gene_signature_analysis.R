# Install required Packages
install.packages(c("dplyr","tidyverse","GEOquery","data.table","ff"))

#Load required packages
library(dplyr)
library(tidyverse)
library(GEOquery)
library(data.table)


# Load clinical data
metadata.modified <- read.csv ("phenodata_clean.csv")
head(metadata.modified)


# Load expression data
library(ff)
expr_data_ff <- read.csv("GSE96058.csv", header = TRUE, stringsAsFactors = FALSE)

# Preview first few rows and dimensions
head(expr_data_ff)
dim(expr_data_ff)
)


# Rename first column to "Gene"
colnames(expr_data_ff)[1] <- "Gene"

# Function to clean sample/column names
clean_names <- function(x) {
  x <- trimws(x)  # remove leading/trailing spaces
  x <- gsub("repl|replicate|rep|_rep|_repl|_R|R$|v2|dup|_old|\\.1|\\.2", "", x, ignore.case = TRUE)  # remove replicate/duplicate indicators
  x <- toupper(x)  # convert to uppercase
  return(x)
}


# Clean sample names in expression data (exclude Gene column)
colnames(expr_data_ff)[-1] <- clean_names(colnames(expr_data_ff)[-1])

# Clean Sample_IDs in clinical data
metadata.modified$Sample_ID <- clean_names(metadata.modified$Sample_ID)


# Keep only common samples in metadata
expr_samples <- colnames(expr_data_ff)[-1]
meta_samples <- metadata.modified$Sample_ID

common_samples <- intersect(expr_samples, meta_samples)
cat("Common samples before NA filtering:", length(common_samples), "\n")

# Filter metadata to keep only samples present in expression data
metadata.modified <- metadata.modified %>%
  filter(Sample_ID %in% common_samples)

# Drop rows with any NA (complete clinical cases only)
metadata.modified <- metadata.modified %>%
  distinct(Sample_ID, .keep_all = TRUE) %>%
  drop_na()  # remove rows with any missing value

cat("Metadata after removing NAs:", nrow(metadata.modified), "rows remain\n")

# Recalculate common samples (after NA filtering)
common_samples <- intersect(colnames(expr_data_ff)[-1], metadata.modified$Sample_ID)
cat("Common samples after NA filtering:", length(common_samples), "\n")

# Numeric order sorting
extract_first_number <- function(x){
  m <- regmatches(x, regexpr("\\d+", x))
  as.integer(ifelse(length(m) == 0 | m == "", NA, m))
}

nums <- extract_first_number(common_samples)
ord <- order(nums, common_samples, na.last = TRUE)
common_samples <- unique(common_samples[ord])

cat("Common samples sorted (head/tail):\n")
print(head(common_samples, 10))
print(tail(common_samples, 10))

# Chunk-wise merge
chunk_size <- 100
chunks <- split(common_samples, ceiling(seq_along(common_samples) / chunk_size))

out_file <- "expr_long_meta_full.csv"
if(file.exists(out_file)) file.remove(out_file)

first_write <- TRUE
for(i in seq_along(chunks)){
  cols_ok <- chunks[[i]]
  
  cat(sprintf("Processing chunk %d/%d : %s -> %s\n", 
              i, length(chunks), cols_ok[1], cols_ok[length(cols_ok)]))
  
  expr_chunk_ff <- expr_data_ff[, c("Gene", cols_ok), drop = FALSE]
  
  expr_long <- pivot_longer(expr_chunk_ff,
                            cols = -Gene,
                            names_to = "Sample_ID",
                            values_to = "Expression")
  
  expr_long_meta <- left_join(expr_long, metadata.modified, by = "Sample_ID") %>%
    filter(!is.na(Sample_ID)) %>%
    distinct()
  
  fwrite(expr_long_meta, out_file, append = !first_write)
  first_write <- FALSE
  
  rm(expr_chunk_ff, expr_long, expr_long_meta)
  gc()
}

cat("\n All chunks processed & saved to", out_file, "\n")

# Verify final result
expr_final <- fread(out_file)
final_unique_samples <- unique(expr_final$Sample_ID)
cat("Final unique samples in merged data:", length(final_unique_samples), "\n")


# Load required libraries
library(R.utils)       # for gzip compression
library(data.table)    # for efficient data manipulation

# Compress merged CSV file
gzip("expr_long_meta_full.csv", 
     destname = "expr_long_meta_full.csv.gz", 
     overwrite = TRUE)
cat("File compressed & saved as expr_long_meta_full.csv.gz\n")

# Load compressed file efficiently
expr_data <- fread("expr_long_meta_full.csv.gz")

# Basic checks
cat("Total rows in merged data:", nrow(expr_data), "\n")
cat("Final unique samples:", length(unique(expr_data$Sample_ID)), "\n")





#  Expression Data Processing

# Load required libraries
library(data.table)   # fast data manipulation
library(dplyr)        # data wrangling
library(ggplot2)      # plotting
library(pheatmap)     # heatmaps
library(reshape2)     # dcast for wide-format conversion

# Load compressed expression CSV
expr_data <- fread("expr_long_meta_full.csv.gz")

# Verify unique samples
cat("Unique samples loaded:", length(unique(expr_data$Sample_ID)), "\n")

# Quick preview
head(expr_data)

# Remove duplicates and NAs
expr_data <- unique(expr_data, by = c("Gene", "Sample_ID"))
expr_data <- expr_data[!is.na(Expression)]

cat("Unique samples after removing NAs/duplicates:", length(unique(expr_data$Sample_ID)), "\n")

# Log2 transformation if needed
q <- quantile(expr_data$Expression, probs = c(0.5, 0.9, 0.99), na.rm = TRUE)
if (q["99%"] > 50) {
  expr_data[, Expression_log2 := log2(Expression + 1)]
  cat("ℹ️  Log2 transformation applied\n")
} else {
  expr_data[, Expression_log2 := Expression]
  cat("ℹ️  Log2 transformation not needed\n")
}

# Filter low-expressed genes
prop_dt <- expr_data[, .(
  prop = mean(Expression_log2 > 1, na.rm = TRUE),
  meanExpr = mean(Expression_log2, na.rm = TRUE)
), by = Gene]

# Keep genes expressed in at least 20% of samples
keep_genes <- prop_dt[prop >= 0.20]$Gene
expr_filt <- expr_data[Gene %in% keep_genes]
cat("Genes after low-expression filter:", length(unique(expr_filt$Gene)), "\n")

# Select top 500 most variable genes
var_dt <- expr_filt[, .(var = var(Expression_log2, na.rm = TRUE)), by = Gene]
top_genes <- var_dt[order(-var)][1:500]$Gene
cat("Top 500 variable genes selected\n")

# Convert to matrix (genes x samples)
tmp <- expr_filt[Gene %in% top_genes, .(Gene, Sample_ID, Expression_log2)]
mat_dt <- dcast(tmp, Gene ~ Sample_ID, value.var = "Expression_log2")
rownames(mat_dt) <- mat_dt$Gene
mat <- as.matrix(mat_dt[, -1, with = FALSE])

cat("Expression matrix created with dimensions:", dim(mat), "\n")

# Save filtered long-format data as compressed CSV
fwrite(expr_filt, "expr_filtered_long.csv.gz", compress = "gzip")
cat("Filtered expression data saved as expr_filtered_long.csv.gz\n")




                                                                    
# Pre-processing & Quality Control

# Load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(reshape2)
library(factoextra)
library(scales)
library(officer)
library(flextable)
library(ggrepel)

# Load filtered dataset
expr_data <- fread("expr_filtered_long.csv.gz")
cat("Loaded filtered dataset with", nrow(expr_data), "rows\n")
str(expr_data)
head(expr_data)
cat("Unique samples loaded:", length(unique(expr_data$Sample_ID)), "\n")

# Prepare clinical table (Distribution Table)
clin_data <- unique(expr_data[, .(
  Sample_ID, age, tumor_size, lymph_node_status,
  ER_status, PR_status, HER2_status, Ki67_status,
  tumor_grade, PAM50_subtype
)])

# Convert binary variables to human-readable
clin_data <- clin_data %>%
  mutate(
    ER_status = ifelse(ER_status == 1, "Positive", "Negative"),
    PR_status = ifelse(PR_status == 1, "Positive", "Negative"),
    HER2_status = ifelse(HER2_status == 1, "Positive", "Negative"),
    Ki67_status = ifelse(Ki67_status == 1, "High", "Low")
  )

# Numeric summary function
num_summary_df <- function(var_name, var) {
  mean_val <- round(mean(var, na.rm = TRUE), 2)
  sd_val   <- round(sd(var, na.rm = TRUE), 2)
  data.frame(
    Variable = var_name,
    Level = "-",
    Count = "-",
    Percent = paste0(mean_val, " ± ", sd_val),
    stringsAsFactors = FALSE
  )
}

# Categorical summary function
cat_summary_df <- function(var_name, var) {
  df <- clin_data %>% count({{var}}) %>% 
    mutate(percent = round(n / sum(n) * 100, 1)) %>%
    rename(Level = {{var}}, Count = n, Percent = percent)
  df$Percent <- paste0(df$Percent, "%")
  df$Variable <- var_name
  df <- df[, c("Variable", "Level", "Count", "Percent")]
  return(df)
}

# Numeric summaries
numeric_tables <- rbind(
  num_summary_df("Age (years)", clin_data$age),
  num_summary_df("Tumor Size (cm)", clin_data$tumor_size)
)

# Categorical summaries
categorical_tables <- rbind(
  cat_summary_df("Tumor Grade", tumor_grade),
  cat_summary_df("Lymph Node Status", lymph_node_status),
  cat_summary_df("ER Status", ER_status),
  cat_summary_df("PR Status", PR_status),
  cat_summary_df("HER2 Status", HER2_status),
  cat_summary_df("Ki67 Status", Ki67_status),
  cat_summary_df("PAM50 Subtype", PAM50_subtype)
)

# Merge all
Table1_final <- rbind(numeric_tables, categorical_tables)
Table1_final$Count <- gsub(",", "", Table1_final$Count)

# Create polished Word table
doc <- read_docx() %>%
  body_add_par("Table 1: Clinical Characteristics", style = "heading 1") %>%
  body_add_flextable(
    flextable(Table1_final) %>%
      autofit() %>%
      bold(i = ~Level == "-", bold = TRUE) %>%              # bold numeric rows
      bg(i = seq(1, nrow(Table1_final), 2), bg = "#F7F7F7") %>% # alternate shading
      border_outer(part = "all", border = fp_border(color = "black", width = 1)) %>%
      border_inner_h(part = "body", border = fp_border(color = "gray", width = 0.5)) %>%
      set_header_labels(
        Variable = "Variable",
        Level = "Level",
        Count = "Count",
        Percent = "Percent"
      )
  )
print(doc, target = "Table1_clinical_final_polished.docx")
cat("Clinical Table saved as Table1_clinical_final_polished.docx\n")

# Preprocess expression data
expr_data <- unique(expr_data, by = c("Gene", "Sample_ID"))
expr_data <- expr_data[!is.na(Expression)]
cat("Unique samples after removing duplicates/NAs:", length(unique(expr_data$Sample_ID)), "\n")

# Log2 transform if needed
q <- quantile(expr_data$Expression, probs = c(0.5, 0.9, 0.99), na.rm = TRUE)
if (q["99%"] > 50) {
  expr_data[, Expression_log2 := log2(Expression + 1)]
  cat("ℹ️ Log2 transformation applied\n")
} else {
  expr_data[, Expression_log2 := Expression]
  cat("ℹ️ Log2 transformation not needed\n")
}

# Filter low-expressed genes
prop_dt <- expr_data[, .(
  prop = mean(Expression_log2 > 1, na.rm = TRUE),
  meanExpr = mean(Expression_log2, na.rm = TRUE)
), by = Gene]
keep_genes <- prop_dt[prop >= 0.20]$Gene
expr_filt <- expr_data[Gene %in% keep_genes]
cat("✅ Genes after low-expression filter:", length(unique(expr_filt$Gene)), "\n")

# Top variable genes
var_dt <- expr_filt[, .(var = var(Expression_log2, na.rm = TRUE)), by = Gene]
top_genes <- var_dt[order(-var)][1:5000]$Gene

# Make expression matrix
tmp <- expr_filt[Gene %in% top_genes, .(Gene, Sample_ID, Expression_log2)]
mat_dt <- dcast(tmp, Gene ~ Sample_ID, value.var = "Expression_log2")
rownames(mat_dt) <- mat_dt$Gene
mat <- as.matrix(mat_dt[, -1, with = FALSE])

# Clean matrix
mat_clean <- mat
mat_clean[!is.finite(mat_clean)] <- NA
mat_clean <- mat_clean[complete.cases(mat_clean), ]



# PCA
pca_res <- prcomp(t(mat_clean), scale. = TRUE, center = TRUE)
var_explained <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)
xlab_pc1 <- paste0("PC1 (", var_explained[1], "%)")
ylab_pc2 <- paste0("PC2 (", var_explained[2], "%)")

pc_df <- data.frame(
  Sample_ID = colnames(mat_clean),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Subtype = expr_data$PAM50_subtype[match(colnames(mat_clean), expr_data$Sample_ID)]
)

# Subtype percentage for legend
subtype_tab <- pc_df %>%
  count(Subtype) %>%
  mutate(
    Percent = round(100 * n / sum(n), 1),
    label = paste0(Subtype, " (", Percent, "%)")
  )
subtype_labels <- subtype_tab$label
names(subtype_labels) <- subtype_tab$Subtype

colors <- brewer.pal(n = length(unique(pc_df$Subtype)), name = "Set2")

# PCA plot
p <- ggplot(pc_df, aes(x = PC1, y = PC2, color = Subtype, fill = Subtype)) +
  geom_point(size = 3, shape = 21, stroke = 0.5) +
  stat_ellipse(aes(group = Subtype), type = "norm", linetype = 2, color = "black") +
  scale_color_manual(values = colors, labels = subtype_labels) +
  scale_fill_manual(values = colors, labels = subtype_labels) +
  labs(x = xlab_pc1, y = ylab_pc2) +
  theme_minimal(base_size = 8) +
  theme(
    plot.title = element_blank(),
    plot.caption = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

# Save PCA EPS
setEPS()
postscript("PCA_plot_advanced.eps", width = 6, height = 6, family = "Helvetica", horizontal = FALSE, onefile = FALSE, paper = "special")
print(p)
dev.off()
cat("✅ PCA plot saved as PCA_plot_advanced.eps\n")



# QC Boxplot by Subtype
expr_plot <- expr_filt %>%
  select(Sample_ID, Gene, Expression_log2) %>%
  mutate(Group = expr_data$PAM50_subtype[match(Sample_ID, expr_data$Sample_ID)]) %>%
  filter(!is.na(Group))

subtypes <- unique(expr_plot$Group)
colors <- setNames(brewer.pal(length(subtypes), "Set2"), subtypes)

# Order groups by median expression
group_order <- expr_plot %>%
  group_by(Group) %>%
  summarize(median_expr = median(Expression_log2, na.rm = TRUE), .groups = "drop") %>%
  arrange(median_expr) %>%
  pull(Group)
expr_plot$Group <- factor(expr_plot$Group, levels = group_order)

summary_stats <- expr_plot %>%
  group_by(Group) %>%
  summarize(median_expr = median(Expression_log2, na.rm = TRUE), n = n(), .groups = "drop")

y_max <- max(expr_plot$Expression_log2, na.rm = TRUE)

# Save Boxplot TIFF
tiff("QC_boxplot_violin_median_sample.tiff", width = 10, height = 5, units = "in", res = 1000, compression = "lzw")
ggplot(expr_plot, aes(x = Group, y = Expression_log2, fill = Group)) +
  geom_violin(alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, color = "black", fatten = 2) +
  geom_jitter(aes(color = Group), width = 0.15, size = 0.5, alpha = 0.6) +
  geom_point(data = summary_stats, aes(x = Group, y = median_expr), color = "lavender", size = 2) +
  geom_text(data = summary_stats, aes(x = Group, y = median_expr, label = round(median_expr, 2)), vjust = -0.7, color = "black", size = 2.5) +
  geom_text(data = summary_stats, aes(x = Group, y = y_max + 0.5, label = paste0("n=", n)), vjust = 0, color = "black", size = 2.5) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_minimal(base_size = 8) +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.title = element_blank(), plot.caption = element_blank(), legend.title = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank()) +
  labs(y = "log2(Expression + 1)")
dev.off()
cat("QC Boxplot saved as QC_boxplot_violin_median_sample.tiff\n")







# QC Boxplot / Violin Plot by PAM50 Subtype

library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Safety check
stopifnot(exists("expr_data"))

# PAM50 Subtype cleanup & ordering
expr_data$PAM50_subtype <- trimws(expr_data$PAM50_subtype)
expr_data$PAM50_subtype <- gsub("_", "-", expr_data$PAM50_subtype)
expr_data$PAM50_subtype <- dplyr::recode(
  expr_data$PAM50_subtype,
  "Basal"  = "Basal-like",
  "Her2"   = "HER2-enriched",
  "LumA"   = "Luminal A",
  "LumB"   = "Luminal B",
  "Normal" = "Normal-like"
)

final_order <- c("Basal-like", "HER2-enriched", "Luminal A", "Luminal B", "Normal-like")
expr_data$PAM50_subtype <- factor(expr_data$PAM50_subtype, levels = final_order)

# Prepare plotting dataframe
expr_filt <- expr_data %>%
  select(Sample_ID, Gene, Expression_log2, PAM50_subtype) %>%
  rename(Group = PAM50_subtype) %>%
  filter(!is.na(Group))

# True sample size per subtype
true_n <- expr_filt %>%
  distinct(Sample_ID, Group) %>%
  count(Group)

# Median expression per group
summary_stats <- expr_filt %>%
  group_by(Group) %>%
  summarize(median_expr = median(Expression_log2, na.rm = TRUE), .groups = "drop")

y_max <- max(expr_filt$Expression_log2, na.rm = TRUE)

# Colors
colors <- setNames(brewer.pal(length(final_order), "Set2"), final_order)

# Save HIGH-RES EPS (optional)
setEPS()
postscript("QC_boxplot_violin_median_sample.eps",
           width = 10, height = 5, horizontal = FALSE,
           onefile = FALSE, paper = "special")

ggplot(expr_filt, aes(x = Group, y = Expression_log2, fill = Group)) +
  geom_violin(alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.1, color = "black", outlier.size = 0.7) +
  geom_jitter(aes(color = Group), width = 0.12, size = 0.5, alpha = 0.5) +
  geom_point(data = summary_stats, aes(x = Group, y = median_expr), size = 3, color = "pink") +
  geom_text(data = summary_stats, aes(x = Group, y = median_expr, label = round(median_expr, 2)),
            vjust = -0.7, size = 3.2, fontface = "bold") +
  geom_text(data = true_n, aes(x = Group, y = y_max + 0.5, label = paste0("N=", n)),
            size = 3.2, fontface = "bold") +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(y = "log2(Expression + 1)") +
  theme_minimal(base_size = 11) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

dev.off()

# Save HIGH-RES TIFF (optional)
tiff("QC_boxplot_violin_median_sample.tiff",
     width = 10, height = 5, units = "in", res = 1000, compression = "lzw")

ggplot(expr_filt, aes(x = Group, y = Expression_log2, fill = Group)) +
  geom_violin(alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.1, color = "black", outlier.size = 0.7) +
  geom_jitter(aes(color = Group), width = 0.12, size = 0.5, alpha = 0.5) +
  geom_point(data = summary_stats, aes(x = Group, y = median_expr), size = 3, color = "pink") +
  geom_text(data = summary_stats, aes(x = Group, y = median_expr, label = round(median_expr, 2)),
            vjust = -0.7, size = 3.2, fontface = "bold") +
  geom_text(data = true_n, aes(x = Group, y = y_max + 0.5, label = paste0("N=", n)),
            size = 3.2, fontface = "bold") +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(y = "log2(Expression + 1)") +
  theme_minimal(base_size = 11) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

dev.off()

cat("QC Boxplot saved as EPS and TIFF\n")







# Clinical Characteristics Barplots

library(data.table)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(gridExtra)
library(grid)
library(scales)

# Load filtered dataset
expr_data <- fread("expr_filtered_long.csv.gz")
stopifnot(exists("expr_data"))

# Prepare clinical data
clin_data <- unique(expr_data[, .(
  Sample_ID,
  tumor_grade,
  lymph_node_status,
  ER_status,
  PAM50_subtype,
  OS_event,
  OS_days
)])

# Clean PAM50 subtypes
clin_data$PAM50_subtype <- trimws(clin_data$PAM50_subtype)
clin_data$PAM50_subtype <- gsub("_", "-", clin_data$PAM50_subtype)
clin_data$PAM50_subtype <- recode(clin_data$PAM50_subtype,
                                  "Basal" = "Basal-like",
                                  "Her2" = "HER2-enriched",
                                  "LumA" = "Luminal A",
                                  "LumB" = "Luminal B",
                                  "Normal" = "Normal-like")
# Factor levels
subtype_order <- c("Basal-like", "HER2-enriched", "Luminal A", "Luminal B", "Normal-like")
clin_data$PAM50_subtype <- factor(clin_data$PAM50_subtype, levels = subtype_order)

# Recode binary variables
clin_data$ER_status <- factor(clin_data$ER_status, levels = c(0,1), labels = c("ER Negative","ER Positive"))
clin_data$OS_event  <- factor(clin_data$OS_event, levels = c(0,1), labels = c("Alive","Death"))

# 3) Color palettes
col_grade  <- c("G1"="#4DBBD5", "G2"="#00A087", "G3"="#E64B35")
col_node   <- c("NodeNegative"="#3C5488", "NodePositive"="#F39B7F")
col_os     <- c("Alive"="#00A087", "Death"="#E64B35")
col_pam50  <- c("Basal-like"="#7E6148", "HER2-enriched"="#DC0000",
                "Luminal A"="#00A087", "Luminal B"="#3C5488",
                "Normal-like"="#B09C85")
col_er     <- c("ER Negative"="#E64B35", "ER Positive"="#4DBBD5")

# 4) Barplot function (% + N)
bar_fun <- function(var, xlab_txt, cols, rotate_x = FALSE) {
  df <- clin_data[, .N, by = var]
  df[, perc := round(N / sum(N) * 100, 1)]
  df[[var]] <- factor(df[[var]], levels = names(cols))
  
  p <- ggplot(df, aes_string(x = var, y = "N", fill = var)) +
    geom_bar(stat = "identity", width = 0.7, color = "black") +
    geom_text(aes(label = paste0(N, "\n(", perc, "%)")),
              position = position_stack(vjust = 0.5), size = 3.5, color = "white") +
    scale_fill_manual(values = cols) +
    theme_classic(base_size = 9) +
    labs(x = xlab_txt, y = "Count") +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10)
    )
  
  if (rotate_x) {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  return(p)
}

# 5) Individual panels
p1 <- bar_fun("tumor_grade", "Tumor grade", col_grade)
p2 <- bar_fun("lymph_node_status", "Lymph node status", col_node)
p3 <- bar_fun("OS_event", "Overall survival event", col_os)
p4 <- bar_fun("PAM50_subtype", "PAM50 subtype", col_pam50, rotate_x = TRUE)
p5 <- bar_fun("ER_status", "ER status", col_er)

p6 <- ggplot(clin_data, aes(x = OS_days)) +
  geom_histogram(bins = 40, fill = "#8491B4", color = "black") +
  theme_classic(base_size = 9) +
  labs(x = "Overall survival (days)", y = "Count")

# 6) Combine panels (grid)
final_plot <- grid.arrange(p1, p2, p3,
                           p4, p5, p6,
                           ncol = 3)

# 7) Export EPS + TIFF (1000 dpi)
# EPS
setEPS()
postscript("Figure1D_Clinical_Distribution_FINAL.eps",
           width = 8, height = 6,
           horizontal = FALSE, onefile = FALSE, paper = "special")
grid.draw(final_plot)
dev.off()

# TIFF
tiff("Figure1D_Clinical_Distribution_FINAL.tiff",
     width = 8, height = 6, units = "in",
     res = 1000, compression = "lzw")
grid.draw(final_plot)
dev.off()

cat("Clinical distribution plots saved as EPS and TIFF\n")













                                                        ########## step:3 : Genome-wide Survival Screening #########



library(data.table)   # Efficient data handling
library(survival)     # Cox PH model
library(dplyr)        # Data manipulation
library(ggplot2)      # Volcano plot
library(pheatmap)     # Heatmap
library(officer)      # Word export
library(flextable)    # Table formatting





# dt is survival ready dataset

dt <- expr_data[, .(
  Gene,
  Sample_ID,
  Expression_log2,
  OS_days,
  OS_event
)]

View(dt)



#Genome-wide Univariate Cox loop

cox_results <- dt[, {
  if(sd(Expression_log2, na.rm = TRUE) == 0) return(NULL)
  fit <- coxph(Surv(OS_days, OS_event) ~ Expression_log2)
  s <- summary(fit)
  data.table(
    HR        = s$coefficients[,"exp(coef)"],
    CI_lower  = s$conf.int[,"lower .95"],
    CI_upper  = s$conf.int[,"upper .95"],
    p_value   = s$coefficients[,"Pr(>|z|)"]
  )
}, by = Gene]




#Risk / Protective classification + ranking


cox_results <- cox_results %>%
  mutate(
    log2HR = log2(HR),
    Risk_type = case_when(
      p_value < 0.05 & HR > 1 ~ "Risk",
      p_value < 0.05 & HR < 1 ~ "Protective",
      TRUE ~ "Non-significant"
    )
  ) %>%
  arrange(p_value)





#Significant genes only (Manuscript)
#Reviewer-safe main table

table2_main <- cox_results %>%
  filter(p_value < 0.05) %>%
  select(Gene, HR, CI_lower, CI_upper, p_value, Risk_type) %>%
  arrange(p_value)

# Save CSV
fwrite(table2_main, "Table2_UnivariateCox_SignificantGenes.csv")

# Save Word
ft <- flextable(table2_main) %>% autofit() %>% theme_box()
doc <- read_docx() %>% body_add_flextable(ft)
print(doc, target = "Table2_UnivariateCox_SignificantGenes.docx")



fwrite(cox_results, "TableS1_UnivariateCox_AllGenes.csv")




                         



                                                    ######## Figure 2A: Cox-Volcano Plot ##########

volcano <- ggplot(cox_results, aes(x = log2HR, y = -log10(p_value))) +
  geom_point(aes(color = Risk_type), size = 1.5) +
  scale_color_manual(values = c(
    "Risk" = "red",
    "Protective" = "blue",
    "Non-significant" = "grey50"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(x = "log2(Hazard Ratio)", y = "-log10(p-value)", color = "Gene type") +
  theme_classic(base_size = 14)

ggsave(
  filename = "Figure2A_CoxVolcano.eps",
  plot = volcano,
  device = cairo_ps,
  width = 6,
  height = 5,
  dpi = 1000
)






#For labeling top 30 genes name in plot


library(ggrepel)  # For labels

# Top 30 genes to label by p-value
top_genes <- cox_results %>%
  filter(p_value < 0.05) %>%
  arrange(p_value) %>%
  slice(1:30) %>%
  pull(Gene)

# Volcano plot with top 30 labels
volcano <- ggplot(cox_results, aes(x = log2HR, y = -log10(p_value))) +
  geom_point(aes(color = Risk_type), size = 1.5) +
  scale_color_manual(values = c(
    "Risk" = "#D55E00",
    "Protective" = "#0072B2",
    "Non-significant" = "#C7C"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = cox_results[Gene %in% top_genes, ],
    aes(label = Gene),
    size = 3,
    box.padding = 0.4,
    max.overlaps = 50
  ) +
  labs(x = "log2(Hazard Ratio)", y = "-log10(p-value)", color = "Gene type") +
  theme_classic(base_size = 14)

# Save EPS
ggsave(
  filename = "Figure2A_CoxVolcano_top30.eps",
  plot = volcano,
  device = cairo_ps,
  width = 6,
  height = 5,
  dpi = 1000
)









 
                                                ######## Figure 2B: Heatmap (Top 20 genes) #######

apply(mat, 1, sd)
sum(is.na(mat))


library(data.table)
library(pheatmap)
library(dplyr)

# 1️⃣ Top 20 significant genes
top20_genes <- cox_results %>%
  filter(p_value < 0.05) %>%
  arrange(p_value) %>%
  slice(1:20) %>%
  pull(Gene)

# 2️⃣ Subset expression and remove genes with zero variance
heat_dt <- expr_data[Gene %in% top20_genes, .(Sample_ID, Gene, Expression_log2, OS_event)]
mat <- dcast(heat_dt, Gene ~ Sample_ID, value.var = "Expression_log2")
mat <- as.data.table(mat)

# Remove zero variance genes
var_check <- apply(mat[,-1, with=FALSE], 1, sd, na.rm=TRUE)
mat <- mat[var_check != 0, ]

# Convert to matrix and set rownames
mat_z <- as.matrix(mat[,-1, with=FALSE])
rownames(mat_z) <- mat$Gene

# 3️⃣ Z-score by gene
mat_z <- t(scale(t(mat_z)))

# 4️⃣ Column annotation
anno_col <- heat_dt[, .(Sample_ID, OS_event)] %>% unique() %>% as.data.frame()
rownames(anno_col) <- anno_col$Sample_ID
anno_col$Sample_ID <- NULL

# Ensure colnames match annotation
mat_z <- mat_z[, rownames(anno_col), drop = FALSE]

# 5️⃣ EPS Heatmap
setEPS()
postscript("Figure2B_Heatmap_Top20SurvivalGenes_fixed.eps",
           width = 8, height = 6, horizontal = FALSE, paper = "special")

pheatmap(
  mat_z,
  annotation_col = anno_col,
  show_colnames = FALSE,
  fontsize_row = 9,
  color = colorRampPalette(c("navy","white","firebrick3"))(100)
)

dev.off()
message("Figure2B Heatmap EPS generated successfully (fixed)!")











                                                     ############### step:4: DEG analysis ################



library(data.table)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

#  Load filtered dataset 
expr_data <- fread("expr_filtered_long.csv.gz")




                                                            ######### substep:4.1 ##########


library(data.table)
setDT(expr_data)

# extract unique clinical info per sample
clin_data <- unique(
  expr_data[, .(
    Sample_ID,
    tumor_grade
  )]
)

# check distribution
table(clin_data$tumor_grade)



#tumor geade recoding
clin_data[, tumor_grade2 :=
            fifelse(
              tumor_grade %in% c("G1", "G2"), "Low",
              fifelse(tumor_grade == "G3", "High", NA_character_)
            )
]

# keep only valid samples
clin_data <- clin_data[!is.na(tumor_grade2)]

# check balance
table(clin_data$tumor_grade2)


#mandatory check
length(unique(clin_data$Sample_ID))

# any missing values?
any(is.na(clin_data$tumor_grade2))

# preview
head(clin_data)





                                                        ######## substep 4.2 ############

cox_results <- read.csv("Significant.csv")


# survival-associated genes
survival_genes <- cox_results$Gene[cox_results$p_value < 0.05]

length(survival_genes)


#Filter expression data
expr_surv <- expr_data[Gene %in% survival_genes]

# keep only Low / High samples
expr_surv <- expr_surv[Sample_ID %in% clin_data$Sample_ID]

# check size
dim(expr_surv)




# cast to Gene x Sample matrix
expr_mat <- dcast(
  expr_surv,
  Gene ~ Sample_ID,
  value.var = "Expression_log2"
)

# convert to matrix
gene_ids <- expr_mat$Gene
expr_mat <- as.matrix(expr_mat[, -1])
rownames(expr_mat) <- gene_ids

# check
dim(expr_mat)[1:2]



#Match group labels to matrix columns
group <- clin_data[
  match(colnames(expr_mat), Sample_ID),
  tumor_grade2
]

group <- factor(group, levels = c("Low", "High"))

# check
table(group)







                                                      ############# substep: 4.3 #################



#packages
library(data.table)
library(limma)



#Survival genes + valid samples filter
deg_dt <- expr_data[
  Gene %in% survival_genes &
    Sample_ID %in% names(group)
]



#Expression matrix (Gene × Sample)
expr_mat <- dcast(
  deg_dt,
  Gene ~ Sample_ID,
  value.var = "Expression_log2"
)

expr_mat <- as.data.frame(expr_mat)
rownames(expr_mat) <- expr_mat$Gene
expr_mat$Gene <- NULL
expr_mat <- as.matrix(expr_mat)




#Group factor (reference = Low grade)
group <- factor(group, levels = c("Low_grade", "High_grade"))



#Design matrix
design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "High_vs_Low")


class(expr_mat)       # should be "matrix" or "data.frame"
str(expr_mat)         # check column types

head(expr_data)
str(expr_data)



library(data.table)

# Pivot to wide format
expr_mat <- dcast(expr_data[, .(Gene, Sample_ID, Expression_log2)],
                  Gene ~ Sample_ID,
                  value.var = "Expression_log2")

# Set rownames and remove Gene column
rownames(expr_mat) <- expr_mat$Gene
expr_mat$Gene <- NULL

# Convert to numeric matrix
expr_mat <- apply(expr_mat, 2, as.numeric)
expr_mat <- t(expr_mat) # apply() returns samples x genes, so transpose back if needed
rownames(expr_mat) <- colnames(expr_data)[2:ncol(expr_data)] # fix rownames




#limma DEG run
fit <- lmFit(expr_mat, design)
fit <- eBayes(fit)

deg_results <- topTable(
  fit,
  coef = "High_vs_Low",
  number = Inf,
  adjust.method = "fdr"
)

deg_results$Gene <- rownames(deg_results)



#DEG classification (Up / Down / NS)
deg_results$DEG_status <- "Non Significaant"

deg_results$DEG_status[
  deg_results$logFC > 1 & deg_results$adj.P.Val < 0.05
] <- "Up regulated"

deg_results$DEG_status[
  deg_results$logFC < -1 & deg_results$adj.P.Val < 0.05
] <- "Down regulated"



#Final sanity check
table(deg_results$DEG_status)
head(deg_results)










                                    ################### substep: 4.4 ############


#packages
library(ggplot2)
library(ggrepel)


#Prepare plotting dataframe
volcano_df <- deg_results
volcano_df$negLogFDR <- -log10(volcano_df$adj.P.Val)



#Select genes to label(Top 10 Up + Top 10 Down only)
label_up <- volcano_df[
  volcano_df$DEG_status == "Up regulated",
][order(volcano_df$adj.P.Val), ][1:10, ]

label_down <- volcano_df[
  volcano_df$DEG_status == "Down regulated",
][order(volcano_df$adj.P.Val), ][1:10, ]

label_genes <- rbind(label_up, label_down)
label_genes$Gene





#Volcano plot (clean & meaningful)
library(ggplot2)
library(ggrepel)


#Base volcano dataframe
volcano_df$Gene <- rownames(volcano_df)



#Base volcano plot (NO labels yet)
p <- ggplot(
  volcano_df,
  aes(x = logFC, y = -log10(adj.P.Val), color = DEG_status)
) +
  geom_point(size = 0.8) +
  scale_color_manual(
    values = c(
      "Up regulated"   = "#D62728",
      "Down regulated" = "#1F77B4",
      "Non Significaant" = "green"
    )
  ) +
  theme_classic(base_size = 8) +
  theme(
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  )





#ADD gene name labels
p <- p +
  geom_text_repel(
    data = label_genes,
    aes(label = Gene),
    size = 2.2,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    segment.size = 0.2,
    show.legend = FALSE
  )


#View plot
print(p)

# Save as TIFF (high-resolution)
ggsave(
  filename = "Figure2C_DEG_Volcano.tiff",
  plot = p,
  device = "tiff",
  width = 4,
  height = 4,
  units = "in",
  dpi = 1000,
  compression = "lzw"
)

#Save EPS
ggsave(
  "Figure2C_DEG_Volcano.eps",
  plot = p,
  device = "eps",
  dpi = 1200,
  width = 4,
  height = 4
)




#for TIFF format
library(ggplot2)
library(ggrepel)

# Ensure volcano_df has a "Gene" column
volcano_df$Gene <- rownames(volcano_df)

# Fix factor levels to match scale_color_manual
volcano_df$DEG_status <- factor(
  volcano_df$DEG_status,
  levels = c("Upregulated", "Downregulated", "Non Significaant")
)

# Base volcano plot
p <- ggplot(
  volcano_df,
  aes(x = logFC, y = -log10(adj.P.Val), color = DEG_status)
) +
  geom_point(size = 0.8) +
  scale_color_manual(
    values = c(
      "Upregulated"      = "#800000",  # Maroon
      "Downregulated"    = "#800080",  # Purple
      "Non Significaant" = "green"
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  )

# Add gene labels
p <- p +
  geom_text_repel(
    data = label_genes,
    aes(label = Gene),
    size = 2.5,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    segment.size = 0.2,
    show.legend = FALSE
  )

# View plot
print(p)

# Save as TIFF (high-resolution)
ggsave(
  filename = "Figure2C_DEG_Volcano.tiff",
  plot = p,
  device = "tiff",
  width = 4,
  height = 4,
  units = "in",
  dpi = 1000,
  compression = "lzw"
)








                                    ############## substep: 4.5 ############
#“Only genes that are both prognostic (survival-associated) and biologically dysregulated were retained.”



#Load survival genes
cox_sig <- read.csv("Significant.csv", stringsAsFactors = FALSE)
surv_genes <- unique(cox_sig$Gene)
length(surv_genes)



#Extract DEG genes (only significant ones)
deg_sig <- volcano_df[
  volcano_df$DEG_status != "Non Significaant",
]

deg_genes <- unique(deg_sig$Gene)
length(deg_genes)



#INTERSECTION (core result)
surv_deg_genes <- intersect(surv_genes, deg_genes)
length(surv_deg_genes)



#Create final survival-relevant DEG table
final_biomarkers <- deg_sig[
  deg_sig$Gene %in% surv_deg_genes,
]

final_biomarkers <- merge(
  final_biomarkers,
  cox_sig,
  by = "Gene"
)

head(final_biomarkers)



#Save table
write.csv(
  final_biomarkers,
  "Table4_Survival_Relevant_DEGs.csv",
  row.names = FALSE
)





                               ########### Venn Diagram #############

#“Overlap between differentially expressed genes and survival-associated genes.”



#Prepare gene sets
deg_genes <- unique(deg_sig$Gene)
surv_genes <- unique(cox_sig$Gene)
length(intersect(deg_genes, surv_genes))


#Draw Venn diagram
library(VennDiagram)
library(grid)

venn.plot <- venn.diagram(
  x = list(
    "DEGs" = deg_genes,
    "Survival genes" = surv_genes
  ),
  filename = NULL,
  fill = c("#4DBBD5", "#E64B35"),  # journal colors
  alpha = 1,                      # ❗ EPS-safe
  cex = 1.5,
  cat.cex = 1.3,
  fontfamily = "Helvetica",
  cat.fontfamily = "Helvetica",
  margin = 0.1
)


postscript(
  "Figure2D_DEG_Survival_Venn.eps",
  width = 4,
  height = 4,
  family = "Helvetica",
  horizontal = FALSE,
  onefile = FALSE,
  paper = "special"
)

grid.draw(venn.plot)
dev.off()




library(dplyr)

top_genes_ordered <- top_genes %>%
  mutate(abs_logFC = abs(logFC)) %>%
  arrange(desc(abs_logFC), adj.P.Val)

top_genes_ordered




# Convert to uppercase
surv_genes <- toupper(trimws(surv_genes))
deg_genes  <- toupper(trimws(deg_genes))

common_genes <- intersect(surv_genes, deg_genes)
common_genes
length(common_genes)

head(surv_genes)
head(deg_genes)







                                            ###  step :5 LASSO Cox Feature Selection (Final Biomarker Construction) 


survival_DEGs <- read.csv("survival_DEGs.csv")

# Save dt dataset as CSV
write.csv(dt, "dt.csv", row.names = FALSE)
dt <- read.csv("dt.csv")


library(data.table)
library(dplyr)
library(tidyr)
library(survival)
library(glmnet)


#Input gene list (198 genes)
genes_use <- unique(trimws(survival_DEGs$Gene))
length(genes_use)   # MUST be 198



#Expression matrix (LONG → WIDE)

expr_long <- dt[, c("Sample_ID", "Gene", "Expression_log2")]

# keep only 198 survival-DEG genes
expr_long <- expr_long %>%
  filter(Gene %in% genes_use)


# create wide
expr_wide <- expr_long %>%
  pivot_wider(
    names_from  = Gene,
    values_from = Expression_log2
  )



#Final expression matrix for glmnet

# Sample_ID rownames
expr_mat <- as.data.frame(expr_wide)
rownames(expr_mat) <- expr_mat$Sample_ID
expr_mat$Sample_ID <- NULL

# convert to matrix
x <- as.matrix(expr_mat)

dim(x)



#Survival object (from SAME dt)
surv_data <- dt[, c("Sample_ID", "OS_days", "OS_event")]
surv_data <- unique(surv_data)

# align with expression matrix
surv_data <- surv_data[match(rownames(x), surv_data$Sample_ID), ]

# create survival object
y <- Surv(
  time  = surv_data$OS_days,
  event = surv_data$OS_event
)



#LASSO Cox regression (CORE STEP )

set.seed(123)

cvfit <- cv.glmnet(
  x = x,
  y = y,
  family = "cox",
  alpha  = 1,
  nfolds = 10
)

lambda_min <- cvfit$lambda.min
lambda_1se <- cvfit$lambda.1se





#Final biomarker gene selection (Table 5)

coef_min <- coef(cvfit, s = "lambda.min")
idx <- which(coef_min != 0)

lasso_table <- data.frame(
  Gene        = rownames(coef_min)[idx],
  Coefficient = as.numeric(coef_min[idx]),
  Risk_type   = ifelse(coef_min[idx] > 0, "Risk", "Protective")
)

lasso_table


#save final Biomarker genes
write.csv(
  lasso_table,
  "Table5_LASSO_signature.csv",
  row.names = FALSE
)

#Figure 3A — CV curve

plot(cvfit)
abline(v = log(lambda_min), col = "red", lty = 2)
abline(v = log(lambda_1se), col = "blue", lty = 2)


#Advanced (Nature Style)
library(ggplot2)

cv_df <- data.frame(
  log_lambda = log(cvfit$lambda),
  cvm = cvfit$cvm,
  cvsd = cvfit$cvsd
)

ggplot(cv_df, aes(x = log_lambda, y = cvm)) +
  geom_line(color = "black", linewidth = 0.9) +
  geom_ribbon(
    aes(ymin = cvm - cvsd, ymax = cvm + cvsd),
    fill = "grey80",
    alpha = 0.6
  ) +
  geom_vline(xintercept = log(lambda_min),
             color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = log(lambda_1se),
             color = "navy", linetype = "dashed", linewidth = 0.8) +
  annotate("text",
           x = log(lambda_min),
           y = min(cv_df$cvm),
           label = "λ.min",
           color = "red",
           angle = 90,
           vjust = -0.5,
           size = 4) +
  annotate("text",
           x = log(lambda_1se),
           y = min(cv_df$cvm),
           label = "λ.1se",
           color = "navy",
           angle = 90,
           vjust = -0.5,
           size = 4) +
  labs(
    x = "log(λ)",
    y = "Partial likelihood deviance",
    title = "(A) LASSO Cox cross-validation curve"
  ) +
  theme_classic(base_size = 12)




#EPS format


postscript(
  file = "Figure_3A_LASSO_CV_Advanced.eps",
  width = 7,
  height = 5.5,
  horizontal = FALSE,
  onefile = FALSE,
  paper = "special"
)

ggplot(cv_df, aes(x = log_lambda, y = cvm)) +
  
  #  Uncertainty band (±1 SE) 
  geom_ribbon(
    aes(ymin = cvm - cvsd, ymax = cvm + cvsd),
    fill = "grey85",
    alpha = 0.8
  ) +
  
  #  Mean CV curve 
  geom_line(
    linewidth = 1.1,
    color = "black"
  ) +
  
  #  lambda.min 
  geom_vline(
    xintercept = log(lambda_min),
    linetype = "dashed",
    linewidth = 0.8,
    color = "red3"
  ) +
  
  #  lambda.1se 
  geom_vline(
    xintercept = log(lambda_1se),
    linetype = "dashed",
    linewidth = 0.8,
    color = "navy"
  ) +
  
  # lambda labels 
  annotate(
    "text",
    x = log(lambda_min),
    y = min(cv_df$cvm),
    label = paste0("lambda.min = ", signif(lambda_min, 3)),
    hjust = -0.1,
    vjust = -1,
    size = 3,
    color = "red3"
  ) +
  
  annotate(
    "text",
    x = log(lambda_1se),
    y = min(cv_df$cvm),
    label = paste0("lambda.1se = ", signif(lambda_1se, 3)),
    hjust = -0.1,
    vjust = -3,
    size = 3,
    color = "navy"
  ) +
  
  # --- Axis & labels ---
  scale_x_reverse() +
  
  labs(
    title = "LASSO Cox regression with 10-fold cross-validation",
    x = "log(Lambda)",
    y = "Partial likelihood deviance"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0, size = 13),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )

dev.off()











#Figure 3B — Coefficient path

plot(cvfit$glmnet.fit, xvar = "lambda", label = TRUE)
abline(v = log(lambda_min), col = "red", lty = 2)
postscript(
  "Figure_3B_LASSO_Coefficient_Path_GeneLabeled.eps",
  width = 7,
  height = 5.5,
  horizontal = FALSE,
  paper = "special"
)





#Advanced (nature style)

library(ggrepel)

fit <- cvfit$glmnet.fit
beta <- as.matrix(fit$beta)
lambda_seq <- log(fit$lambda)

coef_df <- as.data.frame(beta)
coef_df$Gene <- rownames(coef_df)

coef_long <- coef_df %>%
  tidyr::pivot_longer(
    cols = -Gene,
    names_to = "LambdaIndex",
    values_to = "Coefficient"
  )

coef_long$log_lambda <- rep(lambda_seq, each = nrow(beta))

# final selected genes
final_genes <- lasso_table$Gene

coef_long$Group <- ifelse(
  coef_long$Gene %in% final_genes,
  ifelse(lasso_table$Coefficient[match(coef_long$Gene, lasso_table$Gene)] > 0,
         "Risk", "Protective"),
  "Other"
)

ggplot(coef_long, aes(x = log_lambda, y = Coefficient, group = Gene)) +
  geom_line(
    aes(color = Group),
    linewidth = 0.8,
    alpha = 0.8
  ) +
  scale_color_manual(
    values = c(
      "Risk" = "red3",
      "Protective" = "steelblue",
      "Other" = "grey80"
    )
  ) +
  geom_vline(xintercept = log(lambda_min),
             linetype = "dashed", color = "black") +
  geom_text_repel(
    data = subset(coef_long,
                  log_lambda == min(log_lambda) &
                    Gene %in% final_genes),
    aes(label = Gene),
    size = 3,
    box.padding = 0.3,
    max.overlaps = 20
  ) +
  labs(
    x = "log(λ)",
    y = "LASSO coefficient",
    title = "(B) LASSO coefficient path"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank())



#EPS format

library(tidyr)
library(ggrepel)

# glmnet model
fit <- cvfit$glmnet.fit
beta <- as.matrix(fit$beta)
lambda <- fit$lambda

coef_df <- data.frame(Gene = rownames(beta), beta)

coef_long <- coef_df |>
  pivot_longer(
    cols = -Gene,
    names_to = "idx",
    values_to = "Coefficient"
  )

coef_long$idx <- as.numeric(gsub("V", "", coef_long$idx))
coef_long$Lambda <- lambda[coef_long$idx]

# final selected genes
coef_min <- as.matrix(coef(cvfit, s = "lambda.min"))
final_genes <- rownames(coef_min)[coef_min[, 1] != 0]

coef_long$Status <- ifelse(
  coef_long$Gene %in% final_genes,
  "Selected",
  "Other"
)

postscript(
  "Figure_3B_LASSO_Coefficient_Path_GeneLabeled.eps",
  width = 7,
  height = 5.5,
  horizontal = FALSE,
  paper = "special"
)

ggplot(coef_long,
       aes(x = log(Lambda),
           y = Coefficient,
           group = Gene)) +
  
  # background genes
  geom_line(
    data = subset(coef_long, Status == "Other"),
    color = "grey80",
    linewidth = 0.4
  ) +
  
  # selected genes
  geom_line(
    data = subset(coef_long, Status == "Selected"),
    aes(color = Gene),
    linewidth = 1.2,
    show.legend = FALSE
  ) +
  
  # lambda.min line
  geom_vline(
    xintercept = log(cvfit$lambda.min),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  # gene labels (ONLY selected genes)
  geom_text_repel(
    data = subset(
      coef_long,
      Status == "Selected" &
        abs(Lambda - cvfit$lambda.min) ==
        min(abs(Lambda - cvfit$lambda.min))
    ),
    aes(label = Gene),
    size = 3,
    max.overlaps = Inf,
    segment.color = "grey50"
  ) +
  
  scale_x_reverse() +
  
  labs(
    title = "LASSO Cox coefficient paths",
    x = "log(Lambda)",
    y = "Coefficient"
  ) +
  
  theme_classic(base_size = 12)

dev.off()












# Final gene signature Barplot

lasso_table$Gene <- factor(
  lasso_table$Gene,
  levels = lasso_table$Gene[order(lasso_table$Coefficient)]
)

ggplot(lasso_table,
       aes(x = Coefficient, y = Gene, fill = Risk_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(
    values = c("Risk" = "red3", "Protective" = "blue")
  ) +
  labs(
    x = "LASSO coefficient",
    y = "",
    title = " "
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")







#EPS format




library(ggplot2)

# order genes by coefficient
lasso_table$Gene <- factor(
  lasso_table$Gene,
  levels = lasso_table$Gene[order(lasso_table$Coefficient)]
)

postscript(
  "Figure_3C_LASSO_Gene_Coefficients_withValues.eps",
  width = 8,
  height = 6,
  horizontal = FALSE,
  paper = "special"
)

ggplot(lasso_table,
       aes(x = Gene,
           y = Coefficient,
           fill = Risk_type)) +
  
  geom_bar(
    stat = "identity",
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  geom_text(
    aes(label = round(Coefficient, 2)),
    hjust = ifelse(lasso_table$Coefficient > 0, -0.15, 1.15),
    size = 3
  ) +
  
  scale_fill_manual(
    values = c(
      "Risk" = "#D73027",
      "Protective" = "#1A9850"
    )
  ) +
  
  coord_flip() +   # 👈 highly recommended for readability
  
  labs(
    title = "LASSO-selected prognostic gene signature",
    x = "Genes",
    y = "LASSO Cox coefficient"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

dev.off()













                                          ########## step : 6 : Multivariate Cox (Clinical Adjustment) Independent Prognostic Validation






# Age variable merge kora
clin_data <- merge(clin_data, expr_data[, c("Sample_ID", "age")], by = "Sample_ID", all.x = TRUE)
# Unique Sample_ID check
length(unique(clin_data$Sample_ID))

# Unique samples subset
clin_data <- clin_data[!duplicated(clin_data$Sample_ID), ]
nrow(clin_data)  # 1185


#loaded 22 signature gene
signature <- read.csv("signature.csv")
str(signature)







                                                                    #Risk Score Construction


library(data.table)
library(dplyr)

# ensure data.table
expr_deg <- as.data.table(expr_deg)
signature <- as.data.table(signature)

# keep only 22 signature genes
expr_sig <- expr_deg[Gene %in% signature$Gene]

# merge coefficient
expr_sig <- merge(
  expr_sig,
  signature[, .(Gene, Coefficient)],
  by = "Gene"
)

# calculate gene-wise contribution
expr_sig[, weighted_expr := Expression_log2 * Coefficient]

# sum per sample → Risk Score
risk_df <- expr_sig[, .(
  RiskScore = sum(weighted_expr, na.rm = TRUE)
), by = Sample_ID]

# check
head(risk_df)







                                                       # Clinical Merge + Proper Recoding


final_df <- merge(
  clin_data,
  risk_df,
  by = "Sample_ID"
)


#save in csv file
write.csv(
  final_df,
  file = "final_df.csv",
  row.names = FALSE
)


  





                                               # forest plot 


install.packages(c("rms", "survminer", "ggDCA", "foreign"))




library(survival)
library(survminer)
library(dplyr)
library(rms)      
library(ggDCA)    


final_df <- read.csv("final_df.csv")

# 2. Check structure
str(final_df)


# Survival & Age correction
final_df$OS_event <- as.numeric(final_df$OS_event)
final_df$OS_days  <- as.numeric(final_df$OS_days)
final_df$age      <- as.numeric(final_df$age)

#  Tumor Grade (Ref = G1) 
final_df$tumor_grade <- factor(final_df$tumor_grade, levels = c("G1", "G2", "G3"))
final_df$tumor_grade <- relevel(final_df$tumor_grade, ref = "G1")

# ER Status (Ref = Negative) 
# Assuming 0=Negative, 1=Positive in your raw data
final_df$ER_status <- factor(final_df$ER_status, levels = c(0, 1), labels = c("Negative", "Positive"))
final_df$ER_status <- relevel(final_df$ER_status, ref = "Negative")

#  PR Status (Ref = Negative) 
final_df$PR_status <- factor(final_df$PR_status, levels = c(0, 1), labels = c("Negative", "Positive"))
final_df$PR_status <- relevel(final_df$PR_status, ref = "Negative")

#  HER2 Status (Ref = Negative) 
final_df$HER2_status <- factor(final_df$HER2_status, levels = c(0, 1), labels = c("Negative", "Positive"))
final_df$HER2_status <- relevel(final_df$HER2_status, ref = "Negative")

# Ki67 Status (Ref = Negative/Low) 
final_df$Ki67_status <- factor(final_df$Ki67_status, levels = c(0, 1), labels = c("Negative", "Positive"))
final_df$Ki67_status <- relevel(final_df$Ki67_status, ref = "Negative")



# STEP 3: MULTIVARIATE COX REGRESSION

cox_multi <- coxph(
  Surv(OS_days, OS_event) ~ RiskScore + age + tumor_grade + ER_status + PR_status + HER2_status + Ki67_status,
  data = final_df
)

# Check the new result (Is G3 HR > 1 now?)
summary(cox_multi)




uni_grade <- coxph(Surv(OS_days, OS_event) ~ tumor_grade, data = final_df)
summary(uni_grade)


# This creates a professional plot directly from the model
ggforest(
  cox_multi, 
  data = final_df, 
  main = "Hazard Ratio", 
  cpositions = c(0.02, 0.22, 0.4), 
  fontsize = 1.0, 
  refLabel = "Reference", 
  noDigits = 2
)

# Save the plot
ggsave("Figure_Multivariate_Forest.pdf", width = 10, height = 8)









library(survival)
library(broom)
library(dplyr)

covariates <- c("RiskScore", "age", "tumor_grade", "ER_status", "PR_status", "HER2_status", "Ki67_status")

# Run Univariate Analysis Loop (Proof of 1.73)
univ_results <- lapply(covariates, function(x){
  # Create formula dynamically
  f <- as.formula(paste("Surv(OS_days, OS_event) ~", x))
  
  # Run Cox Model
  model <- coxph(f, data = final_df)
  
  # Extract results cleanly
  res <- tidy(model, exponentiate = TRUE, conf.int = TRUE)
  
  # Keep only necessary columns
  res <- res[, c("term", "estimate", "p.value", "conf.low", "conf.high")]
  
  # Add a column to identify variable
  res$Variable <- x
  return(res)
})

# Combine all univariate results
univ_table <- do.call(rbind, univ_results)

# Format Univariate Columns
univ_table$Univariate_HR <- sprintf("%.2f (%.2f-%.2f)", univ_table$estimate, univ_table$conf.low, univ_table$conf.high)
univ_table$Univariate_P  <- sprintf("%.3f", univ_table$p.value)

#  Run Multivariate Analysis (Existing Result)
multi_model <- coxph(Surv(OS_days, OS_event) ~ RiskScore + age + tumor_grade + ER_status + PR_status + HER2_status + Ki67_status, data = final_df)

multi_res <- tidy(multi_model, exponentiate = TRUE, conf.int = TRUE)

# Format Multivariate Columns
multi_res$Multivariate_HR <- sprintf("%.2f (%.2f-%.2f)", multi_res$estimate, multi_res$conf.low, multi_res$conf.high)
multi_res$Multivariate_P  <- sprintf("%.3f", multi_res$p.value)

# 4. Merge Both into One Final Table

# Align rows by 'term'
final_table_2 <- merge(univ_table[, c("term", "Univariate_HR", "Univariate_P")], 
                       multi_res[, c("term", "Multivariate_HR", "Multivariate_P")], 
                       by = "term", all = TRUE)

# Clean up Term names for Paper
final_table_2$term <- gsub("tumor_gradeG2", "Grade (G2 vs G1)", final_table_2$term)
final_table_2$term <- gsub("tumor_gradeG3", "Grade (G3 vs G1)", final_table_2$term)
final_table_2$term <- gsub("RiskScore", "Risk Score", final_table_2$term)

# Save and View
print(final_table_2)

write.csv(final_table_2, "Table_2_Uni_vs_Multi_Analysis.csv", row.names = FALSE)







# generate plot 

library(survival)
library(broom)
library(dplyr)
library(ggplot2)
if(!require(patchwork)) install.packages("patchwork")
library(patchwork)


# data processing 
covariates <- c("RiskScore", "age", "tumor_grade", "ER_status", "PR_status", "HER2_status", "Ki67_status")

#  Univariate 
univ_list <- lapply(covariates, function(x){
  f <- as.formula(paste("Surv(OS_days, OS_event) ~", x))
  model <- coxph(f, data = final_df)
  res <- tidy(model, exponentiate = TRUE, conf.int = TRUE)
  res <- res[, c("term", "estimate", "conf.low", "conf.high", "p.value")]
  res$Type <- "Univariate"
  return(res)
})
univ_df <- do.call(rbind, univ_list)

#  Multivariate 
multi_model <- coxph(Surv(OS_days, OS_event) ~ RiskScore + age + tumor_grade + 
                       ER_status + PR_status + HER2_status + Ki67_status, data = final_df)
multi_df <- tidy(multi_model, exponentiate = TRUE, conf.int = TRUE)
multi_df <- multi_df[, c("term", "estimate", "conf.low", "conf.high", "p.value")]
multi_df$Type <- "Multivariate"

# Combine 
plot_data <- rbind(univ_df, multi_df)




# make label and text

plot_data$term <- recode(plot_data$term,
                         RiskScore = "Risk Score",
                         age = "Age",
                         tumor_gradeG2 = "Grade: G2 vs G1",
                         tumor_gradeG3 = "Grade: G3 vs G1",
                         ER_statusPositive = "ER(+/-)",
                         PR_statusPositive = "PR(+/-)",
                         HER2_statusPositive = "HER2(+/-)",
                         Ki67_statusPositive = "Ki67(+/-)")

plot_data$term <- factor(plot_data$term, levels = rev(unique(plot_data$term)))

# HR label add 
plot_data$hr_text <- sprintf("%.2f (%.2f-%.2f)", plot_data$estimate, plot_data$conf.low, plot_data$conf.high)

# P-value add
plot_data$p_text <- sapply(plot_data$p.value, function(x) {
  if (x < 0.001) sprintf("%.2e", x) else sprintf("%.3f", x)
})



# variables name in left side 
p1 <- ggplot(subset(plot_data, Type == "Univariate"), aes(y = term)) +
  geom_text(aes(x = 0, label = term), hjust = 0, fontface = "bold", size = 3.8) +
  labs(title = "Variable") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0, size = 11),
    plot.margin = margin(t=0, b=0, l=0, r=0)
  ) +
  xlim(0, 1)

#  HR value text
# x = 0 and group = Type 
p2 <- ggplot(plot_data, aes(y = term, label = hr_text, color = Type, group = Type)) +
  geom_text(aes(x = 0), position = position_dodge(width = 0.8), size = 3.5, fontface = "bold") +
  labs(title = "HR (95% CI)") +
  scale_color_manual(values = c("Multivariate" = "#00468B", "Univariate" = "#ED0000")) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    legend.position = "none"
  ) +
  xlim(-0.5, 0.5) #space for text

# Part 3: main graph
p3 <- ggplot(plot_data, aes(y = term, color = Type, group = Type)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray80") +
  geom_pointrange(aes(x = estimate, xmin = conf.low, xmax = conf.high), 
                  position = position_dodge(width = 0.8), size = 0.6, shape = 15) +
  scale_x_log10(breaks = c(0.5, 1, 2, 5)) +
  scale_color_manual(values = c("Multivariate" = "#00468B", "Univariate" = "#ED0000")) +
  labs(title = "Hazard Ratio Plot") +
  theme_bw() +
  theme(
    axis.text.y = element_blank(), 
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    legend.position = "bottom"
  )

#  P-Value (right side)
#  x = 0 and group = Type 
p4 <- ggplot(plot_data, aes(y = term, label = p_text, color = Type, group = Type)) +
  geom_text(aes(x = 0), position = position_dodge(width = 0.8), size = 3.5) +
  labs(title = "P Value") +
  scale_color_manual(values = c("Multivariate" = "#00468B", "Univariate" = "#ED0000")) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
    legend.position = "none"
  ) +
  xlim(-0.5, 0.5)


# Combine (p1 + p2 + p3 + p4)
final_plot <- p1 + p2 + p3 + p4 + 
  plot_layout(nrow = 1, widths = c(1.5, 2, 3, 1))




print(final_plot)
  


# Save
ggsave(
  filename = "Figure_Final_Table_Graph_Fixed.eps",
  plot = final_plot,
  device = "eps",
  width = 15,
  height = 7,
  dpi = 1000

)










                                     # step:7: Risk Stratification & Survival Performance Validation 



dt <- read.csv("dt.csv")



#High-risk vs Low-risk Stratification (Median cut)
# Median cutoff
cutoff <- median(final_df$RiskScore, na.rm = TRUE)

final_df$risk_group <- ifelse(
  final_df$RiskScore > cutoff,
  "High-risk",
  "Low-risk"
)

final_df$risk_group <- factor(
  final_df$risk_group,
  levels = c("Low-risk", "High-risk")
)

table(final_df$risk_group)





   
                                                #Kaplan–Meier Survival Curve (MAIN FIGURE)
 

library(survival)
library(survminer)

surv_obj <- Surv(time = final_df$OS_days, event = final_df$OS_event)

km_fit <- survfit(surv_obj ~ risk_group, data = final_df)

km_plot <- ggsurvplot(
  km_fit,
  data = final_df,
  risk.table = TRUE,             # Add number at risk table below KM curve
  pval = TRUE,                   # Show log-rank p-value
  conf.int = TRUE,               # Show confidence interval
  xlab = "Days",
  ylab = "Overall Survival Probability",
  legend.title = "Risk Group",
  legend.labs = c("Low-risk", "High-risk"),
  palette = c("#1B9E77","#D95F02"),
  risk.table.height = 0.25,      
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,     
  break.time.by = 200,           
  surv.median.line = "hv",       
  tables.theme = theme_cleantable()
)


# Display
km_plot







                                        # ROC curve #

install.packages("timeROC")
# Required library
library(timeROC)

# Run time-dependent ROC
roc_res <- timeROC(
  T = final_df$OS_days,          # survival time
  delta = final_df$OS_event,     # event (1=death, 0=alive)
  marker = final_df$RiskScore,   # continuous risk score
  cause = 1,                     # event of interest
  times = c(365, 1095, 1825),    # 1y, 3y, 5y in days
  iid = TRUE
)

# Plot ROC curves
plot(
  roc_res,
  time = 365,
  col = "#2C7BB6",
  lwd = 2,
  title = FALSE
)
plot(
  roc_res,
  time = 1095,
  col = "#FDAE61",
  lwd = 2,
  add = TRUE
)
plot(
  roc_res,
  time = 1825,
  col = "#D7191C",
  lwd = 2,
  add = TRUE
)

# Add legend with AUC values
legend(
  "bottomright",
  legend = c(
    paste0("1-year AUC = ", round(roc_res$AUC[1], 3)),
    paste0("3-year AUC = ", round(roc_res$AUC[2], 3)),
    paste0("5-year AUC = ", round(roc_res$AUC[3], 3))
  ),
  col = c("#2C7BB6", "#FDAE61", "#D7191C"),
  lwd = 2,
  bty = "n"
)

# Optional: create AUC table for manuscript
roc_table <- data.frame(
  Time = c("1-year", "3-year", "5-year"),
  AUC = round(roc_res$AUC, 3)
)

write.csv(
  roc_table,
  file = "roc_table.csv",
  row.names = FALSE
)






                                      ### Risk Score Distribution Plot
signature <- read.csv("signature.csv")

library(dplyr)
library(ggplot2)

# 1. Order patients by RiskScore
risk_ord <- final_df %>%
  arrange(RiskScore) %>%
  mutate(Order = row_number())

cut_pos <- sum(risk_ord$RiskGroup == "Low-risk")

# 2. Risk score distribution plot
p_A <- ggplot(risk_ord, aes(x = Order, y = RiskScore, color = RiskGroup)) +
  geom_step(linewidth = 1.2) +
  geom_vline(
    xintercept = cut_pos,
    linetype = "dashed",
    linewidth = 0.8,
    color = "black"
  ) +
  scale_color_manual(
    values = c("Low-risk" = "#2C7BB6", "High-risk" = "#D7191C")
  ) +
  labs(
    y = "Risk score",
    x = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
#check
table(risk_ord$RiskGroup)





#Survival Status (Panel B)
# Add survival status labels
risk_ord$Status <- ifelse(risk_ord$OS_event == 1, "Dead", "Alive")

# Survival status plot
p_B <- ggplot(risk_ord, aes(x = Order, y = OS_days, color = Status)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("Alive" = "#2C7BB6", "Dead" = "#D7191C")) +
  labs(
    y = "Survival time (days)",
    x = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )





final_df <- read.csv("final_df.csv")






#Heatmap plot
signature <- read.csv("signature.csv")
final_df <- read.csv("final_df.csv")
str(final_df)
library(dplyr)

risk_ord <- final_df %>%
  arrange(RiskScore) %>%
  mutate(
    Order = row_number(),
    RiskGroup = risk_group   # <-- unify name
  )

write.csv(
  risk_ord,
  file = "risk_ord.csv",
  row.names = FALSE
)


table(risk_ord$RiskGroup)




# main point of Heatmap
library(dplyr)
library(tidyr)
library(pheatmap)
library(tibble)

sig_genes <- signature$Gene

expr_sig <- expr_deg %>%
  filter(Gene %in% sig_genes)

expr_mat_sig <- expr_sig %>%
  pivot_wider(
    names_from  = Sample_ID,
    values_from = Expression_log2
  ) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

expr_mat_sig <- expr_mat_sig[, risk_ord$Sample_ID]

expr_mat_z <- t(scale(t(expr_mat_sig)))
expr_mat_z[is.na(expr_mat_z)] <- 0

annotation_col <- data.frame(
  RiskGroup = risk_ord$RiskGroup,
  row.names = risk_ord$Sample_ID
)

pheatmap(
  expr_mat_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
  annotation_col = annotation_col,
  annotation_colors = list(
    RiskGroup = c("Low-risk" = "#2C7BB6", "High-risk" = "#D7191C")
  ),
  fontsize_row = 8,
  border_color = NA
)







library(ggpubr)
library(grid)
library(gridExtra)



# Draw heatmap and capture as grob
p_C <- pheatmap(
  expr_mat_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
  annotation_col = annotation_col,
  annotation_colors = list(
    RiskGroup = c("Low-risk" = "#2C7BB6", "High-risk" = "#D7191C")
  ),
  fontsize_row = 8,
  border_color = NA,
  silent = TRUE
)

heatmap_grob <- p_C$gtable



#A + B vertical combine
upper_panels <- ggarrange(
  p_A, p_B,
  ncol = 1,
  heights = c(1, 1),
  labels = c("A", "B"),
  font.label = list(size = 14, face = "bold")
)




#(A+B) + C final combine
final_fig <- arrangeGrob(
  upper_panels,
  heatmap_grob,
  ncol = 1,
  heights = c(2, 3)
)



#SAVE AS TIFF (Journal standard)

tiff(
  filename = "Figure4C_RiskScore_Survival_Heatmap.tiff",
  width = 180,
  height = 220,
  units = "mm",
  res = 1000,
  compression = "lzw"
)

grid.draw(final_fig)
dev.off()



library(grid)
cairo_ps(
  filename = "Figure4C_RiskScore_Survival_Heatmap.eps",
  width = 7.1,     # inches (≈ 180 mm)
  height = 8.7,    # inches (≈ 220 mm)
  fallback_resolution = 1200
)

grid.draw(final_fig)
dev.off()



write.csv(
  risk_ord,
  file = "risk_ord.csv",
  row.names = FALSE
)






  










                       #step: 8: Functional Biological Interpretation



str(signature)
str(final_df)
str(risk_ord)
str(expr_wide)



#Biomarker Gene Preparation & Mapping

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

## Step : Map gene symbols to Entrez IDs
gene_map <- bitr(
  signature$Gene,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

## Step 8.0.2: Merge mapping with signature info
sig_entrez <- signature %>%
  inner_join(gene_map, by = c("Gene" = "SYMBOL"))

## Check results
nrow(sig_entrez)        # should be 20
sig_entrez



#save in csv
write.csv(
  sig_entrez,
  file = "sig_entrez.csv",
  row.names = FALSE
)





#STEP 8A — GO Enrichment Analysis (BP / CC / MF)
gene_list <- sig_entrez$ENTREZID




#GO enrichment (BP, CC, MF)
library(clusterProfiler)

## Biological Process
ego_bp <- enrichGO(
  gene          = gene_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)





## Cellular Component
ego_cc <- enrichGO(
  gene          = gene_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

## Molecular Function
ego_mf <- enrichGO(
  gene          = gene_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)




#Combine GO results (for plotting)

go_list <- list(
  BP = as.data.frame(ego_bp),
  CC = as.data.frame(ego_cc),
  MF = as.data.frame(ego_mf)
)

go_list <- read.csv("go_list.csv")
# empty remove and Ontology label add

go_list <- lapply(names(go_list), function(ont) {
  df <- go_list[[ont]]
  if (nrow(df) > 0) {
    df$ONTOLOGY <- ont
    return(df)
  }
})

go_list <- do.call(rbind, go_list)
str(go_list)




                                                   #GO plot 
library(ggplot2)

p_go <- ggplot(
  go_list,
  aes(
    x = GeneRatio,
    y = reorder(Description, GeneRatio),
    color = p.adjust,
    size = Count
  )
) +
  geom_point() +
  facet_grid(ONTOLOGY ~ ., scales = "free_y") +
  theme_bw() +
  labs(
    x = "Gene Ratio",
    y = "GO Term",
    color = "Adjusted P-value",
    size = "Gene Count"
  ) +
  scale_color_gradient(low = "pink", high = "red")  # 🔴 red gradient

p_go




#maroon 
library(ggplot2)

p_go <- ggplot(
  go_list,
  aes(
    x = GeneRatio,
    y = reorder(Description, GeneRatio),
    color = p.adjust,
    size = Count
  )
) +
  geom_point() +
  facet_grid(ONTOLOGY ~ ., scales = "free_y") +
  theme_bw() +
  labs(
    x = "Gene Ratio",
    y = "GO Term",
    color = "Adjusted P-value",
    size = "Gene Count"
  ) +
  scale_color_gradientn(colors = c("red", "darkred", "maroon"))  # 🔴➡️🟤 maroon

p_go



write.csv(
  go_list,
  file = "go_list.csv",
  row.names = FALSE
)


str(go_list)






                                           #KEGG plot 

#KEGG enrichment execution
library(clusterProfiler)
library(org.Hs.eg.db)

# (same vector you used for GO)

write.csv(
  gene_map,
  file = "gene_map.csv",
  row.names = FALSE
)

#Final KEGG input vector
entrez_genes <- unique(gene_map$ENTREZID)
length(entrez_genes)



#KEGG enrichment run

library(clusterProfiler)
library(org.Hs.eg.db)

ekegg <- enrichKEGG(
  gene          = entrez_genes,
  organism      = "hsa",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

ekegg <- setReadable(
  ekegg,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)



#Result table
kegg_df <- as.data.frame(ekegg)
str(kegg_df)
head(kegg_df[, c("ID","Description","GeneRatio","p.adjust","geneID")])




#KEGG Barplot (Red + Maroon theme)
library(ggplot2)

# Prepare plotting data
plot_df <- kegg_df
plot_df$Description <- factor(
  plot_df$Description,
  levels = plot_df$Description
)

# Barplot
p_kegg <- ggplot(plot_df, aes(
  x = Count,
  y = Description
)) +
  geom_bar(
    stat = "identity",
    fill = "#8B0000",   # maroon
    width = 0.6
  ) +
  geom_text(
    aes(label = paste0("Adj.P = ", signif(p.adjust, 3))),
    hjust = -0.1,
    size = 4
  ) +
  labs(
    x = "Gene Count",
    y = NULL,
    title = "KEGG Pathway Enrichment Analysis"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

p_kegg



#advance

library(ggplot2)

ggplot(kegg_df,
       aes(x = Count,
           y = Description)) +
  geom_bar(stat = "identity",
           fill = "#8B0000",   # dark red
           width = 0.6) +
  geom_text(
    aes(label = paste0("Adj.P = ",
                       formatC(p.adjust, format = "e", digits = 2))),
    hjust = -0.1,
    size = 4
  )

  labs(
    title = "KEGG Pathway Enrichment Analysis",
    x = "Gene Count",
    y = NULL,
    caption = "Enriched genes: GSTM5, NAT1, TK1"
  ) +
  theme_classic(base_size = 12) +
  xlim(0, 4)


#eps 
ggsave(
  "Figure_KEGG_Barplot.eps",
  plot = p_kegg,
  width = 7,
  height = 3,
  dpi = 1200,
  device = "eps"
)











                                          #GSEA Analysis


str(final_df)
str(risk_ord)


#expression matrix prepare (Expression matrix + group align)

library(data.table)

# expr_wide → matrix
expr_mat <- as.data.frame(expr_wide)
rownames(expr_mat) <- expr_mat$Gene
expr_mat$Gene <- NULL

# sample order match
samples_use <- intersect(colnames(expr_mat), final_df$Sample_ID)
expr_mat <- expr_mat[, samples_use]

# group factor
group <- factor(
  final_df$risk_group[match(samples_use, final_df$Sample_ID)],
  levels = c("Low-risk", "High-risk")
)

table(group)





#Ranking metric (limma ranking)
library(limma)

design <- model.matrix(~ group)
fit <- lmFit(expr_mat, design)
fit <- eBayes(fit)

deg_all <- topTable(
  fit,
  coef = "groupHigh-risk",
  number = Inf,
  sort.by = "none"
)

# ranking vector
gene_ranking <- deg_all$logFC
names(gene_ranking) <- rownames(deg_all)

# remove NA & sort
gene_ranking <- gene_ranking[!is.na(gene_ranking)]
gene_ranking <- sort(gene_ranking, decreasing = TRUE)

head(gene_ranking)



gene_ranking <- read.csv("gene_ranking.csv")


#SYMBOL → ENTREZ ID mapping (for GSEA)
library(clusterProfiler)
library(org.Hs.eg.db)

gene_map_all <- bitr(
  names(gene_ranking),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

# merge with ranking
gene_ranking_df <- data.frame(
  SYMBOL = names(gene_ranking),
  logFC  = gene_ranking
)

gene_ranking_df <- merge(
  gene_ranking_df,
  gene_map_all,
  by = "SYMBOL"
)

# final ENTREZ ranking vector
gene_ranking_entrez <- gene_ranking_df$logFC
names(gene_ranking_entrez) <- gene_ranking_df$ENTREZID
gene_ranking_entrez <- sort(gene_ranking_entrez, decreasing = TRUE)

length(gene_ranking_entrez)



#KEGG-GSEA
gsea_kegg <- gseKEGG(
  geneList     = gene_ranking_entrez,
  organism     = "hsa",
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

class(gsea_kegg)




#plot
library(dplyr)

# GSEA result table
gsea_kegg@result <- gsea_kegg@result %>%
  mutate(
    pvalue_fmt  = ifelse(pvalue < 0.001, "0.00", sprintf("%.2f", pvalue)),
    padj_fmt    = ifelse(p.adjust < 0.001, "0.00", sprintf("%.2f", p.adjust))
  )


#proteasome pathway
pathway_name <- "Proteasome"

path_id <- gsea_kegg@result$ID[
  gsea_kegg@result$Description == pathway_name
]

path_id



#FINAL GSEA plot
library(enrichplot)
library(ggplot2)

p_gsea <- gseaplot2(
  gsea_kegg,
  geneSetID    = path_id,
  title        = "KEGG Proteasome Pathway",
  color        = "red",
  pvalue_table = TRUE,
  base_size    = 14
)



#save in eps

ggsave(
  filename = "Figure5C_GSEA_Proteasome.eps",
  plot     = p_gsea,
  device   = cairo_ps,   # 🔥 THIS is the key
  width    = 7,
  height   = 5,
  units    = "in",
  dpi      = 1000
)












                                               #step: 10 — Clinical Utility Model



# 10.1.1 Required package


install.packages("remotes")
remotes::install_github("yikeshu0611/ggDCA")



library(survival)
library(rms)
library(ggDCA)
library(ggplot2)





# 10.1.2: Clean & prepare data (NO risk)

step10_df <- final_df[, c(
  "OS_days",
  "OS_event",
  "RiskScore",
  "age",
  "tumor_grade"
)]

# Remove missing values
step10_df <- na.omit(step10_df)

# Convert grade to factor (IMPORTANT)
step10_df$tumor_grade <- factor(
  step10_df$tumor_grade,
  levels = c("G1", "G2", "G3"),
  ordered = TRUE
)







# STEP 10.2 — Multivariate Cox model (Nomogram base)


cox_nomogram <- coxph(
  Surv(OS_days, OS_event) ~ RiskScore + age + tumor_grade,
  data = step10_df
)

summary(cox_nomogram)







# STEP 10.3 — Nomogram construction (Figure 6A)

# Nomogram 



library(survival)
library(rms)
library(ggplot2)


# Clean Data
step10_df <- final_df[, c("OS_days", "OS_event", "RiskScore", "age", "tumor_grade")]
step10_df <- na.omit(step10_df)

# Re-level tumor grade
step10_df$tumor_grade <- factor(step10_df$tumor_grade, 
                                levels = c("G1", "G2", "G3"), 
                                ordered = FALSE)

# Environment Setup for RMS
dd <- datadist(step10_df)
options(datadist = "dd")

# Fit Cox Model (cph)
cph_model <- cph(Surv(OS_days, OS_event) ~ RiskScore + age + tumor_grade,
                 data = step10_df,
                 x = TRUE, y = TRUE, surv = TRUE)

# Create the survival calculator function from the model
surv_calc <- Survival(cph_model)

# Build Nomogram Object
nom <- nomogram(cph_model, 
                fun = list(function(x) surv_calc(365, x),   # 1 Year
                           function(x) surv_calc(1095, x),  # 3 Years
                           function(x) surv_calc(1825, x)), # 5 Years
                funlabel = c("1-Year OS", "3-Year OS", "5-Year OS"),
                lp = FALSE, 
                fun.at = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.1))

#Save
setEPS()
postscript("Figure6A_Nomogram.eps", width = 10, height = 7)

plot(nom, 
     xfrac = .25, 
     cex.axis = 0.85, 
     cex.var = 1, 
     col.grid = gray(c(0.8, 0.95)))

dev.off()

message("Nomogram created successfully!")


















# Calibration Curve


# PRE-CHECK



library(rms)

options(contrasts = c("contr.treatment", "contr.treatment"))
dd <- datadist(step10_df)
options(datadist = "dd")

cph_model <- cph(
  Surv(OS_days, OS_event) ~ RiskScore + age + tumor_grade,
  data = step10_df,
  x = TRUE,
  y = TRUE,
  surv = TRUE,
  time.inc = 365   # ⭐ THIS IS THE KEY FIX
)


# 1-year calibration (365 days)
cal_1y <- calibrate(
  cph_model,
  method = "boot",
  u = 365,
  m = 80,
  B = 300
)

# 3-year calibration
cal_3y <- calibrate(
  cph_model,
  method = "boot",
  u = 1095,
  m = 60,
  B = 300
)

# 5-year calibration (apparent – no resampling)
cal_5y <- calibrate(
  cph_model,
  method = "boot",
  u = 1825,
  B = 0
)




# save


postscript(
  "Figure6B_Calibration_RedBlue.eps",
  width = 7.5,
  height = 7.5,
  horizontal = FALSE,
  onefile = FALSE,
  paper = "special"
)

par(
  mar = c(5,5,3,2),
  bg = "white"
)

# Empty base plot using cal_1y
plot(
  cal_1y,
  lwd = 2.8,
  col = "#E31A1C",
  lty = 1,
  pch = 16,
  xlim = c(0,1),
  ylim = c(0,1),
  subtitles = FALSE,
  xlab = "Predicted Overall Survival Probability",
  ylab = "Observed Overall Survival Probability",
  cex.lab = 1.2,
  cex.axis = 1.1
)

# Add 3-year calibration curve (BLUE)
plot(
  cal_3y,
  add = TRUE,
  lwd = 2.8,
  col = "#1F78B4",
  lty = 1,
  pch = 17
)

# Ideal diagonal reference
abline(0, 1, lty = 2, lwd = 2, col = "gray60")

legend(
  "topleft",
  legend = c(
    "1-year OS (365 days)",
    "3-year OS (1095 days)",
    "Ideal reference"
  ),
  col = c("#E31A1C", "#1F78B4", "gray60"),
  lwd = c(2.8, 2.8, 2),
  pch = c(16, 17, NA),
  lty = c(1, 1, 2),
  bty = "n",
  cex = 1.05
)

dev.off()


















library(survival)
library(rms)
library(ggDCA)
library(ggplot2)


View(step10_df)
# Decision curve analysis  (DCA plot)


install.packages("rmda")

library(rmda)
library(survival)






# ===============================
# Time-specific OS events
# ===============================
step10_df$OS_event_1y <- ifelse(
  step10_df$OS_days <= 365 & step10_df$OS_event == 1, 1, 0
)

step10_df$OS_event_3y <- ifelse(
  step10_df$OS_days <= 1095 & step10_df$OS_event == 1, 1, 0
)

step10_df$OS_event_5y <- ifelse(
  step10_df$OS_days <= 1825 & step10_df$OS_event == 1, 1, 0
)




library(rmda)


library(rmda)
library(survival)

step10_df$OS_1y <- ifelse(step10_df$OS_days <= 365 & step10_df$OS_event == 1, 1, 0)
step10_df$OS_3y <- ifelse(step10_df$OS_days <= 1095 & step10_df$OS_event == 1, 1, 0)
step10_df$OS_5y <- ifelse(step10_df$OS_days <= 1825 & step10_df$OS_event == 1, 1, 0)



library(rmda)

dca_1y <- decision_curve(
  OS_1y ~ RiskScore + age + tumor_grade,
  data = step10_df,
  family = binomial(link = "logit"),
  thresholds = seq(0, 0.5, by = 0.01),
  study.design = "cohort"
)

dca_3y <- decision_curve(
  OS_3y ~ RiskScore + age + tumor_grade,
  data = step10_df,
  family = binomial(link = "logit"),
  thresholds = seq(0, 0.5, by = 0.01),
  study.design = "cohort"
)

dca_5y <- decision_curve(
  OS_5y ~ RiskScore + age + tumor_grade,
  data = step10_df,
  family = binomial(link = "logit"),
  thresholds = seq(0, 0.5, by = 0.01),
  study.design = "cohort"
)



# Save figure 

postscript(
  "Figure6C_DCA_1_3_5_years.eps",
  width = 7.8,
  height = 7.2,
  horizontal = FALSE,
  onefile = FALSE,
  paper = "special"
)

par(
  mar = c(5,5,2,2),
  bg = "white"
)

plot_decision_curve(
  list(dca_1y, dca_3y, dca_5y),
  curve.names = c("1-year OS", "3-year OS", "5-year OS"),
  col = c("#1F78B4", "#E31A1C", "#33A02C"),
  lwd = 2.8,
  confidence.intervals = FALSE,
  xlab = "Threshold Probability",
  ylab = "Net Benefit",
  legend.position = "bottomright"
)

dev.off()

















