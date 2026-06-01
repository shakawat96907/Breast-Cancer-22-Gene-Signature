# Integrative Transcriptomic Analysis Identifies a Novel 22-Gene Signature Driving Breast Cancer Progression via the Proteasome-Chemokine Axis

🔬 **Project Overview**
This repository hosts the computational architecture, reproducible R pipelines, and molecular validation workflows for an integrative transcriptomic profiling study. By screening a discovery cohort of **1,185 patients** and externally validating the model across **1,094 patients**, this study characterizes a novel **"Proliferation-Immune Axis"**—dictated by the antagonistic interplay between proteasome-driven tumor hyper-proliferation and chemokine-mediated immune surveillance.

---

### 📂 Repository Structure & Module Mapping

The core source code and computational pipelines are systematically organized into the following operational directories, matching the chronological analytical execution of the manuscript:

* 📁 **`Main_Internal_Analysis/`**: Holds the initial discovery phase code. Implements background correction, normalization, and genome-wide Univariate Cox regression screening along with Differential Expression Analysis (DEA).
* 📁 **`Cox_Analysis/`**: Contains the scripts for dimensionality reduction and model construction using **LASSO Cox Regression** to extract the final 22-gene signature and calculate individualized patient Risk Scores.
* 📁 **`KM_Curve/`**: Code for Kaplan-Meier survival curves and log-rank test executions to validate risk stratification between high- and low-risk cohorts.
* 📁 **`ROC_Curve/`**: Features time-dependent ROC curve analysis, calculating Area Under the Curve (AUC) values for 1-, 3-, and 5-year overall survival milestones.
* 📁 **`External_Validation/`**: Automates external validation pipelines. Fetches, normalizes, and processes the independent validation cohort (**TCGA-BRCA**).
* 📁 **`C-Index/`**: Code calculating Multivariate Cox proportional hazard models and computing the **Concordance Index (C-index)** to measure the stable discriminative capability of the model.
* 📁 **`Nomogram_Calibration/`**: Features the clinical translation architecture. Utilizes the genetic Risk Score alongside clinical parameters to construct predictive Nomograms, multi-year Calibration Curves, and execute **Decision Curve Analysis (DCA)**.

---

### 📊 Full Methodological Workflow (Point-by-Point)

#### 2.1 Data Acquisition and Pre-processing
* **Discovery Cohort:** Analyzed a comprehensive whole-transcriptome profile and matching clinical annotations of **1,185 breast cancer patients** sourced from the **Gene Expression Omnibus (GEO)** database under accession number **GSE96058** (originating from the **SCAN-B initiative**).
* **Validation Cohort:** Obtained transcriptome and clinical data for the **TCGA-BRCA cohort ($N=1,094$)** from the GDC Data Portal using **`TCGAbiolinks` (v2.28, 2023)** and **`curatedTCGAData`** R packages for independent external validation.
* **Environment:** All computational data processing and statistical analyses were performed using **R (v4.5.2)**.
* **Normalization & Transformation:** Raw count data were background-corrected and normalized utilizing the **`edgeR`** and **`limma`** packages. Gene expression levels were transformed using the **$\log_2(x+1)$** function to stabilize variance.

#### 2.2 Screening for Survival-Associated Candidate Genes
* **Dual-Screening Approach:** Employed a rigorous dual-engine filter to identify highly reliable prognostic indicators.
* **Univariate Cox Filter:** Conducted a genome-wide univariate Cox regression to identify genes substantially linked with overall survival (OS) at a significance threshold of **$p < 0.05$**.
* **Differential Expression Analysis (DEA):** Utilized the **`DESeq2` package (v1.40, 2023)** to identify Differentially Expressed Genes (DEGs) using strict cutoffs of absolute **$|\log_2\text{FC}| > 1$** and Benjamini-Hochberg adjusted **$p < 0.05$**.
* **Intersection:** The final candidate gene set was derived from the strict intersection of survival-associated genes and identified DEGs.

#### 2.3 Construction of the Prognostic Signature via LASSO
* **Multicollinerity Mitigation:** Applied the **Least Absolute Shrinkage and Selection Operator (LASSO) Cox regression** algorithm via the **`glmnet`** R package to eliminate multi-collinearity and build a simplified model.
* **Cross-Validation:** Identified the optimal penalty parameter ($\lambda$) via **10-fold cross-validation**, specifically selecting the $\lambda$ value within one standard error of the minimum value (**$\lambda_{1\text{se}}$**) to maximize model stability.
* **Risk Score Formulation:** Individual patient risk scores were generated using the formula:
  $$\text{Risk Score} = \sum_{i=1}^{n} \beta_i \times \text{Expr}_i$$
  *(Where $\beta_i$ represents the LASSO regression coefficient, and $\text{Expr}_i$ is the $\log_2$-transformed expression value of gene $i$)*.

#### 2.4 Prognostic Validation and Risk Stratification
* **Risk Stratification:** Patients in both discovery (SCAN-B) and validation (TCGA) cohorts were stratified into **high-risk and low-risk groups** based on the **median risk score** threshold.
* **Survival Disparity Evaluation:** Evaluated survival differences using **Kaplan-Meier curves** and the **log-rank test** via the **`survival`** package.
* **Time-Dependent ROC Analysis:** Measured model sensitivity and specificity for predicting **1-, 3-, and 5-year survival** milestones using the **`timeROC`** package.
* **Genomic Landscape Visualization:** Generated comprehensive heatmaps using the **`pheatmap`** package to display gene expression profiles, patient survival status, and clinical features aligned with the ascending risk scores.

#### 2.5 Independent Prognostic Validation and Multivariate Analysis
* **Prognostic Independence:** Built univariate and multivariate Cox proportional hazards models to guarantee that the signature serves as an independent prognostic factor.
* **Clinical Covariates Control:** Controlled for baseline clinicopathological variables: **age** and **pathologic tumor stage** (simplified into **"Early"**: Stages I–II vs. **"Late"**: Stages III–IV).
* **Risk Score Standardization:** The risk score was **standardized (scaled)** to make the interpretation of the Hazard Ratio (HR) straightforward and biologically meaningful.
* **Discriminative Capacity:** Computed the **Concordance Index (C-index)** to measure and compare the stable discriminative performance of the predictive models.

#### 2.6 Functional Enrichment and PPI Network Analysis
* **Functional Elucidation:** Enriched biological processes and molecular mechanisms were mapped using the **`clusterProfiler`** package via **Gene Ontology (GO)** terms and **Kyoto Encyclopedia of Genes and Genomes (KEGG)** pathways.
* **Gene Set Enrichment Analysis (GSEA):** Conducted GSEA to evaluate hallmark pathways operating within the risk groups.
* **PPI Network Modeling:** Constructed a Protein-Protein Interaction (PPI) network using the **STRING database v12.0**, integrating evidence from experiments, databases, and co-expression.
* **Network Tuning:** Applied a minimum required interaction score of **0.150 (low confidence)** to optimize the capture of indirect regulatory links within the *Proliferation-Immune Axis*.
* **Hub Gene Identification:** Used the network’s built-in topological algorithms to pinpoint core **hub genes (TK1, CXCL13)** based on **Node Degree** ranking.

#### 2.7 Protein-level Validation via Immunohistochemistry (IHC)
* **Database Verification:** Sourced clinical protein expression validation data from the **Human Protein Atlas (HPA v23.0, 2023)**.
* **Target Antibodies:** Analyzed representative tissue microscopic images from breast invasive ductal carcinoma patients stained with specific antibodies: **CAB004683 (for TK1)** and **HPA052613 (for CXCL13)**.
* **Validation Parameters:** Validated transcriptomic results by assessing protein **staining intensity** (Strong, Moderate, Weak), **staining quantity** (percentage of stained tumor cells), and exact **cellular localization**.

#### 2.8 Nomogram Development and Clinical Utility
* **Clinical Translation Tool:** Constructed a composite predictive **Nomogram** using the **`rms`** package in R, integrating the genetic Risk Score with critical clinical features (**age, stage, and grade**).
* **Model Calibration:** Generated **Calibration Curves** to evaluate how closely the nomogram’s predicted survival probabilities matched the actual observed survival outcomes.
* **Decision Curve Analysis (DCA):** Performed DCA to quantify the **net clinical benefit** of the nomogram across a wide range of threshold probabilities, directly comparing its utility against traditional "treat-all" or "treat-none" clinical approaches.

---

### 🛠️ Computational Prerequisites
* **Environment:** R (v4.5.2) 
* **Bioconductor & CRAN Libraries:** `TCGAbiolinks`, `curatedTCGAData`, `edgeR`, `limma`, `DESeq2`, `glmnet`, `survival`, `timeROC`, `pheatmap`, `clusterProfiler`, `rms`

---
### 📜 Declarations & Contact
* **Author:** Md. Shakawat Hossain  
* **Affiliation:** Department of Biochemistry and Molecular Biology, Shahjalal University of Science and Technology (SUST), Sylhet-3114, Bangladesh.  
* **Correspondence:** [shakawathossain96907@gmail.com](mailto:shakawathossain96907@gmail.com) | [ORCID Profile](https://orcid.org/0009-0008-5050-2275)
* **Data Citations:** Sourced via NCBI GEO (GSE96058) and GDC TCGA Data Portals. Full reproducible scripts are completely open-source across this repository.
