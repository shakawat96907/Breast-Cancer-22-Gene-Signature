
# Load necessary library
library(survival)

# Simplify Stage for better analysis
final_merged$Stage_Group <- ifelse(grepl("Stage I$|Stage IA|Stage IB|Stage II$|Stage IIA|Stage IIB", 
                                         final_merged$ajcc_pathologic_tumor_stage), "Early", "Late")
# Set Reference levels
final_merged$Stage_Group <- factor(final_merged$Stage_Group, levels = c("Early", "Late"))





# Univariate Cox Regression

# Univariate for Risk Score
uni_risk <- coxph(Surv(time, status) ~ Risk_Score, data = final_merged)
summary(uni_risk)

# Univariate for Age
uni_age <- coxph(Surv(time, status) ~ age_at_diagnosis, data = final_merged)
summary(uni_age)

# Univariate for Stage
uni_stage <- coxph(Surv(time, status) ~ Stage_Group, data = final_merged)
summary(uni_stage)









# Multivariate Cox Regression


# Multivariate Analysis
multi_cox <- coxph(Surv(time, status) ~ Risk_Score + age_at_diagnosis + Stage_Group, data = final_merged)

# Display result
summary(multi_cox)

# Visualizing with a Forest Plot (for your paper)
library(survminer)
ggforest(multi_cox, data = final_merged, main = "Multivariate Hazard Ratio Plot")











# Forest Plot

# 1. Scaling the Risk Score
# This standardizes the score so 1 unit increase is more meaningful
final_merged$Risk_Score_Scaled <- scale(final_merged$Risk_Score)

# 2. Run the Multivariate Model again with the Scaled Score
multi_cox_scaled <- coxph(Surv(time, status) ~ Risk_Score_Scaled + age_at_diagnosis + Stage_Group, 
                          data = final_merged)

# 3. Check the summary to see the new HR (exp(coef))
summary(multi_cox_scaled)

# 4. Save the Scaled Forest Plot as 1000 DPI TIFF
library(survminer)
tiff("Forest_Plot_Scaled_1000DPI.tiff", 
     width = 9, height = 5, units = 'in', res = 1000, compression = "lzw")

ggforest(multi_cox_scaled, 
         data = final_merged, 
         main = "Hazard Ratio (Independent Validation)",
         cpositions = c(0.02, 0.22, 0.4), 
         fontsize = 0.8, 
         refLabel = "Reference",
         noDigits = 3)

dev.off()
