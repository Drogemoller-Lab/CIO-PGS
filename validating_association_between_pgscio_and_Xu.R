
###Logistic regression analyses
library(ggplot2)
library(devtools)
install.packages("devtools")
library(geneticriskR)
library(webr)
library(dplyr)
library(tidyr)
library(multipleROC)
library(data.table)
library(fmsb)

###Load relevant files
phenotype_file <- fread("/home/projects/cisplatin/Xu_CIO_study_data/phenotype/Ototoxicity_phenotype_20130909_v2.txt", head =TRUE)
Xu_PCA <- fread("/home/projects/archive/cisplatin/PGS/results_v4/Xu_LDpruned_pca.eigenvec", header=FALSE)
Xu_kalra_deseq_milo_adjusted_PRS <- fread("/home/projects/archive/cisplatin/PGS/results_v4/filtered_ambiguous_snps/Xu_kalra_deseq_milo_adjusted_PRS.profile")
kalra_PRS <- fread("/home/projects/archive/cisplatin/PGS/results_v4/filtered_ambiguous_snps/no_ambiguous_SNPS_kalra_hg38.raw")

###Tidy datasets
#PCA
Xu_PCA <- subset(Xu_PCA, select = -V2)
new_names <- c("FID","P1","P2","P3","P4","P5","P6","P7","P8","P9","P10")
colnames(Xu_PCA) <- new_names
Xu_PCA <- Xu_PCA %>% 
  mutate(FID = gsub("_Ototoxicity.*", "", FID))
Xu_PCA$FID[237] <- gsub("_PG_0024", "", Xu_PCA$FID[237])
Xu_PCA$FID[238] <- gsub("_TR_0009", "", Xu_PCA$FID[238])

#PRS (deseq_milo_adjusted)
Xu_kalra_deseq_milo_adjusted_PRS <- subset(Xu_kalra_deseq_milo_adjusted_PRS, select = -FID)
Xu_kalra_deseq_milo_adjusted_PRS <- Xu_kalra_deseq_milo_adjusted_PRS %>%
  mutate(IID = gsub("_Ototoxicity.*", "", IID))
Xu_kalra_deseq_milo_adjusted_PRS$IID[237] <- gsub("_PG_0024", "", Xu_kalra_deseq_milo_adjusted_PRS$IID[237])
Xu_kalra_deseq_milo_adjusted_PRS$IID[238] <- gsub("_TR_0009", "", Xu_kalra_deseq_milo_adjusted_PRS$IID[238])

###Creating the functions for logistic regression analyses
simple_logistic_reg_test <- function(df, phenotype, ...) {
  simple_logit_output <- glm(df[[phenotype]] ~ SCORE, data = df, family = "binomial") 
  summary(simple_logit_output)
  return(simple_logit_output)
}

complex_logistic_regwithoutscore <- function(df, PHENO) {
  complex_logit_output <- glm(df[["PHENO"]] ~ Dxage + Protocol + CSI, data = df, family = "binomial")
  summary(complex_logit_output)
  return(complex_logit_output)
}

complex_logistic_reg <- function(df, PHENO) {
  complex_logit_output <- glm(df[["PHENO"]] ~ SCORE * CSI + Dxage + Protocol, data = df, family = "binomial")
  summary(complex_logit_output)
  return(complex_logit_output)
}

###FOR DIFFERENTIALLY EXPRESSED GENES IN ANY CLUSTER, MILO (Padj<0.05)
#Merge
phenotype_pca <- merge(phenotype_file, Xu_PCA, by.x = "Sample_id", by.y="FID") 
phenotype_pca_prs_deseq_milo_adjusted <- merge(phenotype_pca, Xu_kalra_deseq_milo_adjusted_PRS, by.x = "Sample_id", by.y="IID")
phenotype_pca_prs_deseq_milo_adjusted$PHENO <- phenotype_pca_prs_deseq_milo_adjusted$Chang_1a
#phenotype_pca_prs_deseq_adjusted <- phenotype_pca_prs_deseq_adjusted[phenotype_pca_prs_deseq_adjusted$Chang_1a == phenotype_pca_prs_deseq_adjusted$Chang_2a, ]

#Test normality of data using the shapiro test
prs_deseq_milo_adjusted <- Xu_kalra_deseq_milo_adjusted_PRS$SCORE
shapiro_test_deseq_milo_adjusted <- shapiro.test(prs_deseq_milo_adjusted)
print(shapiro_test_deseq_milo_adjusted)

#Plot the distribution of the PRS
prs_deseq_milo_adjusted <- data.frame(PRS = Xu_kalra_deseq_milo_adjusted_PRS$SCORE)
histogram_kalra_deseq_milo_adjusted <- ggplot(data = prs_deseq_milo_adjusted, aes(x = PRS)) +
  geom_histogram(bins =70, fill = "skyblue", color = "black") +
  labs(title = "Distribution of PRS", x = "PRS", y = "Frequency")

print(histogram_kalra_deseq_milo_adjusted)

#ROC analysis
simple_logit_output <- simple_logistic_reg_test(phenotype_pca_prs_deseq_milo_adjusted, "PHENO") # only for the PRS
complex_logi_output <- complex_logistic_reg(phenotype_pca_prs_deseq_milo_adjusted, "PHENO") # all of them
complex_logi_output_without_score <- complex_logistic_regwithoutscore(phenotype_pca_prs_deseq_milo_adjusted, "PHENO")# without PRS

proc_analysis(phenotype_pca_prs_deseq_milo_adjusted,"PHENO", simple_logit_output)
proc_analysis(phenotype_pca_prs_deseq_milo_adjusted, "PHENO", complex_logi_output)
proc_analysis(phenotype_pca_prs_deseq_milo_adjusted, "PHENO", complex_logi_output_without_score)

## Multiple ROC
a_deseq_milo_adjusted=multipleROC(PHENO ~ SCORE, data=phenotype_pca_prs_deseq_milo_adjusted,plot=FALSE)
b_deseq_milo_adjusted=multipleROC(PHENO ~ Dxage + Protocol + CSI + P1 + P2 + P3 + P4 + P5, data=phenotype_pca_prs_deseq_milo_adjusted,plot=FALSE)
c_deseq_milo_adjusted=multipleROC(PHENO ~ Dxage + Protocol + SCORE * CSI + P1 + P2 + P3 + P4 + P5, data=phenotype_pca_prs_deseq_milo_adjusted,plot=FALSE)  
d_deseq_milo_adjusted=multipleROC(PHENO ~ Dxage, data=phenotype_pca_prs_deseq_milo_adjusted,plot=FALSE) 
e_deseq_milo_adjusted=multipleROC(PHENO ~ Protocol, data=phenotype_pca_prs_deseq_milo_adjusted,plot=FALSE) 
f_deseq_milo_adjusted=multipleROC(PHENO ~ CSI, data=phenotype_pca_prs_deseq_milo_adjusted,plot=FALSE) 

plot_ROC(list(a_deseq_milo_adjusted,b_deseq_milo_adjusted,c_deseq_milo_adjusted),show.eta=FALSE,show.sens=FALSE)

###For CS==1/0
phenotype_pca_prs_deseq_milo_adjusted_cs0 <- subset(phenotype_pca_prs_deseq_milo_adjusted, CSI == 0)
phenotype_pca_prs_deseq_milo_adjusted_cs1 <- subset(phenotype_pca_prs_deseq_milo_adjusted, CSI == 1)
# Calculate the null model (model without PRS) using a logistic regression (as only case-control status is available) 
null.model <- glm(PHENO ~ Protocol + Dxage + CSI + P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + P10, 
                  data=phenotype_pca_prs_deseq_milo_adjusted, family = binomial)

null.model.pseudoR2 <- NagelkerkeR2(null.model)

#extract only the R2 value from the list
null.r2 <- as.numeric(null.model.pseudoR2[2])

# Calculate the full model (model with PRS) using a logistic regression (as only case-control status is available) 
model <- glm(PHENO ~ SCORE*CSI + Protocol + Dxage + P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + P10, 
             data = phenotype_pca_prs_deseq_milo_adjusted, family = binomial)
summary(model)
# Calculate the pseudo R2 of the full model (because we are using logistic regression)
model.pseudoR2 <- NagelkerkeR2(model)

#extract only the R2 value from the list
model.r2 <- as.numeric(model.pseudoR2[2])
# R2 of PRS is simply calculated as the model R2 minus the null R2
prs.r2 <- model.r2-null.r2
# Obtain the coeffcient and p-value of association of PRS as follows
prs.coef <- summary(model)$coeff["SCORE",]
prs.beta <- as.numeric(prs.coef[1])
prs.se <- as.numeric(prs.coef[2])
prs.p <- as.numeric(prs.coef[4])
# We can then store the results
prs.result_deseq_milo_adjusted <- rbind(R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se)
View(prs.result_deseq_milo_adjusted)

###Specificity and Sensitivity
# Fit logistic regression model
model <- glm(PHENO ~ SCORE*CSI + Protocol + Dxage + P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9 + P10, 
             data = phenotype_pca_prs_deseq_milo_adjusted, family = binomial)

# Predict probabilities
predicted_probs <- predict(model, type = "response")

# Classify observations
predicted_class <- ifelse(predicted_probs > 0.5, 1, 0)

# Calculate sensitivity and specificity
conf_matrix <- table(Actual = phenotype_pca_prs_deseq_milo_adjusted$PHENO, Predicted = predicted_class)
View(conf_matrix)
sensitivity <- conf_matrix[2, 2] / sum(conf_matrix[2, ])
specificity <- conf_matrix[1, 1] / sum(conf_matrix[1, ])

# Print or store the results
print(paste("Sensitivity:", sensitivity))
print(paste("Specificity:", specificity))
