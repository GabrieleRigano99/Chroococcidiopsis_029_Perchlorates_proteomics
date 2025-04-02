library("limma")
library("ggplot2")
#I upload the already imputated data transformed in log2
proteomics_data <- read.table("limma_input_029.tsv", header = TRUE, sep = "\t")
head(proteomics_data)
#I set the conditions and I put an x before numbers for practicality
conditions <- factor(c(rep("c", 3), rep("x2", 3), rep("x12", 3), rep("x24", 3)))
design <- model.matrix(~ 0 + conditions)
colnames(design) <- levels(conditions)
input <- proteomics_data[,-14]
head(input)
#I use the lmFit function of my data into the design
fit <- lmFit(input, design)
#I build contrasts
contrast_matrix <- makeContrasts(
  x2_vs_c = x2 - c,
  x12_vs_c = x12 - c,
  x24_vs_c = x24 - c,
  x12_vs_x2 = x12 - x2,
  x24_vs_x12 = x24 - x12,
  x24_vs_x2 = x24 - x2,
  levels = design
  
)
#I do the analysis using limma-trend
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE)
#I get the results for every condition
results_2_vs_Control <- topTable(fit2, coef = "x2_vs_c", adjust.method = "BH", number = Inf,genelist = proteomics_data$Protein_Group)
results_12_vs_Control <- topTable(fit2, coef = "x12_vs_c", adjust.method = "BH", number = Inf,genelist = proteomics_data$Protein_Group)
results_24_vs_Control <- topTable(fit2, coef = "x24_vs_c", adjust.method = "BH", number = Inf,genelist = proteomics_data$Protein_Group)
results_12_vs_2 <- topTable(fit2, coef = "x12_vs_x2", adjust.method = "BH", number = Inf,genelist = proteomics_data$Protein_Group)
results_24_vs_12 <- topTable(fit2, coef = "x24_vs_x12", adjust.method = "BH", number = Inf,genelist = proteomics_data$Protein_Group)
results_24_vs_2 <- topTable(fit2, coef = "x24_vs_x2", adjust.method = "BH", number = Inf,genelist = proteomics_data$Protein_Group)

write.csv(results_2_vs_Control, file="results_2_vs_Control.csv")
write.csv(results_12_vs_Control, file="results_12_vs_Control.csv")
write.csv(results_24_vs_Control, file="results_24_vs_Control.csv")
write.csv(results_24_vs_12, file="results_24_vs_12.csv")
write.csv(results_24_vs_2, file="results_24_vs_2.csv")
write.csv(results_12_vs_2, file="results_12_vs_2.csv")