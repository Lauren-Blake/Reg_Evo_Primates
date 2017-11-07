# This script is run the integration_methylation_expression_humans_chimps_only R Markdown file on the cluster

# Load libraries/data

library("edgeR")
library("limma")
library("plyr")

cpm.voom.cyclic <- readRDS("./Methylation_integration/human_chimp_orth_cpm_voom_cyclic.rds")
exp_methyl <- read.table("./Methylation_integration/human_chimp_orth_exp_methyl_7725_hum.txt", header = T, stringsAsFactors = F)
samples <- read.table("./Methylation_integration/human_chimp_orth_new_sample_info.txt", header = T,  stringsAsFactors = F)

species <- as.data.frame(samples[,4])
tissue <- as.data.frame(samples[,5])
RIN <- as.data.frame(samples[,6])

# 1) Human v Chimp Heart at FDR 5%
# First, pick out the genes that are DE (you run the entire analysis only on genes that are DE)
chimp_human_heart <- c(1, 5, 9, 13, 20, 24, 28)
FDR_level <- 0.05

design <- model.matrix(~ as.factor(species[chimp_human_heart,]) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]

## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), 7)) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), 7)) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_human_chimp_heart_FDR_5perc.txt")







# 2) Human v Chimp Heart at FDR 10%
# First, pick out the genes that are DE (you run the entire analysis only on genes that are DE)
chimp_human_heart <- c(1, 5, 9, 13, 20, 24, 28)
FDR_level <- 0.1

design <- model.matrix(~ as.factor(species[chimp_human_heart,]) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]

## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), 7)) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), 7)) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_human_chimp_heart_FDR_10perc.txt")




# 3) Human v Chimp Kidney at FDR 5%
# First, pick out the genes that are DE (you run the entire analysis only on genes that are DE)
chimp_human_heart <- c(2, 6, 10, 14, 17, 21, 25, 29)
FDR_level <- 0.05

design <- model.matrix(~ as.factor(species[chimp_human_heart,]) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]

## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
#for (k in 1:2){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_human_chimp_kidney_FDR_5perc.txt")





# 4) Human v Chimp Kidney at FDR 10%
# First, pick out the genes that are DE (you run the entire analysis only on genes that are DE)
chimp_human_heart <- c(2, 6, 10, 14, 17, 21, 25, 29)
FDR_level <- 0.1

design <- model.matrix(~ as.factor(species[chimp_human_heart,]) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]

## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
  #for (k in 1:2){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_human_chimp_kidney_FDR_10perc.txt")





# 5) Human v Chimp Liver at FDR 5%
# First, pick out the genes that are DE (you run the entire analysis only on genes that are DE)
chimp_human_heart <- c(3, 7, 11, 15, 18, 22, 26, 30)
FDR_level <- 0.05

design <- model.matrix(~ as.factor(species[chimp_human_heart,]) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]

## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
#for (k in 1:10){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_human_chimp_liver_FDR_5perc.txt")




# 6) Human v Chimp Liver at FDR 10%
# First, pick out the genes that are DE (you run the entire analysis only on genes that are DE)
chimp_human_heart <- c(3, 7, 11, 15, 18, 22, 26, 30)
FDR_level <- 0.10

design <- model.matrix(~ as.factor(species[chimp_human_heart,]) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]

## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
#for (k in 1:10){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_human_chimp_liver_FDR_10perc.txt")




# 7) Human v Chimp Lung at FDR 5%
# First, pick out the genes that are DE (you run the entire analysis only on genes that are DE)
chimp_human_heart <- c(4, 8, 12, 16, 19, 23, 27, 31)
FDR_level <- 0.05

design <- model.matrix(~ as.factor(species[chimp_human_heart,]) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]

## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
  #for (k in 1:10){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_human_chimp_lung_FDR_5perc.txt")



# 8) Human v Chimp Lung at FDR 10%
# First, pick out the genes that are DE (you run the entire analysis only on genes that are DE)
chimp_human_heart <- c(4, 8, 12, 16, 19, 23, 27, 31)
FDR_level <- 0.1

design <- model.matrix(~ as.factor(species[chimp_human_heart,]) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]

## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
  #for (k in 1:10){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_human_chimp_lung_FDR_10perc.txt")






# 9) Human Heart versus Kidney at FDR 5%
# First, pick out the genes that are DE (you run the entire analysis only on genes that are DE)
chimp_human_heart <- c(17, 20, 21, 24, 25, 28, 29)
FDR_level <- 0.05
tissue <- tissue[chimp_human_heart,]
tissue_no_extra <- droplevels.factor(tissue, "liver")


design <- model.matrix(~ as.factor(tissue_no_extra) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]

## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
#for (k in 1:2){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_human_only_heart_kidney_FDR_5perc.txt")




# 10) Human Heart versus Liver at FDR 5%
# First, pick out the genes that are DE (you run the entire analysis only on genes that are DE)
chimp_human_heart <- c(18, 20, 22, 24, 26, 28, 30)
FDR_level <- 0.05
tissue <- as.data.frame(samples[,5])
tissue <- tissue[chimp_human_heart,]
tissue_no_extra <- droplevels.factor(tissue, "liver")


design <- model.matrix(~ as.factor(tissue_no_extra) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]
dim(expression_values_only)

## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
  #for (k in 1:2){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_human_only_heart_liver_FDR_5perc.txt")





# 11) Human heart versus lung at FDR 5%
chimp_human_heart <- c(19, 20, 23, 24, 27, 28, 31)
FDR_level <- 0.05
tissue <- as.data.frame(samples[,5])
tissue <- tissue[chimp_human_heart,]
tissue_no_extra <- droplevels.factor(tissue, "liver")


design <- model.matrix(~ as.factor(tissue_no_extra) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]
dim(expression_values_only)
## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
  #for (k in 1:2){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_human_only_heart_lung_FDR_5perc.txt")




# 12) Human Kidney versus Liver at FDR 5%
chimp_human_heart <- c(17, 18, 21, 22, 25, 26, 29, 30)
FDR_level <- 0.05
tissue <- tissue[chimp_human_heart,]
tissue_no_extra <- droplevels.factor(tissue, "lung")


design <- model.matrix(~ as.factor(tissue_no_extra) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]
dim(expression_values_only)
## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
  #for (k in 1:2){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_human_only_kidney_liver_FDR_5perc.txt")







# 13) Human Kidney versus Lung at FDR 5%
chimp_human_heart <- c(17, 19, 21, 23, 25, 27, 29, 31)
FDR_level <- 0.05
tissue <- tissue[chimp_human_heart,]
tissue_no_extra <- droplevels.factor(tissue, "liver")


design <- model.matrix(~ as.factor(tissue_no_extra) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]
dim(expression_values_only)
## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
  #for (k in 1:2){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_human_only_kidney_lung_FDR_5perc.txt")







# 14) Human Liver versus Lung at FDR 5%
chimp_human_heart <- c(18, 19, 22, 23, 26, 27, 30, 31)
FDR_level <- 0.05
tissue <- as.data.frame(samples[,5])
tissue <- tissue[chimp_human_heart,]
tissue_no_extra <- droplevels.factor(tissue, "lung")


design <- model.matrix(~ as.factor(tissue_no_extra) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]

## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
  #for (k in 1:2){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_human_only_liver_lung_FDR_5perc.txt")







# Chimp Heart versus Kidney at FDR 5%
# First, pick out the genes that are DE (you run the entire analysis only on genes that are DE)
chimp_human_heart <- c(1, 2, 5, 6, 9, 10, 13, 14)
FDR_level <- 0.05
tissue <- as.data.frame(samples[,5])
tissue <- tissue[chimp_human_heart,]
tissue_no_extra <- droplevels.factor(tissue, "liver")


design <- model.matrix(~ as.factor(tissue_no_extra) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]

## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
  #for (k in 1:2){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_chimp_only_heart_kidney_FDR_5perc.txt")





# 15) Chimp Heart versus Liver at FDR 5%
# First, pick out the genes that are DE (you run the entire analysis only on genes that are DE)
chimp_human_heart <- c(1, 3, 5, 7, 9, 11, 13, 15)
FDR_level <- 0.05
tissue <- as.data.frame(samples[,5])
tissue <- tissue[chimp_human_heart,]
tissue_no_extra <- droplevels.factor(tissue, "lung")


design <- model.matrix(~ as.factor(tissue_no_extra) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]

## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
  #for (k in 1:2){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_chimp_only_heart_liver_FDR_5perc.txt")


# 17) Chimp Heart versus Lung at FDR 5%
# First, pick out the genes that are DE (you run the entire analysis only on genes that are DE)
chimp_human_heart <- c(1, 4, 5, 8, 9, 12, 13, 16)
FDR_level <- 0.05
tissue <- as.data.frame(samples[,5])
tissue <- tissue[chimp_human_heart,]
tissue_no_extra <- droplevels.factor(tissue, "liver")


design <- model.matrix(~ as.factor(tissue_no_extra) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]

## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
  #for (k in 1:2){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_chimp_only_heart_lung_FDR_5perc.txt")






# 17) Chimp Heart versus Kidney at FDR 5%
# First, pick out the genes that are DE (you run the entire analysis only on genes that are DE)
chimp_human_heart <- c(1, 2, 5, 6, 9, 10, 13, 14)
FDR_level <- 0.05
tissue <- as.data.frame(samples[,5])
tissue <- tissue[chimp_human_heart,]
tissue_no_extra <- droplevels.factor(tissue, "liver")


design <- model.matrix(~ as.factor(tissue_no_extra) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]

## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
  #for (k in 1:2){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_chimp_only_heart_kidney_FDR_5perc.txt")




# 18) Chimp kidney versus liver at FDR 5%
# First, pick out the genes that are DE (you run the entire analysis only on genes that are DE)
chimp_human_heart <- c(2, 3, 6, 7, 10, 11, 14, 15)
FDR_level <- 0.05
tissue <- as.data.frame(samples[,5])
tissue <- tissue[chimp_human_heart,]
tissue_no_extra <- droplevels.factor(tissue, "liver")


design <- model.matrix(~ as.factor(tissue_no_extra) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]

## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
  #for (k in 1:2){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_chimp_only_kidney_liver_FDR_5perc.txt")





# 19) Chimp kidney versus lung at FDR 5%
# First, pick out the genes that are DE (you run the entire analysis only on genes that are DE)
chimp_human_heart <- c(2, 4, 6, 8, 10, 12, 14, 16)
FDR_level <- 0.05
tissue <- as.data.frame(samples[,5])
tissue <- tissue[chimp_human_heart,]
tissue_no_extra <- droplevels.factor(tissue, "liver")


design <- model.matrix(~ as.factor(tissue_no_extra) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]

## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
  #for (k in 1:2){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_chimp_only_kidney_lung_FDR_5perc.txt")



# 19) Chimp liver versus lung at FDR 5%
# First, pick out the genes that are DE (you run the entire analysis only on genes that are DE)
chimp_human_heart <- c(3, 4, 7, 8, 11, 12, 15, 16)
FDR_level <- 0.05
tissue <- as.data.frame(samples[,5])
tissue <- tissue[chimp_human_heart,]
tissue_no_extra <- droplevels.factor(tissue, "liver")


design <- model.matrix(~ as.factor(tissue_no_extra) + RIN[chimp_human_heart,])
fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
fit_all <- eBayes(fit_all)

HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FDR_level), ]
human_chimp_heart <-  rownames(exp_methyl) %in% HvC_Heart_fit_all_5perc$genes
human_chimp_heart <- as.data.frame(human_chimp_heart)
counts_genes_in <- cbind(exp_methyl, human_chimp_heart)
counts_genes_in_cutoff <- subset(counts_genes_in, human_chimp_heart == "TRUE")
counts_genes_in_cutoff <- counts_genes_in_cutoff[,1:79]  
expression_values_only <- counts_genes_in_cutoff[,2:32]
methylation_values_only <- counts_genes_in_cutoff[,49:79]

## Obtain corrected data (regression for scenario 2) by regressing out RIN on a gene-by-gene basis
resid_methyl <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
expression_values <- as.data.frame(expression_values_only)
methylation_values <- as.data.frame(methylation_values_only)
for (i in 1:nrow(expression_values_only)){
  resid_methyl[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ t(methylation_values_only[i,chimp_human_heart]))$resid
}
rownames(resid_methyl) <- rownames(expression_values_only)

# Scenario 1

fit1 <- lmFit(expression_values_only[, chimp_human_heart], design)
fit1 <- eBayes(fit1)
HvC_Heart_fit1 = topTable(fit1, coef=2, adjust="BH", number=Inf, sort.by="none")

# Scenario 2

fit2 <- lmFit(resid_methyl, design)
fit2 <- eBayes(fit2)
HvC_Heart_fit2 = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")

# Bind DE from fit1 and DE from fit2 together

HvC_Heart_fits12 <- as.data.frame(cbind(rownames(HvC_Heart_fit1), HvC_Heart_fit1$adj.P.Val, HvC_Heart_fit2$adj.P.Val), stringsAsFactors = FALSE)
HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")

DE_both <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_both)

DE_before <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
nrow(DE_before)

DE_after <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
nrow(DE_after)

### Permute the methylation values, then re-run

run1000 <- array(0, dim = c(3, 1000)) 

for (k in 1:1000){
  #for (k in 1:2){
  # Run 1 for all genes
  run1 <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  for (i in 1:nrow(expression_values_only)){
    df <- as.data.frame(methylation_values_only[i,chimp_human_heart])
    new_df <- as.data.frame(sample(df))
    run1[i,] <- t(new_df)
  }
  
  
  
  # Make an array to hold the residuals
  
  resid_methyl_perm <- array(0, dim = c(nrow(expression_values_only), length(chimp_human_heart))) 
  
  
  #for (i in 1:2){
  for (i in 1:nrow(expression_values_only)){
    resid_methyl_perm[i,] <- lm(t(expression_values_only[i,chimp_human_heart]) ~ run1[i,])$resid
  }
  
  rownames(resid_methyl_perm) <- rownames(expression_values_only)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl_perm, design)
  fit2 <- eBayes(fit2)
  HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
  
  
  HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR, HvC_Heart_perm$adj.P.Val)
  
  DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
  DE_total_report <- nrow(DE_total)
  
  DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
  DE_before_report <- nrow(DE_before)
  
  DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
  DE_after_report <- nrow(DE_after)
  
  new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
  
  run1000[,k] <- new_numbers
}

median(run1000[1,1:1000])
median(run1000[2,1:1000])
median(run1000[3,1:1000])

write.table(run1000, "./perm_chimp_only_liver_lung_FDR_5perc.txt")