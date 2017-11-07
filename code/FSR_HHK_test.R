#################### SCRIPT FOR FSR TEST ##################

#!/usr/bin/env Rscript

library("edgeR")
library("limma")
library("plyr")
library("ashr")

FSR_level <- 0.05
cpm.voom.cyclic <- readRDS("/mnt/gluster/home/leblake/Methylation/Integration_humans_chimps/Methylation_integration/human_chimp_orth_cpm_voom_cyclic.rds")
exp_methyl <- read.table("/mnt/gluster/home/leblake/Methylation/Integration_humans_chimps/Methylation_integration/human_chimp_orth_exp_methyl_7725_hum.txt", header = T, stringsAsFactors = F)
samples <- read.table("/mnt/gluster/home/leblake/Methylation/Integration_humans_chimps/Methylation_integration/human_chimp_orth_new_sample_info.txt", header = T,  stringsAsFactors = F)

species <- as.data.frame(samples[,4])
tissue <- as.data.frame(samples[,5])
RIN <- as.data.frame(samples[,6])

integration_ash <- function(chimp_human_heart, FSR_level){
  FDR_level <- FSR_level
  tissue <- as.data.frame(samples[,5])
  tissue <- tissue[chimp_human_heart,]
  tissue_no_extra <- droplevels.factor(tissue)
  
  
  design <- model.matrix(~ as.factor(tissue_no_extra) + RIN[chimp_human_heart,])
  fit_all <- lmFit(cpm.voom.cyclic[,chimp_human_heart], design)
  fit_all <- eBayes(fit_all)
  
  
  # Run ASH
  
  # Prepare the data for ASH
  tests <- colnames(fit_all$coefficients)
  results <- vector(length = length(tests), mode = "list")
  names(results) <- tests
  
  # Perform multiple testing correction with adaptive shrinkage (ASH) 
  #
  # x - object MArrayLM from eBayes output
  # coef - coefficient tested by eBayes
  
  run_ash <- function(x, coef){
    #stopifnot(class(x) == "MArrayLM", coef %in% colnames(x$coefficients),
    #             length(unique(x$df.total) == 1))
    result <- ash(betahat = x$coefficients[, coef], sebetahat = x$stdev.unscaled[, coef] * sqrt(x$s2.post), df = x$df.total[1])
    return(result)
  }
  
  get_results <- function(x, number = nrow(x$coefficients), sort.by = "none",
                          ...) {
    # x - object MArrayLM from eBayes output
    # ... - additional arguments passed to topTable
    stopifnot(class(x) == "MArrayLM")
    results <- topTable(x, number = number, sort.by = sort.by, ...)
    return(results)
  }
  
  # Get lfsr, lfdr, s value, q value, and a beta_est value. 
  # Extract limma results
  results <- get_results(fit_all, coef = tests[2])
  # Add mutliple testing correction with ASH
  output_ash <- run_ash(fit_all, coef = tests[2])
  results <- cbind(results, sebetahat = output_ash$data$s, lfsr = output_ash$result$lfsr,
                   lfdr = output_ash$result$lfdr, qvalue = output_ash$result$qvalue,
                   svalue = output_ash$result$svalue, beta_est = output_ash$result$PosteriorMean, se_est =
                     output_ash$result$PosteriorSD)
  
  # Find the genes < FSR
  #HvC_Heart_fit_all = topTable(fit_all, coef=2, adjust="BH", number=Inf, sort.by="none")
  #HvC_Heart_fit_all_5perc <- HvC_Heart_fit_all[which(HvC_Heart_fit_all$adj.P.Val < FSR_level), ]
  #HvC_Heart_fit_all_5perc <- results[which(results$qvalue < FSR_level), ]
  HvC_Heart_fit_all_5perc <- results[which(results$svalue < FSR_level), ]
  
  
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
  
  results <- get_results(fit1, coef = tests[2])
  # Add mutliple testing correction with ASH
  output_ash <- run_ash(fit1, coef = tests[2])
  results_fit1 <- cbind(results, sebetahat = output_ash$data$s, lfsr = output_ash$result$lfsr,
                        lfdr = output_ash$result$lfdr, qvalue = output_ash$result$qvalue,
                        svalue = output_ash$result$svalue, beta_est = output_ash$result$PosteriorMean, se_est =
                          output_ash$result$PosteriorSD)
  
  
  # Scenario 2
  
  fit2 <- lmFit(resid_methyl, design)
  fit2 <- eBayes(fit2)
  
  results <- get_results(fit2, coef = tests[2])
  # Add mutliple testing correction with ASH
  output_ash <- run_ash(fit2, coef = tests[2])
  results_fit2 <- cbind(results, sebetahat = output_ash$data$s, lfsr = output_ash$result$lfsr,
                        lfdr = output_ash$result$lfdr, qvalue = output_ash$result$qvalue,
                        svalue = output_ash$result$svalue, beta_est = output_ash$result$PosteriorMean, se_est =
                          output_ash$result$PosteriorSD)
  
  # Bind DE from fit1 and DE from fit2 together
  
  HvC_Heart_fits12 <- as.data.frame(cbind(rownames(results_fit1), results_fit1$svalue, results_fit2$svalue), stringsAsFactors = FALSE)
  HvC_Heart_fits12[,2] <- as.numeric(HvC_Heart_fits12[,2])
  HvC_Heart_fits12[,3] <- as.numeric(HvC_Heart_fits12[,3])
  colnames(HvC_Heart_fits12) <- c("genes", "fit1_FDR", "fit2_FDR")
  
  DE_both_fit <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
  
  DE_before_fit <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR < FDR_level &  HvC_Heart_fits12$fit2_FDR > FDR_level),]
  
  DE_after_fit <- HvC_Heart_fits12[which(HvC_Heart_fits12$fit1_FDR > FDR_level &  HvC_Heart_fits12$fit2_FDR < FDR_level),]
}

  ### Permute the methylation values, then re-run
  
run1000 <- array(0, dim = c(3, 1000)) 
  
  #for (k in 1:1000){
    for (k in 1:2){
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
    #HvC_Heart_perm = topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="none")
    
    results <- get_results(fit2, coef = tests[2])
    # Add mutliple testing correction with ASH
    output_ash <- run_ash(fit2, coef = tests[2])
    results_fit2 <- cbind(results, sebetahat = output_ash$data$s, lfsr = output_ash$result$lfsr,
                          lfdr = output_ash$result$lfdr, qvalue = output_ash$result$qvalue,
                          svalue = output_ash$result$svalue, beta_est = output_ash$result$PosteriorMean, se_est =
                            output_ash$result$PosteriorSD)
    
    HvCHeart_DE_perm <- cbind(HvC_Heart_fits12$fit1_FDR,  results_fit2$svalue)
    
    DE_total <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] < FDR_level & HvCHeart_DE_perm[,2]  < FDR_level),])
    DE_total_report <- nrow(DE_total)
    
    DE_before <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1]  < FDR_level &  HvCHeart_DE_perm[,2]  > FDR_level),])
    DE_before_report <- nrow(DE_before)
    
    DE_after <- as.data.frame(HvCHeart_DE_perm[which(HvCHeart_DE_perm[,1] > FDR_level &  HvCHeart_DE_perm[,2] < FDR_level),])
    DE_after_report <- nrow(DE_after)
    
    new_numbers <- rbind(DE_total_report, DE_before_report, DE_after_report)
    print(new_numbers)
    run1000[,k] <- new_numbers
  }
  
  row1 <- median(run1000[1,1:1000])
  row2 <- median(run1000[2,1:1000])
  row3 <- median(run1000[3,1:1000])
  

  
  
  important_returns <- rbind(nrow(DE_both_fit), nrow(DE_before_fit), nrow(DE_after_fit), nrow(expression_values_only), row1, row2, row3)
  return(important_returns)
}


results_perm <- integration_ash(c(17, 20, 21, 24, 25, 28, 29), 0.05)
write.table(results_perm, "/mnt/gluster/home/leblake/Methylation/Integration_humans_chimps_robust/results/human_heart_kidney_FSR_with_perm_5perc.txt")