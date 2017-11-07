### Pulling orthologous CpGs for humans and chimps


####### Based on /mnt/lustre/home/jroux/Methylation/bsseq/smooth_data/combine_smoothed_data_commonSites.R
####### Data /mnt/lustre/home/jroux/Methylation/bsseq/smooth_data/combined_samples

library(bsseq)

allFiles <- dir("/mnt/lustre/home/jroux/Methylation/bsseq/smooth_data/combined_samples/separated_samples/", pattern = "_Smoothed.Rda", full.names = TRUE, ignore.case = TRUE)


trans <- NULL
parameters <- NULL
load(allFiles[1])
trans <- getBSseq(data.fit, "trans") ## the function is hard coded so it will be the same or all samples                
parameters <- getBSseq(data.fit, "parameters")
rm(data.fit) ## clean up workspace                                                                                      

load("/mnt/lustre/home/jroux/Methylation/bsseq/smooth_data/combined_samples/gr.RDa")
load("/mnt/lustre/home/jroux/Methylation/bsseq/smooth_data/combined_samples/Cov.RDa")
load("/mnt/lustre/home/jroux/Methylation/bsseq/smooth_data/combined_samples/M.RDa")

# Only want humans and chimps
load("/mnt/lustre/home/jroux/Methylation/bsseq/smooth_data/combined_samples/pData.RDa")
pData <- pData[-(33:48), ]

## Filtering of CpGs: keep sites seen in each of the 2 species                                                          

keepLoci <- which(rowSums(Cov[, grepl("^H", colnames(Cov), perl=T)] >= 1) >=1 & rowSums(Cov[, grepl("^C", colnames(Cov), perl=T)] >= 1) >=1)

M <- M[keepLoci,1:32]
gr <- gr[keepLoci,] # Note: doesn’t have species information so no need for the 1:32
Cov <- Cov[keepLoci,1:32]

## matrix of smoothed data                                                                                              
load("/mnt/lustre/home/jroux/Methylation/bsseq/smooth_data/combined_samples/coef.RDa")
coef <- coef[keepLoci,1:32]

# Save the data
allData.fit.subset.human.chimp <- BSseq(gr = gr, M = M, Cov = Cov, coef = coef, se.coef = NULL, pData = pData, trans = trans, rmZeroCov = FALSE)

save(allData.fit.subset.human.chimp, file = "../combinedSmoothedCommonHumanChimpSites.RDa")

Cov <- getCoverage(allData.fit.subset.human.chimp, type = "Cov")

## Yoav's suggestion for coverage filtering: 

## At least 2 out of 4 individuals per species/tissue with coverage >= 2x                                                                
keepLoci_cov <- which(
  ## Chimp                                                                                                               
  rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[1] & pData$Species == unique(pData$Species)[1]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[2] & pData$Species == unique(pData$Species)[1]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[3] & pData$Species == unique(pData$Species)[1]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[4] & pData$Species == unique(pData$Species)[1]] >= 2) >= 2 &
    ## Human                                                                                                               
    rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[1] & pData$Species == unique(pData$Species)[2]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[2] & pData$Species == unique(pData$Species)[2]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[3] & pData$Species == unique(pData$Species)[2]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[4] & pData$Species == unique(pData$Species)[2]] >= 2) >= 2 &
    
    ##Average coverage between 2x and 10x across all 4 individuals per species/tissue                                                     
    rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[1] & pData$Species == unique(pData$Species)[1]]) >= 2 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[2] & pData$Species == unique(pData$Species)[1]]) >= 2 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[3] & pData$Species == unique(pData$Species)[1]]) >= 2 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[4] & pData$Species == unique(pData$Species)[1]]) >= 2 &
    rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[1] & pData$Species == unique(pData$Species)[1]]) < 10 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[2] & pData$Species == unique(pData$Species)[1]]) < 10 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[3] & pData$Species == unique(pData$Species)[1]]) < 10 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[4] & pData$Species == unique(pData$Species)[1]]) < 10 &
    rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[1] & pData$Species == unique(pData$Species)[2]]) >= 2 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[2] & pData$Species == unique(pData$Species)[2]]) >= 2 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[3] & pData$Species == unique(pData$Species)[2]]) >= 2 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[4] & pData$Species == unique(pData$Species)[2]]) >= 2 &
    rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[1] & pData$Species == unique(pData$Species)[2]]) < 10 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[2] & pData$Species == unique(pData$Species)[2]]) < 10 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[3] & pData$Species == unique(pData$Species)[2]]) < 10 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[4] & pData$Species == unique(pData$Species)[2]]) < 10)

tab <- tab[keepLoci_cov,]
write.table(keepLoci_cov, file = "../combined_human_chimp_Loci_cov_2and10.txt")


#save(allData.fit.subset.filtered, file = "../combinedSmoothedCommonHumanChimpSites_Coverage_filtered_10mil.RDa")


############## 

# Get the smoothed methylation values
tab <- getMeth(allData.fit.subset.human.chimp, type = "smooth")
keepLoci_cov <- read.table("./combined_human_chimp_Loci_cov_2and10.txt")
keepLoci_cov <- as.matrix(keepLoci_cov)
tab_10mil <- tab[keepLoci_cov,]

# Get the values
cpg_250 <- read.table("./chimp_human_250_matrix_grab_250_merged.txt", header = T)

cpg_250_pos <- as.vector(cpg_250[,2])
