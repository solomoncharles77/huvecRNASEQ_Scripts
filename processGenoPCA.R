
# Process PCA file --------------------------------------------------------
# read in the eigenvectors, produced in PLINK
eigenvec <- read.table("../genotypeData/rawGeno/newImpute2025/huvecGeno.eigenvec", header = FALSE, skip=0, sep = ' ')
rownames(eigenvec) <- eigenvec[,2]
eigenvec <- eigenvec[,3:ncol(eigenvec)]
colnames(eigenvec) <- paste("PC", c(1:20), sep = '')
eigenvec <- data.frame(t(eigenvec))
eigenvec <- data.frame(Covariates = rownames(eigenvec), eigenvec)
write.table(eigenvec, "covFiles/huvecCovs.txt", sep = "\t")
