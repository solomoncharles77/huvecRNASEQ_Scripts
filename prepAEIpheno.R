

rnae <- read.csv("RNAEditingIndexer/AEI/huvecRNASeqSUM/EditingIndex.csv")
rnae$Sample <- sub("Aligned.out.sam.sorted", "", rnae$Sample)

# Harmonize names and samples in geno and expr ----------------------------
mapFile <- read.csv("huvecRNASEQ_Scripts/Mappings.csv")
rnae$Sample <- mapFile$supplier_name[match(rnae$Sample, mapFile$sanger_sample_id)]
rnae$Sample <- sub("E", "S", rnae$Sample)
rnae[1:10, c(1,3,5:10)]

rnaeRand <- rnae[rnae$StrandDecidingMethod == "Randomly", ] 
rnaeRef <- rnae[rnae$StrandDecidingMethod == "RefSeqThenMMSites", ] 
rnaeMMS <- rnae[rnae$StrandDecidingMethod == "MMSitesThenRefSeq", ]
write.table(data.frame(rnaeMMS$Sample),
            "../genotypeData/genotypeData_Scripts/huvecRNAseqSamps.txt",
            row.names = F, col.names = F, quote = F)

##########################################################################
fam <- read.table("../genotypeData/rawGeno/newImpute2025/huvecGeno.fam")

rnaeRand <- rnaeRand[, -c(1,2,4)]
# Reorder rnaeRand to match the order of pheno
rnaeRand <- rnaeRand[match(fam$V2, rnaeRand$Sample), ]
rnaeRand[1:10, 1:10]
head(fam)
rnaeRand <- data.frame(fam[, c(1,2)], rnaeRand[, c(2:7)])
colnames(rnaeRand)[1:2] <- c("FID",   "IID")
head(rnaeRand)
write.table(rnaeRand, "phenoFiles/huvecRandomAEI_clean_AssocReady.txt", col.names = F, row.names = F, quote = F, sep = "\t")

##########################################################################################################
# prep covariate
genoPCA <- read.table(paste0("covFiles/huvecCovs.txt"), check.names=FALSE)
genoPCA$Covariates <- NULL
genoPCA <- t(genoPCA[1:3, ])
cov <- data.frame(fam[, c(1,2,5)], genoPCA)
colnames(cov)[1:3] <- c("FID",   "IID", "Sex")

write.table(cov, "covFiles/huvec_Sex_3gpc_covariates.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(cov, "covFiles/huvec_Sex_3gpc_covariates_wtHeader.txt", row.names = F, quote = F, sep = "\t")
########################################################################################################
########################################################################################################


rnaeRef <- rnaeRef[, -c(1,2,4)]
# Reorder rnaeRef to match the order of pheno
rnaeRef <- rnaeRef[match(fam$V2, rnaeRef$Sample), ]
rnaeRef[1:10, 1:10]
head(fam)
rnaeRef <- data.frame(fam[, c(1,2)], rnaeRef[, c(2:7)])
colnames(rnaeRef)[1:2] <- c("FID",   "IID")
head(rnaeRef)
write.table(rnaeRef, "phenoFiles/huvecRefSeqAEI_clean_AssocReady.txt", col.names = F, row.names = F, quote = F, sep = "\t")

###################################################
rnaeMMS <- rnaeMMS[, -c(1,2,4)]
# Reorder rnaeMMS to match the order of pheno
rnaeMMS <- rnaeMMS[match(fam$V2, rnaeMMS$Sample), ]
rnaeMMS[1:10, 1:10]
head(fam)
rnaeMMS <- data.frame(fam[, c(1,2)], rnaeMMS[, c(2:7)])
colnames(rnaeMMS)[1:2] <- c("FID",   "IID")
head(rnaeMMS)
write.table(rnaeMMS, "phenoFiles/huvecMMSitesAEI_clean_AssocReady.txt", col.names = F, row.names = F, quote = F, sep = "\t")
