library(data.table)

knownEdit <- data.frame(fread("phenoFiles/filteredHuvecRediKnown.txt.gz"))
rownames(knownEdit) <- knownEdit$coordID
knownEdit$coordID <- NULL

# Calculate mean editing level per sample (column), ignoring NA
knownEditLevel <- colMeans(knownEdit, na.rm = TRUE)

# Convert to tidy data frame
kelDF <- data.frame(
  sample = names(knownEditLevel),
  overall_editing = round(knownEditLevel, 4),
  n_sites = colSums(!is.na(knownEdit)),  # number of non-NA sites per sample
  stringsAsFactors = FALSE)

head(kelDF)

############################################################
rnae <- data.frame(fread("phenoFiles/huvecRefSeqAEI_clean_AssocReady.txt"))
rnaeNov <- merge(rnae[, c(1,2,4)], kelDF, by.x = "V2", by.y = "sample")
head(rnaeNov)
cor(rnaeNov$V4, rnaeNov$overall_editing)

kelPheno <- rnaeNov[, c(2,1,4)]
fam <- read.table("../genotypeData/rawGeno/newImpute2025/huvecGeno.fam")
kelPheno <- kelPheno[match(fam$V2, kelPheno$V2), ]
colnames(kelPheno) <- c("FID",   "IID", "meanEditLevel")
head(kelPheno)
write.table(kelPheno, "phenoFiles/huvecMeanEditingLevel_clean_AssocReady.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(kelPheno, "phenoFiles/huvecMeanEditingLevel_clean_AssocReady_wtHeader.txt", row.names = F, quote = F, sep = "\t")

###############################################################
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


# normalize the data using rank-based inverse normal transformation-------
intTrans <- function(x){
  y <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(y)
}

kelPheno$meanEditLevel <-  intTrans(kelPheno$meanEditLevel)
kelPhenoNT <- shapiro.test(kelPheno$meanEditLevel)
ifelse(kelPhenoNT$p.value > 0.05, "normal", "not normal")


write.table(kelPheno, "phenoFiles/huvecMEL_Normalized_AssocReady.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(kelPheno, "phenoFiles/huvecMEL_Normalized_AssocReady_wtHeader.txt", row.names = F, quote = F, sep = "\t")
