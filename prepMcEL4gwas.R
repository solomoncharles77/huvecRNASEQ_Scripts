library(data.table)

csEdCoord <- read.table("huvecRNASEQ_Scripts/codingSequenceEditingSitesCoordinates.txt")

# Import raw editing data
freqDF <- data.frame(fread("phenoFiles/rawHuvecRediKnown.txt.gz"))
colnames(freqDF)[1] <- "coordID"
rownames(freqDF) <- freqDF$coordID

cdsDF <- freqDF[freqDF$coordID %in% csEdCoord$hg38_ID, ]
cdsDF$coordID <- NULL

# Calculate mean editing level per sample (column), ignoring NA
cdsEditLevel <- colMeans(cdsDF, na.rm = TRUE)

# Convert to tidy data frame
kelDF <- data.frame(
  sample = names(cdsEditLevel),
  cds_editing = round(cdsEditLevel, 4),
  n_sites = colSums(!is.na(cdsDF)),  # number of non-NA sites per sample
  stringsAsFactors = FALSE)

head(kelDF)

############################################################
rnae <- data.frame(fread("phenoFiles/huvecMMSitesAEI_clean_AssocReady.txt"))
rnaeNov <- merge(rnae[, c(1,2,4)], kelDF, by.x = "V2", by.y = "sample")
head(rnaeNov)
cor(rnaeNov$V4, rnaeNov$cds_editing)

kelPheno <- rnaeNov[, c(2,1,4)]
fam <- read.table("../genotypeData/rawGeno/newImpute2025/huvecGeno.fam")
kelPheno <- kelPheno[match(fam$V2, kelPheno$V2), ]
colnames(kelPheno) <- c("FID",   "IID", "cdsEditLevel")
head(kelPheno)
write.table(kelPheno, "phenoFiles/huvecMcEL_clean_AssocReady.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(kelPheno, "phenoFiles/huvecMcEL_clean_AssocReady_wtHeader.txt", row.names = F, quote = F, sep = "\t")

# normalize the data using rank-based inverse normal transformation-------
intTrans <- function(x){
  y <- qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  return(y)
}

kelPheno$cdsEditLevel <-  intTrans(kelPheno$cdsEditLevel)
kelPhenoNT <- shapiro.test(kelPheno$cdsEditLevel)
ifelse(kelPhenoNT$p.value > 0.05, "normal", "not normal")

write.table(kelPheno, "phenoFiles/huvecMcEL_Normalized_AssocReady.txt", col.names = F, row.names = F, quote = F, sep = "\t")
write.table(kelPheno, "phenoFiles/huvecMcEL_Normalized_AssocReady_wtHeader.txt", row.names = F, quote = F, sep = "\t")
