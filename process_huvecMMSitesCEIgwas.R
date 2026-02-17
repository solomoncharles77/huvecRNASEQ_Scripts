
library(data.table)
library(qqman)
library(tidyverse)

# Load data
gwas <- data.frame(fread("assocRes/huvecMMSitesCEIGWAS.mlma"))

# Rename columns
gwas <- setnames(gwas, c("Chr", "bp", "A2", "b", "se", "p"),
                 c("CHR", "BP", "REF", "BETA", "SE", "P"), skip_absent = T)
# Remove rows with missing p-values
gwas <- gwas[!is.na(gwas$P), ]
gwas$SNP <- sub("_.*", "", gwas$SNP)
gwas <- gwas[order(gwas$P), ]
head(gwas)

png("assocRes/manPlots/huvecMMSitesCEIGWAS_manhattan_plot.png", width=1200, height=600)
manhattan(gwas, main="huvecMMSitesCEIGWAS",
          suggestiveline = -log10(5e-05))
dev.off()
png("assocRes/qqPlots/huvecMMSitesCEIGWAS_qqplot.png", width=400, height=300) 
qq <- qq(gwas$P)
abline(h = -log10(5e-8), col = "red", lty = 2)
dev.off()  
fwrite(gwas, "resFiles/huvecMMSitesCEIGWAS_GWAS_Summary_Statistics.txt.gz")


########################################################################################
# Load data
gwas <- data.frame(fread("assocRes/huvecMMSitesCEIGWAS.mlma"))

# Rename columns
gwas <- setnames(gwas, c("Chr", "bp", "A2", "b", "se", "p"),
                 c("CHR", "BP", "REF", "BETA", "SE", "P"), skip_absent = T)
# Remove rows with missing p-values
gwas <- gwas[!is.na(gwas$P), ]
gwas <- gwas[order(gwas$P), ]
suggestive_hits <- gwas[gwas$P <= 0.00005, ]
head(suggestive_hits)

write.table(suggestive_hits, "assocRes/huvecMMSitesCEIGWAS_suggestiveHits.txt", row.names = F, quote = F, )
clumpHits <- read.table("assocRes/huvecMMSitesCEIGWAS_clumped_suggestive_loci.clumped", header = T)
gwasTop <- merge(clumpHits[, c(3,6)], gwas, by = "SNP")
gwasTop$rsID <- gwasTop$SNP <- sub("_.*", "", gwasTop$SNP)
gwasTop$otID <- paste0(gwasTop$CHR, "_", gwasTop$BP, "_",
                       gwasTop$REF, "_", gwasTop$A1)

# Define the list of rsIDs
rsIDs = gwasTop$rsID
otIDs <- gwasTop$otID

# Map to genes with biomartR
library("biomaRt")
mart <- useMart("ENSEMBL_MART_SNP")
dataset <- useDataset("hsapiens_snp", mart=mart)


# To get the ensembl gene id belonging to the SNPs
annotHit <- getBM(attributes=c('refsnp_id', 'ensembl_gene_stable_id',
                               "ensembl_transcript_stable_id"), 
                  filters = 'snp_filter', 
                  values = rsIDs, 
                  mart = dataset)

annotGWAS <- merge(gwasTop, annotHit, by.x = "rsID", by.y = "refsnp_id")

head(annotGWAS)

write.csv(gwasTop, "assocRes/huvecMMSitesCEIGWAS_gwasTop.csv", row.names = F )


library(otargen)
result <- variantEffectPredictorQuery(variantId = "1_154624272_T_C")
result <- result[, c(6, 9:10, 14:17)]
result <- result[result$target.biotype == "protein_coding", ]
result <- result[order(abs(result$distanceFromTss)), ]
result <- result[1, ]


otAnnot <- lapply(otIDs, function(x){
  cat(x, "\n")
  result <- variantEffectPredictorQuery(variantId = x)
  if (!is.null(result)) {
    result <- result[, c(6, 9:10, 14:17)]
    result <- result[result$target.biotype == "protein_coding", ]
    result <- result[order(abs(result$distanceFromTss)), ]
    result <- result[1, ]
    
  }else{
    result
  }
  
})

otAnnot <- do.call(rbind, otAnnot)
annotGWAS2 <- merge(annotGWAS, otAnnot, by.x = "otID", by.y = "variantId", all.x = T)

annotGWAS2$comGeneID <- ifelse(is.na(annotGWAS2$ensembl_gene_stable_id),
                               annotGWAS2$target.id, annotGWAS2$ensembl_gene_stable_id)
annotGWAS2$comTxID <- ifelse(is.na(annotGWAS2$ensembl_transcript_stable_id),
                             annotGWAS2$transcriptId, annotGWAS2$ensembl_transcript_stable_id)

write.csv(annotGWAS2, "resFiles/huvecMMSitesCEIGWAS_mappedGenes.csv", row.names = F)
write.table(data.frame(unique(annotGWAS2$comGeneID)), "resFiles/huvecMMSitesCEI_nearestGenesList.txt",
            row.names = F, col.names = F, quote = F)
write.table(data.frame(unique(annotGWAS2$comTxID)), "resFiles/huvecMMSitesCEI_nearestTranscriptList.txt",
            row.names = F, col.names = F, quote = F)

############################################################################
annotGWAS2 <- read.csv("resFiles/huvecMMSitesCEIGWAS_mappedGenes.csv")
gwas2 <- merge(annotGWAS2[, c(3,21)], gwas, by = "SNP", all.y = T)
gwas2$SNP <- ifelse(!is.na(gwas2$nearestGene.symbol), gwas2$nearestGene.symbol, gwas2$SNP)

genes_to_annotate <- c("ADAR","ADARB1")

toLabel <- gwas2[grepl("AD", gwas2$SNP), ]
png("assocRes/manPlots/huvecMMSitesCEI_manhattanPlot.png", width=1200, height=600)
manhattan(gwas2, main="",
          suggestiveline = -log10(5e-05),
          highlight = genes_to_annotate, 
          annotateTop = F)

text(
  x = toLabel$CHR,
  y = -log10(toLabel$P),
  labels = toLabel$SNP,
  pos = 3,      # Place text above the point
  cex = 2,    # !!! CUSTOM TEXT SIZE !!!
  font = 2,     # Bold font
  col = "blue4" # Custom color
)
dev.off()


