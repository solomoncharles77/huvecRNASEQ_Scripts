library(tximport)


tx2gene <- read.csv("/scratch/vasccell/cs806/exprPhenoData/tx2geneHomo_sapiens.GRCh38.cdna.all.csv")

samples <- read.table("/scratch/vasccell/cs806/huvecRNASEQ/huvecRNASEQ_Scripts/allSeq3.txt", header = F)
files <- file.path("/scratch/vasccell/cs806/huvecRNASEQ/salmonQuant", samples$V1, "quant.sf")
names(files) <- paste0(samples$V1)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
saveRDS(txi.salmon, file = "/scratch/vasccell/cs806/huvecRNASEQ/salmonQuant/salmonTxi.rds")
