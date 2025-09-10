
qc <- read.delim("/scratch/vasccell/cs806/huvecRNASEQ/fastqcResults/multiqc_data/multiqc_fastqc.txt")
qc <- qc[,  c("Sample", "Filename", "avg_sequence_length", "Total.Sequences", "basic_statistics")]
colnames(qc) <- c("sampleID", "filename", "avg_read_length", "total_reads", "basic_statistics")
write.csv(qc, "/scratch/vasccell/cs806/huvecRNASEQ/fastqcResults/qcSummary.csv", row.names = F)

map <- read.delim("/scratch/vasccell/cs806/huvecRNASEQ/mappedReads/multiqc_data/multiqc_star.txt")
# map <- map[, c("Sample", "total_reads", "uniquely_mapped", "uniquely_mapped_percent", "unmapped_tooshort", "unmapped_tooshort_percent")]
# colnames(map) <- c("Sample", "total_reads", "mapped", "percent_mapped", "unmapped", "percent_unmapped")
map <- map[, c("Sample", "total_reads", "uniquely_mapped", "uniquely_mapped_percent")]
colnames(map) <- c("Sample", "total_reads", "mapped", "percent_mapped")
map <- map[!grepl("_STARpass1", map$Sample), ]
write.csv(map, "/scratch/vasccell/cs806/huvecRNASEQ/mappedReads/mapSummary.csv", row.names = F)

salm <- read.delim("/scratch/vasccell/cs806/huvecRNASEQ/salmonQuant/multiqc_data/multiqc_salmon.txt")
salm <- salm[, c("Sample", "num_processed", "num_mapped", "percent_mapped")]
colnames(salm) <- c("Sample", "total_reads", "mapped", "percent_mapped")
write.csv(salm, "/scratch/vasccell/cs806/huvecRNASEQ/salmonQuant/quantSummary.csv", row.names = F)
