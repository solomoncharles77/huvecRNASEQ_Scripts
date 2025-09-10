
tg <- read.table("/scratch/vasccell/cs806/huvecRNASEQ/huvecRNASEQ_Scripts/target.txt")
tg
tg$group <- as.factor(c("S", "S", "C", "C", "S", "S", "S", "C", "C", "C")) # edit as required.


write.table(tg[, c(1:3)], "/scratch/vasccell/cs806/huvecRNASEQ/huvecRNASEQ_Scripts/target.txt", quote = F, sep="\t")

targetFile <- "/scratch/vasccell/cs806/huvecRNASEQ/huvecRNASEQ_Scripts/target.txt"
target <- read.table(targetFile, header=TRUE, sep="\t", na.strings="")
