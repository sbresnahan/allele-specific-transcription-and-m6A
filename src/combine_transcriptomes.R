# Combine the two parental transcriptomes for each cross

# SRA         Subspecies Sex Cross
# SRR14654186 EHB        F   B
# SRR14654187 AHB        M   B
# SRR14654188 AHB        F   A
# SRR14654189 EHB        M   A

library(Biostrings)

date <- "2022-06-14"

txptome_86 <- readDNAStringSet(paste0("results/parent_transcriptomes/", date, "/SRR14654186_transcriptome.fa"))
txptome_87 <- readDNAStringSet(paste0("results/parent_transcriptomes/", date, "/SRR14654187_transcriptome.fa"))
txptome_88 <- readDNAStringSet(paste0("results/parent_transcriptomes/", date, "/SRR14654188_transcriptome.fa"))
txptome_89 <- readDNAStringSet(paste0("results/parent_transcriptomes/", date, "/SRR14654189_transcriptome.fa"))

mean(txptome_86 == txptome_87) # 28.6% transcripts identical in cross B
mean(txptome_88 == txptome_89) # 28.6% transcripts identical in cross A

names(txptome_86) <- paste0(ifelse(txptome_86 == txptome_87, "Both-", "SRR14654186-"), names(txptome_86))
names(txptome_87) <- paste0(ifelse(txptome_87 == txptome_86, "Both-", "SRR14654187-"), names(txptome_87))

txptome_B <- c(txptome_86, txptome_87[txptome_87 != txptome_86])

names(txptome_88) <- paste0(ifelse(txptome_88 == txptome_89, "Both-", "SRR14654188-"), names(txptome_88))
names(txptome_89) <- paste0(ifelse(txptome_89 == txptome_88, "Both-", "SRR14654189-"), names(txptome_89))

txptome_A <- c(txptome_88, txptome_89[txptome_89 != txptome_88])

# sanity checking
hist(log(width(txptome_86[txptome_86 == txptome_87])))
hist(log(width(txptome_86[txptome_86 != txptome_87]))) # longer in general

# export

writeXStringSet(txptome_A, paste0("results/parent_transcriptomes/", date, "/CrossA_transcriptome.fa"))
writeXStringSet(txptome_B, paste0("results/parent_transcriptomes/", date, "/CrossB_transcriptome.fa"))
