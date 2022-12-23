# Identify all possible m6A (RRACH) sites in the transcriptome.

library(Biostrings)
library(GenomicRanges)

Mine_m6A <- function(file, pat = "RRACH"){
  txptome <- readDNAStringSet(file)
  mtch <- vmatchPattern(pat, txptome, fixed = FALSE)
  gr <- GRanges(mtch)
  return(gr)
}

grA <- Mine_m6A("results/parent_transcriptomes/2022-06-14/CrossA_transcriptome.fa")
grB <- Mine_m6A("results/parent_transcriptomes/2022-06-14/CrossB_transcriptome.fa")

# Narrow down to the central A

grAn <- narrow(grA, start = 3, end = 3)
grAn

grBn <- narrow(grB, start = 3, end = 3)

saveRDS(grAn, file = "results/epinano/CrossA_RRACH_Aonly_2022-07-06.rds")
saveRDS(grBn, file = "results/epinano/CrossB_RRACH_Aonly_2022-07-06.rds")
