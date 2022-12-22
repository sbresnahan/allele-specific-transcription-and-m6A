# Batch script to subset the output of EpinanoVariants to just the RRACH sites.
# Use like:
# Rscript src/epinano/subset_epinano_variants_batch.R \
#   results/minimap/2022-05-23/2C_all_sorted.plus_strand.per.site.csv \
#   results/epinano/CrossB_RRACH.rds \
#   results/epinano/2022-05-23/2C_all_sorted.plus_strand.RRACH.per.site.csv \
#   10

library(GenomicRanges)

infile <- commandArgs(trailingOnly = TRUE)[1]
grfile <- commandArgs(trailingOnly = TRUE)[2]
outfile <- commandArgs(trailingOnly = TRUE)[3]
minreads <- as.integer(commandArgs(trailingOnly = TRUE)[4])

gr <- readRDS(grfile)
tab <- read.csv(infile, check.names = FALSE)

# subset by coverage
tab <- tab[tab[[5]] >= minreads,]

# subset by RRACH positions
gr2 <- GRanges(seqnames = tab[[1]],
               ranges = IRanges(start = tab[[2]], width = 1))

fo <- findOverlaps(gr2, gr, type = "within")

ind <- sort(unique(queryHits(fo)))

write.csv(tab[ind,], file = outfile, row.names = FALSE, quote = FALSE)
