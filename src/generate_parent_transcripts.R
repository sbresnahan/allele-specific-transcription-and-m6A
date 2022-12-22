# Use VCFs, genome annotation, and transcript sequences to build parental transcriptomes

library(GenomicFeatures)
library(VariantAnnotation)
library(ggplot2)

# Get genome annotation ####

#gtf0 <- makeTxDbFromGFF("data/reference/GCF_003254395.2_Amel_HAv3.1+virus_genomic.filled.gtf",
#                        format = "gtf") # error since some CDS do not match to exons

gtf0 <- rtracklayer::import("data/reference/GCF_003254395.2_Amel_HAv3.1+virus_genomic.filled.gtf")
gtf1 <- gtf0[gtf0$type == "exon"]

gtf1

any(is.na(gtf1$transcript_id)) # FALSE; all filled in

# list exons by transcript
txpt_grl <- GRangesList(sapply(unique(gtf1$transcript_id), function(x) gtf1[gtf1$transcript_id == x]))
#saveRDS(txpt_grl, file = "data/reference/exon_by_transcript_grangeslist.rds")
txpt_grl <- readRDS("data/reference/exon_by_transcript_grangeslist.rds")

table(lengths(txpt_grl))

# list transcripts
gtf2 <- gtf0[!is.na(gtf0$transcript_id) & !gtf0$type %in% c("exon", "CDS")] # not as many as in txpt_grl

# list unique exons
temp <- ranges(gtf1)
mcols(temp)$chrom <- seqnames(gtf1)
temp
duplicated(temp)
temp[1:9]

gtf3 <- gtf1[!duplicated(temp)] # eliminate duplicates
gtf3 <- reduce(gtf3) # merge overlapping ranges

# Preview VCF ####
# Four VCF files from Sean Bresnahan
vcffiles <- list.files("data/VCFs", pattern = ".*\\.vcf\\.gz$", full.names = TRUE)
names(vcffiles) <- sub("data/VCFs/", "", sub("\\.vcf\\.gz", "", vcffiles))

# for(v in vcffiles){
#   indexTabix(v, format = "vcf")
# }

test <- VcfFile(vcffiles[1])

hdr <- scanVcfHeader(test)

hdr
meta(hdr)
meta(hdr)$reference$Value # Amel_HAv3.1.fna
meta(hdr)$contig # uses NCBI accessions rather than chromosome numbers (yay)

seqlevels(txpt_grl)

# look at genotypes
View(as.data.frame(info(hdr)))

# Param to just import fields of interest, and just variants in exons
gtf3a <- gtf3[seqnames(gtf3) %in% seqlevels(hdr)]
seqlevels(gtf3a) <- seqlevels(hdr)

identical(as.character(seqnames(gtf3a)),
          as.character(seqnames(gtf3[seqnames(gtf3) %in% seqlevels(hdr)])))

mysvp_test <- ScanVcfParam(fixed = c("ALT", "QUAL"), info = NA, geno = c("GT"), # GQ is all missing data so don't bother
                      which = gtf3a[1:500])

test_vcf_dat <- readVcf(test, genome = seqinfo(hdr), param = mysvp_test)

rowRanges(test_vcf_dat)
geno(test_vcf_dat)$GT
table(geno(test_vcf_dat)$GT)
# 0/0 1/1
# 238 213

ggplot(mapping = aes(x = geno(test_vcf_dat)$GT,
                     y = log(rowRanges(test_vcf_dat)$QUAL))) +
  geom_violin()
# --> Anything with 0/0 genotype has quality below 1, anything with 1/1 has quality above 1

hist(log(rowRanges(test_vcf_dat)$QUAL)[geno(test_vcf_dat)$GT == "1/1"])

# Filter variants
# QUAL = 30 means 1 in 1000 chance of being incorrect
test_keep <- rowRanges(test_vcf_dat)[as.vector(geno(test_vcf_dat)$GT) == "1/1" &
                                       rowRanges(test_vcf_dat)$QUAL > 30]
test_keep

# Test conversion to transcript coordinates

test_keep_txpt <- mapToTranscripts(test_keep, txpt_grl)

test_keep_txpt
anyDuplicated(test_keep_txpt$xHits)

test_keep_txpt$ALT <- test_keep$ALT[test_keep_txpt$xHits]

lengths(test_keep$ALT)
unlist(test_keep$ALT[rep(list(1L), length(test_keep))])

# Test methods for modifying DNAstringset ####

test_dna <- DNAStringSet(c("AAAAAA", "CCCCCCCC", "GGGG", "TTTTT"))
names(test_dna) <- c("seq1", "seq2", "seq3", "seq4")
test_dna

test_rr <- GRanges(seqnames = c("seq1", "seq2", "seq2", "seq3"),
                   ranges = IRanges(start = c(3, 2, 5, 2),
                                    end = c(3, 2, 6, 2)),
                   ALT = DNAStringSet(c("C", "T", "G", "TT")))

test_dna[test_rr]

#test_dna[test_rr] <- DNAStringSet(c("C", "T", "G")) # doesn't work

test_irl <- IRangesList(sapply(seqlevels(test_rr), function(x) ranges(test_rr)[seqnames(test_rr) == x]))

test_irl

test_alt <- DNAStringSetList(sapply(seqlevels(test_rr), function(x) test_rr$ALT[seqnames(test_rr) == x]))

test_alt

test_dna2 <- test_dna
test_dna2[names(test_irl)] <- replaceAt(test_dna2[names(test_irl)], test_irl, test_alt)

test_dna2
test_dna

# Import transcriptome ####
my_transcriptome <- readDNAStringSet("data/reference/GCF_003254395.2_Amel_HAv3.1_rna+virus.fna")

# Function to import, filter, and convert variants and modify transcriptome ####
mysvp <- ScanVcfParam(fixed = c("ALT", "QUAL"), info = NA, geno = c("GT"),
                      which = gtf3a)

# Moved fn to generate_parent_transcript_batch.R to avoid tracking two copies

# bgzip to help indexing issues ####

for(v in vcffiles){
  bg <- bgzip(v)
  indexTabix(bg, format = "vcf4")
}

# test removing out of bound sequences ####
test_trim <- GRanges(c("chr1", "chr1", "chr2"),
                     IRanges(start = c(50, 100, 150),
                             end = c(80, 1600, 180)),
                     seqinfo = Seqinfo(c("chr1", "chr2"),
                                       c(1500, 1000)))
# --> Gives out of bounds warning

test_trim
trim(test_trim) # keeps but trims an out of bounds sequence

# I would rather remove entirely to avoid screwing up the sequence editing
end(test_trim) > seqlengths(seqinfo(test_trim))[as.integer(match(seqnames(test_trim), seqnames(seqinfo(test_trim))))]

# debug disjoint ranges ####
load("results/parent_transcriptomes/ranges_debug.RData")

irl[[32]]
as.data.frame(sort(irl[[32]]))

mychunk <- readVcf(vcffiles[1], genome = seqinfo(gtf3a),
                   param = ScanVcfParam(fixed = c("ALT", "QUAL"), info = NA, geno = c("GT"),
                                        which = GRanges("NC_037638.1",
                                                        IRanges(start = 1159000, end = 1160000))))
mychunk
rowRanges(mychunk)
# NC_037638.1:1159314 only in VCF once despite ending up in irl twice

mytxpt <- names(irl)[32] # XR_407457.3
txpt_grl[mytxpt] # two exons, each listed once

mymap <- mapToTranscripts(rowRanges(mychunk), txpt_grl)

mymap # 1159314 only listed once here

subsetByOverlaps(gtf3a,
                 GRanges("NC_037638.1",
                         IRanges(start = 1159000, end = 1160000)))
# --> two overlapping exons on opposite strands
# --> fixed by adding ignore.strand = TRUE to reduce

# Sanity check on transcript sequence ####
mychunk <- readVcf(vcffiles[1], genome = seqinfo(gtf1),
                   param = ScanVcfParam(fixed = c("ALT", "QUAL"), info = NA, geno = c("GT"),
                                        which = gtf1[gtf1$transcript_id == mytxpt]))
rowRanges(mychunk)

txpt_new <- my_transcriptome
names(txpt_new) <- sub(" .*", "", names(txpt_new))

mymap <- mapToTranscripts(rowRanges(mychunk), txpt_grl)
#mymap <- sort(mymap)
mymap

txpt_new[[mytxpt]][ranges(mymap)]

txpt_new[mymap]

for(i in seq_along(mymap)){
  variant <- names(mymap)[i]
  txpt_seq <- reverseComplement(txpt_new[[mytxpt]][ranges(mymap)[i]])
  gen_seq <- rowRanges(mychunk)[variant]$REF[[1]]
  if(txpt_seq == gen_seq){
    cat(paste0(variant, ": woot!\n"))
  } else {
    print(variant)
    print(txpt_seq)
    print(gen_seq)
  }
}
# --> all match

# Being on reverse strand, does transcript modification work?
irl_test <- irl[[mytxpt]][!duplicated(irl[[mytxpt]])]
irl_alt_test <- irl_alt[[mytxpt]][!duplicated(irl[[mytxpt]])]

txptseq_test <- replaceAt(txpt_new[[mytxpt]], irl_test, irl_alt_test)
txptseq_test
txpt_new[[mytxpt]] # --> No, it was done incorrectly

txptseq_test <- replaceAt(txpt_new[[mytxpt]], irl_test, reverseComplement(irl_alt_test))
txptseq_test
txpt_new[[mytxpt]]

# method to add to fn
mymap$ALT <- unlist(rowRanges(mychunk)$ALT[mymap$xHits][as.list(rep(1, 40))])
mymap
mymap$ALT[strand(mymap) == "-"] <- reverseComplement(mymap$ALT[strand(mymap) == "-"])
mymap

# Additional debugging of disjoint variants ####
load("results/parent_transcriptomes/SRR14654186_ranges_debug.RData")

irl[[2116]]
mytxpt <- names(irl)[2116]

subsetByOverlaps(gtf3a,
                 GRanges("NC_037638.1",
                         IRanges(start = 25046000, end = 25047000)))
# --> large variant that spans 2 exons.

# Spot check parental txptome after generation ####
var_by_txpt <- read.delim("results/parent_transcriptomes/SRR14654186_variants_by_txpt.txt")
head(var_by_txpt)
nrow(var_by_txpt) # 16153 transcripts with variants, out of 27893 transcripts total
plot(table(var_by_txpt$N_variants)) # exponential decrease, elbow at about 15
# this is EHB, so we might get more variants in AHB

txpt_SRR14654186 <- readDNAStringSet("results/parent_transcriptomes/SRR14654186_transcriptome.fa")

load("results/parent_transcriptomes/SRR14654186_ranges_debug.RData")

# grab 10 random transcripts to look at
set.seed(927)
txpts_examine <- sample(var_by_txpt$Transcript, 10)
var_by_txpt[match(txpts_examine, var_by_txpt$Transcript),]

gtf1[gtf1$transcript_id %in% txpts_examine]

# Import VCF in this regions
mychunk <- readVcf(vcffiles[1], genome = seqinfo(gtf1),
                   param = ScanVcfParam(fixed = c("ALT", "QUAL"), info = NA, geno = c("GT"),
                                        which = gtf1[gtf1$transcript_id %in% txpts_examine]))
mychunk <- mychunk[!as.vector(geno(mychunk)$GT) %in% c("0/0", "0") &
                     rowRanges(mychunk)$QUAL > 30,] # 105 variants

sum(var_by_txpt$N_variants[match(txpts_examine, var_by_txpt$Transcript)]) # 105

# Look at one txpt
i <- 6
x <- txpts_examine[i]

gtf1[gtf1$transcript_id == x] # note strand

rowRanges(mychunk)[names(irl[[x]])] # alleles and position

txpt_new[[x]]          # DNA sequences
txpt_SRR14654186[[x]]

txpt_new[[x]][irl[[x]]]
txpt_SRR14654186[[x]][irl[[x]]] # (this code doesn't work right with indels)

width(rowRanges(mychunk)[names(irl[[x]])]$REF)
width(unlist(rowRanges(mychunk)[names(irl[[x]])]$ALT))

aln <- pairwiseAlignment(txpt_new[[x]], txpt_SRR14654186[[x]])
aln
writePairwiseAlignments(aln, file = "results/parent_transcriptomes/test_align.txt")

# Debug non-disjoint ranges problem with male transcriptomes ####
infile <- vcffiles[2]

hdr <- headerTabix(infile)

gtf3a <- gtf3[seqnames(gtf3) %in% hdr$seqnames]
seqlevels(gtf3a) <- hdr$seqnames

load("results/parent_transcriptomes/SRR14654187_ranges_debug.RData")

mytxpt <- names(irl)[4159] # XM_026438858.1
as.data.frame(irl[[mytxpt]])
anyDuplicated(irl[[mytxpt]])

as.data.frame(unname(irl[[mytxpt]]))
anyDuplicated(as.data.frame(unname(irl[[mytxpt]])))
isDisjoint(irl[[mytxpt]])
db <- disjointBins(irl[[mytxpt]])
table(db)
irl[[mytxpt]][db == 2]

irl[[mytxpt]][70:71] # width of one is wrong
irl[[mytxpt]][260]

subsetByOverlaps(gtf1[gtf1$transcript_id == mytxpt],
                 GRanges("NC_037639.1", IRanges(start = 15716000, width = 2000)))
# NC_037639.1:15716973_TGATTTTTTTTTTTTTTTTTTTAGAAAG spanned an intron-exon junction;
# width got trimmed appropriately but of course the allele did not.

# test out fix
mysvp_test <- ScanVcfParam(fixed = c("ALT", "QUAL"), info = NA, geno = c("GT"),
                      which = subsetByOverlaps(gtf3a, GRanges("NC_037639.1", IRanges(start = 15716000, width = 2000))))
vf <- VcfFile(infile)
vcf_dat <- readVcf(vf, genome = seqinfo(gtf3a), param = mysvp_test)
keep <- !as.vector(geno(vcf_dat)$GT) %in% c("0/0", "0") &
  rowRanges(vcf_dat)$QUAL > 30
rr_keep <- rowRanges(vcf_dat)[keep]

"NC_037639.1:15716973_TGATTTTTTTTTTTTTTTTTTTAGAAAG/TGATTTTTTTTTTTTTTTTTAGAAAG" %in% names(rr_keep)

rr_keep <- subsetByOverlaps(rr_keep, gtf1, type = "within")

"NC_037639.1:15716973_TGATTTTTTTTTTTTTTTTTTTAGAAAG/TGATTTTTTTTTTTTTTTTTAGAAAG" %in% names(rr_keep)

# Setup to make parental transcriptomes for isoforms called by flair ####
#date <- "2021-10-20"
#date <- "2021-07-26"
#date <- "2021-10-25"
date <- "2021-12-20"

gtf0 <- rtracklayer::import(paste0("results/flair/", date, "/All_samples_collapse.isoforms.gtf"))
gtf1 <- gtf0[gtf0$type == "exon"]

gtf1

any(is.na(gtf1$transcript_id)) # FALSE; all filled in
table(nchar(gtf1$transcript_id))

sum(grepl("^GeneID", gtf1$gene_id)) # 116570 in Oct 25 dataset. Oct 20 had problems.
# --> 117882 in Dec 20 dataset.

# list exons by transcript
txpt_grl <- GRangesList(sapply(unique(gtf1$transcript_id), function(x) gtf1[gtf1$transcript_id == x]))
#saveRDS(txpt_grl, file = paste0("results/flair/", date, "/exon_by_transcript_grangeslist.rds"))
# --> 27,980 transcripts in Oct 25 dataset
length(unique(gtf1$transcript_id)) # also 27980
# 28283 in Dec2021

# explore data for debugging
gtf1[endsWith(as.character(seqnames(gtf1)), "V")] # 10 novel transcripts from three viral genomes

# Debugging 12/22/2021
txpt_grl <- readRDS(paste0("results/flair/", date, "/exon_by_transcript_grangeslist.rds"))
my_transcriptome <- readDNAStringSet(paste0("results/flair/", date, "/All_samples_collapse.isoforms.fa"))

txpt_new <- my_transcriptome
names(txpt_new) <- sub(" .*", "", names(txpt_new)) # for published transcriptome
names(txpt_new) <- sub("_GeneID:[[:digit:]]+$", "", names(txpt_new)) # Flair format
names(txpt_new) <- sub("_(VD|DW|ABP|CBP)V(:[[:digit:]]+)?$", "", names(txpt_new))
names(txpt_new) <- sub("_NC_[[:digit:]]{6}\\.1:[[:digit:]]+", "_NC", names(txpt_new))
names(txpt_new) <- sub("_NW_[[:digit:]]{9}\\.1:[[:digit:]]+", "_NW", names(txpt_new))
names(txpt_new) <- sub("_l-rRNA", "", names(txpt_new))

setequal(names(txpt_grl), names(txpt_new))
setdiff(names(txpt_grl), names(txpt_new))
setdiff(names(txpt_new), names(txpt_grl))
grep("a984a4fd-6c46-4a17-a1a4-f08abadc3387", names(my_transcriptome), value = TRUE)
