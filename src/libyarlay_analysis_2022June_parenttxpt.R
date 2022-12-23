library(limma)
library(edgeR)
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(dplyr)
library(Glimma)
library(seriation)
library(viridis)
library(RColorBrewer)
library(Biostrings)
library(msa)
library(seqinr)
library(ape)

# Import and initial QC ####

# Sample info
targets <- read.delim("src/targets.txt")
targets <- targets[targets$Use,]
targets <- targets[order(targets$Cross, targets$Aggression),]
targets$Name[c(1:3, 6)] <- c("T3-Ag1:A2C", "T3-Ag1:A7C", "T3-Ag1:A9C", "T3-Ag2:A1T")
targets$Label <- targets$Name
targets$Label <- factor(targets$Label, levels = targets$Label)
targets$Group <- paste0(targets$Cross, targets$Aggression)
targets$FlairName <- paste(targets$Name, targets$Aggression,
                           sub("^[[:digit:]]{1,2}-", "", sub("-[[:digit:]]{2}$", "", targets$Batch)),
                           sep = "_")

colorkey <- brewer.pal(4, "Paired")
names(colorkey) <- unique(targets$Group)

# Import GTF
gtf0 <- readRDS("results/final_report/annotations.rds")
gtf1 <- gtf0[gtf0$type == "gene"]

# Genomic coordinates of isoforms
isoform_anno <- read.table("results/flair/2022-06-14/All_samples_collapse.isoforms.bed",
                           sep = "\t", header = FALSE)

# Import results from flair
countsmatA <- as.matrix(read.table("results/flair/2022-06-14/counts_matrix_parenttxpt_A.tsv",
                                   row.names = 1, header = TRUE, check.names = FALSE))

countsmatB <- as.matrix(read.table("results/flair/2022-06-14/counts_matrix_parenttxpt_B.tsv",
                                   row.names = 1, header = TRUE, check.names = FALSE))

countsmatA[1:10,]
tail(countsmatA)
str(countsmatA)
table(countsmatA %% 1) # all zero because all are integers
table(countsmatB %% 1)

mode(countsmatA) <- "integer"
mode(countsmatB) <- "integer"

setequal(c(colnames(countsmatA), colnames(countsmatB)), targets$FlairName) # TRUE

# Regular expressions for splitting up transcript labels
readidpat <- "^[[:xdigit:]]{8}-([[:xdigit:]]{4}-){3}[[:xdigit:]]{12}(-[[:digit:]])?"
txptpat <- "^[NX][MR]_[[:digit:]]+\\.[[:digit:]]+(-[[:digit:]])?"
pseudogene_pat <- "^gene[[:digit:]]+(-[[:digit:]])?"

unknowngene_pat <- "N[CW]_[[:digit:]]+\\.[[:digit:]]+:[[:digit:]]+$"
knowngene_pat <- "GeneID:[[:digit:]]+$"
viralgene_pat <- "(DW|ABP|IAP|SB|CBP|VD)V(:[[:digit::]+)?$"
rRNA_pat <- "l-rRNA$"

# Function to make a data frame of information about each transcript
makeTxptInfo <- function(countsmat){
  txptinfo <- data.frame(Label = rownames(countsmat),
                         FlairLab = sub("(Both|SRR1465418[6-9])-", "", rownames(countsmat)),
                         Parent = sub("-.*$", "", rownames(countsmat)),
                         Transcript = "",
                         Gene = "",
                         TotDepth = rowSums(countsmat),
                         row.names = NULL)
  
  stopifnot(all(grepl(readidpat, txptinfo$FlairLab) |
                  grepl(txptpat, txptinfo$FlairLab) |
                  grepl(pseudogene_pat, txptinfo$FlairLab)))
  stopifnot(all(grepl(unknowngene_pat, txptinfo$FlairLab) |
                  grepl(knowngene_pat, txptinfo$FlairLab) |
                  grepl(viralgene_pat, txptinfo$FlairLab) |
                  grepl(rRNA_pat, txptinfo$FlairLab)))
  
  # Add transcript and gene ID to table
  temp <- grep(paste0(readidpat, "_", unknowngene_pat), txptinfo$FlairLab)
  txptinfo$Transcript[temp] <- sub(paste0("_", unknowngene_pat), "", txptinfo$FlairLab[temp])
  txptinfo$Gene[temp] <- sub(paste0(readidpat, "_"), "", txptinfo$FlairLab[temp])
  
  temp <- grep(paste0(readidpat, "_", knowngene_pat), txptinfo$FlairLab)
  txptinfo$Transcript[temp] <- sub(paste0("_", knowngene_pat), "", txptinfo$FlairLab[temp])
  txptinfo$Gene[temp] <- sub(paste0(readidpat, "_"), "", txptinfo$FlairLab[temp])
  
  temp <- grep(paste0(txptpat, "_", knowngene_pat), txptinfo$FlairLab)
  txptinfo$Transcript[temp] <- sub(paste0("_", knowngene_pat), "", txptinfo$FlairLab[temp])
  txptinfo$Gene[temp] <- sub(paste0(txptpat, "_"), "", txptinfo$FlairLab[temp])
  
  temp <- grep(paste0(pseudogene_pat, "_", knowngene_pat), txptinfo$FlairLab)
  txptinfo$Transcript[temp] <- sub(paste0("_", knowngene_pat), "", txptinfo$FlairLab[temp])
  txptinfo$Gene[temp] <- sub(paste0(pseudogene_pat, "_"), "", txptinfo$FlairLab[temp])
  
  temp <- grep(paste0(readidpat, "_", viralgene_pat), txptinfo$FlairLab)
  txptinfo$Transcript[temp] <- sub(paste0("_", viralgene_pat), "", txptinfo$FlairLab[temp])
  txptinfo$Gene[temp] <- sub(paste0(readidpat, "_"), "", txptinfo$FlairLab[temp])
  
  temp <- grep(paste0(readidpat, "_", rRNA_pat), txptinfo$FlairLab)
  txptinfo$Transcript[temp] <- sub(paste0("_", rRNA_pat), "", txptinfo$FlairLab[temp])
  txptinfo$Gene[temp] <- sub(paste0(readidpat, "_"), "", txptinfo$FlairLab[temp])
  
  # Sanity check
  stopifnot(!any(txptinfo$Gene == ""))
  stopifnot(!any(txptinfo$Transcript == ""))
  
  # Add gene symbol
  txptinfo$Symbol <- gtf1$gene[match(txptinfo$Gene, gtf1$gene_id)]
  
  temp <- match(txptinfo$FlairLab, isoform_anno$V4)
  stopifnot(!any(is.na(temp)))
  txptinfo$Chromosome <- isoform_anno$V1[temp]
  txptinfo$Start <- isoform_anno$V2[temp]
  txptinfo$End <- isoform_anno$V3[temp]
  
  return(txptinfo)
}

txptinfoA <- makeTxptInfo(countsmatA)
txptinfoB <- makeTxptInfo(countsmatB)

# reads in annotated vs. unannotated txpts
sum(txptinfoA$TotDepth[grepl(readidpat, txptinfoA$Transcript)])  # 1000207 reads in unannotated txpts
sum(txptinfoA$TotDepth[!grepl(readidpat, txptinfoA$Transcript)]) # 1392946 reads in annotated txpts

sum(txptinfoB$TotDepth[grepl(readidpat, txptinfoB$Transcript)])  # 1165883 reads in unannotated txpts
sum(txptinfoB$TotDepth[!grepl(readidpat, txptinfoB$Transcript)]) # 2349497 reads in annotated txpts

# look at counts in viral genomes
viral <- unique(as.character(seqnames(gtf0)[gtf0$source == "VirusChrom"]))
viralrowsA <- sapply(viral, function(x) grep(x, txptinfoA$Label))
viralrowsA <- viralrowsA[lengths(viralrowsA) > 0]
viralrowsB <- sapply(viral, function(x) grep(x, txptinfoB$Label))
viralrowsB <- viralrowsB[lengths(viralrowsB) > 0]

viralrowsA
viralrowsB
txptinfoA[unlist(viralrowsA),]
countsmatA[unlist(viralrowsA),] # T3-Ag1:A7C has high VDV load
txptinfoB[unlist(viralrowsB),]
countsmatB[unlist(viralrowsB),] # T2-Ag1:B2C has high VDV load

# QC stats for samples ###

rowsA <- match(colnames(countsmatA), targets$FlairName)
rowsB <- match(colnames(countsmatB), targets$FlairName)

targets$Reads_aligned <- 0L
targets$Reads_aligned[rowsA] <- colSums(countsmatA)
targets$Reads_aligned[rowsB] <- colSums(countsmatB)

ggplot(targets, aes(x = Label, y = Reads_aligned, fill = Group)) +
  geom_col() + ggtitle("Total aligned reads") # Newer samples have smaller lib size

targets$Reads_known_txpt_nuclear <- 0L
targets$Reads_known_txpt_nuclear[rowsA] <- colSums(countsmatA[grepl(txptpat, txptinfoA$Transcript) &
                                                                txptinfoA$Chromosome != "NC_001566.1",])
targets$Reads_known_txpt_nuclear[rowsB] <- colSums(countsmatB[grepl(txptpat, txptinfoB$Transcript) &
                                                                txptinfoB$Chromosome != "NC_001566.1",])

targets$Reads_viral <- 0L
targets$Reads_viral[rowsA] <- colSums(countsmatA[grep(viralgene_pat, txptinfoA$Gene),])
targets$Reads_viral[rowsB] <- colSums(countsmatB[grep(viralgene_pat, txptinfoB$Gene),])

targets$Reads_known_gene_novel_isoform_nuclear <- 0L
targets$Reads_known_gene_novel_isoform_nuclear[rowsA] <-
  colSums(countsmatA[grepl(knowngene_pat, txptinfoA$Gene) & grepl(readidpat, txptinfoA$Transcript) &
                       txptinfoA$Chromosome != "NC_001566.1",])
targets$Reads_known_gene_novel_isoform_nuclear[rowsB] <-
  colSums(countsmatB[grepl(knowngene_pat, txptinfoB$Gene) & grepl(readidpat, txptinfoB$Transcript) &
                       txptinfoB$Chromosome != "NC_001566.1",])

targets$Reads_pseudogene <- 0L
targets$Reads_pseudogene[rowsA] <- colSums(countsmatA[grep(pseudogene_pat, txptinfoA$Transcript),])
targets$Reads_pseudogene[rowsB] <- colSums(countsmatB[grep(pseudogene_pat, txptinfoB$Transcript),])

# targets$Reads_rRNA <- 0L  ## All mitochondrial
# targets$Reads_rRNA[rowsA] <- colSums(countsmatA[grep(rRNA_pat, txptinfoA$Gene),])
# targets$Reads_rRNA[rowsB] <- colSums(countsmatB[grep(rRNA_pat, txptinfoB$Gene),])

targets$Reads_unknown_gene_nuclear <- 0L
targets$Reads_unknown_gene_nuclear[rowsA] <-
  colSums(countsmatA[grepl(unknowngene_pat, txptinfoA$Gene) &
                       txptinfoA$Chromosome != "NC_001566.1",])
targets$Reads_unknown_gene_nuclear[rowsB] <-
  colSums(countsmatB[grepl(unknowngene_pat, txptinfoB$Gene) &
                       txptinfoB$Chromosome != "NC_001566.1",])

targets$Reads_mitochondrial <- 0L
targets$Reads_mitochondrial[rowsA] <- colSums(countsmatA[txptinfoA$Chromosome == "NC_001566.1",])
targets$Reads_mitochondrial[rowsB] <- colSums(countsmatB[txptinfoB$Chromosome == "NC_001566.1",])

targets$Reads_NC_037641.1_1186000 <- 0L
targets$Reads_NC_037641.1_1186000[rowsA] <- colSums(countsmatA[txptinfoA$Gene == "NC_037641.1:1186000",])
targets$Reads_NC_037641.1_1186000[rowsB] <- colSums(countsmatB[txptinfoB$Gene == "NC_037641.1:1186000",])

targets$Reads_unknown_gene_nuclear <-
  targets$Reads_unknown_gene_nuclear - targets$Reads_NC_037641.1_1186000

all(rowSums(targets[,c("Reads_known_txpt_nuclear", "Reads_viral", "Reads_known_gene_novel_isoform_nuclear",
                       "Reads_pseudogene", "Reads_unknown_gene_nuclear", "Reads_mitochondrial",
                       "Reads_NC_037641.1_1186000")]) ==
      targets$Reads_aligned) # TRUE

ggplot(targets, aes(x = Label, y = Reads_mitochondrial / Reads_aligned * 100, fill = Group)) +
  geom_col() +
  labs(y = "Percentage",
       title = "Percentage of mitochondrial reads in each sample")

# Build longer table for plots
targ_long <- targets %>% select(Label, Group, Batch, starts_with("Reads")) %>%
  pivot_longer(cols = all_of(c("Reads_known_txpt_nuclear", "Reads_viral", "Reads_known_gene_novel_isoform_nuclear",
                               "Reads_unknown_gene_nuclear", "Reads_mitochondrial", "Reads_NC_037641.1_1186000")),
               names_to = "Fate", values_to = "Reads") %>%
  mutate(Fate = factor(sub("txpt", "transcript", sub("Reads_", "", Fate)),
                       levels = c("mitochondrial", "viral", "NC_037641.1_1186000", "unknown_gene_nuclear",
                                  "known_gene_novel_isoform_nuclear", "known_transcript_nuclear")))

ggplot(targ_long, aes(x = Label, y = Reads, fill = Fate)) +
  geom_col() +
  scale_fill_manual(values = dittoSeq::dittoColors(1)) +
  ggtitle("Read fates expressed as number of reads") ## Use this fig for overall library size

ggplot(targ_long, aes(x = Label, y = Reads / Reads_aligned * 100, fill = Fate)) +
  geom_col() +
  scale_fill_manual(values = dittoSeq::dittoColors(1)) +
  ggtitle("Read fates expressed as percentage of reads") +
  labs(y = "Percentage") ## Use this fig for proportions in various categories

# --> Since fixing the Flair bug for bottom strand transcripts, cross A now has
# many more unknown nuclear gene reads, and many fewer mitochondrial reads, than
# cross B.

# Debug mitochondrial vs. nuclear difference in cross A ####

sum(txptinfoA$Chromosome == "NC_001566.1") # 109 mitochondrial txpts
sum(txptinfoB$Chromosome == "NC_001566.1") # 118 mitochondrial txpts

txptinfoA %>%
  filter(Chromosome == "NC_001566.1") %>%
  arrange(desc(TotDepth)) %>%
  head()

txptinfoB %>%
  filter(Chromosome == "NC_001566.1") %>%
  arrange(desc(TotDepth)) %>%
  head()

txptinfoA %>%
  filter(grepl(unknowngene_pat, Label)) %>%
  arrange(desc(TotDepth)) %>%
  head()

txptinfoB %>%
  filter(grepl(unknowngene_pat, Label)) %>%
  arrange(desc(TotDepth)) %>%
  head()

# Transcript SRR14654188-f7f95e2b-d405-4842-a137-20138e4d467e_NC_037641.1:1186000
# is very highly detected in cross A, but equivalent not in cross B.

txptomeA <- readDNAStringSet("results/parent_transcriptomes/2022-06-14/CrossA_transcriptome.fa")
txptomeB <- readDNAStringSet("results/parent_transcriptomes/2022-06-14/CrossB_transcriptome.fa")

# nuclear txpt
huhA <- txptomeA[grep("f7f95e2b-d405-4842-a137-20138e4d467e", names(txptomeA))]
huhB <- txptomeB[grep("f7f95e2b-d405-4842-a137-20138e4d467e", names(txptomeB))]

# mito txpt
bwahA <- txptomeA[grep("21305644-1911-4f7a-a81a-1e77ae9613dc", names(txptomeA))]
bwahB <- txptomeB[grep("21305644-1911-4f7a-a81a-1e77ae9613dc", names(txptomeB))]

align_nuc_mt <- msa(c(huhA, huhB, reverseComplement(c(bwahA, bwahB))))
align_nuc_mt

print(align_nuc_mt, show="complete")

align_nuc_mt2 <- msaConvert(align_nuc_mt, typ = "seqinr::alignment")

dist_nuc_mt <- dist.alignment(align_nuc_mt2, matrix = "identity", gap = FALSE)

nj_nuc_mt <- nj(dist_nuc_mt)

plot(nj_nuc_mt, cex = 0.5) # EHB sorts into nuclear vs. mitochondrial, but AHB grouped by sample.

# Next, check which strand these were on.
isoform_anno[grep("f7f95e2b-d405-4842-a137-20138e4d467e", isoform_anno$V4),] # nuclear on minus
isoform_anno[grep("21305644-1911-4f7a-a81a-1e77ae9613dc", isoform_anno$V4),] # mito on minus

# Sequence to BLAST
as.character(huhA[1]) ## --> All top hits are honeybee mitochondrion, covers 100% of query.
as.character(huhB[1]) ## --> also honeybee mitochondrion

# In genome browser, nothing really annotated in that range for nuclear.
# Spans ND1 and l-rRNA on mitochondrial genome.

gtf_flair <- rtracklayer::import("results/flair/2022-06-14/All_samples_collapse.isoforms.gtf")
gtf_flair[grep("f7f95e2b-d405-4842-a137-20138e4d467e", gtf_flair$transcript_id)]
# --> One long exon

# As a last sanity check, extract this region from the original genome
checkseq <- Rsamtools::scanFa("data/reference/GCF_003254395.2_Amel_HAv3.1+virus_genomic.fna",
                              param = GRanges("NC_037641.1", IRanges(1186491, 1189696), strand = "-"))
reverseComplement(checkseq) # yes, seems to match txptome.

# Confirm mitochondrial reads assigned to mother ####
txptinfoA$Mito <- txptinfoA$Chromosome == "NC_001566.1" | txptinfoA$Gene == "NC_037641.1:1186000"
prop_mito_A <- apply(countsmatA, 2,
                     function(x){
                       tot <- sum(x)
                       sapply(unique(txptinfoA$Parent),
                              function(y) sum(x[txptinfoA$Mito & txptinfoA$Parent == y]) / tot )
                     })
prop_mito_A

txptinfoB$Mito <- txptinfoB$Chromosome == "NC_001566.1" | txptinfoB$Gene == "NC_037641.1:1186000"
prop_mito_B <- apply(countsmatB, 2,
                     function(x){
                       tot <- sum(x)
                       sapply(unique(txptinfoB$Parent),
                              function(y) sum(x[txptinfoB$Mito & txptinfoB$Parent == y]) / tot )
                     })
prop_mito_B

barplot(prop_mito_A, beside = TRUE)

data.frame(Parent = rownames(prop_mito_A),
           prop_mito_A, check.names = FALSE) %>%
  pivot_longer(starts_with("T"), names_to = "Sample",
               values_to = "Proportion_mitochondrial_reads") %>%
  ggplot(aes(x = Sample, fill = Parent, y = Proportion_mitochondrial_reads)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_brewer(palette = "Dark2")

data.frame(Parent = rownames(prop_mito_B),
           prop_mito_B, check.names = FALSE) %>%
  pivot_longer(starts_with("T"), names_to = "Sample",
               values_to = "Proportion_mitochondrial_reads") %>%
  ggplot(aes(x = Sample, fill = Parent, y = Proportion_mitochondrial_reads)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_brewer(palette = "Dark2")

# Summarize to gene ####
allgenes <- sort(unique(c(txptinfoA$Gene, txptinfoB$Gene))) # 19920 genes
length(setdiff(txptinfoA$Gene, txptinfoB$Gene))   #   269 in A but not B
length(setdiff(txptinfoB$Gene, txptinfoA$Gene))   #  1613 in B but not A
length(intersect(txptinfoA$Gene, txptinfoB$Gene)) # 18038 in both

countsmat_gene <- matrix(0L, nrow = length(allgenes), ncol = ncol(countsmatA) + ncol(countsmatB),
                         dimnames = list(allgenes, c(colnames(countsmatA), colnames(countsmatB))))
for(g in allgenes){
  theserows <- which(txptinfoA$Gene == g)
  countsmat_gene[g,colnames(countsmatA)] <- colSums(countsmatA[theserows,,drop = FALSE])
  theserows <- which(txptinfoB$Gene == g)
  countsmat_gene[g,colnames(countsmatB)] <- colSums(countsmatB[theserows,,drop = FALSE])
}

str(countsmat_gene)

#saveRDS(countsmat_gene, "results/stats/gene_counts_2022-06-27.rds")

countsmat_gene <- readRDS("results/stats/gene_counts_2022-06-27.rds")

# write to spreadsheet, include gene metadata
hasbb <- grepl("BEEBASE:", gtf1$Dbxref)
any(hasbb & !grepl("^BEEBASE:", gtf1$Dbxref)) # FALSE, BEEBASE id always first
beebaseID <- rep(NA_character_, length(gtf1))
beebaseID[hasbb] <- sub("^.*BEEBASE:", "", sub(",.*$", "", gtf1$Dbxref[hasbb]))

temp <- match(allgenes, gtf1$gene_id)
# write.table(cbind(Gene = allgenes, Symbol = gtf1$gene[temp], Beebase = beebaseID[temp],
#                   countsmat_gene),
#             row.names = FALSE, sep = "\t",
#             file = "results/stats/gene_counts_2022-06-27.txt")

# QC; how many genes detected per sample?
targets$Genes_detected <- colSums(countsmat_gene > 0)[targets$FlairName]
targets$Known_genes_detected <-
  colSums(countsmat_gene[grep(knowngene_pat, rownames(countsmat_gene)),] > 0)[targets$FlairName]

ggplot(targets, aes(x = Reads_aligned, y = Genes_detected, color = Group)) +
  geom_text(aes(label = Label))

ggplot(targets, aes(x = Reads_known_txpt_nuclear + Reads_known_gene_novel_isoform_nuclear,
                    y = Known_genes_detected, color = Group)) +
  geom_text(aes(label = Label))
# --> Linear relationships although AC samples are at the bottom

# DGEList ####

d <- DGEList(countsmat_gene[,targets$FlairName], samples = targets)

head(d$counts)
d$samples

d$genes <- data.frame(row.names = rownames(countsmat_gene),
                      Symbol = gtf1$gene[match(rownames(countsmat_gene), gtf1$gene_id)],
                      Beebase = beebaseID[match(rownames(countsmat_gene), gtf1$gene_id)])

head(d$genes)

# Quality control from standard RNA-seq ####
#check the highest and lowest counts
range(rowSums(d$counts))
# 1 622114

max(rowSums(d$counts)) / sum(d$counts)
# 0.1052908

hist(log(rowSums(d$counts)))
head(sort(rowSums(d$counts) / sum(d$counts), decreasing = TRUE))

top10 <- rownames(d$counts)[order(rowSums(d$counts), decreasing = TRUE)[1:10]]
d$genes[top10,,drop = FALSE]
# NC_001566.1 = mitochondrium, Melt = honeybee venom, CSP3 = chemosensory protein, Mrjp1 = royal jelly protein

barplot(apply(d$counts, 2, function(x) max(x) / sum(x)), las = 2, cex.names=0.7, 
        names.arg = d$samples$Label, main = "Proportion of highest reads") # somewhat variable by experimental group

rownames(d$counts)[apply(d$counts, 2, which.max)] # it's the nuclear-really-mitochondrial for all of A, mitochondrial for 4/6 B

# Remove mitochondrial genes, since these are impacting TMM factors ####
mt.filt <- txptinfoA$Mito[match(rownames(d), txptinfoA$Gene)] |
  txptinfoB$Mito[match(rownames(d), txptinfoB$Gene)]
mt.filt[is.na(mt.filt)] <- FALSE
sum(!mt.filt) # 19902 genes out of 19920

d$genes[mt.filt,]

d <- d[!mt.filt,, keep.lib.sizes = FALSE]

# Normalization and preliminary clustering ####
d <- calcNormFactors(d, method = "TMM")

barplot(d$samples$norm.factors, names.arg = d$samples$Label, 
        las=2, cex.names=0.8, main="TMM norm factors") # variable, ok considering batch effects
abline(h=1, lty=2)

logCPM <- cpm(d, log=TRUE)

mymds <- plotMDS(logCPM, top = 5000, labels = d$samples$Label)
# First axis is cross A vs. B.
# Second axis separates T3-Ag1:A2C from the rest, somewhat separates BT and BC

glMDSPlot(logCPM, top = 5000, labels = d$samples$Label, groups = d$samples$Group,
          html = "MDSclustering_preFiltering_2022-06-30", folder = "results/stats/glimma/")
# Third axis separates C vs. T

# Filtering ####

hist(rowMeans(logCPM))
summary(logCPM)
# Zeros got converted to 2.338; 2 ^ 2.338 = 5.056

# find a good CPM cutoff

#pdf("results/stats/CPM_cutoff_2022-06-30.pdf")
for(i in 4:12){
  hist(rowSums(logCPM > log2(i)), 100,
       main = i)
}
#dev.off()

# What cutoff would get us 10 reads
10 / 154000 * 1e6 # 64
10 / 400000 * 1e6 # 25

#How many have at least 8 cpm in at least 3 samples? (i.e. detected in at least 3 samples)

min_cpm <- 8
min_samp <- 3

hist(rowSums(logCPM > log2(min_cpm)), 100)
# mostly <3 or all

ggplot(mapping = aes(x = d$samples$Label, y = colSums(logCPM > log2(min_cpm)),
                     fill = d$samples$Group)) +
  geom_col() +
  labs(x = "Sample", y = "Number of genes", fill = "Group",
       title = paste0("Number of genes with >", min_cpm, " log counts per million"))
# reasonably even

i.filter <- rowSums(logCPM > log2(min_cpm)) >= min_samp
sum(i.filter)
# 14462
mean(i.filter) # 0.7266606

sum(grepl(knowngene_pat, rownames(logCPM)[i.filter])) # 8058 annotated genes

#cut down to filtered gene list

d.filt <- d[i.filter, , keep.lib.sizes=F]
dim(d.filt)
# 14462    12

range(d.filt$sample$lib.size / d$samples$lib.size)
#  0.9866894 0.9975208

#re-do TMM factors, just in case:

d.filt <- calcNormFactors(d.filt)
all.equal( d$samples$norm.factors, d.filt$samples$norm.factors)
# "Mean relative difference: 0.01286792"

ggplot(mapping = aes(x = d$samples$norm.factors, y = d.filt$samples$norm.factors,
                     col = d$samples$Group)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_text(aes(label = d$samples$Label)) +
  labs(x = "TMM norm factors before filtering",
       y = "TMM norm factors after filtering",
       col = "Group")
# same general order, follow line pretty well since removing MT reads

# get logCPM after filtering and normalization
logCPM.filt <- cpm(d.filt, log = T, prior.count =2)

plotDensities( logCPM.filt, group = d$samples$Group, col = 1:4, legend="topright" )
# A samples more multimodal than B samples
plotDensities( logCPM.filt, group = d$samples$Batch, col = 1:4, legend="topright" )
# batch effect for evenness/multimodality of distribution

glMDSPlot(logCPM.filt, top = 5000, labels = d.filt$samples$Label, groups = d.filt$samples$Group,
          html = "MDSclustering_postFiltering_2022-06-30", folder = "results/stats/glimma/")
# highly similary to prefiltering

mymds.filt <- plotMDS(logCPM.filt, top = 5000, labels = d$samples$Label, plot = FALSE)

plot(mymds.filt)

# Differential gene expression ####
model0 <- model.matrix(~ 0 + Group,
                       data = d.filt$samples)
colnames(model0) <- sub("Group", "", colnames(model0))

model0

contrasts0 <- makeContrasts(BT - BC,
                            AT - AC,
                            BC - AC,
                            (BT + AT) - (BC + AC),
                            (BT + BC) - (AT + AC),
                            (BT - BC) - (AT - AC),
                            levels = model0)

contrasts0

#put normalized values in EList object type
e <- new("EList", list(E =logCPM.filt,
                       genes=d.filt$genes, targets=d.filt$samples))

# fit model
fit0 <- lmFit(e, design = model0)

fit1 <- eBayes(fit0, trend = TRUE )

summary(decideTests(fit1))
#           AC    AT    BC    BT
# Down       0     0     0     0
# NotSig     3     0     0     0
# Up     14459 14462 14462 14462

fit2 <- contrasts.fit(fit0, contrasts = contrasts0) %>% eBayes(trend = TRUE)

summary(decideTests(fit2))
#        BT - BC AT - AC BC - AC (BT + AT) - (BC + AC) (BT + BC) - (AT + AC) (BT - BC) - (AT - AC)
# Down         0       0     133                     0                   306                     0
# NotSig   14461   14461   13922                 14462                 13465                 14462
# Up           1       1     407                     0                   691                     0

summary(decideTests(fit2, p.value = 0.1))
#        BT - BC AT - AC BC - AC (BT + AT) - (BC + AC) (BT + BC) - (AT + AC) (BT - BC) - (AT - AC)
# Down         0       0     211                     2                   448                     0
# NotSig   14461   14461   13580                 14448                 12910                 14462
# Up           1       1     671                    12                  1104                     0

# export results
res.dge <- data.frame(Gene = row.names(d.filt$genes),
                      d.filt$genes,
                      AveExpr = fit1$Amean)

for(coef in colnames(contrasts0)){
  temp <- topTable(fit2, coef = coef, number = Inf)
  colnames(temp) <- paste(colnames(temp), coef)
  res.dge <- cbind(res.dge, temp[res.dge$Gene, c(3, 6, 7)])
}

temp <- topTable(fit2, number = Inf) # F-statistics
colnames(temp) <- paste(colnames(temp), "ANOVA")
res.dge <- cbind(res.dge, temp[res.dge$Gene, c("P.Value ANOVA", "adj.P.Val ANOVA")])

res.dge$In_heatmap <- res.dge[["adj.P.Val ANOVA"]] < 0.05

# write.table(res.dge, file = "results/stats/differential_gene_expression_results_2022-06-30.txt",
#             sep = "\t", row.names = FALSE)

# Heatmap for DGE ####
siggenes <- rownames(topTable(fit2, number = Inf, p.value = 0.05)) # 619 genes

# not adjusted for tissue effect
heatdata.orig <- logCPM.filt[siggenes,] %>% t() %>% scale() %>% t()
colnames(heatdata.orig) <- d.filt$samples$Label

hmap(heatdata.orig, scale = "none", col = plasma(36),
     ColSideColors = colorkey[d.filt$samples$Group],
     main = "Unadjusted values")

# export adjusted and unadjusted logCPM values
# write.table(cbind(Gene = row.names(d.filt$genes),
#                   d.filt$genes,
#                   logCPM.filt),
#             "results/stats/filtered_logCPM_2022-06-30.txt",
#             sep = "\t", row.names = FALSE)

# Session info ####
sessionInfo()
# R version 4.2.0 (2022-04-22 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
# 
# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] kableExtra_1.3.4     heatmaply_1.3.0      plotly_4.10.0        viridis_0.6.2        viridisLite_0.4.0    seriation_1.3.5     
# [7] ggrepel_0.9.1        RColorBrewer_1.1-3   Glimma_2.6.0         tidyr_1.2.0          dplyr_1.0.9          ggplot2_3.3.6       
# [13] GenomicRanges_1.48.0 GenomeInfoDb_1.32.2  IRanges_2.30.0       S4Vectors_0.34.0     BiocGenerics_0.42.0  edgeR_3.38.1        
# [19] limma_3.52.2        
# 
# loaded via a namespace (and not attached):
# [1] colorspace_2.0-3            rjson_0.2.21                ellipsis_0.3.2              ggridges_0.5.3             
# [5] XVector_0.36.0              rstudioapi_0.13             farver_2.1.0                bit64_4.0.5                
# [9] AnnotationDbi_1.58.0        fansi_1.0.3                 xml2_1.3.3                  codetools_0.2-18           
# [13] splines_4.2.0               cachem_1.0.6                geneplotter_1.74.0          knitr_1.39                 
# [17] jsonlite_1.8.0              Rsamtools_2.12.0            annotate_1.74.0             png_0.1-7                  
# [21] pheatmap_1.0.12             compiler_4.2.0              httr_1.4.3                  assertthat_0.2.1           
# [25] Matrix_1.4-1                fastmap_1.1.0               lazyeval_0.2.2              cli_3.3.0                  
# [29] htmltools_0.5.2             tools_4.2.0                 gtable_0.3.0                glue_1.6.2                 
# [33] GenomeInfoDbData_1.2.8      reshape2_1.4.4              Rcpp_1.0.8.3                Biobase_2.56.0             
# [37] vctrs_0.4.1                 Biostrings_2.64.0           svglite_2.1.0               rtracklayer_1.56.1         
# [41] crosstalk_1.2.0             iterators_1.0.14            xfun_0.31                   stringr_1.4.0              
# [45] rvest_1.0.2                 lifecycle_1.0.1             restfulr_0.0.15             statmod_1.4.36             
# [49] XML_3.99-0.10               dendextend_1.15.2           zlibbioc_1.42.0             scales_1.2.0               
# [53] TSP_1.2-0                   MatrixGenerics_1.8.0        parallel_4.2.0              SummarizedExperiment_1.26.1
# [57] SingleCellExperiment_1.18.0 yaml_2.3.5                  memoise_2.0.1               gridExtra_2.3              
# [61] stringi_1.7.6               RSQLite_2.2.14              genefilter_1.78.0           BiocIO_1.6.0               
# [65] foreach_1.5.2               BiocParallel_1.30.3         systemfonts_1.0.4           rlang_1.0.2                
# [69] pkgconfig_2.0.3             matrixStats_0.62.0          bitops_1.0-7                evaluate_0.15              
# [73] lattice_0.20-45             purrr_0.3.4                 labeling_0.4.2              GenomicAlignments_1.32.0   
# [77] htmlwidgets_1.5.4           cowplot_1.1.1               bit_4.0.4                   tidyselect_1.1.2           
# [81] plyr_1.8.7                  magrittr_2.0.3              DESeq2_1.36.0               R6_2.5.1                   
# [85] generics_0.1.2              DelayedArray_0.22.0         DBI_1.1.3                   pillar_1.7.0               
# [89] withr_2.5.0                 survival_3.3-1              KEGGREST_1.36.2             RCurl_1.98-1.7             
# [93] tibble_3.1.7                crayon_1.5.1                utf8_1.2.2                  dittoSeq_1.8.1             
# [97] rmarkdown_2.14              locfit_1.5-9.5              grid_4.2.0                  data.table_1.14.2          
# [101] blob_1.2.3                  digest_0.6.29               webshot_0.5.3               xtable_1.8-4               
# [105] munsell_0.5.0               registry_0.5-1 
