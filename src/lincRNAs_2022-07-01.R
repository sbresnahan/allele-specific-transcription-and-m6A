# Identify novel lincRNAs from the Flair output

library(GenomicRanges)
library(Biostrings)
library(ORFhunteR)
library(ggplot2)
library(limma)
library(edgeR)
library(Glimma)
library(seriation)
library(viridis)
library(RColorBrewer)
library(magrittr)

# Load genome annotation ####
gene_info <- readRDS("results/final_report/annotations.rds")

table(gene_info$type)

gene_info[gene_info$type == "lnc_RNA"] # 3146 listed

# Transcript info from Flair ####

gtf <- rtracklayer::import("results/flair/2022-06-14/All_samples_collapse.isoforms.gtf")

# correct gene and txpt names
ncrows <- which(endsWith(gtf$transcript_id, "_NC")) # 25412 rows
gtf$gene_id[ncrows] <- paste0("NC_", gtf$gene_id[ncrows])
gtf$transcript_id[ncrows] <- sub("_NC$", "", gtf$transcript_id[ncrows])

nwrows <- which(endsWith(gtf$transcript_id, "_NW")) # 234 rows
gtf$gene_id[nwrows] <- paste0("NW_", gtf$gene_id[nwrows])
gtf$transcript_id[nwrows] <- sub("_NW$", "", gtf$transcript_id[nwrows])

table(gtf$type)
# transcript       exon 
#      31475     137845

# Start table of transcripts and gather info ####
gtf2 <- gtf[gtf$type == "transcript"]

txpts <- data.frame(Transcript = gtf2$transcript_id,
                    Gene = gtf2$gene_id,
                    Mitochondrial = seqnames(gtf2) == "NC_001566.1" | gtf2$gene_id == "NC_037641.1:1186000")

mean(txpts$Mitochondrial) # 0.2 %

annotated_lncRNAs <- unique(gene_info$transcript_id[gene_info$type == "lnc_RNA"])

txpts$Annotated_lncRNA <- sub("-[[:digit:]]$", "", txpts$Transcript) %in% annotated_lncRNAs

mean(txpts$Annotated_lncRNA) # 1.4%

protein_coding_genes <- gene_info[which(gene_info$type == "gene" &
                                    gene_info$gene_biotype == "protein_coding")] # 9935

other_rna_genes <- gene_info[which(gene_info$type == "gene" &
                                     gene_info$gene_biotype %in% c("guide_RNA", "miRNA", "rRNA",
                                                                   "snoRNA", "snRNA", "tRNA"))] # 530

txpts$Other_RNA <- txpts$Gene %in% other_rna_genes$gene_id
mean(txpts$Other_RNA) # 0.04%
# There are just a few of these, so see what they are
other_rna_genes[other_rna_genes$gene_id %in% txpts$Gene[txpts$Other_RNA]]
txpts[txpts$Other_RNA,]

table(other_rna_genes[other_rna_genes$gene_id %in% txpts$Gene[txpts$Other_RNA]]$gene_biotype)

txpts$Protein_coding <- txpts$Gene %in% protein_coding_genes$gene_id
mean(txpts$Protein_coding) # 56.8% 

gene_overlaps <- countOverlaps(gtf2, c(protein_coding_genes, other_rna_genes),
                              type = "any", ignore.strand = TRUE)

txpts$Intergenic <- gene_overlaps == 0
mean(txpts$Intergenic) # 12.4%

# sanity check
any(txpts$Intergenic & txpts$Protein_coding) # FALSE
# txpts[txpts$Intergenic & txpts$Protein_coding,] # two transcripts
# gtf2[gtf2$transcript_id %in% c("188c3fa9-6433-4082-bd83-ba965beb17be", "fb71bd1b-ec8b-4ab5-8101-bdd463991dc5")]
# protein_coding_genes[protein_coding_genes$gene_id %in% c("GeneID:113219001", "GeneID:100577180")]
any(txpts$Intergenic & txpts$Other_RNA) # FALSE

table(txpts$Annotated_lncRNA, txpts$Intergenic)
#       FALSE  TRUE
# FALSE 27281  3743
# TRUE    286   165

# Get transcript length from FASTA ####
myfasta <- readDNAStringSet("results/flair/2022-06-14/All_samples_collapse.isoforms.fa")

setdiff(names(myfasta), paste(txpts$Transcript, txpts$Gene, sep = "_")) # same names but in different order

txpts$Length <- width(myfasta)[match(paste(txpts$Transcript, txpts$Gene, sep = "_"), names(myfasta))]

linc_so_far <- !txpts$Mitochondrial & txpts$Intergenic & txpts$Length >= 200
sum(linc_so_far) # 3806 out of 31475

# Any genes with transcripts within this group and not
intersect(txpts$Gene[linc_so_far], txpts$Gene[!linc_so_far]) # yes, 232

# Predict ORF lengths
linc_so_far_names <- paste(txpts$Transcript[linc_so_far], txpts$Gene[linc_so_far], sep = "_")

myfasta[linc_so_far_names]
head(linc_so_far_names)

orfpos <- lapply(as.character(myfasta[linc_so_far_names]), findORFs)

# Just using most basic prediction, not the machine learning algorithm.

orf_over_100 <- sapply(orfpos, function(x) any(as.integer(x[,"length"]) > 100L))
mean(orf_over_100) # 64%

orf_length <- sapply(orfpos, function(x) max(as.integer(x[,"length"])))
orf_length[is.infinite(orf_length)] <- 0L
summary(orf_length)

txpts$ORF_length[linc_so_far] <- orf_length

sum(!txpts$Mitochondrial & txpts$Intergenic & txpts$Length >= 200 & txpts$ORF_length <= 100) # 1358

#saveRDS(txpts, file = "results/stats/lincRNA_search_2022-07-01.rds")

txpts <- readRDS("results/stats/lincRNA_search_2022-07-01.rds")

summary(txpts$ORF_length[txpts$Annotated_lncRNA]) # even the annotated lncRNAs have orfs using this method

# Export transcripts to analyze with other tools ####
# I'm not going to filter based on ORF length right now, since the coding potential
# calculator should do a much better job of determining if a transcript is coding.

myfasta_sub1 <- myfasta[linc_so_far_names]
names(myfasta_sub1) <- txpts$Transcript[linc_so_far]

writeXStringSet(myfasta_sub1,
                file = "results/lincRNAs/putative_lincRNA_2022-07-01.fasta")

# write in smaller chunks
for(i in seq_len(ceiling(sum(linc_so_far) / 500))){
  i1 <- 500 * (i - 1) + 1
  i2 <- min(c(500 * i, sum(linc_so_far)))
  writeXStringSet(myfasta_sub1[i1:i2],
                  file = paste0("results/lincRNAs/putative_lincRNA_2022-07-01_subset", i, ".fasta"))
}

# load back in
myfasta_sub1 <- readDNAStringSet("results/lincRNAs/putative_lincRNA_2022-07-01.fasta")

# External tools ####
# CPAT at http://lilab.research.bcm.edu/
# Select Fly (dm3, BDGP Release 5)
# --> Problems getting it to work this time

# CPC2 at http://cpc2.gao-lab.org/batch.php
# Select Fruitfly (dm6) (although this might only be for BED/GTF to retrieve sequences)
# --> Worked fine

# Import CPAT results
txpts$CPC2.prob <- rep(NA_real_, nrow(txpts))
txpts$CPC2.label <- rep(NA_character_, nrow(txpts))
cpc2tab <- read.table("results/lincRNAs/result_cpc2_all.txt", sep = "\t", header = FALSE)
theserows <- match(cpc2tab[[1]], txpts$Transcript)
txpts$CPC2.prob[theserows] <- cpc2tab[[6]]
txpts$CPC2.label[theserows] <- cpc2tab[[7]]

table(txpts$CPC2.label[linc_so_far])
#   coding noncoding 
#       61      3745

# BLASTX to UniProt + Swiss-Prot; see script blast_uniprot_sprot_2022-07-01.sh
blastres <- read.csv("results/lincRNAs/blast_uniprot_sprot_2022-07-01.csv",
                     header = FALSE)
colnames(blastres) <- c("qseqid", "qacc", "sseqid", "sacc", "qstart",
                        "qend", "sstart", "send", "evalue", "length", "pident")

nrow(blastres) # 41740 matches
length(unique(blastres$qseqid)) # 972 transcripts
hist(blastres$pident) # most around 30% match

ggplot(blastres, aes(x = pident, y = -log10(evalue), color = length)) +
  geom_point() +
  scale_color_viridis_c()

hist(-log10(blastres$evalue))

length(unique(blastres$qseqid[blastres$evalue < 1e-3])) # 120

txpts$Min_UPSP_evalue <- rep(NA_real_, nrow(txpts))
for(tr in unique(blastres$qseqid)){
  txpts$Min_UPSP_evalue[match(tr, txpts$Transcript)] <-
    min(blastres$evalue[blastres$qseqid == tr])
}

table(txpts$CPC2.label[linc_so_far], is.na(txpts$Min_UPSP_evalue[linc_so_far]))
#           FALSE TRUE
# coding       44   17
# noncoding   928 2817
# --> 61 had coding potential, and most of these aligned to swissprot/uniprot

table(txpts$CPC2.label[linc_so_far],
      !is.na(txpts$Min_UPSP_evalue[linc_so_far]) & txpts$Min_UPSP_evalue[linc_so_far] < 1e-3) # this cutoff used in ms
#           FALSE TRUE
# coding       37   24
# noncoding  3649   96

table(txpts$CPC2.label[linc_so_far],
      !is.na(txpts$Min_UPSP_evalue[linc_so_far]) & txpts$Min_UPSP_evalue[linc_so_far] < 1)
#           FALSE TRUE
# coding       29   32
# noncoding  3509  236

hist(-log10(txpts$Min_UPSP_evalue), breaks = 50)

table(txpts$CPC2.label[linc_so_far],
      !is.na(txpts$Min_UPSP_evalue[linc_so_far]) & txpts$Min_UPSP_evalue[linc_so_far] < 0.05)
#           FALSE TRUE
# coding       36   25
# noncoding  3630  115

# BLAST against other RNAs ####
# tRNAs from http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Amell4.5/apiMel4-mature-tRNAs.fa
# SILVA database on Biocluster at /home/mirrors/silva/138.1/; no FASTAs?
# download https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz
# and https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
# From NONCODE: http://noncode.org/datadownload/NONCODEv6_fruitfly.fa.gz

# Where else can we get snRNAs, snoRNAs, and 7SL/SRP?
# Are any of these annotated in the bee genome?
table(gene_info$type)

# Yes, actually miRNA, tRNA, snoRNA, snRNA, rRNA, and guide_RNA were excluded above based on their annotation (maybe still need SRP)
other_rna_genes # 530

# are they in the transcriptome?
txptome <- readDNAStringSet("data/reference/GCF_003254395.2_Amel_HAv3.1_rna+virus.fna")
head(names(txptome)) # doesn't have gene name

# get txpt IDs from GTF
other_rna_txpts <- gene_info[gene_info$type != "gene" & gene_info$gene_id %in% other_rna_genes$gene_id]
other_rna_txpts
sum(!is.na(other_rna_txpts$Name)) # 334, maybe not everything

# Extract the ncRNA sequences at the gene level
ncRNA_seqs <- Rsamtools::scanFa("data/reference/GCF_003254395.2_Amel_HAv3.1+virus_genomic.fna",
                     param = other_rna_genes)
names(ncRNA_seqs) <- other_rna_genes$Name
ncRNA_seqs

# writeXStringSet(ncRNA_seqs,
#                 file = "data/reference/GCF_003254395.2_Amel_HAv3.1_shortnoncoding_genes.fna")

# 7SL/SRP?
gene_info[grep("signal recognition particle", gene_info$product)] # protein coding
gene_info[grep("7SL", gene_info$product)] # nothing
# Get Drosophila versions from https://www.ncbi.nlm.nih.gov/nuccore/NR_001992.2?report=fasta
# and https://www.ncbi.nlm.nih.gov/nuccore/NR_037753.2?report=fasta
# put in file data/references/Drosophila_7SL_RNA.fa

# From the terminal:
# cd data/reference
# cat GCF_003254395.2_Amel_HAv3.1_shortnoncoding_genes.fna Drosophila_7SL_RNA.fa > noncoding_to_blast_2022-01-27.fa

# Read in BLAST output
blastres2 <- read.csv("results/lincRNAs/blast_ncRNA_2022-07-06.csv",
                     header = FALSE)
colnames(blastres2) <- c("qseqid", "qacc", "sseqid", "sacc", "qstart",
                        "qend", "sstart", "send", "evalue", "length", "pident")
# Two weak miRNA matches, and some other pretty strong matches to some annotated genes
myhits <- unique(blastres2$sseqid[blastres2$evalue < 1e-3]) # six genes

gene_info[which(gene_info$Name %in% myhits)] # one miRNA and six ribosomal RNA

nov_rRNA <- unique(blastres2$qseqid[blastres2$evalue < 1e-3])

gtf2[gtf2$transcript_id %in% nov_rRNA] # not overlapping annotated ones

txpts$Min_ncRNA_evalue <- rep(NA_real_, nrow(txpts))
for(tr in unique(blastres2$qseqid)){
  txpts$Min_ncRNA_evalue[match(tr, txpts$Transcript)] <-
    min(blastres2$evalue[blastres2$qseqid == tr])
}

#saveRDS(txpts, file = "results/stats/lincRNA_search_2022-07-06.rds")

txpts <- readRDS("results/stats/lincRNA_search_2022-07-06.rds")

# How many lincRNA txpts left? ####
readidpat <- "^[[:xdigit:]]{8}-([[:xdigit:]]{4}-){3}[[:xdigit:]]{12}(-[[:digit:]])?"
txptpat <- "^[NX][MR]_[[:digit:]]+\\.[[:digit:]]+(-[[:digit:]])?"
pseudogene_pat <- "^gene[[:digit:]]+(-[[:digit:]])?"

unknowngene_pat <- "N[CW]_[[:digit:]]+\\.[[:digit:]]+:[[:digit:]]+$"
knowngene_pat <- "GeneID:[[:digit:]]+$"
viralgene_pat <- "(DW|ABP|IAP|SB|CBP|VD)V(:[[:digit::]+)?$"
rRNA_pat <- "l-rRNA$"

linc_final <- !(!linc_so_far | (!is.na(txpts$CPC2.label) & txpts$CPC2.label == "coding") |
                  (!is.na(txpts$Min_UPSP_evalue) & txpts$Min_UPSP_evalue < 1e-3) |
                  (!is.na(txpts$Min_ncRNA_evalue) & txpts$Min_ncRNA_evalue < 1e-3) |
                  grepl(viralgene_pat, txpts$Gene))
sum(linc_final) # 3642 transcripts

length(unique(txpts$Gene[linc_final])) # 3819 genes

# Bookkeeping on origins of these
txpts[!grepl(readidpat, txpts$Transcript) & linc_final,]
txpts[!grepl(readidpat, txpts$Transcript) & linc_final & !txpts$Annotated_lncRNA,] # just one pseudogene

txpts[!grepl(knowngene_pat, txpts$Gene) &
        !grepl(unknowngene_pat, txpts$Gene) &
        linc_final,] # no rows once viral removed

# Import read counts for these transcripts ####
mycounts <- read.table("results/flair/2022-06-14/counts_matrix_AvsB_2022-07-06.tsv",
                       sep = "\t", header = TRUE)
colnames(mycounts)[-1] <- sub("_.*$", "", colnames(mycounts)[-1])

# how many of my lincRNAs are in the counts output (a small num of transcripts from the GTF were excluded from counts)
txpts$FlairName <- paste(txpts$Transcript, txpts$Gene, sep = "_")
mean(txpts$FlairName[linc_final] %in% mycounts$ids) # 99.9%
setdiff(txpts$FlairName[linc_final], mycounts$ids) # two novel transcripts

# Total sample read depth
samDepth <- colSums(mycounts[,-1])

# Mitochondrial depth
mitoDepth <- colSums(mycounts[mycounts$id %in% txpts$FlairName[txpts$Mitochondrial],-1])

# Nuclear depth
nucDepth <- colSums(mycounts[mycounts$id %in%
                               txpts$FlairName[!txpts$Mitochondrial &
                                                 !grepl(viralgene_pat, txpts$Gene)],-1])

# subset to lincRNAs
mycounts <- mycounts[mycounts$id %in% txpts$FlairName[linc_final],] # 3640 transcripts

# Export raw counts and metadata
temp <- match(mycounts$ids, txpts$FlairName)
out <- data.frame(ID = mycounts$ids,
                  Transcript = txpts$Transcript[temp],
                  Gene = txpts$Gene[temp],
                  Length = txpts$Length[temp])
identical(temp, match(out$Transcript, gtf2$transcript_id)) # TRUE
out$Chromosome <- as.character(seqnames(gtf2[temp]))
out$Start <- start(gtf2[temp])
out$End <- end(gtf2[temp])
out$Strand <- as.character(strand(gtf2[temp]))
out$Exons <-
  sapply(out$Transcript,
         function(x){
           ex <- gtf[gtf$type == "exon" & gtf$transcript_id == x]
           ex_pos <- paste(start(ex), end(ex), sep = "-")
           paste(ex_pos, collapse = ";")
         })
out$Annotated_lncRNA <- txpts$Annotated_lncRNA[temp]
out$Annotated_pseudogene <- grepl(pseudogene_pat, out$Transcript)
out <- cbind(out, mycounts[,-1])

# write.table(out, sep = "\t", row.names = FALSE,
#             file = "results/lincRNA_results/raw_counts_txpt_info_2022-07-06.txt")

# writeXStringSet(myfasta[out$ID],
#                 file = "results/lincRNA_results/lincRNA_sequence_2022-07-06.fa")

table(out$Strand)
#    -    + 
# 1762 1878

# quick test on improvement in txpt detection
temp <- read.table("results/lincRNA_results/raw_counts_txpt_info_2022-01-27.txt", sep = "\t", header = TRUE)
table(temp$Strand)
#    -    + 
# 1039 2004

# Prep for DGE and QC ####
# Depth just for lincRNAs
lincRNAdepth <- colSums(mycounts[,-1])

# Set up sample metadata table
targets <- data.frame(Sample = colnames(mycounts)[-1])
targets$Cross <- substring(targets$Sample, 8, 8)
targets$Aggression <- substring(targets$Sample, 10, 10)
targets$Group <- paste0(targets$Cross, targets$Aggression)
targets$Label <- factor(targets$Sample,
                        levels = targets$Sample[order(targets$Cross, targets$Aggression)])
targets$Total_depth <- samDepth
targets$Mitochondrial_depth <- mitoDepth
targets$Nuclear_depth <- nucDepth
targets$lincRNA_depth <- lincRNAdepth

# Proportion lincRNA related to aggression?
summary(lm(lincRNA_depth / Nuclear_depth ~ Aggression, data = targets)) # P = 0.0545

# DGEList ####
countsmat <- as.matrix(mycounts[,-1])
rownames(countsmat) <- mycounts$id
d <- DGEList(countsmat, samples = targets)
temp <- match(mycounts$ids, txpts$FlairName)
d$genes <- data.frame(row.names = mycounts$ids,
                      txpts[temp,c("Transcript", "Gene")])

# Normalization and preliminary clustering ####
d <- calcNormFactors(d, method = "TMM")
targets$TMM_prefilt <- d$samples$norm.factors
logCPM <- cpm(d, log=TRUE)

mymds <- plotMDS(logCPM, top = 5000, labels = d$samples$Label)
# First axis separates A vs B, second C vs T to some extent

glMDSPlot(logCPM, top = 5000, labels = d$samples$Label, groups = d$samples$Group,
          html = "MDSclustering_lincRNAs_preFiltering_2022-07-06", folder = "results/stats/glimma/")

# Filtering ####

hist(rowMeans(logCPM))
summary(logCPM)
# Zeros got converted to 7.360; 2 ^ 7.360 = 164

hist(rowSums(logCPM > log2(200)), 100, main = 200)
hist(rowSums(logCPM > log2(300)), 100, main = 300)
hist(rowSums(logCPM > log2(350)), 100, main = 350)
hist(rowSums(logCPM > log2(400)), 100, main = 400)
hist(rowSums(logCPM > log2(500)), 100, main = 500)

#How many have at least 250 cpm in at least 3 samples? (i.e. detected in at least 3 samples)

min_cpm <- 250
min_samp <- 3

hist(rowSums(logCPM > log2(min_cpm)), 100)

ggplot(mapping = aes(x = d$samples$Label, y = colSums(logCPM > log2(min_cpm)),
                     fill = d$samples$Group)) +
  geom_col() +
  labs(x = "Sample", y = "Number of genes", fill = "Group",
       title = paste0("Number of genes with >", min_cpm, " log counts per million"))

i.filter <- rowSums(logCPM > log2(min_cpm)) >= min_samp
sum(i.filter)
# 2415
mean(i.filter) # 0.6634615

sum(grepl(knowngene_pat, rownames(logCPM)[i.filter])) # 192 txpts from known genes

d.filt <- d[i.filter, , keep.lib.sizes=F]
dim(d.filt)
# 2415   12

range(d.filt$sample$lib.size / d$samples$lib.size)
# 0.9579246 0.9768673

#re-do TMM factors, just in case:

d.filt <- calcNormFactors(d.filt)
all.equal( d$samples$norm.factors, d.filt$samples$norm.factors)
# "Mean relative difference: 0.05000326"

ggplot(mapping = aes(x = d$samples$norm.factors, y = d.filt$samples$norm.factors,
                     col = d$samples$Group)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_text(aes(label = d$samples$Label)) +
  labs(x = "TMM norm factors before filtering",
       y = "TMM norm factors after filtering",
       col = "Group")
# same general pattern

# get logCPM after filtering and normalization
logCPM.filt <- cpm(d.filt, log = T, prior.count =2)

plotDensities( logCPM.filt, group = d$samples$Group, col = 1:4, legend="topright" )
# Basically relates to library size and missing data

glMDSPlot(logCPM.filt, top = 5000, labels = d.filt$samples$Label, groups = d.filt$samples$Group,
          html = "MDSclustering_lincRNAs_postFiltering_2022-07-06", folder = "results/stats/glimma/")
# similar to pre-filtering

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
                            BT - AT,
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

fit2 <- contrasts.fit(fit0, contrasts = contrasts0) %>% eBayes(trend = TRUE)

summary(decideTests(fit2))
#        BT - BC AT - AC BC - AC BT - AT (BT + AT) - (BC + AC) (BT + BC) - (AT + AC) (BT - BC) - (AT - AC)
# Down         1       0       1       4                     3                     2                     0
# NotSig    2409    2413    2404    2408                  2401                  2396                  2415
# Up           5       2      10       3                    11                    17                     0

# export results
res.dge <- data.frame(ID = row.names(d.filt$genes),
                      d.filt$genes,
                      AveExpr = fit2$Amean)

for(coef in colnames(contrasts0)){
  temp <- topTable(fit2, coef = coef, number = Inf)
  colnames(temp) <- paste(colnames(temp), coef)
  res.dge <- cbind(res.dge, temp[res.dge$ID, c(3, 6, 7)])
}

temp <- topTable(fit2, number = Inf) # F-statistics
colnames(temp) <- paste(colnames(temp), "ANOVA")
res.dge <- cbind(res.dge, temp[res.dge$ID, c("P.Value ANOVA", "adj.P.Val ANOVA")])

res.dge$In_heatmap <- res.dge[["adj.P.Val ANOVA"]] < 0.05

# write.table(res.dge, file = "results/lincRNA_results/differential_lincRNA_expression_results_2022-07-06.txt",
#             sep = "\t", row.names = FALSE)

# Heatmap for DGE ####
colorkey <- brewer.pal(4, "Paired")
names(colorkey) <- unique(targets$Group)

siggenes <- rownames(topTable(fit2, number = Inf, p.value = 0.05)) # 41 genes

heatdata <- logCPM.filt[siggenes,] %>% t() %>% scale() %>% t()
colnames(heatdata) <- d.filt$samples$Label

hmap(heatdata, scale = "none", col = plasma(36),
     ColSideColors = colorkey[d.filt$samples$Group],
     main = "lincRNAs significant in ANOVA at FDR < 0.05")
