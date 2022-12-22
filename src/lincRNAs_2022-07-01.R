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

# Glimma MA for DGE ####
for(coef in colnames(contrasts0)){
  glimmaMA(fit2, d.filt,
           logCPM.filt,
           groups = d.filt$samples$Group,
           coef = coef, main = coef,
           transform.counts = "none",
           html = paste0("results/stats/glimma/glimma_MA_lncRNA_", coef, "_2022-07-06.html"))
}

# Cross-reference to parent-specific data 1/31/22 ####
countsmatA <- as.matrix(read.table("results/flair/2022-06-14/counts_matrix_parenttxpt_A.tsv",
                                   row.names = 1, header = TRUE, check.names = FALSE))

countsmatB <- as.matrix(read.table("results/flair/2022-06-14/counts_matrix_parenttxpt_B.tsv",
                                   row.names = 1, header = TRUE, check.names = FALSE))

mode(countsmatA) <- "integer"
mode(countsmatB) <- "integer"

txptsA <- data.frame(Parent = sub("-.*$", "", rownames(countsmatA)),
                     FlairName = sub("(Both|SRR1465418[6-9])-", "", rownames(countsmatA)),
                     TotCounts = rowSums(countsmatA))

txptsB <- data.frame(Parent = sub("-.*$", "", rownames(countsmatB)),
                     FlairName = sub("(Both|SRR1465418[6-9])-", "", rownames(countsmatB)),
                     TotCounts = rowSums(countsmatB))

txpts$Parents_distinguished_A <- txptsA$Parent[match(txpts$FlairName, txptsA$FlairName)] != "Both"
txpts$Parents_distinguished_A[is.na(txpts$Parents_distinguished_A)] <- FALSE
mean(txpts$Parents_distinguished_A) # 63.1% of transcripts overall

txpts$Parents_distinguished_B <- txptsB$Parent[match(txpts$FlairName, txptsB$FlairName)] != "Both"
txpts$Parents_distinguished_B[is.na(txpts$Parents_distinguished_B)] <- FALSE
mean(txpts$Parents_distinguished_B) # 69.8% of transcripts overall

txpts$Detected_A <- txpts$FlairName %in% txptsA$FlairName
txpts$Detected_B <- txpts$FlairName %in% txptsB$FlairName

mean(txpts$Parents_distinguished_A & txpts$Parents_distinguished_B) # 50.3% of transcripts overall
mean(txpts$Parents_distinguished_A[linc_final] & txpts$Parents_distinguished_B[linc_final]) # 51.4% of lincRNAs

temp <- txptsA[txptsA$Parent == "SRR14654188",]
txpts$Counts_AHB_F_SRR14654188 <- temp$TotCounts[match(txpts$FlairName, temp$FlairName)]
txpts$Counts_AHB_F_SRR14654188[txpts$Parents_distinguished_A & is.na(txpts$Counts_AHB_F_SRR14654188)] <- 0
temp <- txptsA[txptsA$Parent == "SRR14654189",]
txpts$Counts_EHB_M_SRR14654189 <- temp$TotCounts[match(txpts$FlairName, temp$FlairName)]
txpts$Counts_EHB_M_SRR14654189[txpts$Parents_distinguished_A & is.na(txpts$Counts_EHB_M_SRR14654189)] <- 0
temp <- txptsB[txptsB$Parent == "SRR14654186",]
txpts$Counts_EHB_F_SRR14654186 <- temp$TotCounts[match(txpts$FlairName, temp$FlairName)]
txpts$Counts_EHB_F_SRR14654186[txpts$Parents_distinguished_B & is.na(txpts$Counts_EHB_F_SRR14654186)] <- 0
temp <- txptsB[txptsB$Parent == "SRR14654187",]
txpts$Counts_AHB_M_SRR14654187 <- temp$TotCounts[match(txpts$FlairName, temp$FlairName)]
txpts$Counts_AHB_M_SRR14654187[txpts$Parents_distinguished_B & is.na(txpts$Counts_AHB_M_SRR14654187)] <- 0

View(txpts[linc_final,])

# Read in parent-specific DGE results
res.dge.gxp <- read.table("results/stats/differential_gene_expression_results_geneXparent_2022-06-30.txt",
                          sep = "\t", header = TRUE, check.names = FALSE)
temp <- match(txpts$Gene, res.dge.gxp$Gene)
txpts$Sig_M_vs_F <- res.dge.gxp$`adj.P.Val (ACM + ATM + BCM + BTM) - (ACF + ATF + BCF + BTF)`[temp] < 0.05
txpts$Sig_AHB_vs_EHB <- res.dge.gxp$`adj.P.Val (BCM + BTM - ACM - ATM) - (BCF + BTF - ACF - ATF)`[temp] < 0.05

# Do parent-specific differential expression analysis for lincRNAs ####
parspec_lincRNAs <- txpts$FlairName[(linc_final & txpts$Parents_distinguished_A & txpts$Parents_distinguished_B) |
                                      (linc_final & txpts$Parents_distinguished_A & !txpts$Detected_B) |
                                      (linc_final & txpts$Parents_distinguished_B & !txpts$Detected_A)] # 2266

countsmat_lincXparent <- matrix(0L, nrow = length(parspec_lincRNAs), ncol = (ncol(countsmatA) + ncol(countsmatB)) * 2,
                                dimnames = list(parspec_lincRNAs,
                                                c(paste0(colnames(countsmatA), "_F"),
                                                  paste0(colnames(countsmatA), "_M"),
                                                  paste0(colnames(countsmatB), "_F"),
                                                  paste0(colnames(countsmatB), "_M"))))

SRR14654188_rows <- txptsA$Parent == "SRR14654188"
SRR14654189_rows <- txptsA$Parent == "SRR14654189"
SRR14654186_rows <- txptsB$Parent == "SRR14654186"
SRR14654187_rows <- txptsB$Parent == "SRR14654187"
for(g in parspec_lincRNAs){
  theserows <- txptsA$FlairName == g
  countsmat_lincXparent[g,paste0(colnames(countsmatA), "_F")] <-
    colSums(countsmatA[theserows & SRR14654188_rows,, drop = FALSE])
  countsmat_lincXparent[g,paste0(colnames(countsmatA), "_M")] <-
    colSums(countsmatA[theserows & SRR14654189_rows,, drop = FALSE])
  
  theserows <- txptsB$FlairName == g
  countsmat_lincXparent[g,paste0(colnames(countsmatB), "_F")] <-
    colSums(countsmatB[theserows & SRR14654186_rows,, drop = FALSE])
  countsmat_lincXparent[g,paste0(colnames(countsmatB), "_M")] <-
    colSums(countsmatB[theserows & SRR14654187_rows,, drop = FALSE])
}

# write.table(cbind(ID = parspec_lincRNAs,
#                   txpts[match(parspec_lincRNAs, txpts$FlairName), c("Transcript", "Gene")],
#                   countsmat_lincXparent),
#             sep = "\t", row.names = FALSE,
#             file = "results/lincRNA_results/raw_counts_lincRNAxParent_2022-07-07.txt")

# set up targets
targets.lxp <- data.frame(Sample = colnames(countsmat_lincXparent))
targets.lxp$Cross <- substring(targets.lxp$Sample, 8, 8) 
targets.lxp$Aggression <- substring(targets.lxp$Sample, 10, 10) 
targets.lxp$Parent <- substring(targets.lxp$Sample, 18, 18)
targets.lxp$Group <- paste0(targets.lxp$Cross, targets.lxp$Aggression, targets.lxp$Parent)
targets.lxp$Label <- sub("_[CT]_(Dec|Apr|Nov|Jul)", "", targets.lxp$Sample)
targets.lxp$Label <- factor(targets.lxp$Label,
                            levels = targets.lxp$Label[order(targets.lxp$Cross, targets.lxp$Parent, targets.lxp$Aggression)])
rownames(targets.lxp) <- targets.lxp$Sample

colorkey2 <- brewer.pal(8, "Paired")
names(colorkey2) <- c('ACM', 'ATM', 'BCM', 'BTM', 'ACF', 'ATF', 'BCF', 'BTF')

# stack the counts so we can still normalize properly within sample
temp <- countsmat_lincXparent[,rownames(targets.lxp)]
tempF <- temp[,endsWith(colnames(temp), "F")]
tempM <- temp[,endsWith(colnames(temp), "M")]
colnames(tempF) <- sub("_F$", "", colnames(tempF))
colnames(tempM) <- sub("_M$", "", colnames(tempM))
rownames(tempF) <- paste0(rownames(tempF), "_F")
rownames(tempM) <- paste0(rownames(tempM), "_M")

genesF <- genesM <- d$genes[parspec_lincRNAs,]
rownames(genesF) <- paste0(rownames(genesF), "_F")
rownames(genesM) <- paste0(rownames(genesM), "_M")

d.lxp <- DGEList(rbind(tempF, tempM),
                 samples = data.frame(rownames = colnames(tempF),
                                      Sample = colnames(tempF),
                                      Label = sub("_[CT]_(Dec|Apr|Nov|Jul)", "", colnames(tempF))),
                 genes = rbind(genesF, genesM))

# Normalization and preliminary clustering for gene X parent ####
d.lxp <- calcNormFactors(d.lxp, method = "TMM")

barplot(d.lxp$samples$norm.factors, names.arg = d.lxp$samples$Label, 
        las=2, cex.names=0.8, main="TMM norm factors") # drops off a bit for B
abline(h=1, lty=2)

logCPM.lxp.stacked <- cpm(d.lxp, log=TRUE)
tempF <- logCPM.lxp.stacked[endsWith(rownames(logCPM.lxp.stacked), "_F"),]
tempM <- logCPM.lxp.stacked[endsWith(rownames(logCPM.lxp.stacked), "_M"),]
colnames(tempF) <- paste0(colnames(tempF), "_F")
colnames(tempM) <- paste0(colnames(tempM), "_M")
rownames(tempF) <- sub("_F$", "", rownames(tempF))
rownames(tempM) <- sub("_M$", "", rownames(tempM))
str(tempF)
str(tempM)
logCPM.lxp <- cbind(tempF, tempM)

mymds.lxp <- plotMDS(logCPM.lxp, top = 5000, labels = targets.lxp$Label) # clear grouping by cross and parent

glMDSPlot(logCPM.lxp, top = 5000, labels = d.lxp$samples$Label, groups = targets.lxp$Group,
          html = "MDSclustering_lincRNAXparent_preFiltering_2022-07-07", folder = "results/stats/glimma/")

# Filtering for lincRNA X parent####

hist(rowMeans(logCPM.lxp))
summary(logCPM.lxp)
# Zeros got converted to 8.539; 2 ^ 8.539 = 372

hist(rowSums(logCPM.lxp > 9), 100, main = 2 ^ 9)
hist(rowSums(logCPM.lxp > log2(500)), 100, main = 500)
hist(rowSums(logCPM.lxp > log2(600)), 100, main = 600)
hist(rowSums(logCPM.lxp > log2(650)), 100, main = 650)
hist(rowSums(logCPM.lxp > log2(700)), 100, main = 700)
hist(rowSums(logCPM.lxp > 9.5), 100, main = 2 ^ 9.5)
hist(rowSums(logCPM.lxp > 10), 100, main = 2 ^ 10)

min_cpm.lxp <- 600
min_samp <- 3

hist(rowSums(logCPM.lxp > log2(min_cpm.lxp)), 100)

# d.lxp$samples$Group <- paste0(substring(rownames(d.lxp$samples), 8, 8),
#                               substring(rownames(d.lxp$samples), 10, 10))

ggplot(mapping = aes(x = targets.lxp$Label, y = colSums(logCPM.lxp > log2(min_cpm.lxp)),
                     fill = targets.lxp$Group)) +
  geom_col() +
  labs(x = "Sample", y = "Number of genes", fill = "Group",
       title = paste0("Number of lincRNAs with >", min_cpm.lxp, " log counts per million")) +
  scale_fill_manual(values = colorkey2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

i.filter.lxp <- rowSums(logCPM.lxp > log2(min_cpm.lxp)) >= min_samp
sum(i.filter.lxp)
# 1216
mean(i.filter.lxp) # 53.7%

#cut down to filtered gene list

d.lxp.filt <- d.lxp[rep(i.filter.lxp, times = 2), , keep.lib.sizes=F]
dim(d.lxp.filt)
# 2432   12

range(d.lxp.filt$sample$lib.size / d.lxp$samples$lib.size)
#  0.9280395 0.9566940

#re-do TMM factors, just in case:

d.lxp.filt <- calcNormFactors(d.lxp.filt)
all.equal( d.lxp$samples$norm.factors, d.lxp.filt$samples$norm.factors)
# "Mean relative difference: 0.04942464"

ggplot(mapping = aes(x = d.lxp$samples$norm.factors, y = d.lxp.filt$samples$norm.factors,
                     col = d.lxp$samples$Group)) +
  geom_abline(slope = 1, intercept = 0) +
  geom_text(aes(label = d.lxp$samples$Label)) +
  labs(x = "TMM norm factors before filtering",
       y = "TMM norm factors after filtering",
       col = "Group") +
  scale_color_manual(values = colorkey)
# follows fairly closely

# get logCPM after filtering and normalization
logCPM.lxp.filt.stacked <- cpm(d.lxp.filt, log = T, prior.count =2)
tempF <- logCPM.lxp.filt.stacked[endsWith(rownames(logCPM.lxp.filt.stacked), "_F"),]
tempM <- logCPM.lxp.filt.stacked[endsWith(rownames(logCPM.lxp.filt.stacked), "_M"),]
colnames(tempF) <- paste0(colnames(tempF), "_F")
colnames(tempM) <- paste0(colnames(tempM), "_M")
rownames(tempF) <- sub("_F$", "", rownames(tempF))
rownames(tempM) <- sub("_M$", "", rownames(tempM))
str(tempF)
str(tempM)
logCPM.lxp.filt <- cbind(tempF, tempM)

# write.table(cbind(ID = rownames(logCPM.lxp.filt),
#                   txpts[match(rownames(logCPM.lxp.filt), txpts$FlairName),
#                         c("Transcript", "Gene")],
#                   logCPM.lxp.filt),
# sep = "\t", row.names = FALSE, file = "results/lincRNA_results/filtered_logCPM_lincRNAxParent_2022-07-07.txt")

glMDSPlot(logCPM.lxp.filt, top = 5000, labels = targets.lxp$Label, groups = targets.lxp$Group,
          html = "MDSclustering_lincRNAXparent_postFiltering_2022-07-07", folder = "results/stats/glimma/")

mymds.lxp.filt <- plotMDS(logCPM.lxp.filt, top = 5000, labels = targets.lxp$Label, plot = FALSE)

ggplot(mapping = aes(x = mymds.lxp.filt$x, y = mymds.lxp.filt$y,
                     fill = targets.lxp$Group)) +
  geom_label(aes(label = targets.lxp$Label)) +
  labs(fill = "Group",
       x = paste0("Dim 1 (", round(mymds.lxp.filt$var.exp[1] * 100, 1),"%)"),
       y = paste0("Dim 2 (", round(mymds.lxp.filt$var.exp[2] * 100, 1),"%)"),
       title = "MDS after filtering") +
  scale_fill_manual(values = colorkey2)

# Differential lincRNA expression across parent alleles ####
model0.lxp <- model.matrix(~ 0 + Group,
                           data = targets.lxp)
colnames(model0.lxp) <- sub("Group", "", colnames(model0.lxp))

model0.lxp

contrasts0.lxp <- makeContrasts((ACM + ATM + BCM + BTM) - (ACF + ATF + BCF + BTF), # parent effects
                                (BCM + BCF + BTM + BTF) - (ACM + ACF + ATM + ATF), # cross effects
                                (BCM + BTM - ACM - ATM) - (BCF + BTF - ACF - ATF), # interaction
                                # contrasts to detect effects that depend on MvF and AHBvEHB
                                (ACM + ATM) - (BCM + BTM + ACF + ATF + BCF + BTF) / 3,
                                (BCM + BTM) - (ACM + ATM + ACF + ATF + BCF + BTF) / 3,
                                (ACF + ATF) - (BCM + BTM + ACM + ATM + BCF + BTF) / 3,
                                (BCF + BTF) - (ACM + ATM + ACF + ATF + BCM + BTM) / 3,
                                # allelic effects on aggression
                                (BTM - BCM) - (BTF - BCF),
                                (ATM - ACM) - (ATF - ACF),
                                # Parent effects within crosses, to look at with interaction
                                (BTM + BCM) - (BTF + BCF),
                                (ATM + ACM) - (ATF + ACF),
                                levels = model0.lxp)
contrasts0.lxp

#put normalized values in EList object type
all(rownames(logCPM.lxp.filt) %in% rownames(d$genes)) # FALSE
tempgenes <- txpts[match(rownames(logCPM.lxp.filt), txpts$FlairName),c("Transcript", "Gene")]
rownames(tempgenes) <- rownames(logCPM.lxp.filt)
e.lxp <- new("EList", list(E =logCPM.lxp.filt,
                           genes=tempgenes, targets=targets.lxp))

# Set up for sample as a random effect
targets.lxp$Name <- sub("_[FM]$", "", targets.lxp$Label)
samplecorr <- duplicateCorrelation(e.lxp, design = model0.lxp,
                                   block = targets.lxp$Name)

# fit model
fit0.lxp <- lmFit(e.lxp, design = model0.lxp, block = targets.lxp$Name,
                  correlation = samplecorr$consensus)

fit2.lxp <- contrasts.fit(fit0.lxp, contrasts = contrasts0.lxp) %>% eBayes(trend = TRUE)

summary(decideTests(fit2.lxp))
#        (ACM + ATM + BCM + BTM) - (ACF + ATF + BCF + BTF) (BCM + BCF + BTM + BTF) - (ACM + ACF + ATM + ATF)
# Down                                                  15                                               111
# NotSig                                              1158                                              1057
# Up                                                    43                                                48
#        (BCM + BTM - ACM - ATM) - (BCF + BTF - ACF - ATF) (ACM + ATM) - (BCM + BTM + ACF + ATF + BCF + BTF)/3
# Down                                                 164                                                  22
# NotSig                                              1031                                                 927
# Up                                                    21                                                 267
#        (BCM + BTM) - (ACM + ATM + ACF + ATF + BCF + BTF)/3 (ACF + ATF) - (BCM + BTM + ACM + ATM + BCF + BTF)/3
# Down                                                    83                                                  43
# NotSig                                                1087                                                1154
# Up                                                      46                                                  19
#        (BCF + BTF) - (ACM + ATM + ACF + ATF + BCM + BTM)/3 (BTM - BCM) - (BTF - BCF) (ATM - ACM) - (ATF - ACF)
# Down                                                    26                         0                         0
# NotSig                                                1156                      1216                      1216
# Up                                                      34                         0                         0
#        (BTM + BCM) - (BTF + BCF) (ATM + ACM) - (ATF + ACF)
# Down                          44                        12
# NotSig                      1142                      1016
# Up                            30                       188

# export DGE results
res.dge.lxp <- data.frame(ID = row.names(fit0.lxp$genes),
                          fit0.lxp$genes,
                          AveExpr = fit0.lxp$Amean)

for(coef in colnames(contrasts0.lxp)){
  temp <- topTable(fit2.lxp, coef = coef, number = Inf)
  colnames(temp) <- paste(colnames(temp), coef)
  res.dge.lxp <- cbind(res.dge.lxp, temp[res.dge.lxp$ID, c(3, 6, 7)])
}

temp <- topTable(fit2.lxp, number = Inf) # F-statistics
colnames(temp) <- paste(colnames(temp), "ANOVA")
res.dge.lxp <- cbind(res.dge.lxp, temp[res.dge.lxp$ID, c("P.Value ANOVA", "adj.P.Val ANOVA")])

res.dge.lxp$In_heatmap <- res.dge.lxp[["adj.P.Val ANOVA"]] < 0.05 # 285 lincRNAs

res.dge.lxp[is.na(res.dge.lxp[["adj.P.Val ANOVA"]]),] # should be no rows

# write.table(res.dge.lxp, file = "results/lincRNA_results/differential_lincRNA_expression_results_lincRNAXparent_2022-07-07.txt",
#             sep = "\t", row.names = FALSE)

# Heatmap for gene x parent ####
siggenes.lxp <- rownames(topTable(fit2.lxp, number = Inf, p.value = 0.05))  # 285 lincRNAs

heatdata.lxp <- logCPM.lxp.filt[siggenes.lxp,] %>% t() %>% scale() %>% t()
colnames(heatdata.lxp) <- targets.lxp$Label

# tiff("results/lincRNA_results/heatmap_lincRNAxParent_2022-07-07.tiff",
#      res = 300, height = 9 * 300, width = 7 * 300,
#      compression = "lzw")
hmap(heatdata.lxp, scale = "none", col = plasma(36),
     ColSideColors = colorkey2[targets.lxp$Group],
     main = "285 lincRNAs significant in ANOVA at FDR < 0.05")
# dev.off()

# Glimma MA plots for lincRNA x parent ####

d.lxp.dummy <- DGEList(countsmat_lincXparent[rownames(logCPM.lxp.filt),rownames(targets.lxp)],
                       samples = targets.lxp,
                       genes = d$genes[rownames(logCPM.lxp.filt),])

glimmaMA(fit2.lxp, d.lxp.dummy,
         counts = logCPM.lxp.filt,
         groups = targets.lxp$Group,
         coef = "(ACM + ATM + BCM + BTM) - (ACF + ATF + BCF + BTF)",
         main = "Male vs. female allele",
         transform.counts = "none",
         html = paste0("results/stats/glimma/glimma_MA_lincRNAxParent_parent-effect_unadjusted_2022-07-07.html"))

glimmaMA(fit2.lxp, d.lxp.dummy,
         logCPM.lxp.filt,
         groups = targets.lxp$Group,
         coef = "(BCM + BCF + BTM + BTF) - (ACM + ACF + ATM + ATF)",
         main = "Cross B vs. Cross A",
         transform.counts = "none",
         html = paste0("results/stats/glimma/glimma_MA_lincRNAxParent_cross-effect_unadjusted_2022-07-07.html"))

glimmaMA(fit2.lxp, d.lxp.dummy,
         logCPM.lxp.filt,
         groups = targets.lxp$Group,
         coef = "(BCM + BTM - ACM - ATM) - (BCF + BTF - ACF - ATF)",
         main = "AHB vs. EHB (interaction of allele and cross)",
         transform.counts = "none",
         html = paste0("results/stats/glimma/glimma_MA_lincRNAxParent_interaction-effect_unadjusted_2022-07-07.html"))

# Individual groups vs. others
glimmaMA(fit2.lxp, d.lxp.dummy,
         logCPM.lxp.filt,
         groups = targets.lxp$Group,
         coef = "(ACM + ATM) - (BCM + BTM + ACF + ATF + BCF + BTF)/3",
         main = "AM (EHB male) vs. others",
         transform.counts = "none",
         html = paste0("results/stats/glimma/glimma_MA_lincRNAxParent_AM_unadjusted_2022-07-07.html"))

glimmaMA(fit2.lxp, d.lxp.dummy,
         logCPM.lxp.filt,
         groups = targets.lxp$Group,
         coef = "(BCM + BTM) - (ACM + ATM + ACF + ATF + BCF + BTF)/3",
         main = "BM (AHB male) vs. others",
         transform.counts = "none",
         html = paste0("results/stats/glimma/glimma_MA_lincRNAxParent_BM_unadjusted_2022-07-07.html"))

glimmaMA(fit2.lxp, d.lxp.dummy,
         logCPM.lxp.filt,
         groups = targets.lxp$Group,
         coef = "(ACF + ATF) - (BCM + BTM + ACM + ATM + BCF + BTF)/3",
         main = "AF (AHB female) vs. others",
         transform.counts = "none",
         html = paste0("results/stats/glimma/glimma_MA_lincRNAxParent_AF_unadjusted_2022-07-07.html"))

glimmaMA(fit2.lxp, d.lxp.dummy,
         logCPM.lxp.filt,
         groups = targets.lxp$Group,
         coef = "(BCF + BTF) - (ACM + ATM + ACF + ATF + BCM + BTM)/3",
         main = "BF (EHB female) vs. others",
         transform.counts = "none",
         html = paste0("results/stats/glimma/glimma_MA_lincRNAxParent_BF_unadjusted_2022-07-07.html"))

# C vs. T interaction effects
glimmaMA(fit2.lxp, d.lxp.dummy,
         logCPM.lxp.filt,
         groups = targets.lxp$Group,
         coef = "(BTM - BCM) - (BTF - BCF)",
         main = "Interaction effect of aggression and parental allele in cross B",
         transform.counts = "none",
         html = paste0("results/stats/glimma/glimma_MA_lincRNAxParent_aggression-parent-interaction-crossB_unadjusted_2022-07-07.html"))

glimmaMA(fit2.lxp, d.lxp.dummy,
         logCPM.lxp.filt,
         groups = targets.lxp$Group,
         coef = "(ATM - ACM) - (ATF - ACF)",
         main = "Interaction effect of aggression and parental allele in cross A",
         transform.counts = "none",
         html = paste0("results/stats/glimma/glimma_MA_lincRNAxParent_aggression-parent-interaction-crossA_unadjusted_2022-07-07.html"))

# Search for parental bias found in cross B but not A ####
genesB <- rownames(res.dge.lxp)[res.dge.lxp[["adj.P.Val (BCM + BTM - ACM - ATM) - (BCF + BTF - ACF - ATF)"]] < 0.05 &
                                  res.dge.lxp[["adj.P.Val (BTM + BCM) - (BTF + BCF)"]] < 0.05 &
                                  res.dge.lxp[["adj.P.Val (ATM + ACM) - (ATF + ACF)"]] >= 0.05]

temp <-
  res.dge.lxp[genesB,
              c("ID", "Transcript", "Gene", "AveExpr",
                paste(rep(c("logFC", "P.Value", "adj.P.Val"), times = 3),
                      rep(c("(BCM + BTM - ACM - ATM) - (BCF + BTF - ACF - ATF)",
                            "(BTM + BCM) - (BTF + BCF)", "(ATM + ACM) - (ATF + ACF)"),
                          each = 3)))] # 30 genes

# write.table(temp, row.names = FALSE, sep = "\t",
#             file = "results/lincRNA_results/differential_lincRNA_expression_results_lincRNAxParent_AHBmale_2022-07-07.txt")
