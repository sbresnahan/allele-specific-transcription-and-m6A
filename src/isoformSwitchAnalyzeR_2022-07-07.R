# Try IsoformSwitchAnalyzeR on Flair counts matrix
library(IsoformSwitchAnalyzeR)
library(Biostrings)

# Import counts and transcript info
mycounts <- read.table("results/flair/2022-06-14/counts_matrix_AvsB_2022-07-06.tsv",
                       sep = "\t", header = TRUE)
str(mycounts)

# Regular expressions for splitting up transcript labels
readidpat <- "^[[:xdigit:]]{8}-([[:xdigit:]]{4}-){3}[[:xdigit:]]{12}(-[[:digit:]])?"
txptpat <- "^[NX][MR]_[[:digit:]]+\\.[[:digit:]]+(-[[:digit:]])?"
pseudogene_pat <- "^gene[[:digit:]]+(-[[:digit:]])?"

unknowngene_pat <- "N[CW]_[[:digit:]]+\\.[[:digit:]]+:[[:digit:]]+$"
knowngene_pat <- "GeneID:[[:digit:]]+$"
viralgene_pat <- "(DW|ABP|IAP|SB|CBP|VD)V(:[[:digit::]+)?$"
rRNA_pat <- "l-rRNA$"

# Function to make a data frame of information about each transcript
# (modified from parental allele expression analysis)
makeTxptInfo <- function(ids){
  txptinfo <- data.frame(Label = ids,
                         FlairLab = ids,
                         Transcript = "",
                         Gene = "",
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
  
  return(txptinfo)
}

mytxpts <- makeTxptInfo(mycounts$ids)
anyDuplicated(mytxpts$Transcript) # No duplicate txpt names

countsmat <- mycounts[,-1]
colnames(countsmat) <- sub("_.*", "", colnames(countsmat))
countsmat <- data.frame(isoform_id = mytxpts$Transcript,
                        countsmat)


# Import annotation ####
# Genomic coordinates of isoforms
#isoform_anno <- read.table("results/flair/2022-06-14/All_samples_collapse.isoforms.bed",
#                           sep = "\t", header = FALSE)

# location of exons
gtf <- rtracklayer::import("results/flair/2022-06-14/All_samples_collapse.isoforms.gtf")
gtf1 <- gtf[gtf$type == "exon"]

# correct gene and txpt names
ncrows <- which(endsWith(gtf1$transcript_id, "_NC")) # 12785 rows
gtf1$gene_id[ncrows] <- paste0("NC_", gtf1$gene_id[ncrows])
gtf1$transcript_id[ncrows] <- sub("_NC$", "", gtf1$transcript_id[ncrows])

nwrows <- which(endsWith(gtf1$transcript_id, "_NW")) # 117 rows
gtf1$gene_id[nwrows] <- paste0("NW_", gtf1$gene_id[nwrows])
gtf1$transcript_id[nwrows] <- sub("_NW$", "", gtf1$transcript_id[nwrows])

setdiff(mytxpts$Transcript, gtf1$transcript_id)
setdiff(mytxpts$Gene, gtf1$gene_id)
# --> everything from the counts matrix is in the GTF

# Subset to what was counted
gtf1 <- gtf1[gtf1$transcript_id %in% mytxpts$Transcript] # 137699 exons

# Get the columns that IsoformSwitch wants
mcols(gtf1) <- DataFrame(isoform_id = gtf1$transcript_id,
                         gene_id = gtf1$gene_id)

gene_info <- readRDS("results/final_report/annotations.rds") # NCBI GTF

gtf1$gene_name <- gene_info$gene[match(gtf1$gene_id, gene_info$gene_id)]

# Design matrices ####
# Individual experimental groups
design1 <- data.frame(sampleID = colnames(countsmat)[-1])
design1$condition <-
  paste0(substring(design1$sampleID, 8, 8),
         substring(design1$sampleID, 10, 10))

design1

comparisons1 <- data.frame(condition_1 = c("AC", "AT", "BC", "AC"),
                           condition_2 = c("BC", "BT", "BT", "AT"))

# Cross B vs. Cross A
design2 <- data.frame(sampleID = colnames(countsmat)[-1])
design2$condition <- substring(design2$sampleID, 8, 8)
design2$aggression <- substring(design2$sampleID, 10, 10)

# C vs T
design3 <- data.frame(sampleID = colnames(countsmat)[-1])
design3$condition <- substring(design3$sampleID, 10, 10)
design3$cross <- substring(design3$sampleID, 8, 8)

# Interaction
design4 <- data.frame(sampleID = design1$sampleID,
                      condition = ifelse(design1$condition %in% c("AC", "BT"),
                                         "ACBT", "ATBC"))
design4

# FASTA file with corrected names ####
myfasta <- readDNAStringSet("results/flair/2022-06-14/All_samples_collapse.isoforms.fa")
all(mytxpts$FlairLab %in% names(myfasta)) # TRUE
myfasta <- myfasta[mytxpts$FlairLab]
names(myfasta) <- mytxpts$Transcript

#writeXStringSet(myfasta, "results/IsoformSwitchAnalyzeR/All_samples_collapse_relabeled_2022-07-07.isoforms.fa")

# Create SwitchAnalyzeRlist ####

sarlist1 <-
  importRdata(isoformCountMatrix = countsmat,
              designMatrix = design1,
              isoformExonAnnoation = gtf1,
              isoformNtFasta = "results/IsoformSwitchAnalyzeR/All_samples_collapse_relabeled_2022-07-07.isoforms.fa",
              comparisonsToMake = comparisons1,
              fixStringTieAnnotationProblem = FALSE,
              fixStringTieViaOverlapInMultiGenes = FALSE)
# The GUESSTIMATED number of genes with differential isoform usage are:
#   comparison estimated_genes_with_dtu
# 1   AC vs BC                    0 - 0
# 2   AT vs BT                    0 - 0
# 3   BC vs BT                  32 - 53

sarlist2 <-
  importRdata(isoformCountMatrix = countsmat,
              designMatrix = design2,
              isoformExonAnnoation = gtf1,
              isoformNtFasta = "results/IsoformSwitchAnalyzeR/All_samples_collapse_relabeled_2022-07-07.isoforms.fa",
              fixStringTieAnnotationProblem = FALSE,
              fixStringTieViaOverlapInMultiGenes = FALSE)
#   comparison estimated_genes_with_dtu
# 1     A vs B                    0 - 0

sarlist3 <-
  importRdata(isoformCountMatrix = countsmat,
              designMatrix = design3,
              isoformExonAnnoation = gtf1,
              isoformNtFasta = "results/IsoformSwitchAnalyzeR/All_samples_collapse_relabeled_2022-07-07.isoforms.fa",
              fixStringTieAnnotationProblem = FALSE,
              fixStringTieViaOverlapInMultiGenes = FALSE)
#   comparison estimated_genes_with_dtu
# 1     C vs T                    0 - 0

sarlist4 <-
  importRdata(isoformCountMatrix = countsmat,
              designMatrix = design4,
              isoformExonAnnoation = gtf1,
              isoformNtFasta = "results/IsoformSwitchAnalyzeR/All_samples_collapse_relabeled_2022-07-07.isoforms.fa",
              fixStringTieAnnotationProblem = FALSE,
              fixStringTieViaOverlapInMultiGenes = FALSE)
#     comparison estimated_genes_with_dtu
# 1 ACBT vs ATBC                    0 - 0

# Filtering ####
sarlist1 <- preFilter(sarlist1)
# The filtering removed 16039 ( 51% of ) transcripts. There is now 15411 isoforms left

sarlist2 <- preFilter(sarlist2)
# The filtering removed 16740 ( 53.23% of ) transcripts. There is now 14710 isoforms left

sarlist3 <- preFilter(sarlist3)
# The filtering removed 16560 ( 52.66% of ) transcripts. There is now 14890 isoforms left

sarlist4 <- preFilter(sarlist4)
# The filtering removed 16471 ( 52.37% of ) transcripts. There is now 14979 isoforms left

# Test for differential isoform usage ####
sarlist1_analyzed <- isoformSwitchTestDEXSeq(sarlist1,
                                    reduceToSwitchingGenes = TRUE,
                                    alpha = 0.1)

extractSwitchSummary(sarlist1_analyzed)
#   Comparison nrIsoforms nrSwitches nrGenes
# 1   AC vs AT          0          0       0
# 2   AC vs BC         13          9       8
# 3   AT vs BT         15         16      10
# 4   BC vs BT          8          6       4
# 5   Combined         27         28      15

extractSwitchSummary(sarlist1_analyzed, alpha = 0.1)
#   Comparison nrIsoforms nrSwitches nrGenes
# 1   AC vs AT          0          0       0
# 2   AC vs BC         13          9       8
# 3   AT vs BT         18         18      12
# 4   BC vs BT         11          9       6
# 5   Combined         31         31      18

sarlist2_analyzed <- isoformSwitchTestDEXSeq(sarlist2,
                                             reduceToSwitchingGenes = TRUE,
                                             alpha = 0.1)

extractSwitchSummary(sarlist2_analyzed)
#   Comparison nrIsoforms nrSwitches nrGenes
# 1     A vs B         21         16      13

extractSwitchSummary(sarlist2_analyzed, alpha = 0.1)
#   Comparison nrIsoforms nrSwitches nrGenes
# 1     A vs B         24         18      15

sarlist3_analyzed <- isoformSwitchTestDEXSeq(sarlist3,
                                             reduceToSwitchingGenes = TRUE,
                                             alpha = 0.1)

extractSwitchSummary(sarlist3_analyzed)
#   Comparison nrIsoforms nrSwitches nrGenes
# 1     C vs T          7          5       4

extractSwitchSummary(sarlist3_analyzed, alpha = 0.1)
#   Comparison nrIsoforms nrSwitches nrGenes
# 1     C vs T         14         11       9

sarlist4_analyzed <- isoformSwitchTestDEXSeq(sarlist4,
                                             reduceToSwitchingGenes = TRUE,
                                             alpha = 0.1)
# No genes switching at this cutoff

# Add open reading frames ####
sarlist1_analyzed <- analyzeORF(sarlist1_analyzed)

sarlist2_analyzed <- analyzeORF(sarlist2_analyzed)

sarlist3_analyzed <- analyzeORF(sarlist3_analyzed)

# Write out sequences and import results from external tools ####
sarlist1_analyzed <- extractSequence(sarlist1_analyzed, alpha = 0.1,
                 pathToOutput = "results/IsoformSwitchAnalyzeR/small_group_comparisons_2022-07-07",
                 removeLongAAseq = TRUE,
                 alsoSplitFastaFile = TRUE)
# no sequences filtered for being too long

sarlist2_analyzed <- extractSequence(sarlist2_analyzed, alpha = 0.1,
                                     pathToOutput = "results/IsoformSwitchAnalyzeR/A_vs_B_2022-07-07",
                                     removeLongAAseq = TRUE,
                                     alsoSplitFastaFile = TRUE)
# no sequences filtered for being too long

sarlist3_analyzed <- extractSequence(sarlist3_analyzed, alpha = 0.1,
                                     pathToOutput = "results/IsoformSwitchAnalyzeR/C_vs_T_2022-07-07",
                                     removeLongAAseq = TRUE,
                                     alsoSplitFastaFile = TRUE)
# no sequences filtered for being too long

# save(sarlist1_analyzed, sarlist2_analyzed, sarlist3_analyzed,
#      file = "results/IsoformSwitchAnalyzeR/analysis_2022-07-07.RData")
load("results/IsoformSwitchAnalyzeR/analysis_2022-07-07.RData")

# SignalP 6.0 at https://services.healthtech.dtu.dk/service.php?SignalP
# Had to use SignalP 5.0 at https://services.healthtech.dtu.dk/service.php?SignalP-5.0 since 6.0 is not supported
# Pfam at https://www.ebi.ac.uk/Tools/hmmer/
# NetSurfP 2.0 at https://services.healthtech.dtu.dk/service.php?NetSurfP-2.0
# CPC2 at http://cpc2.gao-lab.org/batch.php

sarlist1_analyzed <- analyzeCPC2(sarlist1_analyzed,
                                 "results/IsoformSwitchAnalyzeR/small_group_comparisons_2022-07-07/result_cpc2.txt",
                                 codingCutoff = 0.7, removeNoncodinORFs = TRUE)
# Added coding potential to 82 (100%) transcripts

sarlist1_analyzed <- analyzePFAM(sarlist1_analyzed,
                                 "results/IsoformSwitchAnalyzeR/small_group_comparisons_2022-07-07/hmmscan_pfam.txt")
# Added domain information to 33 (40.24%) transcripts

sarlist1_analyzed <- analyzeSignalP(sarlist1_analyzed,
                                    "results/IsoformSwitchAnalyzeR/small_group_comparisons_2022-07-07/SignalP5_output_protein_type.txt")
#  Added signal peptide information to 15 (18.29%) transcripts

sarlist1_analyzed <- analyzeNetSurfP2(sarlist1_analyzed,
                                      "results/IsoformSwitchAnalyzeR/small_group_comparisons_2022-07-07/NetSurfP2_62C734E30000796FC67C4B82.csv")
# Added IDR information to 17 (20.73%) transcripts


sarlist2_analyzed <- analyzeCPC2(sarlist2_analyzed,
                                 "results/IsoformSwitchAnalyzeR/A_vs_B_2022-07-07/result_cpc2.txt",
                                 codingCutoff = 0.7, removeNoncodinORFs = TRUE)
# Added coding potential to 71 (100%) transcripts

sarlist2_analyzed <- analyzePFAM(sarlist2_analyzed,
                                 "results/IsoformSwitchAnalyzeR/A_vs_B_2022-07-07/hmmer_pfam.txt")
# Added domain information to 17 (23.94%) transcripts

sarlist2_analyzed <- analyzeSignalP(sarlist2_analyzed,
                                    "results/IsoformSwitchAnalyzeR/A_vs_B_2022-07-07/SignalP5_output_protein_type.txt")
#  Added signal peptide information to 15 (21.13%) transcripts

sarlist2_analyzed <- analyzeNetSurfP2(sarlist2_analyzed,
                                      "results/IsoformSwitchAnalyzeR/A_vs_B_2022-07-07/NetSurfP2_62CB260C000012B09ECDA369.csv")

# Added IDR information to 13 (18.31%) transcripts


sarlist3_analyzed <- analyzeCPC2(sarlist3_analyzed,
                                 "results/IsoformSwitchAnalyzeR/C_vs_T_2022-07-07/result_cpc2.txt",
                                 codingCutoff = 0.7, removeNoncodinORFs = TRUE)
# Added coding potential to 35 (100%) transcripts

sarlist3_analyzed <- analyzePFAM(sarlist3_analyzed,
                                 "results/IsoformSwitchAnalyzeR/C_vs_T_2022-07-07/hmmscan_pfam.txt")
# Added domain information to 32 (91.43%) transcripts

sarlist3_analyzed <- analyzeSignalP(sarlist3_analyzed,
                                    "results/IsoformSwitchAnalyzeR/C_vs_T_2022-07-07/SignalP5_output_protein_type.txt")
#  Added signal peptide information to 3 (8.57%) transcripts 

sarlist3_analyzed <- analyzeNetSurfP2(sarlist3_analyzed,
                                      "results/IsoformSwitchAnalyzeR/C_vs_T_2022-07-07/NetSurfP2_62CB440B000050B79235E116.csv")
# Added IDR information to 7 (20%) transcripts

# Analyze alternative splicing ####
sarlist1_analyzed <- analyzeAlternativeSplicing(sarlist1_analyzed,
                                                alpha = 0.1)
# warning about transcripts not from same strand

sarlist1_analyzed <- analyzeSwitchConsequences(sarlist1_analyzed,
                                               alpha = 0.1,
                                               removeNonConseqSwitches = FALSE)

extractSwitchSummary(sarlist1_analyzed, filterForConsequences = FALSE, alpha = 0.1)
#   Comparison nrIsoforms nrSwitches nrGenes
# 1   AC vs AT          0          0       0
# 2   AC vs BC         13          9       8
# 3   AT vs BT         18         18      12
# 4   BC vs BT         11          9       6
# 5   Combined         31         31      18

extractSwitchSummary(sarlist1_analyzed, filterForConsequences = TRUE, alpha = 0.1)
#   Comparison nrIsoforms nrSwitches nrGenes
# 1   AC vs AT          0          0       0
# 2   AC vs BC          4          3       3
# 3   AT vs BT          6          7       5
# 4   BC vs BT          9          8       5
# 5   Combined         15         16       9

topswitches1 <- extractTopSwitches(sarlist1_analyzed, filterForConsequences = TRUE, alpha = 0.1, n = Inf)
#             gene_ref           gene_id gene_name condition_1 condition_2 gene_switch_q_value
# 1  geneComp_00022113     GeneID:406140     Apid1          AC          BC        5.234770e-54
# 3  geneComp_00062005     GeneID:406140     Apid1          BC          BT        2.299970e-24
# 7  geneComp_00042059     GeneID:406140     Apid1          AT          BT        6.716720e-15
# 8  geneComp_00042046     GeneID:406121     Mrjp3          AT          BT        1.463720e-10
# 10 geneComp_00028837 NC_001566.1:13000      <NA>          AC          BC        5.060549e-05
# 11 geneComp_00022100     GeneID:406121     Mrjp3          AC          BC        1.831075e-03
# 12 geneComp_00067461     GeneID:725215 LOC725215          BC          BT        3.288085e-03
# 14 geneComp_00063185     GeneID:409924 LOC409924          BC          BT        8.016077e-03
# 15 geneComp_00048570     GeneID:727037 LOC727037          AT          BT        1.162371e-02
# 16 geneComp_00043742     GeneID:410806 LOC410806          AT          BT        1.996403e-02
# 17 geneComp_00044940     GeneID:413141 LOC413141          AT          BT        7.090193e-02
# 18 geneComp_00068516     GeneID:727037 LOC727037          BC          BT        8.136820e-02
# 19 geneComp_00065993     GeneID:551806 LOC551806          BC          BT        8.522326e-02
#    switchConsequencesGene Rank
# 1                    TRUE    1
# 3                    TRUE    2
# 7                    TRUE    3
# 8                    TRUE    4
# 10                   TRUE    5
# 11                   TRUE    6
# 12                   TRUE    7
# 14                   TRUE    8
# 15                   TRUE    9
# 16                   TRUE   10
# 17                   TRUE   11
# 18                   TRUE   12
# 19                   TRUE   13

switchPlot(sarlist1_analyzed, gene = "GeneID:406140", condition1 = "AC", condition2 = "BC",
           alphas = c(0.1, 0.001))
switchPlot(sarlist1_analyzed, gene = "GeneID:406140", condition1 = "AT", condition2 = "BT",
           alphas = c(0.1, 0.001))

# A vs B
sarlist2_analyzed <- analyzeAlternativeSplicing(sarlist2_analyzed,
                                                alpha = 0.1)

sarlist2_analyzed <- analyzeSwitchConsequences(sarlist2_analyzed,
                                               alpha = 0.1,
                                               removeNonConseqSwitches = FALSE)

extractSwitchSummary(sarlist2_analyzed, filterForConsequences = FALSE, alpha = 0.1)
#   Comparison nrIsoforms nrSwitches nrGenes
# 1     A vs B         24         18      15

extractSwitchSummary(sarlist2_analyzed, filterForConsequences = TRUE, alpha = 0.1)
#   Comparison nrIsoforms nrSwitches nrGenes
# 1     A vs B          8          7       6

topswitches2 <- extractTopSwitches(sarlist2_analyzed, filterForConsequences = TRUE, alpha = 0.1, n = Inf)
#            gene_ref             gene_id gene_name condition_1 condition_2 gene_switch_q_value
# 1 geneComp_00002167       GeneID:406140     Apid1           A           B        0.0001549395
# 2 geneComp_00002154       GeneID:406121     Mrjp3           A           B        0.0002033795
# 4 geneComp_00014443 NC_037644.1:5664000      <NA>           A           B        0.0025668299
# 5 geneComp_00006951       GeneID:552829 LOC552829           A           B        0.0029525324
# 7 geneComp_00007246       GeneID:724565 LOC724565           A           B        0.0335776368
# 8 geneComp_00005048       GeneID:413141 LOC413141           A           B        0.0553019080
#   switchConsequencesGene Rank
# 1                   TRUE    1
# 2                   TRUE    2
# 4                   TRUE    3
# 5                   TRUE    4
# 7                   TRUE    5
# 8                   TRUE    6

# C vs T
sarlist3_analyzed <- analyzeAlternativeSplicing(sarlist3_analyzed,
                                                alpha = 0.1)

sarlist3_analyzed <- analyzeSwitchConsequences(sarlist3_analyzed,
                                               alpha = 0.1,
                                               removeNonConseqSwitches = FALSE)

extractSwitchSummary(sarlist3_analyzed, filterForConsequences = FALSE, alpha = 0.1)
#   Comparison nrIsoforms nrSwitches nrGenes
# 1     C vs T         14         11       9

extractSwitchSummary(sarlist3_analyzed, filterForConsequences = TRUE, alpha = 0.1)
#   Comparison nrIsoforms nrSwitches nrGenes
# 1     C vs T          8          6       5

topswitches3 <- extractTopSwitches(sarlist3_analyzed, filterForConsequences = TRUE, alpha = 0.1, n = Inf)
#            gene_ref       gene_id gene_name condition_1 condition_2 gene_switch_q_value switchConsequencesGene
# 1 geneComp_00003347 GeneID:409924 LOC409924           C           T        0.0004061671                   TRUE
# 3 geneComp_00002910 GeneID:409299 LOC409299           C           T        0.0008968435                   TRUE
# 5 geneComp_00005010 GeneID:413069 LOC413069           C           T        0.0905798021                   TRUE
# 6 geneComp_00005887 GeneID:551459 LOC551459           C           T        0.0905798021                   TRUE
# 8 geneComp_00006125 GeneID:551763 LOC551763           C           T        0.0905798021                   TRUE
#   Rank
# 1    1
# 3    2
# 5    3
# 6    4
# 8    5

# Results export ####
# cairo_pdf("results/IsoformSwitchAnalyzeR/Significant_switches_small_group_comparisons_2022-07-15.pdf",
#           width = 9, height = 7, onefile = TRUE)
for(i in seq_len(nrow(topswitches1))){
  switchPlot(sarlist1_analyzed, gene = topswitches1$gene_id[i],
             condition1 = topswitches1$condition_1[i],
             condition2 = topswitches1$condition_2[i],
             alphas = c(0.1, 0.001))
}
# dev.off()

head(sarlist1_analyzed$isoformFeatures)

# write.table(sarlist1_analyzed$isoformFeatures[,-(1:2)], sep = "\t", row.names = FALSE,
#             file = "results/IsoformSwitchAnalyzeR/Isoform_results_small_group_comparisons_2022-07-15.txt")

#saveRDS(sarlist1_analyzed, file = "results/IsoformSwitchAnalyzeR/sarlist1_analyzed_2022-07-15.rds")

## Put topswitches tables into report

# A vs B
# cairo_pdf("results/IsoformSwitchAnalyzeR/Significant_switches_A_vs_B_2022-07-15.pdf",
#           width = 9, height = 7, onefile = TRUE)
for(i in seq_len(nrow(topswitches2))){
  switchPlot(sarlist2_analyzed, gene = topswitches2$gene_id[i],
             condition1 = topswitches2$condition_1[i],
             condition2 = topswitches2$condition_2[i],
             alphas = c(0.1, 0.001))
}
# dev.off()

# write.table(sarlist2_analyzed$isoformFeatures[,-(1:2)], sep = "\t", row.names = FALSE,
#             file = "results/IsoformSwitchAnalyzeR/Isoform_results_A_vs_B_2022-07-15.txt")

#saveRDS(sarlist2_analyzed, file = "results/IsoformSwitchAnalyzeR/sarlist2_analyzed_2022-07-15.rds")

# C vs T
# cairo_pdf("results/IsoformSwitchAnalyzeR/Significant_switches_C_vs_T_2022-07-15.pdf",
#           width = 9, height = 7, onefile = TRUE)
for(i in seq_len(nrow(topswitches3))){
  switchPlot(sarlist3_analyzed, gene = topswitches3$gene_id[i],
             condition1 = topswitches3$condition_1[i],
             condition2 = topswitches3$condition_2[i],
             alphas = c(0.1, 0.001))
}
# dev.off()

# write.table(sarlist3_analyzed$isoformFeatures[,-(1:2)], sep = "\t", row.names = FALSE,
#             file = "results/IsoformSwitchAnalyzeR/Isoform_results_C_vs_T_2022-07-15.txt")

#saveRDS(sarlist3_analyzed, file = "results/IsoformSwitchAnalyzeR/sarlist3_analyzed_2022-07-15.rds")

# old code below ####

# Where do the sig genes rank in terms of overall expression? ####
gene_counts <- readRDS("results/stats/gene_counts_2022-01-05.rds")
str(gene_counts)
logCPM <- read.table("results/stats/filtered_logCPM_2022-01-05.txt",
                     header = TRUE, sep = "\t")
logCPM$Mean <- rowMeans(as.matrix(logCPM[,-(1:3)]))

summary(logCPM$Mean)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2.435   2.943   3.611   4.068   4.698  16.631 

hits <- unique(c(topswitches1$gene_id, topswitches2$gene_id))
hits
# Percentiles for expression level
sapply(logCPM$Mean[match(hits, logCPM$Gene)],
       function(x) mean(logCPM$Mean < x))
#0.9844222 0.9988078 0.9480210 0.9210777 0.9535050 0.5652519 0.8457320 0.7416945 0.4232237
