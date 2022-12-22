# Examining the role of transcription and post-transcriptional modifications on mediating intragenomic conflict in honey bees

## Study design

This study uses whole genome resequencing data produced by Illumina DNA sequencing and transcriptome data produced by Oxford Nanopore direct RNA sequencing to understand behavioral differences between non-aggressive and aggressive honey bees. Two crosses between Africanized honey bees (AHB) and European honey bees (EHB) were made. Cross A had AHB as the female parent (SRR14654188) and EHB as the male parent (SRR14654189). Cross B had EHB as the female parent (SRR14654186) and AHB as the male parent (SRR14654187). Three or more aggressive bees and three or more non-aggressive bees (labeled C and T) were collected from each cross. Twelve head and abdomen samples, spanning multiple batches over the course of a year, were ultimately analyzed, with three samples per experimental group. RNA samples were sequenced with Oxford Nanopore Technology direct RNA sequencing on a FLO-MIN106 flowcell.

In this repository we provide scripts for the detection of and relationship between allele-specific transcription and RNA m6A, and the detection of differential gene expression and differential isoform usage. We test two main hypotheses:

1. Aggression is associated with an enrichment for paternal allele-biased transcription,

2. Allele-biased transcription is associated with allele-biased RNA m6A

Additionally, we assess whether either allele-biased transcription or allele-biased RNA m6A are associated with differential gene expression and/or differential isoform usage.

## Methods

### Call parent SNPs on the Amel_HAv3.1 reference genome
##### > [call_parent_SNPs.md](src/call_parent_SNPs.md)

### Create cross transcriptomes
#### Base calling using Guppy
##### > [guppy_for_EpiNano_2022-05-17.sh](src/guppy_for_EpiNano_2022-05-17.sh)
#### Prepare reference genome
##### > [minimap2_index_2022-06-14.sh](src/minimap2_index_2022-06-14.sh)
#### Align with Minimap2 and correct splice junctions & call transcript isoforms using Flair
##### > [concatenate_and_align_2021-07-26.sh](src/concatenate_and_align_2021-07-26.sh)
#### Create parent transcriptomes
##### > [generate_parent_transcripts.R](src/generate_parent_transcripts.R)
#### Combine parent transcriptomes to create cross transcriptomes
##### > [combine_transcriptomes.R](src/combine_transcriptomes.R)

### Align to cross transcriptomes with EpiNano
#### Prepare cross transcriptomes
##### > [txptome_prep_2022-07-06.sh](src/txptome_prep_2022-07-06.sh)
#### Align and predict m6A with EpiNano
#### > [epinano_array_script_2022-07-15.sh](src/epinano_array_script_2022-07-15.sh)
#### Combine EpiNano results into a single table of coverage and methylation probabilities
##### > [combine_results_2022-08-03.R](src/combine_results_2022-08-03.R)

### Detection of lincRNAs
##### > [lincRNAs_2022-07-01.R](src/lincRNAs_2022-07-01.R)

### Detection of isoform switching
##### > [isoformSwitchAnalyzeR_2022-07-07.R](src/isoformSwitchAnalyzeR_2022-07-07.R)

### Differential gene expression analysis
##### > [libyarlay_analysis_2022June_parenttxpt.R](src/libyarlay_analysis_2022June_parenttxpt.R)

### Analysis of parent-of-origin allele-specific transcription and m6A, and their relationship to gene expression and isoform switching
##### > Annotated markdown - [PSm6A.R](src/PSm6A.R)

## Results

### [Identification and differential expression of lincRNAs](https://sbresnahan.github.io/allele-specific-transcription-and-m6A/reports/report_lincRNAs.html)

### [Differential gene expression](https://sbresnahan.github.io/allele-specific-transcription-and-m6A/reports/report_Oxford_RNAseq_QC_DGE.html)

### [Isoform switching](https://sbresnahan.github.io/allele-specific-transcription-and-m6A/reports/report_isoform_switch.html)

### [Parent-of-origin allele-specific transcription and m6A](https://sbresnahan.github.io/allele-specific-transcription-and-m6A/reports/report_PSm6A.html)