#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH -p hpcbio
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lvclark@illinois.edu
#SBATCH -J align
#SBATCH --output=minimap2-%A_%a.out
#SBATCH --array=1-11%3
#SBATCH -D /home/groups/hpcbio/projects/libyarlay/libyarlay-oxfordRNA-2021Jun/src/slurm-out

### This script runs everything that needs to be run on individual samples after base calling.

cd /home/groups/hpcbio/projects/libyarlay/libyarlay-oxfordRNA-2021Jun/

Sample_ID=`head src/sample_list.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 1`
File_folder=`head src/sample_list.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 2`

Date='2021-07-26'

### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

### Load Modules
module load minimap/2.21-IGB-gcc-8.2.0
module load SAMtools/1.12-IGB-gcc-8.2.0
module load FastQC/0.11.8-Java-1.8.0_152
#module load Subread/2.0.0-IGB-gcc-8.2.0
module load flair/1.5-IGB-gcc-8.2.0-Python-3.7.2

### Concatenate FASTQs

cat ${File_folder}/fastq_pass/*.fastq > results/concatenated_fastq/${Date}/${Sample_ID}_all.fastq

### FASTQC

fastqc -o results/fastqc/${Date} -t $SLURM_NTASKS \
  results/concatenated_fastq/${Date}/${Sample_ID}_all.fastq

### Do alignment

minimap2 -ax splice -uf -k14 \
  --junc-bed data/reference/GCF_003254395.2_Amel_HAv3.1+virus_genomic.filled.bed \
  data/reference/GCF_003254395.2_Amel_HAv3.1+virus_genomic.k14.mmi \
  results/concatenated_fastq/${Date}/${Sample_ID}_all.fastq | \
  samtools sort -O bam - > results/minimap/${Date}/${Sample_ID}_all_sorted.bam

samtools index results/minimap/${Date}/${Sample_ID}_all_sorted.bam

### featureCounts for MultiQC # comment out for now since it gives an error

#featureCounts -F GTF -t exon -g gene_id -s 1 -T $SLURM_NTASKS \
#-a data/reference/GCF_003254395.2_Amel_HAv3.1+virus_genomic.filled.gtf \
#-o results/featureCounts/${Date}/${Sample_ID}_counts.txt \
#results/minimap/${Date}/${Sample_ID}_all_sorted.bam

### convert BAM to BED for Flair

bam2Bed12.py -i results/minimap/${Date}/${Sample_ID}_all_sorted.bam > \
  results/flair/${Date}/${Sample_ID}_all_sorted.bed

### correct junctions with Flair

flair.py correct -q results/flair/${Date}/${Sample_ID}_all_sorted.bed \
-g data/reference/GCF_003254395.2_Amel_HAv3.1+virus_genomic.fna \
-f data/reference/GCF_003254395.2_Amel_HAv3.1+virus_genomic.filled.gtf \
-o results/flair/${Date}/${Sample_ID}
