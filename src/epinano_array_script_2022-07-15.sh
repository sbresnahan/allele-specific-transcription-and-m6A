#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=32G
#SBATCH -p hpcbio
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lvclark@illinois.edu
#SBATCH -J epinano
#SBATCH --array=1-12%3
#SBATCH -D /home/groups/hpcbio/projects/libyarlay/libyarlay-oxfordRNA-2021Jun/src/slurm-out

cd /home/groups/hpcbio/projects/libyarlay/libyarlay-oxfordRNA-2021Jun/

Sample_ID=`head src/sample_list_May2022.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 1`
Cross=`head src/sample_list_May2022.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 4`

Guppy=results/guppy/2022-05-17_for_EpiNano
Date=2022-07-06
FlairDate=2022-06-14

Transcriptome=results/parent_transcriptomes/${FlairDate}/Cross${Cross}_transcriptome.fa
RRACH=results/epinano/Cross${Cross}_RRACH_Aonly_2022-07-06.rds

### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

### Modules
module purge
module load minimap/2.21-IGB-gcc-8.2.0
module load SAMtools/1.12-IGB-gcc-8.2.0
module load singularity/3.8.1

### Concatenate FASTQs
mkdir -p results/concatenated_fastq/${Date}
cat ${Guppy}/${Sample_ID}/*.fastq.gz > results/concatenated_fastq/${Date}/${Sample_ID}_all.fastq.gz

### Do alignment
mkdir -p results/minimap/${Date}

minimap2 -ax map-ont -uf --secondary=no -t $SLURM_NTASKS \
  $Transcriptome \
  results/concatenated_fastq/${Date}/${Sample_ID}_all.fastq.gz | \
  samtools sort -O bam - > results/minimap/${Date}/${Sample_ID}_all_sorted.bam

samtools index results/minimap/${Date}/${Sample_ID}_all_sorted.bam

### Run Epinano
mkdir -p results/epinano/${Date}

singularity exec -e epi12_latest.sif python3 /usr/local/bin/EpiNano/Epinano_Variants.py \
  -R $Transcriptome --type t \
  -b results/minimap/${Date}/${Sample_ID}_all_sorted.bam \
  -s /usr/local/bin/EpiNano/misc/sam2tsv.jar -n $SLURM_NTASKS

### Subset Epinano_Variants output to RRACH sites and filter by coverage
module load R/4.1.2-IGB-gcc-8.2.0

Rscript src/epinano/subset_epinano_variants_batch.R \
  results/minimap/${Date}/${Sample_ID}_all_sorted.plus_strand.per.site.csv \
  $RRACH \
  results/epinano/${Date}/${Sample_ID}_all_sorted.plus_strand.RRACH.per.site.csv \
  10

### Predict m6A with EpiNano, using modified version of script with bug fix
module purge
module load singularity/3.8.1
singularity exec -e epi12_latest.sif python3 src/epinano/Epinano_Predict.py \
  --model /usr/local/bin/EpiNano/models/rrach.q3.mis3.del3.linear.dump \
  --predict results/epinano/${Date}/${Sample_ID}_all_sorted.plus_strand.RRACH.per.site.csv \
  --columns 7,9,11 \
  --out_prefix results/epinano/${Date}/${Sample_ID}.modification

### Add genomic coordinates and kmers
module load R/4.1.2-IGB-gcc-8.2.0

Rscript src/epinano/update_coordinates_batch.R \
  results/epinano/${Date}/${Sample_ID}.modification.q_median.mis.del.MODEL.rrach.q3.mis3.del3.linear.dump.csv \
  $FlairDate $Cross
