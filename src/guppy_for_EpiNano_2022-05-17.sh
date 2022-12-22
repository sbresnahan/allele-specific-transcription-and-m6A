#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH -p hpcbio
#SBATCH -n 1
#SBATCH --mem=16g
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH --mail-user=lvclark@illinois.edu
#SBATCH --mail-type=ALL
#SBATCH -J guppy
#SBATCH --array=1-12%1
#SBATCH --output=guppy-%A_%a.out
#SBATCH -D  /home/groups/hpcbio/projects/libyarlay/libyarlay-oxfordRNA-2021Jun/src/slurm-out/
# ----------------Load Modules--------------------
# EpiNano is built on models from Guppy 3.1.5, and therefore needs base calls from that version
module load guppy-gpu/3.1.5
# ----------------Commands------------------------
# bash strict mode (fail if any commands fail)
set -euo pipefail
IFS=$'\n\t'

cd /home/groups/hpcbio/projects/libyarlay/libyarlay-oxfordRNA-2021Jun/

# set variables
Sample_ID=`head src/sample_list_May2022.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 1`

Subfolder=`head src/sample_list_May2022.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 2`

Fast5_base=`head src/sample_list_May2022.txt -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 3`

Fast5_folder=${Fast5_base}${Sample_ID}/${Subfolder}

OUT=results/guppy/2022-05-17_for_EpiNano

echo ${Sample_ID}

mkdir -p $OUT

# run guppy
guppy_basecaller \
-i ${Fast5_folder} --recursive \
-s ${OUT}/${Sample_ID} \
-c rna_r9.4.1_70bps_hac.cfg \
--device cuda:0 \
--compress_fastq --records_per_fastq 0
