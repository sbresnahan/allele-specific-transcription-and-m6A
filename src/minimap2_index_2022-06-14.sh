#!/bin/bash
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=36G
#SBATCH --mail-user=lvclark@illinois.edu
#SBATCH -J minimapindex
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -D /home/groups/hpcbio/projects/libyarlay/libyarlay-oxfordRNA-2021Jun/src/slurm-out
#SBATCH -p hpcbio

module load minimap/2.21-IGB-gcc-8.2.0

cd /home/groups/hpcbio/projects/libyarlay/libyarlay-oxfordRNA-2021Jun/data/reference/

minimap2 -k14 -d GCF_003254395.2_Amel_HAv3.1+virus_genomic_uppercase.k14.mmi \
  GCF_003254395.2_Amel_HAv3.1+virus_genomic_uppercase.fna
