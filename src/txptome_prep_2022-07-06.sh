#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=16G
#SBATCH -p hpcbio
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lvclark@illinois.edu
#SBATCH -J epinano
#SBATCH -D /home/groups/hpcbio/projects/libyarlay/libyarlay-oxfordRNA-2021Jun/src/slurm-out

cd /home/groups/hpcbio/projects/libyarlay/libyarlay-oxfordRNA-2021Jun/

TranscriptomeB=results/parent_transcriptomes/2022-06-14/CrossB_transcriptome.fa
TranscriptomeA=results/parent_transcriptomes/2022-06-14/CrossA_transcriptome.fa

### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

#### Reference sequence indexing for EpiNano

### Modules
module purge
module load SAMtools/1.12-IGB-gcc-8.2.0
module load picard/2.10.1-Java-1.8.0_152

### Index
samtools faidx $TranscriptomeB
samtools faidx $TranscriptomeA

### Sequence dictionary
picard CreateSequenceDictionary R=$TranscriptomeB
picard CreateSequenceDictionary R=$TranscriptomeA

### Rename dictionaries -- EDIT if rerunning
cd /home/groups/hpcbio/projects/libyarlay/libyarlay-oxfordRNA-2021Jun/results/parent_transcriptomes/2022-06-14
mv CrossB_transcriptome.dict CrossB_transcriptome.fa.dict
mv CrossA_transcriptome.dict CrossA_transcriptome.fa.dict
