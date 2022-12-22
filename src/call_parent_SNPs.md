# Calling parent SNPs on the Amel_HAv3.1 reference genome

### Sample metadata

|    SRA     | Cross | Parent | Lineage |
| ---------- | ----- | ------ | ------- |
| SRR3037350 |  875  |    Q   |   EHB   |
| SRR3037351 |  875  |    D   |   AHB   |
| SRR3037352 |  888  |    Q   |   AHB   |
| SRR3037353 |  888  |    D   |   EHB   |
| SRR3037354 |  882  |    Q   |   EHB   |
| SRR3037355 |  882  |    D   |   AHB   |
| SRR3037356 |  894  |    Q   |   AHB   |
| SRR3037357 |  894  |    D   |   EHB   |

## Retrieve F0 DNA-seq libraries [`sra-tools`]

```
cd ${DIR_SRA}

SRA=("SRR3037350" "SRR3037351" "SRR3037352" "SRR3037353" \
     "SRR3037354" "SRR3037355" "SRR3037356" "SRR3037357")

for i in "${SRA[@]}"
do
  prefetch -O ${DIR_SRA} ${i}
  fasterq-dump -O ${DIR_SRA} ${DIR_SRA}/${i}.sra
  rm {i}.sra
done
```

## Trim adapters from parent WGS reads [`fastp`]

```
for i in "${SRA[@]}"
do
  fastp -i ${i}_1.fastq -I ${i}_2.fastq \
  -o ${DIR_TRIM}/${i}_1.fastq -O ${DIR_TRIM}/${i}_2.fastq
done
```

## Generate BWA index for the Amel_HAv3.1 reference genome [`bwa index`]

```
cd ${DIR_INDEX}

wget -O Amel_HAv3.1.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz
gunzip Amel_HAv3.1.fna.gz

bwa index Amel_HAv3.1.fna
```

## Align parent WGS reads [`bwa mem`] convert to BAM and sort [`samtools`]

```
cd ${DIR_SRA}

for i in "${SRA[@]}"
do
  bwa mem ${DIR_INDEX}/Amel_HAv3.1.fna ${i}_1.fastq ${i}_2.fastq > ${DIR_ALIGN}/${i}.sam
done

cd ${DIR_ALIGN}

for i in "${SRA[@]}"
do
  samtools view -O BAM ${i}.sam | samtools sort -@ 8 -O BAM > ${i}.bam
done
```

## Call variants [`freebayes`]

```
DIPLOID=("SRR3037350" "SRR3037352" "SRR3037354" "SRR3037356")
HAPLOID=("SRR3037351" "SRR3037353" "SRR3037355" "SRR3037357")

cd ${DIR_INDEX}

samtools faidx Amel_HAv3.1.fna

cd ${DIR_ALIGN}

for i in "${DIPLOID[@]}"
do
  freebayes -f ${DIR_INDEX}/Amel_HAv3.1.fna ${i}.bam > ${DIR_VARIANTS}/${i}.vcf
done

for i in "${HAPLOID[@]}"
do
  freebayes -f ${DIR_INDEX}/Amel_HAv3.1.fna ${i}.bam -p 1 > ${DIR_VARIANTS}/${i}.vcf
done
```

## Filter heterozygous variants [`bcftools`]

```
cd ${DIR_VARIANTS}

for i in "${SRA[@]}"
do
  bcftools filter -e 'GT="het"' ${i}.vcf > filtered/${i}.vcf
  bgzip -c filtered/${i}.vcf > filtered/${i}.vcf.gz
done
```
