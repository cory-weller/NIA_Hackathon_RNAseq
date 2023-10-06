#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 12
#SBATCH --mem 50G
#SBATCH --time 1:00:00
#SBATCH --partition quick,norm

treatment=${1}
sample=${2}

echo treatment is ${treatment}
echo sample is ${sample}

DATADIR='/data/NIA_Hackathon_2023/RNA_sequencing/data'
FQ1=${DATADIR}/${treatment}_${sample}_R1.clean.fq.gz
FQ2=${DATADIR}/${treatment}_${sample}_R2.clean.fq.gz
GENOME='/data/NIA_Hackathon_2023/RNA_sequencing/STAR_index'
OUTDIR='/data/NIA_Hackathon_2023/RNA_sequencing/STAR_BAM'

mkdir -p ${OUTDIR} 

echo FQ1 is ${FQ1}
echo FQ2 is ${FQ2}

module load STAR

STAR \
    --runThreadN 12 \
    --genomeDir ${GENOME} \
    --sjdbOverhang 149 \
    --readFilesIn ${FQ1} ${FQ2} \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "${OUTDIR}/${treatment}_${sample}_"