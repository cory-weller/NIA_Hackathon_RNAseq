#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 12
#SBATCH --mem 50G
#SBATCH --time 1:00:00
#SBATCH --partition quick,norm

treatment=${1}
sample=${2}

# Print treatment and sample
echo treatment is ${treatment}
echo sample is ${sample}

# Define directories and files
DATADIR='data'
FQ1=${DATADIR}/${treatment}_${sample}_R1.clean.fq.gz
FQ2=${DATADIR}/${treatment}_${sample}_R2.clean.fq.gz
INDEX='STAR_index'
OUTDIR='STAR_BAM'

# Create output directory
mkdir -p ${OUTDIR} 

# Print fastq names
echo FQ1 is ${FQ1}
echo FQ2 is ${FQ2}

# Load STAR module
module load STAR

# Run alignment
STAR \
    --runThreadN 12 \
    --genomeDir ${INDEX} \
    --sjdbOverhang 149 \
    --readFilesIn ${FQ1} ${FQ2} \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "${OUTDIR}/${treatment}_${sample}_"