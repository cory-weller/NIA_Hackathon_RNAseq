#!/usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --time=1:00:00
#SBATCH --mem=32g


module load subread

OUTDIR='FEATURE_COUNTS'
GTF='REF/Sscrofa11.1.gtf'

TISSUE=${1}
SAMPLENO=${2}
BAM=STAR_BAM/${TISSUE}_${SAMPLENO}_Aligned.sortedByCoord.out.bam

mkdir -p ${OUTDIR}

echo ${BAM}

featureCounts -p --countReadPairs -O -T 8 \
              -t gene \
              -g gene_id \
              -a ${GTF} \
              -o ${OUTDIR}/${TISSUE}_${SAMPLENO}.featureCounts.txt \
              ${BAM}
