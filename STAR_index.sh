#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 100G
#SBATCH --time 4:00:00
#SBATCH --partition norm

# load STAR module, so that the STAR command will work
module load STAR

# Create a directory for the genome
INDEX='/data/NIA_Hackathon_2023/RNA_sequencing/STAR_index'
mkdir -p ${INDEX}

GENOME='/data/NIA_Hackathon_2023/RNA_sequencing/ncbi_dataset/data/GCF_000003025.6/GCF_000003025.6_Sscrofa11.1_genomic.fna'
GFF='/data/NIA_Hackathon_2023/RNA_sequencing/ncbi_dataset/data/GCF_000003025.6/genomic.gff'

# Validate arguments are provided and files exist
BADARGS=0
[[ -f ${GENOME} ]]      || { echo "Genome ${GENOME} does not exist!"; ((BADARGS+=1)); }
[[ -f ${GFF} ]]         || { echo "GFF ${GFF} does not exist!"; ((BADARGS+=1)); }
[[ ${BADARGS} == 0 ]]   || { echo "Exiting due to ${BADARGS} missing file(s)"; exit 1; }


STAR \
    --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir ${INDEX} \
    --genomeFastaFiles ${GENOME} \
    --sjdbGTFfile ${GFF} \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbOverhang 149 \
    --genomeChrBinNbits 15
