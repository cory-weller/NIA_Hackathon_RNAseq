#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 100G
#SBATCH --time 4:00:00
#SBATCH --partition norm

# load STAR module, so that the STAR command will work
module load STAR
mkdir -p REF

# Retrieve GTF annotation if necessary
GTF='REF/Sscrofa11.1.gtf'
GTF_URL='https://ftp.ensembl.org/pub/release-110/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.110.gtf.gz'
[[ ! -f ${GTF} ]] && \
    wget -O ${GTF}.gz ${GTF_URL} && \
    gunzip ${GTF}.gz

# Retrieve reference genome if necessary
GENOME='REF/Sscrofa11.1.dna.toplevel.fa'
GENOME_URL='https://ftp.ensembl.org/pub/release-110/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz'
[[ ! -f "${GENOME}" ]] && \
    wget -O ${GENOME}.gz ${GENOME_URL} && \
    gunzip ${GENOME}.gz

INDEX='STAR_index'

mkdir -p ${INDEX}


# Validate arguments are provided and files exist
BADARGS=0
[[ -f ${GENOME} ]]      || { echo "Genome ${GENOME} does not exist!"; ((BADARGS+=1)); }
[[ -f ${GTF} ]]         || { echo "GTF ${GTF} does not exist!"; ((BADARGS+=1)); }
[[ ${BADARGS} == 0 ]]   || { echo "Exiting due to ${BADARGS} missing file(s)"; exit 1; }


STAR \
    --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir ${INDEX} \
    --genomeFastaFiles ${GENOME} \
    --sjdbGTFfile ${GTF} \
    --sjdbOverhang 149 \
    --genomeChrBinNbits 15
