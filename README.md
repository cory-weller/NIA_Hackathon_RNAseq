# NIA RNAseq Hackathon README

## Overview



[Connecting to Biowulf](docs/biowulf.md)
[Mounting Biowulf disks](docs/mount_disk.md)
[Connecting RStudio Server](docs/rstudio_server.md)
[Connecting VScode](docs/vscode.md)


## Overview for hackathon day 1
- get required tools
- connect to biowulf
- set up data mounted drive
- navigate directories
- find our data
- RNA seq broad overview
- specific RNAseq tools and how they differ
- retrieving reference from web
- preparing reference
- aligning to reference overview
- tiny data overview



## RNAseq permutations
**Align with**:
- [Salmon](https://combine-lab.github.io/salmon/) (Transcriptome-based) [Biowulf guide](https://hpc.nih.gov/apps/salmon.html)
- [Kallisto](https://pachterlab.github.io/kallisto/) (Transcriptome-based) [Biowulf guide](https://hpc.nih.gov/apps/kallisto.html)
- [STAR](https://github.com/alexdobin/STAR) (Genome + Features File-based) [Biowulf guide](https://hpc.nih.gov/apps/STAR.html)
- [RSubread](https://bioconductor.org/packages/release/bioc/html/Rsubread.html) (Genome + Features File-based) [Biowulf guide](https://hpc.nih.gov/apps/subread.html)

## Data location
Files are located at `/data/NIA_Hackathon_2023/RNA_sequencing` which contains two subfolders, organized as follows:
```
.
├── data
│   ├── Foliate_1_R1.clean.fq.gz
│   ├── Foliate_1_R2.clean.fq.gz
│   ├── Foliate_2_R1.clean.fq.gz
│   ├── Foliate_2_R2.clean.fq.gz
│   ├── Foliate_3_R1.clean.fq.gz
│   ├── Foliate_3_R2.clean.fq.gz
│   ├── Foliate_4_R1.clean.fq.gz
│   ├── Foliate_4_R2.clean.fq.gz
│   ├── Foliate_5_R1.clean.fq.gz
│   ├── Foliate_5_R2.clean.fq.gz
│   ├── Organoid_1_R1.clean.fq.gz
│   ├── Organoid_1_R2.clean.fq.gz
│   ├── Organoid_2_R1.clean.fq.gz
│   ├── Organoid_2_R2.clean.fq.gz
│   ├── Organoid_3_R1.clean.fq.gz
│   ├── Organoid_3_R2.clean.fq.gz
│   ├── Organoid_4_R1.clean.fq.gz
│   ├── Organoid_4_R2.clean.fq.gz
│   ├── Organoid_5_R1.clean.fq.gz
│   └── Organoid_5_R2.clean.fq.gz
└── tiny
    ├── tinyFoliate_1_R1.clean.fq.gz
    ├── tinyFoliate_1_R2.clean.fq.gz
    ├── tinyFoliate_2_R1.clean.fq.gz
    ├── tinyFoliate_2_R2.clean.fq.gz
    ├── tinyFoliate_3_R1.clean.fq.gz
    ├── tinyFoliate_3_R2.clean.fq.gz
    ├── tinyFoliate_4_R1.clean.fq.gz
    ├── tinyFoliate_4_R2.clean.fq.gz
    ├── tinyFoliate_5_R1.clean.fq.gz
    ├── tinyFoliate_5_R2.clean.fq.gz
    ├── tinyOrganoid_1_R1.clean.fq.gz
    ├── tinyOrganoid_1_R2.clean.fq.gz
    ├── tinyOrganoid_2_R1.clean.fq.gz
    ├── tinyOrganoid_2_R2.clean.fq.gz
    ├── tinyOrganoid_3_R1.clean.fq.gz
    ├── tinyOrganoid_3_R2.clean.fq.gz
    ├── tinyOrganoid_4_R1.clean.fq.gz
    ├── tinyOrganoid_4_R2.clean.fq.gz
    ├── tinyOrganoid_5_R1.clean.fq.gz
    └── tinyOrganoid_5_R2.clean.fq.gz

2 directories, 40 files
```

The `tiny` files are a random subset of reads, representing about 2% of the full data.
The small size files are useful for testing scripts to make sure they work, instead of the big files. Once you know your workflow is correctly built, you can instead use the full-sized files.

## Create STAR genome index
See [STAR_index.sh](STAR_index.sh) which was submitted via `sbatch STAR_index.sh`

## Align `fastq` files to the indexed genome
[`STAR_wrappers.sh`](STAR_wrappers.sh) contains three example functions that will submit all 10 jobs.

```bash
# Load the file with 'source' so you can use the functions
source STAR_wrappers.sh

# submit jobs with the GNU parallel funciton defined in STAR_wrappers.sh
gnu_parallel
```

## Index `bam` files
```bash
cd /data/NIA_Hackathon_2023/RNA_sequencing/STAR_BAM
module load samtools

for file in *.bam; do
    samtools index ${file}
done
```

## Subset `bam` files
Use `samtools view` with `-s` to subset a fraction of reads. `-s 0.02` subets to 2% of all reads.
The parameter `-b` forces the output to `bam` instead of `sam`.
```bash
module load samtools
i=0
for file in *.bam; do
    ((i+=1))
    echo $i
    [[ ! -f "tiny${file}" ]] && { samtools view -b -s 0.02 ${file} > "tiny${file}" &&  samtools index "tiny${file}"; }
done

samtools index tinyFoliate_1_Aligned.sortedByCoord.out.bam
```


## For hackathon day 2:
**Analyze with**
- [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [sleuth](https://pachterlab.github.io/sleuth/about)

**Compare at**
- Transcript-level
- Gene-level

**Visualize**