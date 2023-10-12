# NIA RNAseq Hackathon README

## Overview


Helpful tutorials:
- [Connecting to Biowulf](docs/biowulf.md)
- [Mounting Biowulf disks](docs/mount_disk.md)
- [Connecting RStudio Server](docs/rstudio_server.md)
- [Connecting VScode](docs/vscode.md)


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
See [STAR_index.sh](STAR_index.sh) which was submitted via `sbatch scripts/STAR_index.sh`. The script downloads and unpacks required reference data (genomic `fasta` and annotation `gtf`).

## Align `fastq` files to the indexed genome
[`STAR_wrappers.sh`](STAR_wrappers.sh) contains three example functions that will submit all 10 jobs.

```bash
# Load the file with 'source' so you can use the functions
source STAR_wrappers.sh

# submit jobs with the GNU parallel funciton defined in STAR_wrappers.sh
gnu_parallel
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

```
