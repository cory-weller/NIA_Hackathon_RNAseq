# NIA_Hackathon_RNAseq

## Connecting to Biowulf

You will need two things: a way to connect, and a way to save text files. Connection method will depend on whether you have a mac or windows machine. 

A plain-text editor will also serve you well... not word or google drive, as those save hidden/invisible characters. I recommend VScode. Easiest to use [Web version](https://vscode.dev/) until you can get a 
full version installed with adminstrator privileges by IT.

Take note of your **computing ID**, which is identical to your NIH email alias that ISN'T **lastname.firstname**@nih.gov.
It is typically up to ~8 letters of lastn ame + first initial + middle initial + maybe a number.
For example, user `John A. Smith` could have an nih email `john.smith3@nih.gov` and a computing ID `smithja3`.

## Connecting to Biowulf on windows

For connecting on a Windows machine, I recommend [MobaXterm](https://mobaxterm.mobatek.net/download-home-edition.html). Download the zip file, extract the folder somewhere on your computer, and launch the executable.

Right click in the side pane, then click New session. 
- choose the SSH tab
- for `Remote Host`, enter `biowulf.nih.gov`
- check `Specify username` and put your computing ID

Double click on the session to connect. It will prompt you for your password. If it doesn't, make sure
you are on the NIH network OR connected to NIHVPN.

## Connecting to Biowulf on a mac
For connecting on a Mac, you can use the built-in terminal app. From the terminal, run the following command (replacing your own ID instead of `computingid`):
```bash
ssh computingid@biowulf.nih.gov
```

## Starting an interactive session
When you first connect to biowulf, you are on a shared login node. Running stuff on the shared login node can degrade performance for all users. Before doing any real work, you should request an interactive session with the command `sinteractive`. You will wait until an interactive node is granted to you. The output will look something like this:
```bash
sinteractive
salloc: Pending job allocation 9737355
salloc: job 9737355 queued and waiting for resources
salloc: job 9737355 has been allocated resources
salloc: Granted job allocation 9737355
salloc: Waiting for resource configuration
salloc: Nodes cn4338 are ready for job
```
At which point you're ready to go.


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

## Editing scripts
The simplest way is to mount your `/data` drive following instructions here https://hpc.nih.gov/docs/helixdrive.html and then editing files using VScode.

[VScode web link](https://vscode.dev/)

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


## For hackathon day 2:
**Analyze with**
- [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [sleuth](https://pachterlab.github.io/sleuth/about)

**Compare at**
- Transcript-level
- Gene-level

**Visualize**