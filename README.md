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
- [Salmon](https://combine-lab.github.io/salmon/) (Transcriptome-based)
- [Kallisto](https://pachterlab.github.io/kallisto/) (Transcriptome-based)
- [STAR](https://github.com/alexdobin/STAR) (Genome + Features File-based)
- [RSubread](https://bioconductor.org/packages/release/bioc/html/Rsubread.html) (Genome + Features File-based)

**Analyze with**
- [DEseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- [sleuth](https://pachterlab.github.io/sleuth/about)

**Compare at**
- Transcript-level
- Gene-level