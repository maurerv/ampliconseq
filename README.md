# Process amplicon sequencing data using snakemake

## Overview

This pipeline uses the following tools to align paired end amplicon sequencing data:

1. [FastQC](https://github.com/s-andrews/FastQC), to assess sequencing data quality
2. [Trim Galore](https://github.com/FelixKrueger/TrimGalore) to trim reads using nextera adapters
3. [STAR](https://github.com/alexdobin/STAR) for gapped alignment
4. [samtools](https://github.com/samtools/) for alignment statistics and indexing
5. [MultiQC](https://multiqc.info/) to compile everything into a neat report


## Configuration

All necessary tools can be installed using conda and the provided environment file as per:
```
conda env create -f env.yml
```
The environment specification file does not control the versions of the utilized tools, so its up to the user to keep versions consistent across datasets.


## Running the pipeline

Adapting the pipeline to your use cases requires at least two changes to the config.yml. ```workdir_top```
is the absolute path to the working directory and STAR_reference specifies the absolute path ot the reference
created by the STAR aligner as per the [manual](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).

The working directory needs to contain a folder named 00_fastqc which contains either the fastq files or symbolic links
to them. Only files with suffix "\_R\[0-9\].fastq.gz" will be recognized by the pipeline.


## Result

After running the pipeline, the working directory will have the following structure
```
workdir_top
├── 00_fastqc/            # Contains the fastqc analysis of files in 01_fastq_ra
├── 01_fastq_raw/         # Contains the raw fastqs
├── 02_fastq_trimmed/     # Contains the trimmed fastqs
├── 03_star/              # Contains the alignment files from STAR
├── multiqc_data/         # Contains data used by multiqc to generate the final report
├── multiqc_report.html   # Final per sample quality report
```
