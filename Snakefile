import os
from os import path

configfile: "./config.yml"
workdir: config["workdir_top"]

import sys

def message(mes):
  sys.stderr.write("|--- " + mes + "\n")

import glob
import re

SAMPLES, = glob_wildcards("01_fastq_raw/{sam}_R1.fastq.gz")

rule all:
    input:
        expand("00_fastqc/{sample}_R1_fastqc.zip", sample = SAMPLES),
        expand("03_star/{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand("03_star/{sample}.Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        expand("03_star/{sample}.Aligned.sortedByCoord.out.bam.stats", sample=SAMPLES),
        "multiqc_report.html"

rule fastqc:
    message:
        "Fastqc statistics"
    input:
        fwd = "01_fastq_raw/{sample}_R1.fastq.gz",
        rev = "01_fastq_raw/{sample}_R2.fastq.gz"
    output:
        fwd = "00_fastqc/{sample}_R1_fastqc.zip",
        rev = "00_fastqc/{sample}_R2_fastqc.zip",
    conda:
        "env.yml"
    resources:
        avg_mem  = lambda wildcards, attempt: 800 * attempt,
        mem_mb   = lambda wildcards, attempt: 1000 * attempt,
        walltime = lambda wildcards, attempt: 20 * attempt,
        attempt  = lambda wildcards, attempt: attempt,
    shell:"""
        fastqc {input.fwd} {input.rev} -o 00_fastqc/
        """

rule trim:
    message:
        "Preprocessing of fastq files"
    input:
        fwd = "01_fastq_raw/{sample}_R1.fastq.gz",
        rev = "01_fastq_raw/{sample}_R2.fastq.gz"
    output:
        fwd = "02_fastq_trimmed/{sample}_R1_val_1.fq.gz",
        rev = "02_fastq_trimmed/{sample}_R2_val_2.fq.gz"
    conda:
        "env.yml"
    resources:
        avg_mem  = lambda wildcards, attempt: 1200 * attempt,
        mem_mb   = lambda wildcards, attempt: 1300 * attempt,
        walltime = lambda wildcards, attempt: 20 * attempt,
        attempt  = lambda wildcards, attempt: attempt,
    shell:"""
        trim_galore --paired --nextera  {input.fwd} {input.rev} -o 02_fastq_trimmed/
        touch {output.fwd}
        touch {output.rev}
        """

rule align:
    message:
        "Align trimmed fastq files"
    input:
        fwd = "02_fastq_trimmed/{sample}_R1_val_1.fq.gz",
        rev = "02_fastq_trimmed/{sample}_R2_val_2.fq.gz"
    output:
        "03_star/{sample}.Aligned.sortedByCoord.out.bam"
    conda:
        "env.yml"
    params:
        odir = config["workdir_top"] + "/03_star/" ,
        ref = config["STAR_reference"]
    threads : 4
    resources:
        avg_mem  = lambda wildcards, attempt: 30000 * attempt,
        mem_mb   = lambda wildcards, attempt: 35000 * attempt,
        walltime = lambda wildcards, attempt: 60 * attempt,
        attempt  = lambda wildcards, attempt: attempt,
    shell:"""
        STAR --runThreadN {threads} --genomeDir {params.ref} \
         --readFilesCommand zcat \
         --outFileNamePrefix {params.odir}/{wildcards.sample}. \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --outSAMattributes Standard \
         --readFilesIn {input.fwd} {input.rev}
        """

rule index:
    message:
        "Index bam files"
    input:
        bam = "03_star/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "03_star/{sample}.Aligned.sortedByCoord.out.bam.bai"
    conda:
        "env.yml"
    resources:
        avg_mem  = lambda wildcards, attempt: 1200 * attempt,
        mem_mb   = lambda wildcards, attempt: 1300 * attempt,
        walltime = lambda wildcards, attempt: 20 * attempt,
        attempt  = lambda wildcards, attempt: attempt,
    shell:"""
        samtools index {input.bam}
        """

rule stats:
    message:
        "Running samtools stats"
    input:
        bam = "03_star/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "03_star/{sample}.Aligned.sortedByCoord.out.bam.stats"
    conda:
        "env.yml"
    resources:
        avg_mem  = lambda wildcards, attempt: 200 * attempt,
        mem_mb   = lambda wildcards, attempt: 300 * attempt,
        walltime = lambda wildcards, attempt: 20 * attempt,
        attempt  = lambda wildcards, attempt: attempt,
    shell:"""
        samtools stats {input.bam} > {output}
        """

rule multiqc:
    message:
        "Summarizing using MultiQC"
    input:
        expand(
            '03_star/{sample}.Aligned.sortedByCoord.out.bam.stats',
            sample=SAMPLES
        )
    output:
        "multiqc_report.html"
    conda:
        "env.yml"
    resources:
        avg_mem  = lambda wildcards, attempt: 8000 * attempt,
        mem_mb   = lambda wildcards, attempt: 8200 * attempt,
        walltime = lambda wildcards, attempt: 60 * attempt,
        attempt  = lambda wildcards, attempt: attempt,
    shell:"""
        multiqc 00_fastqc 02_fastq_trimmed 03_star
        """
