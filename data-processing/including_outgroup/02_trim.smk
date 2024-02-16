################################################################################
## Run FASTQC and cutadapt on fastq files
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## To run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 02_trim.smk
################################################################################
################################################################################
import os
import pandas as pd
from snakemake.io import expand

# Define the paths to the raw data and the output directories
data_dir = "/scratch/general/nfs1/u6048240/BOECHERA/caps_outgroup/data"
qc1_dir = f"{data_dir}/FASTQC1"
qc2_dir = f"{data_dir}/FASTQC2"
log_dir = f"{data_dir}/logs"

sample = ["SRR2070914"]

# Define output files for rule all
rule all:
    input:
        expand(f"{qc1_dir}/{{sample}}_1_fastqc.html", sample=samples),
        expand(f"{qc1_dir}/{{sample}}_2_fastqc.html", sample=samples),
        expand(f"{data_dir}/{{sample}}_1_cut.fq.gz", sample=samples),
        expand(f"{data_dir}/{{sample}}_2_cut.fq.gz", sample=samples),
        expand(f"{qc2_dir}/{{sample}}_1_cut_fastqc.html", sample=samples),
        expand(f"{qc2_dir}/{{sample}}_2_cut_fastqc.html", sample=samples)

# Define rule for checking quality of raw fastqs
# FASTQC v 0.12.1 : https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
rule fastqc_raw:
    input:
        fq1=f"{data_dir}/{{sample}}_1.fastq.gz",
        fq2=f"{data_dir}/{{sample}}_2.fastq.gz"
    output:
        fqc1=f"{qc1_dir}/{{sample}}_1_fastqc.html",
        fqc2=f"{qc1_dir}/{{sample}}_2_fastqc.html"
    log:
        log1=f"{log_dir}/{{sample}}_1.log",
        log2=f"{log_dir}/{{sample}}_2.log"
    shell:
        """
        module load FastQC/0.11.9-Java-11
        echo -e "\\n["$(date)"]\\n Run FastQC on fastq file {input.fq1} ...\\n"
        fastqc -o {qc1_dir} --noextract {input.fq1} &> {log.log1} 
        echo -e "\\n["$(date)"]\\n FastQC round 1, read 1 finished ...\\n"

        echo -e "\\n["$(date)"]\\n Run FastQC on fastq file {input.fq2} ...\\n"
        fastqc -o {qc1_dir} --noextract {input.fq2} &> {log.log2}
        echo -e "\\n["$(date)"]\\n FastQC round 1, read 2 finished ...\\n"
        """

# Define rule to trim adapter sequences from raw fastq files
# CutAdapt v 3.5 : https://cutadapt.readthedocs.io/en/stable/
# cutadapt options:
# -a : adapters
# -m : minimum length, removes reads < 50 bps long
rule cutadapt:
    input:
        fq1=f"{data_dir}/{{sample}}_1.fastq.gz"
        fq2=f"{data_dir}/{{sample}}_2.fastq.gz"
    output:
        cut_fq1=f"{data_dir}/{{sample}}_1_cut.fastq.gz"
        cut_fq2=f"{data_dir}/{{sample}}_2_cut.fastq.gz"
    shell:
        """
        echo -e "\n["$(date)"]\n Run cutadapt on fastq file {input} ...\n"
        module load cutadapt/3.5
        cutadapt -a TACACTCTTTCCCTACACGACGCTCTTCCGATCT -a GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -m 50 {input.fq1} | gzip > {output.cut_fq1}
        cutadapt -a TACACTCTTTCCCTACACGACGCTCTTCCGATCT -a GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -m 50 {input.fq2} | gzip > {output.cut_fq2}
        """

# Define rule for checking quality of trimmed fastq files
# FASTQC v 0.12.1 : https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
rule fastqc_trimmed:
    input:
        cut_fq1=f"{data_dir}/{{sample}}_1_cut.fastq.gz",
        cut_fq2=f"{data_dir}/{{sample}}_2_cut.fastq.gz"
    output:
        cut_fqc1=f"{qc2_dir}/{{sample}}_1_cut_fastqc.html",
        cut_fqc2=f"{qc2_dir}/{{sample}}_2_cut_fastqc.html"
    log:
        log1=f"{log_dir}/trim_{{sample}}_1.log",
        log2=f"{log_dir}/trim_{{sample}}_2.log"
    shell:
        """
        module load FastQC/0.11.9-Java-11
        echo -e "\\n["$(date)"]\\n Run FastQC on fastq file {input.trim_fq1} ...\\n"
        fastqc -o {qc2_dir} --noextract {input.cut_fq1} &> {log.log1}
        echo -e "\\n["$(date)"]\\n FastQC round 2, read 1 finished ...\\n"

        echo -e "\\n["$(date)"]\\n Run FastQC on fastq file {input.trim_fq2} ...\\n"
        fastqc -o {qc2_dir} --noextract {input.cut_fq2} &> {log.log2}
        echo -e "\\n["$(date)"]\\n FastQC round 2, read 2 finished ...\\n"
        """
