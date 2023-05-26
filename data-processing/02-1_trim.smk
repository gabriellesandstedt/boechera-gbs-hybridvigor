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
data_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS_smk/data"
scripts_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS_smk/scripts"
qc1_dir = f"{data_dir}/FASTQC1"
qc2_dir = f"{data_dir}/FASTQC2"
log_dir = f"{data_dir}/logs"

# Read the sample and adapter information from the config file into a pandas DataFrame
# config.txt has three columns (col 1| Array ID, col 2| Sample name, col 3| Adapter sequence)
df = pd.read_csv(f"{data_dir}/config.txt", sep="\t")

# Define output files for rule all
rule all:
    input:
        expand(f"{qc1_dir}/{{sample}}_fastqc.html", sample=df['Sample']),
        expand(f"{data_dir}/{{sample}}_cut.fastq.gz", sample=df['Sample']),
        expand(f"{qc2_dir}/{{sample}}_cut_fastqc.html", sample=df['Sample'])

# Define rule for checking quality of raw fastqs
# FASTQC v 0.12.1 : https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
rule fastqc1:
    input:
        fq=f"{scripts_dir}/{{sample}}"
    output:
        zip=f"{qc1_dir}/{{sample}}_fastqc.zip",
        html=f"{qc1_dir}/{{sample}}_fastqc.html"
    log:
        log=f"{log_dir}/fastqc1_{{sample}}.log"
    shell:
        """
        echo -e "\n["$(date)"]\n run FastQC on raw fastq files ...\n"
        module load fastqc/0.12.1
        fastqc -o {qc1_dir} --noextract {input.fq} &> {log.log}
        echo -e "\n["$(date)"]\n FastQC finished for {{sample}} ...\n"
        """

# Define rule to trim adapter sequences from raw fastq files
# CutAdapt v 3.5 : https://cutadapt.readthedocs.io/en/stable/
# cutadapt options:
# -a : adapters
# -m : minimum length, removes reads < 50 bps long
rule cutadapt:
    input:
        fq=f"{scripts_dir}/{{sample}}"
    output:
        cut_fq=f"{data_dir}/{{sample}}_cut.fastq.gz"
    params:
        adapter=lambda wildcards: df[df['Sample']==wildcards.sample]['Adapter'].values[0]
    shell:
        """
        echo -e "\n["$(date)"]\n Run cutadapt on fastq file {input} ...\n"
        module load cutadapt/3.5
        cutadapt -a AGATCGGAAG -a {params.adapter} -m 50 {input.fq} | gzip > {output.cut_fq}
        """

# Define rule for checking quality of trimmed fastq files
# FASTQC v 0.12.1 : https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
rule fastqc2:
    input:
        fq=f"{data_dir}/{{sample}}_cut.fastq.gz"
    output:
        zip=f"{qc2_dir}/{{sample}}_cut_fastqc.zip",
        html=f"{qc2_dir}/{{sample}}_cut_fastqc.html"
    log:
        log=f"{log_dir}/fastqc2_{{sample}}.log"
    shell:
        """
        echo -e "\n["$(date)"]\n Load and run FastQC on trimmed fq files ...\n"
        module load fastqc/0.12.1
        fastqc -o {qc2_dir} --noextract {input.fq} &> {log.log}
        echo -e "\n["$(date)"]\n FastQC finished for {{sample}} ...\n"
        """
