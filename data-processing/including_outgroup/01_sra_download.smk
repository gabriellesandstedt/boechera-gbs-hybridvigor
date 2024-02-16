################################################################################
## download capsella grandiflora dna-sequencing
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## snakemake version: snakemake/7.22.0-foss-2022a
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 01_sra_download.smk
################################################################################
################################################################################
import os

# assign data directory
data_dir = "/scratch/general/nfs1/u6048240/BOECHERA/caps_outgroup/data"

# assign capsella grandiflora sample to be downloaded from NCBI
# SRR2070914

sample = ["SRR2070914"]

# define all output files
rule all:
    input:
        expand(f"{data_dir}/{{sample}}_1.fastq.gz", sample=sample),
        expand(f"{data_dir}/{{sample}}_2.fastq.gz", sample=sample)

# define rule to download fastqs from NCBI
# https://www.ncbi.nlm.nih.gov/sra
rule download_fastq:
    output:
        fq1=f"{data_dir}/{{sample}}_1.fastq.gz",
        fq2=f"{data_dir}/{{sample}}_2.fastq.gz"
    params:
        sra=lambda wildcards: wildcards.sample
    shell:
        """
        module load sra-toolkit/3.0.10
        echo -e "\\n["$(date)"]\\n Download {params.sra} from SRA database \\n"
        prefetch {params.sra} && fastq-dump --split-3 --gzip {params.sra} -O {data_dir}
        sleep 60
        """
