################################################################################
## Combine asex fastqs
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## To run snakemake script: snakemake --rerun-incomplete --latency-wait 60 --cores 4 -s 02-1_combinefastqs.smk
################################################################################
################################################################################
import pandas as pd
from snakemake.io import expand

# Define the paths to the raw data and the output directories
data_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS_DEC23/data"
cut_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS_DEC23/data"

# Read the sample and adapter information from the config file into a pandas DataFrame
# config_retro.txt has two columns (col 1| Array ID, col 2| Sample name without suffix (BrA, BrB))
# config_retro_stricta.txt has two columns (col 1 | Array ID, column 2| Sample name without suffix (BsA, BsB))
df1 = pd.read_csv(f"{data_dir}/config_retro.txt", sep="\t")
df2 = pd.read_csv(f"{data_dir}/config_retro_stricta.txt", sep="\t")

# define output files for rule all
rule all:
    input:
        expand(f"{data_dir}/{{sample_combine1}}_RR_cut.fastq.gz", sample_combine1=df1['Sample']),
        expand(f"{data_dir}/{{sample_combine2}}_SR_cut.fastq.gz", sample_combine2=df2['Sample'])

# combine fastq files that are bioreplicates of retroxretro asexuals
rule combine_retro_asex:
    input:
        fqA=f"{cut_dir}/{{sample_combine1}}_BrA_cut.fastq.gz",
        fqB=f"{cut_dir}/{{sample_combine1}}_BrB_cut.fastq.gz"
    output:
        combined_fq_RR=f"{data_dir}/{{sample_combine1}}_RR_cut.fastq.gz"
    shell:
        """
        zcat {input.fqA} {input.fqB} | gzip -c > {output.combined_fq_RR}
        """

# combine fastq files that are bioreplicates of strictaxretro asexuals
rule combine_stricta_retro_asex:
    input:
        fqA=f"{cut_dir}/{{sample_combine2}}_BsA_cut.fastq.gz",
        fqB=f"{cut_dir}/{{sample_combine2}}_BsB_cut.fastq.gz"
    output:
        combined_fq_SR=f"{data_dir}/{{sample_combine2}}_SR_cut.fastq.gz"
    shell:
        """
        zcat {input.fqA} {input.fqB} | gzip -c > {output.combined_fq_SR}
        """
