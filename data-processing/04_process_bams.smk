################################################################################
## Add readgroups using picardtools
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 04_process_bams.smk
################################################################################
################################################################################
import os
import pandas as pd
from snakemake.io import expand

# Define the paths to data files
data_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS_DEC23/data"

# Read the sample information from the config file into a pandas DataFrame
# config2.txt has two columns (col 1| Array ID, col 2| Sample name)
df = pd.read_csv(f"{data_dir}/config2.txt", sep="\t")

# define output files for rule all
rule all:
    input:
        expand(f"{data_dir}/{{sample}}_sorted_RG.bam", sample=df['Sample']),
        expand(f"{data_dir}/{{sample}}_sorted_RG.bam.bai", sample=df['Sample'])

# define rule to add or replace read groups
# picard v 2.22: https://broadinstitute.github.io/picard/
rule add_or_replace_read_groups:
    input:
        bam=f"{data_dir}/{{sample}}_sorted.bam"
    output:
        RG_bam=f"{data_dir}/{{sample}}_sorted_RG.bam"
    shell:
        """
        module load picard/2.22.0
        echo -e "\\n["$(date)"]\\n Add read groups..\\n"
        java -jar $PICARD AddOrReplaceReadGroups \
            I={input.bam} \
            O={output.RG_bam} \
            RGID={wildcards.sample} \
            RGLB=lib_{wildcards.sample} \
            RGPL=illumina \
            RGPU=unit_{wildcards.sample} \
            RGSM={wildcards.sample}
        """

# define rule to index bam out_files
# Samtools v 1.16 : https://github.com/samtools/samtools
rule index_RG_bam:
    input:
        RG_bam=f"{data_dir}/{{sample}}_sorted_RG.bam"
    output:
        RG_bai=f"{data_dir}/{{sample}}_sorted_RG.bam.bai"
    shell:
        """
        module load samtools/1.16
        echo -e "\\n["$(date)"]\\n Index BAM file...\\n"
        samtools index {input.RG_bam} {output.RG_bai}
        """

# define rule to assess details on bam files
# Samtools v 1.16 : https://github.com/samtools/samtools
rule assess_quality:
    input:
        sorted_RG_bai=f"{data_dir}/{{sample}}_sorted_RG.bam"
    shell:
        """
        module load samtools/1.16
        echo -e "\\n["$(date)"]\\ run samtools flagstat on final bam file...\\n"
        samtools flagstat {input.sorted_RG_bai}
        """
