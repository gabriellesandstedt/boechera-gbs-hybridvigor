################################################################################
## Index reference genome, align trimmed fastq files, remove reads with a map quality score < 29
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 05_genotype.smk
################################################################################
################################################################################
import os
import pandas as pd
from snakemake.io import expand

# Define the paths to data files
data_dir = "/scratch/general/nfs1/u6048240/BOECHERA/caps_outgroup/data"
ref_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS_DEC23/ref_genome"

# reference genome -- B. stricta v 2.2 : https://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Boechera_stricta/latest_assembly_versions/GCA_018361405.1_NTU_Bstr_LTM_2.2/
ref = "GCA_018361405.1_NTU_Bstr_LTM_2.2_genomic.fa"
interval_list = "interval.list"

sample = ["SRR2070914"]

rule all:
    input:
        expand(f"{data_dir}/{{sample}}.g.vcf.gz", sample=sample)
rule hap_caller:
    input:
        ref = f"{ref_dir}/{ref}",
        CS_bam = f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP_CS.bam",
        intervals = f"{data_dir}/{interval_list}"
    output:
        gvcf = f"{data_dir}/{{sample}}.g.vcf.gz"
    shell:
        """
        module load gatk/4.1
        gatk HaplotypeCaller \
            -I {input.CS_bam} \
            -O {output.gvcf} \
            -R {input.ref} \
            -L {input.intervals} \
            -ERC GVCF
        """
