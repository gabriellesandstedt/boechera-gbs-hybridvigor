################################################################################
## joint genotype samples using GATK
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete --latency-wait 60 --cores 4 -s 05_jointgenotyping.smk
################################################################################
################################################################################
import os
import pandas as pd
from snakemake.io import directory
from snakemake.io import expand

# Define the paths to data files
data_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS_DEC23/data"
ref_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS_DEC23/ref_genome"

# Read the sample information from the config file into a pandas DataFrame
# config2.txt has two columns (col 1| Array ID, col 2| Sample name)
# txt files are available in /boechera-gbs-hybridvigor/helper-files
df = pd.read_csv(f"{data_dir}/config2.txt", sep="\t")

# define sample map for all samples (sample_map1) and matrix (sample_map2), this includes two columns: column 1 with sample name and column 2 with path to samples
sample_map1 = "sample_map1.txt"

# define interval list with list of chromosomes
interval_list = "interval.list"

# Define the reference genome filename
ref = "GCA_018361405.1_NTU_Bstr_LTM_2.2_genomic.fa"

# define all output files for rule all 
rule all:
    input:
        expand(f"{ref_dir}/{ref}.dict"),
        expand(f"{data_dir}/{{sample}}_sorted_RG.bam", sample=df['Sample']),
        expand(f"{data_dir}/{{sample}}.g.vcf", sample=df['Sample']),
        "/scratch/general/nfs1/u6048240/BOECHERA/GBS_DEC23/data/DB_allsamples",
        expand(f"{data_dir}/boech_gbs_allsamples.vcf"),

# define rule to index reference with GATK 
# is there an fai file for the reference? if not, need to run : samtools faidx ref.fasta
rule index_reference:
    input:
        ref = f"{ref_dir}/{ref}"
    output:
        index = f"{ref_dir}/{ref}.dict"
    shell:
        """
        module load gatk/4.1
        gatk CreateSequenceDictionary \
            -R GCA_018361405.1_NTU_Bstr_LTM_2.2_genomic.fa \
            -O GCA_018361405.1_NTU_Bstr_LTM_2.2_genomic.dict
        """

# define rule to call potential variants for each sample
rule hap_caller:
    input:
        ref = f"{ref_dir}/{ref}",
        bam = f"{data_dir}/{{sample}}_sorted_RG.bam",
        intervals = f"{data_dir}/{interval_list}"
    output:
        gvcf = f"{data_dir}/{{sample}}.g.vcf"
    shell:
        """
        module load gatk/4.1
        gatk HaplotypeCaller \
            -I {input.bam} \
            -O {output.gvcf} \
            -R {input.ref} \
            -L {input.intervals} \
            -ERC GVCF
        """

# define rule to combine gvcfs of all samples in a database
rule genomicsdb_import_allsamples:
    input:
        gvcf=expand(f"{data_dir}/{{sample}}.g.vcf", sample=df['Sample']),
        interval_list=f"{data_dir}/{interval_list}",
        map_allsamples=f"{data_dir}/{sample_map1}"
    output:
        database=directory("/scratch/general/nfs1/u6048240/BOECHERA/GBS_DEC23/data/DB_allsamples")
    shell:
        """
        module load gatk/4.1
        gatk GenomicsDBImport \
            --genomicsdb-workspace-path {output.database} \
            --batch-size 156 \
            --reader-threads 6 \
            --sample-name-map {input.map_allsamples} \
            --intervals {input.interval_list}
        """

# define rule to joint genotype for all samples at all sites
rule joint_genotype_allsamples:
    input:
        ref=f"{ref_dir}/{ref}",
        intervals=f"{data_dir}/{interval_list}"
    output:
        boech_output=f"{data_dir}/boech_gbs_allsamples.vcf"
    params:
        genomicsdb="gendb://DB_allsamples"
    shell:
        """
        module load gatk/4.1
        gatk GenotypeGVCFs \
            -R {input.ref} \
            -V {params.genomicsdb} \
            -L {input.intervals} \
            --allow-old-rms-mapping-quality-annotation-data \
            --all-sites \
            -O {output.boech_output}
        """
