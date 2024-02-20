23################################################################################
## Index reference genome, align trimmed fastq files, remove reads with a map quality score < 29
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 03_align.smk
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

sample = ["SRR2070914"]

# define all output files from script
rule all:
    input:
        expand(f"{ref_dir}/{ref}.{{suffix}}", suffix=["bwt", "amb", "ann", "pac"]),
        expand(f"{data_dir}/{{sample}}.sam", sample=sample),
        expand(f"{data_dir}/{{sample}}.bam", sample=sample),
        expand(f"{data_dir}/{{sample}}_sorted.bam", sample=sample),
        expand(f"{data_dir}/{{sample}}_sorted.bam.bai", sample=sample)

# define rule to index reference genome using bwa index
# BWA v 2020_03_19 : https://bio-bwa.sourceforge.net
rule bwa_index:
    input:
        ref = f"{ref_dir}/{ref}"
    output:
        bwt = "{ref}.bwt",
        amb = "{ref}.amb",
        ann = "{ref}.ann",
        pac = "{ref}.pac"
    shell:
        """
        module load bwa/2020_03_19
        bwa index -a bwtsw {input.ref}
        echo -e "\n["$(date)"]\n done ...\n"
        """

# define rule to align fastqs to reference genome
# BWA v 2020_03_19 : https://bio-bwa.sourceforge.net
rule bwa_mem:
    input:
        ref = f"{ref_dir}/{ref}",
        cut_fq1 = f"{data_dir}/{{sample}}_1_cut.fastq.gz",
        cut_fq2 = f"{data_dir}/{{sample}}_2_cut.fastq.gz"
    output:
        sam = f"{data_dir}/{{sample}}.sam"
    shell:
        """
        module load bwa/2020_03_19
        echo -e "\\n["$(date)"]\\n run BWA mem..\\n"
        bwa mem {input.ref} {input.cut_fq1} {input.cut_fq2} > {output.sam}
        """

# define rule to convert raw alignments to bam files and removing reads < map quality 29
# Samtools v 1.16 : https://github.com/samtools/samtools
rule sam_bam_q29:
    input:
        sam = f"{data_dir}/{{sample}}.sam"
    output:
        bam = f"{data_dir}/{{sample}}.bam"
    shell:
        """
        module load samtools/1.16
        samtools view -q 29 -b {input.sam} > {output.bam}
        echo -e "\\n["$(date)"]\\n run sam2bamQ29..\\n"
        """

# define rule to sort bam file
# Samtools v 1.16 : https://github.com/samtools/samtools
rule sort_bam:
    input:
        bam = f"{data_dir}/{{sample}}.bam"
    output:
        sorted_bam = f"{data_dir}/{{sample}}_sorted.bam"
    shell:
        """
        module load samtools/1.16
        samtools sort {input.bam} -o {output.sorted_bam}
        echo -e "\\n["$(date)"]\\n bam file is sorted ..\\n"
        """

# define rule to index bam file
# Samtools v 1.16 : https://github.com/samtools/samtools
rule index_bam:
    input:
        sorted_bam=f"{data_dir}/{{sample}}_sorted.bam"
    output:
        index_bam=f"{data_dir}/{{sample}}_sorted.bam.bai"
    shell:
        """
        module load samtools/1.16
        samtools index {input.sorted_bam}
        echo -e "\\n["$(date)"]\\n bam file is indexed ..\\n"
        """
