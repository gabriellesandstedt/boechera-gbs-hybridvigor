23################################################################################
## Index reference genome, align trimmed fastq files, remove reads with a map quality score < 29
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 04_filt_bam.smk
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

# assign all output files to rule all
rule all:
    input:
        expand(f"{data_dir}/{{sample}}_RG.bam", sample=sample),
        expand(f"{data_dir}/{{sample}}_RG_MD.bam", sample=sample),
        expand(f"{data_dir}/{{sample}}_RG_MD_NS.bam", sample=sample),
        expand(f"{data_dir}/{{sample}}_RG_MD_NS_FM.bam", sample=sample),
        expand(f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP.bam", sample=sample),
        expand(f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP_CS.bam", sample=sample)
        
# define rule to add or replace read groups
# picard v. 2.22: https://broadinstitute.github.io/picard/
# samtools v 1.16: https://github.com/samtools/samtools
rule add_or_replace_read_groups:
    input:
        bam=f"{data_dir}/{{sample}}_sorted.bam"
    output:
        RG_bam=f"{data_dir}/{{sample}}_RG.bam",
        RG_bai=f"{data_dir}/{{sample}}_RG.bam.bai"
    shell:
        """
        module load picard/2.22.0
        module load samtools/1.16
        echo -e "\\n["$(date)"]\\n Add read groups..\\n"
        java -jar $PICARD AddOrReplaceReadGroups \
            I={input.bam} \
            O={output.RG_bam} \
            RGID={wildcards.sample} \
            RGLB=lib_{wildcards.sample} \
            RGPL=illumina \
            RGPU=unit_{wildcards.sample} \
            RGSM={wildcards.sample}   
        samtools index {output.RG_bam}
        """

# define rule to mark and remove duplicates 
# picard v. 2.22: https://broadinstitute.github.io/picard/
# samtools v 1.16: https://github.com/samtools/samtools
rule mark_duplicates:
    input:
        RG_bam=f"{data_dir}/{{sample}}_RG.bam"
    output:
        MD_bam=f"{data_dir}/{{sample}}_RG_MD.bam",
        MD_log=f"{data_dir}/{{sample}}_RG_MDlog.txt",
        MD_bai=f"{data_dir}/{{sample}}_RG_MD.bam.bai"
    shell:
        """
        module load picard/2.22.0
        module load samtools/1.16
        echo -e "\\n["$(date)"]\\n mark and remove duplicates..\\n"
        java -jar $PICARD MarkDuplicates I={input.RG_bam} O={output.MD_bam} REMOVE_DUPLICATES=TRUE M={output.MD_log}
        samtools index {output.MD_bam}
        """

# define rule to sort bam files by name
# samtools v 1.16: https://github.com/samtools/samtools
rule namesort:
    input:
        MD_bam=f"{data_dir}/{{sample}}_RG_MD.bam"
    output:
        NS_bam=f"{data_dir}/{{sample}}_RG_MD_NS.bam"
    shell:
        """
        module load samtools/1.16
        echo -e "\\n["$(date)"]\\n sort bam by name...\\n"
        samtools sort -o {output.NS_bam} -n {input.MD_bam}
        """
        
# define rule to correct or adjust fixmate information of bam files
# samtools v 1.16: https://github.com/samtools/samtools
rule fixmate:
    input:
        NS_bam=f"{data_dir}/{{sample}}_RG_MD_NS.bam"
    output:
        FM_bam=f"{data_dir}/{{sample}}_RG_MD_NS_FM.bam"
    shell:
        """
        module load samtools/1.16
        echo -e "\\n["$(date)"]\\n run samtools fixmate...\\n"
        samtools fixmate -r {input.NS_bam} {output.FM_bam}
        """
        
# define rule to ensure that paired end reads map together
# samtools v 1.16: https://github.com/samtools/samtools
rule proper_pair:
    input:
        FM_bam=f"{data_dir}/{{sample}}_RG_MD_NS_FM.bam"
    output:
        PP_bam=f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP.bam"
    shell:
        """
        module load samtools/1.16
        echo -e "\\n["$(date)"]\\n ensure proper pairing...\\n"
        samtools view -b -f 2 -F 2048 {input.FM_bam} > {output.PP_bam}
        """    

# define rule to sort bam files by coordinates
# samtools v 1.16: https://github.com/samtools/samtools
rule coord_sort:
    input:
        PP_bam=f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP.bam"
    output:
        CS_bam=f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP_CS.bam",
        CS_bai=f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP_CS.bam.bai"
    shell:
        """
        module load samtools/1.16
        echo -e "\\n["$(date)"]\\n sort bam by coordinate...\\n"
        samtools sort -o {output.CS_bam} {input.PP_bam}
        samtools index {output.CS_bam}
        """
        
