################################################################################
## downsample vcf for entropy
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s subset_samples_vcf.smk
################################################################################
################################################################################

# Define the paths to data files
data_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS_May23/data"
vcf = "allsamples_allsites_final_snps.vcf"

# define output files for rule all
rule all:
    input:
        f"{data_dir}/allsamples_allsites_final_snps_subset_noSRF2s.vcf",
        f"{data_dir}/allsamples_allsites_final_snps_subset_onlyRetros.vcf"

# define rule to subset vcf to remove SR lab generated F2s
# bcftools: 
rule noSRF2s:
    input:
        bam=f"{data_dir}/{{sample}}_sorted.bam"
    output:
        RG_bam=f"{data_dir}/{{sample}}_sorted_RG.bam"
    shell:
        """
        module bcftools/1.16
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
        
bcftools/1.16


samples_to_remove="sample1,sample2"

# Create a new VCF file without the specified samples
bcftools view -s "^$samples_to_remove" -o output.vcf input.vcf
