################################################################################
## compare all GBS samples with PCA, pi, Dxy
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete --latency-wait 60 --cores 4 -s 07_estimate_dxy_retro_parents.smk
################################################################################
################################################################################
# define directories
data_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS_May23/data"
ref_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS/ref_genomes/data"

# define the reference genome filename
ref_genome = "GCA_018361405.1_NTU_Bstr_LTM_2.2_genomic.fa"

# define rule all statement
rule all:
    input:
        f"{data_dir}/boech_gbs_final_retro_parents_samples.vcf",
        f"{data_dir}/boech_gbs_retro_parents_invariant.vcf",
        f"{data_dir}/boech_gbs_retro_parents_samples_SNPs.vcf",
        f"{data_dir}/boech_gbs_retro_parents_samples_SNPs_filtered.vcf.gz",
        f"{data_dir}/boech_gbs_retro_parents_invariant.vcf.gz",
        f"{data_dir}/boech_gbs_retro_parents_final.vcf.gz",
        f"{data_dir}/pixy_stats.txt"

# select invariant sites from vcf 
rule select_invariant_sites:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        vcf=f"{data_dir}/boech_gbs_allsamples_combined_final.vcf.gz"
    output:
        invariant_vcf=f"{data_dir}/boech_gbs_allsamples_invariant.vcf"
    shell:
        """
        module load gatk/4.1
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            --select-type-to-include NO_VARIATION \
            -O {output.invariant_vcf}
        """

# select biallelic SNPs from combined vcf 
rule select_biallelic_snps:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        vcf=f"{data_dir}/boech_gbs_allsamples_combined_final.vcf.gz"
    output:
        biallelic_vcf=f"{data_dir}/boech_gbs_allsamples_SNPs.vcf"
    shell:
        """
        module load gatk/4.1
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            --select-type-to-include SNP \
            --restrict-alleles-to BIALLELIC \
            -O {output.biallelic_vcf}
        """

# filter SNPs to remove rare alleles
# 156 samples, mac 8/156 samples >0.05% 
# total snps 5316, this file can be used to perform a pca 
rule filter_minor_allele_count:
    input:
        snps_vcf=f"{data_dir}/boech_gbs_allsamples_SNPs.vcf"
    output:
        filtered_mac_vcf_prefix=f"{data_dir}/boech_gbs_allsamples_SNPs_filtered.vcf",
        filtered_mac_vcf=f"{data_dir}/boech_gbs_allsamples_SNPs_filtered.vcf.gz"
    shell:
        """
        module load vcftools/0.1.15-6
        module load htslib/1.18
        vcftools \
            --vcf {input.snps_vcf} \
            --remove-indels \
            --min-alleles 2 \
            --max-alleles 2 \
            --mac 8 \
            --recode \
            --recode-INFO-all \
            --out {output.filtered_mac_vcf_prefix}
        
        bgzip -c {output.filtered_mac_vcf_prefix}.recode.vcf > {output.filtered_mac_vcf}
        tabix -p vcf {output.filtered_mac_vcf}
        """


# bgzip and tabix files
rule vcf_to_gzvcf:
    input:
        inv_vcf=f"{data_dir}/boech_gbs_allsamples_invariant.vcf"
    output:
        gz_invar_vcf=f"{data_dir}/boech_gbs_allsamples_invariant.vcf.gz",
        tabix_invar_vcf=f"{data_dir}/boech_gbs_allsamples_invariant.vcf.gz.tbi"
    shell:
        """
        module load htslib/1.16
        bgzip {input.inv_vcf}
        tabix -p vcf {output.gz_invar_vcf}
        """

rule combine_vcfs:
    input:
       gz_var_vcf=f"{data_dir}/boech_gbs_allsamples_SNPs_filtered.vcf.gz",
       gz_invar_vcf=f"{data_dir}/boech_gbs_allsamples_invariant.vcf.gz"
    output:
       final_vcf=f"{data_dir}/boech_gbs_allsamples_final.vcf.gz",
       tabix_final=f"{data_dir}/boech_gbs_allsamples_final.vcf.gz.tbi"
    shell:
        """
        module load htslib/1.16
        module load bcftools/1.16
        bcftools concat {input.gz_var_vcf} {input.gz_invar_vcf} -a -Oz -o {output.final_vcf}
        tabix -p vcf {output.final_vcf}
        """  
