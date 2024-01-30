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
# 156 samples, mac 16/312 alleles >0.05% 
# total snps 3862, this file can be used to perform a pca 
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
            --mac 16 \
            --max-missing 0.1 \
            --recode \
            --recode-INFO-all \
            --out {output.filtered_mac_vcf_prefix}
        
        bgzip -c {output.filtered_mac_vcf_prefix}.recode.vcf > {output.filtered_mac_vcf}
        tabix -p vcf {output.filtered_mac_vcf}
        """


rule filter_missing_invvcf:
    input:
        inv_vcf=f"{data_dir}/boech_gbs_allsamples_invariant.vcf"
    output:
        filtered_mm_inv_vcf_prefix=f"{data_dir}/boech_gbs_allsamples_invariant_filtered.vcf",
        filtered_mm_inv_vcf=f"{data_dir}/boech_gbs_allsamples_invariant_filtered.vcf.gz"
    shell:
        """
        module load vcftools/0.1.15-6
        module load htslib/1.18
        vcftools \
            --vcf {input.inv_vcf} \
            --max-missing 0.1 \
            --recode \
            --recode-INFO-all \
            --out {output.filtered_mm_invvcf_prefix}
        
        bgzip -c {output.filtered_mm_inv_vcf_prefix}.recode.vcf > {output.filtered_mm_inv_vcf}
        tabix -p vcf {output.filtered_mm_inv_vcf}
        """

rule combine_vcfs:
    input:
       gz_var_vcf=f"{data_dir}/boech_gbs_allsamples_SNPs_filtered.vcf.gz",
       gz_invar_vcf=f"{data_dir}/boech_gbs_allsamples_invariant_filtered.vcf.gz"
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

rule pixy_stats:
    input:
        vcf=f"{data_dir}/boech_gbs_allsamples_final.vcf.gz",
        pop_file=f"{data_dir}/allsamples_pop.txt"
    output:
        pixy_output=f"{data_dir}/pixy_stats.txt"
    params:
        window_size=1000000,
        n_cores=4
    shell:
        """
        ml pixy/1.2.3
        ml htslib/1.16
        pixy --stats pi fst dxy \
        --vcf {input.vcf} \
        --populations {input.pop_file} \
        --window_size {params.window_size} \
        --n_cores {params.n_cores} \
        > {output.pixy_output}
        """
