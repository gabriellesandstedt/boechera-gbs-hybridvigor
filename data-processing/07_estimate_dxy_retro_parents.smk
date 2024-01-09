################################################################################
## estimate Dxy for B. retrofracta parents 
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

# subset samples for dxy among retrofracta parents
rule retro_parent_dxy_samples:
    input:
        vcf=f"{data_dir}/boech_gbs_allsamples_combined_final.vcf.gz"
    output:
        retro_parents_samples_vcf=f"{data_dir}/boech_gbs_final_retro_parents_samples.vcf"
    params:
        samples_to_include="pa1_pool,pa10_pool,pa2_pool,pa3_pool,pa4_pool,pa5_pool,pa6_pool,pa7_pool,pa8_pool,pa9_pool,pw1_pool,pw10_pool,pw11_pool,pw2_pool,pw3_pool,pw4_pool,pw5_pool,pw6b_pool,pw7_pool,pw8b_pool,pw9_pool"
    shell:
        """
        module load bcftools/1.16
        echo -e '\\n['$(date)']\\n subset vcf ..\\n'
        bcftools view -s "{params.samples_to_include}" -o {output.retro_parents_samples_vcf} {input.vcf}
        """
        
# select invariant sites from combined vcf 
rule select_invariant_sites:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        vcf=f"{data_dir}/boech_gbs_final_retro_parents_samples.vcf"
    output:
        invariant_vcf=f"{data_dir}/boech_gbs_retro_parents_invariant.vcf"
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
        vcf=f"{data_dir}/boech_gbs_final_retro_parents_samples.vcf"
    output:
        biallelic_vcf=f"{data_dir}/boech_gbs_retro_parents_samples_SNPs.vcf"
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
# 21 samples, mac 3
# 3/42, > 0.05
rule filter_minor_allele_count:
    input:
        snps_vcf=f"{data_dir}/boech_gbs_retro_parents_samples_SNPs.vcf"
    output:
        filtered_mac_vcf_prefix=f"{data_dir}/boech_gbs_retro_parents_samples_SNPs_filtered.vcf",
        filtered_mac_vcf=f"{data_dir}/boech_gbs_retro_parents_samples_SNPs_filtered.vcf.gz"
    shell:
        """
        module load vcftools/0.1.15-6
        module load htslib/1.18
        vcftools \
            --vcf {input.snps_vcf} \
            --remove-indels \
            --min-alleles 2 \
            --max-alleles 2 \
            --mac 3 \
            --recode \
            --recode-INFO-all \
            --out {output.filtered_mac_vcf_prefix}
        
        bgzip -c {output.filtered_mac_vcf_prefix}.recode.vcf > {output.filtered_mac_vcf}
        tabix -p vcf {output.filtered_mac_vcf}
        """

# bgzip and tabix files
rule vcf_to_gzvcf:
    input:
        inv_vcf=f"{data_dir}/boech_gbs_retro_parents_invariant.vcf"
    output:
        gz_invar_vcf=f"{data_dir}/boech_gbs_retro_parents_invariant.vcf.gz",
        tabix_invar_vcf=f"{data_dir}/boech_gbs_retro_parents_invariant.vcf.gz.tbi"
    shell:
        """
        module load htslib/1.16
        bgzip {input.inv_vcf}
        tabix -p vcf {output.gz_invar_vcf}
        """
      
# combine SNPs and invariant sites 
rule combine_vcfs:
    input:
       gz_var_vcf=f"{data_dir}/boech_gbs_retro_parents_samples_SNPs_filtered.vcf.gz",
       gz_invar_vcf=f"{data_dir}/boech_gbs_retro_parents_invariant.vcf.gz"
    output:
       final_vcf=f"{data_dir}/boech_gbs_retro_parents_final.vcf.gz",
       tabix_final=f"{data_dir}/boech_gbs_retro_parents_final.vcf.gz.tbi"
    shell:
        """
        module load htslib/1.16
        module load bcftools/1.16
        bcftools concat {input.gz_var_vcf} {input.gz_invar_vcf} -a -Oz -o {output.final_vcf}
        tabix -p vcf {output.final_vcf}
        """   


rule pixy_stats:
    input:
        vcf=f"{data_dir}/boech_gbs_retro_parents_final.vcf.gz",
        pop_file=f"{data_dir}/retro_pop.txt"
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
