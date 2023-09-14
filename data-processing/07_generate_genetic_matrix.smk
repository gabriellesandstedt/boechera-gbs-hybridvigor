################################################################################
## generate genetic matrix
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete --latency-wait 60 --cores 4 -s 06_filter_allsites_allsamples_vcf.smk
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
        f"{data_dir}/boech_gbs_final_genetic_matrix_samples.vcf",
        f"{data_dir}/boech_gbs_final_genetic_matrix_samples_SNPs.vcf",
        f"{data_dir}/boech_gbs_final_genetic_matrix_samples_SNPs_filtered.vcf.gz",
        f"{data_dir}/boech_gbs_final_genetic_matrix_samples_SNPs_filtered.vcf.rel"

# subset samples for genetic matrix
rule genetic_matrix_samples:
    input:
        vcf=f"{data_dir}/boech_gbs_allsamples_combined_final.vcf.gz"
    output:
        gen_matrix_samples_vcf=f"{data_dir}/boech_gbs_final_genetic_matrix_samples.vcf"
    params:
        samples_to_include="pa1_pool,pa10_pool,pa2_pool,pa3_pool,pa4_pool,pa5_pool,pa6_pool,pa7_pool,pa8_pool,pa9_pool,pw1_pool,pw10_pool,pw11_pool,pw2_pool,pw3_pool,pw4_pool,pw5_pool,pw6b_pool,pw7_pool,pw8b_pool,pw9_pool,am1_1,am10_1,am11_1,am12_1,am13_1,am14_1,am15_1,am2_1,am3_1,am4_1,am5_1,am6_1,am7_1,am8_1,am9_1,wi1_1,wi2_1,wi4_1,wi7_1,wi8_1,wi9_1"
    shell:
        """
        module load bcftools/1.16
        echo -e '\\n['$(date)']\\n subset vcf ..\\n'
        bcftools view -s "{params.samples_to_include}" -o {output.gen_matrix_samples_vcf} {input.vcf}
        """

# select biallelic SNPs
rule select_biallelic_snps:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        vcf=f"{data_dir}/boech_gbs_final_genetic_matrix_samples.vcf"
    output:
        biallelic_vcf=f"{data_dir}/boech_gbs_final_genetic_matrix_samples_SNPs.vcf"
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

# filter biallelic SNPs for minor allele count
rule filter_minor_allele_count:
    input:
        filtered_hets_vcf=f"{data_dir}/boech_gbs_final_genetic_matrix_samples_SNPs.vcf"
    output:
        filtered_mac_vcf_prefix=f"{data_dir}/boech_gbs_final_genetic_matrix_samples_SNPs_filtered.vcf",
        filtered_mac_vcf=f"{data_dir}/boech_gbs_final_genetic_matrix_samples_SNPs_filtered.vcf.gz"
    shell:
        """
        module load vcftools/0.1.15-6
        module load htslib/1.18
        vcftools \
            --vcf {input.filtered_hets_vcf} \
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

# define rule to convert snps vcf to a bed file         
rule gzvcf_to_bed:
    input:
        final_gzvcf=f"{data_dir}/boech_gbs_final_genetic_matrix_samples_SNPs_filtered.vcf.gz"
    output:
        bed=f"{data_dir}/boech_gbs_final_genetic_matrix_samples_SNPs_filtered.vcf"
    shell:
        """
        module load plink/2.0
        plink2 --vcf {input.final_gzvcf} \
               --set-missing-var-ids @:#[b37] \
               --make-bed \
               --out {output.bed} \
               --freq \
               --allow-extra-chr
        """

# define rule to create a square genetic relatedness matrix 
rule calculate_relatedness:
    input:
        freq=f"{data_dir}/boech_gbs_final_genetic_matrix_samples_SNPs_filtered.vcf.afreq",
    output:
        rel=f"{data_dir}/boech_gbs_final_genetic_matrix_samples_SNPs_filtered.vcf.rel"
    params:
        prefix = "boech_gbs_final_genetic_matrix_samples_SNPs_filtered.vcf"
    shell:
        """
        module load plink/2.0
        plink2 --read-freq {input.freq} \
               --bfile {params.prefix} \
               --make-rel square \
               --allow-extra-chr \
               --out {output.rel}
        """ 
