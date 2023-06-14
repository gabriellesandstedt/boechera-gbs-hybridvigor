################################################################################
## filter vcf for the genetic relatedness matrix 
## most filtering methods were modified from: https://evodify.com/gatk-in-non-model-organism/
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete --latency-wait 60 --cores 4 -s 06_filter_genetic_matrix_vcf.smk
################################################################################
################################################################################
# define directories
data_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS_May23/data"
scripts_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS_May23/scripts"
ref_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS/ref_genomes/data"

# define the reference genome filename
ref_genome = "GCA_018361405.1_NTU_Bstr_LTM_2.2_genomic.fa"

# define rule all statement
rule all:
    input:
        f"{data_dir}/boech_gbs_matrix_biallelic_filter_DP_hets_mac.vcf.recode.vcf.gz",
        f"{data_dir}/boech_gbs_matrix_biallelic_filter_DP_hets_mac.vcf.recode.vcf.gz.tbi",
        f"{data_dir}/boech_gbs_matrix_biallelic_filter_DP_hets_mac.vcf.recode.vcf.bed",
        f"{data_dir}/boech_gbs_matrix_biallelic_filter_DP_hets_mac.vcf.recode.vcf.rel"


# select biallelic SNPs
rule select_biallelic_snps:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        vcf=f"{data_dir}/boechera_gbs_matrix.vcf"
    output:
        biallelic_vcf=f"{data_dir}/boech_gbs_matrix_biallelic.vcf"
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

# create variant table
rule variant_table:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        vcf=f"{data_dir}/boechera_gbs_matrix.vcf",
        rscript=f"{scripts_dir}/filtering_diagnostics.R"
    output:
        table=f"{data_dir}/boech_gbs_matrix_variant.table"
    shell:
        """
        module load gatk/4.1
        module load R
        gatk VariantsToTable \
            -R {input.ref} \
            -V {input.vcf} \
            -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
            -O {output.table}
        
        Rscript {input.rscript}
        """

# rule to filter variants:
rule filter_variants:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        biallelic_vcf=f"{data_dir}/boech_gbs_matrix_biallelic.vcf"
    output:
        filtered_vcf=f"{data_dir}/boech_gbs_matrix_biallelic_filter.vcf"
    shell:
        """
        module load gatk/4.4.0.0
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.biallelic_vcf} \
            --filter-expression "QD < 2.0" --filter-name "QD2" \
            --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
            --filter-expression "SOR > 5.0" --filter-name "SOR5" \
            --filter-expression "FS > 60.0" --filter-name "FS60" \
            --filter-expression "MQ < 40.0" --filter-name "MQ40" \
            --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -O {output.filtered_vcf}
        """

# extract passed variants
rule extract_passed_variants:
    input:
        filtered_vcf=f"{data_dir}/boech_gbs_matrix_biallelic_filter.vcf"
    output:
        filtered_passed_vcf=f"{data_dir}/boech_gbs_matrix_biallelic_filterPASSED.vcf"
    shell:
        """
        grep -E '^#|PASS' {input.filtered_vcf} > {output.filtered_passed_vcf}
        """

# extract depth information from whole-genome VCF
rule table_for_depth:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        vcf=f"{data_dir}/boechera_gbs_matrix.vcf",
        rscript=f"{scripts_dir}/filtering_diagnostics_DP.R"
    output:
        dp_table=f"{data_dir}/boechera_gbs_matrix.DP.table"
    shell:
        """
        module load gatk/4.1
        module load R
        gatk VariantsToTable \
            -R {input.ref} \
            -V {input.vcf} \
            -F CHROM -F POS -GF GT -GF DP \
            -O {output.dp_table}
            
        Rscript {input.rscript}
        """

# filter variants based on DP
rule filter_variants_DP:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        filtered_passed_vcf=f"{data_dir}/boech_gbs_matrix_biallelic_filterPASSED.vcf"
    output:
        filtered_DP_vcf=f"{data_dir}/boech_gbs_matrix_biallelic_filter_DP.vcf"
    shell:
        """
        module load gatk/4.1
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.filtered_passed_vcf} \
            -G-filter "DP < 5 || DP > 140" \
            -G-filter-name "DP_5-140" \
            -O {output.filtered_DP_vcf}
        """

# select variants and set filtered genotypes to no-call
rule select_variants:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        filtered_DP_vcf=f"{data_dir}/boech_gbs_matrix_biallelic_filter_DP.vcf"
    output:
        filtered_noCall_vcf=f"{data_dir}/boech_gbs_matrix_biallelic_filter_DPfilterNoCall.vcf"
    shell:
        """
        module load gatk/4.1
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.filtered_DP_vcf} \
            --set-filtered-gt-to-nocall TRUE \
            -O {output.filtered_noCall_vcf}
        """

# filter heterozygous genotypes
rule filter_heterozygous_genotypes:
    input:
        rscript=f"{scripts_dir}/filter_heterozygous_genotypes.R"
    output:
        f"{data_dir}/boech_gbs_matrix_biallelic_filter_DP_hets.vcf"
    shell:
        """
        module load R
        Rscript {input.rscript}
        """

# filter minor allele count
# mac 3/41 samples, >0.05% 
# final vcf contains 2499 snp positions
rule filter_minor_allele_count:
    input:
        filtered_hets_vcf=f"{data_dir}/boech_gbs_matrix_biallelic_filter_DP_hets.vcf"
    output:
        final_vcf=f"{data_dir}/boech_gbs_matrix_biallelic_filter_DP_hets_mac.vcf"
    shell:
        """
        module load vcftools/0.1.15-6
        vcftools \
            --vcf {input.filtered_hets_vcf} \
            --remove-indels \
            --min-alleles 2 \
            --max-alleles 2 \
            --mac 3 \
            --recode \
            --recode-INFO-all \
            --out {output.final_vcf}
        """

# define rule to bgzip vcf file
rule vcf_to_gzvcf:
    input:
        final_vcf=f"{data_dir}/boech_gbs_matrix_biallelic_filter_DP_hets_mac.vcf.recode.vcf"
    output:
        final_gzvcf=f"{data_dir}/boech_gbs_matrix_biallelic_filter_DP_hets_mac.vcf.recode.vcf.gz",
        tabix_gzvcf=f"{data_dir}/boech_gbs_matrix_biallelic_filter_DP_hets_mac.vcf.recode.vcf.gz.tbi"
    shell:
        """
        module load htslib/1.16
        bgzip {input.final_vcf}
        tabix -p vcf {output.final_gzvcf}
        """

# define rule to convert a vcf to a bed file         
rule vcfgz_to_bed:
    input:
        final_gzvcf=f"{data_dir}/boech_gbs_matrix_biallelic_filter_DP_hets_mac.vcf.recode.vcf.gz"
    output:
        bed=f"{data_dir}/boech_gbs_matrix_biallelic_filter_DP_hets_mac.vcf.recode.vcf"
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
        freq=f"{data_dir}/boech_gbs_matrix_biallelic_filter_DP_hets_mac.vcf.recode.vcf.afreq",
        bed=f"{data_dir}/boech_gbs_matrix_biallelic_filter_DP_hets_mac.vcf.recode.vcf"
    output:
        rel=f"{data_dir}/boech_gbs_matrix_biallelic_filter_DP_hets_mac.vcf.recode.vcf"
    shell:
        """
        module load plink/2.0
        plink2 --read-freq {input.freq} \
               --bfile {input.bed} \
               --make-rel square \
               --allow-extra-chr \
               --out {output.rel}
        """       
