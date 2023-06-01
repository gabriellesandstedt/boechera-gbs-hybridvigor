################################################################################
## filter vcf for the genetic relatedness matrix 
## filter methods were adapted from: https://evodify.com/gatk-in-non-model-organism/
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

# select biallelic SNPs
rule select_biallelic_snps:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        vcf=f"{data_dir}/boechera_gbs_matrix.vcf"
    output:
        biallelic_vcf=f"{data_dir}/boech_gbs_matrix_biallelic.vcf"
    shell:
        """
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
        gatk VariantsToTable \
            -R {input.ref} \
            -V {input.vcf} \
            -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
            -O {output.table}
        
        Rscript {input.rscript}
        """

# filter variants based on plot_densities rule
rule filter_variant_diagnostics:
    output:
        f"{data_dir}/boech_gbs_matrix_variant_plot_densities.pdf"
    shell:
        """
        Rscript {scripts_dir}/filtering_diagnostics.R
        """

# filter variants
rule filter_variants:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        biallelic_vcf=f"{data_dir}/boech_gbs_matrix_biallelic.vcf"
    output:
        filtered_vcf=f"{data_dir}/boech_gbs_matrix_biallelic_filter.vcf"
    shell:
        """
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.biallelic_vcf} \
            --filter-expression "QUAL < 0 || MQ < 40.00 || SOR > 4.000 || QD < 2.00 || FS > 60.000 || MQRankSum < -12.500 || ReadPosRankSum < -10.000 || ReadPosRankSum > 10.000" \
            --filter-name "my_snp_filter" \
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
        java {params.java_opts} -jar {params.gatk_jar} \
            -T VariantsToTable \
            -R {input.ref} \
            -V {input.vcf} \
            -F CHROM -F POS -GF GT -GF DP \
            -o {output.dp_table}
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
        java -jar ~/Programs/GATK/GenomeAnalysisTK.jar \
            -T VariantFiltration \
            -R {input.ref} \
            -V {input.filtered_passed_vcf} \
            -G-filter "DP < 5 || DP > 100" \
            -G-filter-name "DP_5-100" \
            -o {output.filtered_DP_vcf}
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
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.filtered_DP_vcf} \
            --set-filtered-gt-to-nocall \
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
        Rscript {input.rscript}
        """

# filter minor allele count
rule filter_minor_allele_count:
    input:
        filtered_hets_vcf=f"{data_dir}/boech_gbs_matrix_biallelic_filter_DP_hets.vcf"
    output:
        filtered_mac_vcf=f"{data_dir}/boech_gbs_matrix_biallelic_filter_DP_hets_mac.vcf"
    shell:
        """
        vcftools \
            --vcf {input.filtered_hets_vcf} \
            --remove-indels \
            --min-alleles 2 \
            --max-alleles 2 \
            --mac 3 \
            --recode \
            --recode-INFO-all \
            --out {output.filtered_mac_vcf}
        """

