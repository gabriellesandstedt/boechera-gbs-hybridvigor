################################################################################
## filter vcf for the genetic relatedness matrix 
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete --latency-wait 60 --cores 4 -s 06_filter_genetic_matrix_vcf.smk
################################################################################
################################################################################
# define directories
data_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS_May23/data"
ref_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS/ref_genomes/data"

# define data file
input_vcf = "boechera_gbs_matrix.vcf"
# define the reference genome filename
ref_genome = "GCA_018361405.1_NTU_Bstr_LTM_2.2_genomic.fa"

# define the output filenames
biallelic_vcf = "boech_gbs_matrix_biallelic.vcf"
variant_table = "boech_gbs_matrix_variant.table" # indels and snps 
filtered_vcf = "boech_gbs_matrix_biallelic_filter.vcf"
filtered_passed_vcf = "boech_gbs_matrix_biallelic_filterPASSED.vcf"

# select biallelic SNPs
rule select_biallelic_snps:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        vcf=f"{data_dir}/{input_vcf}"
    output:
        biallelic_vcf=f"{data_dir}/{biallelic_vcf}"
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
        vcf=f"{data_dir}/{input_vcf}"
    output:
        table=f"{data_dir}/{variant_table}"
    shell:
        """
        gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
        -O {output.table}
        """
       
rule plot_densities:
    output:
        "boech_gbs_matrix_variant_plot_densities.pdf"
    shell:
        """
        Rscript filtering_diagnostics.R
        """

# filter variants based on plot_densities rule
rule filter_variants:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        biallelic_vcf=f"{data_dir}/{biallelic_vcf}"
    output:
        filtered_vcf=f"{data_dir}/{filtered_vcf}"
    shell:
        """
        gatk VariantFiltration \
        -R {input.ref} \
        -V {input.biallelic_vcf} \
        --filterExpression "QUAL < 0 || MQ < 40.00 || SOR > 4.000 || QD < 2.00 || FS > 60.000 || MQRankSum < -12.500 || ReadPosRankSum < -10.000 || ReadPosRankSum > 10.000" \
        --filterName "my_snp_filter" \
        -o {output.filtered_vcf}
        """

# extract passed variants
rule extract_passed_variants:
    input:
        filtered_vcf=f"{data_dir}/{filtered_vcf}"
    output:
        filtered_passed_vcf=f"{data_dir}/{filtered_passed_vcf}"
    shell:
        """
        grep -E '^#|PASS' {input.filtered_vcf} > {output.filtered_passed_vcf}
        """

