################################################################################
## filter vcf containing all samples and joint genotyped at all sites
## most filtering methods were modified from: https://evodify.com/gatk-in-non-model-organism/
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete --latency-wait 60 --cores 4 -s 06_filter_allsites_allsamples_vcf.smk
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
        f"{data_dir}/boech_gbs_allsites_filter_DP_hets_mac.vcf.recode.vcf.gz",
        f"{data_dir}/boech_gbs_allsites_filter_DP_hets_mac.vcf.recode.vcf.gz.tbi",
        f"{data_dir}/boech_gbs_allsamples_invariant_geno_called_filter_DP.vcf.gz",
        f"{data_dir}/boech_gbs_allsamples_invariant_geno_called_filter_DP.vcf.gz.tbi",
        f"{data_dir}/boech_gbs_allsamples_combined_final.vcf.gz"

# select biallelic SNPs
rule select_biallelic_snps:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        vcf=f"{data_dir}/boech_gbs_allsamples.vcf"
    output:
        biallelic_vcf=f"{data_dir}/boech_gbs_allsamples_biallelic_snps.vcf"
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

# select invariant sites
rule select_invariant_sites:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        vcf=f"{data_dir}/boech_gbs_allsamples.vcf"
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

# define rule to only print positions with called genotypes   
# I ran this rule to speed up following steps
# However, downstream analyses using windows must be specified to start at pos 1
rule select_genotyped_invariant_sites:
    input:
        invariant_vcf=f"{data_dir}/boech_gbs_allsamples_invariant.vcf",
    output: 
        invariant_vcf2=f"{data_dir}/boech_gbs_allsamples_invariant_geno_called.vcf"
    run:
        with open(input.invariant_vcf, 'r') as f:
            with open(output.invariant_vcf2, 'w') as out_file:
                for line in f:
                    if line.startswith('#'):
                        print(line.strip(), file=out_file)
                    else:
                        fields = line.strip().split('\t')
                        if fields[6] == "." and fields[7] == ".":
                            continue
                        else:
                            print(line.strip(), file=out_file)

# create variant table
rule variant_table:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        biallelic_vcf=f"{data_dir}/boech_gbs_allsamples_biallelic_snps.vcf",
        rscript=f"{scripts_dir}/filtering_diagnostics.R"
    output:
        table=f"{data_dir}/boech_gbs_allsamples_variant.table"
    shell:
        """
        module load gatk/4.1
        module load R
        gatk VariantsToTable \
            -R {input.ref} \
            -V {input.biallelic_vcf} \
            -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
            -O {output.table}
        
        Rscript {input.rscript}
        """

# create invariant table
rule invariant_table:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        biallelic_vcf=f"{data_dir}/boech_gbs_allsamples_invariant_geno_called.vcf",
        rscript=f"{scripts_dir}/filtering_diagnostics.R"
    output:
        table=f"{data_dir}/boech_gbs_allsamples_invariant.table"
    shell:
        """
        module load gatk/4.1
        module load R
        gatk VariantsToTable \
            -R {input.ref} \
            -V {input.biallelic_vcf} \
            -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
            -O {output.table}
        
        Rscript {input.rscript}
        """

# rule to filter variants:
rule filter_variants:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        biallelic_vcf=f"{data_dir}/boech_gbs_allsamples_biallelic_snps.vcf"
    output:
        filtered_vcf=f"{data_dir}/boech_gbs_allsamples_biallelic_snps_filter.vcf"
    shell:
        """
        module load gatk/4.1
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.biallelic_vcf} \
            --filter-expression "QD < 2.0" --filter-name "QD2" \
            --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
            --filter-expression "SOR > 3.0" --filter-name "SOR3" \
            --filter-expression "FS > 60.0" --filter-name "FS60" \
            --filter-expression "MQ < 40.0" --filter-name "MQ40" \
            --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -O {output.filtered_vcf}
        """

# extract passed variants
rule extract_passed_variants:
    input:
        filtered_vcf=f"{data_dir}/boech_gbs_allsamples_biallelic_snps_filter.vcf"
    output:
        filtered_passed_vcf=f"{data_dir}/boech_gbs_allsamples_biallelic_snps_filterPASSED.vcf"
    shell:
        """
        grep -E '^#|PASS' {input.filtered_vcf} > {output.filtered_passed_vcf}
        """

# rule to filter invariants:
rule filter_invariants:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        invcf=f"{data_dir}/boech_gbs_allsamples_invariant_geno_called.vcf"
    output:
        filtered_invcf=f"{data_dir}/boech_gbs_allsamples_invariant_filter.vcf"
    shell:
        """
        module load gatk/4.1
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.invcf} \
            --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
            --filter-expression "MQ < 40.0" --filter-name "MQ40" \
            --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -O {output.filtered_invcf}
        """

# extract passed invariants
rule extract_passed_variants:
    input:
        filtered_invcf=f"{data_dir}/boech_gbs_allsamples_invariant_filter.vcf"
    output:
        filtered_passed_invcf=f"{data_dir}/boech_gbs_allsamples_invariant_filterPASSED.vcf"
    shell:
        """
        grep -E '^#|PASS' {input.filtered_invcf} > {output.filtered_passed_invcf}
        """



# extract depth information from whole-genome VCF
rule table_for_depth:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        vcf=f"{data_dir}/boech_gbs_allsamples.vcf"
    output:
        dp_table=f"{data_dir}/boech_gbs_allsamples.DP.table"
    shell:
        """
        module load gatk/4.1
        gatk VariantsToTable \
            -R {input.ref} \
            -V {input.vcf} \
            -F CHROM -F POS -GF GT -GF DP \
            -O {output.dp_table}
        """

# define rule to only print rows of table if there are called genotypes
# this rule removes positions that are missing genotypes from all samples
# then runs R script that assesses the depth to determine filtering thresholds
rule filter_table_for_called_genotypes:
    input:
        dp_table="boech_gbs_allsamples.DP.table"
        rscript=f"{scripts_dir}/filtering_diagnostics_DP.R"
    output:
        filtered_dp_table="filtered_boech_gbs_allsamples.DP.table"
    shell:
        """
        module load R
        awk -F '\t' 'BEGIN {{OFS="\t"}} {{valid=0; for(i=3; i<=314; i++) {{if($i != "./." && $i != "0") {{valid=1; break;}}}} if(valid) {{print}}}}' {input.dp_table}  > {output.filtered_dp_table} 
     
        Rscript {input.rscript}
        """

# filter variants based on DP
rule filter_variants_DP:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        filtered_passed_vcf=f"{data_dir}/boech_gbs_allsamples_biallelic_snps_filterPASSED.vcf"
    output:
        filtered_DP_vcf=f"{data_dir}/boech_gbs_allsamples_biallelic_snps_filterPASSED_DP.vcf"
    shell:
        """
        module load gatk/4.1
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.filtered_passed_vcf} \
            -G-filter "DP < 5 || DP > 200" \
            -G-filter-name "DP_5-200" \
            -O {output.filtered_DP_vcf}
        """

# select variants and set filtered genotypes to no-call
rule select_variants:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        filtered_DP_vcf=f"{data_dir}/boech_gbs_allsamples_biallelic_snps_filterPASSED_DP.vcf"
    output:
        filtered_noCall_vcf=f"{data_dir}/boech_gbs_allsamples_biallelic_snps_filterPASSED_DPfilterNoCall.vcf"
    shell:
        """
        module load gatk/4.1
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.filtered_DP_vcf} \
            --set-filtered-gt-to-nocall TRUE \
            -O {output.filtered_noCall_vcf}
        """
        
# define rule to filter invariants by depth       
rule filter_invariants_DP:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        filtered_passed_vcf=f"{data_dir}/boech_gbs_allsamples_invariant_filterPASSED.vcf"
    output:
        filtered_DP_vcf=f"{data_dir}/boech_gbs_allsamples_invariant_filterPASSED_DP.vcf"
    shell:
        """
        module load gatk/4.1
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.filtered_passed_vcf} \
            -G-filter "DP < 5 || DP > 200" \
            -G-filter-name "DP_5-200" \
            -O {output.filtered_DP_vcf}
        """

# select invariants and set filtered genotypes to no-call
rule select_invariants:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        filtered_DP_vcf=f"{data_dir}/boech_gbs_allsamples_invariant_filterPASSED_DP.vcf"
    output:
        filtered_noCall_vcf=f"{data_dir}/boech_gbs_allsamples_invariant_filterPASSED_DPfilterNoCall.vcf"
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
# input file for R script: boech_gbs_allsamples_biallelic_snps_filterPASSED_DPfilterNoCall.vcf
rule filter_heterozygous_genotypes:
    input:
        rscript=f"{scripts_dir}/filter_heterozygous_genotypes.R"
    output:
        f"{data_dir}/boech_gbs_allsamples_biallelic_snps_filter_DP_hets.vcf"
    shell:
        """
        module load R
        Rscript {input.rscript}
        """

# filter minor allele count
# mac 8/156 samples, >0.05% 
# no. of biallelic snps: 5316
rule filter_minor_allele_count:
    input:
        filtered_hets_vcf=f"{data_dir}/boech_gbs_allsamples_biallelic_snps_filter_DP_hets.vcf"
    output:
        filtered_mac_vcf=f"{data_dir}/boech_gbs_allsamples_biallelic_snps_filter_DP_hets_mac.vcf"
    shell:
        """
        module load vcftools/0.1.15-6
        vcftools \
            --vcf {input.filtered_hets_vcf} \
            --remove-indels \
            --min-alleles 2 \
            --max-alleles 2 \
            --mac 8 \
            --recode \
            --recode-INFO-all \
            --out {output.filtered_mac_vcf}
        """

# define rule to bgzip vcf files
rule vcf_to_gzvcf:
    input:
        var_vcf=f"{data_dir}/boech_gbs_allsamples_biallelic_snps_filter_DP_hets_mac.vcf.recode.vcf",
        inv_vcf=f"{data_dir}/boech_gbs_allsamples_invariant_filterPASSED_DPfilterNoCall.vcf"
    output:
        gz_var_vcf=f"{data_dir}/boech_gbs_allsamples_biallelic_snps_filter_DP_hets_mac.vcf.recode.vcf.gz",
        gz_invar_vcf=f"{data_dir}/boech_gbs_allsamples_invariant_geno_called_filter_DP.vcf.gz",
        tabix_var_vcf=f"{data_dir}/boech_gbs_allsamples_biallelic_snps_filter_DP_hets_mac.vcf.recode.vcf.gz.tbi",
        tabix_invar_vcf=f"{data_dir}/boech_gbs_allsamples_invariant_geno_called_filter_DP.vcf.gz.tbi"
    shell:
        """
        module load htslib/1.16
        bgzip {input.var_vcf}
        bgzip {input.inv_vcf}
        tabix -p vcf {output.gz_var_vcf}
        tabix -p vcf {output.gz_invar_vcf}
        """

rule combine_vcfs:
    input:
       gz_var_vcf=f"{data_dir}/boech_gbs_allsamples_biallelic_snps_filter_DP_hets_mac.vcf.recode.vcf.gz",
       gz_invar_vcf=f"{data_dir}/boech_gbs_allsamples_invariant_filterPASSED_DPfilterNoCall.vcf.gz"
    output:
       final_vcf=f"{data_dir}/boech_gbs_allsamples_combined_final.vcf.gz"
    shell:
        """
        module load bcftools/1.16
        bcftools concat {input.gz_var_vcf} {input.gz_invar_vcf} -a -Oz -o {output.final_vcf}
        """       
        
