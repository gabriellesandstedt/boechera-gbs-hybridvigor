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
        vcf=f"{data_dir}/{vcf}"
    output:
        noSRF2s_vcf=f"{data_dir}/allsamples_allsites_final_snps_subset_noSRF2s.vcf"
    shell:
        """
        module load bcftools/1.16
        echo -e '\\n['$(date)']\\n subset vcf ..\\n'
        samples_to_remove="EH_11,EH_4,EH_9,EJ1_1,EJ1_2,EJ1_4,EJ2_1,EJ2_3,EJ2_7,ET1_12,ET1_6,ET1_9,ET1_ex1,ET2_11,ET2_3,ET2_6,EU4_3,EU4_7,EU4_8,EU5_1,EU5_2,EU5_5"
        bcftools view -s "^$samples_to_remove" -o {output.noSRF2s_vcf} {input.vcf}
        """


# define rule to subset vcf to 
# bcftools: 
rule noSRF2s:
    input:
        vcf=f"{data_dir}/{vcf}"
    output:
        onlyRetros_vcf=f"{data_dir}/allsamples_allsites_final_snps_subset_onlyRetros.vcf"
    shell:
        """
        module load bcftools/1.16
        echo -e '\\n['$(date)']\\n subset vcf ..\\n'
        samples_to_remove="ALDU_SR,BAYH_SR,BLAC_SR,CARM_SR,EH_11,EH_4,EH_9,EJ1_1,EJ1_2,EJ1_4,EJ2_1,EJ2_3,EJ2_7,ET1_12,ET1_6,ET1_9,ET1_ex1,ET2_11,ET2_3,ET2_6,EU4_3,EU4_7,EU4_8,EU5_1,EU5_2,EU5_5,HESS_pool,IRCP_SR,JAM_pool,LPAT_SR,LTM,MEYF_SR,MILL_SR,ODEL_SR,PLCU_SR,QCDH_SR,SD01_SR,SQWL_SR,SVOV_SR,TWMUp_pool,TWTN_SR,Usilv_pool,WILU_SR"
        bcftools view -s "^$samples_to_remove" -o {output.onlyRetros_vcf} {input.vcf}
        """

