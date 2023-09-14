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

# define rule to subset vcf to remove retros hybrids 
# bcftools: 
rule noRhybrids:
    input:
        vcf=f"{data_dir}/{vcf}"
    output:
        noRhybrids_vcf=f"{data_dir}/allsamples_allsites_final_snps_subset_noRhybrids.vcf"
    shell:
        """
        module load bcftools/1.16
        echo -e '\\n['$(date)']\\n subset vcf ..\\n'
        samples_to_remove="BIGC_RR,CAPG_RR,CG05_RR,DUNG_RR,LOBG_RR,LPAT_RR,MEYC_RR,MILL_RR,PLCU_RR,QCDH_RR,SQWC_RR,SVBR_RR,TWTQ_RR,WILC_RR,am1_1,am1_2,am1_3,am10_1,am10_2,am10_3,am11_1,am11_2,am11_3,am12_1,am12_2,am12_3,am13_1,am13_2,am13_3,am14_1,am14_2,am14_3,am15_1,am15_2,am15_3,am2_1,am2_2,am2_3,am3_1,am3_2,am3_3,am4_1,am4_2,am4_3,am5_1,am5_2,am5_3,am6_1,am6_2,am6_3,am7_1,am7_2,am7_3,am8_1,am8_2,am8_3,am9_1,am9_2,am9_3,wi1_1,wi1_2,wi1_3,wi10_1,wi10_2,wi10_3,wi2_1,wi2_2,wi2_3,wi3_1,wi3_2,wi3_3,wi4_1,wi4_2,wi4_3,wi5_1,wi5_2,wi5_3,wi6_1,wi6_2,wi6_3,wi7_1,wi7_2,wi7_3,wi8_1,wi8_2,wi8_3,wi9_1,am14_2b,wi2_1b,wi9_2,wi9_3"
        bcftools view -s "^$samples_to_remove" -o {output.noRhybrids_vcf} {input.vcf}
        """

# define rule to subset vcf to 
# bcftools: 
rule onlyRetros:
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

