################################################################################
## generate mgpl files for entropy ( stricta x retro and retro x retro )
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 07_entropy.smk
################################################################################
################################################################################

# Define the paths to data files
data_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS_May23/data"
ref_genome = "GCA_018361405.1_NTU_Bstr_LTM_2.2_genomic.fa"
vcf = "boech_gbs_allsamples_SNPs_filtered.vcf.gz"

# define rule to subset only retro - stricta parents and hybrids 
rule retro_str_hybrids:
    input:
        vcf=f"{data_dir}/{vcf}"
    output:
        retro_str_hybrids_vcf=f"{data_dir}/boech_gbs_retro_str_entropy_20MAR.vcf"
    shell:
        """
        module load bcftools/1.16
        echo -e '\\n['$(date)']\\n subset vcf ..\\n'
        samples_to_remove="BIGC_RR,CAPG_RR,CG05_RR,DUNG_RR,LOBG_RR,LPAT_RR,MEYC_RR,MILL_RR,PLCU_RR,QCDH_RR,SQWC_RR,SVBR_RR,TWTQ_RR,WILC_RR,am1_1,am1_2,am1_3,am10_1,am10_2,am10_3,am11_1,am11_2,am11_3,am12_1,am12_2,am12_3,am13_1,am13_2,am13_3,am14_1,am14_2,am14_3,am15_1,am15_2,am15_3,am2_1,am2_2,am2_3,am3_1,am3_2,am3_3,am4_1,am4_2,am4_3,am5_1,am5_2,am5_3,am6_1,am6_2,am6_3,am7_1,am7_2,am7_3,am8_1,am8_2,am8_3,am9_1,am9_2,am9_3,wi1_1,wi1_2,wi1_3,wi10_1,wi10_2,wi10_3,wi2_1,wi2_2,wi2_3,wi3_1,wi3_2,wi3_3,wi4_1,wi4_2,wi4_3,wi5_1,wi5_2,wi5_3,wi6_1,wi6_2,wi6_3,wi7_1,wi7_2,wi7_3,wi8_1,wi8_2,wi8_3,wi9_1,am14_2b,wi2_1b,wi9_2,wi9_3"
        bcftools view -s "^$samples_to_remove" -o {output.retro_str_hybrids_vcf} {input.vcf}
        """

rule retroUp_str_hybrids:
    input:
        vcf=f"{data_dir}/{vcf}"
    output:
        retro_str_hybrids_vcf=f"{data_dir}/boech_gbs_retroUP_str_entropy.vcf"
    shell:
        """
        module load bcftools/1.16
        echo -e '\\n['$(date)']\\n subset vcf ..\\n'
        samples_to_remove="BIGC_RR,CAPG_RR,CG05_RR,DUNG_RR,LOBG_RR,LPAT_RR,MEYC_RR,MILL_RR,PLCU_RR,QCDH_RR,SQWC_RR,SVBR_RR,TWTQ_RR,WILC_RR,am1_1,am1_2,am1_3,am10_1,am10_2,am10_3,am11_1,am11_2,am11_3,am12_1,am12_2,am12_3,am13_1,am13_2,am13_3,am14_1,am14_2,am14_3,am15_1,am15_2,am15_3,am2_1,am2_2,am2_3,am3_1,am3_2,am3_3,am4_1,am4_2,am4_3,am5_1,am5_2,am5_3,am6_1,am6_2,am6_3,am7_1,am7_2,am7_3,am8_1,am8_2,am8_3,am9_1,am9_2,am9_3,wi1_1,wi1_2,wi1_3,wi10_1,wi10_2,wi10_3,wi2_1,wi2_2,wi2_3,wi3_1,wi3_2,wi3_3,wi4_1,wi4_2,wi4_3,wi5_1,wi5_2,wi5_3,wi6_1,wi6_2,wi6_3,wi7_1,wi7_2,wi7_3,wi8_1,wi8_2,wi8_3,wi9_1,am14_2b,wi2_1b,wi9_2,wi9_3,pa10_pool,pa1_pool,pa2_pool,pa3_pool,pa4_pool,pa5_pool,pa6_pool,pa8_pool,pw10_pool,pw11_pool,pw1_pool,pw2_pool,pw3_pool,pw8b_pool,pw9_pool,ES913_pool"
        bcftools view -s "^$samples_to_remove" -o {output.retro_str_hybrids_vcf} {input.vcf}
        """

rule retroLow_str_hybrids:
    input:
        vcf=f"{data_dir}/{vcf}"
    output:
        retro_str_hybrids_vcf=f"{data_dir}/boech_gbs_retroLow_str_entropy.vcf"
    shell:
        """
        module load bcftools/1.16
        echo -e '\\n['$(date)']\\n subset vcf ..\\n'
        samples_to_remove="BIGC_RR,CAPG_RR,CG05_RR,DUNG_RR,LOBG_RR,LPAT_RR,MEYC_RR,MILL_RR,PLCU_RR,QCDH_RR,SQWC_RR,SVBR_RR,TWTQ_RR,WILC_RR,am1_1,am1_2,am1_3,am10_1,am10_2,am10_3,am11_1,am11_2,am11_3,am12_1,am12_2,am12_3,am13_1,am13_2,am13_3,am14_1,am14_2,am14_3,am15_1,am15_2,am15_3,am2_1,am2_2,am2_3,am3_1,am3_2,am3_3,am4_1,am4_2,am4_3,am5_1,am5_2,am5_3,am6_1,am6_2,am6_3,am7_1,am7_2,am7_3,am8_1,am8_2,am8_3,am9_1,am9_2,am9_3,wi1_1,wi1_2,wi1_3,wi10_1,wi10_2,wi10_3,wi2_1,wi2_2,wi2_3,wi3_1,wi3_2,wi3_3,wi4_1,wi4_2,wi4_3,wi5_1,wi5_2,wi5_3,wi6_1,wi6_2,wi6_3,wi7_1,wi7_2,wi7_3,wi8_1,wi8_2,wi8_3,wi9_1,am14_2b,wi2_1b,wi9_2,wi9_3,pw7_pool,pw6b_pool,pw5_pool,pw4_pool,pa9_pool,pa7_pool"
        bcftools view -s "^$samples_to_remove" -o {output.retro_str_hybrids_vcf} {input.vcf}
        """


# define rule to select snps from retro - stricta parents and hybrids vcf
rule select_biallelic_snps_retro_str:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        vcf=f"{data_dir}/boech_gbs_retro_str_entropy.vcf"
    output:
        biallelic_vcf=f"{data_dir}/boech_gbs_retro_str_entropy_SNPs.vcf"
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


rule retroUp_str_snps:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        vcf=f"{data_dir}/boech_gbs_retroUP_str_entropy.vcf"
    output:
        biallelic_vcf=f"{data_dir}/boech_gbs_retroUP_str_entropy_SNPs.vcf"
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

rule retroLow_str_snps:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        vcf=f"{data_dir}/boech_gbs_retroLow_str_entropy.vcf"
    output:
        biallelic_vcf=f"{data_dir}/boech_gbs_retroLow_str_entropy_SNPs.vcf"
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

# 65 samples, mac of (7/130) > 0.05 
# max_missing 10%
#4335 loci
rule filter_minor_allele_count_retro_str:
    input:
        snps_vcf=f"{data_dir}/boech_gbs_retro_str_entropy_SNPs.vcf"
    output:
        filtered_mac_vcf=f"{data_dir}/boech_gbs_retro_str_entropy_SNPs_filtered.vcf",
    shell:
        """
        module load vcftools/0.1.15-6
        vcftools \
            --vcf {input.snps_vcf} \
            --remove-indels \
            --min-alleles 2 \
            --max-alleles 2 \
            --mac 7 \
            --max-missing 0.1 \
            --recode \
            --recode-INFO-all \
            --out {output.filtered_mac_vcf}
        """

#3/49
#5496 SNPs
rule filter_minor_allele_count_retroUP_str:
    input:
        snps_vcf=f"{data_dir}/boech_gbs_retroUP_str_entropy_SNPs.vcf"
    output:
        filtered_mac_vcf=f"{data_dir}/boech_gbs_retroUP_str_entropy_SNPs_filtered.vcf",
    shell:
        """
        module load vcftools/0.1.15-6
        vcftools \
            --vcf {input.snps_vcf} \
            --remove-indels \
            --min-alleles 2 \
            --max-alleles 2 \
            --mac 3 \
            --recode \
            --recode-INFO-all \
            --out {output.filtered_mac_vcf}
        """
#3/59
#5740 SNPs
rule filter_minor_allele_count_retroLow_str:
    input:
        snps_vcf=f"{data_dir}/boech_gbs_retroLow_str_entropy_SNPs.vcf"
    output:
        filtered_mac_vcf=f"{data_dir}/boech_gbs_retroLow_str_entropy_SNPs_filtered.vcf",
    shell:
        """
        module load vcftools/0.1.15-6
        vcftools \
            --vcf {input.snps_vcf} \
            --remove-indels \
            --min-alleles 2 \
            --max-alleles 2 \
            --mac 3 \
            --recode \
            --recode-INFO-all \
            --out {output.filtered_mac_vcf}
        """


rule split_retro_str_hybrids_vcf:
    input: 
        filtered_mac_vcf=f"{data_dir}/boech_gbs_retro_str_entropy_20MAR.vcf"
    output:
        split_vcf=f"{data_dir}/boech_gbs_retro_str_entropy_SNPs_filtered_split_20MAR.txt"
    shell:
        """
        echo -e "CHROM\tPOS\t$(head -n 1 {input.filtered_mac_vcf}  | cut -f 10-143)" > {output.split_vcf}
        grep -v '^#' {input.filtered_mac_vcf} | cut -f 10-143 >> {output.split_vcf}
        """

rule split_retroUP_str_hybrids_vcf:
    input: 
        filtered_mac_vcf=f"{data_dir}/boech_gbs_retroUP_str_entropy_SNPs_filtered.vcf.recode.vcf"
    output:
        split_vcf=f"{data_dir}/boech_gbs_retroUP_str_entropy_SNPs_filtered_split.txt"
    shell:
        """
        echo -e "CHROM\tPOS\t$(head -n 1 {input.filtered_mac_vcf}  | cut -f 10-143)" > {output.split_vcf}
        grep -v '^#' {input.filtered_mac_vcf} | cut -f 10-143 >> {output.split_vcf}
        """

rule split_retroLow_str_hybrids_vcf:
    input: 
        filtered_mac_vcf=f"{data_dir}/boech_gbs_retroLow_str_entropy_SNPs_filtered.vcf.recode.vcf"
    output:
        split_vcf=f"{data_dir}/boech_gbs_retroLow_str_entropy_SNPs_filtered_split.txt"
    shell:
        """
        echo -e "CHROM\tPOS\t$(head -n 1 {input.filtered_mac_vcf}  | cut -f 10-143)" > {output.split_vcf}
        grep -v '^#' {input.filtered_mac_vcf} | cut -f 10-143 >> {output.split_vcf}
        """

# generate mgpl file
# total snps 
rule retro_str_hybrids_vcf2mgpl:
    input:
        split_vcf=f"{data_dir}/boech_gbs_retro_str_entropy_SNPs_filtered_split_20MAR.txt"
    output:
        split_mgpl=f"{data_dir}/boech_gbs_retro_str_entropy_SNPs_filtered_20MAR.mgpl"
    shell:
        r"""
        python - <<'EOF'
with open("{input.split_vcf}", "r") as input_file:
    data = input_file.read()

rows = data.strip().split('\n')

output_data = []

for row in rows:
    columns = row.split('\t')
    extracted_values = []

    for col in columns:
        values = col.split(':')
        for v in values:
            comma_values = v.split(',')
            if len(comma_values) == 3:
                extracted_values.extend(comma_values)
                break
        else:
            extracted_values.extend(['0', '0', '0'])

    extracted_values_str = ' '.join(extracted_values)
    output_data.append(extracted_values_str)

with open("{output.split_mgpl}", "w") as output_file:
    for line in output_data:
        output_file.write(line + '\n')
EOF
        """

# extract chromosome and position info 
rule chrpos_retro_str:
    input: 
        retro_str_vcf=f"{data_dir}/boech_gbs_retro_str_entropy_20MAR.vcf"
    output:
        retro_str_chrpos=f"{data_dir}/boech_gbs_retro_str_entropy_SNPs_filtered_chr_pos_20MAR.txt"
    shell:
        """
        echo -e "CHROM:POS\t$(head -n 1 {input.retro_str_vcf} | cut -f 1-2)" > {output.retro_str_chrpos}
        grep -v '^#' {input.retro_str_vcf} | cut -f 1-2 | awk -F '\t' '{{OFS=":"; print $1, $2}}' >> {output.retro_str_chrpos}
        """

# for the final mgpl file, I manually replaced row 1 from the output of this rule. 
# first row has two columns (space delimited). col 1 | number of individuals, col 2| number of loci
# second row are the individuals, space delimited
rule combine_chr_mpgl_retro_str:
    input:
        retro_str_chrpos=f"{data_dir}/boech_gbs_retro_str_entropy_SNPs_filtered_chr_pos_20MAR.txt",
        mgpl=f"{data_dir}/boech_gbs_retro_str_entropy_SNPs_filtered.mgpl"
    output:
        final_mgpl=f"{data_dir}/boech_gbs_retro_str_entropy_final.mgpl"
    shell:
        """
        paste -d ' ' {input.retro_str_chrpos} {input.mgpl} > {output.final_mgpl}
        """


# generate mgpl file
rule retro_str_hybrids_vcf2mgpl:
    input:
        split_vcf=f"{data_dir}/boech_gbs_retroLow_str_entropy_SNPs_filtered_split.txt"
    output:
        split_mgpl=f"{data_dir}/boech_gbs_retroLow_str_entropy_SNPs_filtered.mgpl"
    shell:
        r"""
        python - <<'EOF'
with open("{input.split_vcf}", "r") as input_file:
    data = input_file.read()

rows = data.strip().split('\n')

output_data = []

for row in rows:
    columns = row.split('\t')
    extracted_values = []

    for col in columns:
        values = col.split(':')
        for v in values:
            comma_values = v.split(',')
            if len(comma_values) == 3:
                extracted_values.extend(comma_values)
                break
        else:
            extracted_values.extend(['0', '0', '0'])

    extracted_values_str = ' '.join(extracted_values)
    output_data.append(extracted_values_str)

with open("{output.split_mgpl}", "w") as output_file:
    for line in output_data:
        output_file.write(line + '\n')
EOF
        """

# extract chromosome and position info 
rule chrpos_retroUP_str:
    input: 
        retro_str_vcf=f"{data_dir}/boech_gbs_retroLow_str_entropy_SNPs_filtered.vcf.recode.vcf"
    output:
        retro_str_chrpos=f"{data_dir}/boech_gbs_retroLow_str_entropy_SNPs_filtered_chr_pos.txt"
    shell:
        """
        echo -e "CHROM:POS\t$(head -n 1 {input.retro_str_vcf} | cut -f 1-2)" > {output.retro_str_chrpos}
        grep -v '^#' {input.retro_str_vcf} | cut -f 1-2 | awk -F '\t' '{{OFS=":"; print $1, $2}}' >> {output.retro_str_chrpos}
        """

# for the final mgpl file, I manually replaced row 1 from the output of this rule. 
# first row has two columns (space delimited). col 1 | number of indiviuals, col 2| number of loci
# second row are the individuals, space delimited
rule combine_chr_mpgl_retro_str:
    input:
        retro_str_chrpos=f"{data_dir}/boech_gbs_retroLow_str_entropy_SNPs_filtered_chr_pos.txt",
        mgpl=f"{data_dir}/boech_gbs_retroLow_str_entropy_SNPs_filtered.mgpl"
    output:
        final_mgpl=f"{data_dir}/boech_gbs_retroLow_str_entropy_final.mgpl"
    shell:
        """
        paste -d ' ' {input.retro_str_chrpos} {input.mgpl} > {output.final_mgpl}
        """

# generate mgpl file
# total snps 5486
rule retro_str_hybrids_vcf2mgpl:
    input:
        split_vcf=f"{data_dir}/boech_gbs_retro_str_entropy_SNPs_filtered_split.txt"
    output:
        split_mgpl=f"{data_dir}/boech_gbs_retro_str_entropy_SNPs_filtered.mgpl"
    shell:
        r"""
        python - <<'EOF'
with open("{input.split_vcf}", "r") as input_file:
    data = input_file.read()

rows = data.strip().split('\n')

output_data = []

for row in rows:
    columns = row.split('\t')
    extracted_values = []

    for col in columns:
        values = col.split(':')
        for v in values:
            comma_values = v.split(',')
            if len(comma_values) == 3:
                extracted_values.extend(comma_values)
                break
        else:
            extracted_values.extend(['0', '0', '0'])

    extracted_values_str = ' '.join(extracted_values)
    output_data.append(extracted_values_str)

with open("{output.split_mgpl}", "w") as output_file:
    for line in output_data:
        output_file.write(line + '\n')
EOF
        """

# extract chromosome and position info 
rule chrpos_retro_str:
    input: 
        retro_str_vcf=f"{data_dir}/boech_gbs_retro_str_entropy_SNPs_filtered.vcf.recode.vcf"
    output:
        retro_str_chrpos=f"{data_dir}/boech_gbs_retro_str_entropy_SNPs_filtered_chr_pos.txt"
    shell:
        """
        echo -e "CHROM:POS\t$(head -n 1 {input.retro_str_vcf} | cut -f 1-2)" > {output.retro_str_chrpos}
        grep -v '^#' {input.retro_str_vcf} | cut -f 1-2 | awk -F '\t' '{{OFS=":"; print $1, $2}}' >> {output.retro_str_chrpos}
        """

# for the final mgpl file, I manually replaced row 1 from the output of this rule. 
# first row has two columns (space delimited). col 1 | number of indiviuals, col 2| number of loci
# second row are the individuals, space delimited
#  loci 4335
rule combine_chr_mpgl_retro_str:
    input:
        retro_str_chrpos=f"{data_dir}/boech_gbs_retro_str_entropy_SNPs_filtered_chr_pos.txt",
        mgpl=f"{data_dir}/boech_gbs_retro_str_entropy_SNPs_filtered.mgpl"
    output:
        final_mgpl=f"{data_dir}/boech_gbs_retro_str_entropy_final.mgpl"
    shell:
        """
        paste -d ' ' {input.retro_str_chrpos} {input.mgpl} > {output.final_mgpl}
        """
# define rule to subset only retro - retro parents and hybrids 
rule retro_retro_hybrids:
    input:
        vcf=f"{data_dir}/{vcf}"
    output:
        retro_retro_hybrids_vcf=f"{data_dir}/boech_gbs_retro_retro_entropy.vcf"
    shell:
        """
        module load bcftools/1.16
        echo -e '\\n['$(date)']\\n subset vcf ..\\n'
        samples_to_remove="ALDU_SR,BAYH_SR,BLAC_SR,CARM_SR,EH_11,EH_4,EH_9,EJ1_1,EJ1_2,EJ1_4,EJ2_1,EJ2_3,EJ2_7,ET1_12,ET1_6,ET1_9,ET1_ex1,ET2_11,ET2_3,ET2_6,EU4_3,EU4_7,EU4_8,EU5_1,EU5_2,EU5_5,HESS_pool,IRCP_SR,JAM_pool,LPAT_SR,LTM,MEYF_SR,MILL_SR,ODEL_SR,PLCU_SR,QCDH_SR,SD01_SR,SQWL_SR,SVOV_SR,TWMUp_pool,TWTN_SR,Usilv_pool,WILU_SR"
        bcftools view -s "^$samples_to_remove" -o {output.retro_retro_hybrids_vcf} {input.vcf}
        """

# define rule to filter only retro - retro parents and hybrids vcf 
rule select_biallelic_snps_retro_retro:
    input:
        ref=f"{ref_dir}/{ref_genome}",
        vcf=f"{data_dir}/boech_gbs_retro_retro_entropy.vcf"
    output:
        biallelic_vcf=f"{data_dir}/boech_gbs_retro_retro_entropy_SNPs.vcf"
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

# 113 samples, mac of 12 >0.05 (12/226)
# 2134
rule filter_minor_allele_count_retro_retro:
    input:
        snps_vcf=f"{data_dir}/boech_gbs_retro_retro_entropy_SNPs.vcf"
    output:
        filtered_mac_vcf=f"{data_dir}/boech_gbs_retro_retro_entropy_SNPs_filtered.vcf",
    shell:
        """
        module load vcftools/0.1.15-6
        vcftools \
            --vcf {input.snps_vcf} \
            --remove-indels \
            --min-alleles 2 \
            --max-alleles 2 \
            --max-missing 0.1 \
            --mac 12 \
            --recode \
            --recode-INFO-all \
            --out {output.filtered_mac_vcf}
        """

rule split_retro_retro_hybrids_vcf:
    input: 
        filtered_mac_vcf=f"{data_dir}/boech_gbs_retro_retro_entropy_SNPs_filtered.vcf.recode.vcf"
    output:
        split_vcf=f"{data_dir}/boech_gbs_retro_retro_entropy_SNPs_filtered_split.txt"
    shell:
        """
        echo -e "CHROM\tPOS\t$(head -n 1 {input.filtered_mac_vcf}  | cut -f 10-143)" > {output.split_vcf}
        grep -v '^#' {input.filtered_mac_vcf} | cut -f 10-143 >> {output.split_vcf}
        """

# generate mgpl file
# total snps: 2134
rule retro_retro_hybrids_vcf2mgpl:
    input:
        split_vcf=f"{data_dir}/boech_gbs_retro_retro_entropy_SNPs_filtered_split.txt"
    output:
        split_mgpl=f"{data_dir}/boech_gbs_retro_retro_entropy_SNPs_filtered.mgpl"
    shell:
        r"""
        python - <<'EOF'
with open("{input.split_vcf}", "r") as input_file:
    data = input_file.read()

rows = data.strip().split('\n')

output_data = []

for row in rows:
    columns = row.split('\t')
    extracted_values = []

    for col in columns:
        values = col.split(':')
        for v in values:
            comma_values = v.split(',')
            if len(comma_values) == 3:
                extracted_values.extend(comma_values)
                break
        else:
            extracted_values.extend(['0', '0', '0'])

    extracted_values_str = ' '.join(extracted_values)
    output_data.append(extracted_values_str)

with open("{output.split_mgpl}", "w") as output_file:
    for line in output_data:
        output_file.write(line + '\n')
EOF
        """

# extract chromosome and position info 
rule chrpos_retro_retro:
    input: 
        retro_vcf=f"{data_dir}/boech_gbs_retro_retro_entropy_SNPs_filtered.vcf.recode.vcf"
    output:
        retro_chrpos=f"{data_dir}/boech_gbs_retro_retro_entropy_SNPs_filtered_chr_pos.txt"
    shell:
        """
        echo -e "CHROM:POS\t$(head -n 1 {input.retro_vcf} | cut -f 1-2)" > {output.retro_chrpos}
        grep -v '^#' {input.retro_vcf} | cut -f 1-2 | awk -F '\t' '{{OFS=":"; print $1, $2}}' >> {output.retro_chrpos}
        """

# for the final mgpl file, I manually replaced row 1 from the output of this rule. 
# first row has two columns (space delimited). col 1 | number of indiviuals, col 2| number of loci
# second row are the individuals, space delimited
rule combine_chr_mpgl_retro_retro:
    input:
        retro_chrpos=f"{data_dir}/boech_gbs_retro_retro_entropy_SNPs_filtered_chr_pos.txt",
        mgpl=f"{data_dir}/boech_gbs_retro_retro_entropy_SNPs_filtered.mgpl"
    output:
        final_mgpl=f"{data_dir}/boech_gbs_retro_retro_entropy_final.mgpl"
    shell:
        """
        paste -d ' ' {input.retro_chrpos} {input.mgpl} > {output.final_mgpl}
        """
