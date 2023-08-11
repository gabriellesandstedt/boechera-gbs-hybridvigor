################################################################################
## convert vcf to mpgl (genotype-likelihood) file
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s subset_samples_vcf.smk
################################################################################
################################################################################

# Define the paths to data files
data_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS_May23/data"

rule all:
    input:
        f"{data_dir}/allsamples_allsites_final_snps_subset_noSRF2s.vcf",
        f"{data_dir}/allsamples_allsites_final_snps_subset_noSRF2s_split.mpgl",
        f"{data_dir}/allsamples_allsites_final_snps_subset_noSRF2s_chr_pos.txt",
        f"{data_dir}/noSRF2_final.mgpl",
        f"{data_dir}/allsamples_allsites_final_snps_subset_onlyRetros_split.txt",
        f"{data_dir}/allsamples_allsites_final_snps_subset_onlyRetros_split.mpgl",
        f"{data_dir}/allsamples_allsites_final_snps_subset_onlyRetros_chr_pos.txt",
        f"{data_dir}/onlyRetros_final.mgpl"

rule split_noSRF2_vcf:
    input: 
         noSRF2vcf = f"{data_dir}/allsamples_allsites_final_snps_subset_noSRF2s.vcf"
    output:
         noSRF2_split = f"{data_dir}/allsamples_allsites_final_snps_subset_noSRF2s_split.txt"
    shell:
        """
        echo -e "CHROM\tPOS\t$(head -n 1 {input.noSRF2vcf}  | cut -f 10-143)" > {output.noSRF2_split}
        grep -v '^#' {input.noSRF2vcf} | cut -f 10-143 >> {output.noSRF2_split}
        """

rule noSRF2_vcf2mgpl:
    input:
        noSRF2_split = f"{data_dir}/allsamples_allsites_final_snps_subset_noSRF2s_split.txt"
    output:
        noSRF2_mpgl = f"{data_dir}/allsamples_allsites_final_snps_subset_noSRF2s_split.mpgl"
    shell:
        r"""
        python - <<'EOF'
with open("{input.noSRF2_split}", "r") as input_file:
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

with open("{output.noSRF2_mpgl}", "w") as output_file:
    for line in output_data:
        output_file.write(line + '\n')
EOF
        """

rule chrpos_noSRF2:
    input: 
        noSRF2_vcf=f"{data_dir}/allsamples_allsites_final_snps_subset_noSRF2s.vcf"
    output:
        noSRF2_chrpos=f"{data_dir}/allsamples_allsites_final_snps_subset_noSRF2s_chr_pos.txt"
    shell:
        """
        echo -e "CHROM:POS\t$(head -n 1 {input.noSRF2_vcf} | cut -f 1-2)" > {output.noSRF2_chrpos}
        grep -v '^#' {input.noSRF2_vcf} | cut -f 1-2 | awk -F '\t' '{{OFS=":"; print $1, $2}}' >> {output.noSRF2_chrpos}
        """

# for the final mgpl file, I manually replaced row 1 from the output of this rule. 
# first row has two columns (space delimited). col 1 | number of indiviuals, col 2| number of loci
# second row are the individuals, space delimited
rule combine_chr_mpgl_noSRF2:
    input:
        noSRF2_chrpos=f"{data_dir}/allsamples_allsites_final_snps_subset_noSRF2s_chr_pos.txt",
        noSRF2_mgpl=f"{data_dir}/allsamples_allsites_final_snps_subset_noSRF2s_split.mpgl"
    output:
        noSRF2_final_file=f"{data_dir}/noSRF2_final.mgpl"
    shell:
        """
        paste -d ' ' {input.noSRF2_chrpos} {input.noSRF2_mgpl} > {output.noSRF2_final_file}
        """


rule split_onlyRetros_vcf:
    input: 
        onlyRetros_vcf=f"{data_dir}/allsamples_allsites_final_snps_subset_onlyRetros.vcf"
    output:
        onlyRetros_split=f"{data_dir}/allsamples_allsites_final_snps_subset_onlyRetros_split.txt"
    shell:
        """
        echo -e "CHROM\tPOS\t$(head -n 1 {input.onlyRetros_vcf}  | cut -f 10-122)" > {output.onlyRetros_split}
        grep -v '^#' {input.onlyRetros_vcf} | cut -f 10-122 >> {output.onlyRetros_split}
        """

rule onlyRetros_vcf2mgpl:
    input:
        onlyRetros_split=f"{data_dir}/allsamples_allsites_final_snps_subset_onlyRetros_split.txt"
    output:
        onlyRetros_mpgl = f"{data_dir}/allsamples_allsites_final_snps_subset_onlyRetros_split.mpgl"
    shell:
        r"""
        python - <<'EOF'
with open("{input.onlyRetros_split}", "r") as input_file:
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

with open("{output.onlyRetros_mpgl}", "w") as output_file:
    for line in output_data:
        output_file.write(line + '\n')
EOF
        """

rule chrpos_onlyRetros:
    input: 
        onlyRetros_vcf=f"{data_dir}/allsamples_allsites_final_snps_subset_onlyRetros.vcf"
    output:
        onlyRetros_chrpos=f"{data_dir}/allsamples_allsites_final_snps_subset_onlyRetros_chr_pos.txt"
    shell:
        """
        echo -e "CHROM:POS\t$(head -n 1 {input.onlyRetros_vcf} | cut -f 1-2)" > {output.onlyRetros_chrpos}
        grep -v '^#' {input.onlyRetros_mvcf} | cut -f 1-2 | awk -F '\t' '{{OFS=":"; print $1, $2}}' >> {output.onlyRetros_chrpos}
        """

# for the final mgpl file, I manually replaced row 1 from the output of this rule. 
# first row has two columns (space delimited). col 1 | number of indiviuals, col 2| number of loci
# second row are the individuals, space delimited
rule combine_chr_mpgl_onlyRetros:
    input:
        onlyRetros_chrpos=f"{data_dir}/allsamples_allsites_final_snps_subset_onlyRetros_chr_pos.txt",
        onlyRetros_mgpl=f"{data_dir}/allsamples_allsites_final_snps_subset_onlyRetros_split.mpgl"
    output:
        onlyRetros_final_file=f"{data_dir}/onlyRetros_final.mgpl"
    shell:
        """
        paste -d ' ' {input.onlyRetros_chrpos} {input.onlyRetros_mgpl} > {output.onlyRetros_final_file}
        """
