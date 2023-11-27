################################################################################
## Run Sabre to demultiplex GBS samples
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## Snakemake v 7.25.0 : https://snakemake.github.io
## To run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 01_demulitplex.smk
################################################################################
################################################################################
# assign directories
raw_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS_DEC23/data"
data_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS_DEC23/data1"

# Define rule all output files
rule all:
    input:
        expand(f"{data_dir}/unknown_barcodes_C6NP3ANXX_7"),
        expand(f"{data_dir}/unknown_barcodes_C6NP3ANXX_8")

# Define rule for demultiplexing with SABRE
# note that the fastqs downloaded from HD are labeled *_fastq* vs *.fastq*
# Sabre v 1.00 : https://github.com/najoshi/sabre
rule sabre_demultiplex:
    input:
        FQ1 = f"{raw_dir}/C6NP3ANXX_7.fastq.gz",
        B1 = f"{raw_dir}/barcodes_C6NP3ANXX_7",
        FQ2 = f"{raw_dir}/C6NP3ANXX_8.fastq.gz",
        B2 = f"{raw_dir}/barcodes_C6NP3ANXX_8",
    output:
        unknown_barcodes1 = expand(f"{data_dir}/unknown_barcodes_C6NP3ANXX_7"),
        unknown_barcodes2 = expand(f"{data_dir}/unknown_barcodes_C6NP3ANXX_8")
    shell:
        """
        echo -e "\n["$(date)"]\n demultiplex fastq files ...\n"
            sabre se -f {input.FQ1} -b {input.B1} -u {output.unknown_barcodes1};
            sabre se -f {input.FQ2} -b {input.B2} -u {output.unknown_barcodes2};
        """

################################################################################
## To unlock a snakemake file: snakemake --unlock --snakefile 01_demultiplex.smk
