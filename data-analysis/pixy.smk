################################################################################
## Run pixy to determine Fst, Pi, and Dxy
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete --latency-wait 60 --cores 4 -s pixy.smk
################################################################################
## Pixy documentation: Manuscript: Korunes, K.L. and Samuk, K. (2021), pixy: Unbiased estimation of nucleotide diversity and divergence in the presence of missing data. Molecular Ecology Resources. Accepted Author Manuscript. https://doi.org/10.1111/1755-0998.13326
################################################################################
# define data directory
data_dir = "/scratch/general/nfs1/u6048240/BOECHERA/GBS_May23/data"

# assign rule to run pixy v 1.2.3: https://pixy.readthedocs.io/en/latest/
# outputs are pi, fst, and dxy txt files. 
rule pixy_stats:
    input:
        vcf=f"{data_dir}/boech_gbs_allsamples_combined_final.vcf.gz",
        pop_file=f"{data_dir}/pop_pixy.txt"
    output:
        pixy_output=f"{data_dir}/pixy_stats.txt"
    params:
        window_size=100000,
        n_cores=4
    shell:
        """
        ml pixy/1.2.3
        pixy --stats pi fst dxy \
        --vcf {input.vcf} \
        --populations {input.pop_file} \
        --window_size {params.window_size} \
        --n_cores {params.n_cores} \
        > {output.pixy_output}
        """
