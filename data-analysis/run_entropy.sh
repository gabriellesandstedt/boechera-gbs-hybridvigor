#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=4 #parallel tasks/processes run on separate CPU cores or nodes (coarse-grained internode/intercore parallelization via MPI; requires OpenMPI/impi) Note: sbatch does not launch the tasks, only requests resources and submits the batch script. ntasks tells Slurm that N parallel tasks will be launched and to allocate resources accordingly. Parallel tasks are launched by the script.
#SBATCH --mem=24G
#SBATCH --job-name entropy
#SBATCH --account=rushworth
#SBATCH --partition=kingspeak-shared
#SBATCH --output=ent.out
#SBATCH --error=ent.error

entropy -i boech_gbs_retro_str_entropy_final.mgpl -m 1 -n 2 -k 2 -q retro_stricta_k2.txt -Q 1 -l 15000 -b 5000 -t 5 -s 30 -o mcmc_retro_stricta_m1n2k2Q1l15000b5000t5s30_qfilek2_chain1.hdf5
entropy -i boech_gbs_retro_str_entropy_final.mgpl -m 1 -n 2 -k 2 -q retro_stricta_k2.txt -Q 1 -l 15000 -b 5000 -t 5 -s 30 -o mcmc_retro_stricta_m1n2k2Q1l15000b5000t5s30_qfilek2_chain2.hdf5
entropy -i boech_gbs_retro_str_entropy_final.mgpl -m 1 -n 2 -k 2 -q retro_stricta_k2.txt -Q 1 -l 15000 -b 5000 -t 5 -s 30 -o mcmc_retro_stricta_m1n2k2Q1l15000b5000t5s30_qfilek2_chain3.hdf5

entropy -i boech_gbs_retro_retro_entropy_final.mgpl -m 1 -n 2 -k 2 -q retro_retro_k2.txt -Q 1 -l 15000 -b 5000 -t 5 -s 30 -o mcmc_retro_retro_m1n2k2Q1l15000b5000t5s30_qfilek2_chain1.hdf5
entropy -i boech_gbs_retro_retro_entropy_final.mgpl -m 1 -n 2 -k 2 -q retro_retro_k2.txt -Q 1 -l 15000 -b 5000 -t 5 -s 30 -o mcmc_retro_retro_m1n2k2Q1l15000b5000t5s30_qfilek2_chain2.hdf5
entropy -i boech_gbs_retro_retro_entropy_final.mgpl -m 1 -n 2 -k 2 -q retro_retro_k2.txt -Q 1 -l 15000 -b 5000 -t 5 -s 30 -o mcmc_retro_retro_m1n2k2Q1l15000b5000t5s30_qfilek2_chain3.hdf5


estpost.entropy -p q -s 4mcmc_retro_stricta_m1n2k2Q1l15000b5000t5s30_qfilek2_chain1.hdf5 mcmc_retro_stricta_m1n2k2Q1l15000b5000t5s30_qfilek2_chain2.hdf5 mcmc_retro_stricta_m1n2k2Q1l15000b5000t5s30_qfilek2_chain1.hdf5 -o ent_stricta_retro_convergence.txt


rule all:
    input:
        expand("mcmc_retro_stricta_m1n2k2Q1l15000b5000t5s30_qfilek2_chain{chain}_att2.hdf5", chain=range(1, 4))

rule entropy_retro_stricta:
    input:
        mgpl="boech_gbs_retro_stricta_gpl_entropy_final.mgpl",
        qfile="retro_stricta_k2.txt"
    output:
        "mcmc_retro_stricta_m1n2k2Q1l15000b5000t5s30_qfilek2_chain{chain}_att2.hdf5"
    params:
        m=1,
        n=2,
        k=2,
        Q=1,
        l=15000,
        b=5000,
        t=5,
        s=30
    shell:
        """
        entropy -i {input.mgpl} -m {params.m} -n {params.n} -k {params.k} -q {input.qfile} -Q {params.Q} -l {params.l} -b {params.b} -t {params.t} -s {params.s} -o {output}

