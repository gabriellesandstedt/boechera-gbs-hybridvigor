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


estpost.entropy -p q -s 4 mcmc_retro_stricta_m1n2k2Q1l15000b5000t5s30_qfilek2_chain1.hdf5 mcmc_retro_stricta_m1n2k2Q1l15000b5000t5s30_qfilek2_chain2.hdf5 mcmc_retro_stricta_m1n2k2Q1l15000b5000t5s30_qfilek2_chain1.hdf5 -o ent_stricta_retro_convergence.txt
estpost.entropy -p q -s 0 mcmc_retro_stricta_m1n2k2Q1l15000b5000t5s30_qfilek2_chain1.hdf5 mcmc_retro_stricta_m1n2k2Q1l15000b5000t5s30_qfilek2_chain2.hdf5 mcmc_retro_stricta_m1n2k2Q1l15000b5000t5s30_qfilek2_chain3.hdf5 -o ent_stricta_retro_admix_prop.txt
estpost.entropy -p Q -s 0 mcmc_retro_stricta_m1n2k2Q1l15000b5000t5s30_qfilek2_chain1.hdf5 mcmc_retro_stricta_m1n2k2Q1l15000b5000t5s30_qfilek2_chain2.hdf5 mcmc_retro_stricta_m1n2k2Q1l15000b5000t5s30_qfilek2_chain3.hdf5 -o ent_stricta_retro_inter_anc.txt

estpost.entropy -p q -s 4 mcmc_retro_retro_m1n2k2Q1l15000b5000t5s30_qfilek2_chain1.hdf5 mcmc_retro_retro_m1n2k2Q1l15000b5000t5s30_qfilek2_chain2.hdf5 mcmc_retro_retro_m1n2k2Q1l15000b5000t5s30_qfilek2_chain1.hdf5 -o ent_retro_retro_convergence.txt
estpost.entropy -p q -s 0 mcmc_retro_retro_m1n2k2Q1l15000b5000t5s30_qfilek2_chain1.hdf5 mcmc_retro_retro_m1n2k2Q1l15000b5000t5s30_qfilek2_chain2.hdf5 mcmc_retro_retro_m1n2k2Q1l15000b5000t5s30_qfilek2_chain3.hdf5 -o ent_retro_retro_admix_prop.txt
estpost.entropy -p Q -s 0 mcmc_retro_retro_m1n2k2Q1l15000b5000t5s30_qfilek2_chain1.hdf5 mcmc_retro_retro_m1n2k2Q1l15000b5000t5s30_qfilek2_chain2.hdf5 mcmc_retro_retro_m1n2k2Q1l15000b5000t5s30_qfilek2_chain3.hdf5 -o ent_retro_retro_inter_anc.txt



