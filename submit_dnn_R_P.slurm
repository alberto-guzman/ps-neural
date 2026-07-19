#!/bin/bash
# Production job 1/3: six conventional methods (see code/run_sim_P.R header).
# Measured budget: ~600 core-hours => ~10h on 64 cores; 3-day limit is ample.
#
# !! UPDATE THE MODULE LOADS FOR YOUR CLUSTER before submitting: the versions
# below are from the Pitt CRC circa 2023 and are almost certainly stale. The
# simulation needs R >= 4.5 with the packages listed in code/run_sim_P.R
# (no python/CUDA needed for this job). infra/aws/bootstrap.sh documents a
# complete from-scratch environment build.
#SBATCH --job-name=dnn_n10000_r1000_P_v2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250GB
#SBATCH --cluster=htc
#SBATCH --time=3-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alg223@pitt.edu

module purge
module load gcc/8.2.0
module load r/4.1.0

# keep worker-internal threading at 1: SimDesign already uses all 64 cores
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

Rscript code/run_sim_P.R >& dnn_n10000_r1000_P_v2.log
