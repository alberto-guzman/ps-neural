#!/bin/bash
# Production job 2/3: GBM alone (see code/run_sim_gbm.R header for why it is
# separate). Measured budget under validation-split selection (2026-07-19
# shakedown): ~410 core-hours => ~6.5h on 64 cores; the 3-day limit is ample.
#
# !! UPDATE THE MODULE LOADS FOR YOUR CLUSTER before submitting: the versions
# below are from the Pitt CRC circa 2023 and are almost certainly stale. The
# simulation needs R >= 4.5 with the packages listed in code/run_sim_gbm.R
# (no python/CUDA needed for this job). infra/aws/bootstrap.sh documents a
# complete from-scratch environment build.
#SBATCH --job-name=dnn_n10000_r1000_gbm_v2
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

# keep worker-internal threading at 1: SimDesign already uses all 64 cores,
# and code/02_analyse_fun.R additionally pins gbm's n.cores = 1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

Rscript code/run_sim_gbm.R >& dnn_n10000_r1000_gbm_v2.log
