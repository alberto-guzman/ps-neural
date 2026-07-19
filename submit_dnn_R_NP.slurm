#!/bin/bash
# Production job 3/3: the three keras neural networks (see code/run_sim_NP.R).
#
# !! UPDATE THE MODULE LOADS FOR YOUR CLUSTER before submitting (Pitt CRC
# versions circa 2023 below). This job additionally needs a python with
# TensorFlow reachable by reticulate — see the run_sim_NP.R header and
# infra/aws/bootstrap.sh for the working version pins; TF_USE_LEGACY_KERAS=1
# is required with TensorFlow >= 2.16.
#
# Resource notes from the 2026-07 AWS validation:
#   - The GPU is OPTIONAL: these are small MLPs (~7-14s/replication on one
#     CPU core; ~87 core-hours total). A CPU node works fine.
#   - TF survives SimDesign's fresh-process PSOCK workers, so switching
#     run_sim_NP.R to parallel = TRUE could cut wall time from ~87h serial
#     to a few hours — validate with a 2-replication test first; this job
#     could then move to the CPU cluster with more cores and no GPU.
#SBATCH --job-name=dnn_n10000_r1000_NP_v2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH --cluster=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=4-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alg223@pitt.edu

module purge
module load gcc/8.2.0
module load r/4.1.0
module load python/anaconda3.7-2019.03
module load cuda/11.0

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export TF_USE_LEGACY_KERAS=1

Rscript code/run_sim_NP.R >& n10000_r1000_NP_v2.log
