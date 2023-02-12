#!/bin/bash
#SBATCH --job-name=dnn_n10000_r1000_NP_e_0212
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128GB
#SBATCH --cluster=gpu
#SBATCH --gres=gpu:1
#SBATCH --time=4-10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alg223@pitt.edu

module purge
module load gcc/8.2.0
module load r/4.1.0
module load python/anaconda3.7-2019.03
module load cuda/11.0

Rscript code/run_sim_NP.R >& dnn_n10000_r1000_NP_e_0212.log
