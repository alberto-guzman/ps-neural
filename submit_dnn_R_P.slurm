#!/bin/bash
#SBATCH --job-name=dnn_n10000_r1000_P_lmSE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=250GB
#SBATCH --cluster=htc
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alg223@pitt.edu

module purge
module load gcc/8.2.0
module load r/4.1.0
module load python/anaconda3.7-2019.03
module load cuda/11.0

Rscript code/run_sim_P.R >& dnn_n10000_r1000_P_lmSE.log