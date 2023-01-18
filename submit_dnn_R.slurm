#!/bin/bash
#SBATCH --job-name=dnn_sim
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=500GB
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alg223@pitt.edu

module purge
module load gcc/8.2.0
module load r/4.1.0
module load python/anaconda3.7-2019.03

Rscript code/run_sim.R >& R.log
