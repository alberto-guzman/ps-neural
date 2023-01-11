#!/bin/bash
#SBATCH --job-name=dnn_sim
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=180GB
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alg223@pitt.edu

module purge
module load gcc/8.2.0
module load r/4.1.0

Rscript code/run_sim.R >& R.log

