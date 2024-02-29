#!/bin/bash -e
#SBATCH --account=nesi02659
#SBATCH --job-name=bayescan_starling.sl
#SBATCH --time=12:00:00
#SBATCH --mem=5GB
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --profile task

cd ~/outlier_analysis/analysis/bayescan

# Load bayescan
module purge
module load BayeScan/2.1-GCCcore-7.4.0

# Run bayescan
bayescan_2.1 ./starling_3populations.bs -od ./ -threads 8 -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 10
