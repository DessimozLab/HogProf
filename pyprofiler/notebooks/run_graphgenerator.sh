#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --job-name=python

# load environment, e.g. set virtualenv, environment variables, etc
source /scratch/dmoi/miniconda/etc/profile.d/conda.sh
conda activate ML2
# Run Jupyter

python graph_generator.py
