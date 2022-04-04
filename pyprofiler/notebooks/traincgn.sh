#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=150G
#SBATCH --time=08:00:00
#SBATCH --job-name=traincgn
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --gres-flags=enforce-binding
#SBATCH --mem 150G


# load environment, e.g. set virtualenv, environment variables, etc

source /scratch/dmoi/miniconda/etc/profile.d/conda.sh
conda activate ML2


# Run Jupyter
python train_CGN.py
