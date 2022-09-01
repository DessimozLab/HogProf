#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=50G
#SBATCH --time=08:00:00
#SBATCH --job-name=python


source /work/FAC/FBM/DBC/cdessim2/default/dmoi/condaenvs/etc/profile.d/conda.sh
#conda init bash
conda activate ML2
python lshbuilder.py --OMA /work/FAC/FBM/DBC/cdessim2/default/dmoi/datasets/OMA/sep2022/OmaServer.h5 --dbtype all --name all

