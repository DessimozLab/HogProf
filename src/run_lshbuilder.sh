#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=50G
#SBATCH --time=08:00:00
#SBATCH --job-name=python


source /work/FAC/FBM/DBC/cdessim2/default/dmoi/miniconda3/etc/profile.d/conda.sh 
#conda init bash
conda activate pyprofiler_test

python lshbuilder.py --OMA  /work/FAC/FBM/DBC/cdessim2/default/aaltenho/birds/OmaStandalone_378_v3/OmaServer.h5 --mastertree /work/FAC/FBM/DBC/cdessim2/default/smajidi1/bird_project/birds370_iqtree_treefile_95bootstrap_internal_name_6let_16Nov.nwk 

