#!/bin/sh

##Script assembled by Paul Buckley (University of Oxford)
##
## this is a slurm script to run the TITAN model against test dataset on the HPC

# sbatch commands
#SBATCH --job-name=TRAIN_TITAN_OTB_FINAL_TEST_OTB_POGORELY_MIRA_NOWUHANMUT_100K_NO_OMI_NVE
#SBATCH -p gpu
#SBATCH --gpus=1
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err

module load python-base/3.7.13

cd /project/koohylab/pbuckley/COVID_19/TITAN/TITAN-main/
source venv/bin/activate

python3.7 scripts/flexible_model_eval.py \
TEST_OTB_WUHAN_MUTATED/test.csv \
TEST_OTB_WUHAN_MUTATED/tcr.csv \
TEST_OTB_WUHAN_MUTATED/epitopes.smi \
TRAIN_TITAN_OTB_FINAL_TEST_OTB_POGORELY_MIRA_NOWUHANMUT_100K_NO_OMI_NVE/test_fine_tuned/OTB_POG \
bimodal_mca WUHAN_MUTATED
