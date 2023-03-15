#!/bin/sh

##Script assembled by Paul Buckley (University of Oxford)
##
##

# sbatch commands
#SBATCH --job-name=TRAIN_TITAN_OTB_FINAL_TEST_OTB_POGORELY_MIRA_NOWUHANMUT_100K_NO_OMI_NVE
#SBATCH -p gpu
#SBATCH --gpus=1
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err

module load python-base/3.7.13

cd /project/koohylab/pbuckley/COVID_19/TITAN/TITAN-main/
source venv/bin/activate

python3.7 scripts/semifrozen_finetuning.py \
TRAIN_TITAN_OTB_FINAL_TEST_OTB_POGORELY_MIRA_NOWUHANMUT_100K_NO_OMI_NVE/train.csv \
TRAIN_TITAN_OTB_FINAL_TEST_OTB_POGORELY_MIRA_NOWUHANMUT_100K_NO_OMI_NVE/test.csv TRAIN_TITAN_OTB_FINAL_TEST_OTB_POGORELY_MIRA_NOWUHANMUT_100K_NO_OMI_NVE/tcr.csv TRAIN_TITAN_OTB_FINAL_TEST_OTB_POGORELY_MIRA_NOWUHANMUT_100K_NO_OMI_NVE/epitopes.smi \
public/trained_model TRAIN_TITAN_OTB_FINAL_TEST_OTB_POGORELY_MIRA_NOWUHANMUT_100K_NO_OMI_NVE/test_fine_tuned OTB_POG TRAIN_TITAN_OTB_FINAL_TEST_OTB_POGORELY_MIRA_NOWUHANMUT_100K_NO_OMI_NVE/finetune_params.json bimodal_mca
