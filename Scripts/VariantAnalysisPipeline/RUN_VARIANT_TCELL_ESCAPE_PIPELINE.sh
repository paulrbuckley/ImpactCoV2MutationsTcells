#!/bin/sh

# sbatch commands
#SBATCH --partition=batch
#SBATCH --job-name=VARIANT_TCELL_ESCAPE
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --ntasks-per-node=8
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err

module load R-cbrg/

# CORES must == --ntasks-per-node
CORES=8
HLAS_ANALYSIS_LOCATION=Nersisyan_HLA-I_common_alleles.txt

# Need to change the below to this directory
cd /t1-data/project/koohylab/pbuckley/COVID_19/OMICRON/Variant_Tcell_escape_pipeline/

# 1) Compute RefSeq epitope coverage.
## Parameters
ANALYSIS_NAME=WUHAN
REF_SEQ_LOCATION=SARS-CoV-2-Wuhan.txt
INPUT_EPITOPES_DT=COV2_PEPTIDES_FOR_ANALYSIS.txt
## RUN SCRIPT
#Rscript 1_compute_wuhan_epitope_coverage.R $INPUT_EPITOPES_DT $REF_SEQ_LOCATION $ANALYSIS_NAME

# 2) Compute the mutations for a VOC

EP_COVERAGE_RDS_LOCATION=${ANALYSIS_NAME}_COV2_Ep_Coverage.rds
VARIANT_FASTA=Variants/SARS_CoV2_VOC_BA4.txt
VARIANT_NAME=BA4

Rscript 2_compute_mutations_voc.R $EP_COVERAGE_RDS_LOCATION $VARIANT_FASTA $VARIANT_NAME $CORES

# 3) Compute binding of WT ppeptide and Mutant

module load netMHCpan/4.1b

Rscript 3_predict_antigen_presentation.R ${VARIANT_NAME}_EPITOPE_MUTATIONS.rds $VARIANT_NAME $HLAS_ANALYSIS_LOCATION

# 4) optional: integrate mutation data with known SNPs for a variant. Just need a list: e.g., OMICRON_SNPS.csv



