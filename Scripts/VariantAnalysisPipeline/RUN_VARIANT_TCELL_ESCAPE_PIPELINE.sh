#!/bin/sh

# sbatch commands
#SBATCH --partition=batch
#SBATCH --job-name=VARIANT_TCELL_ESCAPE
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --ntasks-per-node=8
#SBATCH --output=%j_%x.out
#SBATCH --error=%j_%x.err

## Script written by Paul R. Buckley, University of Oxford, 2022
module load R-cbrg/

# CORES must == --ntasks-per-node
CORES=8 # Number of cores to use for parallelization
HLAS_ANALYSIS_LOCATION=Nersisyan_HLA-I_common_alleles.txt # This is a list of HLA alleles to use for the analysis

# Need to change the below to the current directory locating the scripts
cd /scratch/COVID_19/ImpactCoV2MutationsTcells/Scripts/VariantAnalysisPipeline/

# 1) Compute RefSeq epitope coverage.
## Parameters
ANALYSIS_NAME=WUHAN # Name of the analysis
REF_SEQ_LOCATION=SARS-CoV-2-Wuhan.txt # Location of the reference sequence, in this case wuhan
INPUT_EPITOPES_DT=COV2_PEPTIDES_FOR_ANALYSIS.txt # Location of the WT epitopes for comparison
## RUN SCRIPT to determine where the above epitopes are located in RefSeq proteome
# This will output a RDS file with the coverage of the epitopes in the reference sequence
Rscript 1_compute_wuhan_epitope_coverage.R $INPUT_EPITOPES_DT $REF_SEQ_LOCATION $ANALYSIS_NAME

# 2) Compute the mutations for a VOC

EP_COVERAGE_RDS_LOCATION=${ANALYSIS_NAME}_COV2_Ep_Coverage.rds # Finds location of the RDS file from step 1
VARIANT_FASTA=Variants/SARS_CoV2_VOC_BA4.txt #Find the fasta file for the variant i.e BA4
VARIANT_NAME=BA4 # Name of the variant
# This will output a RDS file with the mutations for the variant in epitope positions as above, wrt to RefSeq
Rscript 2_compute_mutations_voc.R $EP_COVERAGE_RDS_LOCATION $VARIANT_FASTA $VARIANT_NAME $CORES

# 3) Compute binding of WT peptide and Mutant

module load netMHCpan/4.1b
# This will output a RDS file with the binding of the WT and Mutant peptides to 62 HLA alleles
Rscript 3_predict_antigen_presentation.R ${VARIANT_NAME}_EPITOPE_MUTATIONS.rds $VARIANT_NAME $HLAS_ANALYSIS_LOCATION



