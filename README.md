# Mutation impact tool
- Script 'RUN_MUTATION_IMPACT_TOOL.R' will run the tool.
-At least three argument must be supplied: 1) Lengths to process: 9 OR 10 OR 9,10 , 2) The analysis to run: either Mutation, ChangeFrom or ChangeTo, 3) The mutation to analyse in style e.g., Y_8_K
- This script will calculate the estimated impact - given our current mutational impact data - of a mutation of interest.
# Variant Analysis Pipeline
- Given a variant proteome sequence (fasta), the pipeline 'Scripts/VariantAnalysisPipeline' will compare known Wuhan T cell epitopes (or another set is provided and modified) with the variant of interest. Mutated epitopes in the variant will be identified and antigen presentation status will be predicted for these epitopes.
- The pipeline is run by running the script 'Scripts/VariantAnalysisPipeline/VariantAnalysisPipeline.sh'.
- This is a slurm batch script and arguments are supplied within the script.

# TRAP is supplied in another manuscript: see references for bioRxiv version.

# Analysis for manuscript

- 1_Compare_CoV2_Tcell_Targets_Omicron - this file curates a set of Wuhan
CD8T epitopes, and compares them to Omicron BA1 strain. This file also for 
BA1 generates WUHAN_PEPTIDE_BINDING_SCORES.rds and 
OMICRON_PETIDE_BINDING_SCORES.rds. These contain antigen presentation 
metrics for wuhan mutated and couinterpart epitopes.
Figure1.rmd - creates Figure 1. Analyses mutation statistics outpit from 
.rmd file above (1_*...rmd)
- 2_Omicron_Vs_Original_BindingRecognition_createAgretopcityDT - reads in 
*_BINDING_SCORES.rds files, integrates them so we have paired WT-MT 
observations with all WT and MT antigen presentaiton metrics from 
netmHCpan etc. Also creates an agretopicity column. Outputs an .rds file 
called 'AgretopicityDT.rds', summarising these data.
- 3_Analyse_AgPres_Subvariants.rmd - after the pipeline in 
Scripts/Variant_TcellEscape/ has been run for BA2, 4, 5, this notebook 
produces the subvariant plots for Figure 2 and Supplementary Figure 2, i.e 
the antigen presentation analysis for omICRON BA1 Subvariants
- 4_Omicron_HLASupertypesAnalysis - creates a suppleentary plot for BA1, 
where supertype is analysed for any associations with changes in binding 
affinity between wuhan and omicron/
Figure2.rmd - creates Figure 2. Analyses antigen presentation for BA1, and 
incorporates subvariant ag pres analysis from .rmd above prefixed 
3_*..rmd.
- 5_TRAPP_TITAN_Analyse_Wuhan_vs_Subvars.rmd - For each subvariant, analyse 
immunogenicity potential using TRAP. Ignore TITAN in name, TITAN not used 
there/
Figure3.rmd - creates Figure 3. Analyses T cell immunogenicity potential 
using TRAP for BA1. Then incorporates subvriant immunogenicity analysis 
from .rmd file prefixed 5_..rmd to create Fig3 and Supp Fig3. - Generates 
RI SCORES and pan HLA RI scores/
- 6_MIRA_TCR_ANALYSIS_V2 AND 7_MIRA_TCR_ANALYSIS_V2_SUBVARS - analysis that 
creates panels for Figure 4. Incorporates pan HLA RI scores with MIRA 
dataset for BA1 mutated epitopes and SUBVARS: for BA2 , 4, 5. 
- 8_GENERATE_MUTAGENESIS_PEPTIDES - Takes 9/10mers of interest from BA1, and 
generates every single amino acid variant. 
- 9_TRAPP_EscapePotential - takes output from 8_..rmd and analyses it. 
Creates neighbor network analysis and all other panels for Fig 5. 
