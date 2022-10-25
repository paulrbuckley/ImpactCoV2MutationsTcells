#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# This script takes a set of epitopes, and compares them to a proteome, to identify originating locatins.
# For example, we take a set of SARS-CoV-2 epitopes and compare them to the Wuhan reference sequence.
# These 'CoV2 ep coverage' locations are used subsequently in the pipeline for comparison to a variant.

# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop("At least two arguments must be supplied: 1) Path to data table with one column listing epitopes 2) path to a fasta proteome sequence to map the epitopes to e.g., Wuhan reference (fasta should contain protein names, preceded by protein=), 3 ) a naming convention for the wild type, e.g., 'Wuhan'.", call.=FALSE)
}

# Packages
library(data.table)
library(dplyr)
library(tidyverse)
library(foreach)
require(Biostrings)

full_peptides = fread(args[1]) %>% pull
#full_peptides = fread("COV2_PEPTIDES_FOR_ANALYSIS.txt")%>% pull

sequence_stringset = readAAStringSet(args[2],use.names = TRUE)
#sequence_stringset = readAAStringSet("SARS-CoV-2-Wuhan.txt",use.names = TRUE)

OUTPUT_NAME = args[3]
#OUTPUT_NAME = "WUHAN"

print(paste0("There are ",length(full_peptides), " epitopes for analysis"))
print(paste0("Reference stringset has ",length(sequence_stringset), " proteins"))

# Map epitopes to proteome sequence [i.e Wuhan isolate proteome]
Ep_Coverage = foreach(i = 1:length(full_peptides),.combine = rbind, .packages = c("dplyr","magrittr","data.table","Biostrings"))%do% {
  peptide_i=AAString(full_peptides[i])
  
  DT=as.data.table(unlist(vmatchPattern(peptide_i,sequence_stringset)))%>% mutate(Peptide = toString(peptide_i))
  DT
  
}

# Cleaning: Extract the protein name 
Ep_Coverage=Ep_Coverage %>% mutate(Protein = gsub(".*protein=(.+)\\].*","\\1",names))%>% mutate(Protein = gsub("\\].*","",Protein))
# Cleaning, rename start and end pos.
Ep_Coverage=Ep_Coverage %>% dplyr::rename(start_pos="start",end_pos="end")

# Compute all positions covered by a particular epitope and add as an additional column.
Protein_positions=Ep_Coverage %>% group_by(r=row_number()) %>% mutate(allpos = list(start_pos:end_pos)) %>% unnest() %>% group_by(start_pos,end_pos) %>% dplyr::summarise(Positions = paste0(allpos,collapse = ","))%>% ungroup
WUHAN_COV2_Ep_Coverage=Ep_Coverage %>% inner_join(Protein_positions)

# Give the reference sequence proteins their individual IDs and add as an additional column.
PROT_IDS=names(sequence_stringset) %>% as.data.table() %>% dplyr::rename(names=".") %>% mutate(RefSeqProteinID = row_number())%>% mutate(Protein = gsub(".*protein=(.+)\\].*","\\1",names))%>% mutate(Protein = gsub("\\].*","",Protein))

WUHAN_COV2_Ep_Coverage=WUHAN_COV2_Ep_Coverage %>% inner_join(PROT_IDS)

print(paste0("Found ", nrow(WUHAN_COV2_Ep_Coverage), " matching positions in the sequence. Please use your judgement to determine whether this is adequate."))

N_PEPTIDES_MAPPED=WUHAN_COV2_Ep_Coverage%>% select(Peptide) %>% distinct()%>% nrow

print(paste0(N_PEPTIDES_MAPPED, " distinct epitopes were successfully mapped to the reference sequence."))
print(paste0(length(full_peptides)-N_PEPTIDES_MAPPED, " epitopes have not been successfully mapped."))

if(N_PEPTIDES_MAPPED < length(full_peptides)){
  warning(paste0(length(full_peptides)-N_PEPTIDES_MAPPED, " epitopes have not been successfully mapped. This may be e.g., that these epitopes are from a different isolate than your reference sequeunce. For SARS-CoV-2 perhaps, this may mean the epitope originated in a different variant to the provided reference sequence."))
}

print(paste0("Outputting data to RDS file.."))
saveRDS(Ep_Coverage, file=paste0(OUTPUT_NAME,"_COV2_Ep_Coverage.rds"))







