#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least 4 argument: if not, return an error
if (length(args)!=4) {
  stop("At least four arguments must be supplied: 1) Output of script 1: an reference sequence 'epitope coverage' RDs file. 2) a fasta file containing protein sequence of the VOC, 3) Name of the VOC, 4 numbers of cores to use for parallelised analysis. ", call.=FALSE)
}

library(data.table)
library(dplyr)
library(tidyr)
library(ggpubr)
library(Biostrings)
library(tidyverse)
library(foreach)
library(stringdist)
library(doParallel)

# important functions
calculateMutation_position_V2 <- function(str1, str2,start_pos) {
  require(Biostrings)
  peptide = str1
  variantalignment = str2
  alignment = pairwiseAlignment(variantalignment,peptide)
  str1=as.character(aligned(pattern(alignment)))
  str2=as.character(aligned(subject(alignment)))
 
  str1vec <- unlist(strsplit(str1, ""))
  str2vec <- unlist(strsplit(str2, ""))
  iMut <- (1:length(str1vec))[str1vec != str2vec]
  MUTPOS = (start_pos-1)+iMut
  return(paste0(MUTPOS))
}

# Read in 'RefSeq' environment
#Ep_Coverage=readRDS(file="WUHAN_COV2_Ep_Coverage.rds")
Ep_Coverage=readRDS(file=args[1])
Ep_Coverage_saved =Ep_Coverage

VARIANT = args[3]
#VARIANT="BA2"

#FASTA_LOCATION = "Variants/SARS_CoV2_VOC_BA2_Genbank_Protein_sequence.fasta.txt"
FASTA_LOCATION=args[2]
print(paste0("Location of sequence for variant ",VARIANT,": ", FASTA_LOCATION))

#cores=args[4]
cores=8

# Read in variant sequence
Variant_StringSet=readAAStringSet(FASTA_LOCATION,use.names = TRUE)

# Compute some IDs for the variant. These are used for mapping between Variant and RefSeq proteins.
Variant_ProteinIDs = names(Variant_StringSet) %>% as.data.table() %>% dplyr::rename(Variantnames=".") %>% mutate(VariantProteinID = row_number())%>% mutate(Protein = gsub(".*protein=(.+)\\].*","\\1",Variantnames))%>% mutate(Protein = gsub("\\].*","",Protein))
Ep_Coverage=Ep_Coverage_saved %>% inner_join(Variant_ProteinIDs)

# Perform local-global sequence alignment to entire variant

VARIANT_PROTEINS = Ep_Coverage %>% select(Protein, VariantProteinID) %>% distinct()
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)

Mutations_Epitopes = foreach(i = 1:nrow(Ep_Coverage),.combine = rbind, .packages = c("dplyr","magrittr","data.table","Biostrings")) %:%
  foreach(ii = 1:length(Variant_StringSet),.combine = rbind, .packages = c("dplyr","magrittr","data.table","Biostrings"))%dopar% {
    peptide_i=Ep_Coverage$Peptide[i]
    AlignmentDT = NULL
    ALIGNMENT=pairwiseAlignment(Variant_StringSet[[ii]], peptide_i, type="local-global")
    
    
    FullAlignmentDT = data.table(Peptide = peptide_i, VariantAlignment =   as.character(aligned(pattern(ALIGNMENT))), score=score(ALIGNMENT), TestedProtein=VARIANT_PROTEINS %>% filter(VariantProteinID == ii) %>% pull(Protein),VariantProteinID=VARIANT_PROTEINS %>% filter(VariantProteinID == ii) %>% pull(VariantProteinID) )
    
  }

stopCluster(cl)
# Take the alignment with the max score.
Mutations_Epitopes = Mutations_Epitopes%>% distinct() %>% group_by(Peptide) %>% slice_max(score, with_ties = TRUE)%>% ungroup
#Clean the data
Mutations_Epitopes=Mutations_Epitopes %>% dplyr::rename(Protein=TestedProtein) %>% select(Peptide, VariantAlignment, VariantProteinID, Protein)
Mutations_Epitopes = Mutations_Epitopes %>% distinct() %>% ungroup
# Calculate hamming distance between the RefSeq 'Peptide' and the 'VariantAlignmet' epitope from the variant sequence.
Mutations_Epitopes=Mutations_Epitopes %>% mutate(Hamming = stringdist::stringdist(Peptide, VariantAlignment, method="hamming"))
# Join the REfSeq Ep Coverage dataset with the VariantAlignment list of mutations
MutationData_RefSeq_vs_Variant=Ep_Coverage %>% inner_join(Mutations_Epitopes)

# To identify only the set with mutations, exclude any exact matches/
MUTATIONDATA=MutationData_RefSeq_vs_Variant %>% filter(! Hamming==0)
# Compute additional information about the mutation.
MUTATIONS = foreach(i = 1:nrow(MUTATIONDATA),.combine = rbind, .packages = c("dplyr","magrittr","data.table","Biostrings"))%do% {
  
  mutations = paste0(calculateMutation_position_V2(MUTATIONDATA$Peptide[i],MUTATIONDATA$VariantAlignment[i],MUTATIONDATA$start_pos[i]),collapse = ",")
  data.table(Peptide = MUTATIONDATA$Peptide[i], VariantAlignment=MUTATIONDATA$VariantAlignment[i], start_pos= MUTATIONDATA$start_pos[i], Protein = MUTATIONDATA$Protein[i], MutationPos = mutations)
  
}

#Below line assumes that mutatations in ORF1a and ORF1AB in the same position with same peptuide, are the same mutation thus are joined to be ORF1AB/ORF1a
MUTATIONS=MUTATIONS %>% group_by(Peptide,VariantAlignment,start_pos,MutationPos) %>% dplyr::summarise(Protein = sort(paste0(Protein,collapse = "/")))
print(paste0(VARIANT, " computation complete."))
print("NOTE: Mutation has not yet been mapped to known mutation list.")
saveRDS(MUTATIONS,paste0(VARIANT,"_EPITOPE_MUTATIONS.rds"))














