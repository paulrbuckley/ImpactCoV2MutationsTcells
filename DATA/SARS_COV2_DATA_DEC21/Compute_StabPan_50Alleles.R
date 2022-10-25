# Title     : TODO
# Objective : TODO
# Created by: paulbuckley
# Created on: 22/12/2021

library(data.table)
library(dplyr)
library(tidyr)
library(Biostrings)
library(foreach)
library(doParallel)
library(stringdist)

# Read in the Full peptide dataset


#FullDataset=readRDS("DEC_2021_Human_Cleaned_SarsCoV2_EpitopeDataset_w_Freq_ONLYPOSITIVE.rds")
FullDataset =readRDS("IEDB_VIPR_MIRA_COV2_DATASET.rds")


# FullDataset=FullDataset %>% mutate(HLA_Allele = replace_na(""))

# for Frequency_of_ImmunogenicityBinary column, it represents the frequency of the observations of immunoenicity bianry. So there may be more than two observations gowever, if the qualititative measure of the peptide has a positive-X observatiobn. But still, the freq column still represents the frequency of pve/nve for said peptide.

#CLASSI.ALLELES=c('HLA-A0101','HLA-A0201','HLA-A0301','HLA-A2402','HLA-B0702','HLA-B0801','HLA-B4001','HLA-C0401','HLA-C0701','HLA-C0702')
#APPEND.ALLELES = readRDS("../SARS_COV2_DATA_DEC21/DEC_2021_Human_Cleaned_SarsCoV2_EpitopeDataset_w_Freq_ONLYPOSITIVE.rds")%>% select(HLA_Allele) %>% separate_rows(HLA_Allele,sep=",")%>%
# distinct() %>% filter(!HLA_Allele %in% c("HLA class I","HLA class II","HLA-DR","-N/A-","","HLA-A2","HLA-A29","HLA-B7","HLA-B35","HLA-A29","HLA-A24"))%>%
# mutate(HLA_Allele = gsub("\\*|\\:","",HLA_Allele))%>% filter(!HLA_Allele %in% grep("HLA-D",HLA_Allele,value=T))%>% unique %>% pull

#APPEND.ALLELES=fread("~/Nexus365/WIMM CCB - Koohy Group - Koohy Group/standardData/standardDataset_200903_PB_Inc_ImmunogenicityEvidence.csv")%>% filter(MHCType=="MHC-I") %>% separate_rows(MHCrestrictions,sep="\\|")%>%
 # select(MHCrestrictions) %>% table %>% data.table() %>% arrange(desc(N))%>% filter(. %in% grep("[0-9]",.,value = T))%>% head(80)%>% select(".") %>% dplyr::rename(HLA_Allele = ".")%>%
  #mutate(HLA_Allele =gsub("\\:",replacement = "",HLA_Allele))%>% mutate(HLA_Allele =gsub("\\*",replacement = "",HLA_Allele)) %>%
  #filter(!HLA_Allele %in% c("HLA class I","HLA class II","HLA-DR","-N/A-","","HLA-A2","HLA-A29","HLA-B7","HLA-B35","HLA-A29","HLA-A24","HLA-A1","HLA-B27","HLA-B51","HLA-A02","HLA-B8","HLA-B44",
   #                         "HLA-A3","HLA-A11","HLA-B57","JLA-A26","HLA-B53","HLA-B40","HLA-B62","HLA-Cw6","HLA-A03","HLA-A33","HLA-B18","HLA-B58","HLA-A68","HLA-A31","HLA-B15","HLA-B60",
    #                        "HLA-A26","HLA-A23","HLA-A25","HLA-B08","HLA-B07","HLA-B37","HLA-Cw4"))%>% pull


#ALLELES.TO.TEST = unique(c(CLASSI.ALLELES,APPEND.ALLELES))
ALLELES.TO.TEST = fread("../COVID_19_OMICRON_VOC_IMMUNE_ESCAPE/Nersisyan_HLA-I_common_alleles.txt", header = FALSE) %>% pull(V1)

PEPTIDES = FullDataset%>% select(Peptide) %>% mutate(Length = width(Peptide)) %>% filter(!Length >14) %>% select(Peptide) %>% unique
PEPTIDES=PEPTIDES %>% mutate(Length = nchar(Peptide))
FullDataset=PEPTIDES

TEST_DATA_LOCATION = "STABPAN/"

# Function to provide a closest match. Used to match HLA Alleles across mixed output styles.
ClosestMatch2 = function(string, stringVector){

  stringVector[amatch(string, stringVector, maxDist=Inf)]

}

ALLOWED_ALLELES=fread("STABPAN_MHC_allele_names.txt",header = F)

ALLELES.TO.TEST=ClosestMatch2(ALLELES.TO.TEST,ALLOWED_ALLELES$V1)

for(allele_i in 1:length(unique(ALLELES.TO.TEST))){


  #HLA_ALLELE_FOR_TESTING = gsub(x=unique(FullDataset$HLA_Allele)[allele_i],pattern=":",replacement = "")
  HLA_ALLELE_FOR_TESTING =unique(ALLELES.TO.TEST)[allele_i]
    HLA_FOR_FILENAME = gsub(":","",HLA_ALLELE_FOR_TESTING)

  DATASET = FullDataset
    for(i in 1:length(unique(DATASET$Length))){

      LENGTH = unique(DATASET$Length)[i]

  testdata=paste0(TEST_DATA_LOCATION,"LENGTH_",LENGTH,"_Allele_",HLA_FOR_FILENAME,"_NetMHC_data.txt")

  write.table(DATASET %>% filter(Length == LENGTH) %>% select(Peptide) %>% pull,file=testdata,sep="\n",col.names = F,row.names = F,quote=F)
  # Run model
  RESULTS_OUTPUT = paste0(TEST_DATA_LOCATION,"LENGTH_",LENGTH,"_Allele_",HLA_FOR_FILENAME,"_NetMHC","_Results.csv")

    system(paste0("/Applications/netMHCstabpan-1.0/netMHCstabpan -p ",testdata," -a ",HLA_ALLELE_FOR_TESTING," -l ",LENGTH," -xls -xlsfile ", RESULTS_OUTPUT))
    }
}
