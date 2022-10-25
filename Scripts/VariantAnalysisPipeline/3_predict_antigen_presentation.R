#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least 1 argument: if not, return an error
if (length(args)!=3) {
  stop("At least 3 arguments must be supplied: 1) A mutations RDs file. 2) Name of the VOC, 3) file with a list of HLA", call.=FALSE)
}

library(data.table)
library(dplyr)
library(tidyr)
library(Biostrings)
library(foreach)
library(doParallel)
library(stringdist)
library(purrr)
# Read in the Full peptide dataset


MUTATIONS_LIST=args[1]
#MUTATIONS_LIST="BA2_EPITOPE_MUTATIONS.rds"

VARIANT = args[2]
#VARIANT="BA2"

#HLA_LOCATION = "Nersisyan_HLA-I_common_alleles.txt"
HLA_LOCATION = args[3]
ALLELES.TO.TEST = fread(HLA_LOCATION,header = FALSE)%>% pull(V1)

FullDataset=readRDS(MUTATIONS_LIST)
FullDataset=FullDataset %>% ungroup()

FILEPATH=paste0(VARIANT,"_WT_NETMHCPAN_CLASSI/")
dir.create(paste0(VARIANT,"_WT_NETMHCPAN_CLASSI/"))
TEST_DATA_LOCATION = FILEPATH

# Compute for WT
PEPTIDES = FullDataset%>% select(Peptide) %>% mutate(Length = width(Peptide)) %>% filter(!Length >14) %>% select(Peptide) %>% unique %>% pull

cores=8
cl <- makeCluster(cores[1]-4) #not to overload your computer
registerDoParallel(cl)

foreach(allele_i = 1:length(unique(ALLELES.TO.TEST)),.packages = c("dplyr","data.table")) %dopar% {
  HLA_ALLELE_FOR_TESTING = ALLELES.TO.TEST[allele_i]
  print(c("ALLELE ",allele_i))
  testdata=paste0(TEST_DATA_LOCATION,"Allele_",HLA_ALLELE_FOR_TESTING,"_NetMHC_data.txt")
  write.table(PEPTIDES,file=testdata,sep="\n",col.names = F,row.names = F,quote=F)
  # Run model
  RESULTS_OUTPUT = paste0(TEST_DATA_LOCATION,"Allele_",HLA_ALLELE_FOR_TESTING,"_NetMHC","_Results.csv")
  system(paste0("netMHCpan -BA -p ",testdata," -a ",HLA_ALLELE_FOR_TESTING," -xls -xlsfile ", RESULTS_OUTPUT))
  
  
}

stopCluster(cl)

# WT read
data_path <- TEST_DATA_LOCATION
files <- dir(data_path, pattern = "NetMHC_Results.csv")

data3 <- data_frame(file = files) %>%
  mutate(file_contents = map(file,
                             ~ fread(file.path(data_path, .),skip = 1))
  )

Netmhcpanres <- unnest(data3)
Netmhcpanres=Netmhcpanres %>% mutate(HLA_Allele = gsub(x=Netmhcpanres$file,pattern="Allele_|_NetMHC_Results.csv",replacement = ""))
Netmhcpanres=Netmhcpanres %>% select(! c(file,Pos,ID,core,icore)) %>% mutate(Binder = ifelse(NB==1,"BINDER","NONBINDER"))
Netmhcpanres=Netmhcpanres %>% mutate(HLA_Allele = gsub(x=HLA_Allele,pattern = "\\:",replacement = ""))

NETMHC_CLASSI=Netmhcpanres
NETMHC_CLASSI=NETMHC_CLASSI %>% mutate(nM_41 = round(50000^(1-`BA-score`) ,digits=5))
NETMHC_CLASSI=NETMHC_CLASSI %>% mutate(Binder = ifelse(BA_Rank>2, "NONBINDER","BINDER"))
NETMHC_CLASSI %>% head
NETMHC_CLASSI= NETMHC_CLASSI %>% group_by(Peptide) %>% dplyr::rename(Predicted_Binding=HLA_Allele)%>% select(!c(NB)) %>% summarise(Predicted_Binding = paste0(unique(Predicted_Binding),collapse = ","),
                                                                                                                                   `EL-score`=paste0(`EL-score`,collapse = ","),
                                                                                                                                   `EL_Rank`=paste0(`EL_Rank`,collapse = ","),
                                                                                                                                   `BA-score`=paste0(`BA-score`,collapse = ","),
                                                                                                                                   `BA_Rank`=paste0(`BA_Rank`,collapse = ","),
                                                                                                                                   `Ave`=paste0(`Ave`,collapse = ","),
                                                                                                                                   `Binder` = paste0(`Binder`,collapse = ","),
                                                                                                                                   nM_41 = paste0(`nM_41`, collapse = ","))

saveRDS(NETMHC_CLASSI,file=paste0(VARIANT,"_WT_BINDING_SCORES_NETMHCPAN.rds"))

# Compute for MT
PEPTIDES = FullDataset%>% select(VariantAlignment)%>%
  dplyr::rename(Peptide=VariantAlignment)%>% mutate(Length = width(Peptide)) %>% filter(!Length >14) %>% select(Peptide) %>% unique %>% pull

FILEPATH=paste0(VARIANT,"_MT_NETMHCPAN_CLASSI/")
dir.create(paste0(VARIANT,"_MT_NETMHCPAN_CLASSI/"))
TEST_DATA_LOCATION = FILEPATH

cores=8
cl <- makeCluster(cores[1]-4) #not to overload your computer
registerDoParallel(cl)

foreach(allele_i = 1:length(unique(ALLELES.TO.TEST)),.packages = c("dplyr","data.table")) %dopar% {
  HLA_ALLELE_FOR_TESTING = ALLELES.TO.TEST[allele_i]
  print(c("ALLELE ",allele_i))
  testdata=paste0(TEST_DATA_LOCATION,"Allele_",HLA_ALLELE_FOR_TESTING,"_NetMHC_data.txt")
  write.table(PEPTIDES,file=testdata,sep="\n",col.names = F,row.names = F,quote=F)
  # Run model
  RESULTS_OUTPUT = paste0(TEST_DATA_LOCATION,"Allele_",HLA_ALLELE_FOR_TESTING,"_NetMHC","_Results.csv")
  system(paste0("netMHCpan -BA -p ",testdata," -a ",HLA_ALLELE_FOR_TESTING," -xls -xlsfile ", RESULTS_OUTPUT))
  
  
}

stopCluster(cl)

data_path <- TEST_DATA_LOCATION
files <- dir(data_path, pattern = "NetMHC_Results.csv")

data3 <- data_frame(file = files) %>%
  mutate(file_contents = map(file,
                             ~ fread(file.path(data_path, .),skip = 1))
  )

Netmhcpanres <- unnest(data3)
Netmhcpanres=Netmhcpanres %>% mutate(HLA_Allele = gsub(x=Netmhcpanres$file,pattern="Allele_|_NetMHC_Results.csv",replacement = ""))
Netmhcpanres=Netmhcpanres %>% select(! c(file,Pos,ID,core,icore)) %>% mutate(Binder = ifelse(NB==1,"BINDER","NONBINDER"))
Netmhcpanres=Netmhcpanres %>% mutate(HLA_Allele = gsub(x=HLA_Allele,pattern = "\\:",replacement = ""))

NETMHC_CLASSI=Netmhcpanres
NETMHC_CLASSI=NETMHC_CLASSI %>% mutate(MT_nM_41 = round(50000^(1-`BA-score`) ,digits=5))
NETMHC_CLASSI=NETMHC_CLASSI %>% mutate(Binder = ifelse(BA_Rank>2, "NONBINDER","BINDER"))

NETMHC_CLASSI= NETMHC_CLASSI %>% dplyr::rename(VariantAlignment = Peptide)%>% group_by(VariantAlignment) %>%
  dplyr::rename(MT_Predicted_Binding=HLA_Allele,
                `MT_EL-score`=`EL-score`,
                `MT_EL_Rank`=`EL_Rank`,
                `MT_BA-score`=`BA-score`,
                `MT_BA_Rank`=`BA_Rank`,
                `MT_Ave`=`Ave`,
                `MT_Binder` = `Binder`) %>% select(!c(NB)) %>%
  summarise(MT_Predicted_Binding = paste0(unique(MT_Predicted_Binding),collapse = ","),
            `MT_EL-score`=paste0(`MT_EL-score`,collapse = ","),
            `MT_EL_Rank`=paste0(`MT_EL_Rank`,collapse = ","),
            `MT_BA-score`=paste0(`MT_BA-score`,collapse = ","),
            `MT_BA_Rank`=paste0(`MT_BA_Rank`,collapse = ","),
            `MT_Ave`=paste0(`MT_Ave`,collapse = ","),
            `MT_Binder` = paste0(`MT_Binder`,collapse = ","),
            MT_nM_41 = paste0(MT_nM_41, collapse = ","))

saveRDS(NETMHC_CLASSI,file=paste0(VARIANT,"_MT_BINDING_SCORES_NETMHCPAN.rds"))





