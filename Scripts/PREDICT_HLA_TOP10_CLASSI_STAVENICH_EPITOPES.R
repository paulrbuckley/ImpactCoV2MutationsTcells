library(data.table)
library(dplyr)
library(tidyr)
library(Biostrings)
library(foreach)
library(doParallel)
# Read in the Full peptide dataset


FullDataset=readRDS("../OUTPUT/stavenichepitopes.rds")
FullDataset=FullDataset %>% select(Peptide)

ALLELES.TO.TEST = "HLA-A01:01"
length(ALLELES.TO.TEST)
PEPTIDES = FullDataset%>% select(Peptide) %>% mutate(Length = width(Peptide)) %>% filter(!Length >14) %>% select(Peptide) %>% unique %>% pull

TEST_DATA_LOCATION = "STAVENICH_EPS/"
dir.create(TEST_DATA_LOCATION)
cores=detectCores()
cl <- makeCluster(cores[1]-4) #not to overload your computer
registerDoParallel(cl)

foreach(allele_i = 1:length(unique(ALLELES.TO.TEST)),.packages = c("dplyr","data.table")) %dopar% {
  HLA_ALLELE_FOR_TESTING = ALLELES.TO.TEST[allele_i]
  HLA_FOR_FILENAME = gsub(":","",HLA_ALLELE_FOR_TESTING)
  print(c("ALLELE ",allele_i))
  testdata=paste0(TEST_DATA_LOCATION,"Allele_",HLA_FOR_FILENAME,"_NetMHC_data.txt")
  write.table(PEPTIDES,file=testdata,sep="\n",col.names = F,row.names = F,quote=F)
  # Run model
  RESULTS_OUTPUT = paste0(TEST_DATA_LOCATION,"Allele_",HLA_FOR_FILENAME,"_NetMHC","_Results.csv")
  system(paste0("/Applications/netMHCpan-4.1/netMHCpan -BA -p ",testdata," -a ",HLA_ALLELE_FOR_TESTING," -xls -xlsfile ", RESULTS_OUTPUT))


}

stopCluster(cl)


