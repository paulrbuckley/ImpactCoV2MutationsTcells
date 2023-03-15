#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least 3 argument: if not, return an error
if (!length(args)==3) {
  stop("At least three argument must be supplied: 1) Lengths to process: 9 OR 10 OR 9,10 , 2) The analysis to run: either Mutation, ChangeFrom or ChangeTo, 3) The mutation to analyse in style e.g., Y_8_K  ", call.=FALSE)
}

### Sample commands
#Rscript RUN_MUTATION_IMPACT_TOOL.R 9,10 Mutation K_4_N
#Rscript RUN_MUTATION_IMPACT_TOOL.R 9 ChangeFrom K_4
#Rscript RUN_MUTATION_IMPACT_TOOL.R 9 ChangeTo 4_N

library(dplyr)
library(stringr)
library(data.table)
library(tidyr)
library(caret)
library(purrr)
library(doParallel)
library(foreach)
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

FULL_PREDICTIONS_DT = readRDS(file="INSILICO_MUTAGENESIS_COMPILED_PREDICTIONS_SUBONLY_w_MUTATIONDATA.rds") # Read in data
FULL_PREDICTIONS_DT = FULL_PREDICTIONS_DT%>% mutate(Mutation=paste0(ChangeFrom,"_",  SeqMutationPos,"_",ChangeTo)) # Create mutation column

# Clean input data
FULL_PREDICTIONS_DT=FULL_PREDICTIONS_DT %>% mutate(logsOdds = log(Prediction/WuhanScore)) # Create log odds column
FULL_PREDICTIONS_DT = FULL_PREDICTIONS_DT %>% mutate(Length=nchar(Peptide)) # Create length column

FULL_PREDICTIONS_DT=FULL_PREDICTIONS_DT %>% mutate(ChangeFrom = paste0(ChangeFrom,"_",SeqMutationPos)) # Create ChangeFrom column
FULL_PREDICTIONS_DT=FULL_PREDICTIONS_DT %>% mutate(ChangeTo = paste0(SeqMutationPos,"_",ChangeTo)) # Create ChangeTo column
FULL_PREDICTIONS_DT = FULL_PREDICTIONS_DT %>% mutate(Length = nchar(Peptide)) # Create length column
FULL_PREDICTIONS_DT=FULL_PREDICTIONS_DT %>%ungroup() # Ungroup data

# ARGUMENTS AND FILTERS
LENGTH_FILTER =as.numeric(str_split(args[1], ",")[[1]]) # Either 9 OR 10 OR 9,10

#ARGS = c("ChangeFrom") # Either 'Mutation', 'ChangeFrom' or 'ChangeTo'
ARGS = args[2] # Either 'Mutation', 'ChangeFrom' or 'ChangeTo'
#MUTATION = "Y_8"
MUTATION = args[3]

if(grepl("_",MUTATION) == FALSE){
  stop("MUTATION DOES NOT HAVE AN UNDERSCORE. PLEASE PUT AN UNDERSCORE BETWEEN MUTATION ELEMENTS E.G., K_2_A")
}

# Check amino acid
AA_CHECK_VECTOR=gsub("\\_|[0-9]","",MUTATION) # Remove underscores and numbers
AA_CHECK_VECTOR=strsplit(AA_CHECK_VECTOR,"")[[1]] # Split into vector
# Check if amino acid is valid
VALIDSEQS=foreach(i=1:length(AA_CHECK_VECTOR),.combine = "rbind", .packages = c("data.table","dplyr","protr"))%do% {
  PEPTIDE=AA_CHECK_VECTOR[i]
  VALID = protr::protcheck(PEPTIDE)
  if(VALID==FALSE){
    stop("LETTER ENTERED IN MUTATION IS NOT AN AMINO ACID")
  }
}

# Apply length filer
FULL_PREDICTIONS_DT=FULL_PREDICTIONS_DT %>% filter(Length %in% LENGTH_FILTER)
# Group by 'ARGS', calculate mean logsOdds, standard error, etc.
SUMMARY_DT=summarySE(FULL_PREDICTIONS_DT, measurevar="logsOdds", groupvars=ARGS)%>% mutate(PveNve = ifelse(logsOdds>0,"Positive","Negative"))

# Round to three decimal places
SUMMARY_DT=SUMMARY_DT %>% mutate_if(is.numeric, round, 3)
# Filter the first column [which is our ARGS column' by the provided mutation
REQUESTED_MUT_RESULTS=SUMMARY_DT %>% filter(.[[1]] == MUTATION)
# Print results to screen

# Arguments
print(paste0("LENGTH(S) ANALYSED: ", paste0(LENGTH_FILTER, collapse = ",")))
print(paste0("Analysis run: ", ARGS))
print(paste0("Mutation analysed ", MUTATION))

if(nrow(REQUESTED_MUT_RESULTS) == 0){
  stop(paste0("MUTATION NOT FOUND. It is likely that this mutation has not been analysed yet at all, or it hasn't been analysed for the particular peptiude length, or the ChangeFrom (X_) and ChangeTo (_X) are the wrong way round"))
}
print(paste0(ARGS, " ", MUTATION, " has a logs odds score of ", REQUESTED_MUT_RESULTS$logsOdds, " with standard error of ",REQUESTED_MUT_RESULTS$se, " from ", REQUESTED_MUT_RESULTS$N, " observations"))
print(paste0(ARGS, " ", MUTATION, " is estimated to have an average ", REQUESTED_MUT_RESULTS$PveNve, " impact on immunogenicity"))




















