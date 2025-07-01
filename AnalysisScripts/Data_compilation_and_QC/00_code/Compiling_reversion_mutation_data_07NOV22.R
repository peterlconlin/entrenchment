## Compiling data on the genome position and identity of reversion mutations at transition and final timepoints
## 07 Novemeber 2022
## Peter Conlin

# Libraries
library(dplyr)
library(vegan)

###-------------------------------------------------------
### COMBINE TRANS DATA FILES
# Set working directory
setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/entrenchment_data/Entrenchment_new_data/004e_new/unis_end/")

#reads in all csv files in current directory
rm(list = ls()) #Clear environment
file_list <- list.files(pattern='detail_uni.csv')

#make empty vectors to save data
propagule_genome <- c()
strategy_rep <- c()

#loop through files to extract propagule genome and ID info
for (i in 1:length(file_list)) {
  temp <- read.csv(file_list[1], header=TRUE, stringsAsFactors=FALSE)
  head(temp)
  (length(temp$uni)/100)/32
  propagule_genome <- c(propagule_genome, unique(temp$propagule_genome[temp$parent_id==0]))
  strategy_rep <- c(strategy_rep, rep(unique(temp$replicate[temp$parent_id==0]), 
                                      length(unique(temp$propagule_genome[temp$parent_id==0]))))
}
  
length(propagule_genome)
length(strategy_rep)


  
  assign(temp[i], read.table(file=temp[i], header=TRUE, stringsAsFactors=FALSE))

#makes a list of all dataframes in the workspace!
vars <- ls()
nvars <- length(vars)
dataframes <-list()
j <- 1
for(i in 1:nvars)
{
  if(class(get(vars[i]))=="data.frame")
  {
    dataframes[[j]] <- get(vars[i])
    names(dataframes)[j] <- vars[i]
    j <- j+1
  }
}

length(dataframes[[1]]$timepoint)
names(dataframes)

for (i in 1:length(dataframes)){
  dataframes[[i]] <- cbind.data.frame(dataframes[[i]], mc_or_uni = rep("uni", length(dataframes[[i]]$timepoint)))
}

head(dataframes[[i]])

details_trans_data <- bind_rows(dataframes, .id = "column_label")
head(details_trans_data)
details_trans_data$timepoint[details_trans_data$timepoint==60] <- 0
unique(details_trans_data$timepoint)
setwd("~/Desktop/Entrenchment_data_05MAY21/")
write.csv(details_trans_data, "details_trans_data_ALL_FILES_MERGED_13SEPT21.csv")

###-------------------------------------------------------
### COMBINE FINAL DATA FILES
# Set working directory
setwd("~/Desktop/Entrenchment_data_05MAY21/details_final")

#reads in all csv files in current directory
rm(list = ls()) #Clear environment
temp <- list.files(pattern='output_v2.txt')
for (i in 1:length(temp)) assign(temp[i], read.table(file=temp[i], header=TRUE, stringsAsFactors=FALSE))

#makes a list of all dataframes in the workspace!
vars <- ls()
nvars <- length(vars)
dataframes <-list()
j <- 1
for(i in 1:nvars)
{
  if(class(get(vars[i]))=="data.frame")
  {
    dataframes[[j]] <- get(vars[i])
    names(dataframes)[j] <- vars[i]
    j <- j+1
  }
}

length(dataframes[[1]]$timepoint)
names(dataframes)

for (i in 1:length(dataframes)){
  dataframes[[i]] <- cbind.data.frame(dataframes[[i]], mc_or_uni = rep("uni", length(dataframes[[i]]$timepoint)))
}

head(dataframes[[i]])

details_final_data <- bind_rows(dataframes, .id = "column_label")
head(details_final_data)
setwd("~/Desktop/Entrenchment_data_05MAY21/")
write.csv(details_final_data, "details_final_data_ALL_FILES_MERGED_13SEPT21.csv")


###-------------------------------------------------------
### READ IN TRANS AND FINAL DATASETS, MERGE
setwd("~/Desktop/Entrenchment_data_05MAY21/")
details_trans_data <- read.csv("details_trans_data_ALL_FILES_MERGED_13SEPT21.csv", header=T, stringsAsFactors=F)
details_final_data <- read.csv("details_final_data_ALL_FILES_MERGED_13SEPT21.csv", header=T, stringsAsFactors=F)
details_full_data <- rbind.data.frame(details_trans_data, details_final_data)
details_full_data <- cbind.data.frame(details_full_data, unique_mutation_id = c(paste(details_full_data$strategy_rep, "_",
                                                                                      details_full_data$timepoint, "_",
                                                                                      details_full_data$genome_location, "_",
                                                                                      details_full_data$mutation,
                                                                                      sep = "")))
head(details_full_data)
write.csv(details_full_data, "details_full_data_trans_and_final_MERGED_13SEPT21.csv")


###-------------------------------------------------------
### GET NUMBER OF UNIQUE REVERSION MUTATIONS
setwd("~/Desktop/Dropbox/Desktop_28NOV21/Entrenchment_processed_data_for_figures/SI_Figure_Reversion_mutation/Data")
details_full_data <- read.csv("details_full_data_trans_and_final_MERGED_13SEPT21.csv", header=T, stringsAsFactors=F)
details_full_data <- details_full_data[,-1]
head(details_full_data)
reversion_mutations <- unique(details_full_data$unique_mutation_id)
head(reversion_mutations)

strategy_rep <- c()
timepoint <- c()
genome_location <- c()
mutation <- c()

for (i in 1:length(reversion_mutations)){
  strategy_rep[length(strategy_rep)+1] <- strsplit(reversion_mutations[i], split="_")[[1]][1]
  timepoint[length(timepoint)+1] <- strsplit(reversion_mutations[i], split="_")[[1]][2]
  genome_location[length(genome_location)+1] <- strsplit(reversion_mutations[i], split="_")[[1]][3]
  mutation[length(mutation)+1] <- strsplit(reversion_mutations[i], split="_")[[1]][4]
}

reversion_mutation_df <- cbind.data.frame(strategy_rep, timepoint, genome_location, mutation)
head(reversion_mutation_df)
write.csv(reversion_mutation_df, "all_reversion_mutation_identities_14SEPT21.csv")


###-------------------------------------------------------

### EXTRACT AND SAVE RELEVANT SUMMARY DATA
df <- reversion_mutation_df
strategy_rep <- c()
timepoint <- c()
num_rev_m <- c()
num_32 <- c()
num_unique_m <- c()
num_unique_p <- c()
mut_div <- c()
pos_div <- c()

for (i in unique(df$strategy_rep)) {
  for (j in unique(df$timepoint)) {
    strategy_rep[length(strategy_rep)+1] <- i
    timepoint[length(timepoint)+1] <- j
    num_rev_m[length(num_rev_m)+1] <- length(df$strategy_rep[df$strategy_rep==i & df$timepoint==j])
    num_32[length(num_32)+1] <- length(df$strategy_rep[df$strategy_rep==i & df$timepoint==j & df$mutation==32])
    num_unique_m[length(num_unique_m)+1] <- length(unique(df$mutation[df$strategy_rep==i & df$timepoint==j]))
    num_unique_p[length(num_unique_p)+1] <- length(unique(df$genome_location[df$strategy_rep==i & df$timepoint==j]))
    mut_div[length(mut_div)+1] <- diversity(as.numeric(df$mutation[df$strategy_rep==i & df$timepoint==j]), index = "shannon")
    pos_div[length(pos_div)+1] <- diversity(as.numeric(df$genome_location[df$strategy_rep==i & df$timepoint==j]), index = "shannon")
  }
}

all_rev_m_summary_data <- cbind.data.frame(strategy_rep, timepoint, num_rev_m, num_32, num_unique_m, num_unique_p, mut_div, pos_div)
head(all_rev_m_summary_data)

write.csv(all_rev_m_summary_data, "all_reversion_mutations_summary_data.csv")


###-------------------------------------------------------
# Calculate the difference between the final and transition timepoint values 
# for all  relevant summary data from all_reversion_mutations dataframe
strategy_rep <- c()
num_rev_m_diff <- c()
num_32_diff <- c()
num_unique_m_diff <- c()
num_unique_p_diff <- c()
mut_div_diff <- c()
pos_div_diff <- c()

for (i in unique(all_rev_m_summary_data$strategy_rep)) {
  strategy_rep[length(strategy_rep)+1] <- i
  num_rev_m_diff[length(num_rev_m_diff)+1] <- as.numeric(all_rev_m_summary_data$num_rev_m[all_rev_m_summary_data$strategy_rep==i & all_rev_m_summary_data$timepoint==1]) - as.numeric(all_rev_m_summary_data$num_rev_m[all_rev_m_summary_data$strategy_rep==i & all_rev_m_summary_data$timepoint==0])
  num_32_diff[length(num_32_diff)+1] <- as.numeric(all_rev_m_summary_data$num_32[all_rev_m_summary_data$strategy_rep==i & all_rev_m_summary_data$timepoint==1]) - as.numeric(all_rev_m_summary_data$num_32[all_rev_m_summary_data$strategy_rep==i & all_rev_m_summary_data$timepoint==0])
  num_unique_m_diff[length(num_unique_m_diff)+1] <- as.numeric(all_rev_m_summary_data$num_unique_m[all_rev_m_summary_data$strategy_rep==i & all_rev_m_summary_data$timepoint==1]) - as.numeric(all_rev_m_summary_data$num_unique_m[all_rev_m_summary_data$strategy_rep==i & all_rev_m_summary_data$timepoint==0])
  num_unique_p_diff[length(num_unique_p_diff)+1] <- as.numeric(all_rev_m_summary_data$num_unique_p[all_rev_m_summary_data$strategy_rep==i & all_rev_m_summary_data$timepoint==1]) - as.numeric(all_rev_m_summary_data$num_unique_p[all_rev_m_summary_data$strategy_rep==i & all_rev_m_summary_data$timepoint==0])
  mut_div_diff[length(mut_div_diff)+1] <- as.numeric(all_rev_m_summary_data$mut_div[all_rev_m_summary_data$strategy_rep==i & all_rev_m_summary_data$timepoint==1]) - as.numeric(all_rev_m_summary_data$mut_div[all_rev_m_summary_data$strategy_rep==i & all_rev_m_summary_data$timepoint==0])
  pos_div_diff[length(pos_div_diff)+1] <- as.numeric(all_rev_m_summary_data$pos_div[all_rev_m_summary_data$strategy_rep==i & all_rev_m_summary_data$timepoint==1]) - as.numeric(all_rev_m_summary_data$pos_div[all_rev_m_summary_data$strategy_rep==i & all_rev_m_summary_data$timepoint==0])
}

diff_df <- cbind.data.frame(strategy_rep, num_rev_m_diff, timepoint=rep("difference", length(num_rev_m_diff)), num_32_diff, num_unique_m_diff, num_unique_p_diff, mut_div_diff, pos_div_diff)
head(diff_df)
write.csv(diff_df, "all_reversion_mutations_difference_data.csv")
