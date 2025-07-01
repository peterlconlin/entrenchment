# Compiling new lod fitness files from lod_fitness_from_detail_files
# 05 February 2024
# Peter Conlin

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot(font_size = 16))
library(dplyr)
library(directlabels)
library(reshape)
library(viridis)
library(ggrepel)
library(plot.matrix)
library(ComplexHeatmap)
library(stringr)
library(aphid)

setwd("/Users/peterconlin/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data/")

filenames <- list.files(pattern='*16AUG23.csv')

for (file in filenames){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("updated_rev_mut_IDs")){
    updated_rev_mut_IDs <- read.csv(file, header=TRUE, stringsAsFactors=F)
  }
  
  # if the merged dataset does exist, append to it
  else if (exists("updated_rev_mut_IDs")){
    temp_dataset <- read.csv(file, header=TRUE, stringsAsFactors=F)
    updated_rev_mut_IDs <- rbind.data.frame(updated_rev_mut_IDs, temp_dataset)
    rm(temp_dataset)
  }
}

updated_rev_mut_IDs <- updated_rev_mut_IDs[,-5]
colnames(updated_rev_mut_IDs)[3] <- "reversion_mut_ID"
head(updated_rev_mut_IDs)


duplicates_to_remove <- updated_rev_mut_IDs[updated_rev_mut_IDs$reversion_mut_ID != updated_rev_mut_IDs$new_rev_mut_ID,]
head(duplicates_to_remove)
length(duplicates_to_remove$strategy_rep)

setwd("/Users/peterconlin/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data/details_trans/lod_fitness_data_trans/")

filenames <- list.files(pattern='*.csv')
length(filenames)


for (file in filenames){
  
  # if the merged dataset doesn't exist, create it
  if (file == filenames[1]){
    lod_fitness_trans_from_detail <- read.csv(file, header=TRUE, stringsAsFactors=F)
  }
  
  # if the merged dataset does exist, append to it
  else if (file != filenames[1]){
    temp_dataset <- read.csv(file, header=TRUE, stringsAsFactors=F)
    lod_fitness_trans_from_detail <- rbind.data.frame(lod_fitness_trans_from_detail, temp_dataset)
    rm(temp_dataset)
  }
}

head(lod_fitness_trans_from_detail)
lod_fit_trans_detail <- left_join(lod_fitness_trans_from_detail, updated_rev_mut_IDs, by=c("strategy_rep", "timepoint", "reversion_mut_ID"))
head(lod_fit_trans_detail)

# lod_fit_trans_filtered <- lod_fit_trans_detail[!(lod_fit_trans_detail$duplicate==TRUE),]



setwd("/Users/peterconlin/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data/details_end/lod_fitness_data_final/")

filenames <- list.files(pattern='*.csv')

for (file in filenames){
  
  # if the merged dataset doesn't exist, create it
  if (file == filenames[1]){
    lod_fitness_final_from_detail <- read.csv(file, header=TRUE, stringsAsFactors=F)
  }
  
  # if the merged dataset does exist, append to it
  else if (file != filenames[1]){
    temp_dataset <- read.csv(file, header=TRUE, stringsAsFactors=F)
    lod_fitness_final_from_detail <- rbind.data.frame(lod_fitness_final_from_detail, temp_dataset)
    rm(temp_dataset)
  }
}

head(lod_fitness_final_from_detail)
lod_fit_final_detail <- left_join(lod_fitness_final_from_detail, updated_rev_mut_IDs, by=c("strategy_rep", "timepoint", "reversion_mut_ID"))



lod_fitness_detail <- rbind.data.frame(lod_fit_trans_detail, lod_fit_final_detail)
lod_fitness_detail_temp <- cbind.data.frame(lod_fitness_detail, index=paste(lod_fitness_detail$strategy_rep, 
                                                                            lod_fitness_detail$reversion_mut_ID, 
                                                                            lod_fitness_detail$timepoint, sep="_"))
head(lod_fitness_detail_temp)

duplicates_to_remove_temp <- cbind.data.frame(duplicates_to_remove, index=paste(duplicates_to_remove$strategy_rep,
                                                                                duplicates_to_remove$reversion_mut_ID,
                                                                                duplicates_to_remove$timepoint, sep="_"))
head(duplicates_to_remove_temp)

lod_fitness_detail_filtered <- lod_fitness_detail_temp[!(lod_fitness_detail_temp$index %in% duplicates_to_remove_temp$index),]


setwd("/Users/peterconlin/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/00_Data_compilation_and_QC")
write.csv(lod_fitness_detail_filtered, "lod_fitness_full_with_duplicates_filtered_24APR24.csv")



# duplicates <- lod_fit_trans_detail[is.na(lod_fit_trans_detail$new_rev_mut_ID)==FALSE & lod_fit_trans_detail$new_rev_mut_ID != lod_fit_trans_detail$reversion_mut_ID,]
# head(duplicates)
# duplicate_rows <- rownames(duplicates)
# head(duplicate_rows)
# 
# lod_fit_trans_detail$duplicate <- sapply(rownames(lod_fit_trans_detail), function(X) {ifelse(X %in% duplicate_rows, TRUE, FALSE)})
# lod_fit_trans_detail[lod_fit_trans_detail$duplicate==TRUE,]
# 
# ggplot(lod_fit_trans_detail[lod_fit_trans_detail$strategy_rep %in% c(4310) & lod_fit_trans_detail$mc_or_uni=="unicell",], aes(x=as.factor(strategy_rep), y=time_to_fill, color=duplicate)) +
#   geom_violin() +
#   scale_y_continuous(trans="log10") +
#   facet_wrap(~strategy_rep, ncol=1, scales="free")
# 
# ggplot(lod_fit_trans_detail[lod_fit_trans_detail$strategy_rep %in% c(4310) & lod_fit_trans_detail$mc_or_uni=="unicell",], aes(x=as.factor(strategy_rep), y=workload, color=duplicate)) +
#   geom_violin() +
#   scale_y_continuous(trans="log10") +
#   facet_wrap(~strategy_rep, ncol=1, scales="free")
# 
# 
# setwd("/Users/peterconlin/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data/details_end/lod_fitness_data_final/")
# 
# filenames <- list.files(pattern='*.csv')
# 
# for (file in filenames){
#   
#   # if the merged dataset doesn't exist, create it
#   if (file == filenames[1]){
#     lod_fitness_final_from_detail <- read.csv(file, header=TRUE, stringsAsFactors=F)
#   }
#   
#   # if the merged dataset does exist, append to it
#   else if (file != filenames[1]){
#     temp_dataset <- read.csv(file, header=TRUE, stringsAsFactors=F)
#     lod_fitness_final_from_detail <- rbind.data.frame(lod_fitness_final_from_detail, temp_dataset)
#     rm(temp_dataset)
#   }
# }
# 
# head(lod_fitness_final_from_detail)
# 
# lod_fit_final_detail <- left_join(lod_fitness_final_from_detail, updated_rev_mut_IDs, by=c("strategy_rep", "timepoint", "reversion_mut_ID"))
# 
# final_duplicates <- lod_fit_final_detail[is.na(lod_fit_final_detail$new_rev_mut_ID)==FALSE & lod_fit_final_detail$new_rev_mut_ID != lod_fit_final_detail$reversion_mut_ID,]
# head(final_duplicates)
# final_duplicate_rows <- rownames(final_duplicates)
# head(final_duplicate_rows)
# 
# lod_fit_final_detail$duplicate <- sapply(rownames(lod_fit_final_detail), function(X) {ifelse(X %in% final_duplicate_rows, TRUE, FALSE)})
# lod_fit_final_detail[lod_fit_final_detail$duplicate==TRUE,]
# head(lod_fit_final_detail)
# 
# 
# ggplot(lod_fit_final_detail[lod_fit_final_detail$strategy_rep %in% c(3038, 4278, 4310) & lod_fit_final_detail$mc_or_uni=="unicell",], aes(x=as.factor(strategy_rep), y=time_to_fill, color=duplicate)) +
#   geom_violin() +
#   scale_y_continuous(trans="log10") +
#   facet_wrap(~strategy_rep, ncol=3, scales="free")
# 
# ggplot(lod_fit_final_detail[lod_fit_final_detail$strategy_rep %in% c(3038, 4278, 4310) & lod_fit_final_detail$mc_or_uni=="unicell",], aes(x=as.factor(strategy_rep), y=workload, color=duplicate)) +
#   geom_violin() +
#   scale_y_continuous(trans="log10") +
#   facet_wrap(~strategy_rep, ncol=3, scales="free")
# 
# 
# lod_fit_full_detail <- rbind.data.frame(lod_fit_trans_detail, lod_fit_final_detail)
# lod_fit_full_detail_filtered <- lod_fit_full_detail[lod_fit_full_detail$duplicate==FALSE,]
# 
# setwd("/Users/peterconlin/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/00_Data_compilation_and_QC")
# write.csv(lod_fit_full_detail_filtered, "lod_fitness_detail_full_FILTERED_29MAR24.csv")

