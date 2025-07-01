# Identifying reversion mutations from Avidian genomes
# Peter Conlin
# 08 March 2023

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


### Functions
###--------------------------------
# 
# # Convert numbers to unique characters
# numbersToCharacters <- function(x){
#   unique_words <- c(0:37)
#   letter_x <- plyr::mapvalues(x,
#                               from = unique_words,
#                               to = c(letters[1:19], LETTERS[1:19]))
#   return(letter_x)
# }
# 
# # Find longest repeated sequence in the genome
# longestRepeatedSubstring <- function(str=str){
#   n <- length(str)
#   LCSRe <-c(rep(0, (n+1)*(n+1)))
#   LCSRe <- matrix(LCSRe, nrow=n+1, ncol=n+1)
#   res <- c()
#   res_length <- 0
#   index <- 1
#   
#   for (i in seq(1, n)) {
#     for (j in seq(i, n)) {
#       if (isTRUE(str[i] == str[j] & LCSRe[i,j] < (j-i))){
#         LCSRe[i+1,j+1] <- LCSRe[i,j]+1
#         if (LCSRe[i+1,j+1] > res_length){
#           res_length <- LCSRe[i+1,j+1]
#           index <- max(i,index)
#         }
#       }
#       else {
#         LCSRe[i,j] <- 0
#       }
#     }
#   }
#   if (res_length > 0){
#     for (i in seq(index - res_length + 1, index)){
#       res <- c(res, str[i])
#     }
#   }
#   return(res)
# }
# 
# 
### Set working directory
###--------------------------------

# Set working directory
# setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data/unis_trans")

#OLD DIRECTORIES....
#setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_identifying_reversion_mutations/SI_identifying_reversion_mutations_data/unis_trans")
#setwd("C:/Users/pconlin3/Dropbox/Desktop_28NOV21/Entrenchment_processed_data_for_figures/Figure_5/Figure_5_data/unis/All_timepoints") #For PC


### Get input data
###--------------------------------

# #reads in all files in the working directory and makes a data frame for each
# filenames <- list.files(pattern='*.csv')
# all_files <- lapply(filenames, function(i){read.csv(i, header=T, stringsAsFactors=F)})
# genomes <- bind_rows(all_files, .id = "column_label")
# colnames(genomes)[1] <- "organism_id"
# colnames(genomes)[2] <- "index"
# colnames(genomes)[14] <- "num_uni"
# colnames(genomes)[15] <- "iteration"
# #Reorder genomes dataframe by timepoint (num_uni) 
# genomes <- genomes[order(as.numeric(genomes$num_uni)),]
# head(genomes)

setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data/details_trans")
filenames <- list.files(pattern='*multicell_detail.dat')

head(uni_genomes[,1:10])
head(multi_genomes)

unique(multi_genomes$update)


for (file in filenames[1:length(filenames)]) {
  
  print(file)
  
  setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data/details_trans")
  genomes <- read.table(file, header=T)
  genomes <- genomes[,!names(genomes) %in% c("count")]
  genomes$strategy_rep <- rep(strsplit(file, split="_")[[1]][1], length(genomes$timepoint))
  genomes$mc_or_uni <- rep(strsplit(file, split="_")[[1]][2], length(genomes$timepoint))
  
  if (!"genome_location" %in% names(genomes)) {
    genomes$genome_location <- rep("NA", length(genomes$timepoint))
  }
  
  if (!"mutation" %in% names(genomes)) {
    genomes$mutation <- rep("NA", length(genomes$timepoint))
  }
  
  genomes <- cbind.data.frame(genomes, reversion_mut_ID=paste(genomes$genome_location, "_", genomes$mutation, sep=""))
  
  #print(head(genomes))
  
  # if (file == filenames[1]) {
  #   genomes <- temp
  # }
  # 
  # else {
  #   genomes <- rbind.data.frame(genomes, temp)
  # }
  
  strategy_rep <- c()
  reversion_mut_ID <- c()
  timepoint <- c()
  mc_or_uni <- c()
  iteration <- c()
  time_to_fill <- c()
  workload <- c()
  num_orgs <- c()
  total_cells <- c()
  
  
  for (i in unique(genomes$reversion_mut_ID)) {
    
    for (j in unique(genomes$iteration[genomes$reversion_mut_ID==i])) {
      
      temp <- genomes[genomes$reversion_mut_ID==i & genomes$iteration==j,]
      print(paste(i, "rep:", j, Sys.time(), sep=" "))
      
      strategy_rep[length(strategy_rep)+1] <- unique(temp$strategy_rep)
      reversion_mut_ID[length(reversion_mut_ID)+1] <- i
      timepoint[length(timepoint)+1] <- unique(temp$timepoint)
      mc_or_uni[length(mc_or_uni)+1] <- unique(temp$mc_or_uni)
      iteration[length(iteration)+1] <- j
      time_to_fill[length(time_to_fill)+1] <- unique(temp$update)
      workload[length(workload)+1] <- sum(temp$workload)
      num_orgs[length(num_orgs)+1] <- max(temp$multicell)+1
      total_cells[length(total_cells)+1] <- length(temp$timepoint)
      

    }
  }
  
  new_lod_fitness_df <- cbind.data.frame(strategy_rep, reversion_mut_ID, timepoint, mc_or_uni, iteration, time_to_fill, workload,  num_orgs, total_cells)
  #head(new_lod_fitness_df)
  
  setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data/details_trans/lod_fitness_data_trans")
  write.csv(new_lod_fitness_df, paste(file, "_lod_fitness_trans_16AUG23.csv", sep=""))
  
  rm(genomes)
  #rm(temp)
  rm(new_lod_fitness_df)
  
}







### DO THE SAME THING FOR THE FINAL TIMEPOINT

setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data/details_end")
filenames <- list.files(pattern='*multicell_detail.dat')


for (file in filenames[1:length(filenames)]) {
  
  print(file)
  
  setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data/details_end")
  genomes <- read.table(file, header=T)
  genomes <- genomes[,!names(genomes) %in% c("count")]
  genomes$strategy_rep <- rep(strsplit(file, split="_")[[1]][1], length(genomes$timepoint))
  genomes$mc_or_uni <- rep(strsplit(file, split="_")[[1]][2], length(genomes$timepoint))
  
  if (!"genome_location" %in% names(genomes)) {
    genomes$genome_location <- rep("NA", length(genomes$timepoint))
  }
  
  if (!"mutation" %in% names(genomes)) {
    genomes$mutation <- rep("NA", length(genomes$timepoint))
  }
  
  genomes <- cbind.data.frame(genomes, reversion_mut_ID=paste(genomes$genome_location, "_", genomes$mutation, sep=""))
  
  #print(head(genomes))
  
  # if (file == filenames[1]) {
  #   genomes <- temp
  # }
  # 
  # else {
  #   genomes <- rbind.data.frame(genomes, temp)
  # }
  
  strategy_rep <- c()
  reversion_mut_ID <- c()
  timepoint <- c()
  mc_or_uni <- c()
  iteration <- c()
  time_to_fill <- c()
  workload <- c()
  num_orgs <- c()
  total_cells <- c()
  
  
  for (i in unique(genomes$reversion_mut_ID)) {
    
    print(paste(i, Sys.time(), sep=" "))
    
    for (j in unique(genomes$iteration[genomes$reversion_mut_ID==i])) {
      
      temp <- genomes[genomes$reversion_mut_ID==i & genomes$iteration==j,]
      
      strategy_rep[length(strategy_rep)+1] <- unique(temp$strategy_rep)
      reversion_mut_ID[length(reversion_mut_ID)+1] <- i
      timepoint[length(timepoint)+1] <- unique(temp$timepoint)
      mc_or_uni[length(mc_or_uni)+1] <- unique(temp$mc_or_uni)
      iteration[length(iteration)+1] <- j
      time_to_fill[length(time_to_fill)+1] <- unique(temp$update)
      workload[length(workload)+1] <- sum(temp$workload)
      num_orgs[length(num_orgs)+1] <- max(temp$multicell)+1
      total_cells[length(total_cells)+1] <- length(temp$timepoint)
      
      
    }
  }
  
  new_lod_fitness_df <- cbind.data.frame(strategy_rep, reversion_mut_ID, timepoint, mc_or_uni, iteration, time_to_fill, workload,  num_orgs, total_cells)
  #head(new_lod_fitness_df)
  
  setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data/details_end/lod_fitness_data_final")
  write.csv(new_lod_fitness_df, paste(file, "_lod_fitness_final_16AUG23.csv", sep=""))
  
  rm(genomes)
  #rm(temp)
  rm(new_lod_fitness_df)
  
}





### DO THE SAME THING FOR THE SPLIT UP PIECES OF 4048 (TRANSITION)

setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data/details_trans/4048_unicell_detail_split")
filenames <- list.files(pattern='Part*')
filenames


for (file in filenames[1:length(filenames)]) {
  
  print(file)
  
  setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data/details_trans/4048_unicell_detail_split")
  genomes <- read.table(file, header=T)
  genomes <- genomes[,!names(genomes) %in% c("count")]
  genomes$strategy_rep <- rep(strsplit(file, split="_")[[1]][1], length(genomes$timepoint))
  genomes$mc_or_uni <- rep(strsplit(file, split="_")[[1]][2], length(genomes$timepoint))
  
  if (!"genome_location" %in% names(genomes)) {
    genomes$genome_location <- rep("NA", length(genomes$timepoint))
  }
  
  if (!"mutation" %in% names(genomes)) {
    genomes$mutation <- rep("NA", length(genomes$timepoint))
  }
  
  genomes <- cbind.data.frame(genomes, reversion_mut_ID=paste(genomes$genome_location, "_", genomes$mutation, sep=""))
  
  #print(head(genomes))
  
  # if (file == filenames[1]) {
  #   genomes <- temp
  # }
  # 
  # else {
  #   genomes <- rbind.data.frame(genomes, temp)
  # }
  
  strategy_rep <- c()
  reversion_mut_ID <- c()
  timepoint <- c()
  mc_or_uni <- c()
  iteration <- c()
  time_to_fill <- c()
  workload <- c()
  num_orgs <- c()
  total_cells <- c()
  
  
  for (i in unique(genomes$reversion_mut_ID)) {
    
    print(paste(i, Sys.time(), sep=" "))
    
    for (j in unique(genomes$iteration[genomes$reversion_mut_ID==i])) {
      
      temp <- genomes[genomes$reversion_mut_ID==i & genomes$iteration==j,]
      
      strategy_rep[length(strategy_rep)+1] <- unique(temp$strategy_rep)
      reversion_mut_ID[length(reversion_mut_ID)+1] <- i
      timepoint[length(timepoint)+1] <- unique(temp$timepoint)
      mc_or_uni[length(mc_or_uni)+1] <- unique(temp$mc_or_uni)
      iteration[length(iteration)+1] <- j
      time_to_fill[length(time_to_fill)+1] <- unique(temp$update)
      workload[length(workload)+1] <- sum(temp$workload)
      num_orgs[length(num_orgs)+1] <- max(temp$multicell)+1
      total_cells[length(total_cells)+1] <- length(temp$timepoint)
      
      
    }
  }
  
  new_lod_fitness_df <- cbind.data.frame(strategy_rep, reversion_mut_ID, timepoint, mc_or_uni, iteration, time_to_fill, workload,  num_orgs, total_cells)
  #head(new_lod_fitness_df)
  
  setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data/details_trans/4048_unicell_detail_split/lod_fitness_data")
  write.csv(new_lod_fitness_df, paste(file, "4048_unicell_detail_lod_fitness_trans_16AUG23.csv", sep=""))
  
  rm(genomes)
  #rm(temp)
  rm(new_lod_fitness_df)
  
}


setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data/details_trans/4048_unicell_detail_split/lod_fitness_data")
filenames <- list.files(pattern='Part*')
filenames

for (file in filenames){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.csv(file, header=TRUE, stringsAsFactors=F)
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <- read.csv(file, header=TRUE, stringsAsFactors=F)
    dataset <- rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
}

head(dataset)
dataset$strategy_rep <- rep(4048, length(dataset$strategy_rep))
dataset$mc_or_uni <- rep("uni", length(dataset$mc_or_uni))
setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data/details_trans/lod_fitness_data_trans/")
write.csv(dataset, "4048_unicell_detail.dat_lod_fitness_trans_16AUG23.csv")


test_4310_df <- read.table("4310_unicell_detail.dat", header=T)
head(test_4310_df[,1:8], 20)

setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data/details_end")
test_3038_df <- read.table("3038_unicell_detail.dat", header=T)
test_4278_df <- read.table("4278_unicell_detail.dat", header=T)
test_4310_df <- read.table("4310_unicell_detail.dat", header=T)


genomes <- cbind.data.frame(test_4310_df, reversion_mut_ID=paste(test_4310_df$genome_location, "_", test_4310_df$mutation, sep=""))
head(genomes)

unique(genomes$multicell)
unique(genomes$iteration)
sort(unique(genomes$organism_id))
sort(unique(genomes$organism_parent_id))
sort(unique(genomes$update))
sort(unique(genomes$reversion_mut_ID))

strategy_rep <- c()
reversion_mut_ID <- c()
timepoint <- c()
mc_or_uni <- c()
iteration <- c()
time_to_fill <- c()
workload <- c()
num_orgs <- c()
total_cells <- c()


for (i in unique(genomes$reversion_mut_ID)) {
  for (j in unique(genomes$iteration[genomes$reversion_mut_ID==i])) {
    temp <- genomes[genomes$reversion_mut_ID==i & genomes$iteration==j,]
    #strategy_rep[length(strategy_rep)+1] <- unique(temp$strategy_rep)
    reversion_mut_ID[length(reversion_mut_ID)+1] <- i
    iteration[length(iteration)+1] <- j
    total_cells[length(total_cells)+1] <- max(temp$organism_id)+1
    num_orgs[length(num_orgs)+1] <- max(temp$multicell)+1
    time_to_fill[length(time_to_fill)+1] <- unique(temp$update)
    workload[length(workload)+1] <- sum(temp$workload)
  }
}

new_lod_fitness_df <- cbind.data.frame(strategy_rep, reversion_mut_ID, timepoint, mc_or_uni, iteration, time_to_fill, workload,  num_orgs, num_cells)
new_lod_fitness_df


# # Create unique ID for each revertant by identifying unique first propagule genomes and assigning each a numerical ID
# unique_revertants <- unique(genomes$genome[genomes$organism_id==0])
# length(unique_revertants)
# uni_id <- rep(NA, length(genomes$genome))
# 
# 
# for (i in 1:length(unique_revertants)) {
#   uni_id[genomes$propagule_genome==unique_revertants[i]] <- i
# }
# 
# head(genomes, 36)
# 
# 
# # unique_revertants_t35 <- unique(genomes$propagule_genome[genomes$num_uni==35 & genomes$parent_id==0])
# # num_revertants_t35 <- length(unique_revertants_t35)
# # uni_id <- rep(NA, length(genomes$propagule_genome))
# # 
# # for (i in 1:length(unique_revertants_t35)) {
# #   uni_id[genomes$propagule_genome==unique_revertants_t35[i]] <- i
# # }
# # 
# # unique_revertants_t2600 <- unique(genomes$propagule_genome[genomes$num_uni==2600 & genomes$parent_id==0])
# # num_revertants_t2600 <- length(unique_revertants_t2600)
# # 
# # for (i in 1:length(unique_revertants_t2600)) {
# #   uni_id[genomes$propagule_genome==unique_revertants_t2600[i]] <- (i+num_revertants_t35)
# # }
# # 
# # unique_revertants_t6026 <- unique(genomes$propagule_genome[genomes$num_uni==6026 & genomes$parent_id==0])
# # num_revertants_t6026 <- length(unique_revertants_t6026)
# # 
# # for (i in 1:length(unique_revertants_t6026)) {
# #   uni_id[genomes$propagule_genome==unique_revertants_t6026[i]] <- (i+num_revertants_t35+num_revertants_t2600)
# # }
# # sort(unique(uni_id))
# # length(unique(uni_id)) == num_revertants_t35+num_revertants_t2600+num_revertants_t6026+1
# # 
# # 
# # length(genomes$organism_id[genomes$num_uni==35 & genomes$iteration==0])
# # length(genomes$organism_id[genomes$num_uni==2600 & genomes$iteration==0])
# 
# 
# # Fill in the gaps with the last non-NA value iterated over
# stored_value <- 0 
# 
# for (i in 1:length(uni_id)) {
#   
#   if (!is.na(uni_id[i])) {
#     
#     if (stored_value == uni_id[i]) {
#       next
#     }
#     
#     else if (stored_value != uni_id[i]) {
#       stored_value <- uni_id[i]
#     }
#   
#   }
#   
#   if (is.na(uni_id[i])) {
#     uni_id[i] <- stored_value
#   }
# 
# }
# 
# 
# genomes <- cbind.data.frame(genomes, uni_id)
# #Reorder genomes dataframe by timepoint (num_uni) and newly created uni_id
# genomes <- genomes[order(as.numeric(genomes$uni_id), as.numeric(genomes$iteration)),]
# head(genomes)


propagule_genomes <- genomes[genomes$organism_id==0,]
propagule_genomes <- cbind.data.frame(propagule_genomes, num_uni=rep(4310, length(propagule_genomes$timepoint)))
head(propagule_genomes)
unique(propagule_genomes$reversion_mut_ID)

# Loop through all propagule genomes to identify reversion mutations
###--------------------------------
setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data")
ancestor_genomes <- read.csv("multicelled_ancestor_genomes_at_transition_and_final_timepoint_29NOV21.csv", header=T, stringsAsFactors=F)

head(ancestor_genomes)

strategy_rep <- c()
timepoint <- c()
genome_length <- c()

for (i in unique(ancestor_genomes$strategy_rep)) {
  strategy_rep[length(strategy_rep)+1] <- i
  timepoint[length(timepoint)+1] <- 0
  genome_length[length(genome_length)+1]  <- length(ancestor_genomes$wt_mut[ancestor_genomes$strategy_rep==i & ancestor_genomes$timepoint==0])
}

for (i in unique(ancestor_genomes$strategy_rep)) {
  strategy_rep[length(strategy_rep)+1] <- i
  timepoint[length(timepoint)+1] <- 1
  genome_length[length(genome_length)+1]  <- length(ancestor_genomes$wt_mut[ancestor_genomes$strategy_rep==i & ancestor_genomes$timepoint==1])
}

genome_length_df <- cbind.data.frame(strategy_rep, timepoint, genome_length)

genome_length_df[genome_length_df$genome_length<100,]
genome_length_df[genome_length_df$genome_length>100,]

length(ancestor_genomes$wt_mut[ancestor_genomes$strategy_rep==4310 & ancestor_genomes$timepoint==0])


# Create empty vectors to populate
multi_parent_id <- c()
uni_id <- c()
iteration <- c()
anc_genome_length <- c()
revertant_genome_length <- c()
mut_position <- c()
wt_allele <- c()
mut_allele <- c()


for (h in unique(propagule_genomes$num_uni[propagule_genomes$num_uni==4310])) {
  
  anc_genome <- ancestor_genomes$wt_mut[ancestor_genomes$strategy_rep==h & ancestor_genomes$timepoint==1]
  length(anc_genome)
  print(Sys.time())
  print(paste("multicellular parental genotype ", h))
  
  for (i in unique(propagule_genomes$reversion_mut_ID[propagule_genomes$num_uni==h & propagule_genomes$timepoint==1])) {
    
    print(paste("unicellular revertant genotype ", i))
    
    for (j in unique(propagule_genomes$iteration[propagule_genomes$num_uni==h & propagule_genomes$timepoint==1 & propagule_genomes$reversion_mut_ID==i])) {
      
      print(paste("replicate ", j))
      org_genome <- as.numeric(unlist(strsplit(trimws(propagule_genomes$genome[propagule_genomes$num_uni==h & propagule_genomes$timepoint==1 & propagule_genomes$reversion_mut_ID==i & propagule_genomes$iteration==j], "b"), "\\s+")))

      glo <- align(list(anc_genome, list(org_genome)), type = "global", seqweights=c(0.5,0.5), progressive=T)
      #print(glo)
      
      temp <- c()
      for (p in 1:length(glo[1,])) {
        
        if (glo[1,p] == glo[2,p]) {
          temp[length(temp)+1] <- TRUE
        }
        
        # else if (glo[2,p] == "-" | glo[1,p] == "-") {
        #   temp[length(temp)+1] <- NA
        #   }
        
        else if (glo[1,p] != glo[2,p]) {
          multi_parent_id[length(multi_parent_id)+1] <- h
          #timepoint[length(timepoint)+1] <- genomes$num_uni[genomes$uni_id==h & genomes$iteration==i & genomes$index==j]
          uni_id[length(uni_id)+1] <- i
          iteration[length(iteration)+1] <- j
          #parent_org[length(parent_org)+1] <- genomes$parent_id[genomes$uni_id==h & genomes$iteration==i & genomes$index==j]
          #offspring_org[length(offspring_org)+1] <- genomes$offspring_id[genomes$uni_id==h & genomes$iteration==i & genomes$index==j]
          anc_genome_length[length(anc_genome_length)+1] <- length(anc_genome)
          revertant_genome_length[length(revertant_genome_length)+1] <- length(org_genome)
          mut_position[length(mut_position)+1] <- length(glo[1,1:p][glo[1,1:p]!="-"])
          wt_allele[length(wt_allele)+1] <- glo[1,p]
          mut_allele[length(mut_allele)+1] <- glo[2,p]
        }
      }
      
    }
    
  }
  
}


reversion_mutations <- cbind.data.frame(multi_parent_id, uni_id, iteration, anc_genome_length, 
                                        revertant_genome_length, mut_position, wt_allele, mut_allele, 
                                        reversion_mut_ID=paste((mut_position-1), mut_allele, sep="_"))
head(reversion_mutations)

reversion_mutations[reversion_mutations$wt_allele!="-" & reversion_mutations$mut_allele!="-",]

#reversion_mutations[reversion_mutations$wt_allele!="-",]

old_rev_mut_ID <- c()
new_rev_mut_ID <- c()
num_revertants <- c()

for (i in unique(reversion_mutations$uni_id)) {
  my_tab <- as.data.frame(table(reversion_mutations$reversion_mut_ID[reversion_mutations$uni_id == i & 
                                                         reversion_mutations$wt_allele!="-" & 
                                                         reversion_mutations$mut_allele!="-"]))
  my_tab_sort2 <- my_tab[order(my_tab$Freq, decreasing = TRUE),]
  rownames(my_tab_sort2) <- NULL
  old_rev_mut_ID[length(old_rev_mut_ID)+1] <- i
  new_rev_mut_ID[length(new_rev_mut_ID)+1] <- paste(my_tab_sort2$Var1[1])
  num_revertants[length(num_revertants)+1] <- my_tab_sort2$Freq[1]
}



updates_to_rev_mut <- cbind.data.frame(old_rev_mut_ID, new_rev_mut_ID, num_revertants)
updates_to_rev_mut


write.csv(updates_to_rev_mut, "4310F_updated_reversion_mutation_identities_16AUG23.csv")









#### Final timepoint
# Set working directory
setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data/unis_end")

#OLD DIRECTORIES....
#setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_identifying_reversion_mutations/SI_identifying_reversion_mutations_data/unis_trans")
#setwd("C:/Users/pconlin3/Dropbox/Desktop_28NOV21/Entrenchment_processed_data_for_figures/Figure_5/Figure_5_data/unis/All_timepoints") #For PC


### Get input data
###--------------------------------

# #reads in all files in the working directory and makes a data frame for each
# filenames <- list.files(pattern='*unicell_detail.dat')
# all_files <- lapply(filenames, function(i){read.csv(i, header=T, stringsAsFactors=F)})
# final_genomes <- bind_rows(all_files, .id = "column_label")
# colnames(final_genomes)[1] <- "organism_id"
# colnames(final_genomes)[2] <- "index"
# colnames(final_genomes)[14] <- "num_uni"
# colnames(final_genomes)[15] <- "iteration"
# #Reorder genomes dataframe by timepoint (num_uni) 
# final_genomes <- final_genomes[order(as.numeric(final_genomes$num_uni)),]
# head(final_genomes)

setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data/details_end")
final_4310_df <- read.table("4310_unicell_detail.dat", header=T)
head(final_4310_df)
head(final_4310_df[,1:8], 20)

unique(final_4310_df$genome_location)
unique(sort(final_4310_df$mutation))
unique(final_4310_df$iteration)

genomes <- cbind.data.frame(final_4310_df, reversion_mut_ID=paste(final_4310_df$genome_location, "_", final_4310_df$mutation, sep=""))
head(genomes)

# # Create unique ID for each revertant by identifying unique first propagule genomes and assigning each a numerical ID
# unique_revertants <- unique(genomes$genome[genomes$organism_id==0])
# length(unique_revertants)
# uni_id <- rep(NA, length(genomes$genome))
# 
# 
# for (i in 1:length(unique_revertants)) {
#   uni_id[genomes$propagule_genome==unique_revertants[i]] <- i
# }
# 
# head(genomes, 36)
# 
# 
# # unique_revertants_t35 <- unique(genomes$propagule_genome[genomes$num_uni==35 & genomes$parent_id==0])
# # num_revertants_t35 <- length(unique_revertants_t35)
# # uni_id <- rep(NA, length(genomes$propagule_genome))
# # 
# # for (i in 1:length(unique_revertants_t35)) {
# #   uni_id[genomes$propagule_genome==unique_revertants_t35[i]] <- i
# # }
# # 
# # unique_revertants_t2600 <- unique(genomes$propagule_genome[genomes$num_uni==2600 & genomes$parent_id==0])
# # num_revertants_t2600 <- length(unique_revertants_t2600)
# # 
# # for (i in 1:length(unique_revertants_t2600)) {
# #   uni_id[genomes$propagule_genome==unique_revertants_t2600[i]] <- (i+num_revertants_t35)
# # }
# # 
# # unique_revertants_t6026 <- unique(genomes$propagule_genome[genomes$num_uni==6026 & genomes$parent_id==0])
# # num_revertants_t6026 <- length(unique_revertants_t6026)
# # 
# # for (i in 1:length(unique_revertants_t6026)) {
# #   uni_id[genomes$propagule_genome==unique_revertants_t6026[i]] <- (i+num_revertants_t35+num_revertants_t2600)
# # }
# # sort(unique(uni_id))
# # length(unique(uni_id)) == num_revertants_t35+num_revertants_t2600+num_revertants_t6026+1
# # 
# # 
# # length(genomes$organism_id[genomes$num_uni==35 & genomes$iteration==0])
# # length(genomes$organism_id[genomes$num_uni==2600 & genomes$iteration==0])
# 
# 
# # Fill in the gaps with the last non-NA value iterated over
# stored_value <- 0 
# 
# for (i in 1:length(uni_id)) {
#   
#   if (!is.na(uni_id[i])) {
#     
#     if (stored_value == uni_id[i]) {
#       next
#     }
#     
#     else if (stored_value != uni_id[i]) {
#       stored_value <- uni_id[i]
#     }
#   
#   }
#   
#   if (is.na(uni_id[i])) {
#     uni_id[i] <- stored_value
#   }
# 
# }
# 
# 
# genomes <- cbind.data.frame(genomes, uni_id)
# #Reorder genomes dataframe by timepoint (num_uni) and newly created uni_id
# genomes <- genomes[order(as.numeric(genomes$uni_id), as.numeric(genomes$iteration)),]
# head(genomes)


propagule_genomes <- genomes[genomes$organism_id==0,]
propagule_genomes <- cbind.data.frame(propagule_genomes, num_uni=rep(4310, length(propagule_genomes$timepoint)))
head(propagule_genomes)


# Loop through all propagule genomes to identify reversion mutations
###--------------------------------
setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Identifying_reversion_mutations/SI_identifying_reversion_mutations_data")
ancestor_genomes <- read.csv("multicelled_ancestor_genomes_at_transition_and_final_timepoint_29NOV21.csv", header=T, stringsAsFactors=F)

for (i in unique(ancestor_genomes$strategy_rep)) {
  print(paste(i, length(ancestor_genomes$wt_mut[ancestor_genomes$strategy_rep==i & ancestor_genomes$timepoint==0])))
}

for (i in unique(ancestor_genomes$strategy_rep)) {
  print(paste(i, length(ancestor_genomes$wt_mut[ancestor_genomes$strategy_rep==i & ancestor_genomes$timepoint==1])))
}

length(ancestor_genomes$wt_mut[ancestor_genomes$strategy_rep==4310 & ancestor_genomes$timepoint==1])

# Create empty vectors to populate
multi_parent_id <- c()
uni_id <- c()
iteration <- c()
anc_genome_length <- c()
revertant_genome_length <- c()
mut_position <- c()
wt_allele <- c()
mut_allele <- c()


for (h in unique(propagule_genomes$num_uni[propagule_genomes$num_uni==4310])) {
  
  anc_genome <- ancestor_genomes$wt_mut[ancestor_genomes$strategy_rep==h & ancestor_genomes$timepoint==1]
  length(anc_genome)
  print(Sys.time())
  print(paste("multicellular parental genotype ", h))
  
  for (i in unique(propagule_genomes$reversion_mut_ID[propagule_genomes$num_uni==h])) {
    
    #print(paste("unicellular revertant genotype ", i))
    
    for (j in unique(propagule_genomes$iteration[propagule_genomes$num_uni==h & propagule_genomes$reversion_mut_ID==i])) {
      
      print(paste("replicate ", j))
      org_genome <- as.numeric(unlist(strsplit(trimws(propagule_genomes$genome[propagule_genomes$num_uni==h & propagule_genomes$reversion_mut_ID==i & propagule_genomes$iteration==j], "b"), "\\s+")))
      
      glo <- align(list(anc_genome, list(org_genome)), type = "global", seqweights=c(0.5,0.5), progressive=T)
      print(glo)
      
      temp <- c()
      for (p in 1:length(glo[1,])) {
        
        if (glo[1,p] == glo[2,p]) {
          temp[length(temp)+1] <- TRUE
        }
        
        # else if (glo[2,p] == "-" | glo[1,p] == "-") {
        #   temp[length(temp)+1] <- NA
        #   }
        
        else if (glo[1,p] != glo[2,p]) {
          multi_parent_id[length(multi_parent_id)+1] <- h
          #timepoint[length(timepoint)+1] <- genomes$num_uni[genomes$uni_id==h & genomes$iteration==i & genomes$index==j]
          uni_id[length(uni_id)+1] <- i
          iteration[length(iteration)+1] <- j
          #parent_org[length(parent_org)+1] <- genomes$parent_id[genomes$uni_id==h & genomes$iteration==i & genomes$index==j]
          #offspring_org[length(offspring_org)+1] <- genomes$offspring_id[genomes$uni_id==h & genomes$iteration==i & genomes$index==j]
          anc_genome_length[length(anc_genome_length)+1] <- length(anc_genome)
          revertant_genome_length[length(revertant_genome_length)+1] <- length(org_genome)
          mut_position[length(mut_position)+1] <- length(glo[1,1:p][glo[1,1:p]!="-"])
          wt_allele[length(wt_allele)+1] <- glo[1,p]
          mut_allele[length(mut_allele)+1] <- glo[2,p]
        }
      }
      
    }
    
  }
  
}


reversion_mutations <- cbind.data.frame(multi_parent_id, uni_id, iteration, anc_genome_length, 
                                        revertant_genome_length, mut_position, wt_allele, mut_allele, 
                                        reversion_mut_ID=paste((mut_position-1),mut_allele, sep="_"))

reversion_mutations[reversion_mutations$wt_allele!="-" & reversion_mutations$mut_allele!="-",]

#reversion_mutations[reversion_mutations$wt_allele!="-",]
my_tab <- table(reversion_mutations$reversion_mut_ID[reversion_mutations$wt_allele!="-" & reversion_mutations$mut_allele!="-"])

my_tab_sort2 <- my_tab[order(my_tab, decreasing = TRUE)]
my_tab_sort2                               # Print ordered table

unique(propagule_genomes$reversion_mut_ID)

length(unique(reversion_mutations$reversion_mut_ID[reversion_mutations$multi_parent_id==4310]))













































for (h in unique(propagule_genomes$uni_id)) {
  
  anc_genome <- as.numeric(unlist(strsplit(genomes$propagule_genome[genomes$uni_id==h & genomes$index==min(genomes$index[genomes$uni_id==h])], "\\s+")))
  anc_genome <- anc_genome[!is.na(anc_genome)]
  print(Sys.time())
  print(paste("genotype ", h))
  
  for (i in unique(genomes$iteration[genomes$uni_id==h])) {
    
    print(paste("replicate ",i))
    
    if (max(genomes$org_size[genomes$uni_id==h & genomes$iteration==i]) == 1) {
      
      for (j in unique(genomes$index[genomes$uni_id==h & genomes$iteration==i])) {
        
        org_genome <- as.numeric(unlist(strsplit(genomes$offspring_genome[genomes$uni_id==h & genomes$iteration==i & genomes$index==j], "\\s+")))
        org_genome <- org_genome[!is.na(org_genome)]
        
        parent_org_genome <- as.numeric(unlist(strsplit(genomes$propagule_genome[genomes$uni_id==h & genomes$iteration==i & genomes$index==j], "\\s+")))
        parent_org_genome <- parent_org_genome[!is.na(parent_org_genome)]
        
        offspring_org_genome <- as.numeric(unlist(strsplit(genomes$offspring_genome[genomes$uni_id==h & genomes$iteration==i & genomes$index==j], "\\s+")))
        offspring_org_genome <- offspring_org_genome[!is.na(offspring_org_genome)]
        
        if (length(parent_org_genome) != 0) {
          
          glo <- align(list(parent_org_genome, list(org_genome)), type = "global", seqweights=c(1,0), progressive=T)
          
          temp <- c()
          for (p in 1:length(glo[1,])) {
            
            if (glo[1,p] == glo[2,p]) {
              temp[length(temp)+1] <- TRUE
            }
            
            # else if (glo[2,p] == "-" | glo[1,p] == "-") {
            #   temp[length(temp)+1] <- NA
            #   }
            
            else if (glo[1,p] != glo[2,p]) {
              uni_id[length(uni_id)+1] <- h
              timepoint[length(timepoint)+1] <- genomes$num_uni[genomes$uni_id==h & genomes$iteration==i & genomes$index==j]
              iteration[length(iteration)+1] <- i
              index[length(index)+1] <- j
              parent_org[length(parent_org)+1] <- genomes$parent_id[genomes$uni_id==h & genomes$iteration==i & genomes$index==j]
              offspring_org[length(offspring_org)+1] <- genomes$offspring_id[genomes$uni_id==h & genomes$iteration==i & genomes$index==j]
              parent_genome_length[length(parent_genome_length)+1] <- length(parent_org_genome)
              offspring_genome_length[length(offspring_genome_length)+1] <- length(offspring_org_genome)
              mut_position[length(mut_position)+1] <- length(glo[1,1:p][glo[1,1:p]!="-"])
              wt_allele[length(wt_allele)+1] <- glo[1,p]
              mut_allele[length(mut_allele)+1] <- glo[2,p]
            }
          }
          # if (length(temp[temp==FALSE])>25) {
          #   print(h)
          # }
          #initial_genome_diffs <- rbind(initial_genome_diffs, temp)
          
          #if (length(anc_genome) == length(org_genome)) {
          
          #if (length(anc_genome) == length(parent_org_genome)) {
          
          # for (p in 1:length(org_genome)) {
          #   
          #   if (is.na(org_genome[p]) == FALSE) {
          #     
          #     if (org_genome[p] != parent_org_genome[p]) {
          #       
          #       uni_id[length(uni_id)+1] <- h
          #       timepoint[length(timepoint)+1] <- genomes$timepoint[genomes$uni_id==h & genomes$iteration==i & genomes$index==j]
          #       iteration[length(iteration)+1] <- i
          #       replication_event[length(replication_event)+1] <- j
          #       parent_org[length(parent_org)+1] <- genomes$parent_id[genomes$uni_id==h & genomes$iteration==i & genomes$index==j]
          #       mut_position[length(mut_position)+1] <- p
          #       wt_allele[length(wt_allele)+1] <- parent_org_genome[p]
          #       mut_allele[length(mut_allele)+1] <- org_genome[p]
          #       
          #     }
          #   }
          # }
        }
        
        #   else {
        #     print("Parent has indel!")
        #     has_indel_iteration[length(has_indel_iteration)+1] <- i
        #     has_indel[length(has_indel)+1] <- j
        #   }
        #   
        # }
        
        #   else {
        #     print("Offspring has indel!")
        #     has_indel_genotype[length(has_indel_genotype)+1] <- h
        #     has_indel_iteration[length(has_indel_iteration)+1] <- i
        #     has_indel[length(has_indel)+1] <- j
        #   }
        #   
        # }
        
        else {
          print("Missing parent organism!")
          has_missing_parent_genotype[length(has_missing_parent)+1] <- h
          has_missing_parent_iteration[length(has_missing_parent_iteration)+1] <- i
          has_missing_parent[length(has_missing_parent)+1] <- j
        }
        
      }
    }
    
    else if (max(genomes$org_size[genomes$uni_id==h & genomes$iteration==i]) > 1){
      print("Multicellularity re-evolved!")
      has_multicell_genotype[length(has_multicell_genotype)+1] <- h
      has_multicell_contaminant[length(has_multicell_contaminant)+1] <- i
    }
  }
}


























#---------------------------------------------------------
# NEED TO CONFIRM THIS LATER!
# Check that each uni_id is used exclusively within a given timepoint 
# (Want to see no identical genomes at different timepoints, recycled uni_ids)
# NEED TO CONFIRM THIS LATER! 



#---------------------------------------------------------
# GET LIST OF REVERSION MUTATIONS
#make empty vectors to save data
revertant_genome <- c()
strategy_rep <- c()

revertant_genomes <- c(unique(genomes$propagule_genome[genomes$parent_id==0]))
uni_id <- c(rep(unique(genomes$uni_id[temp$parent_id==0]), 
                                      length(unique(genomes$propagule_genome[temp$parent_id==0]))))


length(propagule_genome)
length(strategy_rep)













#---------------------------------------------------------


# CALCUATE THE TIME TO FILL FOR EACH POPULATION
head(genomes)

uni_id <- c()
iteration <- c()
timepoint <- c()
time_to_fill <- c()

for (i in unique(genomes$num_uni)) {
  for (j in unique(genomes$uni_id[genomes$num_uni==i])) {
    for (k in unique(genomes$iteration[genomes$num_uni==i & genomes$uni_id==j])) {
      timepoint[length(timepoint)+1] <- i
      uni_id[length(uni_id)+1] <- j
      iteration[length(iteration)+1] <- k
      time_to_fill[length(time_to_fill)+1] <- max(genomes$update[genomes$num_uni==i & genomes$uni_id==j & genomes$iteration==k])
    }
  }
}

time_to_fill_df <- cbind.data.frame(timepoint, uni_id, iteration, time_to_fill)
head(time_to_fill_df)

#setwd("~/Desktop/Dropbox/Desktop_28NOV21/Entrenchment_processed_data_for_figures/Figure_5/Figure_5_data/")
write.csv(time_to_fill_df, "time_to_fill_data_for_3416_case_study_all_timepoints_08DEC21.csv")


# CALCUATING THE NUMBER OF MUTATIONS (m)

# Create empty vectors to populate
uni_id <- c()
timepoint <- c()
iteration <- c()
index <- c()
parent_org <- c()
offspring_org <- c()
parent_genome_length <- c()
offspring_genome_length <- c()
mut_position <- c()
wt_allele <- c()
mut_allele <- c()
has_indel <- c()
has_indel_iteration <- c()
has_indel_genotype <- c()
has_missing_parent_genotype <- c()
has_missing_parent_iteration <- c()
has_missing_parent <- c()
has_multicell_genotype <- c()
has_multicell_contaminant <- c()


# Loop through all 
###--------------------------------

for (h in unique(genomes$uni_id[genomes$uni_id > 463])) {
  
  anc_genome <- as.numeric(unlist(strsplit(genomes$propagule_genome[genomes$uni_id==h & genomes$index==min(genomes$index[genomes$uni_id==h])], "\\s+")))
  anc_genome <- anc_genome[!is.na(anc_genome)]
  print(Sys.time())
  print(paste("genotype ", h))
  
  for (i in unique(genomes$iteration[genomes$uni_id==h])) {
    
    print(paste("replicate ",i))
    
    if (max(genomes$org_size[genomes$uni_id==h & genomes$iteration==i]) == 1) {
      
      for (j in unique(genomes$index[genomes$uni_id==h & genomes$iteration==i])) {
        
        org_genome <- as.numeric(unlist(strsplit(genomes$offspring_genome[genomes$uni_id==h & genomes$iteration==i & genomes$index==j], "\\s+")))
        org_genome <- org_genome[!is.na(org_genome)]
        
        parent_org_genome <- as.numeric(unlist(strsplit(genomes$propagule_genome[genomes$uni_id==h & genomes$iteration==i & genomes$index==j], "\\s+")))
        parent_org_genome <- parent_org_genome[!is.na(parent_org_genome)]
        
        offspring_org_genome <- as.numeric(unlist(strsplit(genomes$offspring_genome[genomes$uni_id==h & genomes$iteration==i & genomes$index==j], "\\s+")))
        offspring_org_genome <- offspring_org_genome[!is.na(offspring_org_genome)]
        
        if (length(parent_org_genome) != 0) {
          
          glo <- align(list(parent_org_genome, list(org_genome)), type = "global", seqweights=c(1,0), progressive=T)
          
          temp <- c()
          for (p in 1:length(glo[1,])) {
            
            if (glo[1,p] == glo[2,p]) {
              temp[length(temp)+1] <- TRUE
            }
            
            # else if (glo[2,p] == "-" | glo[1,p] == "-") {
            #   temp[length(temp)+1] <- NA
            #   }
            
            else if (glo[1,p] != glo[2,p]) {
              uni_id[length(uni_id)+1] <- h
              timepoint[length(timepoint)+1] <- genomes$num_uni[genomes$uni_id==h & genomes$iteration==i & genomes$index==j]
              iteration[length(iteration)+1] <- i
              index[length(index)+1] <- j
              parent_org[length(parent_org)+1] <- genomes$parent_id[genomes$uni_id==h & genomes$iteration==i & genomes$index==j]
              offspring_org[length(offspring_org)+1] <- genomes$offspring_id[genomes$uni_id==h & genomes$iteration==i & genomes$index==j]
              parent_genome_length[length(parent_genome_length)+1] <- length(parent_org_genome)
              offspring_genome_length[length(offspring_genome_length)+1] <- length(offspring_org_genome)
              mut_position[length(mut_position)+1] <- length(glo[1,1:p][glo[1,1:p]!="-"])
              wt_allele[length(wt_allele)+1] <- glo[1,p]
              mut_allele[length(mut_allele)+1] <- glo[2,p]
            }
          }
          # if (length(temp[temp==FALSE])>25) {
          #   print(h)
          # }
          #initial_genome_diffs <- rbind(initial_genome_diffs, temp)
          
          #if (length(anc_genome) == length(org_genome)) {
          
          #if (length(anc_genome) == length(parent_org_genome)) {
          
          # for (p in 1:length(org_genome)) {
          #   
          #   if (is.na(org_genome[p]) == FALSE) {
          #     
          #     if (org_genome[p] != parent_org_genome[p]) {
          #       
          #       uni_id[length(uni_id)+1] <- h
          #       timepoint[length(timepoint)+1] <- genomes$timepoint[genomes$uni_id==h & genomes$iteration==i & genomes$index==j]
          #       iteration[length(iteration)+1] <- i
          #       replication_event[length(replication_event)+1] <- j
          #       parent_org[length(parent_org)+1] <- genomes$parent_id[genomes$uni_id==h & genomes$iteration==i & genomes$index==j]
          #       mut_position[length(mut_position)+1] <- p
          #       wt_allele[length(wt_allele)+1] <- parent_org_genome[p]
          #       mut_allele[length(mut_allele)+1] <- org_genome[p]
          #       
          #     }
          #   }
          # }
        }
        
        #   else {
        #     print("Parent has indel!")
        #     has_indel_iteration[length(has_indel_iteration)+1] <- i
        #     has_indel[length(has_indel)+1] <- j
        #   }
        #   
        # }
        
        #   else {
        #     print("Offspring has indel!")
        #     has_indel_genotype[length(has_indel_genotype)+1] <- h
        #     has_indel_iteration[length(has_indel_iteration)+1] <- i
        #     has_indel[length(has_indel)+1] <- j
        #   }
        #   
        # }
        
        else {
          print("Missing parent organism!")
          has_missing_parent_genotype[length(has_missing_parent)+1] <- h
          has_missing_parent_iteration[length(has_missing_parent_iteration)+1] <- i
          has_missing_parent[length(has_missing_parent)+1] <- j
        }
        
      }
    }
    
    else if (max(genomes$org_size[genomes$uni_id==h & genomes$iteration==i]) > 1){
      print("Multicellularity re-evolved!")
      has_multicell_genotype[length(has_multicell_genotype)+1] <- h
      has_multicell_contaminant[length(has_multicell_contaminant)+1] <- i
    }
  }
}


###--------------------------------

mutation_df_464_and_over <- cbind.data.frame(uni_id, timepoint, iteration, index, parent_org, parent_genome_length, offspring_org, offspring_genome_length, mut_position, wt_allele, mut_allele)
head(mutation_df)
write.csv(mutation_df, "raw_mutation_analysis_output_3416_unicell_revertants_464_and_over_filtered_14DEC21.csv")
`

# indel_df <- cbind.data.frame(has_indel_iteration, has_indel)
# head(indel_df)
# write.csv(indel_df, "indel_mutation_analysis_output_3416_over_time_T35_filtered_21FEB21.csv")

# 
# has_missing_parent_df <- cbind.data.frame(has_missing_parent_genotype, has_missing_parent_iteration, has_missing_parent)
# head(has_missing_parent_df)
# write.csv(has_missing_parent_df, "has_missing_parent_raw_mutation_analysis_3416_T35_filtered_02FEB21.csv")

has_multicell_df <- cbind.data.frame(has_multicell_genotype, has_multicell_contaminant)
write.csv(has_multicell_df, "has_multicell_raw_mutation_analysis_3416_all_timepoints_filtered_08DEC21.csv")




multicell_genomes <- has_multicell_df
colnames(multicell_genomes) <- c("uni_id", "iteration")
filtered_genomes <- anti_join(genomes, multicell_genomes, by=c("uni_id", "iteration"))


# Retain rows with unique 
mutation_df_filtered <- mutation_df %>% distinct(uni_id, iteration, parent_org, mut_position, wt_allele, mut_allele, .keep_all = TRUE)
head(mutation_df_filtered)



# Calculate the number of unique mutations per culture (iteration)
uni_id <- c()
iteration <- c()
m <- c()

for (i in unique(mutation_df_filtered$uni_id)) {
  
  for (j in unique(mutation_df_filtered$iteration[mutation_df_filtered$uni_id==i])) {
    
    #generation[length(generation)+1] <- i
    uni_id[length(uni_id)+1] <- i
    iteration[length(iteration)+1] <- j
    m[length(m)+1] <- length(mutation_df_filtered$mut_allele[mutation_df_filtered$uni_id==i & 
                                                               mutation_df_filtered$iteration==j & 
                                                               mutation_df_filtered$mut_allele!="-" & 
                                                               mutation_df_filtered$wt_allele!="-"])
  }
}

mutation_number_summary <- cbind.data.frame(uni_id, iteration, m)
head(mutation_number_summary)



# Make dataframe with a summary of replication data for each iteration and num_uni combination
timepoint <- c()
uni_id <- c()
iteration <- c()
counts <- c()
avg_workload <- c()


for (i in unique(filtered_genomes$uni_id)) {
  for (j in unique(filtered_genomes$iteration[filtered_genomes$uni_id==i])) {
    timepoint[length(timepoint)+1] <- unique(filtered_genomes$num_uni[filtered_genomes$uni_id==i & filtered_genomes$iteration==j])
    uni_id[length(uni_id)+1] <- i
    iteration[length(iteration)+1] <- j
    counts[length(counts)+1] <- length(filtered_genomes$organism_id[filtered_genomes$uni_id==i & filtered_genomes$iteration==j])
    
    temp <- filtered_genomes[filtered_genomes$uni_id==i & filtered_genomes$iteration==j,]
    max_workload <- c()
    
    for (k in unique(temp$parent_id)) {
      max_workload[length(max_workload)+1] <- max(temp$workload[temp$parent_id==k])
    }
    
    avg_workload[length(avg_workload)+1] <- sum(max_workload)/length(filtered_genomes$organism_id[filtered_genomes$uni_id==i & filtered_genomes$iteration==j])
    
  }
}

replication_counts <-cbind.data.frame(timepoint, uni_id, iteration, counts, avg_workload)
#write.csv(replication_counts, "replication_summary_data_3416_unicell_revertants_T35_T2600_T6026_filtered_23NOV21.csv")





# num_uni <- c()
# iteration <- c()
# total_replications <- c()
# 
# for (i in unique(genomes$num_uni)) {
#   for (j in unique(genomes$iteration[genomes$num_uni==i])) {
#     num_uni[length(num_uni)+1] <- i
#     iteration[length(iteration)+1] <- j
#     total_replications[length(total_replications)+1] <- length(genomes$organism_id[genomes$num_uni==i & genomes$iteration==j])
#   }
# }
# 
# replication_summary_df <- cbind.data.frame(num_uni, iteration, total_replications)
# head(replication_summary_df)
# write.csv(replication_summary_df, "replication_summary_data_3416_unicell_revertants_T35_T2600_T6026_filtered_09NOV21.csv")



# Combine mutation_number_summary w/ replication_summary_df
mutation_summary_df <- full_join(mutation_number_summary, replication_counts, by=c("uni_id", "iteration"))
mutation_summary_complete <- cbind.data.frame(mutation_summary_df, mutation_rate=mutation_summary_df$m/mutation_summary_df$counts)
mutation_summary_complete <- left_join(mutation_summary_complete, time_to_fill_df, by=c("uni_id", "iteration", "timepoint"))
mutation_summary_complete <- cbind.data.frame(mutation_summary_complete, growth_rate = log(mutation_summary_complete$counts/1)/mutation_summary_complete$time_to_fill)
head(mutation_summary_complete)


# Calculate genome size for multicell genome at each timepoint
timepoint <- c()
genome_size <- c()

for (i in unique(genomes$num_uni)) {
  anc_genome <- c(strsplit(unname(unlist(genomes$propagule_genome[genomes$num_uni==i & genomes$iteration==0 & genomes$parent_id==0]))[[1]], split="   "))
  timepoint[length(timepoint)+1] <- i
  genome_size[length(genome_size)+1] <- length(anc_genome[[1]][!is.na(anc_genome)])
}

genome_size_df <- cbind.data.frame(timepoint, genome_size)


# Merge genome_size_df w/ avg_mutation_rate_df
mutation_summary_complete <- left_join(mutation_summary_complete, genome_size_df, by="timepoint")
head(mutation_summary_complete)
# Calculate a per site mutation rate
mutation_summary_complete <- cbind.data.frame(mutation_summary_complete, per_site_mu = mutation_summary_complete$mutation_rate/mutation_summary_complete$genome_size)
head(mutation_summary_complete)


write.csv(mutation_summary_complete, "mutation_summary_data_full_3416_unicell_revertants_T35_T2600_T6026_filtered_06DEC21.csv")

##############################
### START FROM HERE FOR PLOTTING 
##############################

mutation_summary_complete <- read.csv("mutation_summary_data_full_3416_unicell_revertants_T35_T2600_T6026_filtered_06DEC21.csv", header=T, stringsAsFactors=F)
mutation_summary_complete <- mutation_summary_complete[,-1]
head(mutation_summary_complete)


ggplot(data=mutation_summary_complete[mutation_summary_complete$counts>=1,], aes(x=1-per_site_mu, y=avg_workload)) +
  geom_point(alpha=0.8, cex=3, aes(pch=factor(timepoint)), col="black") +
  geom_smooth(method='lm', formula=y~x, col="dodgerblue") +
  scale_shape_discrete(name="Generation") +
  xlab(expression(paste("Fidelity (1"-"µ)"))) +
  ylab("Average workload") +
  ggtitle("")


p1 <- ggplot(data=mutation_summary_complete[mutation_summary_complete$counts>0,], aes(x=1-per_site_mu, y=avg_workload)) +
  geom_point(alpha=0.8, cex=3, aes(pch=factor(timepoint), col=counts)) +
  scale_shape_discrete(name="Generation") +
  scale_color_viridis_c(option="B", name = "# of cell \ndivisions") +
  geom_smooth(method='lm', formula=y~x, col="dodgerblue") +
  xlab(expression(paste("Fidelity (1"-"µ)"))) +
  ylab("Average workload") +
  ggtitle("# cell divisions > 0") +
  theme_cowplot() +
  theme(legend.position="none")

p2 <- ggplot(data=mutation_summary_complete[mutation_summary_complete$counts>1,], aes(x=1-per_site_mu, y=avg_workload)) +
  geom_point(alpha=0.8, cex=3, aes(pch=factor(timepoint), col=counts)) +
  scale_shape_discrete(name="Generation") +
  scale_color_viridis_c(option="B", name = "# of cell \ndivisions") +
  geom_smooth(method='lm', formula=y~x, col="dodgerblue") +
  xlab(expression(paste("Fidelity (1"-"µ)"))) +
  ylab("Average workload") +
  ggtitle("# cell divisions > 1")+
  theme_cowplot() +
  theme(legend.position="none")

p3 <- ggplot(data=mutation_summary_complete[mutation_summary_complete$counts>2,], aes(x=1-per_site_mu, y=avg_workload)) +
  geom_point(alpha=0.8, cex=3, aes(col=counts, pch=factor(timepoint))) +
  scale_color_viridis_c(option="B", name = "# of cell \ndivisions") +
  scale_shape_discrete(name="Generation") +
  geom_smooth(method='lm', formula=y~x, col="dodgerblue") +
  xlab(expression(paste("Fidelity (1"-"µ)"))) +
  ylab("Average workload") +
  ggtitle("# cell divisions > 2") +
  theme_cowplot()

ggdraw() +
  draw_plot(p1, x=0, y=0.0, width=0.3, height=1) +
  draw_plot(p2, x=0.3, y=0.0, width=0.3, height=1) + 
  draw_plot(p3, x=0.6, y=0.0, width=0.4, height=1)


# LOOK AT GROWTH RATE

ggplot(data=mutation_summary_complete[mutation_summary_complete$counts>=1,], aes(x=1-per_site_mu, y=avg_workload)) +
  geom_point(alpha=0.8, cex=3, aes(pch=factor(timepoint), col=growth_rate)) +
  geom_smooth(method='lm', formula=y~x, col="dodgerblue") +
  scale_color_viridis_c(option="B", name = "Growth \nrate") +
  scale_shape_discrete(name="Generation") +
  xlab(expression(paste("Fidelity (1"-"µ)"))) +
  ylab("Average workload") +
  ggtitle("") +
  facet_wrap(~timepoint, ncol=1, scales="free")

p1 <- ggplot(data=mutation_summary_complete[mutation_summary_complete$counts>0,], aes(x=1-per_site_mu, y=avg_workload)) +
  geom_point(alpha=0.8, cex=3, aes(pch=factor(timepoint), col=growth_rate)) +
  scale_shape_discrete(name="Generation") +
  scale_color_viridis_c(option="B", name = "Growth \nrate") +
  geom_smooth(method='lm', formula=y~x, col="dodgerblue") +
  xlab(expression(paste("Fidelity (1"-"µ)"))) +
  ylab("Average workload") +
  ggtitle("# cell divisions > 0") +
  theme_cowplot() +
  theme(legend.position="none") +
  facet_wrap(~timepoint, ncol=1, scales="free")

p2 <- ggplot(data=mutation_summary_complete[mutation_summary_complete$counts>1,], aes(x=1-per_site_mu, y=avg_workload)) +
  geom_point(alpha=0.8, cex=3, aes(pch=factor(timepoint), col=growth_rate)) +
  scale_shape_discrete(name="Generation") +
  scale_color_viridis_c(option="B", name = "Growth \nrate") +
  geom_smooth(method='lm', formula=y~x, col="dodgerblue") +
  xlab(expression(paste("Fidelity (1"-"µ)"))) +
  ylab("Average workload") +
  ggtitle("# cell divisions > 1") +
  theme_cowplot() +
  theme(legend.position="none") +
  facet_wrap(~timepoint, ncol=1, scales="free")

p3 <- ggplot(data=mutation_summary_complete[mutation_summary_complete$counts>2,], aes(x=1-per_site_mu, y=avg_workload)) +
  geom_point(alpha=0.8, cex=3, aes(pch=factor(timepoint), col=growth_rate)) +
  scale_color_viridis_c(option="B", name = "Growth \nrate") +
  scale_shape_discrete(name="Generation") +
  geom_smooth(method='lm', formula=y~x, col="dodgerblue") +
  xlab(expression(paste("Fidelity (1"-"µ)"))) +
  ylab("Average workload") +
  ggtitle("# cell divisions > 2") +
  theme_cowplot() +
  facet_wrap(~timepoint, ncol=1, scales="free")

ggdraw() +
  draw_plot(p1, x=0, y=0.0, width=0.3, height=1) +
  draw_plot(p2, x=0.3, y=0.0, width=0.3, height=1) + 
  draw_plot(p3, x=0.6, y=0.0, width=0.4, height=1)




uni_id <- c()
timepoint <- c()
num_reps <- c()
counts <- c()
avg_workload <- c()
avg_growth_rate <- c()
per_site_mu <- c()

for (i in unique(mutation_summary_complete$uni_id)) {
  uni_id[length(uni_id)+1] <- i
  timepoint[length(timepoint)+1] <- unique(mutation_summary_complete$timepoint[mutation_summary_complete$uni_id==i])
  num_reps[length(num_reps)+1] <- length(mutation_summary_complete$counts[mutation_summary_complete$uni_id==i])
  counts[length(counts)+1] <- mean(mutation_summary_complete$counts[mutation_summary_complete$uni_id==i])
  avg_workload[length(avg_workload)+1] <- mean(mutation_summary_complete$avg_workload[mutation_summary_complete$uni_id==i])
  avg_growth_rate[length(avg_growth_rate)+1] <- mean(mutation_summary_complete$growth_rate[mutation_summary_complete$uni_id==i])
  per_site_mu[length(per_site_mu)+1] <- mean(mutation_summary_complete$per_site_mu[mutation_summary_complete$uni_id==i])
}

final_df <- cbind.data.frame(uni_id, timepoint, num_reps, counts, avg_workload, avg_growth_rate, per_site_mu)

ggplot(data=final_df[final_df$counts>0,], aes(x=1-per_site_mu, y=avg_workload))+
  geom_point(alpha=1, cex=2.5, aes(pch=factor(timepoint), col=avg_growth_rate)) +
  scale_shape_discrete(name="Generation") +
  scale_color_viridis_c(option="B", name = "Revertant\ngrowth\nrate") +
  geom_smooth(method='lm', formula=y~x, col="dodgerblue") +
  xlab(expression(paste("Fidelity (1"-"µ)"))) +
  ylab("Average workload") +
  theme_cowplot(font_size=20) +
  facet_wrap(~timepoint, ncol=1, scales="free")

ggplot(data=final_df[final_df$counts>0,], aes(x=1-per_site_mu, y=avg_workload))+
  geom_point(alpha=0.6, cex=3, aes(pch=factor(timepoint), col=avg_growth_rate)) +
  scale_shape_discrete(name="timepoint") +
  scale_color_viridis_c(option="B", name = "# of cell \ndivisions") +
  xlab(expression(paste("Fidelity (1"-"µ)"))) +
  ylab("Average workload") +
  theme_cowplot(font_size=20) +
  xlim(0.9895, 0.9915) +
  ylim(0, 1) +
  facet_wrap(~timepoint, ncol=3, scales="free")


p1 <- ggplot(data=final_df[final_df$counts>0,], aes(x=1-per_site_mu, y=avg_workload))+
  geom_point(aes(pch=factor(timepoint), cex=counts, col=avg_growth_rate), alpha=0.5) +
  scale_shape_discrete(name="Generation") +
  scale_color_viridis_c(option="B", name = "# of cell \ndivisions") +
  geom_smooth(method='lm', formula=y~x, col="dodgerblue") +
  xlab(expression(paste("Fidelity (1"-"µ)"))) +
  ylab("Average workload") +
  ggtitle("# cell divisions > 0") +
  theme_cowplot(font_size = 20) +
  theme(legend.position="none") +
  facet_wrap(~timepoint, ncol=1, scales="free")

p2 <- ggplot(data=final_df[final_df$counts>1,], aes(x=1-per_site_mu, y=avg_workload))+
  geom_point(aes(pch=factor(timepoint), cex=counts, col=avg_growth_rate), alpha=0.5) +
  scale_shape_discrete(name="Generation") +
  scale_color_viridis_c(option="B", name = "# of cell \ndivisions") +
  geom_smooth(method='lm', formula=y~x, col="dodgerblue") +
  xlab(expression(paste("Fidelity (1"-"µ)"))) +
  ylab("Average workload") +
  ggtitle("# cell divisions > 1") +
  theme_cowplot(font_size=20) +
  theme(legend.position="none") +
  facet_wrap(~timepoint, ncol=1, scales="free")

p3 <- ggplot(data=final_df[final_df$counts>2,], aes(x=1-per_site_mu, y=avg_workload))+
  geom_point(aes(pch=factor(timepoint), cex=counts, col=avg_growth_rate), alpha=0.5) +
  scale_shape_discrete(name="Generation") +
  scale_size_continuous(name="# iterations") +
  scale_color_viridis_c(option="B", name = "# of cell \ndivisions") +
  geom_smooth(method='lm', formula=y~x, col="dodgerblue") +
  xlab(expression(paste("Fidelity (1"-"µ)"))) +
  ylab("Average workload") +
  ggtitle("# cell divisions > 2") +
  theme_cowplot(font_size=20) +
  facet_wrap(~timepoint, ncol=1, scales="free")

ggdraw() +
  draw_plot(p1, x=0, y=0.0, width=0.3, height=1) +
  draw_plot(p2, x=0.3, y=0.0, width=0.3, height=1) + 
  draw_plot(p3, x=0.6, y=0.0, width=0.4, height=1)






super_hi_fidelity_revertants <- final_df[(1 - final_df$per_site_mu) > 0.993,]
revertant_genome <- c()

for (i in super_hi_fidelity_revertants$uni_id) {
  revertant_genome[length(revertant_genome)+1] <- genomes$propagule_genome[genomes$uni_id==i & genomes$iteration==1 & 
                                                                             genomes$update==min(genomes$update[genomes$uni_id==i & genomes$iteration==1])]
}

super_hi_fidelity_revertants_full <- cbind.data.frame(super_hi_fidelity_revertants, revertant_genome)
head(super_hi_fidelity_revertants_full)
setwd("~/Desktop/Dropbox/Desktop_28NOV21/Entrenchment_processed_data_for_figures/Figure_5/Figure_5_data/multis")
mc_genomes_t35 <- read.csv("35_detail_mc.csv", header=T, stringsAsFactors=F)
mc_genomes_t6026 <- read.csv("6026_detail_mc.csv", header=T, stringsAsFactors=F)

head(mc_genomes_t35$replicate)

multi_anc_genomes <- cbind.data.frame(timepoint=c(35, 6026), multi_anc_genome=c(mc_genomes_t35$propagule_genome[1], mc_genomes_t6026$propagule_genome[1]))

complete_data_SHFRs <- left_join(super_hi_fidelity_revertants_full, multi_anc_genomes, by="timepoint")
head(complete_data_SHFRs)


mut_position <- c()
multi_allele <- c()
revertant_allele <- c()

for (i in complete_data_SHFRs$uni_id) {
  rev_gen <- as.numeric(unlist(strsplit(complete_data_SHFRs$revertant_genome[complete_data_SHFRs$uni_id==i], "\\s+")))
  rev_gen <- rev_gen[!is.na(rev_gen)]
  anc_gen <- as.numeric(unlist(strsplit(complete_data_SHFRs$multi_anc_genome[complete_data_SHFRs$uni_id==i], "\\s+")))
  anc_gen <- anc_gen[!is.na(anc_gen)]
  glo <- align(list(list(anc_gen), list(rev_gen)), type = "global", seqweights=c(1,0), progressive=T)
  
  temp <- c()
  for (p in 1:length(glo[1,])) {
    
    if (glo[1,p] == glo[2,p]) {
      temp[length(temp)+1] <- TRUE
    }
    
    # else if (glo[2,p] == "-" | glo[1,p] == "-") {
    #   temp[length(temp)+1] <- NA
    #   }
    
    else if (glo[1,p] != glo[2,p]) {
      #uni_id[length(uni_id)+1] <- h
      #timepoint[length(timepoint)+1] <- genomes$timepoint[genomes$uni_id==h & genomes$iteration==i & genomes$index==j]
      #iteration[length(iteration)+1] <- i
      #index[length(index)+1] <- j
      #parent_org[length(parent_org)+1] <- genomes$parent_id[genomes$uni_id==h & genomes$iteration==i & genomes$index==j]
      #offspring_org[length(offspring_org)+1] <- genomes$offspring_id[genomes$uni_id==h & genomes$iteration==i & genomes$index==j]
      mut_position[length(mut_position)+1] <- length(glo[1,1:p][glo[1,1:p]!="-"])
      multi_allele[length(multi_allele)+1] <- glo[1,p]
      revertant_allele[length(revertant_allele)+1] <- glo[2,p]
    }
  }
  
}

final_SHFR_mutation_df <- cbind.data.frame(complete_data_SHFRs, mut_position, multi_allele, revertant_allele)




check_genome_length <- genomes[genomes$uni_id %in% final_SHFR_mutation_df$uni_id,]
offspring_genome_lengths <- unlist(lapply(check_genome_length$offspring_genome, FUN = function(x) {length(unlist(strsplit(x, "\\s+")))}))
check_genome_length <- cbind.data.frame(check_genome_length, offspring_genome_lengths)

mean_offspring_genome_length <- c()

for (i in unique(check_genome_length$uni_id)) {
  mean_offspring_genome_length[length(mean_offspring_genome_length)+1] <- mean(check_genome_length$offspring_genome_lengths[check_genome_length$uni_id==i])
}

df_to_add <- cbind.data.frame(uni_id = unique(check_genome_length$uni_id), mean_offspring_genome_length)

final_SHFR_mutation_df <- left_join(final_SHFR_mutation_df, df_to_add, by="uni_id")

write.csv(final_SHFR_mutation_df, "super_high_fidelity_mutations_data_07DEC21.csv")
