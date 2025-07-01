# Getting the time at which multicellularity evolved and the end point of the experiment (in generations)
# 10 November 2022
# Peter Conlin

setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/entrenchment_data/2020-04-01_update_lod_fitness/")

trans <- read.csv("entrenchment_data_transition_15JULY20.csv", header=T, stringsAsFactors=F)
unique(trans$rep)
final <- read.csv("entrenchment_data_final_15JULY20.csv", header=T, stringsAsFactors=F)
unique(final$rep)

head(trans)
trans$gen <- as.numeric(trans$generation) - as.numeric(trans$generation_diff)
final$gen <- as.numeric(final$generation) - as.numeric(final$generation_diff)



strategy_rep <- c()
timepoint <- c()
generation <- c()

# Get transition time point data
for (i in unique(trans$rep)) {
  strategy_rep[length(strategy_rep)+1] <- i
  timepoint[length(timepoint)+1] <- 0
  generation[length(generation)+1] <- as.integer(median(trans$gen[trans$rep==i & trans$gen > 0]))
}

# Get final time point data
for (i in unique(final$rep)) {
  strategy_rep[length(strategy_rep)+1] <- i
  timepoint[length(timepoint)+1] <- 1
  generation[length(generation)+1] <- as.integer(median(final$gen[final$rep==i & final$gen > 0]))
}


trans[trans$rep==4178,]

df <- cbind.data.frame(strategy_rep, timepoint, generation)

setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/")
write.csv(df, "Time_of_transition_and_final_timepoint_in_generations_10NOV22.csv")
