# Get multicell genomes from transition and final timepoint
# Get reversion mutations
# 26 November 2021
# Peter Conlin


# Get mutlicell genomes from transition timepoint
setwd("~/Desktop/Dropbox/Desktop 28NOV21/Entrenchment_new_data/004e_new/multis_trans/")

# reads in all files in the working directory and makes a data frame for each
filenames <- list.files(pattern='*mc.csv')
all_files <- lapply(filenames, function(i){read.csv(i, header=T, stringsAsFactors=F)})
genomes <- bind_rows(all_files, .id = "column_label")
head(genomes)

anc_genomes_df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("strategy_rep", "timepoint", "genome_location", "wt_mut"))

for (i in unique(genomes$replicate)) {
  anc_genome <- as.numeric(unlist(strsplit(genomes$propagule_genome[genomes$replicate==i & genomes$X==min(genomes$X[genomes$replicate==i])], "\\s+")))
  anc_genome <- anc_genome[!is.na(anc_genome)]
  anc_genomes_df <- rbind.data.frame(anc_genomes_df, data.frame(strain_id=rep(i, length(anc_genome)),
                                     timepoint=rep(0, length(anc_genome)),
                                     position=c(1:length(anc_genome)),
                                     genome=anc_genome))
}
anc_genomes_df


# Get unicell genomes from transition timepoint
setwd("~/Desktop/Dropbox/Desktop 28NOV21/Entrenchment_new_data/004e_new/unis_trans/")

# reads in all files in the working directory and makes a data frame for each
filenames <- list.files(pattern='*uni.csv')
all_files <- lapply(filenames, function(i){read.csv(i, header=T, stringsAsFactors=F)})
genomes <- bind_rows(all_files, .id = "column_label")
head(genomes)

first_propagules <- genomes[!1:length(genomes$update),]

for (i in unique(genomes$replicate)) {
  for (j in unique(genomes$uni[genomes$replicate==i])) {
    first_propagules <- rbind.data.frame(first_propagules, genomes[genomes$replicate==i & genomes$uni==i])
  }
}


anc_genomes_df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("strategy_rep", "timepoint", "genome_location", "wt_mut"))

for (i in unique(genomes$replicate)) {
  anc_genome <- as.numeric(unlist(strsplit(genomes$propagule_genome[genomes$replicate==i & genomes$X==min(genomes$X[genomes$replicate==i])], "\\s+")))
  anc_genome <- anc_genome[!is.na(anc_genome)]
  anc_genomes_df <- rbind.data.frame(anc_genomes_df, data.frame(strain_id=rep(i, length(anc_genome)),
                                                                timepoint=rep(0, length(anc_genome)),
                                                                position=c(1:length(anc_genome)),
                                                                genome=anc_genome))
}
anc_genomes_df



# Get mutlicell genomes from final timepoint
setwd("~/Desktop/Dropbox/Desktop 28NOV21/Entrenchment_new_data/004e_new/multis_end/")

# reads in all files in the working directory and makes a data frame for each
filenames <- list.files(pattern='*mc.csv')
all_files <- lapply(filenames, function(i){read.csv(i, header=T, stringsAsFactors=F)})
genomes <- bind_rows(all_files, .id = "column_label")

for (i in unique(genomes$replicate)) {
  anc_genome <- as.numeric(unlist(strsplit(genomes$propagule_genome[genomes$replicate==i & genomes$X==min(genomes$X[genomes$replicate==i])], "\\s+")))
  anc_genome <- anc_genome[!is.na(anc_genome)]
  anc_genomes_df <- rbind.data.frame(anc_genomes_df, data.frame(strain_id=rep(i, length(anc_genome)),
                                                                timepoint=rep(1, length(anc_genome)),
                                                                position=c(1:length(anc_genome)),
                                                                genome=anc_genome))
}
head(anc_genomes_df)
colnames(anc_genomes_df) <- c("strategy_rep", "timepoint", "genome_location", "wt_mut")

setwd("~/Desktop/Dropbox/Desktop 28NOV21/Entrenchment_processed_data_for_figures/SI_Figure_Reversion_mutation/Data/")
write.csv(anc_genomes_df, "multicelled_ancestor_genomes_at_transition_and_final_timepoint_29NOV21.csv")
