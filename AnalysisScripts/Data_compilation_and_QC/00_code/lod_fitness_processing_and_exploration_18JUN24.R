#Plotting entrenchment data
#Peter Conlin
#28 September 2022

#----------------------------------------
# PREPARE WORKSPACE
# Load packages
library(ggplot2)
library(gridExtra)
library(GGally)
library(reshape)
library(dplyr)
library(Hmisc)
library(cowplot)
library(forcats)
theme_set(theme_cowplot(font_size = 14))
library(viridis)
library(egg)

one_col <- 2.25
two_col <- 4.75

# Set working directory
setwd("/Users/peterconlin/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/00_Data_compilation_and_QC")
lod_fitness_filtered <- read.csv("lod_fitness_full_with_duplicates_filtered_18JUN24.csv", header=T, stringsAsFactors=F)
lod_fitness_filtered_temp <- lod_fitness_filtered[,-1]
head(lod_fitness_filtered_temp)

length(unique(interaction(lod_fitness_filtered_temp$strategy_rep, lod_fitness_filtered_temp$timepoint, lod_fitness_filtered_temp$mc_or_uni)))

# lod_fitness_filtered_temp_trans <- lod_fitness_filtered_temp[lod_fitness_filtered_temp$timepoint==0,]
# lod_fitness_filtered_trans <- lod_fitness_filtered_temp_trans[!(lod_fitness_filtered_temp_trans$strategy_rep %in% c(4310)),]
# lod_fitness_filtered_temp_final <- lod_fitness_filtered_temp[lod_fitness_filtered_temp$timepoint==1,]
# lod_fitness_filtered_final <- lod_fitness_filtered_temp_final[!(lod_fitness_filtered_temp_final$strategy_rep %in% c(3038, 4278, 4310)),]
# unique(lod_fitness_filtered_final$strategy_rep)
# lod_fitness_filtered_full <- rbind.data.frame(lod_fitness_filtered_temp, lod_fitness_filtered_final)

lod_fitness_trans_reruns <- read.csv("lod_fitness_trans_reruns.csv", header=T, stringsAsFactors=F)
unique(lod_fitness_trans_reruns$strategy_rep)
lod_fitness_final_reruns <- read.csv("lod_fitness_final_reruns.csv", header=T, stringsAsFactors=F)
unique(lod_fitness_final_reruns$strategy_rep)
lod_fitness_reruns <- rbind.data.frame(lod_fitness_trans_reruns, lod_fitness_final_reruns)
head(lod_fitness_reruns)

lod_fitness_reruns_2 <- cbind.data.frame(lod_fitness_reruns, 
                                         new_rev_mut_ID = rep(NA, length(lod_fitness_reruns$X)),
                                         index = rep(NA, length(lod_fitness_reruns$X))
                                         )
                                         
colnames(lod_fitness_reruns_2)[4] <- "reversion_mut_ID"
colnames(lod_fitness_reruns_2)[8] <- "num_orgs"
head(lod_fitness_reruns_2)


# Remove old data for strategy_reps that have been rerun
lod_fitness_filtered_temp2 <- lod_fitness_filtered_temp[!(lod_fitness_filtered_temp$strategy_rep %in% unique(lod_fitness_reruns$strategy_rep)),]
head(lod_fitness_filtered_temp2)


# Merge new data with old
lod_fitness <- rbind.data.frame(lod_fitness_filtered_temp2, lod_fitness_reruns_2)
length(unique(interaction(lod_fitness$strategy_rep,lod_fitness$timepoint)))

head(lod_fitness)
lod_fitness[which(lod_fitness$timepoint==60), ]$timepoint <- 0 # Rewrite mislabeled timepoint for isolate 4178
lod_fitness$mc_or_uni[lod_fitness$mc_or_uni=="multicell"] <- "mc"
lod_fitness$mc_or_uni[lod_fitness$mc_or_uni=="unicell"] <- "uni"
unique(lod_fitness$mc_or_uni)

# #----------------------------------------
# # REFORMATTING SUMMARY DATA ON ENTRENCHMENT
# 
# # Read in summary data 
# lod_fitness_summary <- read.csv("lod_fitness_summary.csv", header=T, stringsAsFactors=F)
# lod_fitness_summary <- lod_fitness_summary[,-1]
# # Rewrite mislabeled timepoint for isolate 4178
# lod_fitness_summary[92,1] <- 0
# colnames(lod_fitness_summary)[3] <- "num_viable_unicells_unfiltered"
# head(lod_fitness_summary)

#----------------------------------------
# PREPARING DATA ON THE DISTRIBUTION OF FITNESS EFFECTS OF REVERSION MUTATIONS

# # Read in lod fitness data
# lod_fitness_final <- read.csv("lod_fitness_end.csv", header=T, stringsAsFactors=F)
# lod_fitness_trans <- read.csv("lod_fitness_trans.csv", header=T, stringsAsFactors=F)
# lod_fitness_trans[which(lod_fitness_trans$timepoint==60), ]$timepoint <- 0 # Rewrite mislabeled timepoint for isolate 4178
# unique(lod_fitness_trans$strategy_rep)
# lod_fitness <- rbind(lod_fitness_trans, lod_fitness_final)
# head(lod_fitness)

# # Check for missing data in lod_fitness_trans
# genotype_ID <- paste(lod_fitness_trans$strategy_rep, lod_fitness_trans$mc_or_uni, lod_fitness_trans$timepoint, lod_fitness_trans$count, sep = "_")
# test_df <- cbind.data.frame(lod_fitness_trans, genotype_ID)
# length(unique(test_df$genotype_ID))
# old_trans_genotype_IDs <- unique(test_df$genotype_ID)
# 
# 
# genotype_ID_temp <- c()
# count_reps <- c()
# for (i in unique(test_df$genotype_ID)) {
#   genotype_ID_temp[length(genotype_ID_temp)+1] <- i
#   count_reps[length(count_reps)+1] <- length(test_df$count[test_df$genotype_ID==i])
# }
# rep_count_data <- cbind.data.frame(genotype_ID=genotype_ID_temp, count_reps)
# missing_data_trans <- rep_count_data[rep_count_data$count_reps!=100,]
# head(missing_data_trans)
# 
# # Check for missing data in lod_fitness_trans_new
# genotype_ID <- paste(lod_fitness_trans_new$strategy_rep, lod_fitness_trans_new$mc_or_uni, lod_fitness_trans_new$timepoint, lod_fitness_trans_new$count, sep = "_")
# test_df <- cbind.data.frame(lod_fitness_trans_new, genotype_ID)
# length(unique(test_df$genotype_ID))
# 
# genotype_ID_temp <- c()
# count_reps <- c()
# for (i in unique(test_df$genotype_ID)) {
#   genotype_ID_temp[length(genotype_ID_temp)+1] <- i
#   count_reps[length(count_reps)+1] <- length(test_df$count[test_df$genotype_ID==i])
# }
# rep_count_data <- cbind.data.frame(genotype_ID=genotype_ID_temp, count_reps)
# missing_data_trans_new <- rep_count_data[rep_count_data$count_reps!=100,]
# 
# 
# write.csv(missing_data, "missing_lod_fitness_data.csv")

# Split mc and uni into separate dataframes
# THIS MC DATAFRAME IS MISSING 17 VALUES?!?! FUCK FUCK FUCK THIS IS SO FUCKING TIRING
# lod_fitness_mc <- lod_fitness[which(lod_fitness$mc_or_uni=="mc"),]
# length(unique(interaction(lod_fitness_mc$strategy_rep,lod_fitness_mc$timepoint)))
# lod_fitness_mc$timepoint_lab <- rep(NA, length(lod_fitness_mc$timepoint))
# lod_fitness_mc$timepoint_lab[lod_fitness_mc$timepoint==0] <- "Transition"
# lod_fitness_mc$timepoint_lab[lod_fitness_mc$timepoint==1] <- "Final"
# lod_fitness_mc$timepoint_lab <- factor(lod_fitness_mc$timepoint_lab, levels=c("Transition", "Final"))

# FOR NOW, WILL USE OLD MC DATA...
setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/Figure_3/Figure_3_data")
old_lod_fitness_final <- read.csv("lod_fitness_end.csv", header=T, stringsAsFactors=F)
old_lod_fitness_trans <- read.csv("lod_fitness_trans.csv", header=T, stringsAsFactors=F)
old_lod_fitness_trans[which(old_lod_fitness_trans$timepoint==60), ]$timepoint <- 0 # Rewrite mislabeled timepoint for isolate 4178
length(unique(old_lod_fitness_trans$strategy_rep))
old_lod_fitness <- rbind(old_lod_fitness_trans, old_lod_fitness_final)
head(old_lod_fitness)
lod_fitness_mc <- old_lod_fitness[old_lod_fitness$mc_or_uni=="mc",]
length(unique(interaction(lod_fitness_mc$strategy_rep,lod_fitness_mc$timepoint)))

lod_fitness_uni <- lod_fitness[which(lod_fitness$mc_or_uni=="uni"),]
length(unique(interaction(lod_fitness_uni$strategy_rep,lod_fitness_uni$timepoint)))

lod_fitness_uni$timepoint_lab <- rep(NA, length(lod_fitness_uni$timepoint))
lod_fitness_uni$timepoint_lab[lod_fitness_uni$timepoint==0] <- "Transition"
lod_fitness_uni$timepoint_lab[lod_fitness_uni$timepoint==1] <- "Final"
lod_fitness_uni$timepoint_lab <- factor(lod_fitness_uni$timepoint_lab, levels=c("Transition", "Final"))

#length(lod_fitness_mc$X) + length(lod_fitness_uni$X) == length(lod_fitness$X)


# Get data on number of replicate measurements (iterations)
strategy_rep <- c()
timepoint <- c()
num_counts <- c()

for (i in unique(lod_fitness_uni$strategy_rep)) {
  for (j in unique(lod_fitness_uni$timepoint[lod_fitness_uni$strategy_rep == i])) {
    strategy_rep[length(strategy_rep)+1] <- i
    timepoint[length(timepoint)+1] <- j
    num_counts[length(num_counts)+1] <- length(unique(lod_fitness_uni$reversion_mut_ID[lod_fitness_uni$strategy_rep == i & lod_fitness_uni$timepoint == j]))
  }
}


lod_fitness_summary <- cbind.data.frame(strategy_rep, timepoint, num_viable_unicells_unfiltered=num_counts)
head(lod_fitness_summary)

#----------------------------------------
# FILTER 1: REMAINING UNICELLULAR

# # Check what proportion of unicellular revertants re-evolve multicellularity during our fitness assay
# S4_inset <- ggplot(data=lod_fitness_uni, aes(x=(total_cells/num_orgs), fill=(total_cells/num_orgs)<=1)) +
#   geom_histogram(binwidth=1, breaks=c(0:25), col="black") +
#   scale_fill_manual(values=c("#78C0A8", "black"), labels=c("Fail", "Pass")) +
#   scale_y_continuous() +
#   guides(fill = guide_legend(reverse=TRUE)) +
#   labs(x="Average organism size", 
#        y="Count", 
#        fill="Data \nfiltering \nresult") +
#   theme(text = element_text(size = 10))
# 
# #+ 
#   #facet_wrap(~ timepoint, ncol = 1, nrow = 2)
# 
# 
# # Plot the proportion of each genotype that re-evolves multicellularity
# setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Figures")
# 
# ggplot(data=lod_fitness_uni, aes(x=factor(strategy_rep), fill=(total_cells/num_orgs)<=1)) +
#   geom_histogram(stat="count", col=NA) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_fill_manual(values=c("#78C0A8", "black"), labels=c("Fail", "Pass")) +
#   guides(fill = guide_legend(reverse=TRUE)) +
#   labs(x="Average organism size \n(Total cells / # of organisms)", 
#        y="Count", 
#        fill="Data \nfiltering \nresult") +
#   facet_wrap(~ timepoint_lab, ncol = 1, nrow = 2)
# 
# #table(lod_fitness_uni$reversion_mut_ID)
# 
# #lod_fitness_uni

genotype_ID_1 <- paste(lod_fitness_uni$strategy_rep, lod_fitness_uni$mc_or_uni, lod_fitness_uni$timepoint, lod_fitness_uni$reversion_mut_ID, sep = "_")
genotype_ID_counts <- as.data.frame(table(genotype_ID_1))
genotype_ID_counts[genotype_ID_counts$Freq>100,]

test_df_1 <- cbind.data.frame(lod_fitness_uni, genotype_ID_1)
#head(test_df_1)
#table(interaction(test_df_1$strategy_rep, test_df_1$reversion_mut_ID))

# S4_A <- ggplot(data=test_df_1[which(test_df_1$timepoint==0),], aes(x=factor(genotype_ID_1), fill=(total_cells/num_orgs)<=1)) +
#   geom_histogram(stat="count", col=NA) +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) +
#   scale_fill_manual(values=c("#78C0A8", "black"), labels=c("Fail", "Pass")) +
#   guides(fill = guide_legend(reverse=TRUE)) +
#   labs(x="Unicellular revertant genotype (transition)", 
#        y="Count", 
#        fill="Data \nfiltering \nresult") +
#   facet_wrap(~ timepoint_lab, ncol = 1, nrow = 2) +
#   theme(legend.position = "none")
# 
# S4_B <- ggplot(data=test_df_1[which(test_df_1$timepoint==1),], aes(x=factor(genotype_ID_1), fill=(total_cells/num_orgs)<=1)) +
#   geom_histogram(stat="count", col=NA) +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) +
#   scale_fill_manual(values=c("#78C0A8", "black"), labels=c("Fail", "Pass")) +
#   guides(fill = guide_legend(reverse=TRUE)) +
#   labs(x="Unicellular revertant genotype (final)", 
#        y="Count", 
#        fill="Data \nfiltering \nresult") +
#   facet_wrap(~ timepoint_lab, ncol = 1, nrow = 2) +
#   theme(legend.position = "none")
# 
# 
# pdf("Figure_S4_Dirty_work_Populations_that_re-evolve_multicellularity_24APR24.pdf", width=one_col*3, height=one_col*2.5)
# ggdraw() +
#   draw_plot(S4_inset, x=0.39, y=0.85, width=0.25, height=0.15) +
#   draw_plot(S4_A, x = 0, y = 0.435, 
#             width = length(unique(test_df_1$genotype_ID_1[test_df_1$timepoint==0])) / max(c(length(unique(test_df_1$genotype_ID_1[test_df_1$timepoint==0])),length(unique(test_df_1$genotype_ID_1[test_df_1$timepoint==1])))), 
#             height = 0.415) +
#   draw_plot(S4_B, x = 0, y = 0, 
#             width = length(unique(test_df_1$genotype_ID_1[test_df_1$timepoint==1])) / max(c(length(unique(test_df_1$genotype_ID_1[test_df_1$timepoint==0])),length(unique(test_df_1$genotype_ID_1[test_df_1$timepoint==1])))), 
#             height = 0.435) 
# dev.off()


# Lenient filtering step 1: remove all replicates where multicellularity re-evolves
lod_fitness_uni_filtered_1 <- lod_fitness_uni[which((lod_fitness_uni$total_cells / lod_fitness_uni$num_orgs) <= 1),]
head(lod_fitness_uni_filtered_1)

# ggplot(data=lod_fitness_uni_filtered_1, aes(x=factor(strategy_rep), fill=(total_cells/num_orgs)<=1)) +
#   geom_histogram(stat="count", col=NA) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   #facet_wrap(~ timepoint, ncol = 1, nrow = 2) +
#   scale_fill_manual(values=c("grey70", "black"), labels=c("Fail", "Pass")) +
#   scale_y_continuous() +
#   guides(fill = guide_legend(reverse=TRUE)) +
#   labs(x="Average organism size \n(Total cells / # of organisms)", 
#        y="Count", 
#        fill="Data \nfiltering \nresult")


test_df_1_mod <- cbind.data.frame(test_df_1, avg_cells_per_org = test_df_1$total_cells/test_df_1$num_orgs)
head(test_df_1_mod)
example <- test_df_1_mod[test_df_1_mod$genotype_ID_1=="3002_uni_0_244",]

genotype_ID_summary <- c()
F1_pass_percentage <- c()
F2_pass_percentage <- c()
CUMULATIVE_pass_percentage <- c()

for (i in unique(test_df_1_mod$genotype_ID_1)) {
  genotype_ID_summary[length(genotype_ID_summary)+1] <- i
  F1_pass_percentage[length(F1_pass_percentage)+1] <- length(test_df_1_mod$avg_cells_per_org[test_df_1_mod$genotype_ID_1==i & test_df_1_mod$avg_cells_per_org==1.])
  F2_pass_percentage[length(F2_pass_percentage)+1] <- length(test_df_1_mod$avg_cells_per_org[test_df_1_mod$genotype_ID_1==i & test_df_1_mod$num_orgs>1])
  CUMULATIVE_pass_percentage[length(CUMULATIVE_pass_percentage)+1] <- length(test_df_1_mod$avg_cells_per_org[test_df_1_mod$genotype_ID_1==i & test_df_1_mod$avg_cells_per_org==1. & test_df_1_mod$num_orgs>1])
}

data_filtering_df <- cbind.data.frame(genotype=c(sapply(genotype_ID_summary, function(X) strsplit(X, split="_")[[1]][1])), 
                                           timepoint=c(sapply(genotype_ID_summary, function(X) strsplit(X, split="_")[[1]][3])),
                                           uni_id=c(sapply(genotype_ID_summary, function(X) strsplit(X, split="_")[[1]][4])),
                                           F1_pass_percentage, F2_pass_percentage, CUMULATIVE_pass_percentage)
head(data_filtering_df, 5)

setwd("/Users/peterconlin/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/00_Data_compilation_and_QC")
write.csv(data_filtering_df, "quality_control_filtering_summary_data_18JUN24.csv")

length(unique(data_filtering_df$genotype))

# ggplot(data_filtering_df, aes(x=as.factor(genotype), y=F1_pass_percentage)) +
#   geom_point(stat="summary", fun.y="mean", color="red") +
#   geom_point(aes(as.factor(genotype), y=F2_pass_percentage), stat="summary", fun.y="mean", color="dodgerblue")

paste(c(strsplit(genotype_ID_summary[1], split="_")[[1]][4:5]), sep="_")

df_mean_std <- data_filtering_df %>%
  group_by(genotype, timepoint) %>%
  summarise_at(vars(F1_pass_percentage, F2_pass_percentage, CUMULATIVE_pass_percentage), list(mean=mean, sd=sd, length=length)) %>% 
  as.data.frame()

head(df_mean_std)
length(unique(df_mean_std$genotype))

df_mean_std <- cbind.data.frame(df_mean_std, F1_pass_percentage_se = (df_mean_std$F1_pass_percentage_sd * (1/sqrt(df_mean_std$F1_pass_percentage_length))), 
                 F2_pass_percentage_se = (df_mean_std$F2_pass_percentage_sd * (1/sqrt(df_mean_std$F2_pass_percentage_length))),
                 CUMULATIVE_pass_percentage_se = (df_mean_std$CUMULATIVE_pass_percentage_sd * (1/sqrt(df_mean_std$CUMULATIVE_pass_percentage_length))))

length(unique(df_mean_std$genotype))

write.csv(df_mean_std, "quality_control_filtering_summary_data_means_18JUN24.csv")


data_filtering_df_NO_ZEROS <- data_filtering_df[data_filtering_df$CUMULATIVE_pass_percentage!=0,]

df_mean_std_NO_ZEROS <- data_filtering_df_NO_ZEROS %>%
  group_by(genotype, timepoint) %>%
  summarise_at(vars(F1_pass_percentage, F2_pass_percentage, CUMULATIVE_pass_percentage), list(mean=mean, sd=sd, length=length)) %>% 
  as.data.frame()

df_mean_std_NO_ZEROS <- cbind.data.frame(df_mean_std_NO_ZEROS, F1_pass_percentage_se = (df_mean_std_NO_ZEROS$F1_pass_percentage_sd * (1/sqrt(df_mean_std_NO_ZEROS$F1_pass_percentage_length))), 
                                F2_pass_percentage_se = (df_mean_std_NO_ZEROS$F2_pass_percentage_sd * (1/sqrt(df_mean_std_NO_ZEROS$F2_pass_percentage_length))),
                                CUMULATIVE_pass_percentage_se = (df_mean_std_NO_ZEROS$CUMULATIVE_pass_percentage_sd * (1/sqrt(df_mean_std_NO_ZEROS$CUMULATIVE_pass_percentage_length))))

write.csv(df_mean_std_NO_ZEROS, "quality_control_filtering_summary_data_means_NO_ZEROS_18JUN24.csv")


# transition_plot <- ggplot(df_mean_std[df_mean_std$timepoint==0,], aes(x=fct_reorder(as.factor(genotype), F2_pass_percentage_mean), y=F1_pass_percentage_mean)) +
#   geom_errorbar(aes(ymin=F1_pass_percentage_mean-F1_pass_percentage_se, ymax=F1_pass_percentage_mean+F1_pass_percentage_se), width=.3, color="red") + #, position = position_nudge(x = -0.1)
#   geom_errorbar(aes(ymin=F2_pass_percentage_mean-F2_pass_percentage_se, ymax=F2_pass_percentage_mean+F2_pass_percentage_se), width=.3, color="dodgerblue") + #, position = position_nudge(x = +0.1)
#   geom_errorbar(aes(ymin=CUMULATIVE_pass_percentage_mean-CUMULATIVE_pass_percentage_se, ymax=CUMULATIVE_pass_percentage_mean+CUMULATIVE_pass_percentage_se), width=.3, color="darkorchid4") + #, position = position_nudge(x = +0.1)
#   geom_point(color="red", pch=19) + #, position = position_nudge(x = -0.1)
#   geom_point(aes(x=fct_reorder(as.factor(genotype), F2_pass_percentage_mean), y=F2_pass_percentage_mean), color="dodgerblue", pch=17)  + #, position = position_nudge(x = +0.1)
#   geom_point(aes(x=fct_reorder(as.factor(genotype), CUMULATIVE_pass_percentage_mean), y=CUMULATIVE_pass_percentage_mean), color="darkorchid4", pch=17)  + #, position = position_nudge(x = +0.1)
#   ylim(0,100) +
#   ylab("% Replicates that pass filter") +
#   xlab("Multicellular progenitor genotype (transition)") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#   
# # SOMETHING IS BUGGED WITH THIS PLOT...
# # transition_plot <- ggplot(df_mean_std[df_mean_std$timepoint==0,], aes(x=fct_reorder(as.factor(genotype), F2_pass_percentage_mean), y=F1_pass_percentage_mean)) +
# #   geom_bar(aes(x=fct_reorder(as.factor(genotype), CUMULATIVE_pass_percentage_mean), y=rep(100, 65)), stat="identity", fill="dodgerblue") +
# #   geom_bar(aes(x=fct_reorder(as.factor(genotype), CUMULATIVE_pass_percentage_mean), y=F2_pass_percentage_mean), stat="identity", fill="red") +
# #   geom_bar(aes(x=fct_reorder(as.factor(genotype), CUMULATIVE_pass_percentage_mean), y=CUMULATIVE_pass_percentage_mean), stat="identity", fill="grey20") +
# #   #geom_errorbar(aes(ymin=F1_pass_percentage_mean-F1_pass_percentage_se, ymax=F1_pass_percentage_mean+F1_pass_percentage_se), width=.3, color="white") + #, position = position_nudge(x = -0.1)
# #   geom_errorbar(aes(ymin=F2_pass_percentage_mean-F2_pass_percentage_se, ymax=F2_pass_percentage_mean+F2_pass_percentage_se), width=.3, color="white") + #, position = position_nudge(x = +0.1)
# #   geom_errorbar(aes(ymin=CUMULATIVE_pass_percentage_mean-CUMULATIVE_pass_percentage_se, ymax=CUMULATIVE_pass_percentage_mean+CUMULATIVE_pass_percentage_se), width=.3, color="white") + #, position = position_nudge(x = +0.1)
# #   ylim(0,100) +
# #   ylab("% Replicates that pass filter") +
# #   xlab("Multicellular progenitor genotype (transition)") +
# #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
#   
# final_plot <- ggplot(df_mean_std[df_mean_std$timepoint==1,], aes(x=fct_reorder(as.factor(genotype), F2_pass_percentage_mean), y=F1_pass_percentage_mean)) +
#   geom_errorbar(aes(ymin=F1_pass_percentage_mean-F1_pass_percentage_se, ymax=F1_pass_percentage_mean+F1_pass_percentage_se), width=.3, color="red") + #, position = position_nudge(x = -0.1)
#   geom_errorbar(aes(ymin=F2_pass_percentage_mean-F2_pass_percentage_se, ymax=F2_pass_percentage_mean+F2_pass_percentage_se), width=.3, color="dodgerblue") + #, position = position_nudge(x = +0.1)
#   geom_errorbar(aes(ymin=CUMULATIVE_pass_percentage_mean-CUMULATIVE_pass_percentage_se, ymax=CUMULATIVE_pass_percentage_mean+CUMULATIVE_pass_percentage_se), width=.3, color="darkorchid4") + #, position = position_nudge(x = +0.1)
#   geom_point(color="red", pch=19) + #, position = position_nudge(x = -0.1)
#   geom_point(aes(x=fct_reorder(as.factor(genotype), F2_pass_percentage_mean), y=F2_pass_percentage_mean), color="dodgerblue", pch=17) + #, position = position_nudge(x = +0.1)
#   geom_point(aes(x=fct_reorder(as.factor(genotype), CUMULATIVE_pass_percentage_mean), y=CUMULATIVE_pass_percentage_mean), color="darkorchid4", pch=17)  + #, position = position_nudge(x = +0.1)
#   ylim(0,100) +
#   ylab("% Replicates that pass filter") +
#   xlab("Multicellular progenitor genotype (final)") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# final_plot <- ggplot(df_mean_std[df_mean_std$timepoint==1,], aes(x=fct_reorder(as.factor(genotype), F2_pass_percentage_mean), y=F1_pass_percentage_mean)) +
#   geom_bar(aes(x=fct_reorder(as.factor(genotype), CUMULATIVE_pass_percentage_mean), y=rep(100, 65)), stat="identity", fill="dodgerblue") +
#   geom_bar(aes(x=fct_reorder(as.factor(genotype), CUMULATIVE_pass_percentage_mean), y=F2_pass_percentage_mean), stat="identity", fill="red") +
#   geom_bar(aes(x=fct_reorder(as.factor(genotype), CUMULATIVE_pass_percentage_mean), y=CUMULATIVE_pass_percentage_mean), stat="identity", fill="grey20") +
#   #geom_errorbar(aes(ymin=F1_pass_percentage_mean-F1_pass_percentage_se, ymax=F1_pass_percentage_mean+F1_pass_percentage_se), width=.3, color="white") + #, position = position_nudge(x = -0.1)
#   geom_errorbar(aes(ymin=F2_pass_percentage_mean-F2_pass_percentage_se, ymax=F2_pass_percentage_mean+F2_pass_percentage_se), width=.3, color="white") + #, position = position_nudge(x = +0.1)
#   geom_errorbar(aes(ymin=CUMULATIVE_pass_percentage_mean-CUMULATIVE_pass_percentage_se, ymax=CUMULATIVE_pass_percentage_mean+CUMULATIVE_pass_percentage_se), width=.3, color="white") + #, position = position_nudge(x = +0.1)
#   ylim(0,100) +
#   ylab("% Replicates that pass filter") +
#   xlab("Multicellular progenitor genotype (final)") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# 
# ggdraw() +
#   draw_plot(transition_plot, x=0.01, y=0.5, width=0.99, height=0.5) +
#   draw_plot(final_plot, x = 0.01, y = 0.0, width =0.99, height = 0.5)



###################
# Set working directory
setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/Figure_6/Figure_6_data")

# Read in main experimental data
avg_lod_fitness_data <- read.csv("avg_lod_fitness_data_filtered_2_14FEB22.csv", header=T, stringsAsFactors=F)
head(avg_lod_fitness_data)
full_entrenchment_data <- read.csv("full_entrenchment_data_filtered_2_14FEB22.csv", header=T, stringsAsFactors=F)
head(full_entrenchment_data)
#lod_fitness_summary <- read.csv("lod_fitness_summary_F2_adjusted_14FEB22.csv", header=T, stringsAsFactors=F)

# Read in division of labor data
avg_dol_summary <- read.csv("avg_dol_summary_14FEB22.csv", header=T, stringsAsFactors=F)
head(avg_dol_summary)
trans_dol <- avg_dol_summary[avg_dol_summary$timepoint==0,]
colnames(trans_dol)[length(trans_dol)] <- "trans_dol"
trans_dol <- trans_dol[,c(2,3,5)]
final_dol <- avg_dol_summary[avg_dol_summary$timepoint==1,]
colnames(final_dol)[length(final_dol)] <- "final_dol"
final_dol <- final_dol[,c(2,3,5)]
dol_values <- full_join(trans_dol, final_dol, by=c("strategy_rep"))
dol_values <- dol_values[,c(1,3,5)]
head(dol_values)

# Reformat DoL data
avg_dol_endpoints <- avg_dol_summary[avg_dol_summary$timepoint==1,]
differentiation_status <- cbind.data.frame(strategy_rep=avg_dol_endpoints$strategy_rep, has_soma=avg_dol_endpoints$has_soma)
head(differentiation_status, 30)
colnames(differentiation_status)[1] <- "genotype"
length(differentiation_status$genotype)


###################

setwd("/Users/peterconlin/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/00_Data_compilation_and_QC")
df_mean_std <- read.csv("quality_control_filtering_summary_data_means_18JUN24.csv", header=T, stringsAsFactors=T)
#df_mean_std <- read.csv("quality_control_filtering_summary_data_means_NO_ZEROS_08MAR23.csv", header=T, stringsAsFactors=T)

df_mean_std <- left_join(df_mean_std, differentiation_status, by="genotype")
head(df_mean_std)

length(unique(interaction(df_mean_std$genotype, df_mean_std$timepoint)))


trans_F1_df <- df_mean_std[df_mean_std$timepoint==0,]
length(unique(trans_F1_df$genotype))

final_F1_df <- df_mean_std[df_mean_std$timepoint==1,]
length(unique(final_F1_df$genotype))

trans_F1_df <- cbind.data.frame(trans_F1_df[order(trans_F1_df$F1_pass_percentage_mean),], index=c(1:65))
final_F1_df <- cbind.data.frame(final_F1_df[order(final_F1_df$F1_pass_percentage_mean),], index=c(1:65))
F1_df <- rbind.data.frame(trans_F1_df, final_F1_df)
head(F1_df)

F1_plot <- ggplot(F1_df, aes(x=index, y=F1_pass_percentage_mean, color=as.factor(timepoint), pch=as.factor(timepoint))) +
  geom_point() +
  geom_errorbar(aes(ymin=F1_pass_percentage_mean-F1_pass_percentage_se, ymax=F1_pass_percentage_mean+F1_pass_percentage_se), width=.3) + #, position = position_nudge(x = -0.1)
  scale_color_manual(values=c("darkgoldenrod1", "green3"), labels=c("Transition", "Final"), name="Time point") +
  scale_shape_manual(values=c(19, 17), labels=c("Transition", "Final"), name="Time point") +
  theme() +
  ylab("% Replicates unicellular") +
  xlab("Rank (filter 1)") +
  ylim(c(0,100))
  
trans_F2_df <- df_mean_std[df_mean_std$timepoint==0,]
trans_F2_df <- cbind.data.frame(trans_F2_df[order(trans_F2_df$F2_pass_percentage_mean),], index=c(1:65))
final_F2_df <- df_mean_std[df_mean_std$timepoint==1,]
final_F2_df <- cbind.data.frame(final_F2_df[order(final_F2_df$F2_pass_percentage_mean),], index=c(1:65))
F2_df <- rbind.data.frame(trans_F2_df, final_F2_df)

F2_plot <- ggplot(F2_df, aes(x=index, y=F2_pass_percentage_mean, color=as.factor(timepoint), pch=as.factor(timepoint))) +
  geom_point() +
  geom_errorbar(aes(ymin=F2_pass_percentage_mean-F2_pass_percentage_se, ymax=F2_pass_percentage_mean+F2_pass_percentage_se), width=.3) + #, position = position_nudge(x = -0.1)
  scale_color_manual(values=c("darkgoldenrod1", "green3"), labels=c("Transition", "Final"), name="Time point") +
  scale_shape_manual(values=c(19, 17), labels=c("Transition", "Final"), name="Time point") +
  ylab("% Replicates viable") +
  xlab("Rank (filter 2)") +
  ylim(c(0,100))


ggplot(F2_df[F2_df$timepoint==1,], aes(x=index, y=F2_pass_percentage_mean, color=as.factor(has_soma), pch=as.factor(timepoint))) +
  geom_point() +
  geom_errorbar(aes(ymin=F2_pass_percentage_mean-F2_pass_percentage_se, ymax=F2_pass_percentage_mean+F2_pass_percentage_se), width=.3) + #, position = position_nudge(x = -0.1)
  scale_color_manual(values=c("dodgerblue", "red"), labels=c("Transition", "Final"), name="Time point") +
  scale_shape_manual(values=c(19, 17), labels=c("Transition", "Final"), name="Time point") +
  ylab("% Replicates viable") +
  xlab("Rank (filter 2)") +
  ylim(c(0,100))


trans_CUMULATIVE_df <- df_mean_std[df_mean_std$timepoint==0,]
trans_CUMULATIVE_df <- cbind.data.frame(trans_CUMULATIVE_df[order(trans_CUMULATIVE_df$CUMULATIVE_pass_percentage_mean),], index=c(1:65))
final_CUMULATIVE_df <- df_mean_std[df_mean_std$timepoint==1,]
final_CUMULATIVE_df <- cbind.data.frame(final_CUMULATIVE_df[order(final_CUMULATIVE_df$CUMULATIVE_pass_percentage_mean),], index=c(1:65))
CUMULATIVE_df <- rbind.data.frame(trans_CUMULATIVE_df, final_CUMULATIVE_df)

CUM_plot <- ggplot(CUMULATIVE_df, aes(x=index, y=CUMULATIVE_pass_percentage_mean, color=as.factor(timepoint), pch=as.factor(timepoint))) +
  geom_point() +
  geom_errorbar(aes(ymin=CUMULATIVE_pass_percentage_mean-CUMULATIVE_pass_percentage_se, ymax=CUMULATIVE_pass_percentage_mean+CUMULATIVE_pass_percentage_se), width=.3) + #, position = position_nudge(x = -0.1)
  scale_color_manual(values=c("darkgoldenrod1", "green3"), labels=c("Transition", "Final"), name="Time point") +
  scale_shape_manual(values=c(19, 17), labels=c("Transition", "Final"), name="Time point") +
  theme(legend.position="blank") +
  ylab("% Replicates (cumulative)") +
  xlab("Rank (filter 1 + filter 2)") +
  ylim(c(0,100))


ggplot(CUMULATIVE_df[CUMULATIVE_df$timepoint==1,], aes(x=index, y=CUMULATIVE_pass_percentage_mean, color=as.factor(has_soma), pch=as.factor(timepoint))) +
  geom_point() +
  geom_errorbar(aes(ymin=CUMULATIVE_pass_percentage_mean-CUMULATIVE_pass_percentage_se, ymax=CUMULATIVE_pass_percentage_mean+CUMULATIVE_pass_percentage_se), width=.3) + #, position = position_nudge(x = -0.1)
  scale_color_manual(values=c("dodgerblue", "red"), labels=c("Transition", "Final"), name="Time point") +
  scale_shape_manual(values=c(19, 17), labels=c("Transition", "Final"), name="Time point") +
  ylab("% Replicates viable") +
  xlab("Rank (filter 2)") +
  ylim(c(0,100))


# Save SI figures
setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Figures")


# Check what proportion of unicellular revertants re-evolve multicellularity during our fitness assay
S4_inset <- ggplot(data=lod_fitness_uni, aes(x=(total_cells/num_orgs), fill=(total_cells/num_orgs)<=1)) +
  geom_histogram(binwidth=1, breaks=c(0:25), col=NA) +
  scale_fill_manual(values=c("grey70", "black"), labels=c("Fail", "Pass")) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE), breaks=c(0,400000,800000)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  labs(x="Average organism size", 
       y="Count", 
       fill="Data \nfiltering \nresult") +
  theme(text = element_text(size = 12), panel.border = element_rect(colour = "black", fill=NA))


pdf("Figure_S4_Dirty_work_Populations_that_re-evolve_multicellularity_18JUN24.pdf", width=one_col*3.25, height=one_col*1.75)
ggdraw() +
  draw_plot(F1_plot, x=0.01, y=0.01, width=0.99, height=0.99) +
  draw_plot(S4_inset, x = 0.4, y = 0.15, width=0.425, height = 0.4)
dev.off()



S5_inset <- ggplot(data=lod_fitness_uni, aes(x=num_orgs, fill=(num_orgs) > 1)) +
  geom_histogram(stat="count", col=NA) +
  #facet_wrap(~ timepoint, ncol = 1, nrow = 2) +
  scale_fill_manual(values=c("grey70", "black"), labels=c("Fail", "Pass")) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  guides(fill = guide_legend(reverse=TRUE)) +
  labs(x="Number of organisms", 
       y="Count", 
       fill="Data \nfiltering \nresult") +
  theme(text = element_text(size = 12), panel.border = element_rect(colour = "black", fill=NA))



pdf("Figure_S5_Dirty_work_Populations_that_fail_to_divide_18JUN24.pdf", width=one_col*3.25, height=one_col*1.75)
ggdraw() +
  draw_plot(F2_plot, x=0.01, y=0.01, width=0.99, height=0.99) +
  draw_plot(S5_inset, x = 0.4, y = 0.15, width=0.425, height = 0.4)
dev.off()




ggarrange(F1_plot, F2_plot, CUM_plot)






### NOT WORKING STARTING HERE
# Update lod_fitness_summary table with number of replicates that pass our filters
post_filtering_lod_fitness_df <- lod_fitness_uni_filtered_1
strategy_rep <- c()
timepoint <- c()
num_counts <- c()

for (i in unique(post_filtering_lod_fitness_df$strategy_rep)) {
  for (j in unique(post_filtering_lod_fitness_df$timepoint[post_filtering_lod_fitness_df$strategy_rep == i])) {
    strategy_rep[length(strategy_rep)+1] <- i
    timepoint[length(timepoint)+1] <- j
    num_counts[length(num_counts)+1] <- length(unique(post_filtering_lod_fitness_df$reversion_mut_ID[post_filtering_lod_fitness_df$strategy_rep == i & post_filtering_lod_fitness_df$timepoint == j]))
  }
}

num_viable_unicells_filtered <- cbind.data.frame(strategy_rep, timepoint, num_counts)
colnames(num_viable_unicells_filtered)[3] <- "num_viable_unicells_F1"
head(num_viable_unicells_filtered)


lod_fitness_summary_F1 <- left_join(lod_fitness_summary, num_viable_unicells_filtered, by=c("strategy_rep", "timepoint"))
head(lod_fitness_summary_F1)


# Get data on number of replicate measurements (iterations) that pass our F1 filter
strategy_rep <- c()
timepoint <- c()
rev_mut_ID <- c()
num_iterations <- c()

for (i in unique(lod_fitness_uni_filtered_1$strategy_rep)) {
  for (j in unique(lod_fitness_uni_filtered_1$timepoint[lod_fitness_uni_filtered_1$strategy_rep == i])) {
    for (k in unique(lod_fitness_uni_filtered_1$reversion_mut_ID[lod_fitness_uni_filtered_1$strategy_rep == i & lod_fitness_uni_filtered_1$timepoint == j])) {
      strategy_rep[length(strategy_rep)+1] <- i
      timepoint[length(timepoint)+1] <- j
      rev_mut_ID[length(rev_mut_ID)+1] <- k
      num_iterations[length(num_iterations)+1] <- length(unique(lod_fitness_uni_filtered_1$iteration[lod_fitness_uni_filtered_1$strategy_rep == i & lod_fitness_uni_filtered_1$timepoint == j & lod_fitness_uni_filtered_1$reversion_mut_ID == k]))
    }
  }
}

F1_filtered_lod_fitness_summary_data <- cbind.data.frame(strategy_rep, timepoint, rev_mut_ID, num_iterations)
colnames(F1_filtered_lod_fitness_summary_data)[4] <- "num_iterations_F1"
head(F1_filtered_lod_fitness_summary_data)


#----------------------------------------
# FILTER 2: SUCCESSFUL REPRODUCTION
# Plot the proportion of each genotype that fails to divide at least once
genotype_ID <- paste(lod_fitness_uni_filtered_1$strategy_rep, lod_fitness_uni_filtered_1$count, sep = "_")
test_df_2 <- cbind.data.frame(lod_fitness_uni_filtered_1, genotype_ID)

# # in aggregate
# S5_inset <- ggplot(data=test_df_2, aes(x=num_orgs, fill=(num_orgs) > 1)) +
#   geom_histogram(stat="count", col="black") +
#   #facet_wrap(~ timepoint, ncol = 1, nrow = 2) +
#   scale_fill_manual(values=c("white", "black"), labels=c("Fail", "Pass")) +
#   scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
#   guides(fill = guide_legend(reverse=TRUE)) +
#   labs(x="Number of organisms",
#        y="Count",
#        fill="Data \nfiltering \nresult") +
#   theme(panel.border = element_rect(colour = "black", fill=NA))
# 
# 
# 
# # by strategy rep
# ggplot(data=test_df_1, aes(x=factor(strategy_rep), fill=(num_orgs) > 1)) +
#   geom_histogram(data=test_df_2, stat="count", col=NA) +
#   facet_wrap(~ timepoint_lab, ncol = 1, nrow = 2) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_fill_manual(values=c("#78C0A8", "black"), labels=c("Fail", "Pass")) +
#   scale_y_continuous() +
#   guides(fill = guide_legend(reverse=TRUE)) +
#   labs(x="Number of organisms",
#        y="Count",
#        fill="Data \nfiltering \nresult")
# 
# 
# # by unicell revertant genotype
# S5_A <- ggplot(data=test_df_1[which(test_df_1$timepoint==0),], aes(x=factor(genotype_ID), fill=(num_orgs) > 1)) +
#   geom_histogram(data=test_df_2, stat="count", col=NA) +
#   facet_wrap(~ timepoint_lab, ncol = 1, nrow = 2) +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) +
#   scale_fill_manual(values=c("#78C0A8", "black"), labels=c("Fail", "Pass")) +
#   guides(fill = guide_legend(reverse=TRUE)) +
#   labs(x="Unicellular revertant genotype (transition)",
#        y="Count",
#        fill="Data \nfiltering \nresult") +
#   facet_wrap(~ timepoint_lab, ncol = 1, nrow = 2) +
#   theme(legend.position = "none")
# 
# 
# S5_B <- ggplot(data=test_df_1[which(test_df_1$timepoint==1),], aes(x=factor(genotype_ID), fill=(num_orgs) > 1)) +
#   geom_histogram(data=test_df_2, stat="count", col=NA) +
#   facet_wrap(~ timepoint_lab, ncol = 1, nrow = 2) +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) +
#   scale_fill_manual(values=c("#78C0A8", "black"), labels=c("Fail", "Pass")) +
#   guides(fill = guide_legend(reverse=TRUE)) +
#   labs(x="Unicellular revertant genotype (final)",
#        y="Count",
#        fill="Data \nfiltering \nresult") +
#   facet_wrap(~ timepoint_lab, ncol = 1, nrow = 2) +
#   theme(legend.position = "none")
# 
# 
# pdf("Figure_S4_Dirty_work_Populations_that_fail_to_divide_12OCT22.pdf", width=one_col*3, height=one_col*2.5)
# ggdraw() +
#   draw_plot(S5_inset, x=0.39, y=0.85, width=0.25, height=0.15) +
#   draw_plot(S5_A, x = 0, y = 0.435, 
#             width = length(unique(test_df$genotype_ID[test_df$timepoint==0])) / max(c(length(unique(test_df$genotype_ID[test_df$timepoint==0])),length(unique(test_df$genotype_ID[test_df$timepoint==1])))), 
#             height = 0.415) +
#   draw_plot(S5_B, x = 0, y = 0, 
#             width = length(unique(test_df$genotype_ID[test_df$timepoint==1])) / max(c(length(unique(test_df$genotype_ID[test_df$timepoint==0])),length(unique(test_df$genotype_ID[test_df$timepoint==1])))), 
#             height = 0.435) 
# dev.off()

# Lenient filtering step 2: remove ALL replicates that failed to divide in the given time
lod_fitness_uni_filtered_2 <- lod_fitness_uni_filtered_1[which(lod_fitness_uni_filtered_1$num_orgs > 1),]
head(lod_fitness_uni_filtered_2)

#----------------------------------------
# # Lenient filtering 1/2 hybrid: remove replicates that UNIFORMLY failed to divide in the given time
# # Get list of sterile revertants
# sterile_revertant_genotypes <- c()
# for (i in unique(test_df$genotype_ID)) {
#   if (mean(test_df$num_orgs[test_df$genotype_ID==i]) <= 1) {
#     sterile_revertant_genotypes[length(sterile_revertant_genotypes)+1] <- i
#   }
# }
# 
# #Remove sterile revertants from full dataset
# for (i in unique(sterile_revertant_genotypes)){
#   test_df <- test_df[! (test_df$genotype_ID == i), ]
# }
# 
# test_df <- test_df[,-11]
# lod_fitness_uni_filtered_2 <- test_df
# 
# # Replot to check that the correct values have been removed
# # by unicell revertant genotype
# P1 <- ggplot(data=test_df[which(test_df$timepoint==0),], aes(x=factor(genotype_ID), fill=(num_orgs) > 1)) +
#   geom_histogram(stat="count", col=NA) +
#   facet_wrap(~ timepoint, ncol = 1, nrow = 2) +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.title.x = element_blank()) +
#   scale_fill_manual(values = c("#F07818", "#78C0A8"))
# P2 <- ggplot(data=test_df[which(test_df$timepoint==1),], aes(x=factor(genotype_ID), fill=(num_orgs) > 1)) +
#   geom_histogram(stat="count", col=NA) +
#   facet_wrap(~ timepoint, ncol = 1, nrow = 2) +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()) +
#   scale_fill_manual(values = c("#F07818", "#78C0A8"))
# ggdraw() +
#   draw_plot(P1, x = 0, y = 0.5, width = 1, height = 0.5) +
#   draw_plot(P2, x = 0, y = 0, width = 1, height = 0.5) 
# 


# #----------------------------------------
# ggplot(data=lod_fitness_uni_filtered_2, aes(x=factor(strategy_rep), fill=(num_orgs) > 1)) +
#   geom_histogram(stat="count", col=NA) +
#   facet_wrap(~ timepoint, ncol = 1, nrow = 2) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_fill_manual(values = c("#78C0A8","#F07818"))


# Update lod_fitness_summary table with number of replicates that pass our filters
post_filtering_lod_fitness_df <- lod_fitness_uni_filtered_2
strategy_rep <- c()
timepoint <- c()
num_counts <- c()

for (i in unique(post_filtering_lod_fitness_df$strategy_rep)) {
  for (j in unique(post_filtering_lod_fitness_df$timepoint[post_filtering_lod_fitness_df$strategy_rep == i])) {
    strategy_rep[length(strategy_rep)+1] <- i
    timepoint[length(timepoint)+1] <- j
    num_counts[length(num_counts)+1] <- length(unique(post_filtering_lod_fitness_df$reversion_mut_ID[post_filtering_lod_fitness_df$strategy_rep == i & post_filtering_lod_fitness_df$timepoint == j]))
  }
}

num_viable_unicells_filtered <- cbind.data.frame(strategy_rep, timepoint, num_counts)
colnames(num_viable_unicells_filtered)[3] <- "num_viable_unicells_F2"
head(num_viable_unicells_filtered)

lod_fitness_summary_F2 <- full_join(lod_fitness_summary_F1, num_viable_unicells_filtered, by=c("strategy_rep", "timepoint"))
#lod_fitness_summary_F1$num_viable_unicells_unfiltered == lod_fitness_summary_F2$num_viable_unicells_F1
head(lod_fitness_summary_F2)

#(length(lod_fitness_uni_filtered_1$X) - length(lod_fitness_uni_filtered_2$X)) / length(lod_fitness_uni_filtered_1$X)


#----------------------------------------
head(lod_fitness_uni)
head(lod_fitness_uni_filtered_1)
head(lod_fitness_uni_filtered_2)

# Get data on number of replicate measurements (iterations) that pass our F2 filter
strategy_rep <- c()
timepoint <- c()
rev_mut_ID <- c()
num_iterations <- c()

for (i in unique(lod_fitness_uni_filtered_2$strategy_rep)) {
  for (j in unique(lod_fitness_uni_filtered_2$timepoint[lod_fitness_uni_filtered_2$strategy_rep == i])) {
    for (k in unique(lod_fitness_uni_filtered_2$reversion_mut_ID[lod_fitness_uni_filtered_2$strategy_rep == i & lod_fitness_uni_filtered_2$timepoint == j])) {
      #print(paste(i, j, k, length(unique(lod_fitness_uni_filtered_2$iteration[lod_fitness_uni_filtered_2$strategy_rep == i & lod_fitness_uni_filtered_2$timepoint == j & lod_fitness_uni_filtered_2$count == k]))))
      strategy_rep[length(strategy_rep)+1] <- i
      timepoint[length(timepoint)+1] <- j
      rev_mut_ID[length(rev_mut_ID)+1] <- k
      num_iterations[length(num_iterations)+1] <- length(unique(lod_fitness_uni_filtered_2$iteration[lod_fitness_uni_filtered_2$strategy_rep == i & lod_fitness_uni_filtered_2$timepoint == j & lod_fitness_uni_filtered_2$reversion_mut_ID == k]))
    }
  }
}

F2_filtered_lod_fitness_summary_data <- cbind.data.frame(strategy_rep, timepoint, rev_mut_ID, num_iterations)
colnames(F2_filtered_lod_fitness_summary_data)[4] <- "num_iterations_F2"
head(F2_filtered_lod_fitness_summary_data)


#----------------------------------------
# Calculate adjusted number of viable unicellular revertants
combined_lod_fitness_summary_iterations <- full_join(F1_filtered_lod_fitness_summary_data, F2_filtered_lod_fitness_summary_data, by=c("strategy_rep", "timepoint", "rev_mut_ID"))
combined_lod_fitness_summary_iterations$proportion_fertile <- combined_lod_fitness_summary_iterations$num_iterations_F2 / combined_lod_fitness_summary_iterations$num_iterations_F1
head(combined_lod_fitness_summary_iterations)


strategy_rep <- c()
timepoint <- c()
adjusted_num_viable_uni <- c()

for (i in unique(combined_lod_fitness_summary_iterations$strategy_rep)) {
  for (j in unique(combined_lod_fitness_summary_iterations$timepoint[combined_lod_fitness_summary_iterations$strategy_rep == i])) {
    strategy_rep[length(strategy_rep)+1] <- i
    timepoint[length(timepoint)+1] <- j
    adjusted_num_viable_uni[length(adjusted_num_viable_uni)+1] <- sum(combined_lod_fitness_summary_iterations$proportion_fertile[combined_lod_fitness_summary_iterations$strategy_rep == i & combined_lod_fitness_summary_iterations$timepoint == j], na.rm=TRUE)
  }
}

adjusted_viable_unicells_df <- cbind.data.frame(strategy_rep, timepoint, adjusted_num_viable_uni)
head(adjusted_viable_unicells_df)
lod_fitness_summary_F2_adjusted <- full_join(lod_fitness_summary_F2, adjusted_viable_unicells_df, by=c("strategy_rep", "timepoint"))
head(lod_fitness_summary_F2_adjusted)


setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/Figure_3/Figure_3_data")
write.csv(lod_fitness_summary_F2_adjusted, "lod_fitness_summary_F2_adjusted_18JUN24.csv")



# #----------------------------------------
# # FILTER 3 (OPTIONAL): SUFFICIENT REPLICATION
# 
# # Plot number of replicates after filtering for each unicell genotype
# plot(x=c(1:length(F2_filtered_lod_fitness_summary_data$num_iterations)), y=F2_filtered_lod_fitness_summary_data$num_iterations, 
#      xlab="Unicell revertant genotype", ylab="# Replicates after filtering", type="h")
# num_iterations_sd <- sd(F2_filtered_lod_fitness_summary_data$num_iterations)
# xx <- c(0, length(F2_filtered_lod_fitness_summary_data$num_iterations), length(F2_filtered_lod_fitness_summary_data$num_iterations), 0)
# yy <- c(mean(F2_filtered_lod_fitness_summary_data$num_iterations)-num_iterations_sd, 
#         mean(F2_filtered_lod_fitness_summary_data$num_iterations)-num_iterations_sd,
#         mean(F2_filtered_lod_fitness_summary_data$num_iterations)+num_iterations_sd,
#         mean(F2_filtered_lod_fitness_summary_data$num_iterations)+num_iterations_sd)
# polygon(xx, yy, col=rgb(1,0,0,0.2), border=NA)
# segments(x0=0, x1=length(F2_filtered_lod_fitness_summary_data$num_iterations), 
#          y0=mean(F2_filtered_lod_fitness_summary_data$num_iterations),
#          y1=mean(F2_filtered_lod_fitness_summary_data$num_iterations),
#          lwd=2, lty=2, col="red")
# 
# 
# # Optional stringent filtering step: remove all genotypes that have fewer replicates than average after applying the lenient filter
# #LENIENT FILTER
# mean_iterations <- mean(F2_filtered_lod_fitness_summary_data$num_iterations)
# sd_iterations <- sd(F2_filtered_lod_fitness_summary_data$num_iterations)
# insufficient_iterations <- F2_filtered_lod_fitness_summary_data[which(F2_filtered_lod_fitness_summary_data$num_iterations < mean_iterations-2*sd_iterations),]
# #length(insufficient_iterations$strategy_rep)
# 
# lod_fitness_uni_filtered_3L <- anti_join(lod_fitness_uni_filtered_2, insufficient_iterations, by = c("strategy_rep", "timepoint", "count"))
# length(lod_fitness_uni_filtered_3L$X)
# length(lod_fitness_uni_filtered_2$X) - length(lod_fitness_uni_filtered_3L$X) == sum(insufficient_iterations$num_iterations)
# 
# 
# # Update lod_fitness_summary table with number of replicates that pass our filters
# post_filtering_lod_fitness_df <- lod_fitness_uni_filtered_3L
# strategy_rep <- c()
# timepoint <- c()
# num_counts <- c()
# 
# for (i in unique(post_filtering_lod_fitness_df$strategy_rep)) {
#   for (j in unique(post_filtering_lod_fitness_df$timepoint[post_filtering_lod_fitness_df$strategy_rep == i])) {
#     strategy_rep[length(strategy_rep)+1] <- i
#     timepoint[length(timepoint)+1] <- j
#     num_counts[length(num_counts)+1] <- length(unique(post_filtering_lod_fitness_df$count[post_filtering_lod_fitness_df$strategy_rep == i & post_filtering_lod_fitness_df$timepoint == j]))
#   }
# }
# 
# num_viable_unicells_filtered <- cbind.data.frame(strategy_rep, timepoint, num_counts)
# colnames(num_viable_unicells_filtered)[3] <- "num_viable_unicells_F3L"
# head(num_viable_unicells_filtered)
# 
# lod_fitness_summary_F3L <- full_join(lod_fitness_summary_F2_adjusted, num_viable_unicells_filtered, by=c("strategy_rep", "timepoint"))
# lod_fitness_summary_F2$num_viable_unicells_unfiltered == lod_fitness_summary_F3L$num_viable_unicells_F2
# head(lod_fitness_summary_F3L)
# 
# 
# #MODERATE FILTER
# mean_iterations <- mean(F2_filtered_lod_fitness_summary_data$num_iterations)
# sd_iterations <- sd(F2_filtered_lod_fitness_summary_data$num_iterations)
# insufficient_iterations <- F2_filtered_lod_fitness_summary_data[which(F2_filtered_lod_fitness_summary_data$num_iterations < mean_iterations-sd_iterations),]
# #length(insufficient_iterations$strategy_rep)
# 
# lod_fitness_uni_filtered_3M <- anti_join(lod_fitness_uni_filtered_2, insufficient_iterations, by = c("strategy_rep", "timepoint", "count"))
# length(lod_fitness_uni_filtered_3M$X)
# length(lod_fitness_uni_filtered_2$X) - length(lod_fitness_uni_filtered_3M$X) == sum(insufficient_iterations$num_iterations)
# 
# 
# # Update lod_fitness_summary table with number of replicates that pass our filters
# post_filtering_lod_fitness_df <- lod_fitness_uni_filtered_3M
# strategy_rep <- c()
# timepoint <- c()
# num_counts <- c()
# 
# for (i in unique(post_filtering_lod_fitness_df$strategy_rep)) {
#   for (j in unique(post_filtering_lod_fitness_df$timepoint[post_filtering_lod_fitness_df$strategy_rep == i])) {
#     strategy_rep[length(strategy_rep)+1] <- i
#     timepoint[length(timepoint)+1] <- j
#     num_counts[length(num_counts)+1] <- length(unique(post_filtering_lod_fitness_df$count[post_filtering_lod_fitness_df$strategy_rep == i & post_filtering_lod_fitness_df$timepoint == j]))
#   }
# }
# 
# num_viable_unicells_filtered <- cbind.data.frame(strategy_rep, timepoint, num_counts)
# colnames(num_viable_unicells_filtered)[3] <- "num_viable_unicells_F3M"
# head(num_viable_unicells_filtered)
# 
# lod_fitness_summary_F3M <- full_join(lod_fitness_summary_F3L, num_viable_unicells_filtered, by=c("strategy_rep", "timepoint"))
# lod_fitness_summary_F3L$num_viable_unicells_unfiltered == lod_fitness_summary_F3M$num_viable_unicells_F3M
# head(lod_fitness_summary_F3M)
# 
# 
# #STRINGENT FILTER
# mean_iterations <- mean(F2_filtered_lod_fitness_summary_data$num_iterations)
# sd_iterations <- sd(F2_filtered_lod_fitness_summary_data$num_iterations)
# insufficient_iterations <- F2_filtered_lod_fitness_summary_data[which(F2_filtered_lod_fitness_summary_data$num_iterations < mean_iterations),]
# #length(insufficient_iterations$strategy_rep)
# 
# lod_fitness_uni_filtered_3S <- anti_join(lod_fitness_uni_filtered_2, insufficient_iterations, by = c("strategy_rep", "timepoint", "count"))
# length(lod_fitness_uni_filtered_3S$X)
# length(lod_fitness_uni_filtered_2$X) - length(lod_fitness_uni_filtered_3S$X) == sum(insufficient_iterations$num_iterations)
# 
# 
# # Update lod_fitness_summary table with number of replicates that pass our filters
# post_filtering_lod_fitness_df <- lod_fitness_uni_filtered_3S
# strategy_rep <- c()
# timepoint <- c()
# num_counts <- c()
# 
# for (i in unique(post_filtering_lod_fitness_df$strategy_rep)) {
#   for (j in unique(post_filtering_lod_fitness_df$timepoint[post_filtering_lod_fitness_df$strategy_rep == i])) {
#     strategy_rep[length(strategy_rep)+1] <- i
#     timepoint[length(timepoint)+1] <- j
#     num_counts[length(num_counts)+1] <- length(unique(post_filtering_lod_fitness_df$count[post_filtering_lod_fitness_df$strategy_rep == i & post_filtering_lod_fitness_df$timepoint == j]))
#   }
# }
# 
# num_viable_unicells_filtered <- cbind.data.frame(strategy_rep, timepoint, num_counts)
# colnames(num_viable_unicells_filtered)[3] <- "num_viable_unicells_F3S"
# head(num_viable_unicells_filtered)
# 
# lod_fitness_summary_F3S <- full_join(lod_fitness_summary_F3M, num_viable_unicells_filtered, by=c("strategy_rep", "timepoint"))
# lod_fitness_summary_F3M$num_viable_unicells_unfiltered == lod_fitness_summary_F3S$num_viable_unicells_F3M
# head(lod_fitness_summary_F3S)
# 
# lod_fitness_summary_F3S[which(lod_fitness_summary_F3S$strategy_rep==3097),]
# 
# 
# ggplot(data=lod_fitness_summary_F3S, aes(x=factor(strategy_rep))) +
#   geom_col(aes(y=num_viable_unicells_unfiltered), fill="#F07818", color=NA, alpha=0.25) +
#   geom_col(aes(y=num_viable_unicells_F1), fill="#F07818", color=NA, alpha=0.5) + 
#   geom_col(aes(y=num_viable_unicells_F2), fill="#F07818", color=NA, alpha=1) + 
#   geom_col(aes(y=num_viable_unicells_F3L), fill="white", color=NA, alpha=1) + 
#   geom_col(aes(y=num_viable_unicells_F3L), fill="#78C0A8", color=NA, alpha=0.25) + 
#   geom_col(aes(y=num_viable_unicells_F3M), fill="#78C0A8", color=NA, alpha=0.5) + 
#   geom_col(aes(y=num_viable_unicells_F3S), fill="#78C0A8", color=NA, alpha=1) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_y_continuous(trans='log10') +
#   xlab("Multicellular progenitor strain") +
#   ylab("Number of viable unicellular revertants") +
#   facet_wrap(~timepoint, ncol=1, nrow=2)


#----------------------------------------
# DATA ANALYSIS: PT 1

# Select which level of filtering to use (1, 2, 3S, 3M, 3L) and rejoin unicell and multicell data
colnames(lod_fitness_mc)[4] <- "reversion_mut_ID"
colnames(lod_fitness_mc)[8] <- "num_orgs"
intersect(colnames(lod_fitness_mc), colnames(lod_fitness_uni_filtered_2))
lod_fitness_mc$new_rev_mut_ID <- rep(NA, length(lod_fitness_mc$X))
lod_fitness_mc$index <- paste(lod_fitness_mc$strategy_rep, lod_fitness_mc$reversion_mut_ID, lod_fitness_mc$timepoint, sep="_")
lod_fitness_mc$timepoint_lab <- ifelse(lod_fitness_mc$timepoint==0, "Transition", "Final")


lod_fitness_filtered_full <- rbind(lod_fitness_mc, lod_fitness_uni_filtered_2)
lod_fitness_filtered_full$growth_rate <- log(lod_fitness_filtered_full$num_orgs / 1) / lod_fitness_filtered_full$time_to_fill
head(lod_fitness_filtered_full)

write.csv(lod_fitness_filtered_full, "lod_fitness_data_full_filtered_2_18JUN24.csv")


# Average avg_time_to_fill and avg_workload across replicate runs of each reversion mutation, extract ID variables
timepoint <- c()
mc_or_uni <- c()
rev_mut_ID <- c()
strategy_rep <- c()
avg_growth_rate <- c()
growth_rate_sd <- c()
avg_workload <- c()
workload_sd <- c()


# First, we average the effects ("time_to_fill" and "workload") of each reversion mutation (unique combination of "replicate" and "id") across all 100 "iterations"
for (i in unique(lod_fitness_filtered_full$strategy_rep)) {
  for (j in unique(lod_fitness_filtered_full$timepoint[lod_fitness_filtered_full$strategy_rep == i])) {
    for (k in unique(lod_fitness_filtered_full$mc_or_uni[lod_fitness_filtered_full$strategy_rep == i & lod_fitness_filtered_full$timepoint == j])) {
      for (l in unique(lod_fitness_filtered_full$reversion_mut_ID[lod_fitness_filtered_full$strategy_rep == i & lod_fitness_filtered_full$timepoint == j & lod_fitness_filtered_full$mc_or_uni == k])) {
        strategy_rep[length(strategy_rep)+1] <- i
        timepoint[length(timepoint)+1] <- j
        mc_or_uni[length(mc_or_uni)+1] <- k
        rev_mut_ID[length(rev_mut_ID)+1] <- l
        avg_growth_rate[length(avg_growth_rate)+1] <- mean(lod_fitness_filtered_full$growth_rate[lod_fitness_filtered_full$strategy_rep == i & lod_fitness_filtered_full$timepoint == j & lod_fitness_filtered_full$mc_or_uni == k & lod_fitness_filtered_full$reversion_mut_ID == l])
        growth_rate_sd[length(growth_rate_sd)+1] <- sd(lod_fitness_filtered_full$growth_rate[lod_fitness_filtered_full$strategy_rep == i & lod_fitness_filtered_full$timepoint == j & lod_fitness_filtered_full$mc_or_uni == k & lod_fitness_filtered_full$reversion_mut_ID == l])
        avg_workload[length(avg_workload)+1] <- mean(lod_fitness_filtered_full$workload[lod_fitness_filtered_full$strategy_rep == i & lod_fitness_filtered_full$timepoint == j & lod_fitness_filtered_full$mc_or_uni == k & lod_fitness_filtered_full$reversion_mut_ID == l])
        workload_sd[length(workload_sd)+1] <- sd(lod_fitness_filtered_full$workload[lod_fitness_filtered_full$strategy_rep == i & lod_fitness_filtered_full$timepoint == j & lod_fitness_filtered_full$mc_or_uni == k & lod_fitness_filtered_full$reversion_mut_ID == l])
      }
    }
  }
}

avg_lod_fitness_data <- cbind.data.frame(strategy_rep, timepoint, mc_or_uni, rev_mut_ID, avg_workload, workload_sd, avg_growth_rate, growth_rate_sd)
head(avg_lod_fitness_data)


# Calculate relative fitness using the averaged avg_time_to_fill values generated above (This should probably be done in the same for loop as above...)
rel_fitness <- c()

for (i in unique(avg_lod_fitness_data$strategy_rep)){
  for (j in unique(avg_lod_fitness_data$timepoint[avg_lod_fitness_data$strategy_rep == i])){
    for (k in avg_lod_fitness_data$avg_growth_rate[avg_lod_fitness_data$strategy_rep == i & avg_lod_fitness_data$timepoint == j]){
      rel_fitness[length(rel_fitness)+1] <- k / avg_lod_fitness_data$avg_growth_rate[avg_lod_fitness_data$strategy_rep == i & avg_lod_fitness_data$timepoint == j & avg_lod_fitness_data$mc_or_uni == "mc"]
    }
  }
}

length(rel_fitness) == length(avg_lod_fitness_data$strategy_rep)


# Merge with avg_lod_fitness_data
avg_lod_fitness_data <- cbind.data.frame(avg_lod_fitness_data, rel_fitness)

# Merge with combined_lod_fitness_summary_iterations (to get adjusted number of unicell revertants)
mc_or_uni <- rep("uni", length(combined_lod_fitness_summary_iterations$strategy_rep))
combined_lod_fitness_summary_iterations <- cbind.data.frame(combined_lod_fitness_summary_iterations, mc_or_uni)
colnames(combined_lod_fitness_summary_iterations)[3] <- "rev_mut_ID"
avg_lod_fitness_data <- left_join(avg_lod_fitness_data, combined_lod_fitness_summary_iterations, by=c("strategy_rep", "timepoint", "rev_mut_ID", "mc_or_uni"))

write.csv(avg_lod_fitness_data, "avg_lod_fitness_data_filtered_2_18JUN24.csv")


# START FROM HERE TO GENERATE PLOTS FROM SAVED DATA
setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/Figure_3/Figure_3_data")
theme_set(theme_cowplot(font_size = 7))

avg_lod_fitness_data <- read.csv("avg_lod_fitness_data_filtered_2_18JUN24.csv", header=T, stringsAsFactors=F)

avg_lod_fitness_data$timepoint <- gsub("0", "transition", avg_lod_fitness_data$timepoint)
avg_lod_fitness_data$timepoint <- gsub("1", "final", avg_lod_fitness_data$timepoint)
head(avg_lod_fitness_data)



# Plot standard deviation of growth rates and workloads as a way of estimating variation across replicates
P1 <- ggplot(data=avg_lod_fitness_data, aes(x=c(1:length(growth_rate_sd)), y=growth_rate_sd, col=mc_or_uni)) +
  geom_point(alpha=0.75) +
  xlab("Index") +
  ylab("Std. dev. (growth rate)") + 
  scale_color_manual(values = c("red", "blue"))
P2 <- ggplot(data=avg_lod_fitness_data, aes(x=c(1:length(workload_sd)), y=workload_sd, col=mc_or_uni)) +
  geom_point(alpha=0.75) + 
  xlab("Index") +
  ylab("Std. dev (workload)") + 
  scale_color_manual(values = c("red", "blue"))
ggdraw() +
  draw_plot(P1, x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(P2, x = 0, y = 0, width = 1, height = 0.5) 



#----------------------------------------
# PLOTTING THE DISTRIBUTION OF FITNESS EFFECTS OF REVERSION MUTATIONS

avg_lod_fitness_data$label <- factor(avg_lod_fitness_data$strategy_rep, levels = unique(avg_lod_fitness_data$strategy_rep), 
                                     labels = paste("DW.", order(unique(avg_lod_fitness_data$strategy_rep)), sep=""))


setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Figures")

# Relative fitness
pdf("Figure_S6_Relative_fitness_distributions_18JUN24.pdf", width=one_col*3, height=one_col*3.75)
ggplot(data=avg_lod_fitness_data[avg_lod_fitness_data$mc_or_uni != "mc", ], aes(x=rel_fitness, color=timepoint)) +
  #Use geom_freqpoly() for traces and geom_area() for filled distributions
  #geom_freqpoly(binwidth=0.025) +
  xlab("relative fitness") +
  #geom_area(aes(y = after_stat(count), fill=timepoint), stat = "bin", binwidth=0.025, position="identity", alpha=0.15) +
  geom_area(aes(y = after_stat(count), fill=timepoint, weight=proportion_fertile), stat = "bin", binwidth=0.025, position="identity", alpha=0.4) +
  facet_wrap(facets=avg_lod_fitness_data[avg_lod_fitness_data$mc_or_uni != "mc", ]$label, 
             scales = "free_y",
             ncol=7) +
  scale_fill_manual(values=c("green3", "darkgoldenrod1")) +
  scale_color_manual(values=c("green3", "darkgoldenrod1")) +
  labs(x="Relative fitness", y="Count", fill="Timepoint", color="Timepoint") +
  theme(legend.position=c(0.9, 0.025))
dev.off()

# Unicellular revertant growth rate
ggplot(data=avg_lod_fitness_data[avg_lod_fitness_data$mc_or_uni != "mc", ], aes(x=avg_growth_rate, color=timepoint)) +
  #xlim(0.08, 0.215) +
  #Use geom_freqpoly() for traces and geom_area() for filled distributions
  #geom_freqpoly(binwidth=0.005) +
  #geom_area(aes(y = after_stat(count), fill=timepoint), stat = "bin", bins = 50, position="identity", alpha=0.15) +
  geom_area(aes(y = after_stat(count), fill=timepoint, weight=proportion_fertile), stat = "bin", bins = 50, position="identity", alpha=0.4) +
  facet_wrap(facets=avg_lod_fitness_data[avg_lod_fitness_data$mc_or_uni != "mc", ]$label, 
             scales = "free_y",
             ncol=7) +
  scale_fill_manual(values=c("green3", "darkgoldenrod1")) +
  scale_color_manual(values=c("green3", "darkgoldenrod1")) +
  labs(x="Revertant growth rate", y="Count", fill="Timepoint", color="Timepoint") +
  theme(legend.position=c(0.9, 0.025))

# Multicellular growth rate
ggplot(data=avg_lod_fitness_data[avg_lod_fitness_data$mc_or_uni == "mc", ], aes(x=avg_growth_rate, color=timepoint)) +
  #xlim(0.08, 0.215) +
  #Use geom_freqpoly() for traces and geom_area() for filled distributions
  #geom_freqpoly(binwidth=0.005) +
  geom_area(aes(y = after_stat(count), fill=timepoint), stat = "bin", bins = 50, position="identity", alpha=0.4) +
  facet_wrap(facets=avg_lod_fitness_data[avg_lod_fitness_data$mc_or_uni == "mc", ]$label, 
             scales = "free_y",
             ncol=7) +
  scale_fill_manual(values=c("green3", "darkgoldenrod1")) +
  scale_color_manual(values=c("green3", "darkgoldenrod1")) +
  labs(x="Multicell growth rate", y="Count", fill="Timepoint", color="Timepoint") +
  theme(legend.position=c(0.9, 0.025))

# Workload
pdf("Figure_S7_Workload_distributions_18JUN24.pdf", width=one_col*3, height=one_col*3.75)
ggplot(data=avg_lod_fitness_data[avg_lod_fitness_data$mc_or_uni != "mc", ], aes(x=avg_workload, color=timepoint)) +
  #geom_freqpoly() +
  #geom_area(aes(y = after_stat(count), fill=timepoint), stat = "bin", binwidth=50, position="identity", alpha=0.15) +
  geom_area(aes(y = after_stat(count), fill=timepoint, weight=proportion_fertile), stat = "bin", binwidth=50, position="identity", alpha=0.4) +
  facet_wrap(facets=avg_lod_fitness_data[avg_lod_fitness_data$mc_or_uni != "mc", ]$label, 
             scales = "free_y",
             ncol=7) +
  scale_fill_manual(values=c("green3", "darkgoldenrod1")) +
  scale_color_manual(values=c("green3", "darkgoldenrod1")) +
  labs(x="Workload", y="Count", fill="Timepoint", color="Timepoint") +
  theme(legend.position=c(0.9, 0.025))
dev.off()


# ggplot(data=avg_lod_fitness_data, aes(x=rel_fitness, y=avg_growth_rate, color=timepoint)) +
#   #xlim(0.08, 0.215) +
#   #Use geom_freqpoly() for traces and geom_area() for filled distributions
#   #geom_freqpoly(binwidth=0.005) +
#   geom_point(aes(fill=timepoint), alpha=0.25) +
#   facet_wrap(facets=avg_lod_fitness_data$mc_or_uni, scales = "free_y") +
#   scale_fill_manual(values=c("green3", "darkgoldenrod1")) +
#   scale_color_manual(values=c("green3", "darkgoldenrod1"))


#Plot for individual isolate with highest entrenchment score
ggplot(data=avg_lod_fitness_data[avg_lod_fitness_data$mc_or_uni != "mc" & avg_lod_fitness_data$strategy_rep==3416,], aes(x=rel_fitness, color=timepoint)) +
  xlim(-0.1, 1.1) +
  xlab("Relative fitness") +
  ylab("Count") +
  #geom_freqpoly(binwidth=500) +
  geom_area(aes(y = ..count.., fill=timepoint), stat = "bin", binwidth=0.01, position="identity", alpha=0.15) +
  geom_area(aes(y = ..count.., fill=timepoint, weight=proportion_fertile), stat = "bin", binwidth=0.01, position="identity", alpha=0.7) +
  geom_rug(sides="b", alpha=0.5) +
  scale_fill_manual(values=c("green3", "darkgoldenrod1")) +
  scale_color_manual(values=c("green3", "darkgoldenrod1")) +
  ggtitle("3416") +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text = element_text(size=14), # Increase font size for axis labels
        axis.title = element_text(size=14))




#----------------------------------------
# SUMMARIZE DFE DATA AND ADD TO SUMMARY DATAFRAME

# Make new dataframe with avg_time_to_fill and relative_fitness
strategy_rep_temp <- c()
timepoint_temp <- c()
uni_avg_avg_growth_rate <- c()
multi_avg_growth_rate <- c()
uni_avg_rel_fitness <- c()
uni_avg_avg_workload <- c()
multi_avg_workload <- c()


for (i in unique(avg_lod_fitness_data$strategy_rep)){
  for (j in unique(avg_lod_fitness_data$timepoint[avg_lod_fitness_data$strategy_rep == i])){
    strategy_rep_temp[length(strategy_rep_temp)+1] <- i
    timepoint_temp[length(timepoint_temp)+1] <- j
    uni_avg_avg_growth_rate[length(uni_avg_avg_growth_rate)+1] <- mean(avg_lod_fitness_data$avg_growth_rate[avg_lod_fitness_data$strategy_rep == i & avg_lod_fitness_data$timepoint == j & avg_lod_fitness_data$mc_or_uni != "mc"])
    multi_avg_growth_rate[length(multi_avg_growth_rate)+1] <- mean(avg_lod_fitness_data$avg_growth_rate[avg_lod_fitness_data$strategy_rep == i & avg_lod_fitness_data$timepoint == j & avg_lod_fitness_data$mc_or_uni == "mc"])
    uni_avg_rel_fitness[length(uni_avg_rel_fitness)+1] <- mean(avg_lod_fitness_data$rel_fitness[avg_lod_fitness_data$strategy_rep == i & avg_lod_fitness_data$timepoint == j & avg_lod_fitness_data$mc_or_uni != "mc"])
    uni_avg_avg_workload[length(uni_avg_avg_workload)+1] <- mean(avg_lod_fitness_data$avg_workload[avg_lod_fitness_data$strategy_rep == i & avg_lod_fitness_data$timepoint == j & avg_lod_fitness_data$mc_or_uni != "mc"])
    multi_avg_workload[length(multi_avg_workload)+1] <- mean(avg_lod_fitness_data$avg_workload[avg_lod_fitness_data$strategy_rep == i & avg_lod_fitness_data$timepoint == j & avg_lod_fitness_data$mc_or_uni == "mc"])
    
  }
}

averaged_dataset <- cbind.data.frame(strategy_rep_temp, timepoint_temp, uni_avg_avg_workload, multi_avg_workload, uni_avg_avg_growth_rate, multi_avg_growth_rate, uni_avg_rel_fitness)
head(averaged_dataset)

strategy_rep_temp_unique <- averaged_dataset$strategy_rep_temp[c(TRUE, FALSE)]
trans_uni_avg_growth_rate <- averaged_dataset$uni_avg_avg_growth_rate[c(TRUE, FALSE)]
final_uni_avg_growth_rate <- averaged_dataset$uni_avg_avg_growth_rate[c(FALSE, TRUE)]
trans_uni_avg_rel_fitness <- averaged_dataset$uni_avg_rel_fitness[c(TRUE, FALSE)]
final_uni_avg_rel_fitness <- averaged_dataset$uni_avg_rel_fitness[c(FALSE, TRUE)]
trans_uni_avg_workload <- averaged_dataset$uni_avg_avg_workload[c(TRUE, FALSE)]
final_uni_avg_workload <- averaged_dataset$uni_avg_avg_workload[c(FALSE, TRUE)]

trans_multi_avg_growth_rate <- averaged_dataset$multi_avg_growth_rate[c(TRUE, FALSE)]
final_multi_avg_growth_rate <- averaged_dataset$multi_avg_growth_rate[c(FALSE, TRUE)]
trans_multi_avg_workload <- averaged_dataset$multi_avg_workload[c(TRUE, FALSE)]
final_multi_avg_workload <- averaged_dataset$multi_avg_workload[c(FALSE, TRUE)]

averaged_dataset_mod <- cbind.data.frame(strategy_rep_temp_unique, 
                                         trans_uni_avg_growth_rate, final_uni_avg_growth_rate, 
                                         trans_uni_avg_rel_fitness, final_uni_avg_rel_fitness, 
                                         trans_uni_avg_workload, final_uni_avg_workload,
                                         trans_multi_avg_growth_rate, final_multi_avg_growth_rate,
                                         trans_multi_avg_workload, final_multi_avg_workload)
head(averaged_dataset_mod)
colnames(averaged_dataset_mod)[1] <- "strategy_rep"


# Extract time as multi and difference in number of viable unicells (final-trans)
strategy_rep <- c()
time_as_mc <- c()
num_viable_unicells_diff <- c()
adjusted_num_viable_unicells_diff <- c()
trans_viable_uni <- c()
adjusted_trans_viable_uni <- c()
final_viable_uni <- c()
adjusted_final_viable_uni <- c()

head(lod_fitness_summary_F2_adjusted)

for (i in unique(lod_fitness_summary_F2_adjusted$strategy_rep)){
  strategy_rep[length(strategy_rep)+1] <- i
  #time_as_mc[length(time_as_mc)+1] <- lod_fitness_summary_F2_adjusted$update[lod_fitness_summary_F2_adjusted$strategy_rep==i & lod_fitness_summary_F2_adjusted$timepoint==1] - lod_fitness_summary_F2_adjusted$update[lod_fitness_summary_F2_adjusted$strategy_rep==i & lod_fitness_summary_F2_adjusted$timepoint==0]
  num_viable_unicells_diff[length(num_viable_unicells_diff)+1] <- lod_fitness_summary_F2_adjusted$num_viable_unicells_F2[lod_fitness_summary_F2_adjusted$strategy_rep==i & lod_fitness_summary_F2_adjusted$timepoint==1] - lod_fitness_summary_F2_adjusted$num_viable_unicells_F2[lod_fitness_summary_F2_adjusted$strategy_rep==i & lod_fitness_summary_F2_adjusted$timepoint==0]
  adjusted_num_viable_unicells_diff[length(adjusted_num_viable_unicells_diff)+1] <- lod_fitness_summary_F2_adjusted$adjusted_num_viable_uni[lod_fitness_summary_F2_adjusted$strategy_rep==i & lod_fitness_summary_F2_adjusted$timepoint==1] - lod_fitness_summary_F2_adjusted$adjusted_num_viable_uni[lod_fitness_summary_F2_adjusted$strategy_rep==i & lod_fitness_summary_F2_adjusted$timepoint==0]
  trans_viable_uni[length(trans_viable_uni)+1] <- lod_fitness_summary_F2_adjusted$num_viable_unicells_F2[lod_fitness_summary_F2_adjusted$strategy_rep==i & lod_fitness_summary_F2_adjusted$timepoint==0]
  final_viable_uni[length(final_viable_uni)+1] <- lod_fitness_summary_F2_adjusted$num_viable_unicells_F2[lod_fitness_summary_F2_adjusted$strategy_rep==i & lod_fitness_summary_F2_adjusted$timepoint==1]
  adjusted_trans_viable_uni[length(adjusted_trans_viable_uni)+1] <- lod_fitness_summary_F2_adjusted$adjusted_num_viable_uni[lod_fitness_summary_F2_adjusted$strategy_rep==i & lod_fitness_summary_F2_adjusted$timepoint==0]
  adjusted_final_viable_uni[length(adjusted_final_viable_uni)+1] <- lod_fitness_summary_F2_adjusted$adjusted_num_viable_uni[lod_fitness_summary_F2_adjusted$strategy_rep==i & lod_fitness_summary_F2_adjusted$timepoint==1]
}


# Save new values as dataframe
lod_fitness_summary_diff_data <- cbind.data.frame(strategy_rep, trans_viable_uni, final_viable_uni, num_viable_unicells_diff, adjusted_trans_viable_uni, adjusted_final_viable_uni, adjusted_num_viable_unicells_diff, rep("difference", length(strategy_rep)))

colnames(lod_fitness_summary_diff_data)[8] <- c("value")
head(lod_fitness_summary_diff_data)
head(averaged_dataset_mod)

full_entrenchment_data <- full_join(lod_fitness_summary_diff_data, averaged_dataset_mod, by="strategy_rep")
head(full_entrenchment_data)

# Save outputs
write.csv(avg_lod_fitness_data, "avg_lod_fitness_data_filtered_2_18JUN24.csv")
write.csv(full_entrenchment_data, "full_entrenchment_data_filtered_2_18JUN24.csv")


