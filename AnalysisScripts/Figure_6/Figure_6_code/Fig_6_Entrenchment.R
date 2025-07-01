# Figure 6, plots
# Peter Conlin
# 18 June 2024

library(ggplot2)
library(ggbeeswarm)
library(gridExtra)
library(GGally)
library(reshape)
library(dplyr)
library(Hmisc)
library(cowplot)
library(circlize)
theme_set(theme_cowplot(font_size = 7, font_family="Helvetica"))
library(viridis)
library(egg)

one_col <- 2.25
two_col <- 4.75

#   Lightest blue: #bbc6cd
#   Second lightest blue: #95aec2
#   Second darkest blue: #8aa8bf
#   Darkest blue: #608fb1
#   Lightest red: #cbb6b5
#   Second lightest red: #be8b86
#   Second darkest red: #b46e68
#   Darkest red: #a0352b


## Set parent working directory 
# If running in RStudio, you can use the following command to setwd() to source file location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')    #Go up a directory
pwd <- getwd()    #Set parent working directory
fwd <- c(strsplit(pwd, split="/"))     #Get a list of directories in file path

## Set working directory to load required data sets
setwd(paste(pwd, "/", fwd[[1]][length(fwd[[1]])], "_data", sep=""))
getwd()

# Read in main experimental data
avg_lod_fitness_data <- read.csv("avg_lod_fitness_data_filtered_2_18JUN24.csv", header=T, stringsAsFactors=F)
head(avg_lod_fitness_data)
full_entrenchment_data <- read.csv("full_entrenchment_data_filtered_2_18JUN24.csv", header=T, stringsAsFactors=F)
head(full_entrenchment_data)
lod_fitness_summary <- read.csv("lod_fitness_summary_F2_adjusted_18JUN24.csv", header=T, stringsAsFactors=F)

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

###################
# THIS IS WHERE I NEED TO CONTINUE FROM
###################

# Reformat DoL data
avg_dol_endpoints <- avg_dol_summary[avg_dol_summary$timepoint==1,]
differentiation_status <- cbind.data.frame(strategy_rep=avg_dol_endpoints$strategy_rep, has_soma=avg_dol_endpoints$has_soma)
head(differentiation_status)
# Merge reformatted dol data w/ summary dataframe
Entrenchment_DoL<- full_join(full_entrenchment_data, differentiation_status, by=c("strategy_rep"))
# Entrenchment_DoL <- full_join(Entrenchment_DoL_temp, dol_values, by=c("strategy_rep"))
head(Entrenchment_DoL)

# Read in entrenchment data
entrenchment_values <- read.csv("Entrenchment_values_from_logistic_regression_14FEB22.csv", header=T, stringsAsFactors=F)
# Reformat entrenchment data
logistic_entrench_trans <- entrenchment_values$inflection_point[entrenchment_values$timepoint==0]
logistic_entrench_final <- entrenchment_values$inflection_point[entrenchment_values$timepoint==1]
entrenchment_score_data <- cbind.data.frame(strategy_rep=entrenchment_values$strategy_rep[entrenchment_values$timepoint==0], logistic_entrench_trans, logistic_entrench_final, value=rep("difference", length(logistic_entrench_trans)))
# Merge reformatted entrenchment data w/ summary dataframe
Entrenchment_DoL_df <- full_join(Entrenchment_DoL, entrenchment_score_data, by=c("strategy_rep"))
head(Entrenchment_DoL_df)
length(Entrenchment_DoL_df)


# Read in distirbuted dirty work control data
DD_avg_lod_fitness_data <- read.csv("Dist_dirt_avg_lod_fitness_data_filtered_2_18JUN24.csv", header=T, stringsAsFactors=F)
DD_full_entrenchment_data <- read.csv("Dist_dirt_full_entrenchment_data_filtered_2_18JUN24.csv", header=T, stringsAsFactors=F)
DD_lod_fitness_summary <- read.csv("Dist_dirt_lod_fitness_summary_F2_adjusted_18JUN24.csv", header=T, stringsAsFactors=F)

# Create dummy division of labor data
DD_differentiation_status <- cbind.data.frame(strategy_rep=DD_full_entrenchment_data$strategy_rep, has_soma=rep("Unknown", length(DD_full_entrenchment_data$strategy_rep)))
# Merge dummy dol data w/ summary dataframe
DD_Entrenchment_DoL <- full_join(DD_full_entrenchment_data, DD_differentiation_status, by=c("strategy_rep"))
head(DD_Entrenchment_DoL)


# Read in entrenchment values
DD_entrenchment_values <- read.csv("Distributed_dirty_work_entrenchment_values_from_logistic_regression_15JULY20.csv", header=T, stringsAsFactors=F)
# Reformat entrenchment values
DD_logistic_entrench_trans <- DD_entrenchment_values$inflection_point[DD_entrenchment_values$timepoint==0]
DD_logistic_entrench_final <- DD_entrenchment_values$inflection_point[DD_entrenchment_values$timepoint==1]
DD_entrenchment_score_data <- cbind.data.frame(strategy_rep=DD_entrenchment_values$strategy_rep[DD_entrenchment_values$timepoint==0], logistic_entrench_trans=DD_logistic_entrench_trans, logistic_entrench_final=DD_logistic_entrench_final, value=rep("difference", length(DD_logistic_entrench_trans)))
# Merge entrenchment values w/ summary dataframe
DD_Entrenchment_DoL_df <- full_join(DD_Entrenchment_DoL, DD_entrenchment_score_data, by=c("strategy_rep"))
head(DD_Entrenchment_DoL_df)
length(DD_Entrenchment_DoL_df)
colnames(DD_Entrenchment_DoL_df)

# Remove unshared columns from Entrenchment_DoL_df
drops <- c("time_as_mc", "num_unicell_revertants", "num_inviable_unicells", "update")
Entrenchment_DoL_df_for_merge <-  Entrenchment_DoL_df[ , !(names(Entrenchment_DoL_df) %in% drops)]
length(Entrenchment_DoL_df_for_merge)
# Add column to specify which experiment data comes from
Entrenchment_DoL_df_for_merge_2 <- cbind.data.frame(Entrenchment_DoL_df_for_merge, experiment=rep("main", length(Entrenchment_DoL_df_for_merge$strategy_rep)))
DD_Entrenchment_DoL_df_for_merge <- cbind.data.frame(DD_Entrenchment_DoL_df, experiment=rep("control", length(DD_Entrenchment_DoL_df$strategy_rep)))

Full_Fig_6_data <- rbind(Entrenchment_DoL_df_for_merge_2, DD_Entrenchment_DoL_df_for_merge)
colnames(Full_Fig_6_data)

Full_Fig_6_data$has_soma <- factor(Full_Fig_6_data$has_soma, levels=c("Unknown", "0", "1"))



Full_Fig_6_data$has_soma <- factor(Full_Fig_6_data$has_soma, levels=c(1,0,"Unknown"))
head(Full_Fig_6_data)

write.csv(Full_Fig_6_data, "Full_Figure_6_data.csv")

x <- seq(from=2, to=100, by=4)
y1 <- rnorm(length(x), mean=10, sd=1.5)
y2 <- (y1/100) * x


df <- cbind.data.frame(x,y1,y2)

#White: #ffffff
#Lightest blue: #bbc6cd
#Second lightest blue: #95aec2
#Second darkest blue: #8aa8bf
#Darkest blue: #608fb1

scale_colour_gradientn(colours = c("#ffffff", "#bbc6cd", "#95aec2", "#8aa8bf", "#608fb1"))


P0B <- ggplot(data=df, aes(x=x, y=y2, fill=y2)) +
  geom_point(pch=21, stroke=0.25, size=1.5, color="black") +
  scale_y_continuous(limits=c(0, 10)) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position="none") +
  scale_fill_gradientn(colours = c("#ffffff", "#bbc6cd", "#95aec2", "#8aa8bf", "#608fb1")) +
  xlab("Workload per cell") +
  ylab("Mutations per cell")

P0D <- ggplot(data=df, aes(x=x, y=y1/2, fill=y2)) +
  geom_point(pch=21, stroke=0.25, size=1.5, color="black") +
  scale_y_continuous(limits=c(0, 10)) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position="none") +
  scale_fill_gradientn(colours = c("#ffffff", "#bbc6cd", "#95aec2", "#8aa8bf", "#608fb1")) +
  xlab("Workload per cell") +
  ylab("Mutations per cell")




# Read in detailed workload analyses
DW_workload_trans <- read.csv("Dirty_work_workload_detail_trans_summary_21AUG20.csv", header=T, stringsAsFactors=F)
DW_workload_final <- read.csv("Dirty_work_workload_detail_final_summary_15AUG20.csv", header=T, stringsAsFactors=F)
DD_workload_trans <- read.csv("Distributed_dirt_workload_detail_trans_summary_13AUG20.csv", header=T, stringsAsFactors=F)
DD_workload_final <- read.csv("Distributed_dirt_workload_detail_final_summary_15AUG20.csv", header=T, stringsAsFactors=F)


# Dirty work experiments
data <- DW_workload_trans
strategy_rep <- c()
avg_total_SD_trans <- c()
avg_total_CV_trans <- c()

for (i in unique(data$strategy_rep)) {
  strategy_rep[length(strategy_rep)+1] <- i
  avg_total_SD_trans[length(avg_total_SD_trans)+1] <- mean(data$total_workload_sd[data$strategy_rep==i], na.rm=TRUE)
  avg_total_CV_trans[length(avg_total_CV_trans)+1] <- mean(data$total_workload_sd[data$strategy_rep==i], na.rm=TRUE) / mean(data$total_workload_mean[data$strategy_rep==i], na.rm=TRUE)
}

DW_workload_trans_avg <- cbind.data.frame(strategy_rep, avg_total_SD_trans, avg_total_CV_trans)

data <- DW_workload_final
strategy_rep <- c()
avg_total_SD_final <- c()
avg_total_CV_final <- c()

for (i in unique(data$strategy_rep)) {
  strategy_rep[length(strategy_rep)+1] <- i
  avg_total_SD_final[length(avg_total_SD_final)+1] <- mean(data$total_workload_sd[data$strategy_rep==i], na.rm=TRUE)
  avg_total_CV_final[length(avg_total_CV_final)+1] <- mean(data$total_workload_sd[data$strategy_rep==i], na.rm=TRUE) / mean(data$total_workload_mean[data$strategy_rep==i], na.rm=TRUE)
}

DW_workload_final_avg <- cbind.data.frame(strategy_rep, avg_total_SD_final, avg_total_CV_final, treatment=rep("DW", length(strategy_rep)))

DW_workload_avg <- full_join(DW_workload_trans_avg, DW_workload_final_avg, by="strategy_rep")
DW_workload_avg_DOL <- full_join(DW_workload_avg, differentiation_status, by="strategy_rep")
head(DW_workload_avg_DOL)


# Distributed dirt experiments
data <- DD_workload_trans
strategy_rep <- c()
avg_total_SD_trans <- c()
avg_total_CV_trans <- c()

for (i in unique(data$strategy_rep)) {
  strategy_rep[length(strategy_rep)+1] <- i
  avg_total_SD_trans[length(avg_total_SD_trans)+1] <- mean(data$total_workload_sd[data$strategy_rep==i], na.rm=TRUE)
  avg_total_CV_trans[length(avg_total_CV_trans)+1] <- mean(data$total_workload_sd[data$strategy_rep==i], na.rm=TRUE) / mean(data$total_workload_mean[data$strategy_rep==i], na.rm=TRUE)
}
DD_workload_trans_avg <- cbind.data.frame(strategy_rep, avg_total_SD_trans, avg_total_CV_trans)

data <- DD_workload_final
strategy_rep <- c()
avg_total_SD_final <- c()
avg_total_CV_final <- c()

for (i in unique(data$strategy_rep)) {
  strategy_rep[length(strategy_rep)+1] <- i
  avg_total_SD_final[length(avg_total_SD_final)+1] <- mean(data$total_workload_sd[data$strategy_rep==i], na.rm=TRUE)
  avg_total_CV_final[length(avg_total_CV_final)+1] <- mean(data$total_workload_sd[data$strategy_rep==i], na.rm=TRUE) / mean(data$total_workload_mean[data$strategy_rep==i], na.rm=TRUE)
}
DD_workload_final_avg <- cbind.data.frame(strategy_rep, avg_total_SD_final, avg_total_CV_final, treatment=rep("DD", length(strategy_rep)))

DD_workload_avg <- full_join(DD_workload_trans_avg, DD_workload_final_avg, by="strategy_rep")
head(DD_workload_avg)

# Create dummy division of labor data
DD_differentiation_status <- cbind.data.frame(strategy_rep=DD_full_entrenchment_data$strategy_rep, has_soma=rep("Unknown", length(DD_full_entrenchment_data$strategy_rep)))
DD_workload_avg_DOL <- full_join(DD_workload_avg, DD_differentiation_status, by="strategy_rep")
head(DD_workload_avg_DOL)
avg_workload_variance <- rbind.data.frame(DW_workload_avg_DOL, DD_workload_avg_DOL)
avg_workload_variance$has_soma <- factor(avg_workload_variance$has_soma, levels=c(1,0,"Unknown"))
head(avg_workload_variance)

# # Save dirty work DOL data to Figure 5 folder
# Fig_5A_data <- Full_Fig_6_data[Full_Fig_6_data$experiment=="main",]
# Fig_5A_data <- left_join(Fig_5A_data, DW_workload_avg_DOL, by=c("strategy_rep"))
# unique(Fig_5A_data$has_soma.y)
# setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/Figure_5/Figure_5_data")
# write.csv(Fig_5A_data, "Figure_5A_DOL_data.csv")

# Standard deviation = DoL
P1 <- ggplot(avg_workload_variance, aes(x=factor(has_soma), y=avg_total_SD_final - avg_total_SD_trans)) +
  geom_violin(aes(fill=factor(has_soma)), scale="width", col=NA, alpha=0.4, trim=FALSE) +
  #stat_summary(geom = "crossbar", width = 1, fatten=0, color="white",
  #             fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  #geom_dotplot(binaxis='y', stackdir='center', fill="white", aes(col=factor(has_soma))) +
  #geom_boxplot(aes(fill=factor(has_soma)), col="black") +
  geom_quasirandom(aes(col=factor(has_soma)), method = "tukeyDense", size=0.1) +
  #geom_jitter(aes(col=factor(has_soma)), height=0, width=0.2, size=0.25) +
  scale_fill_manual(values=c(rgb(147/255, 17/255, 0), rgb(0, 84/255, 148/255), rgb(228/255, 155/255, 54/255))) +
  scale_color_manual(values=c(rgb(147/255, 17/255, 0), rgb(0, 84/255, 148/255), rgb(228/255, 155/255, 54/255))) +
  geom_hline(yintercept=0, lty=2) +
  scale_x_discrete(breaks=c("1","0", "Unknown"), labels=c("Differentiated", "Undifferentiated", "Distributed dirt")) +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank()) +
  xlab("") +
  ylab(expression(paste(Delta, " Division of labor")))


P2 <- ggplot(data=Full_Fig_6_data, 
             aes(x=factor(has_soma),
                 y=logistic_entrench_final - logistic_entrench_trans)) +
  geom_violin(aes(fill=factor(has_soma)), scale="width", col=NA, alpha=0.4, trim=FALSE) +
  #stat_summary(geom = "crossbar", width = 1, fatten=0, color="white",
  #             fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  #geom_dotplot(binaxis='y', stackdir='center', fill="white", aes(col=factor(has_soma))) +
  #geom_boxplot(aes(fill=factor(has_soma)), col="black") +
  geom_quasirandom(aes(col=factor(has_soma)), method = "tukeyDense", size=0.1) +
  #geom_jitter(aes(col=factor(has_soma)), height=0, width=0.2, size=0.25) +
  scale_fill_manual(values=c(rgb(147/255, 17/255, 0), rgb(0, 84/255, 148/255), rgb(228/255, 155/255, 54/255))) +
  scale_color_manual(values=c(rgb(147/255, 17/255, 0), rgb(0, 84/255, 148/255), rgb(228/255, 155/255, 54/255))) +
  scale_x_discrete(breaks=c("1","0", "Unknown"), labels=c("Differentiated", "Undifferentiated", "Distributed dirt")) +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank()) +
  geom_hline(yintercept=0.0, lty=2, col="black", lwd=0.5) +
  xlab("") + 
  ylab(expression(paste("Entrenchment")))

P3 <- ggplot(data=Full_Fig_6_data, 
             aes(x=factor(has_soma),
                 y=adjusted_final_viable_uni - adjusted_trans_viable_uni)) +
  geom_violin(aes(fill=factor(has_soma)), scale="width", col=NA, alpha=0.4, trim=FALSE) +
  #stat_summary(geom = "crossbar", width = 1, fatten=0, color="white",
  #             fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  #geom_dotplot(binaxis='y', stackdir='center', fill="white", aes(col=factor(has_soma))) +
  #geom_boxplot(aes(fill=factor(has_soma)), col="black") +
  geom_quasirandom(aes(col=factor(has_soma)), method = "tukeyDense", size=0.1) +
  #geom_jitter(aes(col=factor(has_soma)), height=0, width=0.2, size=0.25) +
  scale_fill_manual(values=c(rgb(147/255, 17/255, 0), rgb(0, 84/255, 148/255), rgb(228/255, 155/255, 54/255))) +
  scale_color_manual(values=c(rgb(147/255, 17/255, 0), rgb(0, 84/255, 148/255), rgb(228/255, 155/255, 54/255))) +
  scale_x_discrete(breaks=c("1","0", "Unknown"), labels=c("Differentiated", "Undifferentiated", "Distributed dirt")) +
  #theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank()) +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank()) +
  geom_hline(yintercept=0.0, lty=2, col="black", lwd=0.5) +
  xlab("") + 
  ylab(expression(paste(Delta, " # Revertants")))

P4 <- ggplot(data=Full_Fig_6_data, 
             aes(x=factor(has_soma),
                 y=final_uni_avg_rel_fitness - trans_uni_avg_rel_fitness)) +
  geom_violin(aes(fill=factor(has_soma)), scale="width", col=NA, alpha=0.4, trim=FALSE) +
  #stat_summary(geom = "crossbar", width = 1, fatten=0, color="white",
  #             fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  #geom_dotplot(binaxis='y', stackdir='center', fill="white", aes(col=factor(has_soma))) +
  #geom_boxplot(aes(fill=factor(has_soma)), col="black") +
  geom_quasirandom(aes(col=factor(has_soma)), method = "tukeyDense", size=0.1) +
  #geom_jitter(aes(col=factor(has_soma)), height=0, width=0.2, size=0.25) +
  scale_fill_manual(values=c(rgb(147/255, 17/255, 0), rgb(0, 84/255, 148/255), rgb(228/255, 155/255, 54/255))) +
  scale_color_manual(values=c(rgb(147/255, 17/255, 0), rgb(0, 84/255, 148/255), rgb(228/255, 155/255, 54/255))) +
  scale_x_discrete(breaks=c("1","0", "Unknown"), labels=c("Differentiated", "Undifferentiated", "Distributed dirt")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1), legend.position = "none") +
  geom_hline(yintercept=0.0, lty=2, col="black", lwd=0.5) +
  xlab("") + 
  ylab(expression(paste(Delta, " Relative fitness")))

P5 <- ggplot(data=Full_Fig_6_data, 
             aes(x=factor(has_soma),
                 y=final_multi_avg_growth_rate - trans_multi_avg_growth_rate)) +
  geom_violin(aes(fill=factor(has_soma)), scale="width", col=NA, alpha=0.4, trim=FALSE) +
  #stat_summary(geom = "crossbar", width = 1, fatten=0, color="white",
  #             fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  #geom_dotplot(binaxis='y', stackdir='center', fill="white", aes(col=factor(has_soma))) +
  #geom_boxplot(aes(fill=factor(has_soma)), col="black") +
  geom_quasirandom(aes(col=factor(has_soma)), method = "tukeyDense", size=0.1) +
  #geom_jitter(aes(col=factor(has_soma)), height=0, width=0.2, size=0.25) +
  scale_fill_manual(values=c(rgb(147/255, 17/255, 0), rgb(0, 84/255, 148/255), rgb(228/255, 155/255, 54/255))) +
  scale_color_manual(values=c(rgb(147/255, 17/255, 0), rgb(0, 84/255, 148/255), rgb(228/255, 155/255, 54/255))) +
  scale_x_discrete(breaks=c("1","0", "Unknown"), labels=c("Differentiated", "Undifferentiated", "Distributed dirt")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1), legend.position = "none") +
  geom_hline(yintercept=0.0, lty=2, col="black", lwd=0.5) +
  xlab("") + 
  ylab(expression(paste(Delta, " Multicell growth rate")))

P6 <- ggplot(data=Full_Fig_6_data, 
             aes(x=factor(has_soma),
                 y=final_uni_avg_growth_rate - trans_uni_avg_growth_rate)) +
  geom_violin(aes(fill=factor(has_soma)), scale="width", col=NA, alpha=0.4, trim=FALSE) +
  #stat_summary(geom = "crossbar", width = 1, fatten=0, color="white",
  #             fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +
  #geom_dotplot(binaxis='y', stackdir='center', fill="white", aes(col=factor(has_soma))) +
  #geom_boxplot(aes(fill=factor(has_soma)), col="black") +
  geom_quasirandom(aes(col=factor(has_soma)), method = "tukeyDense", size=0.1) +
  #geom_jitter(aes(col=factor(has_soma)), height=0, width=0.2, size=0.25) +
  scale_fill_manual(values=c(rgb(147/255, 17/255, 0), rgb(0, 84/255, 148/255), rgb(228/255, 155/255, 54/255))) +
  scale_color_manual(values=c(rgb(147/255, 17/255, 0), rgb(0, 84/255, 148/255), rgb(228/255, 155/255, 54/255))) +
  scale_x_discrete(breaks=c("1","0", "Unknown"), labels=c("Differentiated", "Undifferentiated", "Distributed dirt")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1), legend.position = "none") +
  geom_hline(yintercept=0.0, lty=2, col="black", lwd=0.5) +
  xlab("") + 
  ylab(expression(paste(Delta, " Unicell growth rate")))


P123456 <- ggarrange(P1, P2, P3, P4, P5, P6, nrow=2)
P1234 <- ggarrange(P1, P2, P3, P4, nrow=2)


# Change working directory
setwd(paste(pwd, "/", fwd[[1]][length(fwd[[1]])], "_drafts", sep=""))


pdf("Figure_6_Entrenchment.pdf", width=two_col, height=two_col)
ggdraw() +
  #draw_plot(P0A, x = 0.0, y = 0.75, width = 0.25, height = 0.25) +
  draw_plot(P0B, x = 0.25, y = 0.75, width = 0.245, height = 0.235) +
  #draw_plot(P0C, x = 0.50, y = 0.75, width = 0.25, height = 0.25) +
  draw_plot(P0D, x = 0.75, y = 0.75, width = 0.245, height = 0.235) +
  draw_plot(P123456, x=0.005, y=0.000, width=1, height=0.75) +
  #draw_plot(P1, x = 0.005, y = 0.41, width = 0.325, height = 0.35) +
  #draw_plot(P2, x = 0.335, y = 0.41, width = 0.325, height = 0.35) +
  #draw_plot(P3, x = 0.665, y = 0.41, width = 0.325, height = 0.35) +
  #draw_plot(P4, x = 0.005, y = 0.000, width = 0.325, height = 0.4) +
  #draw_plot(P5, x = 0.335, y = 0.000, width = 0.325, height = 0.4) +
  #draw_plot(P6, x = 0.665, y = 0.00, width = 0.325, height = 0.4) +
  draw_plot_label(label = c("A", " ","B", " ", "C", "D", "E", "F", "G", "H"), size = 9,
                  family = "Helvetica", fontface = "bold",
                  x = c(0.0, 0.25, 0.5, 0.75, 0.0, 0.33, 0.66,0.0, 0.33, 0.66),
                  y = c(1.0, 1.0, 1.0, 1.0, 0.75, 0.75, 0.75, 0.45, 0.45, 0.45))
dev.off()


