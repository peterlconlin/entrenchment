# Figure 4 re-plotting scripts
# 18 June 2024
# Peter Conlin

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
theme_set(theme_cowplot(font_size = 7, font_family="Helvetica"))
library(viridis)
library(tidyr)

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

# If running in RStudio, you can use the following command to setwd() to source file location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')    #Go up a directory
pwd <- getwd()    #Set parent working directory
fwd <- c(strsplit(pwd, split="/"))     #Get a list of directories in file path

#----------------------------------------
## Set working directory to load required data sets
setwd(paste(pwd, "/", fwd[[1]][length(fwd[[1]])], "_data", sep=""))
getwd()

size_and_workload <- read.csv("3416_over_time_size_workload.csv", header=T, stringsAsFactors=F)
head(size_and_workload)
entrenchment_scores <- read.csv("3416_case_study_over_time_entrenchment_values_from_logistic_regression.csv", header=T, stringsAsFactors=F)
head(entrenchment_scores)
entrenchment_scores_trans_and_final <- read.csv("Entrenchment_values_from_logistic_regression.csv", header=T, stringsAsFactors=F)
entrenchment_scores_trans_and_final <- entrenchment_scores_trans_and_final[entrenchment_scores_trans_and_final$strategy_rep==3416,]
entrenchment_scores_trans_and_final$strategy_rep <- c(35, 6026)
#entrenchment_scores_trans_and_final <- entrenchment_scores_trans_and_final[,-c(1,4)]
head(entrenchment_scores_trans_and_final)

entrenchment_scores <- rbind.data.frame(entrenchment_scores, entrenchment_scores_trans_and_final)
head(entrenchment_scores)

entrenchment_df <- cbind.data.frame(timepoint=unique(entrenchment_scores$strategy_rep), 
                                     entrenchment=entrenchment_scores$inflection_point -
                                       entrenchment_scores$inflection_point[entrenchment_scores$strategy_rep==35])
head(entrenchment_df)


lod_fitness <- read.csv("Case_study_over_time_full_summary_data_filtered_2.csv", header=T, stringsAsFactors=F)
head(lod_fitness)


head(size_and_workload)
temp_1 <- cbind.data.frame(timepoint=size_and_workload$timepoint, count=size_and_workload$germ_count, cell_type=rep("germ", length(size_and_workload$timepoint)))
temp_2 <- cbind.data.frame(timepoint=size_and_workload$timepoint, count=size_and_workload$multicell_size - size_and_workload$germ_count, cell_type=rep("soma", length(size_and_workload$timepoint)))
cell_type_df <- rbind.data.frame(temp_1, temp_2)
cell_type_df$cell_type <- factor(cell_type_df$cell_type, levels=c("soma", "germ"))

head(size_and_workload)

Fig_4A <- ggplot(data=size_and_workload, aes(x=timepoint, y=multicell_size)) +
  geom_area(data=size_and_workload, aes(x=timepoint, y=multicell_size, fill="#608fb1"), lwd=0.00000000001) +
  geom_area(data=size_and_workload, aes(x=timepoint, y=germ_count, fill="#a0352b"), lwd=0.00000000001) +
  geom_line(size=0.1) +
  #scale_x_continuous(limits=c(6,6026), trans="log10") +
  scale_y_continuous() +
  scale_fill_manual(values=c("#a0352b", "#608fb1")) +
  xlab(" ") +
  ylab("Adult size") +
  theme(legend.position = "none", plot.margin = unit(c(7, 1, 0, 4), "pt"), axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

temp_1 <- cbind.data.frame(timepoint=size_and_workload$timepoint, fraction=size_and_workload$germ_workload, cell_type=rep("germ", length(size_and_workload$timepoint)))
temp_2 <- cbind.data.frame(timepoint=size_and_workload$timepoint, fraction=size_and_workload$soma_workload, cell_type=rep("soma", length(size_and_workload$timepoint)))
cell_type_df <- rbind.data.frame(temp_1, temp_2)
cell_type_df$cell_type <- factor(cell_type_df$cell_type, levels=c("soma", "germ"))

Fig_4B <- ggplot(data=size_and_workload, aes(x=timepoint, y=total_workload)) +
  geom_area(data=cell_type_df, aes(x=timepoint, y=fraction, fill=cell_type), lwd=0.00000000001) +
  geom_line(size=0.1) +
  scale_x_continuous(limits=c(0,6026)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  
  scale_fill_manual(values=c("#a0352b", "#608fb1")) +
  xlab(" ") +
  ylab("Workload") +
  theme(legend.position = "none", plot.margin = unit(c(0, 1, 0, 4), "pt"), axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

Fig_4C <- ggplot(entrenchment_df, aes(x=timepoint, y=entrenchment)) +
  geom_line(size=0.25) +
  geom_point(size=0.33) +
  xlab(" ") +
  ylab("Entrenchment") +
  scale_x_continuous(limits=c(0,6026)) +
  theme(plot.margin = unit(c(0, 1, 0, 4), "pt"), axis.ticks.x = element_blank(),
        axis.text.x = element_blank())


head(lod_fitness)
lod_fitness <- lod_fitness[lod_fitness$timepoint!=6000,]
Fig_4D <- ggplot(lod_fitness, aes(x=timepoint, y=adjusted_num_viable_uni)) +
  geom_line(size=0.25) +
  geom_point(size=0.33) +
  scale_x_continuous(limits=c(0,6026)) +
  scale_y_continuous(limits=c(0,max(lod_fitness$adjusted_num_viable_uni)),
                     expand = expansion(mult = c(0, 0.05))) +
  xlab(" ") +
  ylab("# Revertants") +
  theme(plot.margin = unit(c(0, 1, 0, 4), "pt"), axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

Fig_4E <- ggplot(lod_fitness, aes(x=timepoint, y=uni_avg_rel_fitness)) +
  geom_line(size=0.25) +
  geom_point(size=0.33) +
  scale_x_continuous(limits=c(0,6026)) +
  scale_y_continuous(limits=c(0,0.4), breaks=c(0,0.1,0.2,0.3,0.4)) +
  xlab("Time (generations)") +
  ylab("Rel. fitness") +
  theme(plot.margin = unit(c(0, 1, 7, 4), "pt"))


setwd(paste(pwd, "/", fwd[[1]][length(fwd[[1]])], "_drafts", sep=""))

pdf("Figure_4_Entrenchment_FINAL.pdf", width=one_col, height=one_col*1.66)
plot_grid(Fig_4A, Fig_4B, Fig_4C, Fig_4D, Fig_4E, align="v", ncol=1, rel_heights=c(1,1,1,1,1.275))
dev.off()





