library(ggplot2)
library(gridExtra)
library(GGally)
library(reshape)
library(dplyr)
library(Hmisc)
library(cowplot)
theme_set(theme_cowplot(font_size = 14))
library(viridis)
library(ggpubr)


# #Remake Figure S1
# S1_data <- cbind.data.frame(efficiency_level = c("5%", "6%", "7%", "8%"), percent_multi = c(6.5, 7.8, 11.0, 12.8))
# ggplot(S1_data, aes(x=efficiency_level, y=percent_multi)) +
#   geom_bar(stat="identity") +
#   ylab("Percent multicellular") +
#   xlab("Efficiency level")

setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/")
timepoint_data <- read.csv("Time_of_transition_and_final_timepoint_in_generations_10NOV22.csv", header=T, stringsAsFactors=F)
timepoint_data <- timepoint_data[,-1]
head(timepoint_data)

setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/")

# summary_data <- read.csv("lod_fitness_summary.csv", header=T, stringsAsFactors=F)
# summary_data$timepoint[summary_data$timepoint==60] <- 0
# summary_data <- summary_data[,-1]
# head(summary_data)
# time_as_mc_df <- left_join(summary_data, timepoint_data, by=c("strategy_rep", "timepoint"))
# head(time_as_mc_df)
  
strategy_rep <- c()
time_as_mc <- c()

for (i in unique(timepoint_data$strategy_rep)) {
  strategy_rep[length(strategy_rep)+1] <- i
  time_as_mc[length(time_as_mc)+1] <- timepoint_data$generation[timepoint_data$strategy_rep==i & timepoint_data$timepoint==1] - timepoint_data$generation[timepoint_data$strategy_rep==i & timepoint_data$timepoint==0]
}

time_as_mc_df <- cbind.data.frame(strategy_rep, time_as_mc)

setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/Figure_6/Figure_6_data/")
# entrenchment_data <- read.csv("full_entrenchment_data_filtered_2_14FEB22.csv")
# head(entrenchment_data)
Fig_6_data <- read.csv("Full_Figure_6_data_25OCT22.csv", header=T, stringsAsFactors=F)
Fig_6_data <- Fig_6_data[,-c(1,2)]
head(Fig_6_data)

full_df <- left_join(Fig_6_data, time_as_mc_df, by=c("strategy_rep"))
head(full_df)


# Entrenchment
P1 <- ggplot(data=full_df[full_df$experiment=="main",], aes(x=time_as_mc, y=logistic_entrench_final - logistic_entrench_trans)) +
  geom_point(aes(color=has_soma)) +
  #geom_smooth(formula=y~x, method="lm", aes(color=has_soma)) +
  geom_smooth(formula=y~x, method="lm", color="black") +
  #stat_cor(label.x = 3500, label.y = 385, method="spearman", cor.coef.name="rho") +
  geom_hline(yintercept=0, lty=3, color="grey80") +
  # stat_regline_equation( aes(label =  paste(..eq.label.., ..adj.rr.label.., 
  #                                           sep = "~~~~")),
  #                        label.x.npc = "center", label.y.npc = "top", size = 3) +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_color_manual(values=c(rgb(147/255, 17/255, 0), rgb(0, 84/255, 148/255), rgb(228/255, 155/255, 54/255)),
                     breaks=c("1","0", "Unknown"),
                     labels=c("Differentiated", "Non-differentiated", "Distributed dirt")) +  
  xlab(" ") +
  ylab(expression(paste("Entrenchment"))) +
  theme(legend.position="none", plot.margin = margin(0.5,0.4,0,0.3, "cm"),
        axis.text.x = element_text(angle = 20, vjust = 1, hjust=1)) 

# Number of revertants
P2 <- ggplot(data=full_df[full_df$experiment=="main",], aes(x=time_as_mc, y=adjusted_num_viable_unicells_diff)) +
  geom_point(aes(color=has_soma)) +
  #geom_smooth(formula=y~x, method="lm", aes(color=has_soma)) +
  geom_smooth(formula=y~x, method="lm", color="black") +
  #stat_cor(label.x = 3500, label.y = 125, method="spearman", cor.coef.name="rho") +
  geom_hline(yintercept=0, lty=3, color="grey80") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_color_manual(values=c(rgb(147/255, 17/255, 0), rgb(0, 84/255, 148/255), rgb(228/255, 155/255, 54/255)),
                     breaks=c("1","0", "Unknown"),
                     labels=c("Differentiated", "Non-differentiated", "Distributed dirt")) +
  xlab(" ") +
  ylab(expression(paste(Delta, " # Revertants"))) +
  theme(legend.position="none", plot.margin = margin(0.5,0.4,0,0.3, "cm"),
        axis.text.x = element_text(angle = 20, vjust = 1, hjust=1)) 

# Relative fitness
P3 <- ggplot(data=full_df[full_df$experiment=="main",], aes(x=time_as_mc, y=final_uni_avg_rel_fitness - trans_uni_avg_rel_fitness)) +
  geom_point(aes(color=has_soma)) +
  #geom_smooth(formula=y~x, method="lm", aes(color=has_soma)) +
  geom_smooth(formula=y~x, method="lm", color="black") +
  #stat_cor(label.x = 3500, label.y = 0.282, method="spearman", cor.coef.name="rho") +
  geom_hline(yintercept=0, lty=3, color="grey80") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_color_manual(values=c(rgb(147/255, 17/255, 0), rgb(0, 84/255, 148/255), rgb(228/255, 155/255, 54/255)),
                     breaks=c("1","0", "Unknown"),
                     labels=c("Differentiated", "Non-differentiated", "Distributed dirt")) +
  xlab("Time as multicell (generations)") +
  ylab(expression(paste(Delta, " Revertant relative fitness")))  +
  theme(legend.position="none", plot.margin = margin(0.5,0.4,0,0.3, "cm"),
        axis.text.x = element_text(angle = 20, vjust = 1, hjust=1)) 

# Unicell growth rate
P4 <- ggplot(data=full_df[full_df$experiment=="main",], aes(x=time_as_mc, y=final_uni_avg_growth_rate - trans_uni_avg_growth_rate)) +
  geom_point(aes(color=has_soma)) +
  #geom_smooth(formula=y~x, method="lm", aes(color=has_soma)) +
  geom_smooth(formula=y~x, method="lm", color="black") +
  #stat_cor(label.x = 3500, label.y = 0.001175, method="spearman", cor.coef.name="rho") +
  geom_hline(yintercept=0, lty=3, color="grey80") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_color_manual(values=c(rgb(147/255, 17/255, 0), rgb(0, 84/255, 148/255), rgb(228/255, 155/255, 54/255)),
                     breaks=c("1","0", "Unknown"),
                     labels=c("Differentiated", "Non-differentiated", "Distributed dirt")) +
  xlab("Time as multicell (generations)") +
  ylab(expression(paste(Delta, " Revertant growth rate"))) +
  theme(legend.position="none", plot.margin = margin(0.5,0.4,0,0.3, "cm"),
        axis.text.x = element_text(angle = 20, vjust = 1, hjust=1)) 

# Multicell growth rate
P5 <- ggplot(data=full_df[full_df$experiment=="main",], aes(x=time_as_mc, y=final_multi_avg_growth_rate - trans_multi_avg_growth_rate)) +
  geom_point(aes(color=has_soma)) +
  #geom_smooth(formula=y~x, method="lm", aes(color=has_soma)) +
  geom_smooth(formula=y~x, method="lm", color="black") +
  #stat_cor(label.x = 3500, label.y = 0.0275, method="spearman", cor.coef.name="rho") +
  geom_hline(yintercept=0, lty=3, color="grey80") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_color_manual(values=c(rgb(147/255, 17/255, 0), rgb(0, 84/255, 148/255), rgb(228/255, 155/255, 54/255)),
                     breaks=c("1","0", "Unknown"),
                     labels=c("Differentiated", "Non-differentiated", "Distributed dirt")) +
  xlab("Time as multicell (generations)") +
  ylab(expression(paste(Delta, " Multicell growth rate"))) +
  theme(legend.position="none", plot.margin = margin(0.5,0.4,0,0.3, "cm"),
        axis.text.x = element_text(angle = 20, vjust = 1, hjust=1)) 

P6 <- ggplot(data=full_df[full_df$experiment=="main",], aes(x=time_as_mc, y=final_multi_avg_growth_rate - trans_multi_avg_growth_rate)) +
  geom_point(aes(color=has_soma)) +
  #geom_smooth(formula=y~x, method="lm", aes(color=has_soma)) +
  geom_smooth(formula=y~x, method="lm", color="black") +
  #stat_cor(label.x = 3500, label.y = 0.0275, method="spearman", cor.coef.name="rho") +
  geom_hline(yintercept=0, lty=3, color="grey80") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_color_manual(values=c(rgb(147/255, 17/255, 0), rgb(0, 84/255, 148/255), rgb(228/255, 155/255, 54/255)),
                     breaks=c("1","0", "Unknown"),
                     labels=c("Differentiated", "Non-differentiated", "Distributed dirt")) +
  xlab("Time as multicell (generations)") +
  ylab(expression(paste(Delta, " Multicell growth rate"))) +
  labs(color="Cell differentiation status")

leg <- get_legend(P6)

plot_grid(P1, P2, P3, P4, P5, leg, align="hv", ncol=3, axis="l", labels=c("A", "B", "C", "D", "E"))
#ggarrange(P1, P2, P3, P4, P5, leg, labels=c("A", "B", "C", "D", "E"), align="hv", axis="l")
