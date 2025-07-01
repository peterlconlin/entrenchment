# Figure 3, plots
# Peter Conlin
# 18 June 2024

#----------------------------------------
# PREPARE WORKSPACE
# Load packages
library(ggplot2)
library(ggbeeswarm)
library(gridExtra)
library(GGally)
library(reshape)
library(dplyr)
library(Hmisc)
library(cowplot)
theme_set(theme_cowplot(font_size = 7, font_family="Helvetica"))
library(viridis)
library(egg)

one_col <- 2.25
two_col <- 4.75

#Set working directory
# If running in RStudio, you can use the following command to setwd() to source file location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')    #Go up a directory
pwd <- getwd()    #Set parent working directory
fwd <- c(strsplit(pwd, split="/"))     #Get a list of directories in file path

## Set working directory to load required data sets
setwd(paste(pwd, "/", fwd[[1]][length(fwd[[1]])], "_data", sep=""))
getwd()


avg_lod_fitness_data <- read.csv("avg_lod_fitness_data_filtered_2.csv", header=TRUE, stringsAsFactors=FALSE)
head(avg_lod_fitness_data)

entrenchment_score_data <- read.csv("entrenchment.csv", header=T, stringsAsFactors=F)
logistic_entrenchment <- read.csv("Entrenchment_values_from_logistic_regression.csv", header=T, stringsAsFactors=F)
head(entrenchment_score_data)
head(logistic_entrenchment)
logistic_entrench_trans <- logistic_entrenchment$inflection_point[logistic_entrenchment$timepoint==0]
logistic_entrench_final <- logistic_entrenchment$inflection_point[logistic_entrenchment$timepoint==1]
entrenchment_score_data <- cbind.data.frame(entrenchment_score_data, logistic_entrench_trans, logistic_entrench_final, rep("difference", length(entrenchment_score_data$rep)))
colnames(entrenchment_score_data)[length(entrenchment_score_data)] <- "value_2"

full_entrenchment_data <- read.csv("full_entrenchment_data_filtered_2.csv", header=T, stringsAsFactors=F)

## Generate the schematic plots to illustrate possible shifts in the distribution of fitness effects of reversion mutations (DFERMs)
# Set parameters
anc.mean <- 1
anc.sd <- 0.05
evo.mean <- 0.25
evo.sd <- 0.05
x.min <- -0.25
x.max <- 1.25
anc.col <- 'darkgoldenrod1'
evo.col <- 'green3'

# Generate ancestral and evolved distributions
x.values <- c(x.min,seq(x.min,x.max,0.01),x.max) 
y.freq <- c(0,dnorm(seq(x.min,x.max,0.01), mean=anc.mean, sd=anc.sd),0) 
y.anc <- y.freq*2.5
y.both <- c(0,dnorm(seq(x.min,x.max,0.01), mean=evo.mean, sd=evo.sd),0)
y.effect <- y.both*2.5

####
DFEs <- cbind.data.frame(x.values, y.anc, y.freq, y.effect, y.both)


Fig_3A1 <- ggplot(data=DFEs, aes(x=x.values)) +
  geom_area(aes(y=y.anc), alpha=1, fill="darkgoldenrod1", lwd=1) +
  geom_area(aes(y=y.freq), alpha=1, fill="green3") +
  xlab("") +
  scale_x_continuous(limits=c(-0.25,1.25), breaks=c(0.5),
                     expand = expansion(mult = c(0, 0.05))) +
  ylab("Count") +
  scale_y_continuous(limits=c(0,20), breaks=c(0,10,20),
                     expand = expansion(mult = c(0, 0.05))) +
  annotate("text", x=0.35, y=15, label="transition", size=3, angle=0, hjust=0, color="darkgoldenrod1") +
  annotate("text", x=0.5, y=5, label="final", size=3, angle=0, hjust=0, color="green3") +
  theme(panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        #axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  panel_border(colour=NA, size=0.5)

Fig_3A2 <- ggplot(data=DFEs, aes(x=x.values)) +
  geom_area(aes(y=y.anc), alpha=1, fill="darkgoldenrod1") +
  geom_area(aes(y=y.effect), alpha=1, fill="green3") +
  xlab("Relative fitness (hypothetical)") +
  scale_x_continuous(limits=c(-0.25,1.25), breaks=c(0.5),
                     expand = expansion(mult = c(0, 0.05))) +
  ylab("") +
  scale_y_continuous(limits=c(0,20), breaks=c(0,10,20),
                     expand = expansion(mult = c(0, 0.05))) +
  theme(panel.border = element_blank(),
        #axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  panel_border(colour=NA, size=0.5)

Fig_3A3 <- ggplot(data=DFEs, aes(x=x.values)) +
  geom_area(aes(y=y.anc), alpha=1, fill="darkgoldenrod1") +
  geom_area(aes(y=y.both), alpha=1, fill="green3") +
  xlab("") +
  scale_x_continuous(limits=c(-0.25,1.25), breaks=c(0.5),
                     expand = expansion(mult = c(0, 0.05))) +
  ylab("") +
  scale_y_continuous(limits=c(0,20), breaks=c(0,10,20),
                     expand = expansion(mult = c(0, 0.05))) +
  theme(panel.border = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  panel_border(colour=NA, size=0.5)

Fig_3A <- ggarrange(Fig_3A1, Fig_3A2, Fig_3A3, nrow=1)


strategy_rep_labs <- c("Case study: 3416")
names(strategy_rep_labs) <- "3416"

Fig_3B <- ggplot(data=avg_lod_fitness_data[avg_lod_fitness_data$mc_or_uni != "mc" & avg_lod_fitness_data$strategy_rep==3416,], aes(x=rel_fitness)) +
  xlab("Relative fitness (case study)") +
  ylab("Count") +
  #stat_bin(aes(y=after_stat(count), color=timepoint),geom="step", binwidth=0.025) +
  #geom_density(aes(fill=timepoint, color=timepoint), alpha=0.5) +
  #geom_histogram(aes(fill=timepoint, color=timepoint), alpha=0.5, binwidth=0.01) +
  #geom_area(aes(y = after_stat(count), fill=timepoint, color=timepoint), stat = "bin", binwidth=0.025, position="identity", alpha=0.15) +
  geom_area(aes(y = after_stat(count), fill=timepoint, weight=proportion_fertile, color=timepoint), 
            stat = "bin", binwidth=0.025, position="identity", alpha=0.8, lwd=0) +
  geom_rug(aes(color=timepoint), sides="b", alpha=0.3, lwd=0.5, length=unit(0.075, "npc")) +
  #facet_grid(.~ strategy_rep) + #taking the place of title
  scale_x_continuous(breaks=c(0.0, 0.5), limits = c(0, 0.576),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expansion(0.125)) +
  scale_fill_manual(values=c("green3", "darkgoldenrod1"), guide = "none") +
  scale_color_manual(values=c("green3", "darkgoldenrod1"), guide = "none") +
  labs(fill = "Timepoint") +
  #ggtitle("4278") +
  annotate("text", x=0.475, y=30, label="transition", size=3, angle=0, hjust=0, color="darkgoldenrod1") +
  annotate("text", x=0.05, y=30, label="final", size=3, angle=0, hjust=0, color="green3") +
  theme(legend.position = "bottom", 
        legend.justification = c(0.5,0.5),
        panel.border = element_blank())


# Make plot of difference in entrenchment cost

Fig_3C <- ggplot(data=entrenchment_score_data, aes(x=value_2, y=logistic_entrench_final - logistic_entrench_trans)) +
  geom_violin(fill="grey90", scale="count", col=NA, trim=FALSE) +
  geom_quasirandom(data=entrenchment_score_data[entrenchment_score_data$rep!=3416,], method = "tukeyDense", size=0.1, color="grey40") +
  geom_point(data=entrenchment_score_data[entrenchment_score_data$rep==3416,], color="#FF6900", size=2.5, pch=18) +
  # geom_jitter(width=ifelse(entrenchment_score_data$rep==3416, 0, 0.35), 
  #             height=ifelse(entrenchment_score_data$rep==3416, 0, 0.1), 
  #             color=ifelse(entrenchment_score_data$rep==3416,"#FF6900", "grey40"), 
  #             size=ifelse(entrenchment_score_data$rep==3416, 1.5, 0.1),
  #             pch=ifelse(entrenchment_score_data$rep==3416, 18, 19)) +
  geom_hline(yintercept=0, linetype=2, color="black", size=0.25) +
  xlab("") +
  #scale_y_continuous(breaks=c(0,2,4,6,8)) +
  ylab(expression(paste(Delta, " Stability (Entrenchment)"))) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.border = element_blank())


# Make plot of difference in number of reversion mutations
Fig_3D <- ggplot(full_entrenchment_data, aes(x=value, y=adjusted_num_viable_unicells_diff)) +
  geom_violin(fill="grey90", scale="count", col=NA, trim=FALSE) +
  geom_quasirandom(data=full_entrenchment_data[full_entrenchment_data$strategy_rep!=3416,], method = "tukeyDense", size=0.1, color="grey40") +
  geom_point(data=full_entrenchment_data[full_entrenchment_data$strategy_rep==3416,], color="#FF6900", size=2.5, pch=18) +
  geom_hline(yintercept=1, linetype=2, color="black", size=0.25) +
  #geom_dotplot(binaxis='y', stackdir='center', binwidth = 50, fill=rgb(0,0,0,0.5), color=NA) +
  xlab("") +
  ylab(expression(paste(Delta, " # Revertants"))) + 
  #ylim(-max(abs(full_entrenchment_data$num_viable_unicells_diff)), max(abs(full_entrenchment_data$num_viable_unicells_diff))) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.border = element_blank())


# Make plot of difference in average relative fitness effect of reversion mutations
Fig_3E <- ggplot(full_entrenchment_data, aes(x=value, y=final_uni_avg_rel_fitness - trans_uni_avg_rel_fitness)) +
  geom_violin(fill="grey90", scale="count", col=NA, trim=FALSE) + 
  geom_quasirandom(data=full_entrenchment_data[full_entrenchment_data$strategy_rep!=3416,], method = "tukeyDense", size=0.1, color="grey40") +
  geom_point(data=full_entrenchment_data[full_entrenchment_data$strategy_rep==3416,], color="#FF6900", size=2.5, pch=18) +
  geom_hline(yintercept=0, linetype=2, color="black", size=0.25) +
  #geom_dotplot(binaxis='y', stackdir='center', binwidth=0.035, fill=rgb(0,0,0,0.5), color=NA) + 
  xlab("") +
  ylab(expression(paste(Delta, " Relative fitness"))) +  #ylim(-1, 1) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.border = element_blank())





setwd(paste(pwd, "/", fwd[[1]][length(fwd[[1]])], "_drafts", sep=""))

# Figure 3
pdf("Figure_3_Entrenchment_FINAL.pdf", width=two_col, height=two_col*1.2)
ggdraw() +
  draw_plot(Fig_3A, x=0.0475, y=0.75, width=0.98, height=0.25) +
  #draw_plot(Fig_3A1, x = 0.035, y = 0.75, width = 0.33, height = 0.25) +
  #draw_plot(Fig_3A2, x = 0.35, y = 0.75, width = 0.33, height = 0.25) +
  #draw_plot(Fig_3A3, x = 0.66, y = 0.75, width = 0.33, height = 0.25) +
  draw_plot(Fig_3B, x = 0.025, y = 0.42, width = 0.955, height = 0.35) +
  draw_plot(Fig_3C, x = 0.01, y = 0.01, width = 0.33, height = 0.40) +
  draw_plot(Fig_3D, x = 0.335, y = 0.01, width = 0.33, height = 0.40) +
  draw_plot(Fig_3E, x = 0.66, y = 0.01, width = 0.33, height = 0.40) +
  draw_plot_label(label = c("A", " ", " ", "B", "C", "D", "E"), size = 9,
                  family = "Helvetica", fontface = "bold",
                  x = c(0, 0.33, 0.66, 0, 0, 0.33, 0.66),
                  y = c(1, 1, 1, 0.76, 0.425, 0.425, 0.425))
dev.off()


#------------------------------------------------
# STATISTICS

# Wilcoxon signed rank test
# Number of viable revertants @ transition vs. final time point
x <- full_entrenchment_data$adjusted_trans_viable_uni
y <- full_entrenchment_data$adjusted_final_viable_uni

test <- wilcox.test(x=x, y=y, paired=TRUE)
Zstat<-qnorm(test$p.value/2)
test                          #print p-value
abs(Zstat)/sqrt(length(x))    #print effect size


# Wilcoxon signed rank test
# Relative fitness @ transition vs. final time point
x <- full_entrenchment_data$trans_uni_avg_rel_fitness
y <- full_entrenchment_data$final_uni_avg_rel_fitness

test <- wilcox.test(x=x, y=y, paired=TRUE)
Zstat<-qnorm(test$p.value/2)
test                          #print test results
abs(Zstat)/sqrt(length(x))    #print effect size
