# Figure 5, plots
# Peter Conlin
# 11 JULY 20


library(ggplot2)
library(cowplot)
theme_set(theme_cowplot(font_size = 7, font_family="Helvetica"))
library(dplyr)
library(directlabels)
library(reshape)
library(viridis)
library(ggrepel)
library(ComplexHeatmap)
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

# If running in RStudio, you can use the following command to setwd() to source file location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')    #Go up a directory
pwd <- getwd()    #Set parent working directory
fwd <- c(strsplit(pwd, split="/"))     #Get a list of directories in file path

## Set working directory to load required data sets
setwd(paste(pwd, "/", fwd[[1]][length(fwd[[1]])], "_data", sep=""))
getwd()


### FIGURE 5A ###
Fig_5A_data <- read.csv("Figure_5A_DOL_data.csv", header=T, stringsAsFactors=F)

head(Fig_5A_data)

unique(Fig_5A_data$has_soma.y)

head(Fig_5A_data)

ggplot(data=Fig_5A_data, aes(x=avg_total_SD_final - avg_total_SD_trans, y=final_uni_avg_rel_fitness - trans_uni_avg_rel_fitness, color=as.factor(has_soma.x))) +
  geom_point(size=2, alpha=0.8) +
  scale_color_manual(values=c("#608fb1", "#a0352b"), labels=c("Undiff.", "Diff."), name = "Cell diff. status") +
  xlab(expression(paste(Delta, " Division of labor"))) +
  ylab(expression(paste(Delta, " Relative fitness"))) +
  theme(legend.position="top", 
        legend.justification = "left",
        legend.title = element_text(size=14),
        legend.spacing.x = unit(0.01, 'cm'),
        #legend.background = element_rect(fill=NA,
        #                                 size=0.1, linetype="solid", 
        #                                 colour ="black"),
        legend.text=element_text(size=14))



Fig_5A <- ggplot(data=Fig_5A_data, aes(x=avg_total_SD_final - avg_total_SD_trans, y=logistic_entrench_final - logistic_entrench_trans, color=as.factor(has_soma.x))) +
  geom_point(size=0.5, alpha=0.8) +
  scale_color_manual(values=c("#608fb1", "#a0352b"), labels=c("Undifferentiated", "Differentiated"), name = "Differentiation status:") +
  xlab(expression(paste(Delta, " Division of labor"))) +
  ylab(expression(paste("Enternchment"))) +
  guides(color = guide_legend(nrow = 2)) +
  theme(legend.position="top",
        legend.justification = "left",
        legend.direction = "vertical",
        legend.title = element_text(size=6),
        legend.spacing.x = unit(0.001, 'cm'),
        #legend.background = element_rect(fill=NA,
        #                                 size=0.1, linetype="solid", 
        #                                 colour ="black"),
        legend.text=element_text(size=6))


t.test(x=Fig_5A_data$logistic_entrench_final[Fig_5A_data$has_soma.x==1] - Fig_5A_data$logistic_entrench_trans[Fig_5A_data$has_soma.x==1], 
       y=Fig_5A_data$logistic_entrench_final[Fig_5A_data$has_soma.x==0] - Fig_5A_data$logistic_entrench_trans[Fig_5A_data$has_soma.x==0],
       alternative="two.sided")

cor.test(x=Fig_5A_data$avg_total_SD_final - Fig_5A_data$avg_total_SD_trans, 
         y=Fig_5A_data$logistic_entrench_final - Fig_5A_data$logistic_entrench_trans,
         alternative="two.sided",
         method="spearman")


### FIGURES 5D - 5E ###
# Generate prediction for tradeoff between workload and fidelity
genome_length <- 100
work <- c(1:10000)
mutations <- (37/38)*(0.01 + work * 0.00075) * 100
corrected_mutations <- 1 - exp(-mutations/genome_length)
fidelity <- 1 - corrected_mutations
genome_length_100 <- cbind.data.frame(work, fidelity)
head(genome_length_100)

genome_length <- 50
corrected_mutations <- 1 - exp(-mutations/genome_length)
fidelity <- 1 - corrected_mutations
genome_length_50 <- cbind.data.frame(work, fidelity)


### FIGURE 5D ###
mutation_summary_complete_5D <- read.csv("mutation_summary_complete_multicells_at_final_timepoint_filtered.csv", header=T, stringsAsFactors=F)
mutation_summary_complete_5D <- mutation_summary_complete_5D[,-1]
head(mutation_summary_complete_5D)
mutation_summary_complete_5D <- mutation_summary_complete_5D[complete.cases(mutation_summary_complete_5D),]

# Average data by "uni_id"
head(mutation_summary_complete_5D)

uni_id <- c()
timepoint <- c()
num_reps <- c()
counts <- c()
avg_workload <- c()
avg_growth_rate <- c()
avg_OGL <- c()
per_site_mu <- c()

for (i in unique(mutation_summary_complete_5D$uni_id)) {
  uni_id[length(uni_id)+1] <- i
  timepoint[length(timepoint)+1] <- unique(mutation_summary_complete_5D$timepoint[mutation_summary_complete_5D$uni_id==i])
  num_reps[length(num_reps)+1] <- length(mutation_summary_complete_5D$counts[mutation_summary_complete_5D$uni_id==i])
  counts[length(counts)+1] <- mean(mutation_summary_complete_5D$counts[mutation_summary_complete_5D$uni_id==i])
  avg_workload[length(avg_workload)+1] <- mean(mutation_summary_complete_5D$avg_workload[mutation_summary_complete_5D$uni_id==i])
  avg_growth_rate[length(avg_growth_rate)+1] <- mean(mutation_summary_complete_5D$growth_rate[mutation_summary_complete_5D$uni_id==i])
  avg_OGL[length(avg_OGL)+1] <- mean(mutation_summary_complete_5D$avg_OGL[mutation_summary_complete_5D$uni_id==i])
  per_site_mu[length(per_site_mu)+1] <- mean(mutation_summary_complete_5D$per_site_mutation_rate[mutation_summary_complete_5D$uni_id==i])
}

multi_final_df_5D <- cbind.data.frame(uni_id, timepoint, num_reps, counts, avg_workload, avg_growth_rate, avg_OGL, per_site_mu)
head(multi_final_df_5D)

has_soma_df <- cbind.data.frame(timepoint=Fig_5A_data$strategy_rep, has_soma=Fig_5A_data$has_soma.x)
multi_final_df_5D <- left_join(multi_final_df_5D, has_soma_df, by=c("timepoint"))
head(multi_final_df_5D)


### FIGURE 5E DATA###
mutation_summary_complete_5E <- read.csv("mutation_summary_complete_3416_multicells_over_time_filtered.csv", header=T, stringsAsFactors=F)
mutation_summary_complete_5E <- mutation_summary_complete_5E[,-1]
mutation_summary_complete_5E <- mutation_summary_complete_5E[complete.cases(mutation_summary_complete_5E),]
head(mutation_summary_complete_5E)


# AVERAGED BY UNI_ID
head(mutation_summary_complete_5E)

uni_id <- c()
timepoint <- c()
num_reps <- c()
counts <- c()
avg_workload <- c()
avg_growth_rate <- c()
avg_OGL <- c()
per_site_mu <- c()

for (i in unique(mutation_summary_complete_5E$uni_id)) {
  uni_id[length(uni_id)+1] <- i
  timepoint[length(timepoint)+1] <- unique(mutation_summary_complete_5E$timepoint[mutation_summary_complete_5E$uni_id==i])
  num_reps[length(num_reps)+1] <- length(mutation_summary_complete_5E$counts[mutation_summary_complete_5E$uni_id==i])
  counts[length(counts)+1] <- mean(mutation_summary_complete_5E$counts[mutation_summary_complete_5E$uni_id==i])
  avg_workload[length(avg_workload)+1] <- mean(mutation_summary_complete_5E$avg_workload[mutation_summary_complete_5E$uni_id==i])
  avg_growth_rate[length(avg_growth_rate)+1] <- mean(mutation_summary_complete_5E$growth_rate[mutation_summary_complete_5E$uni_id==i])
  avg_OGL[length(avg_OGL)+1] <- mean(mutation_summary_complete_5E$avg_OGL[mutation_summary_complete_5E$uni_id==i])
  per_site_mu[length(per_site_mu)+1] <- mean(mutation_summary_complete_5E$per_site_mutation_rate[mutation_summary_complete_5E$uni_id==i])
}

multi_final_df_5E <- cbind.data.frame(uni_id, timepoint, num_reps, counts, avg_workload, avg_growth_rate, avg_OGL, per_site_mu)
head(multi_final_df_5E)



### FIGURE 5F ###
mutation_summary_complete_5F <- read.csv("mutation_summary_complete_3416_all_timepoints_filtered.csv", header=T, stringsAsFactors=F)
mutation_summary_complete_5F <- mutation_summary_complete_5F[,-1]
mutation_summary_complete_5F <- mutation_summary_complete_5F[complete.cases(mutation_summary_complete_5F),]
head(mutation_summary_complete_5F)

# AVERAGED BY UNI_ID
head(mutation_summary_complete_5F)

uni_id <- c()
timepoint <- c()
num_reps <- c()
counts <- c()
avg_workload <- c()
avg_growth_rate <- c()
avg_OGL <- c()
per_site_mu <- c()

for (i in unique(mutation_summary_complete_5F$uni_id)) {
  uni_id[length(uni_id)+1] <- i
  timepoint[length(timepoint)+1] <- unique(mutation_summary_complete_5F$timepoint[mutation_summary_complete_5F$uni_id==i])
  num_reps[length(num_reps)+1] <- length(mutation_summary_complete_5F$counts[mutation_summary_complete_5F$uni_id==i])
  counts[length(counts)+1] <- mean(mutation_summary_complete_5F$counts[mutation_summary_complete_5F$uni_id==i])
  avg_workload[length(avg_workload)+1] <- mean(mutation_summary_complete_5F$avg_workload[mutation_summary_complete_5F$uni_id==i])
  avg_growth_rate[length(avg_growth_rate)+1] <- mean(mutation_summary_complete_5F$growth_rate[mutation_summary_complete_5F$uni_id==i])
  avg_OGL[length(avg_OGL)+1] <- mean(mutation_summary_complete_5F$avg_OGL[mutation_summary_complete_5F$uni_id==i])
  per_site_mu[length(per_site_mu)+1] <- mean(mutation_summary_complete_5F$per_site_mutation_rate[mutation_summary_complete_5F$uni_id==i])
}

final_df_5F <- cbind.data.frame(uni_id, timepoint, num_reps, counts, avg_workload, avg_growth_rate, avg_OGL, per_site_mu)
head(final_df_5F)


####
# Set shared scale for growth rate in Figure 5D, 5E, and 5F
min_growth_rate <- min(c(multi_final_df_5E$avg_growth_rate, final_df_5F$avg_growth_rate))
max_growth_rate <- max(c(multi_final_df_5E$avg_growth_rate, final_df_5F$avg_growth_rate))
myPalette <- viridis_pal(option="B")
sc <- scale_colour_gradientn(colours = myPalette(1000), limits=c(min_growth_rate, max_growth_rate), name="Average \ngrowth \nrate")
min_y_lim_full <- min(c(multi_final_df_5D$avg_workload, multi_final_df_5E$avg_workload, final_df_5F$avg_workload))
max_y_lim_full <- max(c(multi_final_df_5D$avg_workload, multi_final_df_5E$avg_workload, final_df_5F$avg_workload))
min_y_lim_part <- min(c(multi_final_df_5D$avg_workload, multi_final_df_5E$avg_workload, final_df_5F$avg_workload[final_df_5F$counts>30]))
max_y_lim_part <- max(c(multi_final_df_5D$avg_workload, multi_final_df_5E$avg_workload, final_df_5F$avg_workload[final_df_5F$counts>30]))
min_x_lim_full <- min(c(1-multi_final_df_5D$per_site_mu, 1-multi_final_df_5E$per_site_mu, 1-final_df_5F$per_site_mu))
max_x_lim_full <- max(c(1-multi_final_df_5D$per_site_mu, 1-multi_final_df_5E$per_site_mu, 1-final_df_5F$per_site_mu))
min_x_lim_part <- min(c(1-multi_final_df_5D$per_site_mu, 1-multi_final_df_5E$per_site_mu, 1-final_df_5F$per_site_mu[final_df_5F$counts>30]))
max_x_lim_part <- max(c(1-multi_final_df_5D$per_site_mu, 1-multi_final_df_5E$per_site_mu, 1-final_df_5F$per_site_mu[final_df_5F$counts>30]))


### FIGURE 5D PLOTS ###
# ALL DATA
Fig_5D_log <- ggplot(data=multi_final_df_5D, aes(x=1-per_site_mu, y=avg_workload+1, color=as.factor(has_soma))) +
  scale_color_manual(values=c("#608fb1", "#a0352b")) +  
  scale_y_continuous(trans="log10", limits=c(1, max_y_lim_full+1)) +
  scale_x_continuous(limits=c(min(1-multi_final_df_5D$per_site_mu), 1.0), breaks=c(0.925,0.950,0.975,1.000)) +
  geom_point(size=0.5, alpha=0.8) +
  geom_line(data=genome_length_100, aes(x=fidelity, y=work), col="grey60", lty=1, lwd=0.33) +
  #geom_line(data=genome_length_50[test_df$work<=max(multi_final_df_5D$avg_workload[multi_final_df_5D$counts>0]),], aes(x=fidelity, y=work), col="blue", lty=1, lwd=1) +
  #geom_line(data=test_df[test_df$test_fidelity>=0,], 
  #          aes(x=test_fidelity, y=test_work), col="dodgerblue", lty=1, lwd=1) +
  xlab(expression(paste("Fidelity (1"-"µ)"))) +
  ylab("Average workload + 1") +
  theme(legend.position = "none") 


### FIGURE 5E PLOTS ###

Fig_5E_log <- ggplot(data=multi_final_df_5E, aes(x=1-per_site_mu, y=avg_workload+1)) +
  sc +
  #scale_color_viridis_c(option="B", name = "Revertant\ngrowth\nrate") +
  scale_y_continuous(trans="log10", limits=c(1, max_y_lim_full+1)) +
  scale_x_continuous(limits=c(min(1-multi_final_df_5D$per_site_mu), 1.0), breaks=c(0.925,0.950,0.975,1.000)) +
  geom_point(size=0.5, alpha=0.8, aes(col=avg_growth_rate)) +
  #geom_text_repel(data=multi_final_df_5E,
  #                position="jitter",
  #                label=unique(timepoint), 
  #                aes(color=avg_growth_rate), 
  #                size=2, max_overlap=20) +
  geom_segment(aes(color=avg_growth_rate,
                 xend=c(tail(1-per_site_mu, n=-1), NA), 
                 yend=c(tail(avg_workload, n=-1), NA)),
               lwd=0.25, alpha=0.5) +
  geom_line(data=genome_length_100, aes(x=fidelity, y=work), col="grey60", lty=1, lwd=0.33) +
  #geom_line(data=genome_length_50[test_df$work<=max(multi_final_df_5E$avg_workload[multi_final_df_5E$counts>0]),], aes(x=fidelity, y=work), col="blue", lty=1, lwd=1) +
  #geom_line(data=test_df[test_df$test_fidelity>=0,], 
  #          aes(x=test_fidelity, y=test_work), col="dodgerblue", lty=1, lwd=1) +
  xlab(expression(paste("Fidelity (1"-"µ)"))) +
  ylab("Average workload + 1") +
  theme(legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank()) 


### FIGURE 5F PLOTS ###
Fig_5F_log <- ggplot(data=final_df_5F, aes(x=1-per_site_mu, y=avg_workload+1)) +
  sc +
  #scale_color_viridis_c(option="B", name = "Revertant\ngrowth\nrate") +
  scale_y_continuous(trans="log10", limits=c(1, max_y_lim_full+1)) +
  scale_x_continuous(limits=c(min(1-multi_final_df_5D$per_site_mu), 1.0), breaks=c(0.925,0.950,0.975,1.000)) +
  geom_point(size=0.5, alpha=0.8, aes(col=avg_growth_rate)) +
  geom_line(data=genome_length_100, aes(x=fidelity, y=work), col="grey60", lty=1, lwd=0.33) +
  #geom_rect(xmin=0.95, xmax=1.015, ymin=-0.05, ymax=log(max(final_df_5F$avg_workload[final_df_5F$counts>30]),10), fill=NA, size=0.5, color="orange", alpha=1) +
  #geom_line(data=genome_length_50[test_df$work<=max(final_df_5F$avg_workload[final_df_5F$counts>0]),], aes(x=fidelity, y=work), col="blue", lty=1, lwd=1) +
  #geom_line(data=test_df[test_df$test_fidelity>=0,], 
  #          aes(x=test_fidelity, y=test_work), col="dodgerblue", lty=1, lwd=1) +
  #geom_smooth(method=MASS::rlm, formula=y~x, col="red") +
  xlab(expression(paste("Fidelity (1"-"µ)"))) +
  ylab("Average workload + 1") +
  theme(legend.title = element_text(size=5),
        legend.text=element_text(size=5),
        axis.title.y=element_blank(),
        axis.text.y=element_blank()) 



Fig_5B <- ggplot() + theme_void()
Fig_5C <- ggplot() + theme_void()

Fig_5_all <- ggarrange(Fig_5A, Fig_5B, Fig_5C, Fig_5D_log, Fig_5E_log, Fig_5F_log, nrow=2)

setwd(paste(pwd, "/", fwd[[1]][length(fwd[[1]])], "_drafts", sep=""))

pdf("Figure_5_Entrenchment_log_scale.pdf", width=two_col, height=two_col*0.715)
ggdraw() +
  draw_plot(Fig_5_all, x=0.01, y=0.01, width=0.98, height=0.98) +
  draw_plot_label(label = c("A", "B", "C", "D", "E", "F"), size = 9,
                  x = c(0, 0.315, 0.6, 0, 0.315, 0.6), 
                  y = c(0.9, 0.9, 0.9, 0.46, 0.46, 0.46))
dev.off()


#------------------------------------------------
# STATISTICS
# Spearman's rank order correlation
# ∆ Division of labor vs. ∆ Entrenchment
x <- Fig_5A_data$avg_total_SD_final - Fig_5A_data$avg_total_SD_trans
y <- Fig_5A_data$logistic_entrench_final - Fig_5A_data$logistic_entrench_trans

cor.test(x=x, y=y, method="spearman")


# Wilcoxon rank sum test
# Entrenchment in undifferentiated vs. differentiated isolates
x <- Fig_5A_data$avg_total_SD_trans[Fig_5A_data$has_soma.x==0] - Fig_5A_data$avg_total_SD_final[Fig_5A_data$has_soma.x==0]
y <- Fig_5A_data$avg_total_SD_trans[Fig_5A_data$has_soma.x==1] - Fig_5A_data$avg_total_SD_final[Fig_5A_data$has_soma.x==1]

test <- wilcox.test(x=x, y=y, paired=FALSE)
Zstat<-qnorm(test$p.value/2)
test   #print p-value
abs(Zstat)/sqrt(max(length(x), length(y)))    #print effect size

