# Figure 2, plots
# Peter Conlin
# 17 SEPT 21

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
theme_set(theme_cowplot(font_size = 18, font_family="Helvetica"))
library(viridis)

one_col <- 2.25
two_col <- 4.75

#----------------------------------------

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot(font_size=7))
library(dplyr)
library(directlabels)
library(reshape)
library(viridis)

# If running in RStudio, you can use the following command to setwd() to source file location:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('..')    #Go up a directory
pwd <- getwd()    #Set parent working directory
fwd <- c(strsplit(pwd, split="/"))     #Get a list of directories in file path

## Set working directory to load required data sets
setwd(paste(pwd, "/", fwd[[1]][length(fwd[[1]])], "_data", sep=""))
getwd()
#----------------------------------------

case_study_overtime <- read.csv("3416_over_time_size_workload.csv", header=T, stringsAsFactors=F)
head(case_study_overtime)

Fig_2A <- ggplot(data=case_study_overtime, aes(x=timepoint, y=multicell_size)) +
  xlab("Time (generations)") +
  ylab("Avg. size \nat reproduction") +
  scale_y_continuous(limits=c(0,25), breaks=c(0,5,15,25),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(trans = "log10",
                     expand = expansion(mult = c(0, 0.05))) +
  geom_vline(xintercept = 35, lty=2, color="goldenrod1", linewidth=0.25) +
  geom_vline(xintercept = max(case_study_overtime$timepoint), lty=2, color="green3", linewidth=0.25) +
  annotate(geom = "text", x = 36, y = 25, label = "transition", color = "goldenrod1", angle = 90, hjust=1, vjust=-1, size=2) +
  annotate(geom = "text", x = max(case_study_overtime$timepoint), y = 25, label = "final", color = "green3", angle = 90, hjust=1, vjust=1, size=2) +
  geom_line(alpha=0.9) 


final_full_df <- read.table("3416_final_lod_entrench_all.dat", header=T, stringsAsFactors=F)
trans_full_df <- read.table("3416_transition_lod_entrench_all.dat", header=T, stringsAsFactors=F)
T_df <- cbind.data.frame(trans_full_df, label=rep("trans", length(trans_full_df$cost)))
F_df <- cbind.data.frame(final_full_df, label=rep("final", length(final_full_df$cost)))
full_df <- rbind.data.frame(T_df, F_df)

Fig_2BC <- ggplot(data=full_df[full_df$cost > 32 & full_df$cost < 2048,], aes(x=generation_diff, y=organism_size, group=interaction(label, iteration))) +
  geom_line(aes(x=generation_diff, y=organism_size, color=label), 
            size=0.35, alpha=0.66) +
  xlab("Time (generations)") +
  ylab("Avg. size") +
  scale_y_continuous(limits=c(0,25), breaks=c(0,5,15,25),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(breaks=c(50,100),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_color_manual(values=c("green3", "goldenrod1")) +
  facet_wrap(~cost, nrow=1)  + 
  theme(legend.position = "none",
        strip.background = element_rect(color=NA, 
                                        fill=NA, size=1.5, 
                                        linetype="solid"),
        panel.background = element_rect(fill="white", colour="white")) +
  panel_border(color="grey80", size=0.5)
  



entrenchment_final <- read.csv("entrenchment_data_final.csv", header=TRUE, stringsAsFactors=FALSE)
entrenchment_trans <- read.csv("entrenchment_data_transition.csv", header=TRUE, stringsAsFactors=FALSE)

# Loop through entrechment data and run logisitc regression
tp <- "trans"
data <- entrenchment_trans
trans_logistic_model_predictions <- data.frame(X1_range=numeric(),
                                         res=character(),
                                         rep=numeric(),
                                         timepoint=character())
strategy_rep <- c()
inflection_point <- c()

for (i in unique(data$rep)) {
  
  # Run logistic regression, save predictions
  model <- glm(formula = reverted ~ cost, family = binomial, data = data[data$rep==i,])
  X1_range <- seq(from=min(data[data$rep==i,]$cost), to=max(data[data$rep==i,]$cost), by=1)
  res <- predict(model, list(cost=X1_range), type="response")
  rep <- rep(i, length(X1_range))
  timepoint <- rep(tp, length(X1_range))
  df <- cbind.data.frame(X1_range, res, rep, timepoint)
  trans_logistic_model_predictions <- rbind.data.frame(trans_logistic_model_predictions, df)
  
  # Calculate inflection point, store in a list
  b0 <- model$coef[1]
  X1 <- model$coef[2]
  p <- 0.5
  strategy_rep[length(strategy_rep)+1] <- i
  inflection_point[length(inflection_point)+1] <- (log(p/(1-p)) - b0 / X1)
  
}

transition_entrenchment_values <- cbind.data.frame(strategy_rep, inflection_point, timepoint=rep(0, length(strategy_rep)))

case_study <- 3416
trans_inflection <- transition_entrenchment_values$inflection_point[transition_entrenchment_values$strategy_rep==case_study]


# Loop through entrechment data and run logisitc regression
tp <- "final"
final_data <- entrenchment_final
final_logistic_model_predictions <- data.frame(X1_range=numeric(),
                                         res=character(),
                                         rep=numeric(),
                                         timepoint=character())
strategy_rep <- c()
inflection_point <- c()

for (i in unique(final_data$rep)) {
  
  # Run logistic regression, save predictions
  model <- glm(formula = reverted ~ cost, family = binomial, data = final_data[final_data$rep==i,])
  X1_range <- seq(from=min(final_data[final_data$rep==i,]$cost), to=max(final_data[final_data$rep==i,]$cost), by=1)
  res <- predict(model, list(cost=X1_range), type="response")
  rep <- rep(i, length(X1_range))
  timepoint <- rep(tp, length(X1_range))
  df <- cbind.data.frame(X1_range, res, rep, timepoint)
  final_logistic_model_predictions <- rbind.data.frame(final_logistic_model_predictions, df)
  
  # Calculate inflection point, store in a list
  b0 <- model$coef[1]
  X1 <- model$coef[2]
  p <- 0.5
  strategy_rep[length(strategy_rep)+1] <- i
  inflection_point[length(inflection_point)+1] <- (log(p/(1-p)) - b0 / X1)
  
}


final_entrenchment_values <- cbind.data.frame(strategy_rep, inflection_point, timepoint=rep(1, length(strategy_rep)))

final_inflection <- final_entrenchment_values$inflection_point[final_entrenchment_values$strategy_rep==case_study]


Fig_2D <- ggplot() +
  geom_line(data=final_logistic_model_predictions[final_logistic_model_predictions$rep==case_study & final_logistic_model_predictions$X1_range <= 2048,], 
            aes(x=X1_range, y=res, color=timepoint), alpha=0.75, lwd=0.5, lty=1) +
  geom_line(data=trans_logistic_model_predictions[trans_logistic_model_predictions$rep==case_study & trans_logistic_model_predictions$X1_range <= 2048,], 
            aes(x=X1_range, y=res, color=timepoint), alpha=0.75, lwd=0.5, lty=1) +
  geom_point(data=data[data$rep==case_study & data$cost <= 2048,], aes(x=cost, y=reverted+0.01, color=timepoint), color="goldenrod1", size=0.5, alpha=0.5) +
  geom_point(data=final_data[final_data$rep==case_study & final_data$cost <= 2048,], aes(x=cost, y=reverted-0.01), color="green3", size=0.5, alpha=0.5) +
  scale_color_manual(values=c("green3", "goldenrod1")) +
  scale_x_continuous(trans="log2") +
  xlab("Cost (updates)") +
  ylab("Probability \nof reversion") +
  geom_hline(yintercept=0.5, linetype=2, color="black", size=0.5) +
  geom_point(aes(x=final_inflection, y=0.5),
             pch=21, fill="white", color="green3", size=1, stroke=1.25) +
  geom_point(aes(x=trans_inflection, y=0.5),
             pch=21, fill="white", color="goldenrod1", size=1, stroke=1.25) +
  annotate(geom="text", x=trans_inflection-25, y=0.555, color="darkgoldenrod1", size=2, vjust=0, hjust=1,
           label=format(round(trans_inflection, digits=1), nsmall=1)) +
  annotate(geom="text", x=final_inflection+95, y=0.555, color="green3", size=2, vjust=0, hjust=0,
           label=format(round(final_inflection, digits=1), nsmall=1)) +
  theme(axis.text.x = element_text(angle=0), legend.position = "none") 


plot_data <- full_join(transition_entrenchment_values, final_entrenchment_values, by=c("strategy_rep"))
plot_data_full <- cbind(plot_data, rep("diff", length(transition_entrenchment_values$strategy_rep)))
colnames(plot_data_full)[length(plot_data_full)] <- "value"

Fig_2E <- ggplot(data=plot_data_full, aes(x=inflection_point.x, y=inflection_point.y)) +
  scale_x_continuous(trans="log2", limits=c(1,1024), 
                     expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(trans="log2", limits=c(1,1024), 
                     expand = expansion(mult = c(0, 0.05))) +
  geom_abline(slope=1, intercept=c(0,0), na.rm = FALSE, show.legend = NA, color="grey80", lty=1, linewidth=0.25) +
  geom_hline(yintercept=plot_data_full$inflection_point.y[plot_data_full$strategy_rep==3416], lty=2, linewidth=0.25, color="green3") +
  geom_vline(xintercept=plot_data_full$inflection_point.x[plot_data_full$strategy_rep==3416], lty=2, linewidth=0.25, color="goldenrod1") +
  geom_point(color="grey40", size=0.5, alpha=1) +
   geom_point(color=ifelse(plot_data_full$strategy_rep==3416, "#FF6900", "grey40"),
            pch=ifelse(plot_data_full$strategy_rep==3416, 18, 19),
            size=ifelse(plot_data_full$strategy_rep==3416, 2.5, 0.4)) +
  xlab("Stability (transition)") +
  ylab("Stability (final)")



setwd(paste(pwd, "/", fwd[[1]][length(fwd[[1]])], "_drafts", sep=""))


pdf("Figure_2_Entrenchment_FINAL_v5.pdf", width=two_col, height=two_col*0.9)
ggdraw() +
  draw_plot(Fig_2A, x = 0.035, y=0.775, width = 0.965, height = 0.225) +
  draw_plot(Fig_2BC, x= 0.05, y=0.50, width = 0.95, height = 0.275) +
  draw_plot(Fig_2D, x = 0.01, y = 0.0, width = 0.475, height = 0.475) +
  draw_plot(Fig_2E, x = 0.525, y = 0.0, width = 0.475, height = 0.475) +
  draw_plot_label(label = c("A", " ", "B", " ", "C", "D"), size = 9, 
                  family = "Helvetica", fontface = "bold",
                  x = c(0, 0.3, 0, 0.5, 0, 0.5), 
                  y = c(1, 1, 0.75, 0.75, 0.5, 0.5))
dev.off()





#------------------------------------------------
# STATISTICS

mean(transition_entrenchment_values$inflection_point)
sd(transition_entrenchment_values$inflection_point)
mean(final_entrenchment_values$inflection_point)
sd(final_entrenchment_values$inflection_point)
wilcox.test(transition_entrenchment_values$inflection_point, 
            final_entrenchment_values$inflection_point,
            paired=TRUE, alternative="two.sided", exact=FALSE)

# Wilcoxon signed rank test
# Stability @ transition vs. final timepoint
x <- plot_data_full$inflection_point.x
y <- plot_data_full$inflection_point.y

test <- wilcox.test(x=x, y=y, paired=TRUE)
Zstat<-qnorm(test$p.value/2)
test                          #print p-value
abs(Zstat)/sqrt(length(x))    #print effect size
