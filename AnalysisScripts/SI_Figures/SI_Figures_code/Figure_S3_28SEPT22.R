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


# Set working directory
setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/entrenchment_data/2020-04-01_update_lod_fitness")

entrenchment_final <- read.csv("entrenchment_data_final_15JULY20.csv", header=TRUE, stringsAsFactors=FALSE)
entrenchment_trans <- read.csv("entrenchment_data_transition_15JULY20.csv", header=TRUE, stringsAsFactors=FALSE)


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

ggplot() +
  geom_line(data=trans_logistic_model_predictions[trans_logistic_model_predictions$X1_range <= 2048,], aes(x=X1_range, y=res, color=timepoint), alpha=0.75, lwd=1.5, lty=1) +
  geom_point(data=data[data$cost <= 2048,], aes(x=cost, y=reverted), color=rgb(0,0,0,0.33), size=1) +
  scale_color_manual(values=c("darkgoldenrod1", "green3")) +
  scale_x_continuous(trans="log2") +
  xlab("Cost") +
  ylab("Probability of reversion") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  facet_wrap(~rep)

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


ggplot() +
  geom_line(data=final_logistic_model_predictions[final_logistic_model_predictions$X1_range <= 2048,], aes(x=X1_range, y=res, color=timepoint), alpha=0.75, lwd=1.5, lty=1) +
  geom_point(data=data[data$cost <= 2048,], aes(x=cost, y=reverted), color=rgb(0,0,0,0.33), size=1) +
  scale_color_manual(values=c("green3")) +
  scale_x_continuous(trans="log2") +
  xlab("Cost") +
  ylab("Probability of reversion") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  facet_wrap(~rep)

final_inflection <- final_entrenchment_values$inflection_point[final_entrenchment_values$strategy_rep==case_study]




final_logistic_model_predictions <- left_join(final_logistic_model_predictions, 
                                              cbind.data.frame(rep = unique(final_logistic_model_predictions$rep), 
                                                               strain_ID = c(outer(c("DD."), c(1:length(unique(final_logistic_model_predictions$rep))), FUN=paste, sep=""))
                                                               ), 
                                              by="rep")

final_logistic_model_predictions$strain_ID <- factor(final_logistic_model_predictions$strain_ID, levels=c(outer(c("DD."), c(1:length(unique(final_logistic_model_predictions$rep))), 
                                                                                                                FUN=paste, sep="")))

levels(final_logistic_model_predictions$strain_ID)

trans_logistic_model_predictions <- left_join(trans_logistic_model_predictions, 
                                              cbind.data.frame(rep = unique(trans_logistic_model_predictions$rep), 
                                                               strain_ID = c(outer(c("DD."), c(1:length(unique(trans_logistic_model_predictions$rep))), 
                                                                                   FUN=paste, sep=""))
                                              ), by="rep")



setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Figures")


final_logistic_model_predictions$label <- factor(final_logistic_model_predictions$rep, levels = unique(final_logistic_model_predictions$rep), 
                                     labels = paste("DW.", order(unique(final_logistic_model_predictions$rep)), sep=""))
trans_logistic_model_predictions$label <- factor(trans_logistic_model_predictions$rep, levels = unique(trans_logistic_model_predictions$rep), 
                                                 labels = paste("DW.", order(unique(trans_logistic_model_predictions$rep)), sep=""))
data$label <- factor(data$rep, levels = unique(data$rep), 
                                                 labels = paste("DW.", order(unique(data$rep)), sep=""))
final_data$label <- factor(final_data$rep, levels = unique(final_data$rep), 
                     labels = paste("DW.", order(unique(final_data$rep)), sep=""))




pdf("Figure_S3_Dirty_work_Entrenchment_20DEC22.pdf", width=one_col*3, height=one_col*3.5)
ggplot(data=final_logistic_model_predictions[final_logistic_model_predictions$X1_range <= 2048,]) +
  geom_line(data=final_logistic_model_predictions[final_logistic_model_predictions$X1_range <= 2048,], 
            aes(x=X1_range, y=res, color=timepoint), alpha=0.75, lwd=0.5, lty=1) +
  geom_line(data=trans_logistic_model_predictions[trans_logistic_model_predictions$X1_range <= 2048,], 
            aes(x=X1_range, y=res, color=timepoint), alpha=0.75, lwd=0.5, lty=1) +
  geom_point(data=data, aes(x=cost, y=reverted+0.01, color=timepoint), color="goldenrod1", size=0.5, alpha=0.5) +
  geom_point(data=final_data[final_data$cost <= 2048,], aes(x=cost, y=reverted-0.01), color="green3", size=0.5, alpha=0.5) +
  scale_color_manual(values=c("green3", "goldenrod1"), name = "Timepoint", labels = c("Final", "Transition")) +
  scale_x_continuous(trans="log2") +
  xlab("Cost (updates)") +
  ylab("Probability of reversion") +
  # geom_hline(yintercept=0.5, linetype=2, color="black", size=0.5) +
  # geom_point(aes(x=final_inflection, y=0.5),
  #            pch=21, fill="white", color="green3", size=1, stroke=1.25) +
  # geom_point(aes(x=trans_inflection, y=0.5),
  #            pch=21, fill="white", color="goldenrod1", size=1, stroke=1.25) +
  # annotate(geom="text", x=trans_inflection-25, y=0.555, color="darkgoldenrod1", size=2, vjust=0, hjust=1,
  #          label=format(round(trans_inflection, digits=1), nsmall=1)) +
  # annotate(geom="text", x=final_inflection+95, y=0.555, color="green3", size=2, vjust=0, hjust=0,
  #          label=format(round(final_inflection, digits=1), nsmall=1)) +
  theme(axis.text.x = element_text(angle=0), axis.title=element_text(size=10), legend.position=c(0.9, 0.025)) +
  facet_wrap(~label, ncol=7)
dev.off()





logistic_regression_entrenchment_values <- rbind.data.frame(transition_entrenchment_values, final_entrenchment_values)
setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/Figure_3/Figure_3_data")
write.csv(logistic_regression_entrenchment_values, "Entrenchment_values_from_logistic_regression_14FEB22.csv")


plot_data <- full_join(transition_entrenchment_values, final_entrenchment_values, by=c("strategy_rep"))
plot_data_full <- cbind(plot_data, rep("diff", length(transition_entrenchment_values$strategy_rep)))
colnames(plot_data_full)[length(plot_data_full)] <- "value"


ggplot(data=plot_data_full, aes(x=value, y=inflection_point.y - inflection_point.x)) +
  geom_violin(fill="grey90", scale="count", col=NA, trim=FALSE) +
  geom_hline(yintercept=0, linetype=2, color="black", size=0.5) +
  geom_jitter(width=0.35, 
              height=0.1, 
              color="grey40", 
              size=1,
              pch=19) +
  xlab("") +
  #scale_y_continuous(breaks=c(0,2,4,6,8)) +
  ylab(expression(paste(Delta, " Entrenchment"))) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.border = element_blank())


mean(transition_entrenchment_values$inflection_point)
sd(transition_entrenchment_values$inflection_point)
mean(final_entrenchment_values$inflection_point)
sd(final_entrenchment_values$inflection_point)
wilcox.test(transition_entrenchment_values$inflection_point, 
            final_entrenchment_values$inflection_point,
            paired=TRUE, alternative="two.sided", exact=FALSE)






setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/Figure_2/Figure_2_drafts")

# # Figure 2 new
# pdf("Figure_2_Entrenchment_FINAL_v3.pdf", width=two_col, height=two_col*.7)
# ggdraw() +
#   #draw_plot(Fig_2a, x = 0.01, y = 0.66, width = 0.30, height = 0.30) +
#   draw_plot(Fig_2b_3, x = 0.01, y=0.66, width = 0.99, height = 0.30) +
#   draw_plot(Fig_2b_new, x = 0.01, y = 0.335, width = 0.64, height = 0.295) +
#   draw_plot(Fig_2e_new, x = 0.68, y = 0.33, width = 0.31, height = 0.30) +
#   draw_plot(Fig_2c_new, x = 0.01, y = 0.005, width = 0.64, height = 0.295) +
#   draw_plot(Fig_2f_new, x = 0.68, y = 0.0, width = 0.31, height = 0.30) +
#   draw_plot_label(label = c("A", " ", "B", "C", "D", "E"), size = 9, 
#                   family = "Helvetica", fontface = "bold",
#                   x = c(0, 0.3, 0, 0.67, 0, 0.67), 
#                   y = c(1, 1, 0.65, 0.65, 0.32, 0.32))
# dev.off()






# ADDITIONAL FIGURE PANEL THAT SHOULD BE ADDED TO FIGURE 2

Fig_2E <- ggplot(data=plot_data_full, aes(x=inflection_point.x, y=inflection_point.y)) +
  scale_x_continuous(trans="log2", limits=c(1,1024), 
                     expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(trans="log2", limits=c(1,1024), 
                     expand = expansion(mult = c(0, 0.05))) +
  #geom_abline(slope=1, intercept=c(0,0), na.rm = FALSE, show.legend = NA, color="black", lty=2, lwd=1.15) +
  #geom_hline(yintercept=plot_data_full$inflection_point.y[plot_data_full$strategy_rep==3416], lty=2, size=0.25, color="green3") +
  #geom_vline(xintercept=plot_data_full$inflection_point.x[plot_data_full$strategy_rep==3416], lty=2, size=0.25, color="goldenrod1") +
  geom_point(color="grey40", size=3, alpha=0.8) +
  # geom_point(color=ifelse(plot_data_full$strategy_rep==3416, "#FF6900", "grey40"),
  #          pch=ifelse(plot_data_full$strategy_rep==3416, 18, 19),
  #          size=ifelse(plot_data_full$strategy_rep==3416, 2.5, 0.4)) +
  xlab("Stability (transition)") +
  ylab("Stability (final)")



pdf("Figure_2_Entrenchment_FINAL_v4.pdf", width=two_col, height=two_col*.7)
ggdraw() +
  draw_plot(Fig_2A, x = 0.01, y=0.725, width = 0.99, height = 0.275) +
  draw_plot(Fig_2B, x = 0.01, y = 0.4, width = 0.67, height = 0.325) +
  draw_plot(Fig_2C, x = 0.01, y = 0.03, width = 0.67, height = 0.325) +
  draw_plot(Fig_2D, x = 0.68, y = 0.355, width = 0.31, height = 0.37) +
  draw_plot(Fig_2E, x = 0.70, y = 0.0, width = 0.29, height = 0.37) +
  draw_plot_label(label = c("A", " ", "B", "D", "C", "E"), size = 9, 
                  family = "Helvetica", fontface = "bold",
                  x = c(0, 0.3, 0, 0.68, 0, 0.68), 
                  y = c(1, 1, 0.725, 0.725, 0.375, 0.375))
dev.off()



# RECONFIGURE LAYOUT & COMBINE PLOTS B AND C

pdf("Figure_2_Entrenchment_FINAL_v5.pdf", width=two_col, height=two_col*0.9)
ggdraw() +
  draw_plot(Fig_2A, x = 0.035, y=0.775, width = 0.965, height = 0.225) +
  draw_plot(Fig_2BC, x= 0.05, y=0.50, width = 0.95, height = 0.275) +
  #draw_plot(Fig_2B, x = 0.01, y = 0.4, width = 0.67, height = 0.325) +
  #draw_plot(Fig_2C, x = 0.01, y = 0.03, width = 0.67, height = 0.325) +
  draw_plot(Fig_2D, x = 0.01, y = 0.0, width = 0.475, height = 0.475) +
  draw_plot(Fig_2E, x = 0.525, y = 0.0, width = 0.475, height = 0.475) +
  draw_plot_label(label = c("A", " ", "B", " ", "C", "D"), size = 9, 
                  family = "Helvetica", fontface = "bold",
                  x = c(0, 0.3, 0, 0.5, 0, 0.5), 
                  y = c(1, 1, 0.75, 0.75, 0.5, 0.5))
dev.off()





#------------------------------------------------
# STATISTICS

# Wilcoxon signed rank test
# Stability @ transition vs. final timepoint
x <- plot_data_full$inflection_point.x
y <- plot_data_full$inflection_point.y

test <- wilcox.test(x=x, y=y, paired=TRUE)
Zstat<-qnorm(test$p.value/2)
test                          #print p-value
abs(Zstat)/sqrt(length(x))    #print effect size
