## Supplemental figure examining position and identity of reversion mutations at transition and final timepoints
## 07 September 2021
## Peter Conlin

# Libraries
library(ggplot2)
library(GGally)
library(viridis)
library(cowplot)
theme_set(theme_cowplot(font_size = 9))
# Libraries needed for Sankey diagrams
library(tidyverse)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(networkD3)
library(htmlwidgets)
library(egg)

one_col <- 2.25
two_col <- 4.75

setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Figure_Reversion_mutation/Data/")
#setwd("C:/Users/pconlin3/Dropbox/Desktop 28NOV21/Entrenchment_processed_data_for_figures/SI_Figure_Reversion_mutation/Data/") #For PC

# Read in compiled reversion mutation data
all_reversion_mutations <- read.csv("all_reversion_mutation_identities_14SEPT21.csv", header=T, stringsAsFactors=F)
all_reversion_mutations <- all_reversion_mutations[,-1] #remove indexing column added when data was read in from csv
all_reversion_mutations$genome_location <- all_reversion_mutations$genome_location + 1 #GENOME INDEXES BEGIN AT 0 FOR COMPILED REVERSION MUTATION DATA BUT AT 1 FOR ANCESTOR GENOMES
head(all_reversion_mutations)
147/length(all_reversion_mutations$strategy_rep)
unique(sort(all_reversion_mutations$genome_location))

rev_muts_4310 <- all_reversion_mutations[all_reversion_mutations$strategy_rep==4310 & all_reversion_mutations$timepoint==0,]
rev_muts_3038 <- all_reversion_mutations[all_reversion_mutations$strategy_rep==3038 & all_reversion_mutations$timepoint==1,]
rev_muts_4278 <- all_reversion_mutations[all_reversion_mutations$strategy_rep==4278 & all_reversion_mutations$timepoint==1,]
rev_muts_3416 <- all_reversion_mutations[all_reversion_mutations$strategy_rep==3416 & all_reversion_mutations$timepoint==1,]


# Read in multicell ancestor genomes from transition and final timepoint
multi_genomes <- read.csv("multicelled_ancestor_genomes_at_transition_and_final_timepoint_29NOV21.csv", header=T, stringsAsFactors=F)
multi_genomes <- multi_genomes[,-1]
head(multi_genomes)


unique(sort(multi_genomes$genome_location))

multi_genome_4310 <- multi_genomes[multi_genomes$strategy_rep==4310 & multi_genomes$timepoint==0,]
multi_genome_3038 <- multi_genomes[multi_genomes$strategy_rep==3038 & multi_genomes$timepoint==1,]
multi_genome_4278 <- multi_genomes[multi_genomes$strategy_rep==4278 & multi_genomes$timepoint==1,]


# Merge the two files together using left_join
rev_mut_df <- left_join(x=all_reversion_mutations, y=multi_genomes, by=c("strategy_rep", "timepoint", "genome_location"))
head(rev_mut_df)


# Check that none of the mutations match the wt_mut column
z <- rev_mut_df$mutation==rev_mut_df$wt_mut
length(z[z==TRUE])
# Find the identity of the rows that are problematic
problems <- rev_mut_df[which(rev_mut_df$mutation==rev_mut_df$wt_mut),] # This produced a table w/ all NA values
problems <- rev_mut_df[rowSums(is.na(rev_mut_df)) > 0,] # This produced the expected result, 
                                                        # revealing 147 columns for which no matching multi_genome position existed
write.csv(problems, "unmappable_mutations_18JAN23.csv")


# Read in verbose list of AVIDA instructions
AVIDA_instructions <- read.csv("Avida_instructions_final_20DEC22.csv", header=T, stringsAsFactors=F)
head(AVIDA_instructions)

# Add AVIDA operations for mutation and wt_mut
avida_mut <- AVIDA_instructions
colnames(avida_mut) <- c("mutation", "mut_op")
rev_mut_df <- left_join(rev_mut_df, avida_mut, by=c("mutation"))
avida_wt <- AVIDA_instructions
colnames(avida_wt) <- c("wt_mut", "wt_op")
rev_mut_df <- left_join(rev_mut_df, avida_wt, by=c("wt_mut"))
# Make unique columns for operation position pairs
rev_mut_df$wt_ID <- paste("p", rev_mut_df$genome_location,"; ", rev_mut_df$wt_op, sep="")

example <- rev_mut_df[rev_mut_df$strategy_rep==3416,]
write.csv(example, "3416_data_for_sankey_plot_example_20DEC22.csv")
###-------------------------------------------------------
# MAKING SANKEY DIAGRAMS

# Re-name all_reversion_mutations before modifying
sankey <- rev_mut_df[complete.cases(rev_mut_df),]

# Add color vector to AVIDA instructions for plotting mutation identity in the sankey diagrams
AVIDA_instructions$color <- viridis(n=38, option="D")

for (i in unique(sankey$strategy_rep)) {
  # Get transition data
  data <- sankey[sankey$strategy_rep==i & sankey$timepoint==0,]
  data <- cbind.data.frame(data, value=rep(1, length(data$strategy_rep)))
  data_nodes <- cbind.data.frame(name=c(unique(data$wt_ID), unique(data$mut_op)), empty=rep(NA,length(unique(data$wt_ID)) + length(unique(data$mut_op))))
  data$IDsource=match(data$mut_op, data_nodes$name)-1 
  data$IDtarget=match(data$wt_ID, data_nodes$name)-1 
  
  # prepare D3 colour scale for transition timepoint
  op_colors <- c()
  
  for (j in data_nodes$name) {
    if (j %in% data$wt_ID) {
      temp <- unique(data$wt_ID[data$wt_ID==j])
      temp <- unique(data$wt_op[data$wt_ID==temp])
      op_colors[length(op_colors)+1] <- AVIDA_instructions$color[AVIDA_instructions$operation==temp]
    }
    else if (j %in% data$mut_op) {
      temp <- unique(data$mut_op[data$mut_op==j])
      op_colors[length(op_colors)+1] <- AVIDA_instructions$color[AVIDA_instructions$operation==temp]
    }
  }
  
  Color_scale_trans = paste('d3.scaleOrdinal() .domain([', 
                            gsub("^c\\(|\\)$", "", toString(paste(list(data_nodes$name), sep=""))),
                            ']) .range([',
                            gsub("^c\\(|\\)$", "", toString(list(op_colors))),
                            '])', sep="")
  
  # Make the Network
  sankeyNetwork(Links = data, 
                Nodes = data_nodes, 
                Source = "IDsource", 
                Target = "IDtarget",
                Value = "value", 
                NodeID = "name", 
                sinksRight=FALSE, 
                colourScale=Color_scale_trans, 
                nodeWidth=80, 
                fontSize=18, 
                nodePadding=8, 
                iterations = 0, 
                margin = 200,
                height = 20*max(c(length(unique(data$mut_op)), length(unique(data$wt_ID)))),
                width = 1200) %>%  
    onRender(
      '
  function(el,x){
  // select all our node text
  d3.select(el)
  .selectAll(".node text")
  .filter(function(d) { return d.name.startsWith("<"); })
  .attr("x", x.options.nodeWidth - 90)
  .attr("text-anchor", "end");
  }
  '
    ) %>%
    saveNetwork(file = paste(i,'trans_reversion_mutations_sankey_plot.html'))
  
  # Get final data
  data <- sankey[sankey$strategy_rep==i & sankey$timepoint==1,]
  data <- cbind.data.frame(data, value=rep(1, length(data$strategy_rep)))
  data_nodes <- cbind.data.frame(name=c(unique(data$wt_ID), unique(data$mut_op)), empty=rep(NA,length(unique(data$wt_ID)) + length(unique(data$mut_op))))
  data$IDsource=match(data$mut_op, data_nodes$name)-1 
  data$IDtarget=match(data$wt_ID, data_nodes$name)-1 
  
  # prepare D3 colour scale for final timepoint
  op_colors <- c()
  
  for (j in data_nodes$name) {
    if (j %in% data$wt_ID) {
      temp <- unique(data$wt_ID[data$wt_ID==j])
      temp <- unique(data$wt_op[data$wt_ID==temp])
      op_colors[length(op_colors)+1] <- AVIDA_instructions$color[AVIDA_instructions$operation==temp]
    }
    else if (j %in% data$mut_op) {
      temp <- unique(data$mut_op[data$mut_op==j])
      op_colors[length(op_colors)+1] <- AVIDA_instructions$color[AVIDA_instructions$operation==temp]
    }
  }
  
  Color_scale_final = paste('d3.scaleOrdinal() .domain([', 
                            gsub("^c\\(|\\)$", "", toString(paste(list(data_nodes$name), sep=""))),
                            ']) .range([',
                            gsub("^c\\(|\\)$", "", toString(list(op_colors))),
                            '])', sep="")
  # Make the Network
  sankeyNetwork(Links = data, 
                Nodes = data_nodes,
                Source = "IDsource", 
                Target = "IDtarget",
                Value = "value", 
                NodeID = "name", 
                sinksRight=FALSE, 
                colourScale=Color_scale_final, 
                nodeWidth=80, 
                fontSize=18, 
                nodePadding=8, 
                iterations = 0, 
                margin = 200,
                height = 20*max(c(length(unique(data$mut_op)), length(unique(data$wt_ID)))),
                width = 1200) %>% 
    onRender(
      '
  function(el,x){
  // select all our node text
  d3.select(el)
  .selectAll(".node text")
  .filter(function(d) { return d.name.startsWith("<"); })
  .attr("x", x.options.nodeWidth - 90)
  .attr("text-anchor", "end");
  }
  '
    ) %>%
    saveNetwork(file = paste(i,'final_reversion_mutations_sankey_plot.html'))
  
}



###-------------------------------------------------------
# REVERSION MUTATION HISTOGRAMS

rev_mut_df$wt_op <- factor(rev_mut_df$wt_op, levels=c(AVIDA_instructions$operation))
rev_mut_df$mut_op <- factor(rev_mut_df$mut_op, levels=c(AVIDA_instructions$operation))
rev_mut_df$timepoint[rev_mut_df$timepoint==0] <- "Transition"
rev_mut_df$timepoint[rev_mut_df$timepoint==1] <- "Final"
rev_mut_df$timepoint <- factor(rev_mut_df$timepoint, levels=c("Transition", "Final"))
rev_mut_df_cc <- rev_mut_df[complete.cases(rev_mut_df),]



# Checking on the frequency of nopx operations
multi_nopx_freq_trans <- length(multi_genomes$strategy_rep[multi_genomes$wt_mut==3 & multi_genomes$timepoint==0]) / length(multi_genomes$strategy_rep[multi_genomes$timepoint==0])
multi_nopx_freq_final <- length(multi_genomes$strategy_rep[multi_genomes$wt_mut==3 & multi_genomes$timepoint==1]) / length(multi_genomes$strategy_rep[multi_genomes$timepoint==1])
uni_nopx_freq_trans <- length(rev_mut_df_cc$strategy_rep[rev_mut_df_cc$wt_mut==3 & rev_mut_df_cc$timepoint=="Transition"]) / length(rev_mut_df_cc$strategy_rep[rev_mut_df_cc$timepoint=="Transition"])
uni_nopx_freq_final <- length(rev_mut_df_cc$strategy_rep[rev_mut_df_cc$wt_mut==3 & rev_mut_df_cc$timepoint=="Final"]) / length(rev_mut_df_cc$strategy_rep[rev_mut_df_cc$timepoint=="Final"])

uni_nopx_freq_trans/multi_nopx_freq_trans
uni_nopx_freq_final/multi_nopx_freq_final

# Plot identity of wt_op at transition
wt_ops <- ggplot(data=rev_mut_df_cc, 
                 aes(x=factor(wt_op))) +
  geom_bar(stat="count", width=0.8, aes(fill=factor(timepoint))) +
  facet_wrap(~timepoint, ncol=1) +
  scale_fill_manual(values=c("goldenrod1", "green3"), name = "Timepoint") +
  xlab("Wild type operation") +
  ylab("Count") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


# Plot identity of mut_op at transition
mut_ops <- ggplot(data=rev_mut_df_cc, 
                  aes(x=factor(mut_op))) +
  geom_bar(stat="count", width=0.8, aes(fill=factor(timepoint))) +
  facet_wrap(~timepoint, ncol=1) +
  scale_fill_manual(values=c("goldenrod1", "green3"), name = "Timepoint") +
  xlab("Mutant operation") +
  ylab("Count") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Figures")

pdf("Mutant_operations_at_transition_and_final_time_point_20DEC22.pdf", width=one_col*3, height=one_col*4)
mut_ops
dev.off()

pdf("Wildtype_operations_at_transition_and_final_time_point_20DEC22.pdf", width=one_col*3, height=one_col*4)
wt_ops
dev.off()


# Plot identity of wt_op at transition
wt_ops <- ggplot(data=rev_mut_df_cc, 
                 aes(x=factor(wt_op))) +
  geom_bar(stat="count", width=0.8, aes(fill=factor(timepoint))) +
  facet_wrap(~timepoint, ncol=1) +
  scale_fill_manual(values=c("goldenrod1", "green3"), name = "Timepoint") +
  xlab("Wild type operation") +
  ylab("Count") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 


# Plot identity of mut_op at transition
mut_ops <- ggplot(data=rev_mut_df_cc, 
                  aes(x=factor(mut_op))) +
  geom_bar(stat="count", width=0.8, aes(fill=factor(timepoint))) +
  facet_wrap(~timepoint, ncol=1) +
  scale_fill_manual(values=c("goldenrod1", "green3"), name = "Timepoint") +
  xlab("Mutant operation") +
  ylab("Count") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position="none") 


pdf("Reversion_mutations_at_transition_and_final_time_point_20DEC22.pdf", width=one_col*6, height=one_col*4)
ggarrange(mut_ops, wt_ops, ncol=2, labels=c("A", "B"), label.args=list(gp = grid::gpar(font = 2, cex = 1.2)))
dev.off()

# 
# # Compare histograms of reversion mutations at transition and final timepoint 
# ggplot(data=rev_mut_df, 
#        aes(x=as.numeric(genome_location))) +
#   geom_histogram(data=rev_mut_df[rev_mut_df$timepoint=="Transition",], 
#                  stat="bin", binwidth=1, alpha=0.7, fill="goldenrod1") +
#   geom_histogram(data=rev_mut_df[rev_mut_df$timepoint=="Final",], 
#                  stat="bin", binwidth=1, alpha=0.7, fill="green3") +
#   scale_x_discrete(breaks=seq(from=0, to=100, by=1)) +
#   scale_fill_manual(values = c("dodgerblue", "goldenrod1")) +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   facet_wrap(~factor(strategy_rep))
# 
# # Plot histogram of reversion mutations by position w/ color indicating mutation identity for transition timepoint
# ggplot(data=rev_mut_df[rev_mut_df$timepoint=="Transition",], 
#        aes(x=as.numeric(genome_location), fill=reorder(mutation, as.numeric(mutation)))) +
#   geom_histogram(stat="bin", binwidth=1) +
#   scale_x_discrete(breaks=seq(from=0, to=100, by=1)) +
#   scale_fill_manual(values = viridis(38)) +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   facet_wrap(~factor(strategy_rep))
# 
# # Plot histogram of reversion mutations by position w/ color indicating mutation identity for transition timepoint
# ggplot(data=rev_mut_df[rev_mut_df$timepoint=="Final",],
#        aes(x=as.numeric(genome_location), fill=reorder(mutation, as.numeric(mutation)))) +
#   geom_histogram(stat="bin", binwidth=1) +
#   scale_x_discrete(breaks=seq(from=0, to=100, by=1)) +
#   scale_fill_manual(values = viridis(38)) +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   facet_wrap(~factor(strategy_rep))
# 
# 
# 
# 
# 
# ###-------------------------------------------------------
# # EXPLORING REVERSION SUMMARY DATA AND SUMMARY DIFFERENCE DATA
# 
# # Read in data
# all_rev_m_summary_data <- read.csv("all_reversion_mutations_summary_data.csv", header=T, stringsAsFactors=F)
# diff_df <- read.csv("all_reversion_mutations_difference_data.csv", header=T, stringsAsFactors=F)
# 
# ggplot(data=all_rev_m_summary_data, aes(y=num_rev_m, x=factor(timepoint))) +
#   geom_violin()
# 
# ggplot(data=diff_df, aes(y=num_rev_m_diff, x=factor(timepoint))) +
#   geom_violin()
# ggplot(data=diff_df, aes(y=num_32_diff, x=factor(timepoint))) +
#   geom_violin()
# ggplot(data=diff_df, aes(y=num_unique_m_diff, x=factor(timepoint))) +
#   geom_violin()
# ggplot(data=diff_df, aes(y=num_unique_p_diff, x=factor(timepoint))) +
#   geom_violin()
# ggplot(data=diff_df, aes(y=mut_div_diff, x=factor(timepoint))) +
#   geom_violin()
# ggplot(data=diff_df, aes(y=pos_div_diff, x=factor(timepoint))) +
#   geom_violin()
# 
# 
# ggpairs(data=diff_df[,c(3,5:7)])

head(multi_genomes)
strategy_rep <- c()
timepoint <- c()
num_tiss_acc <- c()

for (i in unique(multi_genomes$strategy_rep)) {
  for (j in unique(multi_genomes$timepoint)) {
    strategy_rep[length(strategy_rep)+1] <- i
    timepoint[length(timepoint)+1] <- j
    num_tiss_acc[length(num_tiss_acc)+1] <- length(multi_genomes$wt_mut[multi_genomes$strategy_rep==i & multi_genomes$timepoint==j & multi_genomes$wt_mut==24])
  }
}

tiss_acc_df <- cbind.data.frame(strategy_rep, timepoint, num_tiss_acc)

tiss_acc_df_1 <- tiss_acc_df[tiss_acc_df$timepoint==1,]
tiss_acc_df$strategy_rep <- factor(tiss_acc_df$strategy_rep, levels=(tiss_acc_df$strategy_rep[tiss_acc_df$timepoint==1])[rev(order(tiss_acc_df$num_tiss_acc[tiss_acc_df$timepoint==1]))], ordered=TRUE)

ggplot() +
  geom_histogram(data=tiss_acc_df[tiss_acc_df$timepoint==1,], aes(x=as.factor(strategy_rep), y=num_tiss_acc), fill="green3", stat="identity", position="identity") +
  geom_histogram(data=tiss_acc_df[tiss_acc_df$timepoint==0,], aes(x=as.factor(strategy_rep), y=num_tiss_acc), fill="darkgoldenrod1", stat="identity", position="identity", alpha=0.8) +
  ylab("Number of tissue accretion instructions") +
  xlab("Multicellular lineage") +
  theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))



reshape_tiss_acc_df <- cbind.data.frame(strategy_rep = unique(tiss_acc_df$strategy_rep), 
                 num_tiss_acc_0 = tiss_acc_df$num_tiss_acc[tiss_acc_df$timepoint==0], 
                 num_tiss_acc_1 = tiss_acc_df$num_tiss_acc[tiss_acc_df$timepoint==1],
                 strain_ID = c(outer(c("DD."), c(1:length(unique(strategy_rep))), FUN=paste, sep=""))
                 )

head(reshape_tiss_acc_df)
reshape_tiss_acc_df$strain_ID <- factor(reshape_tiss_acc_df$strain_ID, levels=(reshape_tiss_acc_df$strain_ID)[rev(order(tiss_acc_df$num_tiss_acc[tiss_acc_df$timepoint==1]))], ordered=TRUE)
reshape_tiss_acc_df$strategy_rep
reshape_tiss_acc_df$strain_ID

P1 <- ggplot() +
  geom_histogram(data=reshape_tiss_acc_df, aes(x=as.factor(strain_ID), y=num_tiss_acc_1), fill="green3", stat="identity", position="identity") +
  geom_histogram(data=reshape_tiss_acc_df, aes(x=as.factor(strain_ID), y=num_tiss_acc_0), fill="darkgoldenrod1", stat="identity", position="identity", alpha=0.8) +
  ylab("# Tissue accretion instructions") +
  xlab("Multicellular lineage") +
  theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))


setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/Figure_6/Figure_6_data/")
# entrenchment_data <- read.csv("full_entrenchment_data_filtered_2_14FEB22.csv")
# head(entrenchment_data)
Fig_6_data <- read.csv("Full_Figure_6_data_25OCT22.csv", header=T, stringsAsFactors=F)
Fig_6_data <- Fig_6_data[,-c(1,2)]
Fig_6_data$strategy_rep <- factor(Fig_6_data$strategy_rep, levels=(Fig_6_data$strategy_rep)[rev(order(tiss_acc_df$num_tiss_acc[tiss_acc_df$timepoint==1]))], ordered=TRUE)

head(Fig_6_data)

final_df <- left_join(reshape_tiss_acc_df, Fig_6_data, by="strategy_rep")

colors <- c("Transition" = "darkgoldenrod1", "Final" = "green3")

P2 <- ggplot() +
  geom_jitter(data=final_df, aes(y=adjusted_trans_viable_uni, x=num_tiss_acc_0, color="Transition"), alpha=0.6, width=0) +
  geom_jitter(data=final_df, aes(y=adjusted_final_viable_uni, x=num_tiss_acc_1, color="Final"), alpha=0.6, width=0) +
  xlab("# Tissue accretion instructions") +
  ylab("# Viable unicellular revertants") +
  labs(color = "Time point") +
  scale_color_manual(values = colors)

setwd("~/Desktop")

pdf("Redundant_tissue_accretion_operations_23DEC22.pdf", width=two_col*1.5, height=two_col*1.3)
ggdraw() +
  draw_plot(P1, x = 0.04, y = 0.60, width = 0.95, height = 0.4) +
  draw_plot(P2, x = 0.02, y = 0.0, width = 0.65, height = 0.6) +
  draw_plot_label(label = c("A","B"), size = 11, 
                  family = "Helvetica", fontface = "bold",
                  x = c(0, 0), 
                  y = c(1, 0.6))
dev.off()


ggarrange(P1, P2, labels = c("A", "B"), label.args = list(gp = grid::gpar(font = 2, cex = 1.2)), nrow=2)
