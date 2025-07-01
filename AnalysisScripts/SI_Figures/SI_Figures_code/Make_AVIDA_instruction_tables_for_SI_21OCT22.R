library("gridExtra")
library("grid")

TSpecial <- ttheme_minimal(
  core=list(bg_params = list(fill = "white", 
                             col="black"),
            fg_params=list(col="black",
                           hjust=0, 
                           x=00.025)),
  colhead=list(bg_params =list(fill = "grey80", 
                               col="black"),
               fg_params=list(col="black", 
                              hjust=0, 
                              x=00.025)),
  rowhead=list(fg_params=list(col=rgb(0,0,0,0), 
                              fontface=4)))


setwd("~/Desktop/Dropbox/2019_Goldsby_et_al_Entrenchment_MS/Entrenchment_processed_data_for_figures/SI_Figures")


no_op <- read.csv("Avida_instructions_No-operation.csv", header=T, stringsAsFactors=F)
for (i in 1:length(no_op$Description)) {
  no_op$Description[i] <- paste(strwrap(no_op$Description[i], width = 60), collapse = "\n")
}

colnames(no_op) <- c("ID #",
                     "Old Instruction                     ",
                     "Instruction                                     ",
                     "Category             ",
                     "Description                                                                               ")


grid.newpage()
pdf("Avida_instructions_No-operation_Table.pdf", width=8, height=3)
grid.table(d=no_op[,c(3,5)],
           theme=TSpecial)
dev.off()


math <- read.csv("Avida_instructions_Math.csv", header=T, stringsAsFactors=F)
for (i in 1:length(math$Description)) {
  math$Description[i] <- paste(strwrap(math$Description[i], width = 60), collapse = "\n")
}

colnames(math) <- c("ID #",
                    "Old Instruction                     ",
                    "Instruction                                     ",
                    "Category             ",
                    "Description                                                                               ")
grid.newpage()
pdf("Avida_instructions_Math_Table_2025.pdf", width=8, height=1.75)
grid.table(d=math[,c(3,5)],
           theme=TSpecial)
dev.off()



control <- read.csv("Avida_instructions_Hardware_control.csv", header=T, stringsAsFactors=F)
for (i in 1:length(control$Description)) {
  control$Description[i] <- paste(strwrap(control$Description[i], width = 60), collapse = "\n")
}

colnames(control) <- c("ID #",
                       "Old Instruction                     ",
                       "Instruction                                     ",
                       "Category             ",
                       "Description                                                                               ")

grid.newpage()
pdf("Avida_instructions_Hardware_control_Table_2025.pdf", width=8, height=8.5)
grid.table(d=control[,c(3,5)],
           theme=TSpecial)
dev.off()



interaction <- read.csv("Avida_instructions_Interaction.csv", header=T, stringsAsFactors=F)
for (i in 1:length(interaction$Description)) {
  interaction$Description[i] <- paste(strwrap(interaction$Description[i], width = 60), collapse = "\n")
}

colnames(interaction) <- c("ID #",
                           "Old Instruction                     ",
                           "Instruction                                     ",
                           "Category             ",
                           "Description                                                                               ")

grid.newpage()
pdf("Avida_instructions_Interaction_Table_2025.pdf", width=8, height=10)
grid.table(d=interaction[,c(3,5)],
           theme=TSpecial)
dev.off()




biological <- read.csv("Avida_instructions_Biological.csv", header=T, stringsAsFactors=F)
for (i in 1:length(biological$Description)) {
  biological$Description[i] <- paste(strwrap(biological$Description[i], width = 60), collapse = "\n")
}

colnames(biological) <- c("ID #",
                          "Old Instruction                     ",
                          "Instruction                                     ",
                          "Category             ",
                          "Description                                                                               ")

grid.newpage()
pdf("Avida_instructions_Biological_Table.pdf", width=8, height=7.75)
grid.table(d=biological[,c(3,5)],
           theme=TSpecial)
dev.off()



