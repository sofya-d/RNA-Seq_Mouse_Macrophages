# setup
library(tidyverse)
library(ggplot2)
library(grid)


## create empty data frame with genes associated with macrophage activation
# create pre-data frame
m1de <- as.data.frame(c("replicate1", "replicate2"))

# create vector of genes associated with macrophage activation
m1genes <- c("Nos2", "Tnf", "Cxcl9", "Cxcl10", "Cxcl11", "Stat1", "Irf1", "Irf8")

# add columns with genes to data frame "m1de"
for (gene in m1genes){
  m1de[[gene]] <- NA
}

# turn first column into row names
rownames(m1de) <- m1de$replicates
m1de <- m1de[, 2:9]


## create list of data frames filled with normalized counts
# create list of empty data frames
m1detpm <- list(Ctr = m1de, IFN_3 = m1de, CPG_3 = m1de, IFN_CPG_3 = m1de, IFN_9 = m1de, CPG_9 = m1de, IFN_CPG_9 = m1de, IFN_15 = m1de, CPG_15 = m1de, IFN_CPG_15 = m1de)

# create vector of names of "m1detpm" objects considering biological replicates
m1detpm_names <- c()

for (i in names(m1detpm)){
  m1detpm_names <- c(m1detpm_names, rep(i, 2))
}

# fill normalized counts in data frames of "m1detpm" list
replicate = 2
colnum = 1
for (name in m1detpm_names){
  replicate = replicate + ((-1) ^ (replicate - 1))
  colnum = colnum + 1
  for (i in 1:length(norm_counts$GeneName)){
    for (j in 1:8){
      if (norm_counts$GeneName[i] == colnames(m1detpm[[name]])[j]){
        m1detpm[[name]][replicate, j] <- norm_counts[i, colnum]
      }
    }
  }
}

# add "sample" column to the data frames in list "m1detpm"
for (name in names(m1detpm)){
  m1detpm[[name]]$sample <- c(name, name)
}

## prepare table for barplot
# concatenate data frames of list "m1detpm"
m1de <- m1detpm[["Ctr"]]
for (name in names(m1detpm)[-1]){
  m1de <- rbind(m1de, m1detpm[[name]])
}

# create data frame grouped by sample and gene
m1g <- gather(m1de, gene, TPM, 1:8)
m1g_grouped <- group_by(m1g, sample, gene)

# create data frame grouped by sample, concatenate replicate columns with mean TPM in cell and add column with standard deviation values
m1g_sumrsd <- summarise(m1g_grouped, mean = mean(TPM), sd = sd(TPM))

# reorder summarized table
m1g_ordrd <- m1g_sumrsd[c(25:32, 41:48, 9:16, 65:72, 49:56, 17:24, 73:80, 33:40, 1:8, 57:64),]
# fix the order
m1g_ordrd$sample <- factor(m1g_ordrd$sample, levels = unique(m1g_ordrd$sample), ordered=TRUE)


#SEPARATE PLOTS
# list of data frames for plots separated by time
m1_plot <- list(hour_3 = m1g_ordrd[1:32,], hour_9 = m1g_ordrd[c(1:8, 33:56),], hour_15 = m1g_ordrd[c(1:8, 57:80),])

# obtain genes that have significant differential expression compared to control
m1express <- list()
for (name in names(lrts_de)){
  vec <- c()
  for (gene in colnames(m1detpm[["Ctr"]])){
    if (gene %in% lrts_de[[name]]$GeneName){
      vec <- c(vec, gene)
m1express[[name]] <- vec
    }
  }
}

# create list of vectors of characters that indicate if gene has significant differential expression. those vectors will be used to mark those genes on barplot with "*" symbol
start <- -2
end <- 0
star <- c()
stars <- list()
for (name in names(m1_plot)){
  start <- start + 3
  end <- end + 3
  for (i in start:end){
    for (gene in m1_plot[["hour_3"]]$gene[1:8]){
      
      if (gene %in% m1express[[i]]){
        star <- c(star, "*")
      } else star <- c(star, " ")
    }
  }
  stars[[name]] <- c(rep(" ", 8), star)
  star <- c()
}


## draw bar plot
pdf("M1_DEG.pdf", width = 8, height = 5)
#PLOT 3h
plot <- ggplot(m1_plot[["hour_3"]], aes(x = sample, y = mean, fill = gene))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(0.9), size = 0.6, width = 0.25)+
  xlab("Samples")+
  ylab("TPM")+
  ggtitle("M1 DEG 3h")+
  geom_text(aes(label = stars[[1]]), position = position_dodge(0.8), size = 6, colour = "red")

grid.draw(plot)

#PLOT 9h
plot <- ggplot(m1_plot[["hour_9"]], aes(x = sample, y = mean, fill = gene))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(0.9), size = 0.6, width = 0.25)+
  xlab("Samples")+
  ylab("TPM")+
  ggtitle("M1 DEG 9h")+
  geom_text(aes(label = stars[[2]]), position = position_dodge(0.8), size = 6, colour = "red")

grid.draw(plot)

#PLOT 15h
plot <- ggplot(m1_plot[["hour_15"]], aes(x = sample, y = mean, fill = gene))+
  geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(0.9), size = 0.6, width = 0.25)+
  xlab("Samples")+
  ylab("TPM")+
  ggtitle("M1 DEG 15h")+
  geom_text(aes(label = stars[[3]]), position = position_dodge(0.8), size = 6, colour = "red")

grid.draw(plot)

dev.off()
