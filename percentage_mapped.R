table <- flagstat_table

table <- table[c(11, 13:20, 1:10, 12),]

samples <- rep(c("Ctrl", "IFN_3h", "CpG_3h", "IFN.CpG_3h", "IFN_9h", "CpG_9h", "IFN.CpG_9h", "IFN_15h", "CpG_15h", "IFN.CpG_15h"), each = 2)

replicate <- 2
for (i in 1:length(samples)){
  replicate <- replicate + ((-1) ^ (replicate - 1))
  samples[i] <- paste(samples[i], ".", replicate, sep = "")
}

table$sample <- samples
rownames(table) <- 1:length(table$sample)

table <- table[, c(2, 3, 1)]

table <- gather(table, condition, percent, 1:2)

table <- group_by(table, sample)
table <-  arrange(table)

table$sample <- factor(table$sample, levels = unique(table$sample), ordered = T)

library(ggplot2)

ggplot(table, aes(fill = condition, y = percent, x = sample)) + 
  geom_bar(position = "stack", stat = "identity") +
  ggtitle("Percentage of mapped reads") +
  guides(x =  guide_axis(angle = 45))

