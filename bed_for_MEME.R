control <- c()

for (name in names(lrts)){
  control <- c(control, subset(rownames(lrts[[name]]), lrts[[name]]$FDR > 0.05))
}

control <- unique(control)

experiment <- c()

for (name in names(lrts_de)){
  experiment <- c(experiment, rownames(lrts_de[[name]]))
}

experiment <- unique(experiment)

cluster_ids <- list(control = control, experiment = experiment)

for (vec in names(cluster_ids)){
  bed <- as.data.frame(cluster_ids[[vec]])
  colnames(bed) <- "name"
  bed$chrom <- NA
  bed$chromStart <- NA
  bed$chromEnd <- NA
  
  for (i in 1:length(counts$id)){
    for (j in 1:length(bed$name)){
      if (bed$name[j] == counts$id[i]){
        bed$chromStart[j] <- counts$pos[i] - 200
        bed$chromEnd[j] <- counts$pos[i] + 200
        bed$chrom[j] <- counts$chrom[i]
      }
    }
  }
  
  bed <- bed[, c(2:4, 1)]
  
  write.table(bed, file = paste(vec, ".bed", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}


######### Top 200 #########
exp_sorted <- sort.int(lrts_de[["IFN.CPGvsCtrl_15h"]]$logFC, index.return = T, decreasing = T)

top_200_ind <- exp_sorted[["ix"]][1:200]

top200_lrt <- rownames(lrts_de[["IFN.CPGvsCtrl_15h"]][top_200_ind, ])

top200 <- subset(counts, counts$id %in% top200_lrt)
top200 <- top200[, c(2, 1, 5, 6)]
top200$chromStart <- top200$pos - 200
top200$chromEnd <- top200$pos + 200
top200$score <- 0
top200 <- top200[, c(1, 5, 6, 2, 7, 3)]
colnames(top200)[4] <- "name"

other <- subset(counts, !(counts$id %in% top200_lrt))
other <- other[, c(2, 1, 5, 6)]
other$chromStart <- other$pos - 200
other$chromEnd <- other$pos + 200
other$score <- 0
other <- other[, c(1, 5, 6, 2, 7, 3)]
colnames(other)[4] <- "name"
  
write.table(top200, file = "top200.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(other, file = "other.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


## experiment separated by different matter
names <- c("exp_Ifn", "exp_CpG", "exp_Ifn.CpG")
cluster_ids <- list()
ind <- c(1, 4, 7)
experiment <- c()

for (num in c(0, 1, 2)){
  ind <- ind + num
  for (name in names(lrts_de[ind])){
    experiment <- c(experiment, rownames(lrts_de[[name]]))
  }
  experiment <- unique(experiment)
  cluster_ids[[names[num + 1]]] <- experiment
  experiment <- c()
}

for (vec in names(cluster_ids)){
  bed <- as.data.frame(cluster_ids[[vec]])
  colnames(bed) <- "name"
  bed$chrom <- NA
  bed$chromStart <- NA
  bed$chromEnd <- NA
  
  for (i in 1:length(counts$id)){
    for (j in 1:length(bed$name)){
      if (bed$name[j] == counts$id[i]){
        bed$chromStart[j] <- counts$pos[i] - 200
        bed$chromEnd[j] <- counts$pos[i] + 200
        bed$chrom[j] <- counts$chrom[i]
      }
    }
  }
  
  bed <- bed[, c(2:4, 1)]
  
  write.table(bed, file = paste(vec, ".bed", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}
