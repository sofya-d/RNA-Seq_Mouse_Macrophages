#PACKAGES
library(factoextra)
library(RColorBrewer)
library(gplots)


# DATA PREPARATION
# vector of DE genes
de_genes <- c()

for (name in names(lrts_de)){
  de_genes <- c(de_genes, rownames(lrts_de[[name]]))
}

de_genes <- unique(de_genes)

# dataframe
de_heatmap <- as.data.frame(de_genes)
rownames(de_heatmap) <- de_genes

de_heatmap$IFN_3h <- NA
de_heatmap$CPG_3h <- NA
de_heatmap$IFNplusCPG_3h <- NA

de_heatmap$IFN_9h <- NA
de_heatmap$CPG_9h <- NA
de_heatmap$IFNplusCPG_9h <- NA

de_heatmap$IFN_15h <- NA
de_heatmap$CPG_15h <- NA
de_heatmap$IFNplusCPG_15h <- NA

de_heatmap <- de_heatmap[, 2:10]

for (i in 1:length(lrts)){
  for (j in 1:nrow(de_heatmap)){
    for (l in 1:nrow(lrts[[i]])){
      if (row.names(de_heatmap)[j] == row.names(lrts[[i]])[l]){
        de_heatmap[j, i] <- (log2(lrts[[i]]$FDR[l]) - mean(log2(lrts[[i]]$FDR))) / sd(log2(lrts[[i]]$FDR))
        break
      }
    }
  }
}


# calculate optimal number of clusters
# ELBOW METHOD (result: ~3)
fviz_nbclust(de_heatmap, kmeans, method = "wss")+
  labs(subtitle = "Elbow method")

# SILHOUETTE METHOD (result: ~2)
fviz_nbclust(de_heatmap, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

# GAP STATISTIC (50 ïðåäóïðåæäåíèé "íå ñîøëîñü çà 10 èòåðàöèé", result: 2)
set.seed(2022)
fviz_nbclust(de_heatmap, kmeans, nstart = 25,  method = "gap_stat", nboot = 500)+
  labs(subtitle = "Gap statistic method")


# K-MEANS
kmeans.res <- kmeans(as.matrix(de_heatmap), 2, nstart = 25)


#HEATMAP
kmeans.DEgenes.clusters <- as.vector(kmeans.res[["cluster"]])

palette  <- c(rep(c("#00D5FF"), 500), colorRampPalette(c("#00D5FF", "#FF00BF"))(n = 100), rep(c("#FF00BF"), 500))

heatmap.2(as.matrix(de_heatmap), main = "DE genes z-score(log2(FDR))", margins=c(9, 4), lhei=c(4, 15), cexCol=0.8, trace = 'none', tracecol = 'white', lwid = c(2, 7), Rowv = kmeans.DEgenes.clusters, col = palette, key.xlab = "z-score(log2(FDR))")


#SAVING DATAFRAME
de_heatmap$id <- rownames(de_heatmap)
de_heatmap$cluster <- kmeans.DEgenes.clusters
de_heatmap <- de_heatmap[, c(10, 11, 1:9)]
de_heatmap <- de_heatmap[, c(1, 11, 2:10)]

write.table(de_heatmap, file = "DE_heatmap.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


##### Less Genes Heatmap
# DATA PREPARATION
# vector of DE genes
less_genes <- c()

for (name in names(lrts_de)){
  less_genes <- c(less_genes, rownames(subset(lrts_de[[name]], abs(lrts_de[[name]]$logFC) > 6)))
}

less_genes <- unique(less_genes)

# dataframe
less_heatmap <- as.data.frame(less_genes)
rownames(less_heatmap) <- less_genes

less_heatmap$IFN_3h <- NA
less_heatmap$CPG_3h <- NA
less_heatmap$IFNplusCPG_3h <- NA

less_heatmap$IFN_9h <- NA
less_heatmap$CPG_9h <- NA
less_heatmap$IFNplusCPG_9h <- NA

less_heatmap$IFN_15h <- NA
less_heatmap$CPG_15h <- NA
less_heatmap$IFNplusCPG_15h <- NA

less_heatmap <- less_heatmap[, 2:10]

# fill values
for (i in 1:length(lrts)){
  for (j in 1:nrow(less_heatmap)){
    for (l in 1:nrow(lrts[[i]])){
      if (row.names(less_heatmap)[j] == row.names(lrts[[i]])[l]){
        less_heatmap[j, i] <- (log2(lrts[[i]]$FDR[l]) - mean(log2(lrts[[i]]$FDR))) / sd(log2(lrts[[i]]$FDR))
        break
      }
    }
  }
}

# bind gene names
less_heatmap$GeneName <- NA

for (i in 1:length(rownames(less_heatmap))){
  for (j in 1:length(rownames(norm_counts))){
    if (rownames(less_heatmap)[i] == rownames(norm_counts)[j]){
      if (norm_counts$GeneName[j] != ""){
        less_heatmap$GeneName[i] <- norm_counts$GeneName[j]
      } else {
        less_heatmap$GeneName[i] <- rownames(norm_counts)[j]
      }
    break
    }
  }
}

less_heatmap <- less_heatmap[, c(10, 1:9)]

# change duplicated gene names
# logical vector of duplicated values
dupl <- duplicated(less_heatmap$GeneName)

# indices of duplicated names
ind <- seq_along(less_heatmap$GeneName)[dupl]

# change duplicated names
for (i in ind){
  less_heatmap$GeneName[i] <- paste(less_heatmap$GeneName[i], ".2", sep = "")
}

# matrix for heatmap
less_heatmapmx <- as.matrix(less_heatmap[, 2:10])
rownames(less_heatmapmx) <- less_heatmap$GeneName


# calculate optimal number of clusters
# ELBOW METHOD (result: 3)
fviz_nbclust(less_heatmapmx, kmeans, method = "wss")+
  labs(subtitle = "less Elbow method")

# SILHOUETTE METHOD (result: 4)
fviz_nbclust(less_heatmapmx, kmeans, method = "silhouette")+
  labs(subtitle = "less Silhouette method")

# GAP STATISTIC (result: 10)
fviz_nbclust(less_heatmapmx, kmeans, nstart = 25,  method = "gap_stat", nboot = 500)+
  labs(subtitle = "less Gap statistic method")


# K-MEANS
less_kmeans.res <- kmeans(less_heatmapmx, centers = 6, nstart = 25)


#HEATMAP
less_kmeans.DEgenes.clusters <- as.vector(less_kmeans.res[["cluster"]])

heatmap.2(less_heatmapmx, main = "DE genes z-score(log2(FDR)) (|logFC| > 6)", margins=c(9, 4), lhei=c(4, 15), cexCol=0.8, trace = 'none', tracecol = 'white', lwid = c(2, 7), Rowv = less_kmeans.DEgenes.clusters, key.xlab = "z-score(log2(FDR))")

#SAVING DATAFRAME
less_heatmap$cluster <- less_kmeans.DEgenes.clusters
less_heatmap <- less_heatmap[, c(1, 11, 2:10)]
write.table(less_heatmap, file = "less_heatmap.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#### less genes heatmap 3h
##### Less Genes Heatmap
# DATA PREPARATION
# vector of DE genes
less3_genes <- c()

for (name in names(lrts_de)[1:3]){
  less3_genes <- c(less3_genes, rownames(subset(lrts_de[[name]], abs(lrts_de[[name]]$logFC) > 6)))
}

less3_genes <- unique(less3_genes)

# dataframe
less3_heatmap <- as.data.frame(less3_genes)
rownames(less3_heatmap) <- less3_genes

less3_heatmap$IFN_3h <- NA
less3_heatmap$CPG_3h <- NA
less3_heatmap$IFNplusCPG_3h <- NA

colnames(less3_heatmap)[1] <- "id"

less3_heatmap <- less3_heatmap[2:4]

# fill values
for (i in 1:3){
  for (j in 1:nrow(less3_heatmap)){
    for (l in 1:nrow(lrts[[i]])){
      if (row.names(less3_heatmap)[j] == row.names(lrts[[i]])[l]){
        less3_heatmap[j, i] <- (log2(lrts[[i]]$FDR[l]) - mean(log2(lrts[[i]]$FDR))) / sd(log2(lrts[[i]]$FDR))
        break
      }
    }
  }
}

# bind gene names
less3_heatmap$GeneName <- NA

for (i in 1:length(rownames(less3_heatmap))){
  for (j in 1:length(rownames(norm_counts))){
    if (rownames(less3_heatmap)[i] == rownames(norm_counts)[j]){
      if (norm_counts$GeneName[j] != ""){
        less3_heatmap$GeneName[i] <- norm_counts$GeneName[j]
      } else {
        less3_heatmap$GeneName[i] <- rownames(norm_counts)[j]
      }
      break
    }
  }
}

less3_heatmap <- less3_heatmap[, c(4, 1:3)]

# change duplicated gene names
# logical vector of duplicated values
dupl <- duplicated(less3_heatmap$GeneName)

# indices of duplicated names
ind <- seq_along(less3_heatmap$GeneName)[dupl]

# change duplicated names
for (i in ind){
  less3_heatmap$GeneName[i] <- paste(less3_heatmap$GeneName[i], ".2", sep = "")
}

# matrix for heatmap
less3_heatmapmx <- as.matrix(less3_heatmap[, 2:4])
rownames(less3_heatmapmx) <- less3_heatmap$GeneName


# calculate optimal number of clusters
# ELBOW METHOD (result: 3)
fviz_nbclust(less3_heatmapmx, kmeans, method = "wss")+
  labs(subtitle = "less 3h Elbow method")

# SILHOUETTE METHOD (result: 3)
fviz_nbclust(less3_heatmapmx, kmeans, method = "silhouette")+
  labs(subtitle = "less 3h Silhouette method")

# GAP STATISTIC (result: 6)
fviz_nbclust(less3_heatmapmx, kmeans, nstart = 25,  method = "gap_stat", nboot = 500)+
  labs(subtitle = "less 3h Gap statistic method")


# K-MEANS
less3_kmeans.res <- kmeans(less3_heatmapmx, centers = 4, nstart = 25)


#HEATMAP
less3_kmeans.DEgenes.clusters <- as.vector(less3_kmeans.res[["cluster"]])

heatmap.2(less3_heatmapmx, main = "DE genes z-score(log2(FDR)) (|logFC| > 6) 3h", margins=c(9, 4), lhei=c(4, 15), cexCol=0.8, trace = 'none', tracecol = 'white', lwid = c(2, 7), Rowv = less3_kmeans.DEgenes.clusters, key.xlab = "z-score(log2(FDR))")

#SAVING DATAFRAME
less3_heatmap$cluster <- less3_kmeans.DEgenes.clusters
less3_heatmap <- less3_heatmap[, c(1, 5, 2:4)]
write.table(less3_heatmap, file = "less3_DE_FDR.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#### less genes heatmap 9h
##### Less Genes Heatmap
# DATA PREPARATION
# vector of DE genes
less9_genes <- c()

for (name in names(lrts_de)[4:6]){
  less9_genes <- c(less9_genes, rownames(subset(lrts_de[[name]], abs(lrts_de[[name]]$logFC) > 6)))
}

less9_genes <- unique(less9_genes)

# dataframe
less9_heatmap <- as.data.frame(less9_genes)
rownames(less9_heatmap) <- less9_genes

less9_heatmap$IFN_9h <- NA
less9_heatmap$CPG_9h <- NA
less9_heatmap$IFNplusCPG_9h <- NA

colnames(less9_heatmap)[1] <- "id"

less9_heatmap <- less9_heatmap[2:4]

# fill values
for (i in 4:6){
  for (j in 1:nrow(less9_heatmap)){
    for (l in 1:nrow(lrts[[i]])){
      if (row.names(less9_heatmap)[j] == row.names(lrts[[i]])[l]){
        less9_heatmap[j, i - 3] <- (log2(lrts[[i]]$FDR[l]) - mean(log2(lrts[[i]]$FDR))) / sd(log2(lrts[[i]]$FDR))
        break
      }
    }
  }
}

# bind gene names
less9_heatmap$GeneName <- NA

for (i in 1:length(rownames(less9_heatmap))){
  for (j in 1:length(rownames(norm_counts))){
    if (rownames(less9_heatmap)[i] == rownames(norm_counts)[j]){
      if (norm_counts$GeneName[j] != ""){
        less9_heatmap$GeneName[i] <- norm_counts$GeneName[j]
      } else {
        less9_heatmap$GeneName[i] <- rownames(norm_counts)[j]
      }
      break
    }
  }
}

less9_heatmap <- less9_heatmap[, c(4, 1:3)]

# change duplicated gene names
# logical vector of duplicated values
dupl <- duplicated(less9_heatmap$GeneName)

# indices of duplicated names
ind <- seq_along(less9_heatmap$GeneName)[dupl]

# change duplicated names
for (i in ind){
  less9_heatmap$GeneName[i] <- paste(less9_heatmap$GeneName[i], ".2", sep = "")
}

# matrix for heatmap
less9_heatmapmx <- as.matrix(less9_heatmap[, 2:4])
rownames(less9_heatmapmx) <- less9_heatmap$GeneName


# calculate optimal number of clusters
# ELBOW METHOD (result: 3)
fviz_nbclust(less9_heatmapmx, kmeans, method = "wss")+
  labs(subtitle = "less 9h Elbow method")

# SILHOUETTE METHOD (result: 3)
fviz_nbclust(less9_heatmapmx, kmeans, method = "silhouette")+
  labs(subtitle = "less 9h Silhouette method")

# GAP STATISTIC (result: 10)
fviz_nbclust(less9_heatmapmx, kmeans, nstart = 25,  method = "gap_stat", nboot = 500)+
  labs(subtitle = "less 9h Gap statistic method")


# K-MEANS
less9_kmeans.res <- kmeans(less9_heatmapmx, centers = 6, nstart = 25)


#HEATMAP
less9_kmeans.DEgenes.clusters <- as.vector(less9_kmeans.res[["cluster"]])

heatmap.2(less9_heatmapmx, main = "DE genes z-score(log2(FDR)) (|logFC| > 6) 9h", margins=c(9, 4), lhei=c(4, 15), cexCol=0.8, trace = 'none', tracecol = 'white', lwid = c(2, 7), Rowv = less9_kmeans.DEgenes.clusters, key.xlab = "z-score(log2(FDR))")

#SAVING DATAFRAME
less9_heatmap$cluster <- less9_kmeans.DEgenes.clusters
less9_heatmap <- less9_heatmap[, c(1, 5, 2:4)]
write.table(less9_heatmap, file = "less9_DE_FDR.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#### less genes heatmap 15h
##### Less Genes Heatmap
# DATA PREPARATION
# vector of DE genes
less15_genes <- c()

for (name in names(lrts_de)[7:9]){
  less15_genes <- c(less15_genes, rownames(subset(lrts_de[[name]], abs(lrts_de[[name]]$logFC) > 6)))
}

less15_genes <- unique(less15_genes)

# dataframe
less15_heatmap <- as.data.frame(less15_genes)
rownames(less15_heatmap) <- less15_genes

less15_heatmap$IFN_15h <- NA
less15_heatmap$CPG_15h <- NA
less15_heatmap$IFNplusCPG_15h <- NA

colnames(less15_heatmap)[1] <- "id"

less15_heatmap <- less15_heatmap[2:4]

# fill values
for (i in 7:9){
  for (j in 1:nrow(less15_heatmap)){
    for (l in 1:nrow(lrts[[i]])){
      if (row.names(less15_heatmap)[j] == row.names(lrts[[i]])[l]){
        less15_heatmap[j, i - 6] <- (log2(lrts[[i]]$FDR[l]) - mean(log2(lrts[[i]]$FDR))) / sd(log2(lrts[[i]]$FDR))
        break
      }
    }
  }
}

# bind gene names
less15_heatmap$GeneName <- NA

for (i in 1:length(rownames(less15_heatmap))){
  for (j in 1:length(rownames(norm_counts))){
    if (rownames(less15_heatmap)[i] == rownames(norm_counts)[j]){
      if (norm_counts$GeneName[j] != ""){
        less15_heatmap$GeneName[i] <- norm_counts$GeneName[j]
      } else {
        less15_heatmap$GeneName[i] <- rownames(norm_counts)[j]
      }
      break
    }
  }
}

less15_heatmap <- less15_heatmap[, c(4, 1:3)]

# change duplicated gene names
# logical vector of duplicated values
dupl <- duplicated(less15_heatmap$GeneName)

# indices of duplicated names
ind <- seq_along(less15_heatmap$GeneName)[dupl]

# change duplicated names
for (i in ind){
  less15_heatmap$GeneName[i] <- paste(less15_heatmap$GeneName[i], ".2", sep = "")
}

# matrix for heatmap
less15_heatmapmx <- as.matrix(less15_heatmap[, 2:4])
rownames(less15_heatmapmx) <- less15_heatmap$GeneName


# calculate optimal number of clusters
# ELBOW METHOD (result: 3)
fviz_nbclust(less15_heatmapmx, kmeans, method = "wss")+
  labs(subtitle = "less 15h Elbow method")

# SILHOUETTE METHOD (result: 3)
fviz_nbclust(less15_heatmapmx, kmeans, method = "silhouette")+
  labs(subtitle = "less 15h Silhouette method")

# GAP STATISTIC (result: 4)
fviz_nbclust(less15_heatmapmx, kmeans, nstart = 25,  method = "gap_stat", nboot = 500)+
  labs(subtitle = "less 15h Gap statistic method")


# K-MEANS
less15_kmeans.res <- kmeans(less15_heatmapmx, centers = 3, nstart = 25)


#HEATMAP
less15_kmeans.DEgenes.clusters <- as.vector(less15_kmeans.res[["cluster"]])

heatmap.2(less15_heatmapmx, main = "DE genes z-score(log2(FDR)) (|logFC| > 6) 15h", margins=c(9, 4), lhei=c(4, 15), cexCol=0.8, trace = 'none', tracecol = 'white', lwid = c(2, 7), Rowv = less15_kmeans.DEgenes.clusters, key.xlab = "z-score(log2(FDR))")

#SAVING DATAFRAME
less15_heatmap$cluster <- less15_kmeans.DEgenes.clusters
less15_heatmap <- less15_heatmap[, c(1, 5, 2:4)]
write.table(less15_heatmap, file = "less15_DE_FDR.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
