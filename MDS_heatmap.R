#SETUP
library(edgeR)
library(heatmap3)

#create design matrix
group <- as.factor(c('ctr', 'ctr', 'h3_i', 'h3_i', 'h3_c', 'h3_c', 'h3_i.c',
                     'h3_i.c', 'h9_i', 'h9_i', 'h9_c', 'h9_c', 'h9_i.c',
                     'h9_i.c', 'h15_i', 'h15_i', 'h15_c', 'h15_c', 'h15_i.c',
                     'h15_i.c'))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

#make contrasts
cont <- list()
#CONTRAST 3h
cont[["IFNvsCtrl_3h"]] <- makeContrasts(h3_i - ctr, levels = design)
cont[["CPGvsCtrl_3h"]] <- makeContrasts(h3_c - ctr, levels = design)
cont[["IFN.CPGvsCtrl_3h"]] <- makeContrasts(h3_i.c - ctr, levels = design)

#CONTRAST 9h
cont[["IFNvsCtrl_9h"]] <- makeContrasts(h9_i - ctr, levels = design)
cont[["CPGvsCtrl_9h"]] <- makeContrasts(h9_c - ctr, levels = design)
cont[["IFN.CPGvsCtrl_9h"]] <- makeContrasts(h9_i.c - ctr, levels = design)

#CONTRAST 15h
cont[["IFNvsCtrl_15h"]] <- makeContrasts(h15_i - ctr, levels = design)
cont[["CPGvsCtrl_15h"]] <- makeContrasts(h15_c - ctr, levels = design)
cont[["IFN.CPGvsCtrl_15h"]] <- makeContrasts(h15_i.c - ctr, levels = design)

# create table with normalized counts only
# level2_mm39 - table with row and normalized counts
norm_counts <- mm10_level2[, 27:46]
norm_counts <- norm_counts[, c(c(11, 13, 14:19, 20, 1, 2:5, 6:10, 12))]
rownames(norm_counts) <- mm10_level2$id

colnames(norm_counts) <- c('ctr1', 'ctr2', 'h3_i1', 'h3_i2', 'h3_c1', 'h3_c2', 'h3_i.c1', 'h3_i.c2', 'h9_i1', 'h9_i2', 'h9_c1', 'h9_c2', 'h9_i.c1', 'h9_i.c2', 'h15_i1', 'h15_i2', 'h15_c1', 'h15_c2', 'h15_i.c1', 'h15_i.c2')


# bind gene names to cluster ids
# an - level2_mm39_annotated, contains gene names
# checked the IDs of genes Tmtc3 and Zfp935
norm_counts$GeneName <- NA
norm_counts <- norm_counts[, c(21, 1:20)]

for (i in 1:nrow(norm_counts)){
  for (j in 1:nrow(mm10_level2_annotated)){
    if (rownames(norm_counts)[i] == mm10_level2_annotated$CLUSTER_ID[j]){
      norm_counts$GeneName[i] <- mm10_level2_annotated$GeneName[j]
    }
  }
}


# create DGEList object
group <- as.factor(c('ctr', 'ctr', 'h3_i', 'h3_i', 'h3_c', 'h3_c', 'h3_i.c',
                     'h3_i.c', 'h9_i', 'h9_i', 'h9_c', 'h9_c', 'h9_i.c',
                     'h9_i.c', 'h15_i', 'h15_i', 'h15_c', 'h15_c', 'h15_i.c',
                     'h15_i.c'))

norm_counts_dge <- DGEList(as.matrix(norm_counts[, 2:21]), group = group, 
                           genes = data.frame(GeneName = norm_counts$GeneName))

# for LRT
norm_counts_dge <- calcNormFactors(norm_counts_dge, method = "TMM")
norm_counts_dge <- estimateGLMTrendedDisp(norm_counts_dge, design)
norm_counts_dge <- estimateGLMTagwiseDisp(norm_counts_dge, design)

fit <- glmFit(norm_counts_dge, design)

# create list of LRTs
lrts <- list()
for (name in names(cont)){
  lrt <- glmLRT(fit, contrast = cont[[name]])
  lrt <- topTags(lrt, n=nrow(lrt$table))
  lrts[[name]] <- lrt$table
}


# subset DE genes
lrts_de <- list()
for (name in names(lrts)){
  lrts_de[[name]] <- subset(lrts[[name]], lrts[[name]]$FDR < 0.05 & abs(lrts[[name]]$logFC) > 1)
}

# rename core tables
counts <- mm10_level2
annotation <- mm10_level2_annotated

# create MDS plot
# changed size of the "Plots" window of RStudio and exported plot as .pdf (it looked well only in portrait orientation in some reason)
group <- factor(colnames(norm_counts)[2:21], levels = (colnames(norm_counts)[2:21]))
points <- c(5,18,1,16,2,17,0,15,1,16,2,17,0,15,1,16,2,17,0,15)
colors <- c(1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4)
plotMDS(norm_counts_dge, col=colors[group], pch=points[group], cex.lab=0.75, cex.axis=0.75, main='MDS')
par(mar=c(10.2, 4.1, 3.1, 4.1), xpd=TRUE)
legend("bottom", legend=group, pch=points, col=colors, ncol=3, 
       xpd=TRUE, inset=c(0, -.8), cex=.6)

# create a heatmap
library(heatmap3)

heatmap3(as.matrix(norm_counts[, 2:21]), main = "Normalized counts (TPM)", scale = "none", Rowv = NA, Colv = NA)

scaled_norm_counts <- norm_counts[, 2:21]
scaled_norm_counts <- abs(log2(scaled_norm_counts))

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

scaled_norm_counts <- t(apply(scaled_norm_counts, 1, cal_z_score))

# create a scaled heat map
heatmap3(scaled_norm_counts, main = "Normalized counts (|log2(TPM)| z-scores)", scale = "none", Rowv = NA)
