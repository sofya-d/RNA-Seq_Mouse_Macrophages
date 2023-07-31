## set up
library(chromVAR)
library(BSgenome)
library(motifmatchr)
library(SummarizedExperiment)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BiocParallel)
library(ggplot2)


## list available genomes. there was no mm10 genome
available.genomes()

## Obtain vector of DE genes IDs
de_genes <- c()
for (name in names(lrts_de)){
  for (rowname in rownames(lrts_de[[name]])){
    de_genes <- c(de_genes, rowname)
  }
}

de_genes <- unique(de_genes)


## Get counts of DE genes
# preparing data frame
norm_counts_de <- as.data.frame(annotation$CLUSTER_ID)
colnames(norm_counts_de) <- "id"
norm_counts_de$chr <- annotation$CHROMOSOME
norm_counts_de$start <- annotation$TSS - 200
norm_counts_de$end <- annotation$TSS + 200
norm_counts_de$strand <- annotation$STRAND
norm_counts_de$h3_i <- NA
norm_counts_de$h3_c <- NA
norm_counts_de$h3_i.c <- NA
norm_counts_de$h9_i <- NA
norm_counts_de$h9_c <- NA
norm_counts_de$h9_i.c <- NA
norm_counts_de$h15_i <- NA
norm_counts_de$h15_c <- NA
norm_counts_de$h15_i.c <- NA

# filling in arithmetic mean of replicates counts
norm_counts_de <- norm_counts_de[match(rownames(norm_counts), norm_counts_de$id),]

for (i in 1:length(rownames(norm_counts))){
  m <- 6
  for (r1 in seq(4, 21, 2)){
    r2 <- r1 + 1
    norm_counts_de[i, m] <- (norm_counts[i, r1] + norm_counts[i, r2]) / 2
    m <- m + 1
  }
}

# remaining only DE genes IDs
norm_counts_de <- subset(norm_counts_de, norm_counts_de$id %in% de_genes)

# removing rows with starts or ends out of cromosome range
seqlen <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)

outrange <- c()
for (i in 1:length(rownames(norm_counts_de))){
  chrom <- norm_counts_de[i, 2]
  if (norm_counts_de$start[i] < 0 | norm_counts_de$end[i] > seqlen[[chrom]]){
    outrange <- c(outrange, i)
  }
}

norm_counts_de <- norm_counts_de[-outrange, ]


#make SummarizedExperiment object
rownames(norm_counts_de) <- norm_counts_de$id
norm_counts_de <- norm_counts_de[, 2:14]

n_c_de <- makeSummarizedExperimentFromDataFrame(norm_counts_de)


## match motifs
motifs <- getJasparMotifs(species = "Mus musculus")

motifs_de <- matchMotifs(motifs, n_c_de, genome = BSgenome.Mmusculus.UCSC.mm10)


## another RangedSummarizedExperiment object
# creating IRanges for GRanges()
ranges <- IRanges(start = norm_counts_de$start, end = norm_counts_de$end, width = rep(401, 4456), names = rownames(norm_counts_de))

# creating granges for SummarizedExperiment()
granges <- GRanges(seqnames = norm_counts_de$chr, ranges = ranges, strand = norm_counts_de$strand)

# creating RangedSummarizedExperiment with counts
n_c_de_counts <- SummarizedExperiment(assays=list(counts=as.matrix(norm_counts_de[, 5:13])), rowRanges=granges, colData = DataFrame(experiment = names(lrts)))


## Compute GC content for peaks
n_c_de_counts <- addGCBias(n_c_de_counts, genome = BSgenome.Mmusculus.UCSC.mm10)


## get a set of background peaks for each peak based on GC content and number 
## of fragments across all samples
bg_peaks <- getBackgroundPeaks(n_c_de_counts)


## Compute deviations in chromatin accessibility across sets of annotations
dev <- computeDeviations(object = n_c_de_counts, annotations = motifs_de, background_peaks = bg_peaks)


## tsne plot (Perplexity (3) given too high; Setting perplexity to 2)
tsne_results <- deviationsTsne(dev, threshold = 1.5, perplexity = 2)

tsne_plot <- plotDeviationsTsne(dev, tsne_results, 
                                annotation_name = NULL, 
                                sample_column = "experiment", 
                                shiny = FALSE)

pdf("tSNE.pdf", width = 10, height = 6)

plotDeviationsTsne(dev, tsne_results, 
                   annotation_name = NULL, 
                   sample_column = "experiment", 
                   shiny = FALSE)

dev.off()
