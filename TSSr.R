## guide https://github.com/Linlab-slu/TSSr

# set up
library(Rsamtools)
library(GenomicRanges)
library(GenomicFeatures)
library(Gviz)
library(rtracklayer)
library(DESeq2)
library(BSgenome)
library(data.table)
library(stringr)
library(TSSr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(readxl)

# create vector of sample labels
sampleLabels <- c('ctr1', 'ctr2', 'h3_i1', 'h3_i2', 'h3_c1', 'h3_c2', 'h3_i.c1', 'h3_i.c2', 'h9_i1', 'h9_i2', 'h9_c1', 'h9_c2', 'h9_i.c1', 'h9_i.c2', 'h15_i1', 'h15_i2', 'h15_c1', 'h15_c2', 'h15_i.c1', 'h15_i.c2')

# input bam files names with path

inputFiles <- c()
for (name in sampleLabels){
  inputFiles <- c(inputFiles, paste(name, ".bam", sep = ""))
}

inputFilesType <- "bam"

# Provide path and file name of genome annotation file
refSource <- "mm10.ncbiRefSeq.gtf"

# create new tssr object
myTSSr <- new("TSSr", genomeName = "BSgenome.Mmusculus.UCSC.mm10"
              ,inputFiles = inputFiles
              ,inputFilesType = inputFilesType
              ,sampleLabels = sampleLabels
              ,sampleLabelsMerged = c("Ctrl","IFN_3h", "CPG_3h", "IFN_CPG_3h", "IFN_9h", "CPG_9h", "IFN_CPG_9h", "IFN_15h", "CPG_15h", "IFN_CPG_15h")
              ,mergeIndex = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10)
              ,refSource = refSource
              ,organismName = "Mus musculus")

# The "getTSS" function identifies genomic coordinates of all TSSs and calculates read counts supporting each TSS from bam file and return values to the slot of TSSrawMatrix
getTSS(myTSSr)
# check
# myTSSr@TSSrawMatrix

# Merging samples (biological replicates)
mergeSamples(myTSSr)
# check
# myTSSr@TSSprocessedMatrix

# normalize raw counts
normalizeTSS(myTSSr)

# Clustering TSSs to infer core promoters with clusterTSS() function. Default settings
clusterTSS(myTSSr, method = "peakclu", peakDistance=100, extensionDistance=30
           , localThreshold = 0.02,clusterThreshold = 1
           , useMultiCore=FALSE, numCores=NULL)
# check
# myTSSr@tagClusters

# Aggregating consensus TSS clusters. Default settings
consensusCluster(myTSSr, dis = 50)

# Annotation (Assigning TCs to genes). Default settings
annotateCluster(myTSSr, clusters = "consensusClusters",filterCluster = TRUE,
                filterClusterThreshold = 0.02, annotationType = "genes"
                ,upstream=1000, upstreamOverlap = 500, downstream = 0)
# Annotating...
# Import genomic features from the file as a GRanges object ... OK
# Prepare the 'metadata' data frame ... OK
# Make the TxDb object ... OK
# 3 genes were dropped because they have exons located on both strands of the
# same reference sequence or on more than one reference sequence, so cannot be
# represented by a single genomic range.
# Use 'single.strand.genes.only=FALSE' to get all the genes in a GRangesList
# object, or use suppressMessages() to suppress this message.
# Ïðåäóïðåæäåíèå:
#   Â .get_cds_IDX(mcols0$type, mcols0$phase) :
#   The "phase" metadata column contains non-NA values for features of type
# stop_codon. This information was ignored.

# check
# myTSSr@assignedClusters

# import macrophage GO terms
library(readxl)
GO_term_summary_20220224_025538 <- read_excel("GO_term_summary_20220224_025538.xlsx")
macgo <- GO_term_summary_20220224_025538


# list with promoter shifts
shifts <- list()

# Core promoter shifts. Default settings
######################################
shiftPromoter(myTSSr, comparePairs = list(c("Ctrl","IFN_CPG_3h")), pval = 0.01)
shifts[["Ctrl_VS_IFN_CPG_3h"]][["all"]] <- myTSSr@PromoterShift[["Ctrl_VS_IFN_CPG_3h"]]

shiftPromoter(myTSSr, comparePairs = list(c("IFN_3h","IFN_CPG_3h")), pval = 0.01)
shifts[["IFN_3h_VS_IFN_CPG_3h"]][["all"]] <- myTSSr@PromoterShift[["IFN_3h_VS_IFN_CPG_3h"]]

shiftPromoter(myTSSr, comparePairs = list(c("CPG_3h","IFN_CPG_3h")), pval = 0.01)
shifts[["CPG_3h_VS_IFN_CPG_3h"]][["all"]] <- myTSSr@PromoterShift[["CPG_3h_VS_IFN_CPG_3h"]]

#####################################
shiftPromoter(myTSSr, comparePairs = list(c("Ctrl","IFN_CPG_9h")), pval = 0.01)
shifts[["Ctrl_VS_IFN_CPG_9h"]][["all"]] <- myTSSr@PromoterShift[["Ctrl_VS_IFN_CPG_9h"]]

shiftPromoter(myTSSr, comparePairs = list(c("IFN_9h","IFN_CPG_9h")), pval = 0.01)
shifts[["IFN_9h_VS_IFN_CPG_9h"]][["all"]] <- myTSSr@PromoterShift[["IFN_9h_VS_IFN_CPG_9h"]]

shiftPromoter(myTSSr, comparePairs = list(c("CPG_9h","IFN_CPG_9h")), pval = 0.01)
shifts[["CPG_9h_VS_IFN_CPG_9h"]][["all"]] <- myTSSr@PromoterShift[["CPG_9h_VS_IFN_CPG_9h"]]

####################################
shiftPromoter(myTSSr, comparePairs = list(c("Ctrl","IFN_CPG_15h")), pval = 0.01)
shifts[["Ctrl_VS_IFN_CPG_15h"]][["all"]] <- myTSSr@PromoterShift[["Ctrl_VS_IFN_CPG_15h"]]

shiftPromoter(myTSSr, comparePairs = list(c("IFN_15h","IFN_CPG_15h")), pval = 0.01)
shifts[["IFN_15h_VS_IFN_CPG_15h"]][["all"]] <- myTSSr@PromoterShift[["IFN_15h_VS_IFN_CPG_15h"]]

shiftPromoter(myTSSr, comparePairs = list(c("CPG_15h","IFN_CPG_15h")), pval = 0.01)
shifts[["CPG_15h_VS_IFN_CPG_15h"]][["all"]] <- myTSSr@PromoterShift[["CPG_15h_VS_IFN_CPG_15h"]]

#Calculating core promoter shifts...
#???? 50 ??? ????? ?????????????? (??????? warnings() ????? ??????????? ?????? 50)
#> warnings()
#??????????????:
#  1: In chisq.test(data[, c("tags.x", "tags.y")]) :
#  Chi-squared approximation may be incorrect

# match genes which promoters have shifts and genes associated with macrophage activation
matched <- c()
for (name in names(shifts)){
  for (i in 1:length(shifts[[name]][["all"]])){
    gene <- shifts[[name]][["all"]]$gene[i]
    if (gene %in% macgo$Symbol){
      matched <- c(matched, i)
    }
  }
  shifts[[name]][["macrophage"]] <- shifts[[name]][["all"]][matched,]
  matched <- c()
}


# ANNOTATE GENES THAT HAVE PROMOTER SHIFTS
library(biomaRt)

# see the different datasets available within a biomaRt for "dataset" argument of "useMart" function
mart = useMart("ensembl")
listDatasets(mart)
# mmusculus_gene_ensembl

# create Mart object for "mart" argument of "getGene" function
mm_ensembl = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

# possible values for "attributes" argument of "getBM" function
listAttributes()

shifts_descr <- list()
for (comparison in names(shifts)){
  genes <- shifts[[comparison]][["all"]]$gene
  descriptions <- getBM(attributes = c("external_gene_name", "description"), filters = "external_gene_name", values = genes, mart = mm_ensembl)
  shifts_descr[[comparison]] <- descriptions
}

for (comparison in names(shifts)){
  shifts[[comparison]][["all"]]$description <- ""
}

for (comparison in names(shifts)){
  inds <- match(shifts[[comparison]][["all"]]$gene, shifts_descr[[comparison]]$external_gene_name)
  for (i in 1:length(inds)){
    if (!is.na(inds[i])){
      shifts[[comparison]][["all"]]$description[i] <- shifts_descr[[comparison]]$description[inds[i]]
    }
  }
}

for (name in names(shifts)){
  write.table(shifts[[name]][["all"]], file = paste("~/tssr/promoter_shifts/", name, ".txt", sep = ""), quote = F, sep = "\t", row.names = F, col.names = T)
}
