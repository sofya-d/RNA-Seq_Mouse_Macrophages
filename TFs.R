colnames(venn)[1] <- "TF"

# ANNOTATE GENES
library(biomaRt)

# see the different datasets available within a biomaRt for "dataset" argument of "useMart" function
mart = useMart("ensembl")
listDatasets(mart)
# mmusculus_gene_ensembl

# create Mart object for "mart" argument of "getGene" function
mm_ensembl = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

# possible values for "attributes" argument of "getBM" function
listAttributes()

TFs <- venn$TF
descriptions <- getBM(attributes = c("external_gene_name", "description"), filters = "external_gene_name", values = TFs, mart = mm_ensembl)

venn$description <- ""

inds <- match(venn$TF, toupper(descriptions$external_gene_name))
for (i in 1:length(inds)){
  if (!is.na(inds[i])){
    venn$description[i] <- descriptions$description[inds[i]]
  }
}


# VENN ARRAYS SEPARATION
venn$TF <- tolower(venn$TF)

capitalize <- function(chr){
  paste(toupper(substring(chr, 1, 1)), substring(chr, 2), sep = "")
}

venn$TF <- capitalize(venn$TF)

id <- bitr(venn$TF, fromType = "SYMBOL", toType = "ENTREZID",  OrgDb = org.Mm.eg.db, drop = TRUE)

venn$ENTREZID <- id$ENTREZID
venn <- venn[, c(1, 6, 2:5)]

ifn <- c()
cpg <- c()
ifn.cpg <- c()

for (i in 1:length(rownames(venn))){
  if (venn$IFN[i] == 1){
    ifn <- c(ifn, venn$ENTREZID[i])
  }
  
  if (venn$CPG[i] == 1){
    cpg <- c(cpg, venn$ENTREZID[i])
  }
  
  if (venn$CPG.IFN[i] == 1){
    ifn.cpg <- c(ifn.cpg, venn$ENTREZID[i])
  }
}

IDs <- list(ifn = ifn, cpg = cpg, ifn.cpg = ifn.cpg)


# GO ANNOTATION
library(clusterProfiler)
library(org.Mm.eg.db)

go <- list()

for (name in names(IDs)){
  ids <- IDs[[name]]
  enrich <- enrichGO(gene = ids, OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "BP")
  # default settings
  simplify <- simplify(enrich, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
  go[[name]] <- simplify@result
}

for (name in names(go)){
  rownames(go[[name]]) <- 1:length(rownames(go[[name]]))
}

go[["ifn.cpg"]]$geneID[203]
venn$description[16]


# COMPARISON WITH MOTIFS FROM ARTICLE
ifn <- c()
cpg <- c()
ifn.cpg <- c()

for (i in 1:length(rownames(venn))){
  if (venn$IFN[i] == 1){
    ifn <- c(ifn, venn$TF[i])
  }
  
  if (venn$CPG[i] == 1){
    cpg <- c(cpg, venn$TF[i])
  }
  
  if (venn$CPG.IFN[i] == 1){
    ifn.cpg <- c(ifn.cpg, venn$TF[i])
  }
}

colnames(art_ifn) <- "tf"
art.ifn <- art_ifn$tf
art.ifn <- tolower(art.ifn)
art.ifn <- capitalize(art.ifn)

TFs <- list(ifn = ifn, cpg = cpg, ifn.cpg = ifn.cpg, art.ifn = art.ifn)

# DRAW PLOT
library(RColorBrewer)
library(grid)
library(VennDiagram)

myCol <- brewer.pal(4, "Accent")

pdf("Venn.pdf")

grid.newpage()
venn_plot <- venn.diagram(
  x = TFs,
  category.names = c("Ifn", "CpG", "Ifn + CpG", "Ifn Roy et al."),
  filename = NULL,
  output=TRUE,
  # Main
  main = "TFs of enriched motifs including motifs from Roy et al.",
  main.cex = 1.5,
  main.fontface = "bold",
  main.fontfamily = "sans",
  # Circles
  fill = myCol,
  # Numbers
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.fontfamily = "sans")
  grid.draw(venn_plot)

dev.off()

# GET INTERSECTIONS
library(gplots)

go_venn <- list()

for (name in names(go_de_entrezid)){
  venn_plot <- venn(go_de_entrezid[[name]], show.plot = FALSE)
  go_venn[[name]] <- attr(venn_plot, "intersections")
}
