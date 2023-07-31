# lrts_de contains genes with |logFC| > 1 and FDR < 0.05

image <- Image.venn
rm(Image.venn)
colnames(image)[1] <- "TF"

image$TF <- tolower(image$TF)

capitalize <- function(chr){
  paste(toupper(substring(chr, 1, 1)), substring(chr, 2), sep = "")
}

image$TF <- capitalize(image$TF)

ifn <- c()
cpg <- c()
ifn.cpg <- c()

for (i in 1:length(rownames(image))){
  if (image$IFN[i] == 1){
    ifn <- c(ifn, image$TF[i])
  }
  
  if (image$CPG[i] == 1){
    cpg <- c(cpg, image$TF[i])
  }
  
  if (image$CPG.IFN[i] == 1){
    ifn.cpg <- c(ifn.cpg, image$TF[i])
  }
}

TFs <- list(ifn = ifn, cpg = cpg, ifn.cpg = ifn.cpg)
DE_TFs <- list(ifn = c(), cpg = c(), ifn.cpg = c())

count = -2
de_tfs <- c()
for (group in names(TFs)){
  count = count + 3
  for (sample in names(lrts_de)[count:(count + 2)]){
    inds <- na.omit(match(TFs[[group]], lrts_de[[sample]]$GeneName))
    de_tfs <- c(de_tfs, lrts_de[[sample]]$GeneName[inds])
  }
  DE_TFs[[group]] <- unique(de_tfs)
  de_tfs <- c()
}


colnames(DE_TF_Roy) <- "TF"
DE_TF_Roy <- DE_TF_Roy$TF
DE_TF_Roy <- na.omit(DE_TF_Roy)
DE_TF_Roy <- DE_TF_Roy[-27]


DE_TFs[["ifn_Roy"]] <- DE_TF_Roy

# DRAW PLOT
library(RColorBrewer)
library(grid)
library(VennDiagram)

myCol <- brewer.pal(4, "Accent")

pdf("DE_Venn.pdf")

grid.newpage()
venn_plot <- venn.diagram(
  x = DE_TFs,
  category.names = c("Ifn", "CpG", "Ifn + CpG", "Ifn Roy et al."),
  filename = NULL,
  output=TRUE,
  # Main
  main = "DE TFs of enriched motifs including motifs from Roy et al.",
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

