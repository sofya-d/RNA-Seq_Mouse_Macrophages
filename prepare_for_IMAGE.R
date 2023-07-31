## links
# IMAGE
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5793788/
# https://github.com/JesperGrud/IMAGE

## prepare table with gene expression data with raw counts
# structure preparetion
exp <- mm10_level2[, c(1:5, 17, 19:26, 7:16, 18)]
colnames(exp) <- c("Symbol", "Chr", "Start", "End", "Strand", "ctrl1", "ctrl2", "h3_i1", "h3_i2", "h3_c1", "h3_c2", "h3_i.c1", "h3_i.c2", "h9_i1", "h9_i2", "h9_c1", "h9_c2", "h9_i.c1", "h9_i.c2", "h15_i1", "h15_i2", "h15_c1", "h15_c2", "h15_i.c1", "h15_i.c2")

exp$RefSeq <- ""
exp$Blank_Column <- ""
exp$Gene_Length_For_RPKM_Calculation <- 1
exp <- exp[, c(26, 1:5, 27, 28, 6:25)]

# fill in "Symbol" column
for (i in 1:length(mm10_level2_annotated$GeneName)){
  if (mm10_level2_annotated$GeneName[i] == ""){
    mm10_level2_annotated$GeneName[i] = mm10_level2_annotated$CLUSTER_ID[i]
  }
}

for (i in 1:length(exp$Symbol)){
  for (j in 1:length(mm10_level2_annotated$GeneName)){
    if (exp$Symbol[i] == mm10_level2_annotated$CLUSTER_ID[j]){
      exp$Symbol[i] = mm10_level2_annotated$GeneName[j]
      break
    }
  }
}

# find indeces of duplicated gene names (except the first appearance)
dupl <- duplicated(exp$Symbol)
ind <- c()

for (i in 1:length(dupl)){
  if (dupl[i] == T){
    ind <- c(ind, i)
  }
}

# obtain names of duplicated genes
gen <- unique(exp$Symbol[ind])
duplgen <- subset(exp, exp$Symbol %in% gen)
rownames(duplgen) <- 1:length(duplgen$Symbol)


# create a list of vector containing indeces of duplicated genes (including first appearance). if duplicated gene name belongs to rows with different "Chr" value, only an index of latter appearance of gene name is included to the list. So a length of indeces of that gene name vector equals 1.
duplinds <- list()
inds <- c()
for (i in 1:length(duplgen$Symbol)){
  if (i == length(duplgen$Symbol)){
    inds <- c(inds, i)
    duplinds[[duplgen$Symbol[i]]] <- inds
    break
  }
  if (duplgen$Symbol[i] == duplgen$Symbol[i + 1]){
    inds <- c(inds, i)
  }
  if (duplgen$Symbol[i] != duplgen$Symbol[i + 1]){
    inds <- c(inds, i)
    duplinds[[duplgen$Symbol[i]]] <- inds
    inds <- c()
  }
}

# make data frame with summed counts of raws with duplicataed gene names
# prepare row structure
RefSeq <- ""
summed <- as.data.frame(RefSeq)

for (name in colnames(duplgen)[-1]){
  summed[[name]] <- NA
}

# fill in data frame
raw <- summed
for (name in names(duplinds)){
  if (length(duplinds[[name]]) > 1){
    tmp_raw <- raw
    inds <- duplinds[[name]]
    tmp_raw$Chr <- duplgen$Chr[inds[1]]
    tmp_raw$Start <- duplgen$Start[inds[1]]
    tmp_raw$End <- duplgen$End[inds[length(inds)]]
    tmp_raw$Strand <- duplgen$Strand[inds[length(inds)]]
    tmp_raw$Symbol <- name
    for (sample in colnames(duplgen)[9:length(colnames(duplgen))]){
      summed_counts <- 0
      for (i in inds){
        summed_counts <- summed_counts + duplgen[[sample]][i]
      }
      tmp_raw[[sample]] <- summed_counts
    }
    summed <- rbind(summed, tmp_raw)
  }
}

summed <- summed[-1,]
rownames(summed) <- 1:length(summed$RefSeq)
summed$Blank_Column <- ""
summed$Gene_Length_For_RPKM_Calculation <- 1

# get rid of duplications
exp <- subset(exp, !(exp$Symbol %in% gen))

# bind summed counts to the enhancer table
exp <- rbind(exp, summed)

#fill in NCBI RefSeq ID in exp$RefSeq
for (gene_exp in 1:length(exp$Symbol)){
  for (gene_ann in 1:length(mm10_level2_annotated$GeneName)){
    if (exp$Symbol[gene_exp] == mm10_level2_annotated$GeneName[gene_ann]){
      if (mm10_level2_annotated$ncbiRefSeq_50bp[gene_ann] != ""){
      exp$RefSeq[gene_exp] <- mm10_level2_annotated$ncbiRefSeq_50bp[gene_ann]
      }
      if (mm10_level2_annotated$ncbiRefSeq_500bp[gene_ann] != ""){
        exp$RefSeq[gene_exp] <- mm10_level2_annotated$ncbiRefSeq_500bp[gene_ann]
      }
      if (mm10_level2_annotated$ncbiRefSeq_1000bp[gene_ann] != ""){
        exp$RefSeq[gene_exp] <- mm10_level2_annotated$ncbiRefSeq_1000bp[gene_ann]
      }
    }
  }
}

# check if all gene names get their refseq id (except rows where value of Symbol column is cluster id)
View(subset(exp, exp$RefSeq == "" & !(grepl('L2__chr', exp$Symbol))))
View(subset(exp, exp$RefSeq == ""))


## prepare enchancer table
# structure preparation
en <- mm10_level2[, c(2, 6, 1, 37, 39:46, 27:36, 38)]
en$Start <- en$pos - 250
en$End <- en$pos + 250
en <- en[, c(1, 24, 25, 3:23)]
colnames(en) <- c("Chr", "Start", "End", "PeakID", "ctrl1", "ctrl2", "h3_i1", "h3_i2", "h3_c1", "h3_c2", "h3_i.c1", "h3_i.c2", "h9_i1", "h9_i2", "h9_c1", "h9_c2", "h9_i.c1", "h9_i.c2", "h15_i1", "h15_i2", "h15_c1", "h15_c2", "h15_i.c1", "h15_i.c2")

for (i in 1:length(en$PeakID)){
  for (j in 1:length(mm10_level2_annotated$GeneName)){
    if (en$PeakID[i] == mm10_level2_annotated$CLUSTER_ID[j]){
      en$PeakID[i] = mm10_level2_annotated$GeneName[j]
      break
    }
  }
}

# find indeces of duplicated gene names (except the first appearance)
dupl <- duplicated(en$PeakID)
ind <- c()

for (i in 1:length(dupl)){
  if (dupl[i] == T){
    ind <- c(ind, i)
  }
}

# obtain names of duplicated genes
gen <- unique(en$PeakID[ind])
duplgen <- subset(en, en$PeakID %in% gen)
rownames(duplgen) <- 1:length(duplgen$Chr)

# create a list of vector containing indeces of duplicated genes (including first appearance). if duplicated gene name belongs to rows with different "Chr" value, only an index of latter appearance of gene name is included to the list. So a length of indeces of that gene name vector equals 1.
duplinds <- list()
inds <- c()
for (i in 1:length(duplgen$PeakID)){
  if (i == length(duplgen$PeakID)){
    inds <- c(inds, i)
    duplinds[[duplgen$PeakID[i]]] <- inds
    break
  }
  if (duplgen$PeakID[i] == duplgen$PeakID[i + 1]){
    inds <- c(inds, i)
  }
  if (duplgen$PeakID[i] != duplgen$PeakID[i + 1]){
    inds <- c(inds, i)
    duplinds[[duplgen$PeakID[i]]] <- inds
    inds <- c()
  }
}

# make data frame with summed counts of raws with duplicataed gene names
# prepare row structure
Chr <- NA
summed <- as.data.frame(Chr)

for (name in colnames(duplgen)[-1]){
  summed[[name]] <- NA
}

# fill in data frame
raw <- summed
for (name in names(duplinds)){
  if (length(duplinds[[name]]) > 1){
    tmp_raw <- raw
    inds <- duplinds[[name]]
    tmp_raw$Chr <- duplgen$Chr[inds[1]]
    tmp_raw$Start <- duplgen$Start[inds[1]]
    tmp_raw$End <- duplgen$End[inds[length(inds)]]
    tmp_raw$PeakID <- name
    for (sample in colnames(duplgen)[5:length(colnames(duplgen))]){
      summed_counts <- 0
      for (i in inds){
        summed_counts <- summed_counts + duplgen[[sample]][i]
      }
      tmp_raw[[sample]] <- summed_counts
    }
    summed <- rbind(summed, tmp_raw)
  }
}

summed <- summed[-1,]
rownames(summed) <- 1:length(summed$Chr)

# get rid of duplications
en <- subset(en, !(en$PeakID %in% gen))

# bind summed counts to the enhancer table
en <- rbind(en, summed)

# count how many cluster from different chromosomes annotated as the same gene
apart <- c()
for (name in names(duplinds)){
  if (length(duplinds[[name]]) == 1){
    apart <- c(apart, name)
  }
}

for (i in c(2, 3, 5:24)){
  en[, i] <- as.numeric(format(en[, i], scientific = FALSE))
}

for (i in c(4, 5, 9:28)){
  exp[, i] <- as.numeric(format(exp[, i], scientific = FALSE))
}

options(scipen = 10000)

# write expression and enchancer tables
write.table(en, file = "enchancer.txt", sep = "\t", row.names = FALSE, col.names = F, quote = FALSE)

write.table(exp[, c(1:6, 8:28)], file = "expression.txt", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)

write.table(exp, file = "~/9_IMAGE/expression_blank.txt", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
