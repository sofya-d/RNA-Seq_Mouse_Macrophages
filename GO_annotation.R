## CREATE VENN DIAGRAMS OF DE GENES
# PREPARE INPUT
library(edgeR)

venn <- list(venn_3h = list(), venn_9h = list(), venn_15h = list())

for (i in c(1, 4, 7)){
  for (name in names(lrts_de)[i:(i+2)]){
    venn[[names(venn)[(i+2) / 3]]][[name]] <- rownames(lrts_de[[name]])
  }
}

# DRAW PLOT
library(RColorBrewer)
library(grid)
library(VennDiagram)

myCol <- brewer.pal(3, "Pastel2")
hours <- c("3", "9", "15")

pdf("Venn_diagrams.pdf")

for(name in names(venn)){
  grid.newpage()
  venn_plot <- venn.diagram(
    x = venn[[name]],
    category.names = names(venn[[name]]),
    filename = NULL,
    output=TRUE,
    # Main
    main = paste("DE clusters, ", hours[match(name, names(venn))], " hours", sep = ''),
    main.fontface = "bold",
    main.fontfamily = "sans",
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol,
    # Numbers
    cex = 1,
    fontface = "bold",
    fontfamily = "sans",
    # Set names
    cat.cex = 1,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    rotation = 1)
  grid.draw(venn_plot)
}

dev.off()




### GO ANNOTATION
## PREPARE INPUT
# GET GENE NAMES
go_de_genename <- list(venn_3h = list(), venn_9h = list(), venn_15h = list())

for (i in c(1, 4, 7)){
  for (name in names(lrts_de)[i:(i+2)]){
    vec <- lrts_de[[name]]$GeneName
    go_de_genename[[names(go_de_genename)[(i+2) / 3]]][[name]] <- subset(vec, vec != "")
  }
}

# GET GENE ID
library(clusterProfiler)
library(org.Mm.eg.db)

go_de_entrezid <- list()

for (list in names(go_de_genename)){
  for (name in names(go_de_genename[[list]])){
    id <- bitr(go_de_genename[[list]][[name]], fromType = "SYMBOL", toType = "ENTREZID",  OrgDb = org.Mm.eg.db, drop = TRUE)
    go_de_entrezid[[list]][[name]] <- id$ENTREZID
  }
}


# GET INTERSECTIONS
library(gplots)

go_venn <- list()

for (name in names(go_de_entrezid)){
  venn_plot <- venn(go_de_entrezid[[name]], show.plot = FALSE)
  go_venn[[name]] <- attr(venn_plot, "intersections")
}




# ENRICH GO ANALYSIS
go_enrich <- list()

for (list in names(go_venn)){
  for (name in names(go_venn[[list]])[c(1:3, 7)]){
    go_enrich[[list]][[name]] <- enrichGO(gene = as.numeric(go_venn[[list]][[name]]), OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "BP")
  }
}

# SUMMARY TABLE OF ENRICH GO TEGS
sample = c()
terms_num <- c()
macr_terms_num <- c()
macr_terms <- c()

for (time in names(go_enrich)){
  for (name in names(go_enrich[[time]])){
    sample <- c(sample, name)
    
    result <- go_enrich[[time]][[name]]@result
    
    terms_num <- c(terms_num, length(rownames(result)))
    
    macr_vec <- subset(result$Description, grepl("macrophage", result$Description, fixed = TRUE))
    macr_num <- length(macr_vec)
    macr_terms_num <- c(macr_terms_num, macr_num)
    
    terms <- paste(macr_vec, collapse = "/")
    macr_terms <- c(macr_terms, terms)
  }
}

enrich_sum <- data.frame(sample)
enrich_sum$terms_num <- terms_num
enrich_sum$macr_terms_num <- macr_terms_num

per <- (macr_terms_num / terms_num) * 100
enrich_sum$macr_terms_per <- round(per, digits = 3)

enrich_sum$macr_terms <- macr_terms

rownames(enrich_sum) <- enrich_sum$sample

write.csv(enrich_sum, "enrich_sum.csv", quote = F, row.names = F)

library(heatmaply)

mx <- as.matrix(enrich_sum$macr_terms_per)
rownames(mx) <- enrich_sum$sample
colnames(mx) <- "macrophage GO terms percent"

heatmaply(mx, Rowv = F, Colv = F, column_text_angle = 0, main = "Percent of GO terms related to macrophage")




# GO SIMPLIFY (default settings)
go_simplify <- list()

for (list in names(go_enrich)){
  for (name in names(go_enrich[[list]])){
    go_simplify[[list]][[name]] <- simplify(go_enrich[[list]][[name]], cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
  }
}

# SUMMARY TABLE OF simplify GO TEGS
sample = c()
terms_num <- c()
macr_terms_num <- c()
macr_terms <- c()

for (time in names(go_simplify)){
  for (name in names(go_simplify[[time]])){
    sample <- c(sample, name)
    
    result <- go_simplify[[time]][[name]]@result
    
    terms_num <- c(terms_num, length(rownames(result)))
    
    macr_vec <- subset(result$Description, grepl("macrophage", result$Description, fixed = TRUE))
    macr_num <- length(macr_vec)
    macr_terms_num <- c(macr_terms_num, macr_num)
    
    terms <- paste(macr_vec, collapse = "/")
    macr_terms <- c(macr_terms, terms)
  }
}

simplify_sum <- data.frame(sample)
simplify_sum$terms_num <- terms_num
simplify_sum$macr_terms_num <- macr_terms_num

per <- (macr_terms_num / terms_num) * 100
simplify_sum$macr_terms_per <- round(per, digits = 3)

simplify_sum$macr_terms <- macr_terms

rownames(simplify_sum) <- simplify_sum$sample

write.csv(simplify_sum, "simplify_sum.csv", quote = F, row.names = F)

mx <- as.matrix(enrich_sum$macr_terms_per)
mx <- cbind(mx, simplify_sum$macr_terms_per)
rownames(mx) <- simplify_sum$sample
colnames(mx) <- c("enrich", "simplify")

library(heatmaply)

heatmaply(mx, Rowv = F, Colv = F, column_text_angle = 0, main = "Percent of GO terms related to macrophage")




# BARPLOT
library(grid)
library(ggplot2)

pdf("GO_barplots.pdf")

for (list in names(go_simplify)){
  for (name in names(go_simplify[[list]])){
    result <- go_simplify[[list]][[name]]@result
    
    if (nrow(result) > 20){
      num <- 20
    } else num <- nrow(result)
    
    result <- result[1:num,]
    
    y = paste(result$ID, ": ", substr(result$Description, 1, 40), sep = "")
    y = factor(y, levels = unique(y), ordered = T)
    
    plot <- ggplot(data = result, aes(y = y, x = Count, fill = p.adjust))+
      geom_col()+
      coord_cartesian()+
      # coord_flip()+
      ggtitle(name)+
      scale_fill_gradient2(low = "#4500BD", 
                           high = "#BD0084",
                           midpoint = median(result$p.adjust))+
      labs(x = "Gene count", y = "GO term")
    grid.draw(plot)
  }
}

dev.off()




# barplot of GO terms related to macrophages
# create dataframe
go_macrophage <- list()

for (list in names(go_enrich)){
  for (name in names(go_enrich[[list]])){
    result <- go_enrich[[list]][[name]]@result
    go_macrophage[[name]] <- result[grepl("macrophage", result$Description),]
  }
}

# plot
pdf("GO_macrophages_barplots.pdf")

for (name in names(go_macrophage)){
  result <- go_macrophage[[name]]
  
  y = paste(result$ID, ": ", substr(result$Description, 1, 40), sep = "")
  y = factor(y, levels = unique(y), ordered = T)
  
  plot <- ggplot(data = go_macrophage[[name]], aes(y = y, x = Count, fill = p.adjust))+
    geom_col()+
    coord_cartesian()+
    # coord_flip()+
    ggtitle(name)+
    scale_fill_gradient2(low = "#4500BD", 
                         high = "#BD0084",
                         midpoint = median(result$p.adjust))+
    labs(x = "Gene count", y = "GO term")
    
  grid.draw(plot)
}

dev.off()




# write table of GO terms
names <- c("CPG_3h", "IFN_3h", "IFNplusCPG_3h", "common_3h", "CPG_9h", "IFN_9h", "IFNplusCPG_9h", "common_9h", "CPG_15h", "IFN_15h", "IFNplusCPG_15h", "common_15h")

count <- 0
for (list in names(go_enrich)){
  for (name in names(go_enrich[[list]])){
    count <- count + 1
    file <- go_enrich[[list]][[name]]@result
    write.table(file, file = paste("~/3_GO/term_tables/", names[[count]], ".txt", sep = ''), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
}




## GO TERMS HEATMAP
library(org.Mm.eg.db)
library(clusterProfiler)
go_hmy <- list()
for (list in names(go_de_entrezid)){
  for (name in names(go_de_entrezid[[list]])){
    go_hmy[[list]][[name]] <- enrichGO(gene = as.numeric(go_de_entrezid[[list]][[name]]), OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "BP")
  }
}

# check
# dotplot(go_hmy[["venn_3h"]][["CPGvsCtrl_3h"]])
# head(go_hmy[[1]]@result)

tabs <- list()
for (venn in names(go_hmy)){
  for (sample in names(go_hmy[[venn]])){
    tabs[[sample]] <- go_hmy[[venn]][[sample]]@result[go_hmy[[venn]][[sample]]@result$p.adjust < 0.05, ]
  }
}

tab_hmy <- merge(tabs[["IFNvsCtrl_3h"]], tabs[["CPGvsCtrl_3h"]], by.x="ID", by.y="ID", all.x=T)
for (sample in names(tabs)[3:9]){
  tab_hmy <- merge(tab_hmy, tabs[[sample]], by.x="ID", by.y="ID", all.x=T)
}

tab_hmy <- na.omit(tab_hmy)
go_names <- as.character(paste(tab_hmy$ID, substr(tab_hmy$Description.x, 1, 40)))
# check
# head(go_names)

tab_hmy <- tab_hmy[,grep("p.adjust", colnames(tab_hmy))]
colnames(tab_hmy) <- c("IFNvsCtrl_3h", "CPGvsCtrl_3h", "IFN.CPGvsCtrl_3h", "IFNvsCtrl_9h", "CPGvsCtrl_9h", "IFN.CPGvsCtrl_9h", "IFNvsCtrl_15h", "CPGvsCtrl_15h", "IFN.CPGvsCtrl_15h")
rownames(tab_hmy) <- go_names

tab_hmy <- abs(log2(tab_hmy))

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

tab_hmy <- t(apply(tab_hmy, 1, cal_z_score))

library(heatmaply)
library(viridis)

colors_hm <- viridis(256, alpha = 1, begin = 0, end = 1, direction = 1, option = "viridis")
heatmaply(tab_hmy, k_row = 9, colors = colors_hm, showticklabels = c(T,F), main = "GO terms of DEG z-score(|log2(p.adjusted|))", margins = c(50, 70, 35, 0))


# heatmap
go_hmy <- list()
for (list in names(go_de_entrezid)){
  for (name in names(go_de_entrezid[[list]])){
    go_hmy[[name]] <- enrichGO(gene = as.numeric(go_de_entrezid[[list]][[name]]), OrgDb = org.Mm.eg.db, keyType = 'ENTREZID', ont = "BP")
  }
}

## FIND GO TERMS INTERSECTIONS
# CREATE LIST OF VECTORS OF GO TERMS
go_terms <- list()
for (time in names(go_hmy)){
  for (sample in names(go_hmy[[time]])){
    result <- go_hmy[[time]][[sample]]
    terms <- as.character(result$ID)
    go_terms[[time]][[sample]] <- terms
  }
}

# FIND INTERSECTIONS
library(gplots)

go_intersect <- list()
for (time in names(go_terms)){
  intersections <- venn(go_terms[[time]], show.plot = FALSE)
  go_intersect[[time]] <- attr(intersections, "intersections")
}

# SUBSTRACT GO ENRICH RESULTS
go_int <- list()
for (time in names(go_hmy)){
  for (sample in names(go_hmy[[time]])){
    result <- go_hmy[[time]][[sample]]
    terms <- go_intersect[[time]][[sample]]
    go_int[[sample]] <- subset(result, result$ID %in% terms)
  }
}





# CREATE BARPLOT
library(ggplot2)
library(grid)

pdf("unique_GO_barplots.pdf")

for (name in names(go_int)){
  result <- go_int[[name]][1:20,]
  y = paste(result$ID, ": ", substr(result$Description, 1, 40), sep = "")
  y = factor(y, levels = y, ordered = T)
  plot <- ggplot(data = result, aes(y = y, x = Count, fill = p.adjust))+
    geom_col()+
    coord_cartesian()+
    # coord_flip()+
    ggtitle(name)+
    scale_fill_gradient2(low = "#4500BD", 
                         high = "#BD0084",
                         midpoint = median(result$p.adjust))+
    labs(x = "Gene count", y = "GO term")
  grid.draw(plot)
}

dev.off()
