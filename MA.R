# setup
library(RColorBrewer)
library(wordcloud)
library(scales)

regulation <- c("up", "down")
up_down <- as.data.frame(regulation)

#WRITING PLOT
pdf("MA_plots.pdf", width = 6, height = 6)
#LOOP
for (i in names(cont)){
        lrt <- glmLRT(fit, contrast = cont[[i]])
        
        #òàáëèöà äëÿ ãðàôèêà
        sorted <- topTags(lrt, n=nrow(lrt$table))
        final <- sorted$table
        final$cluster_ID <- rownames(final)
        final <- final[, c(7, c(1:6))]
        #logPval
        final$logPval <- abs(log10(final$PValue))
        #colour
        final$Colour<-"black"
        final$Colour[final$FDR<0.05 & final$logFC> 1] <- brewer.pal(8, "Paired")[2]
        final$Colour[final$FDR<0.05 & final$logFC< -1] <- brewer.pal(8, "Paired")[4]
        #symbol
        final$cex <- 0.5
        final$cex[final$FDR<0.05 & abs(final$logFC)>1] <- 0.2
        #transparancy
        final$alpha<-0.01
        final$alpha[final$FDR<0.05 & abs(final$logFC)>1]<-0.2
        
        write.table(final, file = paste(i,".txt", sep = ''), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        
        #GRAPHICAL PARAMETERS
        par(mfrow=c(1,1),mai=c(1,1, 1, 1), oma=c(.1,0.1,0.1,.1), mgp = c(0, 1, 0))
        #FOR PLOT FUNCTION
        len <- abs(min(final$logFC))+abs(max(final$logFC))
        add <- len*10/100
        #FOR LEGEND FUNCTION
        up <- nrow(subset(final, logFC > 1 & FDR < 0.05))
        down <- nrow(subset(final, logFC < -1 & FDR < 0.05))
        
        up_down[[i]] <- c(up, down)
        
        #PLOT
        plot(final$logFC~c(final$logCPM), xlab="", ylim=c(min(final$logFC)-add, max(final$logFC)+add), ylab="", bg=alpha(final$Colour, final$alpha/2), col=alpha(final$Colour,final$alpha), pch=21, cex.axis = 1, cex.lab= 1, cex=final$cex, main = i)
        #LINES ON PLOT
        abline(h=c(-1,1), col="gold2", lty=2)
        #LEGEND
        legend("topleft", legend=c(paste("Up regulated: ", up, sep = "")), ncol=1, cex = 1, col = "gray80", bty = "n")
        legend("bottomleft", legend=c(paste("Down regulated: ", down, sep = "")), ncol=1, cex = 1, col = "gray80", bty = "n")
        #PLOT BOX
        box()
        mtext('logFC', side = 2, line = 2, cex = 1)
        mtext('logTPM', side = 1, line = 2, cex = 1)
}

#CLOSE PLOT
dev.off()

rownames(up_down) <- up_down$regulation
up_down <- as.matrix(up_down[, 2:10])

pheatmap(mat = up_down, kmeans_k = NA, border_color = NA, cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = F, cluster_cols = F, legend = T, show_rownames = T, show_colnames = T, main = "Number of up and down regulated DEC", angle_col = 315, display_numbers = TRUE, fontsize_number = 12)
