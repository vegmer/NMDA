DiffFR_ByBin <- function(data=Diff_DSNS_NotDown_35bins, cueExcOnly=FALSE, color=colindx[1], ymin=-2, ymax=18,
                         graphFolder=MixedGraphFolder, experiment="Exp 1 Not Down", points=FALSE, comparisons=NA){
        
        uniqueBins <- unique(data$bin)
        nBins <- 1:length(uniqueBins)
        nTrials <- uniqueBins[2]-uniqueBins[1]
        
        xmin=0; xmax=nBins+1
        
        if(cueExcOnly==TRUE){unitSel <- "CueExc Only"} else {unitSel <- "All units"}
        
        filename <- paste(graphFolder, experiment, "DIFF FR trial from CP", unitSel, ".pdf", sep="_")
        pdf(file = filename)
        
        plot.new()
        plot.window(xlim=c(xmin, max(xmax)), ylim=c(ymin, ymax))
        
        byBinDiff <- lapply(seq(1, length(uniqueBins)), function(i){
                
                perBin <- data[data$bin==uniqueBins[i], ]
                
                if(cueExcOnly==TRUE){
                        perBin <- perBin[perBin$CueExcited==TRUE, ]
                }
                
                summ <- summary(perBin$DiffFR)
                
                MIN <- summ[1]; MAX  <- summ[6]
                Q1 <- summ[2]; Q3 <- summ[5]; IQR <- Q3-Q1
                MEAN <- summ[4]; MEDIAN <- summ[3]
                
                #Plot Interquartile range with a box
                rect(xleft = i-0.3, xright = i+0.3, 
                     ybottom = Q1, ytop=Q3, col=color, border = "white")
                
                #Plot mean with a black line
                segments(x0=i-0.3, x1=i+0.3, y0=MEAN, y1=MEAN, col="black", lwd=2)
                
                #Plot median with a white line
                segments(x0=i-0.3, x1=i+0.3, y0=MEDIAN, y1=MEDIAN, col="white", lwd=2)
                
                #Plot data of individual neurons as points over the boxplot
                if(points==TRUE){
                        points(x=rep(i, nrow(perBin)), y=perBin$DiffFR, pch=19, cex=0.5)
                }
                
                perBin$DiffFR
           
        })
        
        atVals <- (1:(length(uniqueBins)+1))-0.5
        labVals <- c(uniqueBins, max(uniqueBins)+(uniqueBins[2]-uniqueBins[1]))
        
        axis(side=1, at=atVals, labels=labVals, cex.axis=1.4)
        abline(v=atVals[which(labVals==0)], lwd=2)
        axis(side=2, at=seq(ymin, ymax, by=2), cex.axis=1.4, las=2)
        abline(h=0, lty=3)
        
        dev.off()
        
        ### ANALYSIS: COMPARE NON-CONSECUTIVE BINS TO GUARANTEE THAT DIFFERENT NEURONS ARE COMPARED
        forComp1 <- byBinDiff[c(!is.even(x=nBins))] #Bins 1, 3, 5, 7, etc.
        forComp2 <- byBinDiff[c(is.even(x=nBins))] #Bins 2, 4, 6, 8, etc.
        
        Comp1 <- do.call("rbind", lapply(seq(1, length(forComp1)-1), function(n){
                w <- wilcox.test(x=forComp1[[n]], y=forComp1[[n+1]], paired=FALSE, alternative = "less")
                bin1Idx <- nBins[!is.even(nBins)][n]
                bin1Name1 <- uniqueBins[bin1Idx]
                bin1Name2 <- bin1Name1+nTrials
                
                bin2Idx <- nBins[!is.even(nBins)][n+1]
                bin2Name1 <- uniqueBins[bin2Idx]
                bin2Name2 <- bin2Name1+nTrials
                
                data.frame(idx=bin1Idx, bin=paste(bin1Idx, " vs. ", bin2Idx, sep=""), 
                           comp=paste(bin1Name1, " to ", bin1Name2, " vs. ", bin2Name1, " to ", bin2Name2, sep=""), 
                           W=w$statistic, p.val=w$p.value)
        })
        )
        
        Comp2 <- do.call("rbind", lapply(seq(1, length(forComp2)-1), function(n){
                w <- wilcox.test(x=forComp2[[n]], y=forComp2[[n+1]], paired=FALSE, alternative = "less")
                bin1Idx <- nBins[is.even(nBins)][n]
                bin1Name1 <- uniqueBins[bin1Idx]
                bin1Name2 <- bin1Name1+nTrials
                
                bin2Idx <- nBins[is.even(nBins)][n+1]
                bin2Name1 <- uniqueBins[bin2Idx]
                bin2Name2 <- bin2Name1+nTrials
                
                data.frame(idx=bin1Idx, bin=paste(bin1Idx, " vs. ", bin2Idx, sep=""), 
                           comp=paste(bin1Name1, " to ", bin1Name2, " vs. ", bin2Name1, " to ", bin2Name2, sep=""), 
                           W=w$statistic, p.val=w$p.value)
        })
        )
        
        Comp <- rbind(Comp1, Comp2)
        
        #Sort by first bin of the pair being compared
        Comp <- Comp[match((1:nrow(Comp)), Comp$idx), ]
        
        Comp$p.adj <- p.adjust(Comp$p.val, method="holm")
        
        Comp$sig <- giveStars(Comp$p.adj)
        
        results <- Comp
        
        
        #########################################
        ## JUST COMPARE THE BIN BEFORE CHANGE POINT WITH THE FIRST BIN
        # This is in case I want to define which bins I want to compare. In that case, define it in the "comparisons" parameter as a list of comparisons.
        # e.g. comparisons = list(c(1, 4), c(5, 8)) to compare bins 1 vs. 4 and 5 vs. 8
        
        if(!is.na(comparisons[1])){
                
                select.comp <- do.call("rbind", lapply(seq(1, length(comparisons)), function(x){
                        tocomp1<- unlist(byBinDiff[comparisons[[x]][1]])
                        tocomp2 <- unlist(byBinDiff[comparisons[[x]][2]])
                        
                        tocomp1Name1 <- uniqueBins[comparisons[[x]][1]]; tocomp1Name2 <- tocomp1Name1+nTrials
                        tocomp2Name1 <- uniqueBins[comparisons[[x]][2]]; tocomp2Name2 <- tocomp2Name1+nTrials
                        
                        wtest <- wilcox.test(x=tocomp1, y=tocomp2, paired=FALSE, alternative="less")
                        
                        data.frame(bins=paste(comparisons[[x]][1], " vs. ", comparisons[[x]][2], sep=""), 
                                   comparison=paste(tocomp1Name1, " to ", tocomp1Name2, " vs. ", tocomp2Name1, " to ", tocomp2Name2, sep=''),
                                   W=wtest$statistic, p.val=wtest$p.value)
                })
                )
                
                select.comp$p.adj <- p.adjust(p=select.comp$p.val, method="holm")
                select.comp$sig <- giveStars(select.comp$p.adj)
                
                results <- list(Comp, select.comp)
        }
        
        return(results)
        
}



save(DiffFR_ByBin, file="E:/Dropbox/NMDA/R Functions/DiffFR_ByBin.R")
